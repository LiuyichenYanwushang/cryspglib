#!/usr/bin/env python3
"""Derive the exact ISO-IR direct-coordinate data-Hall frames.

This is deliberately a pure in-memory derivation.  It identifies the Hall
setting already represented by the pinned PIR/CIR source coordinates, whose
affine frame is the singleton ``P = I, p = 0``.  It is not a general affine
setting reconstruction: there is no basis/origin search, 24-grid search,
score, tie-break, label matching, or runtime-data fallback, and neither ISO-IR
header declares a Hall number.

The only production inputs are the public exact ISO-IR source loader and the
public committed spglib magnetic-provenance loader.  The centering generators
below are transcribed from the official ``centeringmatrix`` and
``centeringmatrix_denom`` DATA statements in both pinned Fortran source
members.  The Fortran reader selects one of the seven structures from the
first symbol character; exact additive closure modulo conventional integer
translations produces the cosets used here.  Tests independently read both
members and verify that provenance.

This module intentionally does not write an artifact.  It reports the raw
candidate set, the centering-filtered Hall, ordered source/Hall maps, and the
full derivation census so a later sidecar step can freeze those values.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from typing import Optional, Tuple

try:  # ``scripts`` is normally imported as a namespace package.
    from . import iso_irrep_exact
    from . import spglib_magnetic_provenance
except ImportError:  # pragma: no cover - useful when run from scripts/.
    import iso_irrep_exact
    import spglib_magnetic_provenance


TRANSLATION_DENOMINATOR = 12
IDENTITY_ROTATION = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
ZERO3 = (Fraction(0), Fraction(0), Fraction(0))

EXPECTED_PIR_RECORDS = 10_294
EXPECTED_CIR_RECORDS = 11_202
EXPECTED_SOURCE_REPRESENTATIVES = 2_609
EXPECTED_RAW_UNIQUE = 220
EXPECTED_RAW_AMBIGUOUS = 10
EXPECTED_RAW_MISSING = 0
EXPECTED_RAW_AMBIGUOUS_SPACEGROUPS = (
    5, 8, 9, 12, 15, 21, 38, 39, 65, 67,
)
EXPECTED_FILTERED_UNIQUE = 230
EXPECTED_FILTERED_AMBIGUOUS = 0
EXPECTED_FILTERED_MISSING = 0
EXPECTED_SELECTED_HALL_OPERATIONS = 4_425
EXPECTED_SOURCE_TO_HALL = 2_609
EXPECTED_SOURCE_TO_HALL_NONZERO = 0
EXPECTED_HALL_TO_SOURCE = 4_425
EXPECTED_HALL_TO_SOURCE_NONZERO = 1_816
EXPECTED_EXPANDED_NORMALIZATION_NONZERO = 410
EXPECTED_CENTERING_COUNTS = (
    ("P", 149), ("A", 4), ("B", 0), ("C", 16),
    ("F", 16), ("I", 38), ("R", 7),
)


class IsoIrrepDataHallError(ValueError):
    """Base class for typed exact data-Hall derivation failures."""


class DataHallDerivationError(IsoIrrepDataHallError):
    """A unique direct-coordinate Hall frame cannot be derived."""


class DataHallInvariantError(IsoIrrepDataHallError):
    """A source, operation, mapping, or census invariant failed."""


class DataHallLookupError(IsoIrrepDataHallError):
    """A public data-Hall frame lookup has an invalid target."""


class DataHallSchemaError(IsoIrrepDataHallError):
    """A synthetic typed input has an invalid shape or scalar type."""


Int3 = Tuple[int, int, int]
Fraction3 = Tuple[Fraction, Fraction, Fraction]
Rotation3 = Tuple[Int3, Int3, Int3]
Distribution = Tuple[Tuple[Int3, int], ...]


# Official Fortran values from the multiplication-table routine.  Fortran's
# DATA initialization of centeringmatrix(3,3,7) is column-major: each group of
# nine values contains three translation-generator columns.
CENTERING_MATRIX_DATA = (
    ("P", 1, (1, 0, 0, 0, 1, 0, 0, 0, 1)),
    ("A", 2, (2, 0, 0, 0, 1, -1, 0, 1, 1)),
    ("B", 2, (1, 0, -1, 0, 2, 0, 1, 0, 1)),
    ("C", 2, (1, -1, 0, 1, 1, 0, 0, 0, 2)),
    ("F", 2, (0, 1, 1, 1, 0, 1, 1, 1, 0)),
    ("I", 2, (-1, 1, 1, 1, -1, 1, 1, 1, -1)),
    ("R", 3, (2, 1, 1, -1, 1, 1, -1, -2, 1)),
)


def _mod1(value: Fraction) -> Fraction:
    return value - (value.numerator // value.denominator)


def _mod1_vector(vector: Fraction3) -> Fraction3:
    return tuple(_mod1(value) for value in vector)  # type: ignore[return-value]


def _centering_cosets(
    matrix_values: Tuple[int, ...], denominator: int
) -> Tuple[Fraction3, ...]:
    if type(denominator) is not int or denominator <= 0:
        raise DataHallInvariantError("centering denominator is invalid")
    if type(matrix_values) is not tuple or len(matrix_values) != 9:
        raise DataHallInvariantError("centering matrix must contain nine integers")
    if any(type(value) is not int for value in matrix_values):
        raise DataHallInvariantError("centering matrix entries must be integers")
    generators = tuple(
        tuple(Fraction(matrix_values[3 * column + row], denominator)
              for row in range(3))
        for column in range(3)
    )
    cosets = {(Fraction(0), Fraction(0), Fraction(0))}
    frontier = list(cosets)
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = _mod1_vector(tuple(
                current[index] + generator[index] for index in range(3)
            ))
            if candidate not in cosets:
                cosets.add(candidate)
                frontier.append(candidate)
    return tuple(sorted(cosets))


CENTERING_COSETS = tuple(
    (name, _centering_cosets(matrix, denominator))
    for name, denominator, matrix in CENTERING_MATRIX_DATA
)


def _centering_cosets_for(name: str) -> Tuple[Fraction3, ...]:
    if type(name) is not str:
        raise DataHallSchemaError("centering name must be an exact string")
    for centering, cosets in CENTERING_COSETS:
        if centering == name:
            return cosets
    raise DataHallInvariantError(f"unknown centering {name!r}")


def _require_exact_int(value, context: str) -> int:
    if type(value) is not int:
        raise DataHallSchemaError(f"{context} must be an exact integer")
    return value


def _require_int3(value, context: str) -> Int3:
    if type(value) is not tuple or len(value) != 3:
        raise DataHallSchemaError(f"{context} must be a three-integer tuple")
    if any(type(item) is not int for item in value):
        raise DataHallSchemaError(f"{context} must contain exact integers")
    return value  # type: ignore[return-value]


def _require_fraction3(value, context: str) -> Fraction3:
    if type(value) is not tuple or len(value) != 3:
        raise DataHallSchemaError(f"{context} must be a three-Fraction tuple")
    if any(type(item) is not Fraction for item in value):
        raise DataHallSchemaError(f"{context} must contain exact Fractions")
    return value  # type: ignore[return-value]


def _validate_rotation(value, context: str) -> Rotation3:
    if type(value) is not tuple or len(value) != 3:
        raise DataHallSchemaError(f"{context} must be a 3x3 tuple")
    if any(type(row) is not tuple or len(row) != 3 for row in value):
        raise DataHallSchemaError(f"{context} must be a 3x3 tuple")
    if any(type(item) is not int for row in value for item in row):
        raise DataHallSchemaError(f"{context} entries must be exact integers")
    determinant = _rotation_determinant(value)
    if determinant not in (-1, 1):
        raise DataHallInvariantError(f"{context} is not GL(3,Z)")
    return value  # type: ignore[return-value]


def _rotation_determinant(rotation: Rotation3) -> int:
    a, b, c = rotation
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _rotation_product(left: Rotation3, right: Rotation3) -> Rotation3:
    return tuple(
        tuple(sum(left[row][index] * right[index][column] for index in range(3))
              for column in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def _fraction_difference_is_integer(left: Fraction3, right: Fraction3) -> bool:
    return all((left[index] - right[index]).denominator == 1 for index in range(3))


def _numerators_over_12(vector: Fraction3, context: str) -> Int3:
    values = []
    for index, value in enumerate(_require_fraction3(vector, context)):
        scaled = value * TRANSLATION_DENOMINATOR
        if scaled.denominator != 1:
            raise DataHallDerivationError(
                f"{context}[{index}] cannot be represented exactly over 12"
            )
        values.append(scaled.numerator)
    return tuple(values)  # type: ignore[return-value]


def _distribution(values: Counter) -> Distribution:
    return tuple(sorted((tuple(key), count) for key, count in values.items()))


@dataclass(frozen=True)
class SourceToHall:
    """Ordered source-operation to selected-Hall-operation mapping."""

    __slots__ = ("source_operation_index", "hall_operation_index", "shift_numerator")

    source_operation_index: int
    hall_operation_index: int
    shift_numerator: Int3

    def __post_init__(self):
        source_index = _require_exact_int(
            self.source_operation_index, "SourceToHall.source_operation_index"
        )
        hall_index = _require_exact_int(
            self.hall_operation_index, "SourceToHall.hall_operation_index"
        )
        if source_index < 0 or hall_index < 0:
            raise DataHallInvariantError("SourceToHall operation index is negative")
        _require_int3(self.shift_numerator, "SourceToHall.shift_numerator")


@dataclass(frozen=True)
class HallToSource:
    """Ordered selected-Hall-operation to source-operation mapping."""

    __slots__ = ("hall_operation_index", "source_operation_index", "shift_numerator")

    hall_operation_index: int
    source_operation_index: int
    shift_numerator: Int3

    def __post_init__(self):
        hall_index = _require_exact_int(
            self.hall_operation_index, "HallToSource.hall_operation_index"
        )
        source_index = _require_exact_int(
            self.source_operation_index, "HallToSource.source_operation_index"
        )
        if hall_index < 0 or source_index < 0:
            raise DataHallInvariantError("HallToSource operation index is negative")
        _require_int3(self.shift_numerator, "HallToSource.shift_numerator")


def _validate_symbol(value, centering: str, context: str) -> str:
    if type(value) is not str or not value:
        raise DataHallSchemaError(f"{context} must be a nonempty string")
    if any(not 0x20 <= ord(character) <= 0x7E for character in value):
        raise DataHallSchemaError(f"{context} must be printable ASCII")
    if value[0] != centering:
        raise DataHallInvariantError(f"{context} does not select {centering!r}")
    return value


@dataclass(frozen=True)
class ExactDataHallFrame:
    """One SG's exact source-coordinate Hall frame and ordered maps."""

    __slots__ = (
        "spacegroup", "source_symbol", "centering", "raw_candidate_halls",
        "data_hall", "source_operation_count", "hall_operation_count",
        "source_to_hall", "hall_to_source",
    )

    spacegroup: int
    source_symbol: str
    centering: str
    raw_candidate_halls: Tuple[int, ...]
    data_hall: int
    source_operation_count: int
    hall_operation_count: int
    source_to_hall: Tuple[SourceToHall, ...]
    hall_to_source: Tuple[HallToSource, ...]

    def __post_init__(self):
        spacegroup = _require_exact_int(self.spacegroup, "ExactDataHallFrame.spacegroup")
        if not 1 <= spacegroup <= 230:
            raise DataHallInvariantError("ExactDataHallFrame.spacegroup outside 1..230")
        if type(self.centering) is not str or self.centering not in {
            name for name, _ in CENTERING_COSETS
        }:
            raise DataHallSchemaError("ExactDataHallFrame.centering is unknown")
        _validate_symbol(self.source_symbol, self.centering, "ExactDataHallFrame.source_symbol")
        if type(self.raw_candidate_halls) is not tuple or not self.raw_candidate_halls:
            raise DataHallSchemaError("ExactDataHallFrame.raw_candidate_halls is malformed")
        if any(type(hall) is not int or not 1 <= hall <= 530
               for hall in self.raw_candidate_halls):
            raise DataHallSchemaError("ExactDataHallFrame.raw_candidate_halls is malformed")
        if tuple(sorted(set(self.raw_candidate_halls))) != self.raw_candidate_halls:
            raise DataHallInvariantError("ExactDataHallFrame.raw_candidate_halls is not ordered")
        data_hall = _require_exact_int(self.data_hall, "ExactDataHallFrame.data_hall")
        if data_hall not in self.raw_candidate_halls:
            raise DataHallInvariantError("selected Hall is not a raw candidate")
        source_count = _require_exact_int(
            self.source_operation_count, "ExactDataHallFrame.source_operation_count"
        )
        hall_count = _require_exact_int(
            self.hall_operation_count, "ExactDataHallFrame.hall_operation_count"
        )
        if source_count <= 0 or hall_count <= 0:
            raise DataHallInvariantError("operation counts must be positive")
        if type(self.source_to_hall) is not tuple or len(self.source_to_hall) != source_count:
            raise DataHallInvariantError("source_to_hall cardinality mismatch")
        if type(self.hall_to_source) is not tuple or len(self.hall_to_source) != hall_count:
            raise DataHallInvariantError("hall_to_source cardinality mismatch")
        if any(not isinstance(mapping, SourceToHall)
               for mapping in self.source_to_hall):
            raise DataHallSchemaError("source_to_hall contains a wrong type")
        if any(not isinstance(mapping, HallToSource)
               for mapping in self.hall_to_source):
            raise DataHallSchemaError("hall_to_source contains a wrong type")
        if any(mapping.source_operation_index != index
               for index, mapping in enumerate(self.source_to_hall)):
            raise DataHallInvariantError("source_to_hall is not source ordered")
        if any(mapping.hall_operation_index != index
               for index, mapping in enumerate(self.hall_to_source)):
            raise DataHallInvariantError("hall_to_source is not Hall ordered")
        if any(mapping.hall_operation_index >= hall_count
               for mapping in self.source_to_hall):
            raise DataHallInvariantError("source_to_hall Hall index is out of range")
        if len({mapping.hall_operation_index for mapping in self.source_to_hall}) != source_count:
            raise DataHallInvariantError("source_to_hall Hall indices are not unique")
        if any(mapping.source_operation_index >= source_count
               for mapping in self.hall_to_source):
            raise DataHallInvariantError("hall_to_source source index is out of range")
        if len({mapping.source_operation_index for mapping in self.hall_to_source}) < source_count:
            raise DataHallInvariantError("hall_to_source does not cover source operations")
        for mapping in self.source_to_hall:
            inverse = self.hall_to_source[mapping.hall_operation_index]
            if inverse.source_operation_index != mapping.source_operation_index:
                raise DataHallInvariantError("source/Hall mapping direction disagrees")


@dataclass(frozen=True)
class DerivationCensus:
    """Immutable aggregate census, including both distinct shift metrics."""

    __slots__ = (
        "pir_records", "cir_records", "source_representatives", "raw_unique",
        "raw_ambiguous", "raw_missing", "raw_ambiguous_spacegroups",
        "filtered_unique", "filtered_ambiguous", "filtered_missing",
        "selected_hall_operations", "source_to_hall", "source_to_hall_nonzero",
        "hall_to_source", "hall_to_source_nonzero", "hall_to_source_shifts",
        "hall_to_source_cosets", "expanded_normalization_nonzero",
        "expanded_normalization_shifts", "centering_counts",
    )

    pir_records: int
    cir_records: int
    source_representatives: int
    raw_unique: int
    raw_ambiguous: int
    raw_missing: int
    raw_ambiguous_spacegroups: Tuple[int, ...]
    filtered_unique: int
    filtered_ambiguous: int
    filtered_missing: int
    selected_hall_operations: int
    source_to_hall: int
    source_to_hall_nonzero: int
    hall_to_source: int
    hall_to_source_nonzero: int
    hall_to_source_shifts: Distribution
    hall_to_source_cosets: Distribution
    expanded_normalization_nonzero: int
    expanded_normalization_shifts: Distribution
    centering_counts: Tuple[Tuple[str, int], ...]

    def __post_init__(self):
        integer_fields = (
            "pir_records", "cir_records", "source_representatives", "raw_unique",
            "raw_ambiguous", "raw_missing", "filtered_unique",
            "filtered_ambiguous", "filtered_missing", "selected_hall_operations",
            "source_to_hall", "source_to_hall_nonzero", "hall_to_source",
            "hall_to_source_nonzero", "expanded_normalization_nonzero",
        )
        for field in integer_fields:
            value = _require_exact_int(getattr(self, field), f"DerivationCensus.{field}")
            if value < 0:
                raise DataHallInvariantError(f"DerivationCensus.{field} is negative")
        if type(self.raw_ambiguous_spacegroups) is not tuple:
            raise DataHallSchemaError("raw ambiguous SG census must be a tuple")
        if any(type(value) is not int or not 1 <= value <= 230
               for value in self.raw_ambiguous_spacegroups):
            raise DataHallSchemaError("raw ambiguous SG census is malformed")
        if tuple(sorted(set(self.raw_ambiguous_spacegroups))) != self.raw_ambiguous_spacegroups:
            raise DataHallInvariantError("raw ambiguous SG census is not ordered")
        for field in (
            "hall_to_source_shifts", "hall_to_source_cosets",
            "expanded_normalization_shifts",
        ):
            value = getattr(self, field)
            if type(value) is not tuple:
                raise DataHallSchemaError(f"DerivationCensus.{field} must be a tuple")
            if tuple(sorted(value)) != value:
                raise DataHallInvariantError(f"DerivationCensus.{field} is not ordered")
            if any(
                type(row) is not tuple or len(row) != 2
                or type(row[0]) is not tuple or len(row[0]) != 3
                or any(type(component) is not int for component in row[0])
                or type(row[1]) is not int or row[1] <= 0
                for row in value
            ):
                raise DataHallSchemaError(f"DerivationCensus.{field} is malformed")
        if type(self.centering_counts) is not tuple:
            raise DataHallSchemaError("centering census must be a tuple")
        if any(
            type(row) is not tuple or len(row) != 2 or type(row[0]) is not str
            or type(row[1]) is not int or row[1] < 0
            for row in self.centering_counts
        ):
            raise DataHallSchemaError("centering census is malformed")


@dataclass(frozen=True)
class ExactDataHallDatabase:
    """Immutable ordered result returned by :func:`derive_data_hall_frames`."""

    __slots__ = ("frames", "census")

    frames: Tuple[ExactDataHallFrame, ...]
    census: DerivationCensus

    def __post_init__(self):
        if type(self.frames) is not tuple or len(self.frames) != 230:
            raise DataHallInvariantError("data-Hall result must contain 230 frames")
        if any(not isinstance(frame, ExactDataHallFrame)
               for frame in self.frames):
            raise DataHallSchemaError("data-Hall result contains a wrong frame type")
        if any(frame.spacegroup != index for index, frame in enumerate(self.frames, 1)):
            raise DataHallInvariantError("data-Hall frames are not ordered 1..230")
        if not isinstance(self.census, DerivationCensus):
            raise DataHallSchemaError("data-Hall result census has a wrong type")

    def source_frame(self, spacegroup: int) -> ExactDataHallFrame:
        if type(spacegroup) is not int or not 1 <= spacegroup <= 230:
            raise DataHallLookupError("spacegroup must be an exact int in 1..230")
        return self.frames[spacegroup - 1]

    # These names make the immutable result convenient to consume without
    # exposing any alternate mutable representation.
    @property
    def spacegroups(self) -> Tuple[ExactDataHallFrame, ...]:
        return self.frames

    def source_universe(self, spacegroup: int) -> ExactDataHallFrame:
        return self.source_frame(spacegroup)

    def __iter__(self):
        """Permit ``frames, census = derive_data_hall_frames()`` ergonomics."""

        yield self.frames
        yield self.census


# Names used by an earlier draft are harmless aliases, but the canonical API
# for this pure stage is the shorter pair above.
SourceToHallMapping = SourceToHall
HallToSourceMapping = HallToSource


def _source_spacegroup(source) -> int:
    value = _require_exact_int(source.spacegroup, "source.spacegroup")
    if not 1 <= value <= 230:
        raise DataHallInvariantError("source.spacegroup outside 1..230")
    return value


def _source_centering(source) -> Tuple[str, str]:
    try:
        symbol = source.space_group_symbol
        centering = source.centering.value
    except Exception as error:
        raise DataHallSchemaError("source lacks symbol or centering") from error
    if type(centering) is not str:
        raise DataHallSchemaError("source centering is not an exact name")
    _centering_cosets_for(centering)
    _validate_symbol(symbol, centering, "source.space_group_symbol")
    return symbol, centering


def _source_operations(source):
    try:
        operations = source.operations
    except Exception as error:
        raise DataHallSchemaError("source lacks operations") from error
    if type(operations) is not tuple or not operations:
        raise DataHallSchemaError("source operations must be a nonempty tuple")
    checked = []
    for index, operation in enumerate(operations):
        try:
            rotation = _validate_rotation(operation.rotation, f"source operation {index} rotation")
            translation = _require_fraction3(
                operation.translation, f"source operation {index} translation"
            )
        except AttributeError as error:
            raise DataHallSchemaError("source operation lacks exact fields") from error
        checked.append((rotation, translation))
    rotations = tuple(item[0] for item in checked)
    if len(set(rotations)) != len(rotations):
        raise DataHallInvariantError("source operation rotations are not unique")
    if rotations[0] != IDENTITY_ROTATION or checked[0][1] != ZERO3:
        raise DataHallInvariantError("source operation slot 0 is not exact identity")
    return tuple(checked)


def _hall_operations(provenance, hall: int):
    try:
        operations = tuple(provenance.spg_operations(hall))
    except Exception as error:
        raise DataHallDerivationError(f"unable to read SPG Hall {hall}") from error
    if not operations:
        raise DataHallDerivationError(f"SPG Hall {hall} has no operations")
    checked = []
    for index, operation in enumerate(operations):
        try:
            rotation = _validate_rotation(operation.rotation, f"Hall {hall} operation {index} rotation")
            translation = _require_fraction3(
                operation.translation, f"Hall {hall} operation {index} translation"
            )
        except AttributeError as error:
            raise DataHallSchemaError("Hall operation lacks exact fields") from error
        checked.append((rotation, translation))
    return tuple(checked)


def _raw_candidates(source, provenance):
    """Return every raw direct-coordinate candidate, preserving Hall order."""

    spacegroup = _source_spacegroup(source)
    source_operations = _source_operations(source)
    source_rotation_set = {rotation for rotation, _ in source_operations}
    candidates = []
    details = {}
    for hall in range(1, 531):
        try:
            parent = provenance.spacegroup_number_for_hall(hall)
        except Exception as error:
            raise DataHallDerivationError(
                f"unable to resolve parent SG for Hall {hall}"
            ) from error
        if type(parent) is not int:
            raise DataHallSchemaError(f"Hall {hall} parent spacegroup is not an int")
        if parent != spacegroup:
            continue
        hall_operations = _hall_operations(provenance, hall)
        if {rotation for rotation, _ in hall_operations} != source_rotation_set:
            continue
        source_matches = []
        for source_index, (source_rotation, source_translation) in enumerate(source_operations):
            matches = [
                hall_index
                for hall_index, (hall_rotation, hall_translation) in enumerate(hall_operations)
                if hall_rotation == source_rotation
                and _fraction_difference_is_integer(source_translation, hall_translation)
            ]
            if len(matches) > 1:
                raise DataHallDerivationError(
                    f"Hall {hall} has duplicate same-rotation match for source operation {source_index}"
                )
            if len(matches) != 1:
                source_matches = []
                break
            source_matches.append(matches[0])
        if len(source_matches) == len(source_operations):
            candidates.append(hall)
            details[hall] = (hall_operations, tuple(source_matches))
    return tuple(candidates), tuple((hall, details[hall]) for hall in candidates)


def _details_for(details, hall: int):
    for candidate, value in details:
        if candidate == hall:
            return value
    raise DataHallDerivationError(f"missing raw Hall details for Hall {hall}")


def _pure_translation_cosets(hall_operations) -> Tuple[Fraction3, ...]:
    return tuple(sorted({
        _mod1_vector(translation)
        for rotation, translation in hall_operations
        if rotation == IDENTITY_ROTATION
    }))


def _filtered_candidates(candidates, details, centering_cosets):
    if type(centering_cosets) is not tuple or not centering_cosets:
        raise DataHallDerivationError("centering cosets must be a nonempty tuple")
    expected = tuple(sorted(centering_cosets))
    return tuple(
        hall for hall in candidates
        if _pure_translation_cosets(_details_for(details, hall)[0]) == expected
    )


def _mapping_for_selected_hall(source_operations, hall_operations, centering_cosets, context):
    source_to_hall = []
    used_hall_indices = set()
    for source_index, (source_rotation, source_translation) in enumerate(source_operations):
        matches = [
            hall_index
            for hall_index, (hall_rotation, hall_translation) in enumerate(hall_operations)
            if hall_rotation == source_rotation
            and _fraction_difference_is_integer(source_translation, hall_translation)
        ]
        if len(matches) != 1:
            raise DataHallDerivationError(f"{context} source-to-Hall mapping is not unique")
        hall_index = matches[0]
        if hall_index in used_hall_indices:
            raise DataHallDerivationError(f"{context} source-to-Hall Hall index repeats")
        used_hall_indices.add(hall_index)
        hall_translation = hall_operations[hall_index][1]
        shift = tuple(
            source_translation[index] - hall_translation[index]
            for index in range(3)
        )
        source_to_hall.append(SourceToHall(
            source_operation_index=source_index,
            hall_operation_index=hall_index,
            shift_numerator=_numerators_over_12(shift, f"{context} source-to-Hall shift"),
        ))

    source_by_rotation = {
        rotation: index for index, (rotation, _) in enumerate(source_operations)
    }
    hall_to_source = []
    for hall_index, (hall_rotation, hall_translation) in enumerate(hall_operations):
        source_index = source_by_rotation.get(hall_rotation)
        if source_index is None:
            raise DataHallDerivationError(f"{context} Hall-to-source rotation is missing")
        source_translation = source_operations[source_index][1]
        shift = tuple(
            hall_translation[index] - source_translation[index]
            for index in range(3)
        )
        shift_numerator = _numerators_over_12(shift, f"{context} Hall-to-source shift")
        if _mod1_vector(shift) not in centering_cosets:
            raise DataHallInvariantError(
                f"{context} Hall-to-source shift is outside centering cosets"
            )
        hall_to_source.append(HallToSource(
            hall_operation_index=hall_index,
            source_operation_index=source_index,
            shift_numerator=shift_numerator,
        ))
    if len(hall_to_source) != len(hall_operations):
        raise DataHallDerivationError(f"{context} Hall-to-source mapping is incomplete")
    return tuple(source_to_hall), tuple(hall_to_source)


def _source_quotient_check(source_operations, centering_cosets, context: str) -> None:
    rotation_index = {rotation: index for index, (rotation, _) in enumerate(source_operations)}
    for left_index, (left_rotation, left_translation) in enumerate(source_operations):
        for right_index, (right_rotation, right_translation) in enumerate(source_operations):
            product_rotation = _rotation_product(left_rotation, right_rotation)
            product_index = rotation_index.get(product_rotation)
            if product_index is None:
                raise DataHallInvariantError(
                    f"{context} source quotient rotation product is missing"
                )
            target_translation = source_operations[product_index][1]
            delta = tuple(
                left_translation[row]
                + sum(left_rotation[row][column] * right_translation[column]
                      for column in range(3))
                - target_translation[row]
                for row in range(3)
            )
            if _mod1_vector(delta) not in centering_cosets:
                raise DataHallInvariantError(
                    f"{context} source quotient translation is outside centering lattice"
                )


def _expanded_normalization(source_operations, hall_operations, centering_cosets, context):
    hall_by_rotation_translation = {}
    for hall_index, (rotation, translation) in enumerate(hall_operations):
        key = (rotation, translation)
        if key in hall_by_rotation_translation:
            raise DataHallInvariantError(f"{context} Hall operation key is duplicated")
        hall_by_rotation_translation[key] = hall_index
    shift_counts = Counter()
    for source_index, (source_rotation, source_translation) in enumerate(source_operations):
        for coset in centering_cosets:
            expanded = tuple(
                source_translation[index] + coset[index] for index in range(3)
            )
            canonical = _mod1_vector(expanded)
            shift = tuple(expanded[index] - canonical[index] for index in range(3))
            shift_numerator = _numerators_over_12(
                shift, f"{context} expanded normalization shift"
            )
            if (source_rotation, canonical) not in hall_by_rotation_translation:
                raise DataHallInvariantError(
                    f"{context} expanded source operation has no selected Hall representative"
                )
            shift_counts[shift_numerator] += 1
    return _distribution(shift_counts)


def _derive_one(source, provenance, *, candidates=None, details=None):
    spacegroup = _source_spacegroup(source)
    source_symbol, centering = _source_centering(source)
    source_operations = _source_operations(source)
    if candidates is None or details is None:
        candidates, details = _raw_candidates(source, provenance)
    if not candidates:
        raise DataHallDerivationError(f"SG{spacegroup} has no raw Hall candidate")
    centering_cosets = _centering_cosets_for(centering)
    filtered = _filtered_candidates(candidates, details, centering_cosets)
    if len(filtered) != 1:
        raise DataHallDerivationError(
            f"SG{spacegroup} centering filter is not unique: {filtered}"
        )
    data_hall = filtered[0]
    hall_operations = _details_for(details, data_hall)[0]
    source_to_hall, hall_to_source = _mapping_for_selected_hall(
        source_operations, hall_operations, centering_cosets, f"SG{spacegroup}"
    )
    _source_quotient_check(source_operations, centering_cosets, f"SG{spacegroup}")
    expanded_shifts = _expanded_normalization(
        source_operations, hall_operations, centering_cosets, f"SG{spacegroup}"
    )
    source_nonzero = sum(
        mapping.shift_numerator != (0, 0, 0) for mapping in source_to_hall
    )
    hall_shift_counts = Counter(tuple(mapping.shift_numerator)
                                for mapping in hall_to_source)
    hall_coset_counts = Counter(
        tuple(mapping.shift_numerator[index] % TRANSLATION_DENOMINATOR
              for index in range(3))
        for mapping in hall_to_source
    )
    hall_nonzero = sum(
        count for shift, count in hall_shift_counts.items()
        if shift != (0, 0, 0)
    )
    frame = ExactDataHallFrame(
        spacegroup=spacegroup,
        source_symbol=source_symbol,
        centering=centering,
        raw_candidate_halls=tuple(candidates),
        data_hall=data_hall,
        source_operation_count=len(source_operations),
        hall_operation_count=len(hall_operations),
        source_to_hall=source_to_hall,
        hall_to_source=hall_to_source,
    )
    return frame, source_nonzero, hall_nonzero, _distribution(hall_shift_counts), _distribution(hall_coset_counts), expanded_shifts


def _records_and_count(source_database):
    try:
        pir_records = tuple(source_database.pir_records)
        cir_records = tuple(source_database.cir_records)
        universes = tuple(source_database.universes)
    except Exception as error:
        raise DataHallDerivationError("exact source database is not typed") from error
    if len(universes) != 231 or any(universe is None for universe in universes[1:]):
        raise DataHallInvariantError("exact source database does not contain 230 universes")
    if universes[0] is not None:
        raise DataHallInvariantError("exact source universe slot zero is not None")
    return pir_records, cir_records


def _derive_from_databases(source_database, provenance, *, enforce_census: bool):
    pir_records, cir_records = _records_and_count(source_database)
    if enforce_census and (
        len(pir_records) != EXPECTED_PIR_RECORDS
        or len(cir_records) != EXPECTED_CIR_RECORDS
    ):
        raise DataHallInvariantError("exact source record census mismatch")

    frames = []
    raw_counts = Counter()
    filtered_counts = Counter()
    raw_ambiguous = []
    centering_counts = Counter()
    hall_shift_counts = Counter()
    hall_coset_counts = Counter()
    expanded_shift_counts = Counter()
    source_representatives = 0
    selected_hall_operations = 0
    source_to_hall_total = 0
    source_to_hall_nonzero = 0
    hall_to_source_total = 0
    hall_to_source_nonzero = 0

    for spacegroup in range(1, 231):
        try:
            source = source_database.source_universe(spacegroup)
        except Exception as error:
            raise DataHallDerivationError(
                f"unable to read exact source universe SG{spacegroup}"
            ) from error
        if _source_spacegroup(source) != spacegroup:
            raise DataHallInvariantError(f"source universe slot mismatch at SG{spacegroup}")
        candidates, details = _raw_candidates(source, provenance)
        raw_counts[len(candidates)] += 1
        if len(candidates) > 1:
            raw_ambiguous.append(spacegroup)
        if not candidates:
            raise DataHallDerivationError(f"SG{spacegroup} has no raw Hall candidate")
        _, centering = _source_centering(source)
        filtered = _filtered_candidates(
            candidates, details, _centering_cosets_for(centering)
        )
        filtered_counts[len(filtered)] += 1
        if len(filtered) != 1:
            raise DataHallDerivationError(
                f"SG{spacegroup} centering filter is not unique: {filtered}"
            )
        frame, source_nonzero, hall_nonzero, hall_shifts, hall_cosets, expanded_shifts = _derive_one(
            source, provenance, candidates=candidates, details=details
        )
        frames.append(frame)
        source_representatives += frame.source_operation_count
        selected_hall_operations += frame.hall_operation_count
        source_to_hall_total += len(frame.source_to_hall)
        source_to_hall_nonzero += source_nonzero
        hall_to_source_total += len(frame.hall_to_source)
        hall_to_source_nonzero += hall_nonzero
        for shift, count in hall_shifts:
            hall_shift_counts[shift] += count
        for shift, count in hall_cosets:
            hall_coset_counts[shift] += count
        for shift, count in expanded_shifts:
            expanded_shift_counts[shift] += count
        centering_counts[centering] += 1

    if enforce_census:
        if raw_counts[1] != EXPECTED_RAW_UNIQUE:
            raise DataHallInvariantError("raw unique Hall census mismatch")
        if tuple(raw_ambiguous) != EXPECTED_RAW_AMBIGUOUS_SPACEGROUPS:
            raise DataHallInvariantError("raw ambiguous Hall census mismatch")
        if raw_counts[0] != EXPECTED_RAW_MISSING:
            raise DataHallInvariantError("raw missing Hall census mismatch")
        if filtered_counts[1] != EXPECTED_FILTERED_UNIQUE:
            raise DataHallInvariantError("filtered unique Hall census mismatch")
        if filtered_counts[2] != EXPECTED_FILTERED_AMBIGUOUS:
            raise DataHallInvariantError("filtered ambiguous Hall census mismatch")
        if filtered_counts[0] != EXPECTED_FILTERED_MISSING:
            raise DataHallInvariantError("filtered missing Hall census mismatch")
        if source_representatives != EXPECTED_SOURCE_REPRESENTATIVES:
            raise DataHallInvariantError("source representative census mismatch")
        if selected_hall_operations != EXPECTED_SELECTED_HALL_OPERATIONS:
            raise DataHallInvariantError("selected Hall operation census mismatch")
        if source_to_hall_total != EXPECTED_SOURCE_TO_HALL:
            raise DataHallInvariantError("source-to-Hall census mismatch")
        if hall_to_source_total != EXPECTED_HALL_TO_SOURCE:
            raise DataHallInvariantError("Hall-to-source census mismatch")
        if source_to_hall_nonzero != EXPECTED_SOURCE_TO_HALL_NONZERO:
            raise DataHallInvariantError("source-to-Hall nonzero census mismatch")
        if hall_to_source_nonzero != EXPECTED_HALL_TO_SOURCE_NONZERO:
            raise DataHallInvariantError("Hall-to-source nonzero census mismatch")
        expanded_nonzero = sum(
            count for shift, count in expanded_shift_counts.items()
            if shift != (0, 0, 0)
        )
        if expanded_nonzero != EXPECTED_EXPANDED_NORMALIZATION_NONZERO:
            raise DataHallInvariantError("expanded normalization census mismatch")
        observed_centering = tuple(
            (name, centering_counts.get(name, 0))
            for name, _ in EXPECTED_CENTERING_COUNTS
        )
        if observed_centering != EXPECTED_CENTERING_COUNTS:
            raise DataHallInvariantError("centering census mismatch")
    else:
        expanded_nonzero = sum(
            count for shift, count in expanded_shift_counts.items()
            if shift != (0, 0, 0)
        )

    census = DerivationCensus(
        pir_records=len(pir_records),
        cir_records=len(cir_records),
        source_representatives=source_representatives,
        raw_unique=raw_counts[1],
        raw_ambiguous=len(raw_ambiguous),
        raw_missing=raw_counts[0],
        raw_ambiguous_spacegroups=tuple(raw_ambiguous),
        filtered_unique=filtered_counts[1],
        filtered_ambiguous=filtered_counts[2],
        filtered_missing=filtered_counts[0],
        selected_hall_operations=selected_hall_operations,
        source_to_hall=source_to_hall_total,
        source_to_hall_nonzero=source_to_hall_nonzero,
        hall_to_source=hall_to_source_total,
        hall_to_source_nonzero=hall_to_source_nonzero,
        hall_to_source_shifts=_distribution(hall_shift_counts),
        hall_to_source_cosets=_distribution(hall_coset_counts),
        expanded_normalization_nonzero=expanded_nonzero,
        expanded_normalization_shifts=_distribution(expanded_shift_counts),
        centering_counts=tuple(
            (name, centering_counts.get(name, 0))
            for name, _ in EXPECTED_CENTERING_COUNTS
        ),
    )
    return ExactDataHallDatabase(tuple(frames), census)


def derive_data_hall_frames(source_db=None, spg_db=None) -> ExactDataHallDatabase:
    """Derive all exact direct-source Hall frames in memory.

    With no arguments, only the public committed exact ISO-IR and spglib
    provenance loaders are called and the complete pinned census is enforced.
    Supplying both databases is a test seam: all operation/mapping invariants
    still apply, while a synthetic complete fixture may omit the production
    census.  Supplying only one database is an error.
    """

    if (source_db is None) != (spg_db is None):
        raise TypeError("source_db and spg_db must be supplied together")
    enforce_census = source_db is None
    if source_db is None:
        try:
            source_db = iso_irrep_exact.load_exact_iso_irrep_sources()
            spg_db = spglib_magnetic_provenance.load_committed_provenance()
        except IsoIrrepDataHallError:
            raise
        except (ValueError, OSError, TypeError) as error:
            raise DataHallDerivationError("authoritative input loader failed") from error
    return _derive_from_databases(source_db, spg_db, enforce_census=enforce_census)


# Internal spelling retained for focused tests and downstream callers that
# want to make the production-vs-synthetic census choice explicit.
derive_data_hall_authority = derive_data_hall_frames


__all__ = [
    "CENTERING_COSETS", "CENTERING_MATRIX_DATA", "DataHallDerivationError",
    "DataHallInvariantError", "DataHallLookupError", "DataHallSchemaError",
    "DerivationCensus", "ExactDataHallDatabase", "ExactDataHallFrame",
    "HallToSource", "HallToSourceMapping", "IsoIrrepDataHallError",
    "SourceToHall", "SourceToHallMapping", "TRANSLATION_DENOMINATOR",
    "derive_data_hall_authority", "derive_data_hall_frames",
]
