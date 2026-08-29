#!/usr/bin/env python3
"""Load the committed exact ISO--IR data--Hall provenance sidecar.

The loader is intentionally an independent, fixed-file provenance boundary.
It reads the two committed JSON files exactly once on a cold load, verifies
their fixed bytes and SHA-256 values before parsing, then validates the full
canonical JSON and immutable typed graph.  The stored frame is the direct
source frame ``P=I, p=0``; a source-to-Hall shift means
``source = hall + shift`` and a Hall-to-source shift means
``hall = source + shift``.  Hashes are file-integrity checks only, never a
search or selection mechanism.

Dataclass constructors are value constructors: a caller may construct a
semantically valid value for testing, but that value is not the pinned
authority.  Only :func:`load_committed_data_hall_provenance` establishes the
fixed committed authority.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import threading
from typing import Optional


_NATIVE_DATACLASS_SLOTS = False
try:
    _NATIVE_DATACLASS_SLOTS = "slots" in __import__(
        "inspect"
    ).signature(dataclass).parameters
except (TypeError, ValueError):
    pass
_DATACLASS_OPTIONS = {"frozen": True}
if _NATIVE_DATACLASS_SLOTS:
    _DATACLASS_OPTIONS["slots"] = True


_MODULE_DIR = Path(__file__).resolve().parent
_DATA_DIR = _MODULE_DIR / "data"
_ARTIFACT_NAME = "iso_irrep_data_hall_v1.json"
_MANIFEST_NAME = "iso_irrep_data_hall_v1.manifest.json"
_ARTIFACT_PATH = _DATA_DIR / _ARTIFACT_NAME
_MANIFEST_PATH = _DATA_DIR / _MANIFEST_NAME

SCHEMA = "cryspglib-iso-irrep-data-hall-v1"
MANIFEST_SCHEMA = "cryspglib-iso-irrep-data-hall-manifest-v1"
TRANSLATION_DENOMINATOR = 12
FRAME_SEMANTICS = "direct-source-frame-P=I-p=0"
_MAPPING_SEMANTICS = {
    "source_to_hall": "source=hall+shift",
    "hall_to_source": "hall=source+shift",
}

ARTIFACT_BYTE_LENGTH = 697_730
ARTIFACT_SHA256 = (
    "35bcb00958021eb6fc5a330f8dbf85a80be78ccec324f441e6138cdba4b617e0"
)
MANIFEST_BYTE_LENGTH = 869
MANIFEST_SHA256 = (
    "bc6aa7a94d698f2193e7cb623b16dded2dd8e0307d502cf28b68c554f364d7e2"
)

_CENTERING_ORDER = ("P", "A", "B", "C", "F", "I", "R")
_CENTERING_RESIDUES = {
    "P": ((0, 0, 0),),
    "A": ((0, 0, 0), (0, 6, 6)),
    "B": ((0, 0, 0), (6, 0, 6)),
    "C": ((0, 0, 0), (6, 6, 0)),
    "F": ((0, 0, 0), (0, 6, 6), (6, 0, 6), (6, 6, 0)),
    "I": ((0, 0, 0), (6, 6, 6)),
    "R": ((0, 0, 0), (4, 8, 8), (8, 4, 4)),
}
_IDENTITY_BASIS = (1, 0, 0, 0, 1, 0, 0, 0, 1)
_ZERO_ORIGIN = (0, 0, 0)

_INPUT_SPECS = (
    (
        "pir_zip", "isotropy_subgroup/PIR_data.zip", 1_235_319,
        "e909a4f0121688b0590ccaec10b0276171bc24619cf7eb562ba441268c01e121",
    ),
    (
        "cir_zip", "isotropy_subgroup/CIR_data.zip", 1_555_153,
        "f4edcb2852b83a86d1b58f29fb862d9124a227cfc90f9e1ae17d2c97585264e6",
    ),
    (
        "spglib_artifact", "scripts/data/spglib_magnetic_provenance_v1.json",
        1_537_875,
        "933a52a6696e7f6a1a2e426825ad92c377c6e96330e18c5c045d659798d740b9",
    ),
    (
        "spglib_manifest",
        "scripts/data/spglib_magnetic_provenance_v1.manifest.json", 570,
        "6a9e1b64c190c30a556d63e51e5b896b967d33e8821714beb745ae699fab84bf",
    ),
)

_EXPECTED_CENSUS = {
    "pir_records": 10_294,
    "cir_records": 11_202,
    "source_representatives": 2_609,
    "raw_unique": 220,
    "raw_ambiguous": 10,
    "raw_missing": 0,
    "raw_ambiguous_spacegroups": (
        5, 8, 9, 12, 15, 21, 38, 39, 65, 67,
    ),
    "filtered_unique": 230,
    "filtered_ambiguous": 0,
    "filtered_missing": 0,
    "selected_hall_operations": 4_425,
    "source_to_hall": 2_609,
    "source_to_hall_nonzero": 0,
    "hall_to_source": 4_425,
    "hall_to_source_nonzero": 1_816,
    "hall_to_source_shifts": (
        ((-6, -6, -6), 1),
        ((-6, -6, 0), 48),
        ((-6, -6, 6), 4),
        ((-6, 0, -6), 46),
        ((-6, 0, 6), 10),
        ((-6, 6, -6), 3),
        ((-6, 6, 0), 20),
        ((-6, 6, 6), 44),
        ((0, -6, -6), 46),
        ((0, -6, 6), 10),
        ((0, 0, 0), 2_609),
        ((0, 6, -6), 14),
        ((0, 6, 6), 322),
        ((4, 8, -4), 9),
        ((4, 8, 8), 42),
        ((6, -6, -6), 2),
        ((6, -6, 0), 12),
        ((6, -6, 6), 67),
        ((6, 0, -6), 10),
        ((6, 0, 6), 310),
        ((6, 6, -6), 64),
        ((6, 6, 0), 378),
        ((6, 6, 6), 303),
        ((8, 4, 4), 51),
    ),
    "hall_to_source_cosets": (
        ((0, 0, 0), 2_609),
        ((0, 6, 6), 392),
        ((4, 8, 8), 51),
        ((6, 0, 6), 376),
        ((6, 6, 0), 458),
        ((6, 6, 6), 488),
        ((8, 4, 4), 51),
    ),
    "expanded_normalization_nonzero": 410,
    "expanded_normalization_shifts": (
        ((0, 0, 0), 4_015),
        ((0, 0, 12), 97),
        ((0, 12, 0), 89),
        ((0, 12, 12), 48),
        ((12, 0, 0), 74),
        ((12, 0, 12), 49),
        ((12, 12, 0), 52),
        ((12, 12, 12), 1),
    ),
    "centering_counts": (
        ("P", 149), ("A", 4), ("B", 0), ("C", 16),
        ("F", 16), ("I", 38), ("R", 7),
    ),
}


class DataHallProvenanceError(ValueError):
    """Base class for fixed data--Hall provenance failures."""


class DataHallIntegrityError(DataHallProvenanceError):
    """A committed file or byte-level commitment is not trustworthy."""


class DataHallSchemaError(DataHallProvenanceError):
    """A JSON value or typed value has the wrong exact schema."""


class DataHallInvariantError(DataHallProvenanceError):
    """A schema-valid value violates the frozen graph or census laws."""


class DataHallLookupError(DataHallProvenanceError):
    """A public space-group lookup target is invalid."""


def _schema(context: str, message: str):
    raise DataHallSchemaError(f"{context}: {message}")


def _invariant(context: str, message: str):
    raise DataHallInvariantError(f"{context}: {message}")


def _exact_int(value, context: str, *, minimum: Optional[int] = None,
               maximum: Optional[int] = None) -> int:
    if type(value) is not int:
        _schema(context, "must be an exact integer")
    if minimum is not None and value < minimum:
        _invariant(context, f"must be at least {minimum}")
    if maximum is not None and value > maximum:
        _invariant(context, f"must be at most {maximum}")
    return value


def _exact_str(value, context: str, *, nonempty: bool = False) -> str:
    if type(value) is not str:
        _schema(context, "must be an exact string")
    if any(ord(char) < 0x20 or ord(char) > 0x7E for char in value):
        _schema(context, "must contain printable ASCII only")
    if nonempty and not value:
        _invariant(context, "must not be empty")
    return value


def _exact_tuple(value, context: str, length: Optional[int] = None):
    if type(value) is not tuple:
        _schema(context, "must be an exact tuple")
    if length is not None and len(value) != length:
        _schema(context, f"must contain exactly {length} items")
    return value


def _exact_shift(value, context: str):
    shift = _exact_tuple(value, context, 3)
    return tuple(
        _exact_int(component, f"{context}[{index}]")
        for index, component in enumerate(shift)
    )


@dataclass(**_DATACLASS_OPTIONS)
class SourceToHallMapping:
    if not _NATIVE_DATACLASS_SLOTS:
        __slots__ = (
            "source_operation_index", "hall_operation_index", "shift_numerator"
        )
    source_operation_index: int
    hall_operation_index: int
    shift_numerator: tuple

    def __post_init__(self):
        _exact_int(self.source_operation_index, "source_operation_index", minimum=0)
        _exact_int(self.hall_operation_index, "hall_operation_index", minimum=0)
        _exact_shift(self.shift_numerator, "shift_numerator")


@dataclass(**_DATACLASS_OPTIONS)
class HallToSourceMapping:
    if not _NATIVE_DATACLASS_SLOTS:
        __slots__ = (
            "hall_operation_index", "source_operation_index", "shift_numerator"
        )
    hall_operation_index: int
    source_operation_index: int
    shift_numerator: tuple

    def __post_init__(self):
        _exact_int(self.hall_operation_index, "hall_operation_index", minimum=0)
        _exact_int(self.source_operation_index, "source_operation_index", minimum=0)
        _exact_shift(self.shift_numerator, "shift_numerator")


@dataclass(**_DATACLASS_OPTIONS)
class DataHallFrame:
    if not _NATIVE_DATACLASS_SLOTS:
        __slots__ = (
            "spacegroup", "source_symbol", "centering",
            "pir_anchor_irnumber", "cir_anchor_irnumber",
            "raw_candidate_halls", "data_hall", "basis", "origin_numerator",
            "source_operation_count", "hall_operation_count", "source_to_hall",
            "hall_to_source",
        )
    spacegroup: int
    source_symbol: str
    centering: str
    pir_anchor_irnumber: int
    cir_anchor_irnumber: int
    raw_candidate_halls: tuple
    data_hall: int
    basis: tuple
    origin_numerator: tuple
    source_operation_count: int
    hall_operation_count: int
    source_to_hall: tuple
    hall_to_source: tuple

    def __post_init__(self):
        _validate_frame(self)


@dataclass(**_DATACLASS_OPTIONS)
class DataHallCensus:
    if not _NATIVE_DATACLASS_SLOTS:
        __slots__ = (
            "pir_records", "cir_records", "source_representatives",
            "raw_unique", "raw_ambiguous", "raw_missing",
            "raw_ambiguous_spacegroups", "filtered_unique",
            "filtered_ambiguous", "filtered_missing", "selected_hall_operations",
            "source_to_hall", "source_to_hall_nonzero", "hall_to_source",
            "hall_to_source_nonzero", "hall_to_source_shifts",
            "hall_to_source_cosets", "expanded_normalization_nonzero",
            "expanded_normalization_shifts", "centering_counts",
        )
    pir_records: int
    cir_records: int
    source_representatives: int
    raw_unique: int
    raw_ambiguous: int
    raw_missing: int
    raw_ambiguous_spacegroups: tuple
    filtered_unique: int
    filtered_ambiguous: int
    filtered_missing: int
    selected_hall_operations: int
    source_to_hall: int
    source_to_hall_nonzero: int
    hall_to_source: int
    hall_to_source_nonzero: int
    hall_to_source_shifts: tuple
    hall_to_source_cosets: tuple
    expanded_normalization_nonzero: int
    expanded_normalization_shifts: tuple
    centering_counts: tuple

    def __post_init__(self):
        _validate_census_fields(self)


@dataclass(**_DATACLASS_OPTIONS)
class DataHallProvenanceDatabase:
    if not _NATIVE_DATACLASS_SLOTS:
        __slots__ = ("frames", "census")
    frames: tuple
    census: DataHallCensus

    def __post_init__(self):
        _validate_database_graph(self.frames, self.census)


_CENSUS_SCALAR_FIELDS = (
    "pir_records", "cir_records", "source_representatives", "raw_unique",
    "raw_ambiguous", "raw_missing", "filtered_unique",
    "filtered_ambiguous", "filtered_missing", "selected_hall_operations",
    "source_to_hall", "source_to_hall_nonzero", "hall_to_source",
    "hall_to_source_nonzero", "expanded_normalization_nonzero",
)


def _validate_mapping_leaf(mapping, expected_type, context: str):
    if type(mapping) is not expected_type:
        _schema(context, "must be the exact mapping type")
    if expected_type is SourceToHallMapping:
        _exact_int(mapping.source_operation_index, f"{context}.source_operation_index", minimum=0)
        _exact_int(mapping.hall_operation_index, f"{context}.hall_operation_index", minimum=0)
    else:
        _exact_int(mapping.hall_operation_index, f"{context}.hall_operation_index", minimum=0)
        _exact_int(mapping.source_operation_index, f"{context}.source_operation_index", minimum=0)
    _exact_shift(mapping.shift_numerator, f"{context}.shift_numerator")


def _validate_frame(frame: DataHallFrame):
    context = "frame"
    spacegroup = _exact_int(frame.spacegroup, f"{context}.spacegroup", minimum=1, maximum=230)
    symbol = _exact_str(frame.source_symbol, f"{context}.source_symbol", nonempty=True)
    centering = _exact_str(frame.centering, f"{context}.centering", nonempty=True)
    if centering not in _CENTERING_RESIDUES:
        _invariant(f"{context}.centering", "is unknown")
    if symbol[0] != centering:
        _invariant(context, "source symbol and centering disagree")
    _exact_int(frame.pir_anchor_irnumber, f"{context}.pir_anchor_irnumber", minimum=1)
    _exact_int(frame.cir_anchor_irnumber, f"{context}.cir_anchor_irnumber", minimum=1)

    raw_halls = _exact_tuple(frame.raw_candidate_halls, f"{context}.raw_candidate_halls")
    if not raw_halls:
        _invariant(f"{context}.raw_candidate_halls", "must not be empty")
    for index, hall in enumerate(raw_halls):
        _exact_int(hall, f"{context}.raw_candidate_halls[{index}]", minimum=1, maximum=530)
    if raw_halls != tuple(sorted(set(raw_halls))):
        _invariant(f"{context}.raw_candidate_halls", "must be sorted and unique")
    data_hall = _exact_int(frame.data_hall, f"{context}.data_hall", minimum=1, maximum=530)
    if data_hall not in raw_halls:
        _invariant(context, "selected Hall is not a raw candidate")

    basis = _exact_tuple(frame.basis, f"{context}.basis", 9)
    for index, value in enumerate(basis):
        _exact_int(value, f"{context}.basis[{index}]")
    if basis != _IDENTITY_BASIS:
        _invariant(f"{context}.basis", "must be identity")
    origin = _exact_tuple(frame.origin_numerator, f"{context}.origin_numerator", 3)
    for index, value in enumerate(origin):
        _exact_int(value, f"{context}.origin_numerator[{index}]")
    if origin != _ZERO_ORIGIN:
        _invariant(f"{context}.origin_numerator", "must be zero")

    source_count = _exact_int(
        frame.source_operation_count, f"{context}.source_operation_count",
        minimum=1, maximum=48,
    )
    hall_count = _exact_int(
        frame.hall_operation_count, f"{context}.hall_operation_count",
        minimum=1, maximum=192,
    )
    residues = _CENTERING_RESIDUES[centering]
    if hall_count != source_count * len(residues):
        _invariant(context, "operation counts disagree with centering")

    source_maps = _exact_tuple(frame.source_to_hall, f"{context}.source_to_hall")
    hall_maps = _exact_tuple(frame.hall_to_source, f"{context}.hall_to_source")
    if len(source_maps) != source_count:
        _invariant(f"{context}.source_to_hall", "cardinality disagrees")
    if len(hall_maps) != hall_count:
        _invariant(f"{context}.hall_to_source", "cardinality disagrees")

    source_hall_indices = []
    for index, mapping in enumerate(source_maps):
        _validate_mapping_leaf(
            mapping, SourceToHallMapping, f"{context}.source_to_hall[{index}]"
        )
        if mapping.source_operation_index != index:
            _invariant(f"{context}.source_to_hall[{index}]", "source order changed")
        if mapping.hall_operation_index >= hall_count:
            _invariant(f"{context}.source_to_hall[{index}]", "Hall index out of range")
        if any(value % TRANSLATION_DENOMINATOR for value in mapping.shift_numerator):
            _invariant(f"{context}.source_to_hall[{index}]", "shift is not integral")
        source_hall_indices.append(mapping.hall_operation_index)
    if len(set(source_hall_indices)) != source_count:
        _invariant(f"{context}.source_to_hall", "Hall indices are not unique")

    residues_by_source = [[] for _ in range(source_count)]
    for index, mapping in enumerate(hall_maps):
        _validate_mapping_leaf(
            mapping, HallToSourceMapping, f"{context}.hall_to_source[{index}]"
        )
        if mapping.hall_operation_index != index:
            _invariant(f"{context}.hall_to_source[{index}]", "Hall order changed")
        if mapping.source_operation_index >= source_count:
            _invariant(f"{context}.hall_to_source[{index}]", "source index out of range")
        residue = tuple(
            value % TRANSLATION_DENOMINATOR
            for value in mapping.shift_numerator
        )
        residues_by_source[mapping.source_operation_index].append(residue)

    expected_residues = set(residues)
    for source_index, source_residues in enumerate(residues_by_source):
        if len(source_residues) != len(expected_residues):
            _invariant(
                f"{context}.hall_to_source[{source_index}]",
                "centering residue cardinality disagrees",
            )
        if set(source_residues) != expected_residues:
            _invariant(
                f"{context}.hall_to_source[{source_index}]",
                "centering residues are incomplete",
            )

    for source_index, mapping in enumerate(source_maps):
        inverse = hall_maps[mapping.hall_operation_index]
        if inverse.source_operation_index != source_index:
            _invariant(context, "source/Hall mapping indices disagree")
        if inverse.shift_numerator != tuple(-value for value in mapping.shift_numerator):
            _invariant(context, "source/Hall mapping shifts disagree")


def _validate_distribution(value, context: str):
    rows = _exact_tuple(value, context)
    parsed = []
    for index, row in enumerate(rows):
        row = _exact_tuple(row, f"{context}[{index}]", 2)
        key = _exact_tuple(row[0], f"{context}[{index}].key", 3)
        key = tuple(
            _exact_int(component, f"{context}[{index}].key[{component_index}]")
            for component_index, component in enumerate(key)
        )
        count = _exact_int(row[1], f"{context}[{index}].count", minimum=1)
        parsed.append((key, count))
    if tuple(parsed) != tuple(sorted(parsed)):
        _invariant(context, "distribution rows are not canonically ordered")
    if len({key for key, _ in parsed}) != len(parsed):
        _invariant(context, "distribution keys are duplicated")


def _validate_centering_counts(value, context: str):
    rows = _exact_tuple(value, context, len(_CENTERING_ORDER))
    parsed = []
    for index, row in enumerate(rows):
        row = _exact_tuple(row, f"{context}[{index}]", 2)
        name = _exact_str(row[0], f"{context}[{index}].name")
        count = _exact_int(row[1], f"{context}[{index}].count", minimum=0)
        if name != _CENTERING_ORDER[index]:
            _invariant(context, "centering order changed")
        parsed.append((name, count))
    return tuple(parsed)


def _validate_census_fields(census: DataHallCensus):
    for field in _CENSUS_SCALAR_FIELDS:
        _exact_int(getattr(census, field), f"census.{field}", minimum=0)
    ambiguous = _exact_tuple(
        census.raw_ambiguous_spacegroups,
        "census.raw_ambiguous_spacegroups",
    )
    for index, spacegroup in enumerate(ambiguous):
        _exact_int(
            spacegroup, f"census.raw_ambiguous_spacegroups[{index}]",
            minimum=1, maximum=230,
        )
    if ambiguous != tuple(sorted(set(ambiguous))):
        _invariant("census.raw_ambiguous_spacegroups", "must be sorted and unique")
    for name in (
        "hall_to_source_shifts", "hall_to_source_cosets",
        "expanded_normalization_shifts",
    ):
        _validate_distribution(getattr(census, name), f"census.{name}")
    _validate_centering_counts(census.centering_counts, "census.centering_counts")


def _compute_census(frames, pir_records: int, cir_records: int):
    raw_unique = 0
    raw_ambiguous = 0
    raw_missing = 0
    ambiguous_spacegroups = []
    centerings = Counter()
    source_representatives = 0
    selected_hall_operations = 0
    source_to_hall = 0
    source_to_hall_nonzero = 0
    hall_to_source = 0
    hall_to_source_nonzero = 0
    hall_shifts = Counter()
    hall_cosets = Counter()
    expanded_shifts = Counter()

    for frame in frames:
        _validate_frame(frame)
        candidate_count = len(frame.raw_candidate_halls)
        if candidate_count == 0:
            raw_missing += 1
        elif candidate_count == 1:
            raw_unique += 1
        else:
            raw_ambiguous += 1
            ambiguous_spacegroups.append(frame.spacegroup)
        centerings[frame.centering] += 1
        source_representatives += frame.source_operation_count
        selected_hall_operations += frame.hall_operation_count
        source_to_hall += len(frame.source_to_hall)
        source_to_hall_nonzero += sum(
            shift.shift_numerator != (0, 0, 0)
            for shift in frame.source_to_hall
        )
        hall_to_source += len(frame.hall_to_source)
        for mapping in frame.hall_to_source:
            shift = mapping.shift_numerator
            hall_shifts[shift] += 1
            if shift != (0, 0, 0):
                hall_to_source_nonzero += 1
            residue = tuple(
                value % TRANSLATION_DENOMINATOR for value in shift
            )
            hall_cosets[residue] += 1
            expanded_shifts[tuple(
                residue[index] - shift[index] for index in range(3)
            )] += 1

    def ordered_distribution(counter):
        return tuple(sorted(counter.items()))

    return DataHallCensus(
        pir_records=pir_records,
        cir_records=cir_records,
        source_representatives=source_representatives,
        raw_unique=raw_unique,
        raw_ambiguous=raw_ambiguous,
        raw_missing=raw_missing,
        raw_ambiguous_spacegroups=tuple(ambiguous_spacegroups),
        filtered_unique=len(frames),
        filtered_ambiguous=0,
        filtered_missing=0,
        selected_hall_operations=selected_hall_operations,
        source_to_hall=source_to_hall,
        source_to_hall_nonzero=source_to_hall_nonzero,
        hall_to_source=hall_to_source,
        hall_to_source_nonzero=hall_to_source_nonzero,
        hall_to_source_shifts=ordered_distribution(hall_shifts),
        hall_to_source_cosets=ordered_distribution(hall_cosets),
        expanded_normalization_nonzero=sum(
            count for shift, count in expanded_shifts.items()
            if shift != (0, 0, 0)
        ),
        expanded_normalization_shifts=ordered_distribution(expanded_shifts),
        centering_counts=tuple(
            (name, centerings.get(name, 0)) for name in _CENTERING_ORDER
        ),
    )


def _validate_database_graph(frames, census):
    if type(frames) is not tuple:
        _schema("database.frames", "must be an exact tuple")
    if len(frames) != 230:
        _invariant("database.frames", "must contain exactly 230 frames")
    for index, frame in enumerate(frames, 1):
        if type(frame) is not DataHallFrame:
            _schema(f"database.frames[{index - 1}]", "has wrong exact type")
        _validate_frame(frame)
        if frame.spacegroup != index:
            _invariant("database.frames", "spacegroups must be ordered 1..230")
    if type(census) is not DataHallCensus:
        _schema("database.census", "has wrong exact type")
    _validate_census_fields(census)
    computed = _compute_census(frames, census.pir_records, census.cir_records)
    if computed != census:
        _invariant("database.census", "does not match frame mappings")
    for field, expected in _EXPECTED_CENSUS.items():
        if getattr(census, field) != expected:
            _invariant(f"database.census.{field}", "does not match frozen census")


def _require_keys(value, keys, context: str):
    if type(value) is not dict:
        _schema(context, "must be an exact object")
    if set(value) != set(keys):
        _schema(context, "has missing or extra keys")


def _parse_canonical_json(data: bytes, context: str):
    if type(data) is not bytes:
        _schema(context, "must be exact bytes")
    if not data.endswith(b"\n"):
        _schema(context, "must have exactly one final LF")
    if any(byte >= 0x80 for byte in data):
        _schema(context, "must contain ASCII UTF-8 only")
    try:
        text = data[:-1].decode("utf-8", errors="strict")
        value = json.loads(
            text,
            object_pairs_hook=_pairs_without_duplicates,
            parse_float=lambda value: (_ for _ in ()).throw(
                DataHallSchemaError(f"{context}: floats are forbidden")
            ),
            parse_constant=lambda value: (_ for _ in ()).throw(
                DataHallSchemaError(f"{context}: {value} is forbidden")
            ),
        )
    except DataHallProvenanceError:
        raise
    except RecursionError as error:
        raise DataHallSchemaError(f"{context}: JSON nesting is too deep") from error
    except (UnicodeError, ValueError, TypeError) as error:
        raise DataHallSchemaError(f"{context}: invalid JSON") from error
    try:
        _validate_json_tree(value, context)
        if _canonical_json(value) != data:
            _schema(context, "is not canonical JSON")
    except DataHallProvenanceError:
        raise
    except RecursionError as error:
        raise DataHallSchemaError(f"{context}: JSON nesting is too deep") from error
    return value


def _pairs_without_duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise DataHallSchemaError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _validate_json_tree(value, context: str = "$"):
    if type(value) is dict:
        for key, child in value.items():
            if type(key) is not str:
                _schema(context, "object keys must be exact strings")
            if any(ord(char) < 0x20 or ord(char) > 0x7E for char in key):
                _schema(context, "object keys must be printable ASCII")
            _validate_json_tree(child, f"{context}.{key}")
        return
    if type(value) is list:
        for index, child in enumerate(value):
            _validate_json_tree(child, f"{context}[{index}]")
        return
    if type(value) is str:
        if any(ord(char) < 0x20 or ord(char) > 0x7E for char in value):
            _schema(context, "strings must be printable ASCII")
        return
    if type(value) is int:
        return
    _schema(context, "null, bool, float, tuple, and non-built-in values are forbidden")


def _canonical_json(value) -> bytes:
    try:
        _validate_json_tree(value)
        encoded = json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except DataHallProvenanceError:
        raise
    except RecursionError as error:
        raise DataHallSchemaError("JSON value is too deeply nested") from error
    except (TypeError, ValueError, UnicodeError) as error:
        raise DataHallSchemaError("unable to encode canonical JSON") from error
    return encoded + b"\n"


def _parse_int_list(value, context: str, *, length: Optional[int] = None):
    if type(value) is not list:
        _schema(context, "must be an exact array")
    if length is not None and len(value) != length:
        _schema(context, f"must contain exactly {length} items")
    return tuple(
        _exact_int(item, f"{context}[{index}]")
        for index, item in enumerate(value)
    )


def _parse_mapping(value, expected_keys, expected_type, context: str):
    _require_keys(value, expected_keys, context)
    first = _exact_int(value[expected_keys[0]], f"{context}.{expected_keys[0]}", minimum=0)
    second = _exact_int(value[expected_keys[1]], f"{context}.{expected_keys[1]}", minimum=0)
    shift = _parse_int_list(
        value["lattice_shift_numerator"],
        f"{context}.lattice_shift_numerator",
        length=3,
    )
    if expected_type is SourceToHallMapping:
        return SourceToHallMapping(first, second, shift)
    return HallToSourceMapping(first, second, shift)


def _parse_frame(value, context: str):
    keys = (
        "spacegroup", "source_symbol", "centering", "pir_anchor_irnumber",
        "cir_anchor_irnumber", "raw_candidate_halls", "data_hall", "basis",
        "origin_numerator", "source_operation_count", "hall_operation_count",
        "source_to_hall", "hall_to_source",
    )
    _require_keys(value, keys, context)
    source_maps_raw = value["source_to_hall"]
    hall_maps_raw = value["hall_to_source"]
    if type(source_maps_raw) is not list:
        _schema(f"{context}.source_to_hall", "must be an exact array")
    if type(hall_maps_raw) is not list:
        _schema(f"{context}.hall_to_source", "must be an exact array")
    source_maps = tuple(
        _parse_mapping(
            mapping,
            ("source_operation_index", "hall_operation_index", "lattice_shift_numerator"),
            SourceToHallMapping,
            f"{context}.source_to_hall[{index}]",
        )
        for index, mapping in enumerate(source_maps_raw)
    )
    hall_maps = tuple(
        _parse_mapping(
            mapping,
            ("hall_operation_index", "source_operation_index", "lattice_shift_numerator"),
            HallToSourceMapping,
            f"{context}.hall_to_source[{index}]",
        )
        for index, mapping in enumerate(hall_maps_raw)
    )
    return DataHallFrame(
        spacegroup=_exact_int(value["spacegroup"], f"{context}.spacegroup"),
        source_symbol=_exact_str(value["source_symbol"], f"{context}.source_symbol"),
        centering=_exact_str(value["centering"], f"{context}.centering"),
        pir_anchor_irnumber=_exact_int(
            value["pir_anchor_irnumber"], f"{context}.pir_anchor_irnumber"
        ),
        cir_anchor_irnumber=_exact_int(
            value["cir_anchor_irnumber"], f"{context}.cir_anchor_irnumber"
        ),
        raw_candidate_halls=_parse_int_list(
            value["raw_candidate_halls"], f"{context}.raw_candidate_halls"
        ),
        data_hall=_exact_int(value["data_hall"], f"{context}.data_hall"),
        basis=_parse_int_list(value["basis"], f"{context}.basis", length=9),
        origin_numerator=_parse_int_list(
            value["origin_numerator"], f"{context}.origin_numerator", length=3
        ),
        source_operation_count=_exact_int(
            value["source_operation_count"], f"{context}.source_operation_count"
        ),
        hall_operation_count=_exact_int(
            value["hall_operation_count"], f"{context}.hall_operation_count"
        ),
        source_to_hall=source_maps,
        hall_to_source=hall_maps,
    )


def _parse_distribution(value, context: str):
    if type(value) is not list:
        _schema(context, "must be an exact array")
    rows = []
    for index, row in enumerate(value):
        if type(row) is not list or len(row) != 2:
            _schema(f"{context}[{index}]", "must be a two-item array")
        key = _parse_int_list(row[0], f"{context}[{index}][0]", length=3)
        count = _exact_int(row[1], f"{context}[{index}][1]", minimum=1)
        rows.append((key, count))
    return tuple(rows)


def _parse_centering_counts(value, context: str):
    if type(value) is not list:
        _schema(context, "must be an exact array")
    rows = []
    for index, row in enumerate(value):
        if type(row) is not list or len(row) != 2:
            _schema(f"{context}[{index}]", "must be a two-item array")
        rows.append((
            _exact_str(row[0], f"{context}[{index}][0]"),
            _exact_int(row[1], f"{context}[{index}][1]", minimum=0),
        ))
    return tuple(rows)


def _parse_census(value, context: str):
    keys = (
        "pir_records", "cir_records", "source_representatives", "raw_unique",
        "raw_ambiguous", "raw_missing", "raw_ambiguous_spacegroups",
        "filtered_unique", "filtered_ambiguous", "filtered_missing",
        "selected_hall_operations", "source_to_hall", "source_to_hall_nonzero",
        "hall_to_source", "hall_to_source_nonzero", "hall_to_source_shifts",
        "hall_to_source_cosets", "expanded_normalization_nonzero",
        "expanded_normalization_shifts", "centering_counts",
    )
    _require_keys(value, keys, context)
    scalar_values = {
        name: _exact_int(value[name], f"{context}.{name}", minimum=0)
        for name in _CENSUS_SCALAR_FIELDS
    }
    ambiguous = _parse_int_list(
        value["raw_ambiguous_spacegroups"],
        f"{context}.raw_ambiguous_spacegroups",
    )
    return DataHallCensus(
        **scalar_values,
        raw_ambiguous_spacegroups=ambiguous,
        hall_to_source_shifts=_parse_distribution(
            value["hall_to_source_shifts"], f"{context}.hall_to_source_shifts"
        ),
        hall_to_source_cosets=_parse_distribution(
            value["hall_to_source_cosets"], f"{context}.hall_to_source_cosets"
        ),
        expanded_normalization_shifts=_parse_distribution(
            value["expanded_normalization_shifts"],
            f"{context}.expanded_normalization_shifts",
        ),
        centering_counts=_parse_centering_counts(
            value["centering_counts"], f"{context}.centering_counts"
        ),
    )


def _input_descriptors():
    return {
        name: {"path": path, "bytes": byte_count, "sha256": digest}
        for name, path, byte_count, digest in _INPUT_SPECS
    }


def _validate_inputs(value, context: str):
    expected = _input_descriptors()
    _require_keys(value, tuple(expected), context)
    for name, descriptor in value.items():
        _require_keys(descriptor, ("path", "bytes", "sha256"), f"{context}.{name}")
        path = _exact_str(descriptor["path"], f"{context}.{name}.path", nonempty=True)
        byte_count = _exact_int(
            descriptor["bytes"], f"{context}.{name}.bytes", minimum=0
        )
        digest = _exact_str(descriptor["sha256"], f"{context}.{name}.sha256")
        if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
            _schema(f"{context}.{name}.sha256", "must be lowercase hexadecimal SHA-256")
        if (path, byte_count, digest) != (
            expected[name]["path"], expected[name]["bytes"], expected[name]["sha256"]
        ):
            raise DataHallIntegrityError(f"{context}.{name}: commitment mismatch")


def _parse_artifact(value):
    keys = (
        "schema", "translation_denominator", "frame_semantics",
        "mapping_semantics", "inputs", "census", "spacegroups",
    )
    _require_keys(value, keys, "artifact")
    if _exact_str(value["schema"], "artifact.schema") != SCHEMA:
        _schema("artifact.schema", "schema mismatch")
    if _exact_int(value["translation_denominator"], "artifact.translation_denominator") != TRANSLATION_DENOMINATOR:
        _invariant("artifact.translation_denominator", "denominator mismatch")
    if _exact_str(value["frame_semantics"], "artifact.frame_semantics") != FRAME_SEMANTICS:
        _invariant("artifact.frame_semantics", "frame semantics mismatch")
    _require_keys(value["mapping_semantics"], ("source_to_hall", "hall_to_source"), "artifact.mapping_semantics")
    mapping_semantics = {
        name: _exact_str(
            value["mapping_semantics"][name],
            f"artifact.mapping_semantics.{name}",
        )
        for name in ("source_to_hall", "hall_to_source")
    }
    if mapping_semantics != _MAPPING_SEMANTICS:
        _invariant("artifact.mapping_semantics", "mapping semantics mismatch")
    _validate_inputs(value["inputs"], "artifact.inputs")

    raw_frames = value["spacegroups"]
    if type(raw_frames) is not list or len(raw_frames) != 230:
        _schema("artifact.spacegroups", "must be a 230-entry array")
    frames = tuple(
        _parse_frame(record, f"artifact.spacegroups[{index}]")
        for index, record in enumerate(raw_frames)
    )
    census = _parse_census(value["census"], "artifact.census")
    return DataHallProvenanceDatabase(frames, census)


def _parse_manifest(value, artifact_bytes: bytes):
    keys = ("schema", "generator_schema_version", "inputs", "artifact")
    _require_keys(value, keys, "manifest")
    if _exact_str(value["schema"], "manifest.schema") != MANIFEST_SCHEMA:
        _schema("manifest.schema", "schema mismatch")
    if _exact_str(value["generator_schema_version"], "manifest.generator_schema_version") != "1":
        _schema("manifest.generator_schema_version", "version mismatch")
    _validate_inputs(value["inputs"], "manifest.inputs")
    if type(value["artifact"]) is not dict:
        _schema("manifest.artifact", "must be an exact object")
    _require_keys(value["artifact"], ("path", "bytes", "sha256"), "manifest.artifact")
    path = _exact_str(value["artifact"]["path"], "manifest.artifact.path")
    if path != _ARTIFACT_NAME:
        _invariant("manifest.artifact.path", "artifact name mismatch")
    byte_count = _exact_int(value["artifact"]["bytes"], "manifest.artifact.bytes", minimum=0)
    digest = _exact_str(value["artifact"]["sha256"], "manifest.artifact.sha256")
    if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
        _schema("manifest.artifact.sha256", "must be lowercase hexadecimal SHA-256")
    if value["inputs"] != _input_descriptors():
        raise DataHallIntegrityError("manifest input commitments differ")
    if byte_count != len(artifact_bytes):
        raise DataHallIntegrityError("manifest artifact byte length differs")
    if digest != hashlib.sha256(artifact_bytes).hexdigest():
        raise DataHallIntegrityError("manifest artifact SHA-256 differs")


def _frame_to_json(frame: DataHallFrame):
    return {
        "spacegroup": frame.spacegroup,
        "source_symbol": frame.source_symbol,
        "centering": frame.centering,
        "pir_anchor_irnumber": frame.pir_anchor_irnumber,
        "cir_anchor_irnumber": frame.cir_anchor_irnumber,
        "raw_candidate_halls": list(frame.raw_candidate_halls),
        "data_hall": frame.data_hall,
        "basis": list(frame.basis),
        "origin_numerator": list(frame.origin_numerator),
        "source_operation_count": frame.source_operation_count,
        "hall_operation_count": frame.hall_operation_count,
        "source_to_hall": [
            {
                "source_operation_index": mapping.source_operation_index,
                "hall_operation_index": mapping.hall_operation_index,
                "lattice_shift_numerator": list(mapping.shift_numerator),
            }
            for mapping in frame.source_to_hall
        ],
        "hall_to_source": [
            {
                "hall_operation_index": mapping.hall_operation_index,
                "source_operation_index": mapping.source_operation_index,
                "lattice_shift_numerator": list(mapping.shift_numerator),
            }
            for mapping in frame.hall_to_source
        ],
    }


def _census_to_json(census: DataHallCensus):
    def distribution(rows):
        return [[list(key), count] for key, count in rows]

    return {
        "pir_records": census.pir_records,
        "cir_records": census.cir_records,
        "source_representatives": census.source_representatives,
        "raw_unique": census.raw_unique,
        "raw_ambiguous": census.raw_ambiguous,
        "raw_missing": census.raw_missing,
        "raw_ambiguous_spacegroups": list(census.raw_ambiguous_spacegroups),
        "filtered_unique": census.filtered_unique,
        "filtered_ambiguous": census.filtered_ambiguous,
        "filtered_missing": census.filtered_missing,
        "selected_hall_operations": census.selected_hall_operations,
        "source_to_hall": census.source_to_hall,
        "source_to_hall_nonzero": census.source_to_hall_nonzero,
        "hall_to_source": census.hall_to_source,
        "hall_to_source_nonzero": census.hall_to_source_nonzero,
        "hall_to_source_shifts": distribution(census.hall_to_source_shifts),
        "hall_to_source_cosets": distribution(census.hall_to_source_cosets),
        "expanded_normalization_nonzero": census.expanded_normalization_nonzero,
        "expanded_normalization_shifts": distribution(census.expanded_normalization_shifts),
        "centering_counts": [list(row) for row in census.centering_counts],
    }


def _project_database(database: DataHallProvenanceDatabase):
    return {
        "schema": SCHEMA,
        "translation_denominator": TRANSLATION_DENOMINATOR,
        "frame_semantics": FRAME_SEMANTICS,
        "mapping_semantics": dict(_MAPPING_SEMANTICS),
        "inputs": _input_descriptors(),
        "census": _census_to_json(database.census),
        "spacegroups": [_frame_to_json(frame) for frame in database.frames],
    }


def _parse_and_validate_pair(artifact_bytes: bytes, manifest_bytes: bytes):
    """Parse a canonical pair for the fixed loader's internal graph checks."""

    artifact = _parse_canonical_json(artifact_bytes, "artifact")
    manifest = _parse_canonical_json(manifest_bytes, "manifest")
    database = _parse_artifact(artifact)
    _parse_manifest(manifest, artifact_bytes)
    if manifest["inputs"] != artifact["inputs"]:
        raise DataHallInvariantError("artifact and manifest inputs differ")
    if _project_database(database) != artifact:
        raise DataHallInvariantError("typed graph projection differs from artifact")
    return database


def _read_verified(path: Path, expected_length: int, expected_sha256: str,
                   context: str) -> bytes:
    try:
        payload = path.read_bytes()
    except (OSError, RuntimeError, ValueError, UnicodeError) as error:
        raise DataHallIntegrityError(f"unable to read {context}") from error
    if type(payload) is not bytes:
        raise DataHallIntegrityError(f"{context} did not return exact bytes")
    if len(payload) != expected_length:
        raise DataHallIntegrityError(f"{context} byte length mismatch")
    if hashlib.sha256(payload).hexdigest() != expected_sha256:
        raise DataHallIntegrityError(f"{context} SHA-256 mismatch")
    return payload


def _lookup_index(spacegroup) -> int:
    if type(spacegroup) is not int:
        raise DataHallLookupError("spacegroup must be an exact integer")
    if not 1 <= spacegroup <= 230:
        raise DataHallLookupError("spacegroup must be in 1..230")
    return spacegroup - 1


_DATABASE = None
_CACHE_LOCK = threading.Lock()


def load_committed_data_hall_provenance() -> DataHallProvenanceDatabase:
    """Load the one fixed sidecar with single-flight caching."""

    global _DATABASE
    cached = _DATABASE
    if cached is not None:
        return cached
    with _CACHE_LOCK:
        cached = _DATABASE
        if cached is not None:
            return cached
        try:
            # Both fixed byte gates complete before either JSON parser runs.
            artifact_bytes = _read_verified(
                _ARTIFACT_PATH, ARTIFACT_BYTE_LENGTH, ARTIFACT_SHA256,
                "artifact",
            )
            manifest_bytes = _read_verified(
                _MANIFEST_PATH, MANIFEST_BYTE_LENGTH, MANIFEST_SHA256,
                "manifest",
            )
            database = _parse_and_validate_pair(artifact_bytes, manifest_bytes)
        except DataHallProvenanceError:
            raise
        except Exception as error:
            raise DataHallInvariantError("committed sidecar validation failed") from error
        _DATABASE = database
        return database


def data_hall_frame(spacegroup) -> DataHallFrame:
    """Return the immutable pinned frame for one space-group number."""

    index = _lookup_index(spacegroup)
    return load_committed_data_hall_provenance().frames[index]


def data_hall_for_spacegroup(spacegroup) -> int:
    """Return the pinned Hall number for one space-group number."""

    return data_hall_frame(spacegroup).data_hall


def _reset_cache_for_test():
    """Clear the singleton for deterministic isolated tests."""

    global _DATABASE
    with _CACHE_LOCK:
        _DATABASE = None


__all__ = [
    "ARTIFACT_BYTE_LENGTH", "ARTIFACT_SHA256", "DataHallCensus",
    "DataHallFrame", "DataHallIntegrityError", "DataHallInvariantError",
    "DataHallLookupError", "DataHallProvenanceDatabase",
    "DataHallProvenanceError", "DataHallSchemaError", "FRAME_SEMANTICS",
    "HallToSourceMapping", "MANIFEST_BYTE_LENGTH", "MANIFEST_SHA256",
    "SCHEMA", "SourceToHallMapping", "data_hall_for_spacegroup",
    "data_hall_frame", "load_committed_data_hall_provenance",
]
