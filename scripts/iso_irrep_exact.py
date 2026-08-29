#!/usr/bin/env python3
"""Strict, immutable loader for the pinned ISO-IR PIR/CIR source frames.

This module is deliberately a source-frame loader rather than an irrep
generator.  It reads the two pinned ZIP members, retains exact rational Seitz
and k-vector data, and validates (but does not materialise) the archived
matrix tokens.  In particular, no Hall choice, affine search, or phase
correction is performed here.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
import hashlib
import io
from pathlib import Path
import re
import threading
from typing import Iterable, Optional, Sequence, Tuple
import zipfile
import zlib

try:  # ``scripts`` is normally imported as a namespace package.
    from .generate_irrep_data import (
        CIR_COMPLEX_TOKEN_SPELLINGS,
        PIR_MATRIX_TOKEN_SPELLINGS,
    )
except ImportError:  # pragma: no cover - useful when run from scripts/.
    from generate_irrep_data import (
        CIR_COMPLEX_TOKEN_SPELLINGS,
        PIR_MATRIX_TOKEN_SPELLINGS,
    )


# Resolve the source root exactly once at import time.  The production entry
# point never consults an extracted directory or a process-current directory.
_MODULE_DIR = Path(__file__).resolve().parent
_REPOSITORY_ROOT = _MODULE_DIR.parent
_ISO_DIR = _REPOSITORY_ROOT / "isotropy_subgroup"
PIR_ARCHIVE_PATH = _ISO_DIR / "PIR_data.zip"
CIR_ARCHIVE_PATH = _ISO_DIR / "CIR_data.zip"

PIR_ARCHIVE_BYTES = 1_235_319
CIR_ARCHIVE_BYTES = 1_555_153
PIR_ARCHIVE_SHA256 = (
    "e909a4f0121688b0590ccaec10b0276171bc24619cf7eb562ba441268c01e121"
)
CIR_ARCHIVE_SHA256 = (
    "f4edcb2852b83a86d1b58f29fb862d9124a227cfc90f9e1ae17d2c97585264e6"
)

PIR_RECORD_COUNT = 10_294
CIR_RECORD_COUNT = 11_202
PIR_IRTRANSLATION_ROW_COUNT = 64_588
CIR_IRTRANSLATION_ROW_COUNT = 68_612
EXPECTED_DENOMINATORS = frozenset({1, 2, 3, 4, 6})
EXPECTED_CENTERING_COUNTS = {
    "P": 149,
    "A": 4,
    "B": 0,
    "C": 16,
    "F": 16,
    "I": 38,
    "R": 7,
}

_PIR_TITLES = (
    "ISO-IR: Physically Irreducible Representations of the 230 Crystallographic Space Groups",
    "2011 Version",
    "Harold T. Stokes and Branton J. Campbell, 2022",
)
_CIR_TITLES = (
    "ISO-IR: Complex Irreducible Representations of the 230 Crystallographic Space Groups",
    "2011 Version",
    "Harold T. Stokes and Branton J. Campbell, 2022",
)


class IsoIrrepExactError(ValueError):
    """Base class for strict source loading failures."""


class ArchiveIntegrityError(IsoIrrepExactError):
    """A pinned archive or its authoritative member is not trustworthy."""


class SourceSchemaError(IsoIrrepExactError):
    """A source member does not obey the official textual grammar."""


class SourceInvariantError(IsoIrrepExactError):
    """A well-formed source violates a cross-record invariant."""


class SourceLookupError(IsoIrrepExactError):
    """A source-universe lookup has an invalid type, range, or target."""


class SourceArchive(Enum):
    PIR = "PIR"
    CIR = "CIR"


class Centering(Enum):
    P = "P"
    A = "A"
    B = "B"
    C = "C"
    F = "F"
    I = "I"
    R = "R"


Int3 = Tuple[int, int, int]
Fraction3 = Tuple[Fraction, Fraction, Fraction]
Rotation3 = Tuple[Int3, Int3, Int3]
RawAugmented = Tuple[int, int, int, int, int, int, int, int,
                     int, int, int, int, int, int, int, int]
RawKVector = RawAugmented
RawIrTranslation = Tuple[int, int, int, int]
OptionalFraction3 = Optional[Fraction3]


def _require_tuple(value, length: Optional[int], context: str) -> None:
    if type(value) is not tuple:
        raise TypeError(f"{context} must be an immutable tuple")
    if length is not None and len(value) != length:
        raise TypeError(f"{context} must have length {length}")


@dataclass(frozen=True)
class ExactSeitz:
    __slots__ = ("rotation", "translation", "raw_augmented")

    rotation: Rotation3
    translation: Fraction3
    raw_augmented: RawAugmented

    def __post_init__(self):
        _require_tuple(self.rotation, 3, "ExactSeitz.rotation")
        if any(type(row) is not tuple or len(row) != 3 for row in self.rotation):
            raise TypeError("ExactSeitz.rotation must be a 3x3 tuple")
        if any(type(value) is not int for row in self.rotation for value in row):
            raise TypeError("ExactSeitz.rotation entries must be int")
        _require_tuple(self.translation, 3, "ExactSeitz.translation")
        if any(type(value) is not Fraction for value in self.translation):
            raise TypeError("ExactSeitz.translation entries must be Fraction")
        _require_tuple(self.raw_augmented, 16, "ExactSeitz.raw_augmented")
        if any(type(value) is not int for value in self.raw_augmented):
            raise TypeError("ExactSeitz.raw_augmented entries must be int")


@dataclass(frozen=True)
class ExactKArm:
    __slots__ = ("constant", "parameters", "raw_augmented")

    constant: Fraction3
    parameters: tuple[OptionalFraction3, OptionalFraction3, OptionalFraction3]
    raw_augmented: RawKVector

    def __post_init__(self):
        _require_tuple(self.constant, 3, "ExactKArm.constant")
        if any(type(value) is not Fraction for value in self.constant):
            raise TypeError("ExactKArm.constant entries must be Fraction")
        _require_tuple(self.parameters, 3, "ExactKArm.parameters")
        for parameter in self.parameters:
            if parameter is not None:
                _require_tuple(parameter, 3, "ExactKArm parameter")
                if any(type(value) is not Fraction for value in parameter):
                    raise TypeError("ExactKArm parameter entries must be Fraction")
        _require_tuple(self.raw_augmented, 16, "ExactKArm.raw_augmented")
        if any(type(value) is not int for value in self.raw_augmented):
            raise TypeError("ExactKArm.raw_augmented entries must be int")


@dataclass(frozen=True)
class ExactIrTranslation:
    __slots__ = ("vector", "raw")

    vector: Fraction3
    raw: RawIrTranslation

    def __post_init__(self):
        _require_tuple(self.vector, 3, "ExactIrTranslation.vector")
        if any(type(value) is not Fraction for value in self.vector):
            raise TypeError("ExactIrTranslation.vector entries must be Fraction")
        _require_tuple(self.raw, 4, "ExactIrTranslation.raw")
        if any(type(value) is not int for value in self.raw):
            raise TypeError("ExactIrTranslation.raw entries must be int")


@dataclass(frozen=True)
class ExactSourceRecord:
    __slots__ = (
        "archive", "irnumber", "spacegroup", "space_group_symbol",
        "centering", "irrep_label", "dimension", "irtype", "kcount",
        "pmkcount", "k_arms", "operations", "irtranslations",
    )

    archive: SourceArchive
    irnumber: int
    spacegroup: int
    space_group_symbol: str
    centering: Centering
    irrep_label: str
    dimension: int
    irtype: int
    kcount: int
    pmkcount: int
    k_arms: tuple[ExactKArm, ...]
    operations: tuple[ExactSeitz, ...]
    irtranslations: tuple[Optional[ExactIrTranslation], ...]

    def __post_init__(self):
        if not isinstance(self.archive, SourceArchive):
            raise TypeError("ExactSourceRecord.archive must be SourceArchive")
        if type(self.irnumber) is not int or type(self.spacegroup) is not int:
            raise TypeError("ExactSourceRecord source numbers must be int")
        if type(self.space_group_symbol) is not str or type(self.irrep_label) is not str:
            raise TypeError("ExactSourceRecord symbols and labels must be str")
        if not isinstance(self.centering, Centering):
            raise TypeError("ExactSourceRecord.centering must be Centering")
        for name in ("dimension", "irtype", "kcount", "pmkcount"):
            if type(getattr(self, name)) is not int:
                raise TypeError(f"ExactSourceRecord.{name} must be int")
        _require_tuple(self.k_arms, None, "ExactSourceRecord.k_arms")
        _require_tuple(self.operations, None, "ExactSourceRecord.operations")
        _require_tuple(self.irtranslations, None, "ExactSourceRecord.irtranslations")
        if not self.k_arms or not self.operations:
            raise ValueError("ExactSourceRecord must contain k arms and operations")
        if len(self.irtranslations) != len(self.operations):
            raise ValueError("ExactSourceRecord translation/op count mismatch")
        expected_arms = self.pmkcount if self.archive is SourceArchive.PIR else self.kcount
        if len(self.k_arms) != expected_arms:
            raise ValueError("ExactSourceRecord k-arm count mismatch")
        if any(not isinstance(arm, ExactKArm) for arm in self.k_arms):
            raise TypeError("ExactSourceRecord.k_arms entries must be ExactKArm")
        if any(not isinstance(operation, ExactSeitz) for operation in self.operations):
            raise TypeError("ExactSourceRecord.operations entries must be ExactSeitz")
        if any(
            translation is not None and not isinstance(translation, ExactIrTranslation)
            for translation in self.irtranslations
        ):
            raise TypeError("ExactSourceRecord.irtranslations entries must be ExactIrTranslation or None")

    @property
    def opcount(self) -> int:
        """The parsed ninth header field, represented by operation cardinality."""

        return len(self.operations)

    @property
    def special(self) -> bool:
        """Whether the official first-arm ``kspecial`` predicate is true."""

        return all(
            self.k_arms[0].raw_augmented[offset] == 0
            for offset in (4, 5, 6, 8, 9, 10, 12, 13, 14)
        )


@dataclass(frozen=True)
class ExactSpaceGroupUniverse:
    __slots__ = (
        "spacegroup", "space_group_symbol", "centering", "operations",
        "pir_irnumbers", "cir_irnumbers",
    )

    spacegroup: int
    space_group_symbol: str
    centering: Centering
    operations: tuple[ExactSeitz, ...]
    pir_irnumbers: tuple[int, ...]
    cir_irnumbers: tuple[int, ...]

    def __post_init__(self):
        if type(self.spacegroup) is not int:
            raise TypeError("ExactSpaceGroupUniverse.spacegroup must be int")
        if type(self.space_group_symbol) is not str:
            raise TypeError("ExactSpaceGroupUniverse.space_group_symbol must be str")
        if not isinstance(self.centering, Centering):
            raise TypeError("ExactSpaceGroupUniverse.centering must be Centering")
        _require_tuple(self.operations, None, "ExactSpaceGroupUniverse.operations")
        _require_tuple(self.pir_irnumbers, None, "ExactSpaceGroupUniverse.pir_irnumbers")
        _require_tuple(self.cir_irnumbers, None, "ExactSpaceGroupUniverse.cir_irnumbers")
        if any(not isinstance(operation, ExactSeitz) for operation in self.operations):
            raise TypeError("ExactSpaceGroupUniverse.operations entries must be ExactSeitz")
        if any(type(number) is not int for number in self.pir_irnumbers + self.cir_irnumbers):
            raise TypeError("ExactSpaceGroupUniverse irnumbers must be int")


@dataclass(frozen=True)
class ExactIsoIrrepDatabase:
    __slots__ = ("pir_records", "cir_records", "universes")

    pir_records: tuple[ExactSourceRecord, ...]
    cir_records: tuple[ExactSourceRecord, ...]
    universes: tuple[Optional[ExactSpaceGroupUniverse], ...]

    def __post_init__(self):
        _require_tuple(self.pir_records, None, "ExactIsoIrrepDatabase.pir_records")
        _require_tuple(self.cir_records, None, "ExactIsoIrrepDatabase.cir_records")
        _require_tuple(self.universes, 231, "ExactIsoIrrepDatabase.universes")
        if any(
            not isinstance(record, ExactSourceRecord)
            or record.archive is not SourceArchive.PIR
            for record in self.pir_records
        ):
            raise TypeError("ExactIsoIrrepDatabase.pir_records entries must be PIR records")
        if any(
            not isinstance(record, ExactSourceRecord)
            or record.archive is not SourceArchive.CIR
            for record in self.cir_records
        ):
            raise TypeError("ExactIsoIrrepDatabase.cir_records entries must be CIR records")
        if self.universes[0] is not None or any(
            universe is not None and not isinstance(universe, ExactSpaceGroupUniverse)
            for universe in self.universes
        ):
            raise TypeError("ExactIsoIrrepDatabase.universes entries are malformed")

    def source_universe(self, spacegroup: int) -> ExactSpaceGroupUniverse:
        return _lookup_universe(self.universes, spacegroup)


_SIGNED_INTEGER_RE = re.compile(r"(?:0|-[1-9][0-9]*|[1-9][0-9]*)\Z", re.ASCII)
_UNSIGNED_INTEGER_RE = re.compile(r"(?:0|[1-9][0-9]*)\Z", re.ASCII)
_CIR_COMPLEX_TOKEN_RE = re.compile(r"\(([^,]+),([^\)]+)\)\Z", re.ASCII)
_IDENTITY_ROTATION = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
_ZERO_TRANSLATION = (Fraction(0), Fraction(0), Fraction(0))


def _error(error_type: type[IsoIrrepExactError], message: str) -> IsoIrrepExactError:
    return error_type(message)


def _parse_integer(token: str, *, context: str, unsigned: bool = False) -> int:
    """Parse one canonical ASCII integer and nothing accepted by ``int`` more."""

    pattern = _UNSIGNED_INTEGER_RE if unsigned else _SIGNED_INTEGER_RE
    if pattern.fullmatch(token) is None:
        raise SourceSchemaError(f"non-canonical integer {token!r} for {context}")
    # The regular expression has already excluded Unicode digits and signs
    # that Python's int() would otherwise accept.
    try:
        return int(token)
    except ValueError as error:
        # Keep even a syntactically canonical but excessively long integer
        # inside the typed source-schema error boundary (Python may reject it
        # because of its interpreter-wide integer conversion limit).
        raise SourceSchemaError(f"invalid integer {token!r} for {context}") from error


def _validate_ascii_text(text: str, context: str, *, require_final_lf: bool) -> None:
    """Require source text to be printable ASCII records separated by LF."""

    if type(text) is not str:
        raise TypeError(f"{context} must be str")
    if require_final_lf and not text.endswith("\n"):
        raise SourceSchemaError(f"{context} is missing its final LF")
    for position, character in enumerate(text):
        if character == "\n":
            continue
        codepoint = ord(character)
        if codepoint < 0x20 or codepoint > 0x7E:
            raise SourceSchemaError(
                f"{context} contains non-printable/non-ASCII character "
                f"U+{codepoint:04X} at offset {position}"
            )


def _source_lines(text: str, context: str) -> tuple[str, ...]:
    _validate_ascii_text(text, context, require_final_lf=True)
    # ``splitlines`` would silently accept CRLF and Unicode line separators;
    # this explicit split is part of the official source grammar.
    pieces = text.split("\n")
    if not pieces or pieces[-1] != "":  # defensive after the final-LF check
        raise SourceSchemaError(f"{context} is missing its final LF")
    return tuple(pieces[:-1])


def _parse_fixed_unsigned(field: str, width: int, name: str, archive: SourceArchive, line_number: int) -> int:
    if len(field) != width:
        raise SourceSchemaError(
            f"{archive.value} header {name} is not width {width} at line {line_number}"
        )
    digits = field.lstrip(" ")
    value = _parse_integer(
        digits,
        context=f"{archive.value} header {name} at line {line_number}",
        unsigned=True,
    )
    if field != str(value).rjust(width, " "):
        raise SourceSchemaError(
            f"non-canonical {archive.value} header {name} at line {line_number}"
        )
    return value


def _parse_padded_field(field: str, width: int, name: str, archive: SourceArchive, line_number: int) -> str:
    if len(field) != width:
        raise SourceSchemaError(
            f"{archive.value} header {name} is not width {width} at line {line_number}"
        )
    if any(not (0x20 <= ord(character) <= 0x7E) for character in field):
        raise SourceSchemaError(
            f"{archive.value} header {name} is not printable ASCII at line {line_number}"
        )
    semantic = field.rstrip(" ")
    if not semantic or field[0] == " " or '"' in semantic:
        raise SourceSchemaError(
            f"invalid {archive.value} header {name} padding/content at line {line_number}"
        )
    if field != semantic + " " * (width - len(semantic)):
        raise SourceSchemaError(
            f"invalid {archive.value} header {name} right-padding at line {line_number}"
        )
    return semantic


def _parse_header(line: str, archive: SourceArchive, line_number: int) -> tuple[int, int, str, str, int, int, int, int, int]:
    """Parse the exact ``(i5,i4,a,5i3)`` official writer output."""

    if type(line) is not str:
        raise TypeError("header line must be str")
    expected_length = 48 if archive is SourceArchive.PIR else 44
    if len(line) != expected_length:
        raise SourceSchemaError(
            f"malformed {archive.value} header at line {line_number}: "
            f"expected {expected_length} characters, got {len(line)}"
        )
    _validate_ascii_text(line, f"{archive.value} header line {line_number}", require_final_lf=False)
    irnumber = _parse_fixed_unsigned(line[0:5], 5, "irnumber", archive, line_number)
    spacegroup = _parse_fixed_unsigned(line[5:9], 4, "spacegroup", archive, line_number)
    label_width = 8 if archive is SourceArchive.PIR else 4
    if line[9:11] != ' "' or line[21:24] != '" "':
        raise SourceSchemaError(f"malformed {archive.value} header literals at line {line_number}")
    raw_symbol = line[11:21]
    symbol = _parse_padded_field(raw_symbol, 10, "space-group symbol", archive, line_number)
    label_start = 24
    label_end = label_start + label_width
    if line[label_end] != '"':
        raise SourceSchemaError(f"malformed {archive.value} header label quote at line {line_number}")
    raw_label = line[label_start:label_end]
    label = _parse_padded_field(raw_label, label_width, "irrep label", archive, line_number)
    numeric_start = label_end + 1
    values = []
    for offset, field in enumerate(("dimension", "irtype", "kcount", "pmkcount", "opcount")):
        start = numeric_start + offset * 3
        values.append(_parse_fixed_unsigned(line[start:start + 3], 3, field, archive, line_number))
    if line != (
        f"{irnumber:5d}{spacegroup:4d} \"{raw_symbol}\" \"{raw_label}\""
        + "".join(str(value).rjust(3, " ") for value in values)
    ):
        raise SourceSchemaError(f"non-canonical {archive.value} header at line {line_number}")
    dimension, irtype, kcount, pmkcount, opcount = values

    if not 1 <= irnumber:
        raise SourceSchemaError(f"irnumber {irnumber} must be positive at line {line_number}")
    try:
        Centering(raw_symbol[0])
    except ValueError as error:
        raise SourceSchemaError(
            f"unknown centering {raw_symbol[0]!r} at line {line_number}"
        ) from error
    if not 1 <= spacegroup <= 230:
        raise SourceSchemaError(f"spacegroup {spacegroup} outside 1..230 at line {line_number}")
    if not 1 <= dimension <= 48:
        raise SourceSchemaError(f"dimension {dimension} outside 1..48 at line {line_number}")
    if irtype not in (1, 2, 3):
        raise SourceSchemaError(f"irtype {irtype} outside {{1,2,3}} at line {line_number}")
    for field, value in (("kcount", kcount), ("pmkcount", pmkcount), ("opcount", opcount)):
        if not 1 <= value <= 48:
            raise SourceSchemaError(f"{field} {value} outside 1..48 at line {line_number}")
    if kcount not in (pmkcount, 2 * pmkcount):
        raise SourceSchemaError(
            f"kcount={kcount}, pmkcount={pmkcount} violates star relation at line {line_number}"
        )
    return (irnumber, spacegroup, symbol, label, dimension, irtype, kcount, pmkcount, opcount)


def _line_tokens(lines: Sequence[str], index: int, *, context: str) -> list[str]:
    if index >= len(lines):
        raise SourceSchemaError(f"truncated {context}")
    line = lines[index]
    if not line.strip():
        raise SourceSchemaError(f"blank line {index + 1} in {context}")
    _validate_ascii_text(line, f"{context} line {index + 1}", require_final_lf=False)
    # The Fortran writer separates values with ordinary ASCII spaces.  Do not
    # let str.split() silently accept tabs, NBSP, or other Unicode separators.
    return [token for token in line.split(" ") if token]


def _read_exact_block(
    lines: Sequence[str],
    start: int,
    count: int,
    *,
    context: str,
    parser,
) -> tuple[tuple, int]:
    """Read exactly ``count`` tokens, rejecting terminal-line overrun/blank."""

    values = []
    index = start
    while len(values) < count:
        tokens = _line_tokens(lines, index, context=context)
        remaining = count - len(values)
        if len(tokens) > remaining:
            raise SourceSchemaError(
                f"extra tokens at line {index + 1} in {context}: "
                f"expected {remaining}, got {len(tokens)}"
            )
        values.extend(parser(token, index + 1, context) for token in tokens)
        index += 1
    return tuple(values), index


def _skip_exact_block(
    lines: Sequence[str],
    start: int,
    count: int,
    *,
    context: str,
    parser,
) -> int:
    """Validate and skip exactly ``count`` tokens without retaining a block."""

    consumed = 0
    index = start
    while consumed < count:
        tokens = _line_tokens(lines, index, context=context)
        remaining = count - consumed
        if len(tokens) > remaining:
            raise SourceSchemaError(
                f"extra tokens at line {index + 1} in {context}: "
                f"expected {remaining}, got {len(tokens)}"
            )
        for token in tokens:
            parser(token, index + 1, context)
        consumed += len(tokens)
        index += 1
    return index


def _read_exact_row(
    lines: Sequence[str],
    start: int,
    count: int,
    *,
    context: str,
    parser,
) -> tuple[tuple, int]:
    tokens = _line_tokens(lines, start, context=context)
    if len(tokens) != count:
        raise SourceSchemaError(
            f"{context} line {start + 1} has {len(tokens)} tokens, expected {count}"
        )
    return tuple(parser(token, start + 1, context) for token in tokens), start + 1


def _parse_payload_integer(token: str, line_number: int, context: str) -> int:
    return _parse_integer(token, context=f"{context} line {line_number}")


def _parse_matrix_token(token: str, archive: SourceArchive, line_number: int, context: str) -> str:
    if archive is SourceArchive.PIR:
        if token not in PIR_MATRIX_TOKEN_SPELLINGS:
            raise SourceSchemaError(
                f"unknown PIR matrix token {token!r} at line {line_number} for {context}"
            )
        return token
    if token not in CIR_COMPLEX_TOKEN_SPELLINGS:
        raise SourceSchemaError(
            f"unknown CIR complex token {token!r} at line {line_number} for {context}"
        )
    match = _CIR_COMPLEX_TOKEN_RE.fullmatch(token)
    if match is None:
        raise SourceSchemaError(
            f"malformed CIR complex token {token!r} at line {line_number} for {context}"
        )
    real, imaginary = match.groups()
    if real not in PIR_MATRIX_TOKEN_SPELLINGS or imaginary not in PIR_MATRIX_TOKEN_SPELLINGS:
        raise SourceSchemaError(
            f"unknown CIR complex component in {token!r} at line {line_number} for {context}"
        )
    return token


def _parse_k_arms(raw_values: tuple[int, ...], arm_count: int, archive: SourceArchive, context: str) -> tuple[ExactKArm, ...]:
    arms = []
    for arm_index in range(arm_count):
        base = arm_index * 16
        raw = tuple(raw_values[base:base + 16])
        constant_denominator = raw[3]
        if constant_denominator <= 0:
            raise SourceSchemaError(
                f"invalid {archive.value} constant k denominator for {context} arm {arm_index}"
            )
        constant = tuple(
            Fraction(raw[offset], constant_denominator) for offset in (0, 1, 2)
        )
        parameters = []
        for numerator_offset, denominator_offset in ((4, 7), (8, 11), (12, 15)):
            denominator = raw[denominator_offset]
            numerators = raw[numerator_offset:numerator_offset + 3]
            if denominator == 0:
                if any(numerators):
                    raise SourceSchemaError(
                        f"zero {archive.value} parameter denominator with nonzero numerator "
                        f"for {context} arm {arm_index}"
                    )
                parameters.append(None)
            elif denominator > 0:
                parameters.append(tuple(Fraction(n, denominator) for n in numerators))
            else:
                raise SourceSchemaError(
                    f"negative {archive.value} parameter denominator for {context} arm {arm_index}"
                )
        arms.append(
            ExactKArm(
                constant=constant,
                parameters=tuple(parameters),  # type: ignore[arg-type]
                raw_augmented=raw,  # type: ignore[arg-type]
            )
        )
    return tuple(arms)


def _parse_operation(raw: tuple[int, ...], archive: SourceArchive, context: str) -> ExactSeitz:
    denominator = raw[15]
    if denominator <= 0:
        raise SourceSchemaError(f"invalid {archive.value} operation denominator for {context}")
    rotation_offsets = (0, 1, 2, 4, 5, 6, 8, 9, 10)
    if any(raw[offset] % denominator for offset in rotation_offsets):
        raise SourceSchemaError(f"rotation is not divisible by denominator for {context}")
    if raw[12:15] != (0, 0, 0):
        raise SourceSchemaError(f"invalid augmented bottom row for {context}")
    rotation = (
        tuple(raw[offset] // denominator for offset in (0, 1, 2)),
        tuple(raw[offset] // denominator for offset in (4, 5, 6)),
        tuple(raw[offset] // denominator for offset in (8, 9, 10)),
    )
    translation = tuple(Fraction(raw[offset], denominator) for offset in (3, 7, 11))
    return ExactSeitz(
        rotation=rotation,  # type: ignore[arg-type]
        translation=translation,  # type: ignore[arg-type]
        raw_augmented=raw,  # type: ignore[arg-type]
    )


def _parse_irtranslation(raw: tuple[int, ...], archive: SourceArchive, context: str) -> ExactIrTranslation:
    denominator = raw[3]
    if denominator <= 0:
        raise SourceSchemaError(f"invalid {archive.value} irtranslation denominator for {context}")
    vector = tuple(Fraction(raw[offset], denominator) for offset in (0, 1, 2))
    return ExactIrTranslation(vector=vector, raw=raw)  # type: ignore[arg-type]


def _rotation_determinant(rotation: Rotation3) -> int:
    a, b, c = rotation
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _rotation_product(left: Rotation3, right: Rotation3) -> Rotation3:
    return tuple(
        tuple(sum(left[row][index] * right[index][column] for index in range(3)) for column in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def _validate_operations(
    operations: Sequence[ExactSeitz],
    spacegroup: int,
    label: str,
    *,
    check_closure: bool,
) -> None:
    """Validate GL(3,Z), ordered identity, uniqueness, and rotation closure."""

    context = f"SG{spacegroup} {label!r}"
    if not operations:
        raise SourceInvariantError(f"{context} has no Seitz operations")
    for operation in operations:
        determinant = _rotation_determinant(operation.rotation)
        if determinant not in (-1, 1):
            raise SourceInvariantError(
                f"{context} has a non-GL(3,Z) rotation determinant {determinant}"
            )
    rotations = tuple(operation.rotation for operation in operations)
    rotation_set = set(rotations)
    if len(rotation_set) != len(rotations):
        raise SourceInvariantError(f"{context} contains duplicate rotations")
    first = operations[0]
    if first.rotation != _IDENTITY_ROTATION or first.translation != _ZERO_TRANSLATION:
        raise SourceInvariantError(f"{context} operation slot 0 is not exact {{I|0}}")
    if not check_closure:
        return
    for left in rotation_set:
        for right in rotation_set:
            if _rotation_product(left, right) not in rotation_set:
                raise SourceInvariantError(f"{context} rotation set is not multiplication-closed")


def _parse_records_from_lines(
    lines: Sequence[str],
    archive: SourceArchive,
    *,
    validate_census: bool,
) -> tuple[ExactSourceRecord, ...]:
    titles = _PIR_TITLES if archive is SourceArchive.PIR else _CIR_TITLES
    if len(lines) < 3:
        raise SourceSchemaError(f"truncated {archive.value} title block")
    for index, expected in enumerate(titles):
        if lines[index] != expected:
            raise SourceSchemaError(
                f"unexpected {archive.value} title line {index + 1}: {lines[index]!r}"
            )

    records = []
    seen_keys = set()
    # Cache closure by the ordered rotation tuple, rather than by SG number:
    # a malformed same-SG record with a different operation set must still be
    # checked before the later universe-folding comparison.
    closure_checked = set()
    expected_irnumber = 1
    index = 3
    while index < len(lines):
        if not lines[index].strip():
            raise SourceSchemaError(f"blank line {index + 1} between {archive.value} records")
        (
            irnumber,
            spacegroup,
            symbol,
            label,
            dimension,
            irtype,
            kcount,
            pmkcount,
            opcount,
        ) = _parse_header(lines[index], archive, index + 1)
        if irnumber != expected_irnumber:
            raise SourceInvariantError(
                f"unexpected {archive.value} irnumber {irnumber} at line {index + 1}; "
                f"expected {expected_irnumber}"
            )
        expected_irnumber += 1
        key = (spacegroup, label)
        if key in seen_keys:
            raise SourceInvariantError(f"duplicate {archive.value} source key {key!r}")
        seen_keys.add(key)
        try:
            centering = Centering(symbol[0])
        except ValueError as error:
            raise SourceSchemaError(
                f"unknown centering {symbol[0]!r} in {archive.value} symbol {symbol!r}"
            ) from error

        index += 1
        k_payload_count = pmkcount if archive is SourceArchive.PIR else kcount
        divisibility_count = pmkcount if archive is SourceArchive.PIR else kcount
        if dimension % divisibility_count:
            raise SourceInvariantError(
                f"{archive.value} SG{spacegroup} {label!r} dimension {dimension} "
                f"is not divisible by {('pmkcount' if archive is SourceArchive.PIR else 'kcount')} "
                f"{divisibility_count}"
            )
        raw_k, index = _read_exact_block(
            lines,
            index,
            16 * k_payload_count,
            context=f"{archive.value} SG{spacegroup} {label!r} k payload",
            parser=_parse_payload_integer,
        )
        k_arms = _parse_k_arms(
            raw_k,
            k_payload_count,
            archive,
            f"SG{spacegroup} {label!r}",
        )
        special = all(
            raw_k[offset] == 0
            for offset in (4, 5, 6, 8, 9, 10, 12, 13, 14)
        )

        operations = []
        irtranslations = []
        for operation_index in range(opcount):
            context = f"{archive.value} SG{spacegroup} {label!r} operation {operation_index}"
            raw_operation, index = _read_exact_row(
                lines,
                index,
                16,
                context=context,
                parser=_parse_payload_integer,
            )
            operations.append(_parse_operation(raw_operation, archive, context))
            if special:
                irtranslations.append(None)
            else:
                raw_translation, index = _read_exact_row(
                    lines,
                    index,
                    4,
                    context=f"{context} irtranslation",
                    parser=_parse_payload_integer,
                )
                irtranslations.append(_parse_irtranslation(raw_translation, archive, context))
            # Matrix data is an exact source token block, intentionally not
            # materialised in the public model.
            index = _skip_exact_block(
                lines,
                index,
                dimension * dimension,
                context=f"{context} matrix",
                parser=lambda token, line_number, matrix_context: _parse_matrix_token(
                    token, archive, line_number, matrix_context
                ),
            )

        rotation_key = tuple(operation.rotation for operation in operations)
        _validate_operations(
            operations,
            spacegroup,
            label,
            check_closure=rotation_key not in closure_checked,
        )
        closure_checked.add(rotation_key)

        records.append(
            ExactSourceRecord(
                archive=archive,
                irnumber=irnumber,
                spacegroup=spacegroup,
                space_group_symbol=symbol,
                centering=centering,
                irrep_label=label,
                dimension=dimension,
                irtype=irtype,
                kcount=kcount,
                pmkcount=pmkcount,
                k_arms=k_arms,
                operations=tuple(operations),
                irtranslations=tuple(irtranslations),
            )
        )

    if index != len(lines):
        raise SourceSchemaError(f"{archive.value} parser did not reach exact EOF")
    if validate_census:
        expected_count = PIR_RECORD_COUNT if archive is SourceArchive.PIR else CIR_RECORD_COUNT
        if len(records) != expected_count:
            raise SourceInvariantError(
                f"{archive.value} record census mismatch: got {len(records)}, expected {expected_count}"
            )
    return tuple(records)


def parse_exact_source_text(
    text: str,
    archive: SourceArchive,
    *,
    validate_census: bool = False,
) -> tuple[ExactSourceRecord, ...]:
    """Parse source text through the strict record seam.

    ``validate_census=False`` is intended for small synthetic adversarial
    fixtures.  Grammar, token, line, and cross-record checks remain enabled;
    only the pinned full-file record count is relaxed.
    """

    if not isinstance(archive, SourceArchive):
        raise TypeError("archive must be a SourceArchive")
    if type(text) is not str:
        raise TypeError("source text must be str")
    return _parse_records_from_lines(
        _source_lines(text, f"{archive.value} source text"),
        archive,
        validate_census=validate_census,
    )


def parse_exact_source_lines(
    lines: Sequence[str],
    archive: SourceArchive,
    *,
    validate_census: bool = False,
) -> tuple[ExactSourceRecord, ...]:
    """Parse an already split source through the same strict seam."""

    if not isinstance(archive, SourceArchive):
        raise TypeError("archive must be a SourceArchive")
    if not isinstance(lines, (tuple, list)):
        raise TypeError("source lines must be a list or tuple of str")
    if any(type(line) is not str for line in lines):
        raise TypeError("source lines must contain only str")
    for index, line in enumerate(lines):
        _validate_ascii_text(line, f"{archive.value} source line {index + 1}", require_final_lf=False)
        if "\n" in line:
            raise SourceSchemaError(
                f"{archive.value} source line {index + 1} contains an embedded LF"
            )
    return _parse_records_from_lines(lines, archive, validate_census=validate_census)


def _parse_archive_text(
    text: str,
    archive: SourceArchive,
    *,
    validate_census: bool = False,
) -> tuple[ExactSourceRecord, ...]:
    """Compatibility seam used by focused tests; see parse_exact_source_text."""

    return parse_exact_source_text(text, archive, validate_census=validate_census)


def _zip_basename(name: str) -> str:
    normalized = name.replace("\\", "/")
    return normalized.rstrip("/").rsplit("/", 1)[-1]


def _unsafe_zip_name(name: str) -> bool:
    """Reject traversal/absolute archive names before selecting a member."""

    normalized = name.replace("\\", "/")
    if normalized.startswith("/") or normalized.startswith("\\"):
        return True
    if re.match(r"^[A-Za-z]:", normalized):
        return True
    return any(part == ".." for part in normalized.split("/"))


def _read_verified_member(
    archive: SourceArchive,
    *,
    path: Optional[Path] = None,
    expected_bytes: Optional[int] = None,
    expected_sha256: Optional[str] = None,
) -> str:
    if not isinstance(archive, SourceArchive):
        raise TypeError("archive must be a SourceArchive")
    if path is None:
        path = PIR_ARCHIVE_PATH if archive is SourceArchive.PIR else CIR_ARCHIVE_PATH
    elif not hasattr(path, "read_bytes"):
        path = Path(path)
    expected_bytes = (
        PIR_ARCHIVE_BYTES if archive is SourceArchive.PIR else CIR_ARCHIVE_BYTES
    ) if expected_bytes is None else expected_bytes
    expected_sha256 = (
        PIR_ARCHIVE_SHA256 if archive is SourceArchive.PIR else CIR_ARCHIVE_SHA256
    ) if expected_sha256 is None else expected_sha256
    try:
        # This is intentionally the sole path read.  Verification and ZIP
        # parsing both consume this immutable snapshot, avoiding a TOCTOU
        # window between stat/hash and opening the archive again.
        archive_bytes = path.read_bytes()
    except (
        OSError,
        IOError,
        NotImplementedError,
        RuntimeError,
        EOFError,
        ValueError,
        zipfile.BadZipFile,
        zipfile.LargeZipFile,
        KeyError,
        zlib.error,
    ) as error:
        raise ArchiveIntegrityError(f"unable to read {archive.value} archive {path}") from error
    if type(archive_bytes) is not bytes:
        raise ArchiveIntegrityError(f"{archive.value} archive read did not return bytes")
    actual_bytes = len(archive_bytes)
    if actual_bytes != expected_bytes:
        raise ArchiveIntegrityError(
            f"{archive.value} archive byte-size mismatch: expected {expected_bytes}, got {actual_bytes}"
        )
    actual_sha256 = hashlib.sha256(archive_bytes).hexdigest()
    if actual_sha256 != expected_sha256:
        raise ArchiveIntegrityError(
            f"{archive.value} archive SHA-256 mismatch: expected {expected_sha256}, got {actual_sha256}"
        )
    member_name = "PIR_data.txt" if archive is SourceArchive.PIR else "CIR_data.txt"
    try:
        with zipfile.ZipFile(io.BytesIO(archive_bytes)) as zip_file:
            infos = tuple(zip_file.infolist())
            names = tuple(info.filename for info in infos)
            if any(_unsafe_zip_name(name) for name in names):
                raise ArchiveIntegrityError(
                    f"{archive.value} archive contains an unsafe member path"
                )
            if any(
                _zip_basename(name) == member_name and name != member_name
                for name in names
            ):
                raise ArchiveIntegrityError(
                    f"{archive.value} archive contains a non-root {member_name!r} collision"
                )
            matches = tuple(info for info in infos if info.filename == member_name)
            if len(matches) != 1:
                raise ArchiveIntegrityError(
                    f"{archive.value} archive member {member_name!r} matched {len(matches)} entries"
                )
            member_bytes = zip_file.read(matches[0])
    except ArchiveIntegrityError:
        raise
    except (OSError, zipfile.BadZipFile, zipfile.LargeZipFile, KeyError,
            NotImplementedError, RuntimeError, EOFError, ValueError, zlib.error) as error:
        raise ArchiveIntegrityError(
            f"unable to read authoritative {member_name} from {archive.value} archive"
        ) from error
    try:
        return member_bytes.decode("ascii", errors="strict")
    except (UnicodeDecodeError, UnicodeError) as error:
        raise SourceSchemaError(
            f"authoritative {member_name} member is not strict ASCII"
        ) from error


def _record_denominators(records: Iterable[ExactSourceRecord]) -> tuple[set[int], set[int]]:
    operation_denominators = set()
    translation_denominators = set()
    for record in records:
        operation_denominators.update(operation.raw_augmented[15] for operation in record.operations)
        translation_denominators.update(
            translation.raw[3]
            for translation in record.irtranslations
            if translation is not None
        )
    return operation_denominators, translation_denominators


def _build_universes(
    pir_records: Sequence[ExactSourceRecord],
    cir_records: Sequence[ExactSourceRecord],
    *,
    validate_census: bool,
) -> tuple[Optional[ExactSpaceGroupUniverse], ...]:
    grouped: dict[int, dict[SourceArchive, list[ExactSourceRecord]]] = defaultdict(
        lambda: {SourceArchive.PIR: [], SourceArchive.CIR: []}
    )
    for record in tuple(pir_records) + tuple(cir_records):
        grouped[record.spacegroup][record.archive].append(record)

    universes: list[Optional[ExactSpaceGroupUniverse]] = [None]
    for spacegroup in range(1, 231):
        by_archive = grouped.get(spacegroup)
        if by_archive is None:
            if validate_census:
                raise SourceInvariantError(f"missing source records for SG{spacegroup}")
            universes.append(None)
            continue
        all_records = by_archive[SourceArchive.PIR] + by_archive[SourceArchive.CIR]
        baseline = all_records[0] if all_records else None
        if baseline is None:  # defensive; the branch above makes this unreachable.
            universes.append(None)
            continue
        _validate_operations(
            baseline.operations,
            spacegroup,
            baseline.irrep_label,
            check_closure=True,
        )
        for record in all_records[1:]:
            if record.space_group_symbol != baseline.space_group_symbol:
                raise SourceInvariantError(
                    f"SG{spacegroup} source symbols differ: "
                    f"{baseline.space_group_symbol!r} vs {record.space_group_symbol!r}"
                )
            if record.centering is not baseline.centering:
                raise SourceInvariantError(f"SG{spacegroup} source centerings differ")
            if record.operations != baseline.operations:
                raise SourceInvariantError(f"SG{spacegroup} ordered Seitz operations differ")
        pir = by_archive[SourceArchive.PIR]
        cir = by_archive[SourceArchive.CIR]
        if validate_census and (not pir or not cir):
            raise SourceInvariantError(f"SG{spacegroup} lacks a PIR/CIR universe counterpart")
        if pir and cir:
            if pir[0].space_group_symbol != cir[0].space_group_symbol:
                raise SourceInvariantError(f"SG{spacegroup} PIR/CIR symbols differ")
            if pir[0].centering is not cir[0].centering:
                raise SourceInvariantError(f"SG{spacegroup} PIR/CIR centerings differ")
            if pir[0].operations != cir[0].operations:
                raise SourceInvariantError(f"SG{spacegroup} PIR/CIR ordered operations differ")
        universes.append(
            ExactSpaceGroupUniverse(
                spacegroup=spacegroup,
                space_group_symbol=baseline.space_group_symbol,
                centering=baseline.centering,
                operations=baseline.operations,
                pir_irnumbers=tuple(record.irnumber for record in pir),
                cir_irnumbers=tuple(record.irnumber for record in cir),
            )
        )
    return tuple(universes)


def _assemble_database(
    pir_records: Sequence[ExactSourceRecord],
    cir_records: Sequence[ExactSourceRecord],
    *,
    validate_census: bool,
) -> ExactIsoIrrepDatabase:
    pir_records = tuple(pir_records)
    cir_records = tuple(cir_records)
    if any(record.archive is not SourceArchive.PIR for record in pir_records):
        raise SourceInvariantError("pir_records contains a non-PIR source record")
    if any(record.archive is not SourceArchive.CIR for record in cir_records):
        raise SourceInvariantError("cir_records contains a non-CIR source record")
    universes = _build_universes(pir_records, cir_records, validate_census=validate_census)
    if validate_census:
        if len(pir_records) != PIR_RECORD_COUNT or len(cir_records) != CIR_RECORD_COUNT:
            raise SourceInvariantError("source record census mismatch")
        if sum(universe is not None for universe in universes) != 230:
            raise SourceInvariantError("source universe census is not exactly 230")
        if sum(len(universe.operations) for universe in universes[1:] if universe is not None) != 2_609:
            raise SourceInvariantError("source representative operation census is not 2609")
        operation_denominators_pir, translation_denominators_pir = _record_denominators(pir_records)
        operation_denominators_cir, translation_denominators_cir = _record_denominators(cir_records)
        if operation_denominators_pir != EXPECTED_DENOMINATORS or operation_denominators_cir != EXPECTED_DENOMINATORS:
            raise SourceInvariantError("operation denominator census mismatch")
        if translation_denominators_pir != EXPECTED_DENOMINATORS or translation_denominators_cir != EXPECTED_DENOMINATORS:
            raise SourceInvariantError("irtranslation denominator census mismatch")
        pir_translation_count = sum(
            translation is not None
            for record in pir_records
            for translation in record.irtranslations
        )
        cir_translation_count = sum(
            translation is not None
            for record in cir_records
            for translation in record.irtranslations
        )
        if pir_translation_count != PIR_IRTRANSLATION_ROW_COUNT:
            raise SourceInvariantError("PIR irtranslation row census mismatch")
        if cir_translation_count != CIR_IRTRANSLATION_ROW_COUNT:
            raise SourceInvariantError("CIR irtranslation row census mismatch")
        special_count = sum(record.special for record in pir_records + cir_records)
        if special_count != 10_073 or len(pir_records) + len(cir_records) - special_count != 11_423:
            raise SourceInvariantError("special/parameterized source census mismatch")
        centering_counts = Counter(universe.centering.value for universe in universes[1:] if universe is not None)
        observed_centering_counts = {
            centering: centering_counts.get(centering, 0)
            for centering in EXPECTED_CENTERING_COUNTS
        }
        if observed_centering_counts != EXPECTED_CENTERING_COUNTS:
            raise SourceInvariantError(
                f"centering census mismatch: got {observed_centering_counts!r}"
            )
    return ExactIsoIrrepDatabase(
        pir_records=pir_records,
        cir_records=cir_records,
        universes=universes,
    )


def _load_uncached() -> ExactIsoIrrepDatabase:
    pir_text = _read_verified_member(SourceArchive.PIR)
    pir_records = _parse_records_from_lines(
        _source_lines(pir_text, "PIR source text"),
        SourceArchive.PIR,
        validate_census=True,
    )
    # Release the decoded member before reading the larger CIR source.
    del pir_text
    cir_text = _read_verified_member(SourceArchive.CIR)
    cir_records = _parse_records_from_lines(
        _source_lines(cir_text, "CIR source text"),
        SourceArchive.CIR,
        validate_census=True,
    )
    del cir_text
    return _assemble_database(pir_records, cir_records, validate_census=True)


_DATABASE_LOCK = threading.Lock()
_DATABASE: Optional[ExactIsoIrrepDatabase] = None


def load_exact_iso_irrep_sources() -> ExactIsoIrrepDatabase:
    """Load the verified PIR/CIR source frames once and return one identity."""

    global _DATABASE
    cached = _DATABASE
    if cached is not None:
        return cached
    with _DATABASE_LOCK:
        cached = _DATABASE
        if cached is None:
            cached = _load_uncached()
            _DATABASE = cached
        return cached


def _lookup_universe(
    universes: tuple[Optional[ExactSpaceGroupUniverse], ...],
    spacegroup: int,
) -> ExactSpaceGroupUniverse:
    _validate_spacegroup(spacegroup)
    if type(universes) is not tuple or len(universes) != 231:
        raise SourceLookupError("source universes must be an immutable 231-entry tuple")
    universe = universes[spacegroup]
    if universe is None:
        raise SourceLookupError(f"no source universe for SG{spacegroup}")
    return universe


def _validate_spacegroup(spacegroup: int) -> None:
    if type(spacegroup) is not int:
        raise SourceLookupError("spacegroup lookup requires an exact int")
    if not 1 <= spacegroup <= 230:
        raise SourceLookupError(f"spacegroup lookup {spacegroup!r} outside 1..230")


def source_universe(spacegroup: int) -> ExactSpaceGroupUniverse:
    """Look up one strict source universe by SG number."""

    _validate_spacegroup(spacegroup)
    return _lookup_universe(load_exact_iso_irrep_sources().universes, spacegroup)


__all__ = [
    "ArchiveIntegrityError",
    "Centering",
    "CIR_ARCHIVE_BYTES",
    "CIR_ARCHIVE_PATH",
    "CIR_ARCHIVE_SHA256",
    "ExactIrTranslation",
    "ExactIsoIrrepDatabase",
    "ExactKArm",
    "ExactSeitz",
    "ExactSourceRecord",
    "ExactSpaceGroupUniverse",
    "IsoIrrepExactError",
    "PIR_ARCHIVE_BYTES",
    "PIR_ARCHIVE_PATH",
    "PIR_ARCHIVE_SHA256",
    "SourceArchive",
    "SourceInvariantError",
    "SourceLookupError",
    "SourceSchemaError",
    "load_exact_iso_irrep_sources",
    "parse_exact_source_lines",
    "parse_exact_source_text",
    "source_universe",
]
