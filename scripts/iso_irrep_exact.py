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
from pathlib import Path
import re
import threading
from typing import Iterable, Optional, Sequence, Tuple
import zipfile

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


@dataclass(frozen=True)
class ExactSeitz:
    __slots__ = ("rotation", "translation", "raw_augmented")

    rotation: Rotation3
    translation: Fraction3
    raw_augmented: RawAugmented


@dataclass(frozen=True)
class ExactKArm:
    __slots__ = ("constant", "parameters", "raw_augmented")

    constant: Fraction3
    parameters: tuple[OptionalFraction3, OptionalFraction3, OptionalFraction3]
    raw_augmented: RawKVector


@dataclass(frozen=True)
class ExactIrTranslation:
    __slots__ = ("vector", "raw")

    vector: Fraction3
    raw: RawIrTranslation


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


@dataclass(frozen=True)
class ExactIsoIrrepDatabase:
    __slots__ = ("pir_records", "cir_records", "universes")

    pir_records: tuple[ExactSourceRecord, ...]
    cir_records: tuple[ExactSourceRecord, ...]
    universes: tuple[Optional[ExactSpaceGroupUniverse], ...]

    def source_universe(self, spacegroup: int) -> ExactSpaceGroupUniverse:
        return _lookup_universe(self.universes, spacegroup)


_SIGNED_INTEGER_RE = re.compile(r"(?:0|-[1-9][0-9]*|[1-9][0-9]*)\Z", re.ASCII)
_UNSIGNED_INTEGER_RE = re.compile(r"(?:0|[1-9][0-9]*)\Z", re.ASCII)
_CIR_COMPLEX_TOKEN_RE = re.compile(r"\(([^,]+),([^\)]+)\)\Z", re.ASCII)


def _error(error_type: type[IsoIrrepExactError], message: str) -> IsoIrrepExactError:
    return error_type(message)


def _parse_integer(token: str, *, context: str, unsigned: bool = False) -> int:
    """Parse one canonical ASCII integer and nothing accepted by ``int`` more."""

    pattern = _UNSIGNED_INTEGER_RE if unsigned else _SIGNED_INTEGER_RE
    if pattern.fullmatch(token) is None:
        raise SourceSchemaError(f"non-canonical integer {token!r} for {context}")
    # The regular expression has already excluded Unicode digits and signs
    # that Python's int() would otherwise accept.
    return int(token)


def _is_ascii_space(char: str) -> bool:
    return char in " \t"


def _parse_header(line: str, archive: SourceArchive, line_number: int) -> tuple[int, int, str, str, int, int, int, int, int]:
    """Parse the exact nine-field quoted header without a permissive regex."""

    position = 0
    length = len(line)

    def skip_space() -> None:
        nonlocal position
        while position < length and _is_ascii_space(line[position]):
            position += 1

    def read_unsigned(field: str) -> int:
        nonlocal position
        start = position
        while position < length and not _is_ascii_space(line[position]):
            position += 1
        if start == position:
            raise SourceSchemaError(
                f"missing {archive.value} header field {field} at line {line_number}"
            )
        return _parse_integer(line[start:position], context=f"{archive.value} header {field}", unsigned=True)

    def read_quoted(field: str) -> str:
        nonlocal position
        if position >= length or line[position] != '"':
            raise SourceSchemaError(
                f"missing quoted {archive.value} header field {field} at line {line_number}"
            )
        position += 1
        start = position
        end = line.find('"', position)
        if end < 0:
            raise SourceSchemaError(
                f"unterminated quoted {archive.value} header field {field} at line {line_number}"
            )
        value = line[start:end]
        position = end + 1
        if position < length and not _is_ascii_space(line[position]):
            raise SourceSchemaError(
                f"unexpected character after quoted {field} at line {line_number}"
            )
        return value.strip()

    skip_space()
    irnumber = read_unsigned("irnumber")
    skip_space()
    spacegroup = read_unsigned("spacegroup")
    skip_space()
    raw_symbol = read_quoted("space-group symbol")
    skip_space()
    raw_label = read_quoted("irrep label")
    symbol = raw_symbol.strip()
    label = raw_label.strip()
    values = []
    for field in ("dimension", "irtype", "kcount", "pmkcount", "opcount"):
        skip_space()
        values.append(read_unsigned(field))
    skip_space()
    if position != length:
        raise SourceSchemaError(
            f"extra {archive.value} header fields at line {line_number}"
        )
    dimension, irtype, kcount, pmkcount, opcount = values

    if not raw_symbol:
        raise SourceSchemaError(f"empty space-group symbol at line {line_number}")
    try:
        Centering(raw_symbol[0])
    except ValueError as error:
        raise SourceSchemaError(
            f"unknown centering {raw_symbol[0]!r} at line {line_number}"
        ) from error
    if not label:
        raise SourceSchemaError(f"empty irrep label at line {line_number}")
    if irnumber <= 0:
        raise SourceSchemaError(f"irnumber {irnumber} must be positive at line {line_number}")
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
    # Source rows are ASCII-token records.  split() is only used for the
    # separators; every integer token still passes the canonical ASCII gate.
    return line.split()


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
    return _parse_records_from_lines(text.splitlines(), archive, validate_census=validate_census)


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
    return _parse_records_from_lines(lines, archive, validate_census=validate_census)


def _parse_archive_text(
    text: str,
    archive: SourceArchive,
    *,
    validate_census: bool = False,
) -> tuple[ExactSourceRecord, ...]:
    """Compatibility seam used by focused tests; see parse_exact_source_text."""

    return parse_exact_source_text(text, archive, validate_census=validate_census)


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
    else:
        path = Path(path)
    expected_bytes = (
        PIR_ARCHIVE_BYTES if archive is SourceArchive.PIR else CIR_ARCHIVE_BYTES
    ) if expected_bytes is None else expected_bytes
    expected_sha256 = (
        PIR_ARCHIVE_SHA256 if archive is SourceArchive.PIR else CIR_ARCHIVE_SHA256
    ) if expected_sha256 is None else expected_sha256
    try:
        actual_bytes = path.stat().st_size
        if actual_bytes != expected_bytes:
            raise ArchiveIntegrityError(
                f"{archive.value} archive byte-size mismatch: expected {expected_bytes}, got {actual_bytes}"
            )
        archive_bytes = path.read_bytes()
    except ArchiveIntegrityError:
        raise
    except OSError as error:
        raise ArchiveIntegrityError(f"unable to read {archive.value} archive {path}") from error
    actual_sha256 = hashlib.sha256(archive_bytes).hexdigest()
    if actual_sha256 != expected_sha256:
        raise ArchiveIntegrityError(
            f"{archive.value} archive SHA-256 mismatch: expected {expected_sha256}, got {actual_sha256}"
        )
    member_name = "PIR_data.txt" if archive is SourceArchive.PIR else "CIR_data.txt"
    try:
        with zipfile.ZipFile(path) as zip_file:
            matches = tuple(
                name
                for name in zip_file.namelist()
                if name == member_name or name.endswith("/" + member_name)
            )
            if len(matches) != 1:
                raise ArchiveIntegrityError(
                    f"{archive.value} archive member {member_name!r} matched {len(matches)} entries"
                )
            member_bytes = zip_file.read(matches[0])
    except ArchiveIntegrityError:
        raise
    except (OSError, zipfile.BadZipFile, KeyError) as error:
        raise ArchiveIntegrityError(
            f"unable to read authoritative {member_name} from {archive.value} archive"
        ) from error
    try:
        return member_bytes.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise SourceSchemaError(
            f"authoritative {member_name} member is not strict UTF-8"
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
        pir_text.splitlines(), SourceArchive.PIR, validate_census=True
    )
    # Release the decoded member before reading the larger CIR source.
    del pir_text
    cir_text = _read_verified_member(SourceArchive.CIR)
    cir_records = _parse_records_from_lines(
        cir_text.splitlines(), SourceArchive.CIR, validate_census=True
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
