"""Typed, exact access to the committed spglib magnetic provenance artifact.

The public loader has one fixed repository-relative input.  It verifies the
committed bytes before parsing them, then converts the already validated raw
tables into immutable dataclasses and independently checks their integer
group laws.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from fractions import Fraction
import hashlib
from pathlib import Path
import threading
from typing import Optional, Tuple

try:
    from . import extract_spglib_magnetic_provenance as _extractor
except ImportError:
    import extract_spglib_magnetic_provenance as _extractor


SCHEMA = "cryspglib-spglib-magnetic-v1"
TRANSLATION_DENOMINATOR = 12
ROTATION_RADIX = 3
ROTATION_DIGITS = 9
ROTATION_PAYLOAD = ROTATION_RADIX ** ROTATION_DIGITS
TRANSLATION_DIGITS = 3
TRANSLATION_PAYLOAD = TRANSLATION_DENOMINATOR ** TRANSLATION_DIGITS
SPACE_OPERATION_SCALE = ROTATION_PAYLOAD * TRANSLATION_PAYLOAD
MAGNETIC_OPERATION_ENCODING_LIMIT = 2 * SPACE_OPERATION_SCALE

SPG_HALL_COUNT = 531
SPG_HALL_SETTINGS = 530
SPG_OPERATION_COUNT = 8_147
SPG_STANDARD_OPERATION_END = 7_389
SPG_STANDARD_OPERATION_COUNT = 7_388
SPG_LAYER_OPERATION_COUNT = 758
MSG_UNI_COUNT = 1_652
MSG_HALL_SLOTS = 18
MSG_OPERATION_COUNT = 76_683
MSG_ACTIVE_SPAN_COUNT = 4_479
ALTERNATIVE_TRANSFORMATION_VALUE_COUNT = 536

ARTIFACT_BYTE_LENGTH = 1_537_875
ARTIFACT_SHA256 = (
    "933a52a6696e7f6a1a2e426825ad92c377c6e96330e18c5c045d659798d740b9"
)
MANIFEST_BYTE_LENGTH = 570
MANIFEST_SHA256 = (
    "6a9e1b64c190c30a556d63e51e5b896b967d33e8821714beb745ae699fab84bf"
)
_ARTIFACT_NAME = "spglib_magnetic_provenance_v1.json"
_MANIFEST_NAME = "spglib_magnetic_provenance_v1.manifest.json"

_IDENTITY_ROTATION = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
_IDENTITY_ROTATION_FLAT = (1, 0, 0, 0, 1, 0, 0, 0, 1)
_IDENTITY_KEY = (_IDENTITY_ROTATION_FLAT, (0, 0, 0), 0)
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


class MagneticProvenanceError(ValueError):
    """Base class for all typed loader failures."""


class MagneticProvenanceIntegrityError(MagneticProvenanceError):
    """The fixed bytes or their cryptographic commitments are invalid."""


def _resolve_data_dir():
    try:
        module_path = Path(__file__).resolve()
    except (OSError, RuntimeError) as error:
        raise MagneticProvenanceIntegrityError(
            "unable to resolve magnetic provenance module path"
        ) from error
    return module_path.parent / "data"


# Resolve this while importing the module.  The public loader must not depend
# on the caller's later working directory or on a mutable relative __file__.
_DATA_DIR = _resolve_data_dir()


class MagneticProvenanceSchemaError(MagneticProvenanceError):
    """The parsed artifact cannot satisfy the committed schema."""


class MagneticProvenanceDecodeError(MagneticProvenanceError):
    """An operation encoding cannot be decoded and re-encoded exactly."""


class MagneticProvenanceInvariantError(MagneticProvenanceError):
    """A cross-table, algebraic, or closure invariant failed."""


class MagneticProvenanceLookupError(MagneticProvenanceError):
    """A query used an invalid one-based index or inactive Hall setting."""


class ArtifactLookupError(MagneticProvenanceLookupError):
    """Public alias for structured invalid-query failures."""


class MagneticKind(IntEnum):
    TYPE_I = 1
    GREY = 2
    TYPE_II = 2
    TYPE_III = 3
    ANTI_TRANSLATION = 4
    TYPE_IV = 4


class TimeReversal(IntEnum):
    UNITARY = 0
    NONE = 0
    ANTIUNITARY = 1
    PRIMED = 1


@dataclass(**_DATACLASS_OPTIONS)
class OperationSpan:
    if not _NATIVE_DATACLASS_SLOTS:
        __slots__ = ("order", "offset")
    order: int
    offset: int


@dataclass(**_DATACLASS_OPTIONS)
class ExactSeitzOperation:
    if not _NATIVE_DATACLASS_SLOTS:
        __slots__ = (
            "encoded", "rotation", "translation_numerator", "time_reversal"
        )
    encoded: int
    rotation: Tuple[
        Tuple[int, int, int], Tuple[int, int, int], Tuple[int, int, int]
    ]
    translation_numerator: Tuple[int, int, int]
    time_reversal: TimeReversal

    @property
    def translation(self) -> Tuple[Fraction, Fraction, Fraction]:
        return tuple(
            Fraction(value, TRANSLATION_DENOMINATOR)
            for value in self.translation_numerator
        )


@dataclass(**_DATACLASS_OPTIONS)
class MagneticGroupMetadata:
    if not _NATIVE_DATACLASS_SLOTS:
        __slots__ = (
            "uni", "litvin", "bns", "og", "parent_spacegroup", "kind"
        )
    uni: int
    litvin: int
    bns: str
    og: str
    parent_spacegroup: int
    kind: MagneticKind


@dataclass(**_DATACLASS_OPTIONS)
class SpgProvenance:
    if not _NATIVE_DATACLASS_SLOTS:
        __slots__ = (
            "spacegroup_numbers", "operation_index", "raw_operation_codes",
            "decoded_operations"
        )
    spacegroup_numbers: Tuple[int, ...]
    operation_index: Tuple[OperationSpan, ...]
    raw_operation_codes: Tuple[int, ...]
    decoded_operations: Tuple[Optional[ExactSeitzOperation], ...]

    def _hall(self, hall: int) -> int:
        if type(hall) is not int or not 1 <= hall < SPG_HALL_COUNT:
            raise ArtifactLookupError("Hall number must be in 1..530")
        return hall

    def spacegroup_number_for_hall(self, hall: int) -> int:
        return self.spacegroup_numbers[self._hall(hall)]

    def spg_operation_span(self, hall: int) -> OperationSpan:
        return self.operation_index[self._hall(hall)]

    def spg_operations(self, hall: int) -> Tuple[ExactSeitzOperation, ...]:
        span = self.spg_operation_span(hall)
        operations = self.decoded_operations[span.offset:span.offset + span.order]
        if any(operation is None for operation in operations):
            raise MagneticProvenanceInvariantError("SPG query crossed sentinel")
        return operations  # type: ignore[return-value]


@dataclass(**_DATACLASS_OPTIONS)
class MsgProvenance:
    if not _NATIVE_DATACLASS_SLOTS:
        __slots__ = (
            "metadata", "uni_mapping", "derived_hall_mapping", "operation_index",
            "raw_operation_codes", "decoded_operations", "alternative_codes"
        )
    metadata: Tuple[Optional[MagneticGroupMetadata], ...]
    uni_mapping: Tuple[Tuple[int, int], ...]
    derived_hall_mapping: Tuple[Tuple[int, int], ...]
    operation_index: Tuple[Tuple[OperationSpan, ...], ...]
    raw_operation_codes: Tuple[int, ...]
    decoded_operations: Tuple[Optional[ExactSeitzOperation], ...]
    alternative_codes: Tuple[Tuple[Tuple[int, ...], ...], ...]

    def _uni(self, uni: int) -> int:
        if type(uni) is not int or not 1 <= uni < MSG_UNI_COUNT:
            raise ArtifactLookupError("UNI number must be in 1..1651")
        return uni

    def _hall(self, hall: int) -> int:
        if type(hall) is not int or not 1 <= hall < SPG_HALL_COUNT:
            raise ArtifactLookupError("Hall number must be in 1..530")
        return hall

    def magnetic_metadata(self, uni: int) -> MagneticGroupMetadata:
        value = self.metadata[self._uni(uni)]
        if value is None:
            raise MagneticProvenanceInvariantError("metadata sentinel in active UNI")
        return value

    def halls_for_uni(self, uni: int) -> Tuple[int, ...]:
        count, first = self.uni_mapping[self._uni(uni)]
        return tuple(range(first, first + count))

    def unis_for_hall(self, hall: int) -> Tuple[int, ...]:
        smallest, largest = self.derived_hall_mapping[self._hall(hall)]
        return tuple(range(smallest, largest + 1))

    def _slot(self, uni: int, hall: int) -> Tuple[int, int]:
        uni = self._uni(uni)
        hall = self._hall(hall)
        count, first = self.uni_mapping[uni]
        slot = hall - first
        if not 0 <= slot < count:
            raise ArtifactLookupError("Hall is inactive for this UNI")
        return uni, slot

    def magnetic_operation_span(self, uni: int, hall: int) -> OperationSpan:
        uni, slot = self._slot(uni, hall)
        return self.operation_index[uni][slot]

    def magnetic_operations(
        self, uni: int, hall: int
    ) -> Tuple[ExactSeitzOperation, ...]:
        span = self.magnetic_operation_span(uni, hall)
        operations = self.decoded_operations[span.offset:span.offset + span.order]
        if any(operation is None for operation in operations):
            raise MagneticProvenanceInvariantError("MSG query crossed sentinel")
        return operations  # type: ignore[return-value]

    def raw_alternative_codes(self, uni: int, hall: int) -> Tuple[int, ...]:
        uni, slot = self._slot(uni, hall)
        return self.alternative_codes[uni][slot]

    def std_transformations(
        self, uni: int, hall: int
    ) -> Tuple[ExactSeitzOperation, ...]:
        codes = self.raw_alternative_codes(uni, hall)
        return (_IDENTITY_OPERATION,) + tuple(_decode_operation(code) for code in codes)


@dataclass(**_DATACLASS_OPTIONS)
class MagneticProvenanceDatabase:
    if not _NATIVE_DATACLASS_SLOTS:
        __slots__ = ("spg", "msg")
    spg: SpgProvenance
    msg: MsgProvenance

    def spacegroup_number_for_hall(self, hall: int) -> int:
        return self.spg.spacegroup_number_for_hall(hall)

    def spg_operation_span(self, hall: int) -> OperationSpan:
        return self.spg.spg_operation_span(hall)

    def spg_operations(self, hall: int) -> Tuple[ExactSeitzOperation, ...]:
        return self.spg.spg_operations(hall)

    def magnetic_metadata(self, uni: int) -> MagneticGroupMetadata:
        return self.msg.magnetic_metadata(uni)

    def halls_for_uni(self, uni: int) -> Tuple[int, ...]:
        return self.msg.halls_for_uni(uni)

    def unis_for_hall(self, hall: int) -> Tuple[int, ...]:
        return self.msg.unis_for_hall(hall)

    def magnetic_operation_span(self, uni: int, hall: int) -> OperationSpan:
        return self.msg.magnetic_operation_span(uni, hall)

    def magnetic_operations(
        self, uni: int, hall: int
    ) -> Tuple[ExactSeitzOperation, ...]:
        return self.msg.magnetic_operations(uni, hall)

    def raw_alternative_codes(self, uni: int, hall: int) -> Tuple[int, ...]:
        return self.msg.raw_alternative_codes(uni, hall)

    def std_transformations(
        self, uni: int, hall: int
    ) -> Tuple[ExactSeitzOperation, ...]:
        return self.msg.std_transformations(uni, hall)


_CACHE_LOCK = threading.Lock()
_CACHED_DATABASE = None  # type: Optional[MagneticProvenanceDatabase]


def _decode_operation(encoded: int) -> ExactSeitzOperation:
    if (type(encoded) is not int
            or not 0 < encoded < MAGNETIC_OPERATION_ENCODING_LIMIT):
        raise MagneticProvenanceDecodeError("operation encoding is out of range")
    time_value = encoded // SPACE_OPERATION_SCALE
    spatial = encoded % SPACE_OPERATION_SCALE
    if time_value not in (0, 1):
        raise MagneticProvenanceDecodeError("operation time-reversal bit is invalid")
    rotation_payload = spatial % ROTATION_PAYLOAD
    translation_payload = spatial // ROTATION_PAYLOAD
    rotation_digits = tuple(
        (rotation_payload // (ROTATION_RADIX ** exponent)) % ROTATION_RADIX
        for exponent in range(ROTATION_DIGITS - 1, -1, -1)
    )
    translation = tuple(
        (translation_payload // (TRANSLATION_DENOMINATOR ** exponent))
        % TRANSLATION_DENOMINATOR
        for exponent in range(TRANSLATION_DIGITS - 1, -1, -1)
    )
    rotation = tuple(
        tuple(rotation_digits[3 * row + column] - 1 for column in range(3))
        for row in range(3)
    )
    reconstructed_rotation = sum(
        (value + 1) * ROTATION_RADIX ** (ROTATION_DIGITS - 1 - index)
        for index, value in enumerate(
            value for row in rotation for value in row
        )
    )
    reconstructed_translation = sum(
        value * TRANSLATION_DENOMINATOR ** (TRANSLATION_DIGITS - 1 - index)
        for index, value in enumerate(translation)
    )
    reconstructed = (
        time_value * SPACE_OPERATION_SCALE
        + reconstructed_rotation
        + ROTATION_PAYLOAD * reconstructed_translation
    )
    if reconstructed != encoded:
        raise MagneticProvenanceDecodeError("operation round-trip mismatch")
    try:
        time_reversal = TimeReversal(time_value)
    except ValueError as error:
        raise MagneticProvenanceDecodeError("operation time-reversal value invalid") from error
    return ExactSeitzOperation(
        encoded=encoded,
        rotation=rotation,
        translation_numerator=translation,
        time_reversal=time_reversal,
    )


def _encode_key(rotation, translation, time_reversal):
    rotation_payload = sum(
        (value + 1) * ROTATION_RADIX ** (ROTATION_DIGITS - 1 - index)
        for index, value in enumerate(
            value for row in rotation for value in row
        )
    )
    translation_payload = sum(
        value * TRANSLATION_DENOMINATOR ** (TRANSLATION_DIGITS - 1 - index)
        for index, value in enumerate(translation)
    )
    return (
        int(time_reversal) * SPACE_OPERATION_SCALE
        + rotation_payload
        + ROTATION_PAYLOAD * translation_payload
    )


def _operation_key(operation: ExactSeitzOperation):
    return (
        tuple(value for row in operation.rotation for value in row),
        operation.translation_numerator,
        int(operation.time_reversal),
    )


def _determinant(rotation):
    a, b, c, d, e, f, g, h, i = rotation
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def _matrix_multiply(left, right):
    return (
        left[0] * right[0] + left[1] * right[3] + left[2] * right[6],
        left[0] * right[1] + left[1] * right[4] + left[2] * right[7],
        left[0] * right[2] + left[1] * right[5] + left[2] * right[8],
        left[3] * right[0] + left[4] * right[3] + left[5] * right[6],
        left[3] * right[1] + left[4] * right[4] + left[5] * right[7],
        left[3] * right[2] + left[4] * right[5] + left[5] * right[8],
        left[6] * right[0] + left[7] * right[3] + left[8] * right[6],
        left[6] * right[1] + left[7] * right[4] + left[8] * right[7],
        left[6] * right[2] + left[7] * right[5] + left[8] * right[8],
    )


def _matrix_inverse(rotation):
    a, b, c, d, e, f, g, h, i = rotation
    determinant = _determinant(rotation)
    if determinant not in (-1, 1):
        raise MagneticProvenanceInvariantError("rotation determinant is not ±1")
    adjugate = (
        e * i - f * h, c * h - b * i, b * f - c * e,
        f * g - d * i, a * i - c * g, c * d - a * f,
        d * h - e * g, b * g - a * h, a * e - b * d,
    )
    return tuple(value // determinant for value in adjugate)


def _compose_keys(left, right):
    rotation_left, translation_left, time_left = left
    rotation_right, translation_right, time_right = right
    rotation = _matrix_multiply(rotation_left, rotation_right)
    translation = tuple(
        (translation_left[row]
         + sum(rotation_left[3 * row + k] * translation_right[k]
               for k in range(3))) % TRANSLATION_DENOMINATOR
        for row in range(3)
    )
    return rotation, translation, time_left ^ time_right


def _inverse_key(operation_key):
    rotation, translation, time_reversal = operation_key
    inverse_rotation = _matrix_inverse(rotation)
    inverse_translation = tuple(
        (-sum(inverse_rotation[3 * row + k] * translation[k]
              for k in range(3))) % TRANSLATION_DENOMINATOR
        for row in range(3)
    )
    return inverse_rotation, inverse_translation, time_reversal


_IDENTITY_OPERATION = ExactSeitzOperation(
    encoded=16484,
    rotation=_IDENTITY_ROTATION,
    translation_numerator=(0, 0, 0),
    time_reversal=TimeReversal.UNITARY,
)


__all__ = [
    "ArtifactLookupError", "ExactSeitzOperation", "MagneticGroupMetadata",
    "MagneticKind", "MagneticProvenanceDatabase", "MagneticProvenanceDecodeError",
    "MagneticProvenanceError", "MagneticProvenanceIntegrityError",
    "MagneticProvenanceInvariantError", "MagneticProvenanceLookupError",
    "MagneticProvenanceSchemaError", "MsgProvenance", "OperationSpan",
    "SpgProvenance", "TimeReversal", "load_committed_provenance",
]


def _schema_value(value, label):
    if type(value) is not int:
        raise MagneticProvenanceSchemaError(f"{label} must be an integer")
    return value


def _schema_pair(value, label):
    if type(value) is not list or len(value) != 2:
        raise MagneticProvenanceSchemaError(f"{label} must be an integer pair")
    return (
        _schema_value(value[0], f"{label}[0]"),
        _schema_value(value[1], f"{label}[1]"),
    )


def _schema_span(value, label):
    order, offset = _schema_pair(value, label)
    return OperationSpan(order=order, offset=offset)


def _validate_operation(operation, label, allow_antiunitary):
    flat_rotation = tuple(value for row in operation.rotation for value in row)
    if _determinant(flat_rotation) not in (-1, 1):
        raise MagneticProvenanceInvariantError(
            f"{label} rotation determinant is not ±1"
        )
    if any(type(value) is not int or not 0 <= value < TRANSLATION_DENOMINATOR
           for value in operation.translation_numerator):
        raise MagneticProvenanceInvariantError(
            f"{label} translation numerator is out of range"
        )
    if not allow_antiunitary and operation.time_reversal is not TimeReversal.UNITARY:
        raise MagneticProvenanceInvariantError(
            f"{label} unexpectedly has time reversal"
        )
    if _encode_key(
        operation.rotation,
        operation.translation_numerator,
        operation.time_reversal,
    ) != operation.encoded:
        raise MagneticProvenanceDecodeError(f"{label} does not re-encode exactly")


def _validate_group(operations, label, allow_antiunitary):
    if not operations:
        raise MagneticProvenanceInvariantError(f"{label} is empty")
    for index, operation in enumerate(operations):
        _validate_operation(operation, f"{label}[{index}]", allow_antiunitary)
    operation_keys = tuple(_operation_key(operation) for operation in operations)
    keys = set(operation_keys)
    if len(keys) != len(operation_keys):
        raise MagneticProvenanceInvariantError(f"{label} contains duplicate operations")
    identity_count = sum(key == _IDENTITY_KEY for key in operation_keys)
    if identity_count != 1:
        raise MagneticProvenanceInvariantError(
            f"{label} does not contain exactly one {{I|0}}"
        )
    for left_key in operation_keys:
        for right_key in operation_keys:
            if _compose_keys(left_key, right_key) not in keys:
                raise MagneticProvenanceInvariantError(f"{label} is not closed")
    for operation_key in operation_keys:
        if _inverse_key(operation_key) not in keys:
            raise MagneticProvenanceInvariantError(f"{label} lacks an inverse")


def _convert_spg(data):
    try:
        spg = data["spg"]
        numbers_data = spg["spacegroup_number"]
        index_data = spg["symmetry_operation_index"]
        raw_data = spg["symmetry_operations"]
    except (KeyError, TypeError) as error:
        raise MagneticProvenanceSchemaError("SPG artifact fields are missing") from error
    if type(numbers_data) is not list or len(numbers_data) != SPG_HALL_COUNT:
        raise MagneticProvenanceSchemaError("SPG spacegroup census mismatch")
    if type(index_data) is not list or len(index_data) != SPG_HALL_COUNT:
        raise MagneticProvenanceSchemaError("SPG operation-index census mismatch")
    if type(raw_data) is not list or len(raw_data) != SPG_OPERATION_COUNT:
        raise MagneticProvenanceSchemaError("SPG operation census mismatch")

    numbers = tuple(_schema_value(value, f"SPG spacegroup_number[{index}]")
                    for index, value in enumerate(numbers_data))
    operation_index = tuple(
        _schema_span(value, f"SPG operation_index[{index}]")
        for index, value in enumerate(index_data)
    )
    raw_codes = tuple(
        _schema_value(value, f"SPG operation[{index}]")
        for index, value in enumerate(raw_data)
    )
    if numbers[0] != 0 or any(not 1 <= value <= 230 for value in numbers[1:]):
        raise MagneticProvenanceInvariantError("SPG spacegroup number range mismatch")
    if operation_index[0] != OperationSpan(0, 0) or raw_codes[0] != 0:
        raise MagneticProvenanceInvariantError("SPG sentinel mismatch")
    if any(not 0 < value < SPACE_OPERATION_SCALE for value in raw_codes[1:]):
        raise MagneticProvenanceInvariantError("SPG operation encoding range mismatch")

    decoded = [None]
    for index, raw_code in enumerate(raw_codes[1:], 1):
        operation = _decode_operation(raw_code)
        _validate_operation(operation, f"SPG operation[{index}]", False)
        decoded.append(operation)
    decoded_operations = tuple(decoded)

    previous_end = 1
    for hall, span in enumerate(operation_index[1:], 1):
        if (span.order <= 0 or span.offset != previous_end
                or span.offset + span.order > SPG_STANDARD_OPERATION_END):
            raise MagneticProvenanceInvariantError(
                f"SPG operation span {hall} is not in the standard layer"
            )
        previous_end = span.offset + span.order
    if previous_end != SPG_STANDARD_OPERATION_END:
        raise MagneticProvenanceInvariantError("SPG standard span boundary mismatch")
    if len(raw_codes) - SPG_STANDARD_OPERATION_END != SPG_LAYER_OPERATION_COUNT:
        raise MagneticProvenanceInvariantError("SPG layer-tail census mismatch")
    for hall, span in enumerate(operation_index[1:], 1):
        group = decoded_operations[span.offset:span.offset + span.order]
        _validate_group(group, f"SPG Hall {hall}", False)

    result = SpgProvenance(
        spacegroup_numbers=numbers,
        operation_index=operation_index,
        raw_operation_codes=raw_codes,
        decoded_operations=decoded_operations,
    )
    return result


def _convert_metadata(data):
    if type(data) is not list or len(data) != MSG_UNI_COUNT:
        raise MagneticProvenanceSchemaError("MSG metadata census mismatch")
    sentinel = data[0]
    if (type(sentinel) is not dict
            or sentinel != {
                "uni": 0, "litvin": 0, "bns": "", "og": "",
                "parent_spacegroup": 0, "type": 0,
            }):
        raise MagneticProvenanceInvariantError("MSG metadata sentinel mismatch")
    metadata = [None]
    for uni, row in enumerate(data[1:], 1):
        if type(row) is not dict:
            raise MagneticProvenanceSchemaError(f"MSG metadata {uni} is not an object")
        required = {"uni", "litvin", "bns", "og", "parent_spacegroup", "type"}
        if set(row) != required:
            raise MagneticProvenanceSchemaError(f"MSG metadata {uni} keys mismatch")
        if type(row["bns"]) is not str or type(row["og"]) is not str:
            raise MagneticProvenanceSchemaError(f"MSG metadata {uni} labels are not text")
        row_uni = _schema_value(row["uni"], f"MSG metadata {uni}.uni")
        if row_uni != uni:
            raise MagneticProvenanceInvariantError(f"MSG metadata {uni} identity mismatch")
        litvin = _schema_value(row["litvin"], f"MSG metadata {uni}.litvin")
        parent = _schema_value(
            row["parent_spacegroup"], f"MSG metadata {uni}.parent_spacegroup"
        )
        if not 1 <= litvin < MSG_UNI_COUNT:
            raise MagneticProvenanceInvariantError(
                f"MSG metadata {uni} Litvin number is out of range"
            )
        if not 1 <= parent <= 230:
            raise MagneticProvenanceInvariantError(
                f"MSG metadata {uni} parent spacegroup is out of range"
            )
        kind_value = _schema_value(row["type"], f"MSG metadata {uni}.type")
        try:
            kind = MagneticKind(kind_value)
        except ValueError as error:
            raise MagneticProvenanceSchemaError(
                f"MSG metadata {uni} kind is not 1..4"
            ) from error
        metadata.append(MagneticGroupMetadata(
            uni=row_uni,
            litvin=litvin,
            bns=row["bns"],
            og=row["og"],
            parent_spacegroup=parent,
            kind=kind,
        ))
    return tuple(metadata)


def _convert_mapping(data):
    if type(data) is not list or len(data) != MSG_UNI_COUNT:
        raise MagneticProvenanceSchemaError("MSG UNI mapping census mismatch")
    mapping = tuple(_schema_pair(value, f"MSG UNI mapping[{index}]")
                    for index, value in enumerate(data))
    if mapping[0] != (0, 0):
        raise MagneticProvenanceInvariantError("MSG UNI mapping sentinel mismatch")
    for uni, (count, first) in enumerate(mapping[1:], 1):
        if not 1 <= count <= MSG_HALL_SLOTS:
            raise MagneticProvenanceInvariantError(f"UNI {uni} Hall count is invalid")
        if not 1 <= first <= SPG_HALL_SETTINGS or first + count - 1 > SPG_HALL_SETTINGS:
            raise MagneticProvenanceInvariantError(f"UNI {uni} Hall range is invalid")
    return mapping


def _derive_hall_mapping(mapping, spacegroup_numbers, metadata):
    if len(spacegroup_numbers) != SPG_HALL_COUNT:
        raise MagneticProvenanceInvariantError("SPG Hall census is not 531")
    derived = [(0, 0)]
    for hall in range(1, SPG_HALL_COUNT):
        unis = [
            uni for uni, (count, first) in enumerate(mapping[1:], 1)
            if first <= hall < first + count
        ]
        if not unis or unis != list(range(unis[0], unis[-1] + 1)):
            raise MagneticProvenanceInvariantError(
                f"Hall {hall} UNI inverse range is not continuous"
            )
        for uni in unis:
            group_metadata = metadata[uni]
            if (group_metadata is None
                    or spacegroup_numbers[hall] != group_metadata.parent_spacegroup):
                raise MagneticProvenanceInvariantError(
                    f"Hall {hall} parent spacegroup disagrees with UNI {uni}"
                )
        derived.append((unis[0], unis[-1]))
    return tuple(derived)


def _convert_operation_index(data, rows, label):
    if type(data) is not list or len(data) != rows:
        raise MagneticProvenanceSchemaError(f"{label} census mismatch")
    result = []
    for row, row_data in enumerate(data):
        if type(row_data) is not list or len(row_data) != MSG_HALL_SLOTS:
            raise MagneticProvenanceSchemaError(
                f"{label}[{row}] must have 18 entries"
            )
        result.append(tuple(
            _schema_span(value, f"{label}[{row}][{slot}]")
            for slot, value in enumerate(row_data)
        ))
    return tuple(result)


def _convert_raw_codes(data, count, label):
    if type(data) is not list or len(data) != count:
        raise MagneticProvenanceSchemaError(f"{label} census mismatch")
    return tuple(_schema_value(value, f"{label}[{index}]")
                 for index, value in enumerate(data))


def _convert_alternatives(data):
    if type(data) is not list or len(data) != MSG_UNI_COUNT:
        raise MagneticProvenanceSchemaError("alternative transformation census mismatch")
    result = []
    for uni, row in enumerate(data):
        if type(row) is not list or len(row) != MSG_HALL_SLOTS:
            raise MagneticProvenanceSchemaError(
                f"alternative_transformations[{uni}] must have 18 entries"
            )
        converted_row = []
        for slot, values in enumerate(row):
            if type(values) is not list or len(values) != 7:
                raise MagneticProvenanceSchemaError(
                    f"alternative_transformations[{uni}][{slot}] must have 7 entries"
                )
            converted_row.append(tuple(
                _schema_value(value, f"alternative_transformations[{uni}][{slot}][{index}]")
                for index, value in enumerate(values)
            ))
        result.append(tuple(converted_row))
    return tuple(result)


def _validate_alternatives(alternative_rows, mapping):
    value_count = 0
    result = []
    for uni, row in enumerate(alternative_rows):
        converted_row = []
        count = 0 if uni == 0 else mapping[uni][0]
        for slot, values in enumerate(row):
            first_zero = next(
                (index for index, value in enumerate(values) if value == 0),
                len(values),
            )
            if slot >= count:
                if any(value != 0 for value in values):
                    raise MagneticProvenanceInvariantError(
                        f"inactive alternative slot {uni}/{slot} is nonzero"
                    )
                converted_row.append(())
                continue
            if first_zero == len(values):
                raise MagneticProvenanceInvariantError(
                    f"alternative slot {uni}/{slot} lacks a terminator"
                )
            if any(value != 0 for value in values[first_zero:]):
                raise MagneticProvenanceInvariantError(
                    f"alternative slot {uni}/{slot} has a nonzero tail"
                )
            codes = values[:first_zero]
            if any(not 0 < value < SPACE_OPERATION_SCALE for value in codes):
                raise MagneticProvenanceInvariantError(
                    f"alternative slot {uni}/{slot} has an invalid encoding"
                )
            for index, code in enumerate(codes):
                operation = _decode_operation(code)
                _validate_operation(
                    operation, f"alternative operation {uni}/{slot}/{index}", False
                )
            value_count += len(codes)
            converted_row.append(tuple(codes))
        result.append(tuple(converted_row))
    if value_count != ALTERNATIVE_TRANSFORMATION_VALUE_COUNT:
        raise MagneticProvenanceInvariantError("alternative transformation census mismatch")
    return tuple(result)


def _convert_msg(data, spg):
    try:
        msg = data["msg"]
        metadata_data = msg["magnetic_spacegroup_types"]
        mapping_data = msg["magnetic_spacegroup_uni_mapping"]
        index_data = msg["magnetic_spacegroup_operation_index"]
        raw_data = msg["magnetic_symmetry_operations"]
        alternative_data = msg["alternative_transformations"]
    except (KeyError, TypeError) as error:
        raise MagneticProvenanceSchemaError("MSG artifact fields are missing") from error

    metadata = _convert_metadata(metadata_data)
    mapping = _convert_mapping(mapping_data)
    operation_index = _convert_operation_index(
        index_data, MSG_UNI_COUNT, "MSG operation_index"
    )
    raw_codes = _convert_raw_codes(raw_data, MSG_OPERATION_COUNT,
                                   "MSG operation")
    alternative_rows = _validate_alternatives(
        _convert_alternatives(alternative_data), mapping
    )
    derived_hall_mapping = _derive_hall_mapping(
        mapping, spg.spacegroup_numbers, metadata
    )
    if metadata[0] is not None:
        raise MagneticProvenanceInvariantError("MSG metadata sentinel is not empty")
    if raw_codes[0] != 0:
        raise MagneticProvenanceInvariantError("MSG operation sentinel mismatch")
    if any(not 0 < value < MAGNETIC_OPERATION_ENCODING_LIMIT
           for value in raw_codes[1:]):
        raise MagneticProvenanceInvariantError("MSG operation encoding range mismatch")

    decoded = [None]
    for index, raw_code in enumerate(raw_codes[1:], 1):
        operation = _decode_operation(raw_code)
        _validate_operation(operation, f"MSG operation[{index}]", True)
        decoded.append(operation)
    decoded_operations = tuple(decoded)

    if operation_index[0] != tuple(OperationSpan(0, 0) for _ in range(MSG_HALL_SLOTS)):
        raise MagneticProvenanceInvariantError("MSG operation-index sentinel mismatch")
    if alternative_rows[0] != tuple(() for _ in range(MSG_HALL_SLOTS)):
        raise MagneticProvenanceInvariantError("MSG alternative sentinel mismatch")

    active_spans = []
    for uni in range(1, MSG_UNI_COUNT):
        count, first = mapping[uni]
        for slot, span in enumerate(operation_index[uni]):
            if slot < count:
                if (span.order <= 0 or span.offset < 1
                        or span.offset + span.order > MSG_OPERATION_COUNT):
                    raise MagneticProvenanceInvariantError(
                        f"MSG operation span {uni}/{slot} is out of range"
                    )
                active_spans.append((span.offset, span.offset + span.order))
            elif span != OperationSpan(0, 0):
                raise MagneticProvenanceInvariantError(
                    f"MSG inactive operation span {uni}/{slot} is nonzero"
                )
        group_metadata = metadata[uni]
        if group_metadata is None:
            raise MagneticProvenanceInvariantError("MSG active metadata is empty")
        if spg.spacegroup_numbers[first] != group_metadata.parent_spacegroup:
            raise MagneticProvenanceInvariantError(
                f"UNI {uni} parent spacegroup disagrees with SPG"
            )

    if len(active_spans) != MSG_ACTIVE_SPAN_COUNT:
        raise MagneticProvenanceInvariantError("MSG active span census mismatch")
    previous_end = 1
    for start, end in sorted(active_spans):
        if start != previous_end or end <= start:
            raise MagneticProvenanceInvariantError("MSG spans are not contiguous")
        previous_end = end
    if previous_end != MSG_OPERATION_COUNT:
        raise MagneticProvenanceInvariantError("MSG span boundary mismatch")

    for uni in range(1, MSG_UNI_COUNT):
        count, first = mapping[uni]
        for slot in range(count):
            span = operation_index[uni][slot]
            group = decoded_operations[span.offset:span.offset + span.order]
            _validate_group(group, f"MSG UNI {uni} Hall {first + slot}", True)

    type_counts = {kind: 0 for kind in (1, 2, 3, 4)}
    for value in metadata[1:]:
        if value is None:
            raise MagneticProvenanceInvariantError("MSG metadata unexpectedly empty")
        type_counts[int(value.kind)] += 1
    if type_counts != {1: 230, 2: 230, 3: 674, 4: 517}:
        raise MagneticProvenanceInvariantError("MSG type census mismatch")

    return MsgProvenance(
        metadata=metadata,
        uni_mapping=mapping,
        derived_hall_mapping=derived_hall_mapping,
        operation_index=operation_index,
        raw_operation_codes=raw_codes,
        decoded_operations=decoded_operations,
        alternative_codes=alternative_rows,
    )


def _validate_witnesses(database):
    msg = database.msg
    metadata = msg.metadata[7]
    if metadata != MagneticGroupMetadata(
        uni=7, litvin=7, bns="2.7", og="2.4.7", parent_spacegroup=2,
        kind=MagneticKind.ANTI_TRANSLATION,
    ):
        raise MagneticProvenanceInvariantError("UNI7 metadata witness mismatch")
    if msg.halls_for_uni(7) != (2,):
        raise MagneticProvenanceInvariantError("UNI7 Hall witness mismatch")
    if msg.magnetic_operation_span(7, 2) != OperationSpan(4, 14):
        raise MagneticProvenanceInvariantError("UNI7 span witness mismatch")
    if msg.raw_operation_codes[14:18] != tuple(
        (16484, 3198, 34146806, 34133520)
    ):
        raise MagneticProvenanceInvariantError("UNI7 raw witness mismatch")
    uni7_operations = msg.magnetic_operations(7, 2)
    if tuple(operation.rotation for operation in uni7_operations) != (
        _IDENTITY_ROTATION,
        ((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
        _IDENTITY_ROTATION,
        ((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
    ) or tuple(operation.translation_numerator for operation in uni7_operations) != (
        (0, 0, 0), (0, 0, 0), (0, 0, 6), (0, 0, 6)
    ) or tuple(operation.time_reversal for operation in uni7_operations) != (
        TimeReversal.UNITARY, TimeReversal.UNITARY,
        TimeReversal.ANTIUNITARY, TimeReversal.ANTIUNITARY,
    ):
        raise MagneticProvenanceInvariantError("UNI7 decoded witness mismatch")
    if msg.raw_alternative_codes(7, 2) != (30, 90, 111, 810, 2301, 6831):
        raise MagneticProvenanceInvariantError("UNI7 alternative witness mismatch")
    if len(msg.std_transformations(7, 2)) != 7:
        raise MagneticProvenanceInvariantError("UNI7 transformation witness mismatch")
    if msg.std_transformations(7, 2)[0] != _IDENTITY_OPERATION:
        raise MagneticProvenanceInvariantError("UNI7 identity transformation mismatch")

    metadata = msg.metadata[9]
    if metadata != MagneticGroupMetadata(
        uni=9, litvin=9, bns="3.2", og="3.2.9", parent_spacegroup=3,
        kind=MagneticKind.GREY,
    ):
        raise MagneticProvenanceInvariantError("UNI9 metadata witness mismatch")
    if msg.halls_for_uni(9) != (3, 4, 5):
        raise MagneticProvenanceInvariantError("UNI9 Hall witness mismatch")
    expected_spans = (OperationSpan(4, 20), OperationSpan(4, 40),
                      OperationSpan(4, 60))
    expected_raw_groups = (
        (16484, 3360, 34028708, 34015584),
        (16484, 3200, 34028708, 34015424),
        (16484, 16320, 34028708, 34028544),
    )
    if tuple(msg.magnetic_operation_span(9, hall) for hall in (3, 4, 5)) != expected_spans:
        raise MagneticProvenanceInvariantError("UNI9 span witness mismatch")
    for hall, expected_raw in zip((3, 4, 5), expected_raw_groups):
        operations = msg.magnetic_operations(9, hall)
        if tuple(operation.encoded for operation in operations) != expected_raw:
            raise MagneticProvenanceInvariantError("UNI9 raw witness mismatch")
        if any(operation.translation_numerator != (0, 0, 0)
               for operation in operations):
            raise MagneticProvenanceInvariantError("UNI9 translation witness mismatch")
        if sum(operation.time_reversal is TimeReversal.UNITARY
               for operation in operations) != 2:
            raise MagneticProvenanceInvariantError("UNI9 unitary witness mismatch")
        if sum(operation.time_reversal is TimeReversal.ANTIUNITARY
               for operation in operations) != 2:
            raise MagneticProvenanceInvariantError("UNI9 antiunitary witness mismatch")
        if ExactSeitzOperation(
            encoded=34028708,
            rotation=_IDENTITY_ROTATION,
            translation_numerator=(0, 0, 0),
            time_reversal=TimeReversal.ANTIUNITARY,
        ) not in operations:
            raise MagneticProvenanceInvariantError("UNI9 pure-theta witness mismatch")
        if msg.raw_alternative_codes(9, hall) != ():
            raise MagneticProvenanceInvariantError("UNI9 alternatives witness mismatch")
        if msg.std_transformations(9, hall) != (_IDENTITY_OPERATION,):
            raise MagneticProvenanceInvariantError("UNI9 transformation witness mismatch")


def _build_database(data):
    try:
        spg = _convert_spg(data)
        msg = _convert_msg(data, spg)
        database = MagneticProvenanceDatabase(spg=spg, msg=msg)
        _validate_witnesses(database)
        return database
    except MagneticProvenanceError:
        raise
    except (KeyError, IndexError, TypeError, ValueError) as error:
        raise MagneticProvenanceSchemaError("typed artifact conversion failed") from error


def _verify_fixed_bytes(data, length, digest, label):
    if type(data) is not bytes:
        raise MagneticProvenanceIntegrityError(f"{label} is not bytes")
    if len(data) != length:
        raise MagneticProvenanceIntegrityError(
            f"{label} length mismatch: {len(data)} != {length}"
        )
    actual = hashlib.sha256(data).hexdigest()
    if actual != digest:
        raise MagneticProvenanceIntegrityError(
            f"{label} SHA256 mismatch: {actual}"
        )


def _from_bytes(artifact_bytes, manifest_bytes, artifact_name=_ARTIFACT_NAME):
    if artifact_name != _ARTIFACT_NAME:
        raise MagneticProvenanceIntegrityError("artifact name is not committed")
    _verify_fixed_bytes(artifact_bytes, ARTIFACT_BYTE_LENGTH,
                        ARTIFACT_SHA256, "artifact")
    _verify_fixed_bytes(manifest_bytes, MANIFEST_BYTE_LENGTH,
                        MANIFEST_SHA256, "manifest")
    return _from_pair_bytes(artifact_bytes, manifest_bytes, artifact_name)


def _from_pair_bytes(artifact_bytes, manifest_bytes, artifact_name):
    try:
        artifact = _extractor.parse_and_validate_committed_pair(
            artifact_bytes, manifest_bytes, artifact_name
        )
    except _extractor.ExtractionError as error:
        raise MagneticProvenanceSchemaError(
            "committed artifact/manifest schema validation failed"
        ) from error
    return _build_database(artifact)


def _from_uncommitted_pair_for_test(
    artifact_bytes, manifest_bytes, artifact_name=_ARTIFACT_NAME
):
    """Build a typed database from a test-mutated, synchronized byte pair."""
    if artifact_name != _ARTIFACT_NAME:
        raise MagneticProvenanceIntegrityError("artifact name is not committed")
    return _from_pair_bytes(artifact_bytes, manifest_bytes, artifact_name)


def _reset_cache_for_test():
    """Clear the singleton for isolated tests; never part of the public API."""
    global _CACHED_DATABASE
    with _CACHE_LOCK:
        _CACHED_DATABASE = None


def load_committed_provenance() -> MagneticProvenanceDatabase:
    """Load and validate only the committed repository artifact pair."""
    global _CACHED_DATABASE
    cached = _CACHED_DATABASE
    if cached is not None:
        return cached
    with _CACHE_LOCK:
        cached = _CACHED_DATABASE
        if cached is not None:
            return cached
        artifact_path = _DATA_DIR / _ARTIFACT_NAME
        manifest_path = _DATA_DIR / _MANIFEST_NAME
        try:
            artifact_bytes = artifact_path.read_bytes()
        except OSError as error:
            raise MagneticProvenanceIntegrityError(
                f"unable to read committed artifact {artifact_path}"
            ) from error
        _verify_fixed_bytes(artifact_bytes, ARTIFACT_BYTE_LENGTH,
                            ARTIFACT_SHA256, "artifact")
        try:
            manifest_bytes = manifest_path.read_bytes()
        except OSError as error:
            raise MagneticProvenanceIntegrityError(
                f"unable to read committed manifest {manifest_path}"
            ) from error
        _verify_fixed_bytes(manifest_bytes, MANIFEST_BYTE_LENGTH,
                            MANIFEST_SHA256, "manifest")
        database = _from_bytes(artifact_bytes, manifest_bytes, _ARTIFACT_NAME)
        _CACHED_DATABASE = database
        return database
