"""Canonical runtime-parity frame for the committed magnetic database.

The expected frame is deliberately built only from the typed, immutable
provenance loader.  The Rust side emits the same wire format from its public
raw tables and runtime accessors; neither side reads the other's source.
"""

from __future__ import annotations

import struct
from typing import Iterable

try:
    from . import spglib_magnetic_provenance as provenance
except ImportError:
    import spglib_magnetic_provenance as provenance


MAGIC = b"CRYMDBP\0"
VERSION = 1
SECTION_COUNT = 13
SECTION_TAGS = (
    "SGNO", "SGIX", "SGRW", "MTYP", "MUNI", "MHLL", "MIDX",
    "MRAW", "MALT", "SDEC", "SAPI", "MAPI", "TAPI",
)
SECTION_COUNTS = (
    531, 531, 8147, 1652, 1652, 531, 1652 * 18, 76683,
    1652 * 18, 8146, 530, 4479, 4479,
)
SPG_HALL_COUNT = 531
SPG_HALL_SETTINGS = 530
SPG_OPERATION_COUNT = 8147
MSG_UNI_COUNT = 1652
MSG_HALL_SLOTS = 18
MSG_OPERATION_COUNT = 76683
MSG_ACTIVE_SPAN_COUNT = 4479
ALT_VALUE_COUNT = 536
TRANSLATION_DENOMINATOR = 12

GOLDEN_FRAME_LENGTH = 2_862_892
GOLDEN_FRAME_SHA256 = (
    "4c0c0e488c5826369502631799f92bae868427c82d46155e55728e9f65ba14b3"
)
GOLDEN_SECTION_PAYLOAD_SHA256 = (
    "5a2128fb9295556c5e6b8128d6de36cc610e11ee233207ffffe3f278766bea79",
    "2e98578749f3b305429bd10292326637c603d779de17e36a7a1ab9832a80d8ad",
    "7636d8394ce8a5f0d08bca4547cf8adb8b77d2d1be34e95cb94613d83b0abf1f",
    "205eb282403885779ad029a853efbc832893edbd9e8d0b777c59bfa3749eb26d",
    "3747baa65d7ceb373744adbdd1baf7c6c830436e0cac993bf0b50dc9274afd98",
    "a1e942d246a88eeb142d2ac7013bec7b20ee6d14e118e30cb8ac2ca0528e41aa",
    "5353c68604f6e69e0f069a08a9d3aa510bba814c81eb3a22d75b31796591869a",
    "07ec2870830f577687e0ddfc577d59df8bd93fe4e9f1fd518d31c663d3010df6",
    "6e0e4ffd1746d0cdf1ec16efd3f3a64cd177eb4195d3aae5669048d277d2a28d",
    "87f7d2c424073f1940f578075977b4898925b776e094b5446deda2f25b016dd9",
    "634f51515caf52e8f154089522834ea4094673fc92719a41db9b7e30ea21fb93",
    "9ca6b204df245cac6198e55e9f500aa53dc21c291ce27df70ebc36502c31e5c0",
    "ee5353329cbef42d66eda407018bb52b74b51c9766e6cc1658c70a6f34d00685",
)


class FrameError(ValueError):
    """Malformed canonical frame or an impossible typed value."""


def _pack_i32(value: int, label: str) -> bytes:
    if type(value) is not int or not -(1 << 31) <= value < (1 << 31):
        raise FrameError(f"{label} is not a signed i32")
    return struct.pack("<i", value)


def _pack_u8(value: int, label: str) -> bytes:
    if type(value) is not int or not 0 <= value <= 0xFF:
        raise FrameError(f"{label} is not a u8")
    return struct.pack("<B", value)


def _pack_u16(value: int, label: str) -> bytes:
    if type(value) is not int or not 0 <= value <= 0xFFFF:
        raise FrameError(f"{label} is not a u16")
    return struct.pack("<H", value)


def _pack_u32(value: int, label: str) -> bytes:
    if type(value) is not int or not 0 <= value <= 0xFFFFFFFF:
        raise FrameError(f"{label} is not a u32")
    return struct.pack("<I", value)


def _pack_u64(value: int, label: str) -> bytes:
    if type(value) is not int or not 0 <= value <= 0xFFFFFFFFFFFFFFFF:
        raise FrameError(f"{label} is not a u64")
    return struct.pack("<Q", value)


def _pack_string(value: str, label: str) -> bytes:
    if type(value) is not str:
        raise FrameError(f"{label} is not a string")
    if "\0" in value:
        raise FrameError(f"{label} contains NUL")
    try:
        encoded = value.encode("utf-8", "strict")
    except UnicodeError as error:
        raise FrameError(f"{label} is not valid UTF-8") from error
    return _pack_u32(len(encoded), f"{label} length") + encoded


def _flat_rotation(operation, label: str) -> tuple[int, ...]:
    rotation = operation.rotation
    if type(rotation) is not tuple or len(rotation) != 3:
        raise FrameError(f"{label} rotation shape mismatch")
    flat = []
    for row_index, row in enumerate(rotation):
        if type(row) is not tuple or len(row) != 3:
            raise FrameError(f"{label} rotation row shape mismatch")
        for column_index, value in enumerate(row):
            if type(value) is not int or not -128 <= value <= 127:
                raise FrameError(
                    f"{label} rotation[{row_index}][{column_index}] is not i8"
                )
            flat.append(value)
    return tuple(flat)


def _operation_payload(operation, magnetic: bool, label: str) -> bytes:
    flat = _flat_rotation(operation, label)
    translation = operation.translation_numerator
    if (type(translation) is not tuple or len(translation) != 3
            or any(type(value) is not int or not 0 <= value < TRANSLATION_DENOMINATOR
                   for value in translation)):
        raise FrameError(f"{label} translation is not three u8 numerators")
    result = bytearray(struct.pack("<9b", *flat))
    result.extend(struct.pack("<3B", *translation))
    if magnetic:
        time_reversal = int(operation.time_reversal)
        if time_reversal not in (0, 1):
            raise FrameError(f"{label} time reversal is not 0/1")
        result.extend(_pack_u8(time_reversal, f"{label} time reversal"))
    return bytes(result)


def _section(tag: str, record_count: int, payload: bytes) -> bytes:
    if tag not in SECTION_TAGS:
        raise FrameError(f"unknown section tag {tag!r}")
    if type(payload) is not bytes:
        raise FrameError(f"{tag} payload is not bytes")
    return tag.encode("ascii") + _pack_u64(record_count, f"{tag} record count") + _pack_u64(
        len(payload), f"{tag} payload length"
    ) + payload


def _span_payload(span, label: str) -> bytes:
    return _pack_i32(span.order, f"{label}.order") + _pack_i32(
        span.offset, f"{label}.offset"
    )


def _append_spans(payload: bytearray, rows: Iterable, label: str) -> None:
    for index, span in enumerate(rows):
        payload.extend(_span_payload(span, f"{label}[{index}]"))


def _build_sections(database) -> tuple[bytes, ...]:
    spg = database.spg
    msg = database.msg

    sgno = bytearray()
    for hall, number in enumerate(spg.spacegroup_numbers):
        sgno.extend(_pack_i32(number, f"SGNO[{hall}]"))

    sgix = bytearray()
    _append_spans(sgix, spg.operation_index, "SGIX")

    sgrw = bytearray()
    for index, code in enumerate(spg.raw_operation_codes):
        sgrw.extend(_pack_i32(code, f"SGRW[{index}]"))

    mtyp = bytearray()
    for uni, metadata in enumerate(msg.metadata):
        if metadata is None:
            mtyp.extend(struct.pack("<4i", 0, 0, 0, 0))
            mtyp.extend(_pack_string("", f"MTYP[{uni}].bns"))
            mtyp.extend(_pack_string("", f"MTYP[{uni}].og"))
            continue
        mtyp.extend(_pack_i32(metadata.uni, f"MTYP[{uni}].uni"))
        mtyp.extend(_pack_i32(metadata.litvin, f"MTYP[{uni}].litvin"))
        mtyp.extend(_pack_i32(
            metadata.parent_spacegroup, f"MTYP[{uni}].parent_spacegroup"
        ))
        mtyp.extend(_pack_i32(int(metadata.kind), f"MTYP[{uni}].type"))
        mtyp.extend(_pack_string(metadata.bns, f"MTYP[{uni}].bns"))
        mtyp.extend(_pack_string(metadata.og, f"MTYP[{uni}].og"))

    muni = bytearray()
    for index, (count, first) in enumerate(msg.uni_mapping):
        muni.extend(_pack_i32(count, f"MUNI[{index}].count"))
        muni.extend(_pack_i32(first, f"MUNI[{index}].first"))

    mhll = bytearray()
    for hall, (smallest, largest) in enumerate(msg.derived_hall_mapping):
        mhll.extend(_pack_i32(smallest, f"MHLL[{hall}].smallest"))
        mhll.extend(_pack_i32(largest, f"MHLL[{hall}].largest"))

    midx = bytearray()
    for uni, row in enumerate(msg.operation_index):
        for slot, span in enumerate(row):
            midx.extend(_span_payload(span, f"MIDX[{uni}][{slot}]"))

    mraw = bytearray()
    for index, code in enumerate(msg.raw_operation_codes):
        mraw.extend(_pack_i32(code, f"MRAW[{index}]"))

    malt = bytearray()
    for uni, row in enumerate(msg.alternative_codes):
        for slot, prefix in enumerate(row):
            if len(prefix) > 6:
                raise FrameError(f"MALT[{uni}][{slot}] has no unique padding")
            values = tuple(prefix) + (0,) * (7 - len(prefix))
            for index, code in enumerate(values):
                malt.extend(_pack_i32(code, f"MALT[{uni}][{slot}][{index}]"))

    sdec = bytearray()
    for index in range(1, len(spg.raw_operation_codes)):
        operation = spg.decoded_operations[index]
        if operation is None:
            raise FrameError(f"SDEC[{index}] crossed sentinel")
        sdec.extend(_pack_u32(index, f"SDEC[{index}].raw_index"))
        sdec.extend(_pack_i32(spg.raw_operation_codes[index], f"SDEC[{index}].encoded"))
        sdec.extend(_operation_payload(operation, False, f"SDEC[{index}]"))

    sapi = bytearray()
    for hall in range(1, SPG_HALL_SETTINGS + 1):
        operations = database.spg_operations(hall)
        sapi.extend(_pack_u16(hall, f"SAPI[{hall}].hall"))
        sapi.extend(_pack_u16(len(operations), f"SAPI[{hall}].count"))
        for index, operation in enumerate(operations):
            sapi.extend(_operation_payload(operation, False, f"SAPI[{hall}][{index}]"))

    mapi = bytearray()
    for uni in range(1, MSG_UNI_COUNT):
        for hall in database.halls_for_uni(uni):
            operations = database.magnetic_operations(uni, hall)
            mapi.extend(_pack_u16(uni, f"MAPI[{uni},{hall}].uni"))
            mapi.extend(_pack_u16(hall, f"MAPI[{uni},{hall}].hall"))
            mapi.extend(_pack_u16(len(operations), f"MAPI[{uni},{hall}].count"))
            for index, operation in enumerate(operations):
                mapi.extend(_operation_payload(
                    operation, True, f"MAPI[{uni},{hall}][{index}]"
                ))

    tapi = bytearray()
    alternative_occurrences = 0
    for uni in range(1, MSG_UNI_COUNT):
        for hall in database.halls_for_uni(uni):
            transformations = database.std_transformations(uni, hall)
            prefix = database.raw_alternative_codes(uni, hall)
            if len(transformations) != len(prefix) + 1:
                raise FrameError(f"TAPI[{uni},{hall}] transformation count mismatch")
            alternative_occurrences += len(prefix)
            tapi.extend(_pack_u16(uni, f"TAPI[{uni},{hall}].uni"))
            tapi.extend(_pack_u16(hall, f"TAPI[{uni},{hall}].hall"))
            tapi.extend(_pack_u16(
                len(transformations), f"TAPI[{uni},{hall}].count"
            ))
            for index, operation in enumerate(transformations):
                if int(operation.time_reversal) != 0:
                    raise FrameError(f"TAPI[{uni},{hall}][{index}] is antiunitary")
                tapi.extend(_operation_payload(
                    operation, False, f"TAPI[{uni},{hall}][{index}]"
                ))
    if alternative_occurrences != ALT_VALUE_COUNT:
        raise FrameError(
            f"TAPI raw alternative occurrence count {alternative_occurrences}"
            f" != {ALT_VALUE_COUNT}"
        )

    return (
        _section("SGNO", SPG_HALL_COUNT, bytes(sgno)),
        _section("SGIX", SPG_HALL_COUNT, bytes(sgix)),
        _section("SGRW", SPG_OPERATION_COUNT, bytes(sgrw)),
        _section("MTYP", MSG_UNI_COUNT, bytes(mtyp)),
        _section("MUNI", MSG_UNI_COUNT, bytes(muni)),
        _section("MHLL", SPG_HALL_COUNT, bytes(mhll)),
        _section("MIDX", MSG_UNI_COUNT * MSG_HALL_SLOTS, bytes(midx)),
        _section("MRAW", MSG_OPERATION_COUNT, bytes(mraw)),
        _section("MALT", MSG_UNI_COUNT * MSG_HALL_SLOTS, bytes(malt)),
        _section("SDEC", SPG_OPERATION_COUNT - 1, bytes(sdec)),
        _section("SAPI", SPG_HALL_SETTINGS, bytes(sapi)),
        _section("MAPI", MSG_ACTIVE_SPAN_COUNT, bytes(mapi)),
        _section("TAPI", MSG_ACTIVE_SPAN_COUNT, bytes(tapi)),
    )


def build_expected_frame(database=None) -> bytes:
    """Build the canonical expected frame from the immutable typed loader."""
    if database is None:
        database = provenance.load_committed_provenance()
    sections = _build_sections(database)
    return MAGIC + struct.pack("<II", VERSION, SECTION_COUNT) + b"".join(sections)


def parse_frame(frame: bytes) -> tuple[tuple[str, int, bytes], ...]:
    """Strictly parse a canonical frame and return `(tag, count, payload)` rows."""
    if type(frame) is not bytes:
        raise FrameError("frame is not bytes")
    header_length = len(MAGIC) + 8
    if len(frame) < header_length or frame[:len(MAGIC)] != MAGIC:
        raise FrameError("frame magic mismatch")
    version, section_count = struct.unpack_from("<II", frame, len(MAGIC))
    if version != VERSION or section_count != SECTION_COUNT:
        raise FrameError("frame version/section count mismatch")
    offset = header_length
    result = []
    for index, expected_tag in enumerate(SECTION_TAGS):
        if offset + 20 > len(frame):
            raise FrameError(f"section {index} header is truncated")
        tag_bytes = frame[offset:offset + 4]
        try:
            tag = tag_bytes.decode("ascii")
        except UnicodeDecodeError as error:
            raise FrameError(f"section {index} tag is not ASCII") from error
        if tag != expected_tag:
            raise FrameError(f"section {index} tag {tag!r} != {expected_tag!r}")
        count, payload_length = struct.unpack_from("<QQ", frame, offset + 4)
        offset += 20
        end = offset + payload_length
        if end > len(frame):
            raise FrameError(f"section {tag} payload is truncated")
        payload = frame[offset:end]
        result.append((tag, count, payload))
        offset = end
    if offset != len(frame):
        raise FrameError("frame has trailing bytes")
    return tuple(result)


def expected_section_metadata(database=None):
    """Return parsed sections, useful to tests and first-run diagnostics."""
    return parse_frame(build_expected_frame(database))


__all__ = [
    "ALT_VALUE_COUNT", "FrameError", "GOLDEN_FRAME_LENGTH",
    "GOLDEN_FRAME_SHA256", "GOLDEN_SECTION_PAYLOAD_SHA256", "MAGIC",
    "SECTION_COUNT", "SECTION_COUNTS", "SECTION_TAGS", "VERSION",
    "build_expected_frame", "expected_section_metadata", "parse_frame",
]
