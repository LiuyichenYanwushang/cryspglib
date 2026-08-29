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
_FIXED_SECTION_WIDTHS = {
    "SGNO": 4,
    "SGIX": 8,
    "SGRW": 4,
    "MUNI": 8,
    "MHLL": 8,
    "MIDX": 8,
    "MRAW": 4,
    "MALT": 28,
    "SDEC": 20,
}
SPG_HALL_COUNT = 531
SPG_HALL_SETTINGS = 530
SPG_OPERATION_COUNT = 8147
SPG_STANDARD_OPERATION_END = 7389
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


def _require_payload(payload: bytes, offset: int, length: int, label: str) -> int:
    end = offset + length
    if end > len(payload):
        raise FrameError(f"{label} is truncated")
    return end


def _parse_mtyp_records(payload: bytes, record_count: int):
    offset = 0
    records = []
    for record in range(record_count):
        offset = _require_payload(payload, offset, 16, f"MTYP[{record}] header")
        values = _unpack_from(payload, "<4i", offset - 16, f"MTYP[{record}] header")
        strings = []
        for field in ("bns", "og"):
            offset = _require_payload(payload, offset, 4, f"MTYP[{record}].{field} length")
            length = struct.unpack_from("<I", payload, offset - 4)[0]
            value_offset = offset
            end = _require_payload(payload, value_offset, length, f"MTYP[{record}].{field}")
            value = payload[value_offset:end]
            if b"\0" in value:
                raise FrameError(f"MTYP[{record}].{field} contains NUL")
            try:
                value.decode("utf-8", "strict")
            except UnicodeDecodeError as error:
                raise FrameError(f"MTYP[{record}].{field} is not UTF-8") from error
            strings.append(value.decode("utf-8", "strict"))
            offset = end
        records.append((*values, *strings))
    if offset != len(payload):
        raise FrameError("MTYP has trailing bytes")
    return tuple(records)


def _validate_mtyp_payload(payload: bytes, record_count: int) -> None:
    _parse_mtyp_records(payload, record_count)


def _validate_operation_payload(
    payload: bytes, offset: int, count: int, magnetic: bool, label: str
) -> int:
    width = 13 if magnetic else 12
    for index in range(count):
        end = _require_payload(payload, offset, width, f"{label}[{index}]")
        rotation = struct.unpack_from("<9b", payload, offset)
        if any(value not in (-1, 0, 1) for value in rotation):
            raise FrameError(f"{label}[{index}] rotation is outside -1..1")
        translation = payload[offset + 9:offset + 12]
        if any(value >= TRANSLATION_DENOMINATOR for value in translation):
            raise FrameError(f"{label}[{index}] translation is out of range")
        if magnetic and payload[offset + 12] not in (0, 1):
            raise FrameError(f"{label}[{index}] time reversal is out of range")
        offset = end
    return offset


def _validate_sapi_payload(payload: bytes, record_count: int) -> None:
    offset = 0
    for record in range(record_count):
        offset = _require_payload(payload, offset, 4, f"SAPI[{record}] header")
        hall, count = struct.unpack_from("<HH", payload, offset - 4)
        if hall != record + 1 or count == 0:
            raise FrameError(f"SAPI[{record}] header is not canonical")
        offset = _validate_operation_payload(payload, offset, count, False, f"SAPI[{hall}]")
    if offset != len(payload):
        raise FrameError("SAPI has trailing bytes")


def _validate_msg_variable_payload(
    payload: bytes, record_count: int, magnetic: bool, label: str
) -> None:
    offset = 0
    previous = None
    for record in range(record_count):
        offset = _require_payload(payload, offset, 6, f"{label}[{record}] header")
        uni, hall, count = struct.unpack_from("<HHH", payload, offset - 6)
        if not 1 <= uni < MSG_UNI_COUNT or not 1 <= hall < SPG_HALL_COUNT or count == 0:
            raise FrameError(f"{label}[{record}] header is out of range")
        if previous is None:
            if (uni, hall) != (1, 1):
                raise FrameError(f"{label}[{record}] does not start at UNI1/Hall1")
        elif uni == previous[0]:
            if hall != previous[1] + 1:
                raise FrameError(f"{label}[{record}] Hall sequence is not contiguous")
        elif uni == previous[0] + 1:
            # A new UNI starts at its first active Hall; the exact mapping is
            # checked by the typed loader/Rust parity comparison.
            pass
        else:
            raise FrameError(f"{label}[{record}] UNI sequence is not contiguous")
        offset = _validate_operation_payload(payload, offset, count, magnetic, f"{label}[{uni},{hall}]")
        previous = (uni, hall)
    if offset != len(payload):
        raise FrameError(f"{label} has trailing bytes")


def _validate_section_payload(tag: str, count: int, payload: bytes) -> None:
    if tag in _FIXED_SECTION_WIDTHS:
        expected_length = count * _FIXED_SECTION_WIDTHS[tag]
        if len(payload) != expected_length:
            raise FrameError(
                f"{tag} payload length {len(payload)} != {expected_length}"
            )
        return
    if tag == "MTYP":
        _validate_mtyp_payload(payload, count)
    elif tag == "SAPI":
        _validate_sapi_payload(payload, count)
    elif tag == "MAPI":
        _validate_msg_variable_payload(payload, count, True, tag)
    elif tag == "TAPI":
        _validate_msg_variable_payload(payload, count, False, tag)
    else:
        raise FrameError(f"no parser for section {tag!r}")


def _unpack_from(payload: bytes, format_string: str, offset: int, label: str):
    try:
        return struct.unpack_from(format_string, payload, offset)
    except (IndexError, struct.error) as error:
        raise FrameError(f"{label} is truncated") from error


def _parse_i32_rows(payload: bytes, row_count: int, width: int, label: str):
    rows = []
    for row in range(row_count):
        offset = row * width * 4
        values = _unpack_from(
            payload, "<" + "i" * width, offset, f"{label}[{row}]"
        )
        rows.append(tuple(values))
    return tuple(rows)


def _parse_operation_payload(
    payload: bytes, offset: int, magnetic: bool, label: str
):
    width = 13 if magnetic else 12
    end = _require_payload(payload, offset, width, label)
    rotation = tuple(_unpack_from(payload, "<9b", offset, f"{label}.rotation"))
    if any(value not in (-1, 0, 1) for value in rotation):
        raise FrameError(f"{label} rotation is outside -1..1")
    translation = tuple(payload[offset + 9:offset + 12])
    if any(value >= TRANSLATION_DENOMINATOR for value in translation):
        raise FrameError(f"{label} translation is out of range")
    time_reversal = payload[offset + 12] if magnetic else None
    if magnetic and time_reversal not in (0, 1):
        raise FrameError(f"{label} time reversal is out of range")
    determinant = (
        rotation[0] * (rotation[4] * rotation[8] - rotation[5] * rotation[7])
        - rotation[1] * (rotation[3] * rotation[8] - rotation[5] * rotation[6])
        + rotation[2] * (rotation[3] * rotation[7] - rotation[4] * rotation[6])
    )
    if determinant not in (-1, 1):
        raise FrameError(f"{label} rotation determinant is not ±1")
    return rotation, translation, time_reversal, end


def _decode_exact_operation(encoded: int, label: str):
    try:
        return provenance._decode_operation(encoded)
    except provenance.MagneticProvenanceError as error:
        raise FrameError(f"{label} cannot be decoded exactly") from error


def _compare_wire_operation(wire, expected, magnetic: bool, label: str) -> None:
    rotation, translation, time_reversal = wire
    expected_rotation = tuple(
        value for row in expected.rotation for value in row
    )
    if rotation != expected_rotation:
        raise FrameError(f"{label} rotation does not match encoded operation")
    if translation != expected.translation_numerator:
        raise FrameError(f"{label} translation does not match encoded operation")
    expected_time = int(expected.time_reversal)
    if magnetic:
        if time_reversal != expected_time:
            raise FrameError(f"{label} time reversal does not match encoded operation")
    elif time_reversal is not None or expected_time != 0:
        raise FrameError(f"{label} unexpectedly has time reversal")


def _parse_sdec_relations(payloads) -> None:
    sgrw = _parse_i32_rows(payloads["SGRW"], SPG_OPERATION_COUNT, 1, "SGRW")
    raw_codes = tuple(row[0] for row in sgrw)
    if raw_codes[0] != 0:
        raise FrameError("SGRW sentinel mismatch")
    offset = 0
    for index in range(1, SPG_OPERATION_COUNT):
        raw_index = _unpack_from(payloads["SDEC"], "<I", offset, f"SDEC[{index}].index")[0]
        encoded = _unpack_from(payloads["SDEC"], "<i", offset + 4, f"SDEC[{index}].encoded")[0]
        if raw_index != index:
            raise FrameError(f"SDEC[{index}] raw index is not canonical")
        if encoded != raw_codes[index]:
            raise FrameError(f"SDEC[{index}] does not match SGRW[{index}]")
        wire = _parse_operation_payload(
            payloads["SDEC"], offset + 8, False, f"SDEC[{index}]"
        )
        expected = _decode_exact_operation(encoded, f"SDEC[{index}]")
        _compare_wire_operation(wire[:3], expected, False, f"SDEC[{index}]")
        offset += 20
    if offset != len(payloads["SDEC"]):
        raise FrameError("SDEC has trailing bytes")


def _checked_spg_spans(payloads):
    rows = _parse_i32_rows(payloads["SGIX"], SPG_HALL_COUNT, 2, "SGIX")
    if rows[0] != (0, 0):
        raise FrameError("SGIX sentinel mismatch")
    spans = []
    previous_end = 1
    for hall in range(1, SPG_HALL_COUNT):
        order, offset = rows[hall]
        if order <= 0 or offset < 1 or offset > SPG_STANDARD_OPERATION_END:
            raise FrameError(f"SGIX Hall {hall} span is out of range")
        if order > SPG_STANDARD_OPERATION_END - offset:
            raise FrameError(f"SGIX Hall {hall} span overflows standard table")
        if offset != previous_end:
            raise FrameError(f"SGIX Hall {hall} span is not adjacent")
        previous_end = offset + order
        spans.append((order, offset))
    if previous_end != SPG_STANDARD_OPERATION_END:
        raise FrameError("SGIX standard span boundary mismatch")
    return tuple(spans)


def _parse_sapi_relations(payloads, spg_spans) -> None:
    sgrw = _parse_i32_rows(payloads["SGRW"], SPG_OPERATION_COUNT, 1, "SGRW")
    raw_codes = tuple(row[0] for row in sgrw)
    offset = 0
    for hall, (order, raw_offset) in enumerate(spg_spans, 1):
        record_hall, record_count = _unpack_from(
            payloads["SAPI"], "<HH", offset, f"SAPI[{hall}] header"
        )
        if record_hall != hall or record_count != order:
            raise FrameError(f"SAPI[{hall}] header does not match SGIX")
        offset += 4
        for index in range(order):
            wire = _parse_operation_payload(
                payloads["SAPI"], offset, False, f"SAPI[{hall}][{index}]"
            )
            encoded = raw_codes[raw_offset + index]
            expected = _decode_exact_operation(encoded, f"SAPI[{hall}][{index}]")
            _compare_wire_operation(
                wire[:3], expected, False, f"SAPI[{hall}][{index}]"
            )
            offset += 12
    if offset != len(payloads["SAPI"]):
        raise FrameError("SAPI has trailing bytes")


def _checked_msg_spans(payloads):
    muni_rows = _parse_i32_rows(payloads["MUNI"], MSG_UNI_COUNT, 2, "MUNI")
    if muni_rows[0] != (0, 0):
        raise FrameError("MUNI sentinel mismatch")
    midx_rows = _parse_i32_rows(
        payloads["MIDX"], MSG_UNI_COUNT * MSG_HALL_SLOTS, 2, "MIDX"
    )
    if any(midx_rows[slot] != (0, 0) for slot in range(MSG_HALL_SLOTS)):
        raise FrameError("MIDX sentinel row is nonzero")
    active = []
    seen = set()
    for uni in range(1, MSG_UNI_COUNT):
        count, first = muni_rows[uni]
        if not 1 <= count <= MSG_HALL_SLOTS:
            raise FrameError(f"MUNI UNI {uni} Hall count is out of range")
        if not 1 <= first <= SPG_HALL_SETTINGS or first + count - 1 > SPG_HALL_SETTINGS:
            raise FrameError(f"MUNI UNI {uni} Hall range is out of range")
        for slot in range(MSG_HALL_SLOTS):
            order, offset = midx_rows[uni * MSG_HALL_SLOTS + slot]
            if slot >= count:
                if (order, offset) != (0, 0):
                    raise FrameError(f"MIDX UNI {uni} inactive slot {slot} is nonzero")
                continue
            if order <= 0 or offset < 1 or offset > MSG_OPERATION_COUNT:
                raise FrameError(f"MIDX UNI {uni} slot {slot} span is out of range")
            if order > MSG_OPERATION_COUNT - offset:
                raise FrameError(f"MIDX UNI {uni} slot {slot} span overflows table")
            for raw_index in range(offset, offset + order):
                if raw_index in seen:
                    raise FrameError(f"MIDX raw index {raw_index} is duplicated")
                seen.add(raw_index)
            active.append((uni, first + slot, order, offset))
    if len(active) != MSG_ACTIVE_SPAN_COUNT:
        raise FrameError("MIDX active span census mismatch")
    if seen != set(range(1, MSG_OPERATION_COUNT)):
        raise FrameError("MIDX spans do not cover MRAW exactly once")
    return tuple(active), muni_rows


def _validate_metadata_and_hall_relations(payloads, muni_rows) -> None:
    sgno_rows = _parse_i32_rows(payloads["SGNO"], SPG_HALL_COUNT, 1, "SGNO")
    spacegroup_numbers = tuple(row[0] for row in sgno_rows)
    if spacegroup_numbers[0] != 0 or any(
        not 1 <= value <= 230 for value in spacegroup_numbers[1:]
    ):
        raise FrameError("SGNO sentinel or range mismatch")

    metadata = _parse_mtyp_records(payloads["MTYP"], MSG_UNI_COUNT)
    if metadata[0] != (0, 0, 0, 0, "", ""):
        raise FrameError("MTYP sentinel mismatch")
    type_counts = {1: 0, 2: 0, 3: 0, 4: 0}
    for uni in range(1, MSG_UNI_COUNT):
        row_uni, litvin, parent, kind, _, _ = metadata[uni]
        if row_uni != uni or not 1 <= litvin < MSG_UNI_COUNT:
            raise FrameError(f"MTYP UNI {uni} identity/range mismatch")
        if not 1 <= parent <= 230 or kind not in type_counts:
            raise FrameError(f"MTYP UNI {uni} parent/type mismatch")
        type_counts[kind] += 1
        count, first = muni_rows[uni]
        for hall in range(first, first + count):
            if spacegroup_numbers[hall] != parent:
                raise FrameError(f"MTYP UNI {uni} parent disagrees at Hall {hall}")
    if type_counts != {1: 230, 2: 230, 3: 674, 4: 517}:
        raise FrameError("MTYP type census mismatch")

    mhll_rows = _parse_i32_rows(payloads["MHLL"], SPG_HALL_COUNT, 2, "MHLL")
    if mhll_rows[0] != (0, 0):
        raise FrameError("MHLL sentinel mismatch")
    for hall in range(1, SPG_HALL_SETTINGS + 1):
        expected_unis = tuple(
            uni for uni in range(1, MSG_UNI_COUNT)
            if muni_rows[uni][1] <= hall < muni_rows[uni][1] + muni_rows[uni][0]
        )
        if not expected_unis or expected_unis != tuple(
            range(expected_unis[0], expected_unis[-1] + 1)
        ):
            raise FrameError(f"MHLL Hall {hall} inverse range is not continuous")
        if mhll_rows[hall] != (expected_unis[0], expected_unis[-1]):
            raise FrameError(f"MHLL Hall {hall} does not match MUNI")


def _parse_mapi_relations(payloads, active_spans):
    mraw_rows = _parse_i32_rows(payloads["MRAW"], MSG_OPERATION_COUNT, 1, "MRAW")
    raw_codes = tuple(row[0] for row in mraw_rows)
    if raw_codes[0] != 0:
        raise FrameError("MRAW sentinel mismatch")
    offset = 0
    for record, (uni, hall, order, raw_offset) in enumerate(active_spans):
        record_uni, record_hall, record_count = _unpack_from(
            payloads["MAPI"], "<HHH", offset, f"MAPI[{record}] header"
        )
        if (record_uni, record_hall, record_count) != (uni, hall, order):
            raise FrameError(f"MAPI[{record}] header does not match MUNI/MIDX")
        offset += 6
        for index in range(order):
            wire = _parse_operation_payload(
                payloads["MAPI"], offset, True, f"MAPI[{uni},{hall}][{index}]"
            )
            encoded = raw_codes[raw_offset + index]
            expected = _decode_exact_operation(encoded, f"MAPI[{uni},{hall}][{index}]")
            _compare_wire_operation(
                wire[:3], expected, True, f"MAPI[{uni},{hall}][{index}]"
            )
            offset += 13
    if offset != len(payloads["MAPI"]):
        raise FrameError("MAPI has trailing bytes")


def _checked_alternative_prefixes(payloads, muni_rows):
    malt_rows = _parse_i32_rows(
        payloads["MALT"], MSG_UNI_COUNT * MSG_HALL_SLOTS, 7, "MALT"
    )
    prefixes = {}
    occurrences = 0
    for uni in range(MSG_UNI_COUNT):
        active_count = 0 if uni == 0 else muni_rows[uni][0]
        for slot in range(MSG_HALL_SLOTS):
            values = malt_rows[uni * MSG_HALL_SLOTS + slot]
            first_zero = next(
                (index for index, value in enumerate(values) if value == 0),
                len(values),
            )
            if slot >= active_count:
                if any(value != 0 for value in values):
                    raise FrameError(f"MALT UNI {uni} inactive slot {slot} is nonzero")
                continue
            if first_zero > 6 or any(value != 0 for value in values[first_zero:]):
                raise FrameError(f"MALT UNI {uni} slot {slot} terminator/tail is invalid")
            prefix = values[:first_zero]
            if any(not 0 < value < provenance.SPACE_OPERATION_SCALE for value in prefix):
                raise FrameError(f"MALT UNI {uni} slot {slot} encoding is invalid")
            prefixes[(uni, slot)] = prefix
            occurrences += len(prefix)
    if occurrences != ALT_VALUE_COUNT:
        raise FrameError("MALT alternative occurrence census mismatch")
    return prefixes


def _parse_tapi_relations(payloads, active_spans, muni_rows):
    prefixes = _checked_alternative_prefixes(payloads, muni_rows)
    offset = 0
    for record, (uni, hall, expected_order, _) in enumerate(active_spans):
        slot = hall - muni_rows[uni][1]
        prefix = prefixes[(uni, slot)]
        record_uni, record_hall, record_count = _unpack_from(
            payloads["TAPI"], "<HHH", offset, f"TAPI[{record}] header"
        )
        if (record_uni, record_hall) != (uni, hall):
            raise FrameError(f"TAPI[{record}] header does not match MUNI/MIDX")
        if record_count != len(prefix) + 1:
            raise FrameError(f"TAPI[{uni},{hall}] count does not match MALT")
        if expected_order <= 0:
            raise FrameError(f"TAPI[{uni},{hall}] has invalid source span")
        offset += 6
        for index in range(record_count):
            wire = _parse_operation_payload(
                payloads["TAPI"], offset, False, f"TAPI[{uni},{hall}][{index}]"
            )
            if index == 0:
                if wire[:2] != (
                    (1, 0, 0, 0, 1, 0, 0, 0, 1), (0, 0, 0)
                ):
                    raise FrameError(f"TAPI[{uni},{hall}] first operation is not identity")
                expected = None
            else:
                encoded = prefix[index - 1]
                expected = _decode_exact_operation(
                    encoded, f"TAPI[{uni},{hall}][{index}]"
                )
                _compare_wire_operation(
                    wire[:3], expected, False, f"TAPI[{uni},{hall}][{index}]"
                )
            offset += 12
    if offset != len(payloads["TAPI"]):
        raise FrameError("TAPI has trailing bytes")


def _validate_cross_section_relations(payloads) -> None:
    _parse_sdec_relations(payloads)
    spg_spans = _checked_spg_spans(payloads)
    _parse_sapi_relations(payloads, spg_spans)
    active_spans, muni_rows = _checked_msg_spans(payloads)
    _validate_metadata_and_hall_relations(payloads, muni_rows)
    _parse_mapi_relations(payloads, active_spans)
    _parse_tapi_relations(payloads, active_spans, muni_rows)


def parse_frame(frame: bytes) -> tuple[tuple[str, int, bytes], ...]:
    """Strictly parse a canonical frame and return `(tag, count, payload)` rows."""
    try:
        if type(frame) is not bytes:
            raise FrameError("frame is not bytes")
        header_length = len(MAGIC) + 8
        if len(frame) < header_length or frame[:len(MAGIC)] != MAGIC:
            raise FrameError("frame magic mismatch")
        version, section_count = _unpack_from(
            frame, "<II", len(MAGIC), "frame header"
        )
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
            count, payload_length = _unpack_from(
                frame, "<QQ", offset + 4, f"section {tag} header"
            )
            expected_count = SECTION_COUNTS[index]
            if count != expected_count:
                raise FrameError(
                    f"section {tag} record count {count} != {expected_count}"
                )
            offset += 20
            end = offset + payload_length
            if end > len(frame):
                raise FrameError(f"section {tag} payload is truncated")
            payload = frame[offset:end]
            result.append((tag, count, payload))
            offset = end
        if offset != len(frame):
            raise FrameError("frame has trailing bytes")
        for tag, count, payload in result:
            _validate_section_payload(tag, count, payload)
        payloads = {tag: payload for tag, _, payload in result}
        _validate_cross_section_relations(payloads)
        return tuple(result)
    except FrameError:
        raise
    except (IndexError, MemoryError, OverflowError, struct.error, ValueError) as error:
        raise FrameError("frame validation failed") from error


def expected_section_metadata(database=None):
    """Return parsed sections, useful to tests and first-run diagnostics."""
    return parse_frame(build_expected_frame(database))


__all__ = [
    "ALT_VALUE_COUNT", "FrameError", "GOLDEN_FRAME_LENGTH",
    "GOLDEN_FRAME_SHA256", "GOLDEN_SECTION_PAYLOAD_SHA256", "MAGIC",
    "SECTION_COUNT", "SECTION_COUNTS", "SECTION_TAGS", "VERSION",
    "build_expected_frame", "expected_section_metadata", "parse_frame",
]
