#!/usr/bin/env python3
"""Extract pinned spglib magnetic database tables into a raw JSON artifact.

This is deliberately a small, strict C-initializer reader.  It consumes only
the two pinned upstream C files and never imports the generated Rust tables.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


SCHEMA = "cryspglib-spglib-magnetic-v1"
MANIFEST_SCHEMA = "cryspglib-spglib-magnetic-manifest-v1"
EXTRACTOR_VERSION = "1"
TRANSLATION_DENOMINATOR = 12
ROTATION_RADIX = 3
ROTATION_DIGITS = 9
ROTATION_PAYLOAD = ROTATION_RADIX ** ROTATION_DIGITS
TRANSLATION_DIGITS = 3
TRANSLATION_PAYLOAD = TRANSLATION_DENOMINATOR ** TRANSLATION_DIGITS
MSG_OPERATION_SCALE = 34_012_224
MAGNETIC_OPERATION_ENCODING_LIMIT = 2 * MSG_OPERATION_SCALE
SPG_OPERATION_COUNT = 8_147
SPG_STANDARD_OPERATION_END = 7_389
MSG_OPERATION_COUNT = 76_683
MSG_ACTIVE_SPAN_COUNT = 4_479
SPG_HALL_COUNT = 531
SPG_HALL_SETTINGS = SPG_HALL_COUNT - 1
MSG_UNI_COUNT = 1_652
MSG_HALL_SLOTS = 18

EXPECTED_SOURCES = {
    "msg_database.c": "2c45cf2e3827b48bab149ca59951f6b242f3754cc7d559585dfe107aa2a94288",
    "spg_database.c": "1ad4d3fd4ee2b39d43bf102bd6464bc6d1e07bc825dc4ac40df0235501823896",
}
EXPECTED_TYPE_COUNTS = {"1": 230, "2": 230, "3": 674, "4": 517}


class ExtractionError(ValueError):
    """Raised for any malformed or unexpected upstream initializer."""


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def _initializer_text(source, name):
    declaration = re.compile(r"\b" + re.escape(name) + r"\b(?:\s*\[[^\]]*\])*\s*=")
    matches = list(declaration.finditer(source))
    if not matches:
        raise ExtractionError("missing initializer " + name)
    if len(matches) != 1:
        raise ExtractionError("duplicate initializer " + name)
    match = matches[0]
    start = source.find("{", match.end())
    if start < 0:
        raise ExtractionError("missing opening brace for " + name)
    depth = 0
    quote = False
    escaped = False
    for index in range(start, len(source)):
        char = source[index]
        if quote:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                quote = False
            continue
        if char == '"':
            quote = True
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return source[start:index + 1]
            if depth < 0:
                break
    raise ExtractionError("unbalanced initializer " + name)


_TOKEN = re.compile(r"\s*(\{|\}|,|[-+]?[0-9]+|\"(?:\\.|[^\"])*\"|[A-Za-z_][A-Za-z_0-9]*)")


def _tokens(text):
    result = []
    position = 0
    while position < len(text):
        match = _TOKEN.match(text, position)
        if match is None:
            if text[position:].strip():
                raise ExtractionError("unexpected initializer text near " + repr(text[position:position + 30]))
            break
        result.append(match.group(1))
        position = match.end()
    return result


def _parse_initializer(text):
    tokens = _tokens(text)
    cursor = 0

    def parse_value():
        nonlocal cursor
        if cursor >= len(tokens):
            raise ExtractionError("truncated initializer")
        token = tokens[cursor]
        if token == "{":
            cursor += 1
            values = []
            if cursor < len(tokens) and tokens[cursor] == "}":
                cursor += 1
                return values
            while True:
                values.append(parse_value())
                if cursor >= len(tokens):
                    raise ExtractionError("unterminated initializer list")
                if tokens[cursor] == ",":
                    cursor += 1
                    if cursor < len(tokens) and tokens[cursor] == "}":
                        cursor += 1
                        return values
                    continue
                if tokens[cursor] == "}":
                    cursor += 1
                    return values
                raise ExtractionError("expected comma or closing brace")
        if token in ("{", "}", ","):
            raise ExtractionError("unexpected delimiter")
        cursor += 1
        if re.fullmatch(r"[-+]?[0-9]+", token):
            return int(token)
        if token.startswith('"'):
            try:
                value = json.loads(token)
            except json.JSONDecodeError as error:
                raise ExtractionError("invalid C string") from error
            if not isinstance(value, str):
                raise ExtractionError("initializer string is not text")
            return value
        return token

    value = parse_value()
    if cursor != len(tokens):
        raise ExtractionError("trailing initializer tokens")
    return value


def _ints(value, name):
    if not isinstance(value, list) or not all(isinstance(item, int) for item in value):
        raise ExtractionError(name + " must be an integer array")
    return value


def _normalize_rows(value, rows, columns, name):
    if not isinstance(value, list) or len(value) > rows:
        raise ExtractionError(name + " has an invalid row count")
    result = []
    for row in value:
        if not isinstance(row, list) or len(row) > columns:
            raise ExtractionError(name + " has an invalid row width")
        if not all(isinstance(item, int) for item in row):
            raise ExtractionError(name + " contains non-integer entries")
        result.append(row + [0] * (columns - len(row)))
    result.extend([[0] * columns for _ in range(rows - len(result))])
    return result


def _normalize_3d(value, rows, columns, width, name):
    if not isinstance(value, list) or len(value) > rows:
        raise ExtractionError(name + " has an invalid outer count")
    result = []
    for row in value:
        result.append(_normalize_rows(row, columns, width, name))
    result.extend([[[0] * width for _ in range(columns)] for _ in range(rows - len(result))])
    return result


def _fixed_width_digits(value, radix, width, name):
    limit = radix ** width
    if (not isinstance(value, int) or isinstance(value, bool)
            or not 0 <= value < limit):
        raise ExtractionError(f"{name} payload out of range")
    return [(value // (radix ** exponent)) % radix
            for exponent in range(width - 1, -1, -1)]


def _decode_magnetic_operation(encoded):
    if (not isinstance(encoded, int) or isinstance(encoded, bool)
            or not 0 <= encoded < MAGNETIC_OPERATION_ENCODING_LIMIT):
        raise ExtractionError("magnetic operation encoding is out of range")
    time_reversal, spatial = divmod(encoded, MSG_OPERATION_SCALE)
    if time_reversal not in (0, 1):
        raise ExtractionError("magnetic operation has invalid time-reversal bit")

    # This mirrors spgdb_decode_symmetry in the pinned C source.  Each digit
    # is selected from the original payload; consuming a quotient/remainder
    # pair on every iteration changes the positional-base interpretation.
    rotation_payload = spatial % ROTATION_PAYLOAD
    translation_payload = spatial // ROTATION_PAYLOAD
    rotation = [digit - 1 for digit in
                _fixed_width_digits(rotation_payload, ROTATION_RADIX,
                                    ROTATION_DIGITS, "rotation")]
    translation = _fixed_width_digits(translation_payload,
                                      TRANSLATION_DENOMINATOR,
                                      TRANSLATION_DIGITS, "translation")

    if any(item not in (-1, 0, 1) for item in rotation):
        raise ExtractionError("decoded rotation trit out of range")
    if any(item < 0 or item >= TRANSLATION_DENOMINATOR for item in translation):
        raise ExtractionError("decoded translation digit out of range")
    if sum((item + 1) * ROTATION_RADIX ** (ROTATION_DIGITS - 1 - i)
           for i, item in enumerate(rotation)) != rotation_payload:
        raise ExtractionError("rotation payload has a nonzero remainder")
    if sum(item * TRANSLATION_DENOMINATOR ** (TRANSLATION_DIGITS - 1 - i)
           for i, item in enumerate(translation)) != translation_payload:
        raise ExtractionError("translation payload has a nonzero remainder")
    reconstructed = (time_reversal * MSG_OPERATION_SCALE
                     + sum((item + 1) * ROTATION_RADIX ** (ROTATION_DIGITS - 1 - i)
                           for i, item in enumerate(rotation))
                     + ROTATION_PAYLOAD * sum(
                         item * TRANSLATION_DENOMINATOR ** (TRANSLATION_DIGITS - 1 - i)
                         for i, item in enumerate(translation)))
    if reconstructed != encoded:
        raise ExtractionError("magnetic operation encoding round-trip mismatch")
    return {"rotation": rotation, "translation_numerator": translation,
            "time_reversal": time_reversal}


def _source(path, expected_hash):
    path = Path(path)
    actual = _sha256(path)
    if actual != expected_hash:
        raise ExtractionError(f"{path}: SHA256 mismatch: {actual}")
    return _strip_comments(path.read_text(encoding="utf-8")), actual


def extract(upstream):
    upstream = Path(upstream)
    msg_path = upstream / "src/msg_database.c"
    spg_path = upstream / "src/spg_database.c"
    msg, msg_hash = _source(msg_path, EXPECTED_SOURCES["msg_database.c"])
    spg, spg_hash = _source(spg_path, EXPECTED_SOURCES["spg_database.c"])

    spg_types = _parse_initializer(_initializer_text(spg, "spacegroup_types"))
    spg_index = _parse_initializer(_initializer_text(spg, "symmetry_operation_index"))
    spg_operations = _ints(_parse_initializer(_initializer_text(spg, "symmetry_operations")), "symmetry_operations")
    if len(spg_types) != SPG_HALL_COUNT or len(spg_index) != SPG_HALL_COUNT:
        raise ExtractionError("spg Hall census is not 531 including sentinel")
    spg_numbers = []
    for entry in spg_types:
        if (not isinstance(entry, list) or len(entry) != 9
                or not isinstance(entry[0], int)
                or not all(isinstance(item, str) for item in entry[1:7])
                or entry[7] not in {"CENTERING_ERROR", "PRIMITIVE", "C_FACE", "A_FACE", "BODY", "FACE", "R_CENTER"}
                or not isinstance(entry[8], int)):
            raise ExtractionError("spacegroup_types entry has invalid grammar")
        spg_numbers.append(entry[0])
    spg_operation_index = []
    for entry in spg_index:
        if not isinstance(entry, list) or len(entry) != 2 or not all(isinstance(x, int) for x in entry):
            raise ExtractionError("spg operation index entry must have width 2")
        spg_operation_index.append(entry)
    if len(spg_operations) != SPG_OPERATION_COUNT or spg_operations[0] != 0:
        raise ExtractionError("spg operation census/sentinel mismatch")
    if any(not isinstance(item, int) or isinstance(item, bool)
           or not 0 < item < MSG_OPERATION_SCALE for item in spg_operations[1:]):
        raise ExtractionError("spg operation encoding out of range")
    if spg_operation_index[0] != [0, 0]:
        raise ExtractionError("spg operation index dummy mismatch")
    previous_end = 1
    for hall_number, entry in enumerate(spg_operation_index[1:], 1):
        order, offset = entry
        if order <= 0 or offset < 1 or offset + order > len(spg_operations):
            raise ExtractionError("spg operation span out of range")
        if offset + order > SPG_STANDARD_OPERATION_END or offset != previous_end:
            raise ExtractionError("spg operation spans are not contiguous")
        previous_end = offset + order
    if previous_end != SPG_STANDARD_OPERATION_END:
        raise ExtractionError("spg operation standard boundary mismatch")
    if not SPG_STANDARD_OPERATION_END < len(spg_operations):
        raise ExtractionError("spg operation layer tail is missing")

    msg_types = _parse_initializer(_initializer_text(msg, "magnetic_spacegroup_types"))
    uni_mapping = _parse_initializer(_initializer_text(msg, "magnetic_spacegroup_uni_mapping"))
    operation_index = _parse_initializer(_initializer_text(msg, "magnetic_spacegroup_operation_index"))
    magnetic_operations = _ints(_parse_initializer(_initializer_text(msg, "magnetic_symmetry_operations")), "magnetic_symmetry_operations")
    transformations = _parse_initializer(_initializer_text(msg, "alternative_transformations"))
    if not all(isinstance(row, list) and len(row) == 6
               and isinstance(row[0], int) and isinstance(row[1], int)
               and isinstance(row[2], str) and isinstance(row[3], str)
               and isinstance(row[4], int) and isinstance(row[5], int)
               for row in msg_types):
        raise ExtractionError("magnetic spacegroup type row width mismatch")
    if len(msg_types) != MSG_UNI_COUNT or len(uni_mapping) != MSG_UNI_COUNT:
        raise ExtractionError("UNI census is not 1652 including sentinel")
    if len(magnetic_operations) != MSG_OPERATION_COUNT or magnetic_operations[0] != 0:
        raise ExtractionError("magnetic operation census/sentinel mismatch")
    if any(not isinstance(item, int) or isinstance(item, bool)
           or not 0 < item < MAGNETIC_OPERATION_ENCODING_LIMIT
           for item in magnetic_operations[1:]):
        raise ExtractionError("magnetic operation encoding out of range")
    uni_mapping = [
        [int(row[0]), int(row[1])] if isinstance(row, list) and len(row) == 2 and all(isinstance(x, int) for x in row)
        else (_ for _ in ()).throw(ExtractionError("UNI mapping row width mismatch"))
        for row in uni_mapping
    ]
    operation_index = _normalize_3d(operation_index, MSG_UNI_COUNT, MSG_HALL_SLOTS, 2, "magnetic_spacegroup_operation_index")
    transformations = _normalize_3d(transformations, MSG_UNI_COUNT, MSG_HALL_SLOTS, 7, "alternative_transformations")
    if msg_types[0] != [0, 0, "", "", 0, 0]:
        raise ExtractionError("magnetic type sentinel mismatch")
    for uni, row in enumerate(msg_types[1:], 1):
        if row[0] != uni or not isinstance(row[1], int) or not isinstance(row[4], int) or row[5] not in (1, 2, 3, 4):
            raise ExtractionError("magnetic type row identity/range mismatch")
    type_counts = {str(kind): sum(row[5] == kind for row in msg_types[1:]) for kind in (1, 2, 3, 4)}
    if type_counts != EXPECTED_TYPE_COUNTS:
        raise ExtractionError(f"magnetic type census mismatch: {type_counts}")
    decoded_operations = [_decode_magnetic_operation(item) for item in magnetic_operations]

    if uni_mapping[0] != [0, 0] or operation_index[0][0] != [0, 0]:
        raise ExtractionError("magnetic database dummy mapping mismatch")
    if any(entry != [0, 0] for entry in operation_index[0][1:]):
        raise ExtractionError("magnetic database dummy operation tail mismatch")
    if any(entry != [0] * 7 for entry in transformations[0]):
        raise ExtractionError("magnetic database dummy transformation mismatch")

    active_spans = []
    for uni in range(1, MSG_UNI_COUNT):
        hall_count, first_hall = uni_mapping[uni]
        if hall_count < 1 or hall_count > MSG_HALL_SLOTS:
            raise ExtractionError("UNI Hall count out of range")
        if (first_hall < 1 or first_hall > SPG_HALL_SETTINGS
                or first_hall + hall_count - 1 > SPG_HALL_SETTINGS):
            raise ExtractionError("UNI mapping range mismatch")
        for slot, (order, offset) in enumerate(operation_index[uni]):
            if slot < hall_count:
                if (order <= 0 or offset < 1
                        or offset + order > len(magnetic_operations)):
                    raise ExtractionError("magnetic operation span out of range")
                active_spans.append((offset, offset + order))
            elif [order, offset] != [0, 0]:
                raise ExtractionError("nonzero operation index beyond Hall count")
        for slot, values in enumerate(transformations[uni]):
            if slot >= hall_count:
                if values != [0] * 7:
                    raise ExtractionError("nonzero transformation beyond Hall count")
                continue
            first_zero = next((index for index, value in enumerate(values)
                               if value == 0), len(values))
            if any(value == 0 for value in values[:first_zero]):
                raise ExtractionError("invalid transformation zero terminator")
            if any(value != 0 for value in values[first_zero:]):
                raise ExtractionError("nonzero transformation tail")
            if any(not isinstance(value, int) or isinstance(value, bool)
                   or not 0 < value < MSG_OPERATION_SCALE
                   for value in values[:first_zero]):
                raise ExtractionError("transformation encoding out of range")

    if len(active_spans) != MSG_ACTIVE_SPAN_COUNT:
        raise ExtractionError("magnetic active span census mismatch")
    previous_end = 1
    for start, end in sorted(active_spans):
        if start != previous_end or end <= start:
            raise ExtractionError("magnetic operation spans are not contiguous")
        previous_end = end
    if previous_end != len(magnetic_operations):
        raise ExtractionError("magnetic operation span boundary mismatch")

    expected_witnesses = {
        16484: {"rotation": [1, 0, 0, 0, 1, 0, 0, 0, 1],
                "translation_numerator": [0, 0, 0], "time_reversal": 0},
        34146806: {"rotation": [1, 0, 0, 0, 1, 0, 0, 0, 1],
                   "translation_numerator": [0, 0, 6], "time_reversal": 1},
        3198: {"rotation": [-1, 0, 0, 0, -1, 0, 0, 0, -1],
               "translation_numerator": [0, 0, 0], "time_reversal": 0},
        34133520: {"rotation": [-1, 0, 0, 0, -1, 0, 0, 0, -1],
                   "translation_numerator": [0, 0, 6], "time_reversal": 1},
        3360: {"rotation": [-1, 0, 0, 0, 1, 0, 0, 0, -1],
               "translation_numerator": [0, 0, 0], "time_reversal": 0},
        34028708: {"rotation": [1, 0, 0, 0, 1, 0, 0, 0, 1],
                   "translation_numerator": [0, 0, 0], "time_reversal": 1},
        34015584: {"rotation": [-1, 0, 0, 0, 1, 0, 0, 0, -1],
                   "translation_numerator": [0, 0, 0], "time_reversal": 1},
        3200: {"rotation": [-1, 0, 0, 0, -1, 0, 0, 0, 1],
               "translation_numerator": [0, 0, 0], "time_reversal": 0},
        34015424: {"rotation": [-1, 0, 0, 0, -1, 0, 0, 0, 1],
                   "translation_numerator": [0, 0, 0], "time_reversal": 1},
        16320: {"rotation": [1, 0, 0, 0, -1, 0, 0, 0, -1],
                "translation_numerator": [0, 0, 0], "time_reversal": 0},
        34028544: {"rotation": [1, 0, 0, 0, -1, 0, 0, 0, -1],
                   "translation_numerator": [0, 0, 0], "time_reversal": 1},
    }
    for encoded, expected in expected_witnesses.items():
        if _decode_magnetic_operation(encoded) != expected:
            raise ExtractionError(f"magnetic decoder witness mismatch: {encoded}")

    artifact = {
        "schema": SCHEMA,
        "translation_denominator": TRANSLATION_DENOMINATOR,
        "spg": {
            "spacegroup_number": spg_numbers,
            "symmetry_operation_index": spg_operation_index,
            "symmetry_operations": spg_operations,
        },
        "msg": {
            "magnetic_spacegroup_types": [
                {"uni": row[0], "litvin": row[1], "bns": row[2], "og": row[3],
                 "parent_spacegroup": row[4], "type": row[5]}
                for row in msg_types
            ],
            "magnetic_spacegroup_uni_mapping": uni_mapping,
            "magnetic_spacegroup_operation_index": operation_index,
            "magnetic_symmetry_operations": magnetic_operations,
            "alternative_transformations": transformations,
        },
    }
    return artifact, {
        "msg_database.c": msg_hash,
        "spg_database.c": spg_hash,
        "type_counts": type_counts,
        "decoded_magnetic_operations": decoded_operations,
    }


def canonical_json(value):
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def write_outputs(upstream, output, manifest):
    artifact, details = extract(upstream)
    output = Path(output)
    manifest = Path(manifest)
    artifact_bytes = canonical_json(artifact)
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(artifact_bytes)
    manifest_value = {
        "schema": MANIFEST_SCHEMA,
        "repository": "https://github.com/spglib/spglib",
        "tag": "v2.5.0",
        "commit": "e4531bb49371dce3e807c2095a4d9d9b7245c524",
        "sources": [
            {"path": "src/msg_database.c", "sha256": details["msg_database.c"]},
            {"path": "src/spg_database.c", "sha256": details["spg_database.c"]},
        ],
        "extractor_schema_version": EXTRACTOR_VERSION,
        "artifact": {"path": str(output.name), "bytes": len(artifact_bytes),
                      "sha256": hashlib.sha256(artifact_bytes).hexdigest()},
    }
    manifest.write_bytes(canonical_json(manifest_value))
    return artifact, details, manifest_value


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args(argv)
    write_outputs(args.upstream, args.output, args.manifest)


if __name__ == "__main__":
    main()
