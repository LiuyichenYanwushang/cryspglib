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
MSG_OPERATION_SCALE = 34_012_224
MSG_OPERATION_COUNT = 76_683
SPG_HALL_COUNT = 531
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


def _decode_magnetic_operation(encoded):
    if not isinstance(encoded, int) or encoded < 0:
        raise ExtractionError("magnetic operation encoding must be nonnegative")
    timerev, payload = divmod(encoded, MSG_OPERATION_SCALE)
    if timerev not in (0, 1):
        raise ExtractionError("magnetic operation has invalid time-reversal bit")
    rotation_payload, translation_payload = payload % 19_683, payload // 19_683
    rotation = []
    power = 6_561
    for _ in range(9):
        digit, rotation_payload = divmod(rotation_payload, power * 3)
        rotation.append(digit - 1)
        power //= 3
    translation = []
    power = 144
    for _ in range(3):
        digit, translation_payload = divmod(translation_payload, power * 12)
        translation.append(digit)
        power //= 12
    if any(item not in (-1, 0, 1) for item in rotation):
        raise ExtractionError("decoded rotation trit out of range")
    if any(item < 0 or item >= TRANSLATION_DENOMINATOR for item in translation):
        raise ExtractionError("decoded translation digit out of range")
    return {"rotation": rotation, "translation_numerator": translation, "time_reversal": timerev}


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
    if len(spg_operations) <= 0 or spg_operations[0] != 0:
        raise ExtractionError("spg operation encoding sentinel missing")
    for start, count in spg_operation_index:
        if start < 0 or count < 0 or start + count > len(spg_operations):
            raise ExtractionError("spg operation span out of range")

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
    for uni in range(1, MSG_UNI_COUNT):
        hall_count, first_hall = uni_mapping[uni]
        if hall_count < 1 or first_hall < 1 or first_hall > SPG_HALL_COUNT - 1:
            raise ExtractionError("UNI mapping range mismatch")
        for slot, (order, offset) in enumerate(operation_index[uni]):
            if slot < hall_count:
                if order <= 0 or offset <= 0 or offset + order > len(magnetic_operations):
                    raise ExtractionError("magnetic operation span out of range")
            elif [order, offset] != [0, 0]:
                raise ExtractionError("nonzero operation index beyond Hall count")
        for slot, values in enumerate(transformations[uni]):
            if slot >= hall_count and values != [0] * 7:
                raise ExtractionError("nonzero transformation beyond Hall count")

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
