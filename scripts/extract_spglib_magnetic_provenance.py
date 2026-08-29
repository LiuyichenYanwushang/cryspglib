#!/usr/bin/env python3
"""Extract pinned spglib magnetic database tables into a raw JSON artifact.

This is deliberately a small, strict C-initializer reader.  It consumes only
the two pinned upstream C files and never imports the generated Rust tables.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import os
import re
import subprocess
import tempfile
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
# Names matching the upstream decoder's terminology are kept alongside the
# descriptive names above for callers that want to state the C contract.
ROT_RADIX = ROTATION_RADIX
ROT_WIDTH = ROTATION_DIGITS
ROT_SCALE = ROTATION_PAYLOAD
TRANS_RADIX = TRANSLATION_DENOMINATOR
TRANS_WIDTH = TRANSLATION_DIGITS
TRANS_SCALE = TRANSLATION_PAYLOAD
SPACE_OPERATION_SCALE = MSG_OPERATION_SCALE
MAGNETIC_OPERATION_ENCODING_LIMIT = 2 * MSG_OPERATION_SCALE
SPG_OPERATION_COUNT = 8_147
SPG_STANDARD_OPERATION_END = 7_389
MSG_OPERATION_COUNT = 76_683
MSG_ACTIVE_SPAN_COUNT = 4_479
ALTERNATIVE_TRANSFORMATION_VALUE_COUNT = 536
SPG_HALL_COUNT = 531
SPG_HALL_SETTINGS = SPG_HALL_COUNT - 1
MSG_UNI_COUNT = 1_652
MSG_HALL_SLOTS = 18

EXPECTED_SOURCES = {
    "msg_database.c": "2c45cf2e3827b48bab149ca59951f6b242f3754cc7d559585dfe107aa2a94288",
    "spg_database.c": "1ad4d3fd4ee2b39d43bf102bd6464bc6d1e07bc825dc4ac40df0235501823896",
}
EXPECTED_TYPE_COUNTS = {"1": 230, "2": 230, "3": 674, "4": 517}
EXPECTED_UPSTREAM_TAG = "v2.5.0"
EXPECTED_UPSTREAM_COMMIT = "e4531bb49371dce3e807c2095a4d9d9b7245c524"

PINNED_DECLARATIONS = {
    "spacegroup_types": (
        r"^[ \t]*static[ \t]+SpacegroupType[ \t]+const[ \t]+"
        r"spacegroup_types[ \t]*\[[ \t]*\][ \t]*="
    ),
    "symmetry_operation_index": (
        r"^[ \t]*static[ \t]+int[ \t]+const[ \t]+"
        r"symmetry_operation_index[ \t]*\[[ \t]*\][ \t]*\[[ \t]*2[ \t]*\][ \t]*="
    ),
    "symmetry_operations": (
        r"^[ \t]*static[ \t]+int[ \t]+const[ \t]+"
        r"symmetry_operations[ \t]*\[[ \t]*\][ \t]*="
    ),
    "magnetic_spacegroup_types": (
        r"^[ \t]*static[ \t]+const[ \t]+MagneticSpacegroupType[ \t]+"
        r"magnetic_spacegroup_types[ \t]*\[[ \t]*\][ \t]*="
    ),
    "magnetic_spacegroup_uni_mapping": (
        r"^[ \t]*static[ \t]+const[ \t]+int[ \t]+"
        r"magnetic_spacegroup_uni_mapping[ \t]*\[[ \t]*\][ \t]*\[[ \t]*2[ \t]*\][ \t]*="
    ),
    "magnetic_spacegroup_hall_mapping": (
        r"^[ \t]*static[ \t]+const[ \t]+int[ \t]+"
        r"magnetic_spacegroup_hall_mapping[ \t]*\[[ \t]*\][ \t]*\[[ \t]*2[ \t]*\][ \t]*="
    ),
    "magnetic_spacegroup_operation_index": (
        r"^[ \t]*static[ \t]+const[ \t]+int[ \t]+"
        r"magnetic_spacegroup_operation_index[ \t]*\[[ \t]*\][ \t]*\[[ \t]*18[ \t]*\]"
        r"[ \t]*\[[ \t]*2[ \t]*\][ \t]*="
    ),
    "magnetic_symmetry_operations": (
        r"^[ \t]*static[ \t]+const[ \t]+int[ \t]+"
        r"magnetic_symmetry_operations[ \t]*\[[ \t]*\][ \t]*="
    ),
    "alternative_transformations": (
        r"^[ \t]*static[ \t]+const[ \t]+int[ \t]+"
        r"alternative_transformations[ \t]*\[[ \t]*\][ \t]*\[[ \t]*18[ \t]*\]"
        r"[ \t]*\[[ \t]*7[ \t]*\][ \t]*="
    ),
}

MAGNETIC_DECODER_WITNESSES = {
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
UNI7_RAW_WITNESS = [16484, 3198, 34146806, 34133520]
UNI9_RAW_WITNESSES = [
    [16484, 3360, 34028708, 34015584],
    [16484, 3200, 34028708, 34015424],
    [16484, 16320, 34028708, 34028544],
]


class ExtractionError(ValueError):
    """Raised for any malformed or unexpected upstream initializer."""


class CString(str):
    """A decoded C string literal, distinct from a bare identifier."""


class Identifier(str):
    """A bare C identifier used only by the pinned enum-like fields."""


@dataclass(frozen=True)
class _Token:
    kind: str
    value: str


def _strip_comments(text):
    """Remove C comments without changing token boundaries or strings."""
    result = []
    state = "normal"
    index = 0
    while index < len(text):
        char = text[index]
        if char == "\\" and index + 1 < len(text) and text[index + 1] in "\r\n":
            raise ExtractionError("backslash-newline line splice is not supported")
        if state == "normal":
            if char == '"':
                result.append(char)
                state = "string"
                index += 1
            elif char == "'":
                result.append(char)
                state = "char"
                index += 1
            elif char == "/" and index + 1 < len(text) and text[index + 1] == "/":
                result.append(" ")
                state = "line-comment"
                index += 2
            elif char == "/" and index + 1 < len(text) and text[index + 1] == "*":
                result.append(" ")
                state = "block-comment"
                index += 2
            else:
                result.append(char)
                index += 1
        elif state == "string":
            result.append(char)
            if char == "\\":
                if index + 1 >= len(text):
                    raise ExtractionError("unterminated C string escape")
                result.append(text[index + 1])
                index += 2
            elif char == '"':
                state = "normal"
                index += 1
            elif char in "\r\n":
                raise ExtractionError("newline in C string literal")
            else:
                index += 1
        elif state == "char":
            result.append(char)
            if char == "\\":
                if index + 1 >= len(text):
                    raise ExtractionError("unterminated C character escape")
                result.append(text[index + 1])
                index += 2
            elif char == "'":
                state = "normal"
                index += 1
            elif char in "\r\n":
                raise ExtractionError("newline in C character literal")
            else:
                index += 1
        elif state == "line-comment":
            if char in "\r\n":
                result.append(char)
                state = "normal"
            index += 1
        else:
            if char == "*" and index + 1 < len(text) and text[index + 1] == "/":
                result.append(" ")
                state = "normal"
                index += 2
            else:
                if char in "\r\n":
                    result.append(char)
                index += 1
    if state == "block-comment":
        raise ExtractionError("unterminated C block comment")
    if state == "string":
        raise ExtractionError("unterminated C string literal")
    if state == "char":
        raise ExtractionError("unterminated C character literal")
    return "".join(result)


def _initializer_text(source, name):
    source = _strip_comments(source)
    if name not in PINNED_DECLARATIONS:
        raise ExtractionError("unrecognized pinned initializer " + name)
    declaration = re.compile(PINNED_DECLARATIONS[name], flags=re.MULTILINE)
    matches = []
    for candidate in declaration.finditer(source):
        before = source[:candidate.start()]
        escaped = False
        quote = None
        for char in before:
            if quote is not None:
                if escaped:
                    escaped = False
                elif char == "\\":
                    escaped = True
                elif char == quote:
                    quote = None
            elif char in ('"', "'"):
                quote = char
        if quote is None:
            matches.append(candidate)
    if not matches:
        raise ExtractionError("missing initializer " + name)
    if len(matches) != 1:
        raise ExtractionError("duplicate initializer " + name)
    match = matches[0]
    opening = re.match(r"\s*(\{)", source[match.end():])
    if opening is None:
        raise ExtractionError("missing opening brace for " + name)
    start = match.end() + opening.start(1)
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
                trailer = source[index + 1:]
                if re.match(r"\s*;", trailer) is None:
                    raise ExtractionError("trailing initializer tokens")
                return source[start:index + 1]
            if depth < 0:
                break
    raise ExtractionError("unbalanced initializer " + name)


def _tokens(text):
    text = _strip_comments(text)
    result = []
    position = 0
    while position < len(text):
        char = text[position]
        if char.isspace():
            position += 1
            continue
        if char in "{},":
            result.append(_Token({"{": "lbrace", "}": "rbrace", ",": "comma"}[char], char))
            position += 1
            continue
        if char == "-":
            result.append(_Token("minus", char))
            position += 1
            continue
        if char == "+":
            raise ExtractionError("explicit plus is not allowed in pinned C integers")
        if char == '"':
            start = position
            position += 1
            escaped = False
            while position < len(text):
                current = text[position]
                if current in "\r\n":
                    raise ExtractionError("newline in C string literal")
                if escaped:
                    escaped = False
                elif current == "\\":
                    escaped = True
                elif current == '"':
                    position += 1
                    break
                position += 1
            else:
                raise ExtractionError("unterminated C string literal")
            if escaped:
                raise ExtractionError("unterminated C string escape")
            result.append(_Token("cstring", text[start:position]))
            continue
        if char in "0123456789":
            start = position
            while position < len(text) and text[position] in "0123456789":
                position += 1
            spelling = text[start:position]
            if len(spelling) > 1 and spelling.startswith("0"):
                raise ExtractionError("non-canonical or invalid-octal integer: " + spelling)
            if (position < len(text)
                    and (("A" <= text[position] <= "Z")
                         or ("a" <= text[position] <= "z")
                         or text[position] in "._")):
                raise ExtractionError("malformed integer spelling near " + spelling)
            result.append(_Token("integer", spelling))
            continue
        if (("A" <= char <= "Z") or ("a" <= char <= "z") or char == "_"):
            start = position
            position += 1
            while position < len(text) and (
                    ("A" <= text[position] <= "Z")
                    or ("a" <= text[position] <= "z")
                    or text[position] in "0123456789_"):
                position += 1
            result.append(_Token("identifier", text[start:position]))
            continue
        raise ExtractionError("unexpected initializer text near " + repr(text[position:position + 30]))
    return result


def _c_escape_character(codepoint, token, escape_position):
    if (not isinstance(codepoint, int) or not 0 <= codepoint <= 0x10FFFF
            or 0xD800 <= codepoint <= 0xDFFF):
        raise ExtractionError(
            f"invalid Unicode scalar in C string {token!r} at offset "
            f"{escape_position}: U+{codepoint:X}"
        )
    try:
        return chr(codepoint)
    except (ValueError, UnicodeError) as error:
        raise ExtractionError(
            f"invalid C string escape in {token!r} at offset {escape_position}"
        ) from error


def _decode_c_string(token):
    if not (isinstance(token, str) and len(token) >= 2
            and token[0] == '"' and token[-1] == '"'):
        raise ExtractionError("invalid C string")
    value = []
    position = 1
    end = len(token) - 1
    simple_escapes = {
        "a": "\a", "b": "\b", "f": "\f", "n": "\n",
        "r": "\r", "t": "\t", "v": "\v", "\\": "\\",
        '"': '"', "?": "?", "'": "'",
    }
    while position < end:
        char = token[position]
        if char != "\\":
            value.append(char)
            position += 1
            continue
        escape_position = position
        position += 1
        if position >= end:
            raise ExtractionError("unterminated C string escape")
        escaped = token[position]
        if escaped in simple_escapes:
            value.append(simple_escapes[escaped])
            position += 1
            continue
        if escaped in "01234567":
            start = position
            position += 1
            while position < end and position - start < 3 and token[position] in "01234567":
                position += 1
            try:
                codepoint = int(token[start:position], 8)
            except ValueError as error:
                raise ExtractionError(
                    f"invalid octal escape in {token!r} at offset {escape_position}"
                ) from error
            value.append(_c_escape_character(codepoint, token, escape_position))
            continue
        if escaped == "x":
            position += 1
            start = position
            while position < end and token[position] in "0123456789abcdefABCDEF":
                position += 1
            if start == position:
                raise ExtractionError("C hexadecimal escape requires digits")
            try:
                codepoint = int(token[start:position], 16)
            except ValueError as error:
                raise ExtractionError(
                    f"invalid hexadecimal escape in {token!r} at offset {escape_position}"
                ) from error
            value.append(_c_escape_character(codepoint, token, escape_position))
            continue
        raise ExtractionError("unsupported/non-C string escape")
    decoded = "".join(value)
    try:
        decoded.encode("utf-8")
    except UnicodeEncodeError as error:
        raise ExtractionError(
            f"C string {token!r} is not valid UTF-8"
        ) from error
    return CString(decoded)


def _parse_initializer(text):
    tokens = _tokens(text)
    cursor = 0

    def parse_value():
        nonlocal cursor
        if cursor >= len(tokens):
            raise ExtractionError("truncated initializer")
        token = tokens[cursor]
        if token.kind == "lbrace":
            cursor += 1
            values = []
            if cursor < len(tokens) and tokens[cursor].kind == "rbrace":
                cursor += 1
                return values
            while True:
                values.append(parse_value())
                if cursor >= len(tokens):
                    raise ExtractionError("unterminated initializer list")
                if tokens[cursor].kind == "comma":
                    cursor += 1
                    if cursor < len(tokens) and tokens[cursor].kind == "rbrace":
                        cursor += 1
                        return values
                    continue
                if tokens[cursor].kind == "rbrace":
                    cursor += 1
                    return values
                raise ExtractionError("expected comma or closing brace")
        if token.kind in ("lbrace", "rbrace", "comma"):
            raise ExtractionError("unexpected delimiter")
        if token.kind == "minus":
            cursor += 1
            if cursor >= len(tokens) or tokens[cursor].kind != "integer":
                raise ExtractionError("minus must precede an integer")
            token_value = tokens[cursor].value
            cursor += 1
            return -int(token_value)
        cursor += 1
        if token.kind == "integer":
            return int(token.value)
        if token.kind == "cstring":
            return _decode_c_string(token.value)
        if token.kind == "identifier":
            return Identifier(token.value)
        raise ExtractionError("unknown initializer token kind")

    value = parse_value()
    if cursor != len(tokens):
        raise ExtractionError("trailing initializer tokens")
    return value


def _ints(value, name):
    if not isinstance(value, list) or not all(type(item) is int for item in value):
        raise ExtractionError(name + " must be an integer array")
    return value


def _normalize_rows(value, rows, columns, name):
    if not isinstance(value, list) or len(value) > rows:
        raise ExtractionError(name + " has an invalid row count")
    result = []
    for row in value:
        if not isinstance(row, list) or len(row) > columns:
            raise ExtractionError(name + " has an invalid row width")
        if not all(type(item) is int for item in row):
            raise ExtractionError(name + " contains non-integer entries")
        result.append(row + [0] * (columns - len(row)))
    result.extend([[0] * columns for _ in range(rows - len(result))])
    return result


def _normalize_3d(value, rows, columns, width, name):
    if not isinstance(value, list) or len(value) != rows:
        raise ExtractionError(name + " has an invalid outer count")
    result = []
    for row in value:
        result.append(_normalize_rows(row, columns, width, name))
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


def _validate_hall_mapping(hall_mapping, uni_mapping, msg_types, spg_numbers):
    """Check the bidirectional Hall↔UNI index and parent SG relation."""
    if not isinstance(hall_mapping, list) or len(hall_mapping) != SPG_HALL_COUNT:
        raise ExtractionError("magnetic Hall mapping census mismatch")
    for hall, pair in enumerate(hall_mapping):
        if (not isinstance(pair, list) or len(pair) != 2
                or not all(type(value) is int for value in pair)):
            raise ExtractionError(f"magnetic Hall mapping row {hall} is invalid")
    if hall_mapping[0] != [0, 0] or uni_mapping[0] != [0, 0]:
        raise ExtractionError("magnetic Hall mapping dummy mismatch")
    if uni_mapping[1] != [1, 1]:
        raise ExtractionError("UNI1 mapping witness mismatch")
    for hall, (smallest_uni, largest_uni) in enumerate(hall_mapping[1:], 1):
        if not 1 <= smallest_uni <= largest_uni < MSG_UNI_COUNT:
            raise ExtractionError("magnetic Hall mapping UNI range mismatch")
        expected_unis = []
        for uni in range(1, MSG_UNI_COUNT):
            hall_count, first_hall = uni_mapping[uni]
            if first_hall <= hall < first_hall + hall_count:
                expected_unis.append(uni)
        if expected_unis != list(range(smallest_uni, largest_uni + 1)):
            raise ExtractionError("magnetic Hall mapping is not bidirectional")
    for uni in range(1, MSG_UNI_COUNT):
        hall_count, first_hall = uni_mapping[uni]
        parent_spacegroup = msg_types[uni][4]
        for hall in range(first_hall, first_hall + hall_count):
            smallest_uni, largest_uni = hall_mapping[hall]
            if not smallest_uni <= uni <= largest_uni:
                raise ExtractionError("UNI Hall mapping inverse mismatch")
            if spg_numbers[hall] != parent_spacegroup:
                raise ExtractionError("UNI parent spacegroup mismatch")


def _git_output(upstream, arguments):
    try:
        result = subprocess.run(
            ["git", "-C", str(upstream), *arguments],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="strict",
        )
    except (OSError, UnicodeError, subprocess.CalledProcessError) as error:
        raise ExtractionError(f"unable to verify git provenance for {upstream}") from error
    return result.stdout.strip()


def _verify_upstream_provenance(upstream):
    upstream = Path(upstream)
    if _git_output(upstream, ["rev-parse", "--is-inside-work-tree"]) != "true":
        raise ExtractionError(f"{upstream}: not a git work tree")
    head = _git_output(upstream, ["rev-parse", "HEAD"])
    if head != EXPECTED_UPSTREAM_COMMIT:
        raise ExtractionError(
            f"{upstream}: git HEAD {head!r} is not {EXPECTED_UPSTREAM_COMMIT}"
        )
    tags = _git_output(upstream, ["tag", "--points-at", "HEAD"]).splitlines()
    if EXPECTED_UPSTREAM_TAG not in tags:
        raise ExtractionError(
            f"{upstream}: HEAD is not tagged exactly {EXPECTED_UPSTREAM_TAG}"
        )


def _source_from_bytes(label, source_bytes, expected_hash):
    if type(source_bytes) is not bytes:
        raise ExtractionError(f"{label}: git blob did not return bytes")
    actual = hashlib.sha256(source_bytes).hexdigest()
    if actual != expected_hash:
        raise ExtractionError(f"{label}: SHA256 mismatch: {actual}")
    try:
        source_text = source_bytes.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise ExtractionError(f"{label}: source is not strict UTF-8") from error
    return _strip_comments(source_text), actual


def _source(path, expected_hash):
    path = Path(path)
    try:
        source_bytes = path.read_bytes()
    except OSError as error:
        raise ExtractionError(f"unable to read source {path}") from error
    return _source_from_bytes(path, source_bytes, expected_hash)


def _git_blob(upstream, relative_path):
    object_name = f"{EXPECTED_UPSTREAM_COMMIT}:{relative_path}"
    try:
        result = subprocess.run(
            ["git", "-C", str(upstream), "cat-file", "blob", object_name],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError as error:
        raise ExtractionError(f"unable to run git cat-file for {relative_path}") from error
    except subprocess.CalledProcessError as error:
        stderr = error.stderr
        if isinstance(stderr, bytes):
            stderr = stderr.decode("utf-8", errors="replace").strip()
        raise ExtractionError(
            f"git cat-file failed for {relative_path}: {stderr or 'unknown error'}"
        ) from error
    if type(result.stdout) is not bytes:
        raise ExtractionError(f"git cat-file returned non-binary output for {relative_path}")
    return result.stdout


def _git_source(upstream, relative_path, expected_hash):
    source_bytes = _git_blob(upstream, relative_path)
    return _source_from_bytes(f"{upstream}:{relative_path}", source_bytes, expected_hash)


def extract(upstream):
    upstream = Path(upstream)
    _verify_upstream_provenance(upstream)
    msg, msg_hash = _git_source(
        upstream, "src/msg_database.c", EXPECTED_SOURCES["msg_database.c"]
    )
    spg, spg_hash = _git_source(
        upstream, "src/spg_database.c", EXPECTED_SOURCES["spg_database.c"]
    )

    spg_types = _parse_initializer(_initializer_text(spg, "spacegroup_types"))
    spg_index = _parse_initializer(_initializer_text(spg, "symmetry_operation_index"))
    spg_operations = _ints(_parse_initializer(_initializer_text(spg, "symmetry_operations")), "symmetry_operations")
    if len(spg_types) != SPG_HALL_COUNT or len(spg_index) != SPG_HALL_COUNT:
        raise ExtractionError("spg Hall census is not 531 including sentinel")
    spg_numbers = []
    for entry in spg_types:
        if (not isinstance(entry, list) or len(entry) != 9
                or type(entry[0]) is not int
                or not all(type(item) is CString for item in entry[1:7])
                or type(entry[7]) is not Identifier
                or entry[7] not in {"CENTERING_ERROR", "PRIMITIVE", "C_FACE", "A_FACE", "BODY", "FACE", "R_CENTER"}
                or type(entry[8]) is not int):
            raise ExtractionError("spacegroup_types entry has invalid grammar")
        spg_numbers.append(entry[0])
    if (spg_types[0][0] != 0 or spg_types[0][7] != "CENTERING_ERROR"
            or spg_types[0][8] != 0):
        raise ExtractionError("spacegroup type sentinel mismatch")
    spg_operation_index = []
    for entry in spg_index:
        if (not isinstance(entry, list) or len(entry) != 2
                or not all(type(x) is int for x in entry)):
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
    hall_mapping = _parse_initializer(_initializer_text(msg, "magnetic_spacegroup_hall_mapping"))
    operation_index = _parse_initializer(_initializer_text(msg, "magnetic_spacegroup_operation_index"))
    magnetic_operations = _ints(_parse_initializer(_initializer_text(msg, "magnetic_symmetry_operations")), "magnetic_symmetry_operations")
    transformations = _parse_initializer(_initializer_text(msg, "alternative_transformations"))
    if not all(isinstance(row, list) and len(row) == 6
               and type(row[0]) is int and type(row[1]) is int
               and type(row[2]) is CString and type(row[3]) is CString
               and type(row[4]) is int and type(row[5]) is int
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
        [row[0], row[1]] if (isinstance(row, list) and len(row) == 2
                              and all(type(x) is int for x in row))
        else (_ for _ in ()).throw(ExtractionError("UNI mapping row width mismatch"))
        for row in uni_mapping
    ]
    operation_index = _normalize_3d(operation_index, MSG_UNI_COUNT, MSG_HALL_SLOTS, 2, "magnetic_spacegroup_operation_index")
    transformations = _normalize_3d(transformations, MSG_UNI_COUNT, MSG_HALL_SLOTS, 7, "alternative_transformations")
    if msg_types[0] != [0, 0, "", "", 0, 0]:
        raise ExtractionError("magnetic type sentinel mismatch")
    for uni, row in enumerate(msg_types[1:], 1):
        if (row[0] != uni or type(row[1]) is not int
                or type(row[4]) is not int or row[5] not in (1, 2, 3, 4)):
            raise ExtractionError("magnetic type row identity/range mismatch")
    type_counts = {str(kind): sum(row[5] == kind for row in msg_types[1:]) for kind in (1, 2, 3, 4)}
    if type_counts != EXPECTED_TYPE_COUNTS:
        raise ExtractionError(f"magnetic type census mismatch: {type_counts}")
    _validate_hall_mapping(hall_mapping, uni_mapping, msg_types, spg_numbers)
    decoded_operations = [_decode_magnetic_operation(item) for item in magnetic_operations]

    if uni_mapping[0] != [0, 0] or operation_index[0][0] != [0, 0]:
        raise ExtractionError("magnetic database dummy mapping mismatch")
    if any(entry != [0, 0] for entry in operation_index[0][1:]):
        raise ExtractionError("magnetic database dummy operation tail mismatch")
    if any(entry != [0] * 7 for entry in transformations[0]):
        raise ExtractionError("magnetic database dummy transformation mismatch")

    active_spans = []
    transformation_value_count = 0
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
            if first_zero == len(values):
                raise ExtractionError("alternative transformation terminator missing")
            if any(value != 0 for value in values[first_zero:]):
                raise ExtractionError("nonzero transformation tail")
            if any(not isinstance(value, int) or isinstance(value, bool)
                   or not 0 < value < MSG_OPERATION_SCALE
                   for value in values[:first_zero]):
                raise ExtractionError("transformation encoding out of range")
            transformation_value_count += first_zero

    if len(active_spans) != MSG_ACTIVE_SPAN_COUNT:
        raise ExtractionError("magnetic active span census mismatch")
    previous_end = 1
    for start, end in sorted(active_spans):
        if start != previous_end or end <= start:
            raise ExtractionError("magnetic operation spans are not contiguous")
        previous_end = end
    if previous_end != len(magnetic_operations):
        raise ExtractionError("magnetic operation span boundary mismatch")
    if transformation_value_count != ALTERNATIVE_TRANSFORMATION_VALUE_COUNT:
        raise ExtractionError("alternative transformation value census mismatch")

    for encoded, expected in MAGNETIC_DECODER_WITNESSES.items():
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
                {"uni": row[0], "litvin": row[1], "bns": str(row[2]), "og": str(row[3]),
                 "parent_spacegroup": row[4], "type": row[5]}
                for row in msg_types
            ],
            "magnetic_spacegroup_uni_mapping": uni_mapping,
            "magnetic_spacegroup_operation_index": operation_index,
            "magnetic_symmetry_operations": magnetic_operations,
            "alternative_transformations": transformations,
        },
    }
    validate_artifact(artifact)
    return artifact, {
        "msg_database.c": msg_hash,
        "spg_database.c": spg_hash,
        "type_counts": type_counts,
        "decoded_magnetic_operations": decoded_operations,
        "magnetic_spacegroup_hall_mapping": hall_mapping,
    }


def _reject_json_float(value):
    raise ExtractionError("JSON floating-point values are not allowed")


def _reject_json_constant(value):
    raise ExtractionError("JSON non-finite numeric constant is not allowed")


def _validate_json_value(value, path="$" ):
    if value is None or isinstance(value, (bool, int, str)):
        return
    if isinstance(value, float):
        raise ExtractionError(f"JSON floating-point value at {path}")
    if isinstance(value, list):
        for index, item in enumerate(value):
            _validate_json_value(item, f"{path}[{index}]")
        return
    if isinstance(value, dict):
        for key, item in value.items():
            if not isinstance(key, str):
                raise ExtractionError(f"JSON object key at {path} is not text")
            _validate_json_value(item, f"{path}.{key}")
        return
    raise ExtractionError(f"unsupported JSON value at {path}")


def _validate_artifact_tree(value, path="$" ):
    """Reject non-schema JSON types, including bool masquerading as int."""
    if type(value) is int or type(value) is str:
        return
    if type(value) is list:
        for index, item in enumerate(value):
            _validate_artifact_tree(item, f"{path}[{index}]")
        return
    if type(value) is dict:
        for key, item in value.items():
            if type(key) is not str:
                raise ExtractionError(f"artifact key at {path} is not text")
            _validate_artifact_tree(item, f"{path}.{key}")
        return
    raise ExtractionError(f"unsupported artifact value at {path}")


def _artifact_keys(value, expected, name):
    if type(value) is not dict or set(value) != set(expected):
        raise ExtractionError(f"{name} keys mismatch")


def _artifact_list(value, length, name):
    if type(value) is not list or len(value) != length:
        raise ExtractionError(f"{name} length mismatch")


def _artifact_int(value, name):
    if type(value) is not int:
        raise ExtractionError(f"{name} must be an integer")


def _artifact_int_pair(value, name):
    if type(value) is not list or len(value) != 2:
        raise ExtractionError(f"{name} must be a pair")
    if type(value[0]) is not int or type(value[1]) is not int:
        raise ExtractionError(f"{name} must contain integers")


def _validate_artifact_spg(spg):
    _artifact_keys(spg, {"spacegroup_number", "symmetry_operation_index",
                         "symmetry_operations"}, "artifact spg")
    numbers = spg["spacegroup_number"]
    _artifact_list(numbers, SPG_HALL_COUNT, "spacegroup_number")
    if any(type(value) is not int or not 0 <= value <= 230 for value in numbers):
        raise ExtractionError("spacegroup number range/type mismatch")
    if numbers[0] != 0:
        raise ExtractionError("spacegroup number sentinel mismatch")

    operation_index = spg["symmetry_operation_index"]
    _artifact_list(operation_index, SPG_HALL_COUNT, "symmetry_operation_index")
    for hall_number, pair in enumerate(operation_index):
        _artifact_int_pair(pair, f"symmetry_operation_index[{hall_number}]")
    if operation_index[0] != [0, 0]:
        raise ExtractionError("spg operation index dummy mismatch")

    operations = spg["symmetry_operations"]
    _artifact_list(operations, SPG_OPERATION_COUNT, "symmetry_operations")
    if operations[0] != 0:
        raise ExtractionError("spg operation sentinel mismatch")
    if any(type(value) is not int or not 0 < value < MSG_OPERATION_SCALE
           for value in operations[1:]):
        raise ExtractionError("spg operation encoding out of range")
    previous_end = 1
    for hall_number, pair in enumerate(operation_index[1:], 1):
        order, offset = pair
        if order <= 0 or offset < 1 or offset + order > len(operations):
            raise ExtractionError("spg operation span out of range")
        if offset + order > SPG_STANDARD_OPERATION_END or offset != previous_end:
            raise ExtractionError("spg operation spans are not contiguous")
        previous_end = offset + order
    if previous_end != SPG_STANDARD_OPERATION_END:
        raise ExtractionError("spg operation standard boundary mismatch")
    for value in operations[1:]:
        if _decode_magnetic_operation(value)["time_reversal"] != 0:
            raise ExtractionError("spg operation has time reversal")


def _validate_artifact_msg(msg):
    _artifact_keys(
        msg,
        {"magnetic_spacegroup_types", "magnetic_spacegroup_uni_mapping",
         "magnetic_spacegroup_operation_index", "magnetic_symmetry_operations",
         "alternative_transformations"},
        "artifact msg",
    )
    types = msg["magnetic_spacegroup_types"]
    _artifact_list(types, MSG_UNI_COUNT, "magnetic_spacegroup_types")
    type_keys = {"uni", "litvin", "bns", "og", "parent_spacegroup", "type"}
    for uni, row in enumerate(types):
        _artifact_keys(row, type_keys, f"magnetic_spacegroup_types[{uni}]")
        for field in ("uni", "litvin", "parent_spacegroup", "type"):
            _artifact_int(row[field], f"magnetic_spacegroup_types[{uni}].{field}")
        if type(row["bns"]) is not str or type(row["og"]) is not str:
            raise ExtractionError("magnetic type labels must be strings")
        if uni == 0:
            if row != {"uni": 0, "litvin": 0, "bns": "", "og": "",
                       "parent_spacegroup": 0, "type": 0}:
                raise ExtractionError("magnetic type sentinel mismatch")
        elif (row["uni"] != uni or not 1 <= row["litvin"] <= MSG_UNI_COUNT - 1
              or not 1 <= row["parent_spacegroup"] <= 230
              or row["type"] not in (1, 2, 3, 4)):
            raise ExtractionError("magnetic type identity/range mismatch")
    type_counts = {
        str(kind): sum(row["type"] == kind for row in types[1:])
        for kind in (1, 2, 3, 4)
    }
    if type_counts != EXPECTED_TYPE_COUNTS:
        raise ExtractionError(f"magnetic type census mismatch: {type_counts}")

    mapping = msg["magnetic_spacegroup_uni_mapping"]
    _artifact_list(mapping, MSG_UNI_COUNT, "magnetic_spacegroup_uni_mapping")
    for uni, pair in enumerate(mapping):
        _artifact_int_pair(pair, f"magnetic_spacegroup_uni_mapping[{uni}]")
    if mapping[0] != [0, 0]:
        raise ExtractionError("magnetic mapping dummy mismatch")

    operation_index = msg["magnetic_spacegroup_operation_index"]
    _artifact_list(operation_index, MSG_UNI_COUNT,
                   "magnetic_spacegroup_operation_index")
    for uni, row in enumerate(operation_index):
        _artifact_list(row, MSG_HALL_SLOTS,
                       f"magnetic_spacegroup_operation_index[{uni}]")
        for slot, pair in enumerate(row):
            _artifact_int_pair(pair,
                               f"magnetic_spacegroup_operation_index[{uni}][{slot}]")
    if any(pair != [0, 0] for pair in operation_index[0]):
        raise ExtractionError("magnetic operation index dummy mismatch")

    operations = msg["magnetic_symmetry_operations"]
    _artifact_list(operations, MSG_OPERATION_COUNT, "magnetic_symmetry_operations")
    if operations[0] != 0:
        raise ExtractionError("magnetic operation sentinel mismatch")
    if any(type(value) is not int or not 0 < value < MAGNETIC_OPERATION_ENCODING_LIMIT
           for value in operations[1:]):
        raise ExtractionError("magnetic operation encoding out of range")
    for value in operations:
        _decode_magnetic_operation(value)

    active_spans = []
    for uni in range(1, MSG_UNI_COUNT):
        hall_count, first_hall = mapping[uni]
        if hall_count < 1 or hall_count > MSG_HALL_SLOTS:
            raise ExtractionError("UNI Hall count out of range")
        if (first_hall < 1 or first_hall > SPG_HALL_SETTINGS
                or first_hall + hall_count - 1 > SPG_HALL_SETTINGS):
            raise ExtractionError("UNI mapping range mismatch")
        for slot, (order, offset) in enumerate(operation_index[uni]):
            if slot < hall_count:
                if order <= 0 or offset < 1 or offset + order > len(operations):
                    raise ExtractionError("magnetic operation span out of range")
                active_spans.append((offset, offset + order))
            elif [order, offset] != [0, 0]:
                raise ExtractionError("nonzero operation index beyond Hall count")
    if len(active_spans) != MSG_ACTIVE_SPAN_COUNT:
        raise ExtractionError("magnetic active span census mismatch")
    previous_end = 1
    for start, end in sorted(active_spans):
        if start != previous_end or end <= start:
            raise ExtractionError("magnetic operation spans are not contiguous")
        previous_end = end
    if previous_end != len(operations):
        raise ExtractionError("magnetic operation span boundary mismatch")

    transformations = msg["alternative_transformations"]
    _artifact_list(transformations, MSG_UNI_COUNT, "alternative_transformations")
    transformation_value_count = 0
    for uni, row in enumerate(transformations):
        _artifact_list(row, MSG_HALL_SLOTS,
                       f"alternative_transformations[{uni}]")
        for slot, values in enumerate(row):
            _artifact_list(values, 7,
                           f"alternative_transformations[{uni}][{slot}]")
            if uni == 0 or slot >= mapping[uni][0]:
                if any(value != 0 for value in values):
                    raise ExtractionError("nonzero transformation in inactive slot")
                continue
            first_zero = next((index for index, value in enumerate(values)
                               if value == 0), len(values))
            if first_zero == len(values):
                raise ExtractionError("alternative transformation terminator missing")
            if any(value != 0 and not 0 < value < MSG_OPERATION_SCALE
                   for value in values[:first_zero]):
                raise ExtractionError("transformation encoding out of range")
            if any(value != 0 for value in values[first_zero:]):
                raise ExtractionError("nonzero transformation tail")
            transformation_value_count += first_zero
    if transformation_value_count != ALTERNATIVE_TRANSFORMATION_VALUE_COUNT:
        raise ExtractionError("alternative transformation value census mismatch")

    def decoded(raw):
        return _decode_magnetic_operation(raw)

    for raw, expected in MAGNETIC_DECODER_WITNESSES.items():
        if decoded(raw) != expected:
            raise ExtractionError(f"magnetic decoder witness mismatch: {raw}")

    def span(uni, slot):
        order, offset = operation_index[uni][slot]
        return operations[offset:offset + order]

    if span(7, 0) != UNI7_RAW_WITNESS:
        raise ExtractionError("UNI7 raw witness mismatch")
    for slot, expected in enumerate(UNI9_RAW_WITNESSES):
        if span(9, slot) != expected:
            raise ExtractionError(f"UNI9 raw witness mismatch at slot {slot}")


def validate_artifact(data):
    """Validate a parsed artifact independently of its manifest."""
    _validate_artifact_tree(data)
    _artifact_keys(data, {"schema", "translation_denominator", "spg", "msg"},
                   "artifact")
    if data["schema"] != SCHEMA:
        raise ExtractionError("artifact schema mismatch")
    if data["translation_denominator"] != TRANSLATION_DENOMINATOR:
        raise ExtractionError("artifact translation denominator mismatch")
    _validate_artifact_spg(data["spg"])
    _validate_artifact_msg(data["msg"])


def _parse_json_bytes(data, name="JSON"):
    if not isinstance(data, bytes):
        raise ExtractionError(f"{name}: expected bytes")
    try:
        value = json.loads(
            data.decode("utf-8", errors="strict"),
            parse_float=_reject_json_float,
            parse_constant=_reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ExtractionError) as error:
        if isinstance(error, ExtractionError):
            raise
        raise ExtractionError(f"{name}: invalid strict JSON") from error
    _validate_json_value(value)
    return value


def _load_json(path):
    path = Path(path)
    try:
        data = path.read_bytes()
    except OSError as error:
        raise ExtractionError(f"unable to read JSON {path}") from error
    return _parse_json_bytes(data, str(path))


def validate_manifest(manifest, artifact_bytes, artifact_name):
    """Validate a manifest's closed provenance and artifact commitment."""
    if not isinstance(artifact_bytes, bytes):
        raise ExtractionError("artifact commitment must be bytes")
    _validate_json_value(manifest)
    if not isinstance(manifest, dict):
        raise ExtractionError("manifest must be a JSON object")
    expected_keys = {
        "schema", "repository", "tag", "commit", "sources",
        "extractor_schema_version", "artifact",
    }
    if set(manifest) != expected_keys:
        raise ExtractionError("manifest schema keys mismatch")
    if manifest["schema"] != MANIFEST_SCHEMA:
        raise ExtractionError("manifest schema mismatch")
    if manifest["repository"] != "https://github.com/spglib/spglib":
        raise ExtractionError("manifest repository mismatch")
    if manifest["tag"] != EXPECTED_UPSTREAM_TAG:
        raise ExtractionError("manifest tag mismatch")
    if manifest["commit"] != EXPECTED_UPSTREAM_COMMIT:
        raise ExtractionError("manifest commit mismatch")
    if manifest["extractor_schema_version"] != EXTRACTOR_VERSION:
        raise ExtractionError("manifest extractor version mismatch")
    expected_sources = [
        {"path": "src/msg_database.c", "sha256": EXPECTED_SOURCES["msg_database.c"]},
        {"path": "src/spg_database.c", "sha256": EXPECTED_SOURCES["spg_database.c"]},
    ]
    if manifest["sources"] != expected_sources:
        raise ExtractionError("manifest source commitments mismatch")
    expected_artifact = {
        "path": str(artifact_name),
        "bytes": len(artifact_bytes),
        "sha256": hashlib.sha256(artifact_bytes).hexdigest(),
    }
    if manifest["artifact"] != expected_artifact:
        raise ExtractionError("manifest artifact commitment mismatch")


def canonical_json(value):
    _validate_json_value(value)
    try:
        encoded = json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        return (encoded + "\n").encode("utf-8")
    except (TypeError, ValueError, UnicodeError) as error:
        raise ExtractionError("value cannot be represented as strict JSON") from error


def parse_and_validate_committed_pair(artifact_bytes: bytes,
                                      manifest_bytes: bytes,
                                      artifact_name: str) -> dict:
    """Parse and close a caller-provided artifact/manifest byte pair."""
    if type(artifact_name) is not str:
        raise ExtractionError("artifact name must be text")
    artifact = _parse_json_bytes(artifact_bytes, "artifact")
    manifest = _parse_json_bytes(manifest_bytes, "manifest")
    validate_artifact(artifact)
    validate_manifest(manifest, artifact_bytes, artifact_name)
    if canonical_json(artifact) != artifact_bytes:
        raise ExtractionError("artifact is not canonical JSON")
    if canonical_json(manifest) != manifest_bytes:
        raise ExtractionError("manifest is not canonical JSON")
    return artifact


def _validate_output_targets(output, manifest):
    output = Path(output)
    manifest = Path(manifest)
    try:
        if output.resolve(strict=False) == manifest.resolve(strict=False):
            raise ExtractionError("artifact and manifest paths must differ")
        output_stat = os.stat(output)
    except FileNotFoundError:
        output_stat = None
    except OSError as error:
        raise ExtractionError("unable to inspect artifact/manifest paths") from error
    try:
        manifest_stat = os.stat(manifest)
    except FileNotFoundError:
        manifest_stat = None
    except OSError as error:
        raise ExtractionError("unable to inspect artifact/manifest paths") from error
    if (output_stat is not None and manifest_stat is not None
            and output_stat.st_dev == manifest_stat.st_dev
            and output_stat.st_ino == manifest_stat.st_ino):
        raise ExtractionError("artifact and manifest must not share an inode")


def _atomic_write(path, data):
    path = Path(path)
    temporary_path = None
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
        )
        temporary_path = Path(temporary_name)
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
        temporary_path = None
    except OSError as error:
        raise ExtractionError(f"unable to atomically write {path}") from error
    finally:
        if temporary_path is not None:
            try:
                temporary_path.unlink()
            except FileNotFoundError:
                pass
            except OSError:
                pass


def write_outputs(upstream, output, manifest):
    output = Path(output)
    manifest = Path(manifest)
    _validate_output_targets(output, manifest)
    artifact, details = extract(upstream)
    artifact_bytes = canonical_json(artifact)
    manifest_value = {
        "schema": MANIFEST_SCHEMA,
        "repository": "https://github.com/spglib/spglib",
        "tag": EXPECTED_UPSTREAM_TAG,
        "commit": EXPECTED_UPSTREAM_COMMIT,
        "sources": [
            {"path": "src/msg_database.c", "sha256": details["msg_database.c"]},
            {"path": "src/spg_database.c", "sha256": details["spg_database.c"]},
        ],
        "extractor_schema_version": EXTRACTOR_VERSION,
        "artifact": {"path": str(output.name), "bytes": len(artifact_bytes),
                      "sha256": hashlib.sha256(artifact_bytes).hexdigest()},
    }
    _atomic_write(output, artifact_bytes)
    _atomic_write(manifest, canonical_json(manifest_value))
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
