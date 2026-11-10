#!/usr/bin/env python3
"""Generate the pinned collinear spin-space-group table.

Reads the verified FindSpinGroup download manifest (target/, not committed),
extracts the collinear quotient signature ``(spatial operation, channel)``
from each record and emits ``src/spin_group/database/ssg_collinear_data.rs``.
The generated file records the database release id/name, the source SHA-256,
the manifest SHA-256 and the identifier-list SHA-256.  It is deterministic:
records are sorted by their four-index number and operations by
(rotation, translation, channel).  The committed table is post-processed with
`cargo fmt`; run `cargo fmt -p cryspglib` after regeneration so a fresh checkout
and a regenerated table stay byte-identical.

This script is a data-import tool only; all acceptance checks of the generated
table live in Rust integration tests.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RECORDS = ROOT / "target" / "ssg-collinear-records-v1"
DEFAULT_OUTPUT = ROOT / "src" / "spin_group" / "database" / "ssg_collinear_data.rs"

NUMBER_RE = re.compile(r"^(\d+)\.(\d+)\.(\d+)\.(\d+)\.([A-Za-z])$")
RATIONAL_RE = re.compile(r"-?\d+(?:/\d+)?")


def rust_string(value: str) -> str:
    out = ['"']
    for ch in value:
        if ch == "\\":
            out.append("\\\\")
        elif ch == '"':
            out.append('\\"')
        elif ch == "\n":
            out.append("\\n")
        elif ch == "\t":
            out.append("\\t")
        elif ch == "\r":
            out.append("\\r")
        elif ord(ch) < 0x20:
            out.append(f"\\u{{{ord(ch):x}}}")
        else:
            out.append(ch)
    out.append('"')
    return "".join(out)


def plain_symbol(notation: str) -> str:
    text = notation.replace("$", "")
    text = text.replace("\\ce{", "").replace("\\ce", "")
    text = text.replace("\\infty", "∞")
    text = text.replace("\\text{-}", "-")
    text = text.replace("\\text{", "")
    text = text.replace("\\overline{", "")
    text = text.replace("\\", "")
    text = text.replace("{", "").replace("}", "").replace("^", "")
    return "".join(text.split())


def rationals(text: str) -> list[Fraction]:
    return [Fraction(token) for token in RATIONAL_RE.findall(text)]


def symbolic_spin_matrix(expression: str) -> list[list[str]]:
    """Normalize a symbolic spin expression into a 3x3 cell grid."""
    cleaned = expression.replace(" ", "")
    cleaned = cleaned.replace("cos(\\phi)", "C").replace("cos(phi)", "C")
    cleaned = cleaned.replace("sin(\\phi)", "S").replace("sin(phi)", "S")
    cleaned = cleaned.strip()
    if not (cleaned.startswith("[") and cleaned.endswith("]")):
        raise ValueError(f"unrecognised symbolic spin expression {expression!r}")
    cells = cleaned[1:-1].split(",")
    if len(cells) != 9:
        raise ValueError(f"symbolic spin expression is not 3x3: {expression!r}")
    return [cells[0:3], cells[3:6], cells[6:9]]


def infinity_axis(expression: str) -> int:
    cells = symbolic_spin_matrix(expression)
    for axis in range(3):
        unit_column = [cells[row][axis] for row in range(3)]
        unit_row = [cells[axis][column] for column in range(3)]
        if unit_column == ["0", "0", "0"]:
            continue
        if unit_column[axis] != "1" or unit_row[axis] != "1":
            continue
        if all(cells[row][axis] == "0" for row in range(3) if row != axis) and all(
            cells[axis][column] == "0" for column in range(3) if column != axis
        ):
            return axis
    raise ValueError(f"cannot determine rotation axis of {expression!r}")


PERMUTATIONS = {
    0: ((0, 0, 1), (0, 1, 0), (1, 0, 0)),
    1: ((1, 0, 0), (0, 0, 1), (0, 1, 0)),
    2: ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
}


def mat_mul(left, right):
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )


def mat_transpose(matrix):
    return tuple(tuple(matrix[j][i] for j in range(3)) for i in range(3))


def spin_eta(expression: str, permutation) -> int:
    """Return the net action on the canonical z axis after rotating the spin frame."""
    if "cos" in expression or "sin" in expression:
        # A rotation about the (already detected) collinear axis maps z to z.
        return 1
    values = rationals(expression)
    if len(values) != 9:
        raise ValueError(f"spin expression is not 3x3: {expression!r}")
    matrix = tuple(tuple(values[3 * i + j] for j in range(3)) for i in range(3))
    # Orthogonality and proper/improper character are checked as data invariants.
    product = mat_mul(mat_transpose(matrix), matrix)
    if product != ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        raise ValueError(f"spin matrix is not orthogonal: {expression!r}")
    conjugated = mat_mul(mat_mul(permutation, matrix), mat_transpose(permutation))
    column = (conjugated[0][2], conjugated[1][2], conjugated[2][2])
    if column == (0, 0, 1):
        return 1
    if column == (0, 0, -1):
        return -1
    raise ValueError(f"spin matrix does not preserve the collinear axis: {expression!r}")


def normalise_translation(value: Fraction) -> Fraction:
    result = value % 1
    return Fraction(0) if result == 1 else result


def extract_signature(document: dict, identifier: str):
    operations = document.get("operations")
    if not operations:
        raise ValueError(f"{identifier}: no operations")
    infinity_axes = set()
    for operation in operations:
        if "cos" in operation["spin_expr"] or "sin" in operation["spin_expr"]:
            infinity_axes.add(infinity_axis(operation["spin_expr"]))
    if len(infinity_axes) > 1:
        raise ValueError(f"{identifier}: multiple continuous spin axes {sorted(infinity_axes)}")
    if infinity_axes:
        axis = infinity_axes.pop()
    else:
        # Fall back to the axis preserved by every discrete spin matrix.
        candidates = []
        for candidate in range(3):
            try:
                if all(
                    spin_eta(operation["spin_expr"], PERMUTATIONS[candidate]) is not None
                    for operation in operations
                ):
                    candidates.append(candidate)
            except ValueError:
                continue
        if len(candidates) != 1:
            raise ValueError(f"{identifier}: cannot determine collinear spin axis")
        axis = candidates[0]
    permutation = PERMUTATIONS[axis]
    signature = {}
    for operation in operations:
        space = rationals(operation["space_expr"])
        if len(space) != 12:
            raise ValueError(f"{identifier}: space expression is not 3x3+3: {operation['space_expr']!r}")
        rotation = tuple(
            tuple(int(space[3 * row + column]) for column in range(3)) for row in range(3)
        )
        translation = tuple(normalise_translation(space[9 + axis]) for axis in range(3))
        eta = spin_eta(operation["spin_expr"], permutation)
        key = (rotation, translation)
        if key in signature and signature[key] != eta:
            raise ValueError(f"{identifier}: conflicting channel for spatial operation {key}")
        signature[key] = eta
    # Identity must be present and channel-preserving.
    identity_key = (((1, 0, 0), (0, 1, 0), (0, 0, 1)), (0, 0, 0))
    if signature.get(identity_key) != 1:
        raise ValueError(f"{identifier}: no channel-preserving identity")
    # The full closure/homomorphism check is performed by the Rust tests on
    # the generated table (the raw spatial closure check here is quadratic in
    # the largest 384-element records and is not needed for data extraction).
    ordered = sorted(
        ((rotation, translation, eta) for (rotation, translation), eta in signature.items()),
        key=lambda entry: (entry[0], entry[1], entry[2]),
    )
    return ordered


def number_parts(identifier: str):
    match = NUMBER_RE.match(identifier)
    if match is None:
        raise ValueError(f"bad identifier {identifier!r}")
    g0, l0, k_index, map_index, letter = match.groups()
    return int(g0), int(l0), int(k_index), int(map_index), letter


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args()

    manifest_path = args.records / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("problems") and not args.allow_partial:
        raise SystemExit(f"manifest has problems: {manifest['problems'][:3]}")
    expected = manifest["count"]
    if len(manifest["records"]) != expected and not args.allow_partial:
        raise SystemExit(f"manifest records {len(manifest['records'])} != {expected}")

    rows = []
    signature_index = {}
    total_ops = 0
    for position, record in enumerate(manifest["records"], 1):
        if position % 50 == 0:
            print(f'processed {position}/{len(manifest["records"])}', flush=True)
        identifier = record["identifier"]
        document = json.loads((args.records / (identifier + ".json")).read_text())
        if document["path"]["index_label"] != identifier:
            raise SystemExit(f"{identifier}: path label mismatch")
        if document["path"]["configuration"] != "Collinear":
            raise SystemExit(f"{identifier}: not a collinear record")
        g0, l0, k_index, map_index, letter = number_parts(identifier)
        signature = extract_signature(document, identifier)
        total_ops += len(signature)
        meta = document["group"]
        rows.append(
            {
                "g0": g0,
                "l0": l0,
                "k_index": k_index,
                "map_index": map_index,
                "letter": letter,
                "identifier": identifier,
                "international_symbol": plain_symbol(meta["international_notation"]),
                "international_notation_latex": meta["international_notation"],
                "group_type": document["path"]["group_type"],
                "basis": meta.get("basis", ""),
                "spin_point_group": meta["spin_part_point_group"],
                "tspin": int(meta["tspin"]),
                "gs": int(meta["gs"]),
                "gso": meta["gso"],
                "num_all_elements": int(meta["num_all_elements"]),
                "operations": signature,
            }
        )
        hashable = tuple(
            (rotation, tuple(translation), eta) for rotation, translation, eta in signature
        )
        if hashable in signature_index:
            raise SystemExit(
                f"{identifier}: duplicate quotient signature with {signature_index[hashable]}"
            )
        signature_index[hashable] = identifier

    rows.sort(
        key=lambda row: (
            row["g0"],
            row["l0"],
            row["k_index"],
            row["map_index"],
            row["letter"],
        )
    )
    if len(rows) != expected and not args.allow_partial:
        raise SystemExit(f"selected records {len(rows)} != {expected}")

    release = manifest.get("release") or {}
    lines = [
        "//! Generated by `scripts/generate_ssg_collinear_data.py`; do not edit.",
        "//!",
        "//! Source: FindSpinGroup database (<https://database.findspingroup.com>),",
        "//! filtered by `configuration=Collinear`; the quotient signature stores",
        "//! only spatial operations and their net channel action.",
        "",
        "use super::{SpinSpaceGroupNumber, SpinSpaceGroupOperation, SpinSpaceGroupRecord};",
        "use crate::spin_group::SpinGroupKind;",
        "use crate::spin_group::collinear::CollinearAction;",
        "",
        f"pub(crate) const RECORD_COUNT: usize = {len(rows)};",
        f"pub(crate) const RELEASE_ID: u32 = {int(release.get('id') or 0)};",
        f"pub(crate) const RELEASE_NAME: &str = {rust_string(str(release.get('name') or ''))};",
        f"pub(crate) const SOURCE_SHA256: &str = {rust_string(str(release.get('source_sha256') or ''))};",
        f"pub(crate) const MANIFEST_SHA256: &str = {rust_string(manifest.get('manifest_sha256') or hashlib.sha256(manifest_path.read_bytes()).hexdigest())};",
        f"pub(crate) const IDENTIFIERS_SHA256: &str = {rust_string(str(manifest.get('identifiers_sha256') or ''))};",
        "",
        "#[rustfmt::skip]",
        "pub(crate) static RECORDS: &[SpinSpaceGroupRecord] = &[",
    ]
    for row in rows:
        lines.append("    SpinSpaceGroupRecord {")
        lines.append(
            "        number: SpinSpaceGroupNumber::from_parts_unchecked("
            f"{row['g0']}, {row['l0']}, {row['k_index']}, {row['map_index']}, b'{row['letter']}'),"
        )
        lines.append(f"        identifier: {rust_string(row['identifier'])},")
        lines.append(f"        international_symbol: {rust_string(row['international_symbol'])},")
        lines.append(
            f"        international_notation_latex: {rust_string(row['international_notation_latex'])},"
        )
        lines.append("        configuration: SpinGroupKind::Collinear,")
        lines.append(f"        group_type: {rust_string(row['group_type'])},")
        lines.append(f"        basis: {rust_string(row['basis'])},")
        lines.append(f"        spin_point_group: {rust_string(row['spin_point_group'])},")
        lines.append(f"        tspin: {row['tspin']},")
        lines.append(f"        gs: {row['gs']},")
        lines.append(f"        gso: {rust_string(row['gso'])},")
        lines.append(f"        num_all_elements: {row['num_all_elements']},")
        lines.append("        operations: &[")
        for rotation, translation, eta in row["operations"]:
            rotation_literal = (
                "[["
                + "], [".join(
                    ", ".join(str(rotation[row][column]) for column in range(3))
                    for row in range(3)
                )
                + "]]"
            )
            translation_literal = (
                "["
                + ", ".join(f"{float(value)!r}" for value in translation)
                + "]"
            )
            channel = "CollinearAction::Preserve" if eta == 1 else "CollinearAction::Swap"
            lines.append(
                "            SpinSpaceGroupOperation { "
                f"rotation: {rotation_literal}, translation: {translation_literal}, channel: {channel} }},"
            )
        lines.append("        ],")
        lines.append("    },")
    lines.append("];")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n")
    print(
        f"wrote {args.output} records={len(rows)} operations={total_ops} "
        f"release={release.get('id')} source={release.get('source_sha256')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
