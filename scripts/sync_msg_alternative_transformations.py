#!/usr/bin/env python3
"""Sync MSG alternative setting transformations from upstream spglib.

The C table uses partial initializers such as ``{66459, 0}`` for its
``[][18][7]`` array.  Missing values are zero-initialized by C, so the
converter must pad both dimensions instead of discarding short rows.

Usage:
    python3 scripts/sync_msg_alternative_transformations.py \
        /path/to/spglib/src/msg_database.c
"""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
import re


DECLARATION = (
    "pub static ALTERNATIVE_TRANSFORMATIONS: [[[i32; 7]; 18]; 1652] = ["
)


def extract_c_initializer(source: str, array_name: str) -> str:
    """Return the balanced outer initializer for a named C array."""
    name_pos = source.find(array_name)
    if name_pos < 0:
        raise ValueError(f"array {array_name!r} was not found")

    equals_pos = source.find("=", name_pos)
    start = source.find("{", equals_pos)
    if equals_pos < 0 or start < 0:
        raise ValueError(f"initializer for {array_name!r} was not found")

    depth = 0
    for end in range(start, len(source)):
        if source[end] == "{":
            depth += 1
        elif source[end] == "}":
            depth -= 1
            if depth == 0:
                return source[start : end + 1]

    raise ValueError(f"initializer for {array_name!r} is unbalanced")


def parse_c_table(source: str) -> list[list[list[int]]]:
    initializer = extract_c_initializer(source, "alternative_transformations")
    initializer = re.sub(r"/\*.*?\*/", "", initializer, flags=re.DOTALL)
    python_literal = initializer.replace("{", "[").replace("}", "]")
    table = ast.literal_eval(python_literal)

    if len(table) != 1652:
        raise ValueError(f"expected 1652 UNI entries, found {len(table)}")

    normalized = []
    for uni_number, settings in enumerate(table):
        if len(settings) > 18:
            raise ValueError(
                f"UNI {uni_number} has {len(settings)} settings; maximum is 18"
            )

        normalized_settings = []
        for setting in settings:
            if len(setting) > 7:
                raise ValueError(
                    f"UNI {uni_number} has a transformation row longer than 7"
                )
            normalized_settings.append(setting + [0] * (7 - len(setting)))

        while len(normalized_settings) < 18:
            normalized_settings.append([0] * 7)
        normalized.append(normalized_settings)

    return normalized


def render_rust_table(table: list[list[list[int]]]) -> str:
    lines = [DECLARATION]
    for settings in table:
        rows = ["[" + ", ".join(map(str, row)) + "]" for row in settings]
        lines.append("[" + ", ".join(rows) + "],")
    lines.append("];")
    return "\n".join(lines) + "\n"


def replace_rust_table(source: str, rendered_table: str) -> str:
    pattern = re.compile(
        r"pub static ALTERNATIVE_TRANSFORMATIONS:"
        r" \[\[\[i32; 7\]; 18\]; 1652\] = \[.*?\n\];\s*\Z",
        flags=re.DOTALL,
    )
    updated, count = pattern.subn(rendered_table, source)
    if count != 1:
        raise ValueError(
            "expected exactly one trailing ALTERNATIVE_TRANSFORMATIONS table"
        )
    return updated


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("upstream_msg_database_c", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "src/msg_database_gen.rs",
    )
    args = parser.parse_args()

    upstream_source = args.upstream_msg_database_c.read_text(encoding="utf-8")
    rust_source = args.output.read_text(encoding="utf-8")
    table = parse_c_table(upstream_source)
    updated = replace_rust_table(rust_source, render_rust_table(table))
    args.output.write_text(updated, encoding="utf-8")

    nontrivial = sum(any(row) for settings in table for row in settings)
    print(
        f"updated {args.output} with {nontrivial} nontrivial setting rows"
    )


if __name__ == "__main__":
    main()
