#!/usr/bin/env python3
"""Assemble the full frozen embedding table from the derived evidence.

Inputs
------
* the pinned OR-1 census (``U`` for every record the official program prints),
* the per-record child-shift solution (engine-validated),
* the compact-spelling evidence for the 112 parent irreps whose legacy compound
  spelling makes the official program print an empty table.

Output is the Rust module the engine reads.  Every entry is a convention derived
from the official program; the child shift is accepted only where the engine's
own validation accepted it.

This assembler is the **only** writer of
``src/irrep/subduction_settings_data.rs`` (the older
``scripts/generate_subduction_settings.py`` derives conventions for its 75-record
task-8 subset and refuses to write the module).  It validates before it writes:

* a repeated ordinal with a different payload is a conflict, not an overwrite,
* a repeated ordinal with the same payload is still a duplicated input line,
* the assembled ordinal set must equal the set the output module already carries
  (or an explicit ``--expected`` list), so a truncated census cannot silently
  produce a shorter table,
* unmatched census records are reported by status,
* the file is written through a temporary file and ``os.replace``.

``--check`` assembles in memory and compares with the committed module instead of
writing; ``--partial`` allows an incomplete table for experiments and says so.
"""
import argparse
import json
import os
import sys
import tempfile
from fractions import Fraction
from math import gcd

HEADER = '''//! Generated per-ordinal embedding settings for the full-subduction engine.
//!
//! `DO NOT EDIT`: this module is assembled by the task-9 settings pipeline
//! (`scripts/task9/README.md`); `scripts/task9/build_table.py --check`
//! re-assembles it and compares byte for byte, and the older
//! `scripts/generate_subduction_settings.py` refuses to write it.
//!
//! Every entry is derived from the pinned ISOTROPY archive
//! `isotropy_subgroup/iso.zip`
//! (SHA-256 `568667bfc8027095537d642297b319c872d00016b868143c666f90d5931d9f7b`),
//! official `iso` 9.6.1:
//!
//! * `U` is the exact change of basis
//!   `W . P_parent . (P_sub . B_oracle)^-1` the printed basis implies, stored as
//!   `U numerator / U denominator` because a few monoclinic records reach the
//!   official cell only through a fractional change of basis;
//! * `child_shift` is the origin correction between the recorded ITA setting of
//!   the subgroup and the shipped `SG_DATA_HALL` row of the subgroup, in the
//!   subgroup's own conventional cell.  It is present only where the engine's
//!   validation accepted it.
//!
//! A record whose legacy compound Miller-Love spelling (`W1W1`) makes the
//! official program answer with an empty table is addressed through its compact
//! `little_irr_full_label` spelling (`W1WA1`); the mapping is proved per record
//! by requiring the printed table to reproduce the stored rows exactly.

use crate::mathfunc::Mat3I;

/// Identity setting transform.
const IDENTITY: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// Cyclic axis permutation used by the rhombohedral and C-centred records.
const CYCLIC: Mat3I = [[0, 0, 1], [1, 0, 0], [0, 1, 0]];

/// No child origin shift: the isotropy record already uses the Hall frame.
const NO_SHIFT: [i32; 4] = [0, 0, 0, 1];

/// `(ordinal, parent SG, subgroup SG, U numerator, U denominator, child shift)`.
pub type FrozenEmbeddingSetting = (usize, u8, u8, Mat3I, i32, [i32; 4]);

/// One entry per covered isotropy record ordinal.
pub static FROZEN_EMBEDDING_SETTINGS: &[FrozenEmbeddingSetting] = &[
'''


def rationalise(entries):
    values = [[Fraction(v) for v in row] for row in entries]
    denominator = 1
    for row in values:
        for value in row:
            denominator = denominator * value.denominator // gcd(denominator, value.denominator)
    return [[int(value * denominator) for value in row] for row in values], denominator


def reduce_shift(shift):
    """Reduce `(x, y, z, d)`; the same convention must have one spelling."""
    x, y, z, d = (int(v) for v in shift)
    if d <= 0:
        raise SystemExit(f"non-positive shift denominator {d}")
    divisor = gcd(gcd(abs(x), abs(y)), gcd(abs(z), d))
    if divisor == 0:
        return [0, 0, 0, 1]
    return [x // divisor, y // divisor, z // divisor, d // divisor]


def rust_matrix(numerator, denominator):
    if numerator == [[1, 0, 0], [0, 1, 0], [0, 0, 1]] and denominator == 1:
        return "IDENTITY"
    if numerator == [[0, 0, 1], [1, 0, 0], [0, 1, 0]] and denominator == 1:
        return "CYCLIC"
    return "[" + ", ".join(
        "[" + ", ".join(str(v) for v in row) + "]" for row in numerator
    ) + "]"


def load_shifts(path):
    """Child shifts from any of the pipeline's evidence shapes.

    The derivation tools write three shapes: ``{"solved": {ordinal: shift}}``,
    ``{"accepted": [...], "shifts": {ordinal: shift}}`` (apply only the accepted
    ones) and a bare ``{ordinal: shift}`` map.  A shift may also be wrapped as
    ``{"shift": [...]}``.  Accepting all of them lets one command express the
    multi-pass chain the table was actually frozen by, instead of re-running the
    assembler once per fix file.
    """
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    out = {}
    if isinstance(data, dict) and "shifts" in data:
        accepted = data.get("accepted")
        for key, value in data["shifts"].items():
            ordinal = int(key)
            if accepted is not None and ordinal not in accepted:
                continue
            if isinstance(value, dict):
                value = value.get("shift")
            if isinstance(value, list):
                out[ordinal] = [int(v) for v in value]
        return out
    source = data.get("solved", data) if isinstance(data, dict) else {}
    for key, value in source.items():
        if isinstance(value, dict):
            value = value.get("shift")
        if isinstance(value, list):
            out[int(key)] = [int(v) for v in value]
    return out


def load_jsonl(path, only_records):
    """JSON lines in file order.

    The census carries `type` tags and only its `"record"` lines describe
    ordinals; the compact-label evidence (`empty_fixed.jsonl`) has no `type` at
    all, so filtering it the same way would silently drop all 195 lines and the
    assembler would report them as uncovered.
    """
    records = []
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                entry = json.loads(line)
            except json.JSONDecodeError as error:
                raise SystemExit(f"{path}:{number}: {error}") from error
            if only_records and entry.get("type") != "record":
                continue
            records.append(entry)
    return records


def split_top_level(text):
    """Split on commas that are not inside `[...]`."""
    fields = []
    depth = 0
    current = ""
    for character in text:
        if character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
        if character == "," and depth == 0:
            fields.append(current.strip())
            current = ""
        else:
            current += character
    fields.append(current.strip())
    return fields


def parse_committed(path):
    """The `(ordinal, parent, child)` identity of every row of a module.

    Used both for `--check` and to learn the ordinal set the shipped table must
    keep covering, so a truncated input cannot pass as a complete table.  The
    matrix field is not decoded here: `--check` compares the rendered text, and
    the expected set only needs the ordinal.
    """
    entries = {}
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            stripped = line.split("//")[0].strip()
            if not stripped.startswith("("):
                continue
            fields = split_top_level(stripped.rstrip(",").strip("()"))
            if len(fields) != 6:
                raise SystemExit(f"{path}: cannot parse entry {line!r}")
            ordinal = int(fields[0])
            if ordinal in entries:
                raise SystemExit(f"{path}: ordinal {ordinal} appears twice")
            entries[ordinal] = (int(fields[1]), int(fields[2]))
    return entries


def assemble(census, shifts, empty):
    """Build the ordinal -> row map, refusing duplicates with any payload."""
    entries = {}
    duplicates = []
    for source, records in (("census", census), ("compact-label", empty)):
        for record in records:
            if source == "census" and record.get("u") is None:
                continue
            ordinal = record["ordinal"]
            numerator, denominator = rationalise(record["u"])
            shift = reduce_shift(shifts.get(ordinal, [0, 0, 0, 1]))
            row = (
                record["parent"], record["child"], numerator, denominator, shift,
            )
            if ordinal in entries:
                duplicates.append((ordinal, source, entries[ordinal], row))
                continue
            entries[ordinal] = row
    return entries, duplicates


def render(entries):
    lines = [HEADER]
    for ordinal in sorted(entries):
        parent, child, numerator, denominator, shift = entries[ordinal]
        if shift == [0, 0, 0, 1]:
            shift_text = "NO_SHIFT"
        else:
            shift_text = "[" + ", ".join(str(v) for v in shift) + "]"
        lines.append(
            f"    ({ordinal}, {parent}, {child},"
            f" {rust_matrix(numerator, denominator)}, {denominator}, {shift_text}),\n"
        )
    lines.append("];\n")
    return "".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", required=True)
    parser.add_argument("--shifts", required=True)
    parser.add_argument(
        "--derived",
        action="append",
        default=[],
        help="oracle-derived shifts (preferred); repeatable, later files win",
    )
    parser.add_argument("--legacy", help="previously frozen, independently verified shifts")
    parser.add_argument("--empty", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--expected",
        help="JSON list of every ordinal the table must cover; defaults to the "
             "ordinal set of the existing --out module",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="assemble in memory and compare with --out instead of writing",
    )
    parser.add_argument(
        "--partial",
        action="store_true",
        help="allow a table that does not cover the expected ordinals (experiments only)",
    )
    args = parser.parse_args(argv)

    census = load_jsonl(args.census, only_records=True)
    shifts = load_shifts(args.shifts)
    # Precedence: the engine-validated search base, then the oracle-derived
    # shifts (later files win), then the independently verified legacy pin.
    for path in args.derived:
        layer = load_shifts(path)
        shifts.update(layer)
        print(f"oracle-derived shifts from {os.path.basename(path)}: {len(layer)}")
    if args.legacy:
        legacy = load_shifts(args.legacy)
        shifts.update(legacy)
        print(f"legacy verified shifts pinned: {len(legacy)}")
    empty = load_jsonl(args.empty, only_records=False)

    committed = parse_committed(args.out)
    if args.expected:
        with open(args.expected, encoding="utf-8") as handle:
            expected = {int(value) for value in json.load(handle)}
    elif committed is not None:
        expected = set(committed)
    else:
        expected = None

    entries, duplicates = assemble(census, shifts, empty)
    if duplicates:
        for ordinal, source, first, second in duplicates[:5]:
            print(f"ordinal {ordinal}: duplicated by {source}: {first} vs {second}",
                  file=sys.stderr)
        raise SystemExit(f"{len(duplicates)} duplicated ordinals in the inputs")

    missing = [record["ordinal"] for record in census if record["ordinal"] not in entries]
    kinds = {}
    for record in census:
        if record["ordinal"] in missing:
            kinds[record["status"]] = kinds.get(record["status"], 0) + 1

    text = render(entries)
    print(f"entries: {len(entries)}")
    print(f"records still without a frozen convention: {len(missing)}")
    if missing:
        print("  by census status:", kinds)

    if expected is not None:
        uncovered = sorted(expected - set(entries))
        extra = sorted(set(entries) - expected)
        if uncovered or extra:
            message = (
                f"the assembled table does not cover the expected ordinals: "
                f"{len(uncovered)} missing (first {uncovered[:5]}), "
                f"{len(extra)} unexpected (first {extra[:5]})"
            )
            if not args.partial:
                raise SystemExit(message)
            print(f"warning: --partial: {message}")

    if args.check:
        if not os.path.exists(args.out):
            raise SystemExit(f"--check: {args.out} does not exist")
        with open(args.out, encoding="utf-8") as handle:
            current = handle.read()
        if current != text:
            raise SystemExit(
                f"{args.out} differs from the freshly assembled table "
                "(rebuild without --check to replace it)"
            )
        print(f"{args.out}: identical to the freshly assembled table")
        return

    directory = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(directory, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(dir=directory, suffix=".tmp")
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(text)
        os.replace(temporary, args.out)
    except BaseException:
        os.unlink(temporary)
        raise
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
