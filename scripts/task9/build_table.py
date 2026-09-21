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
"""
import argparse
import json
from fractions import Fraction
from math import gcd

HEADER = '''//! Generated per-ordinal embedding settings for the full-subduction engine.
//!
//! `DO NOT EDIT`: regenerate with the task-9 settings pipeline
//! (`scripts/generate_subduction_settings.py`).
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", required=True)
    parser.add_argument("--shifts", required=True)
    parser.add_argument("--derived", help="oracle-derived shifts (preferred)")
    parser.add_argument("--legacy", help="previously frozen, independently verified shifts")
    parser.add_argument("--empty", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    census = [
        json.loads(line) for line in open(args.census, encoding="utf-8")
        if json.loads(line).get("type") == "record"
    ]
    shifts = json.load(open(args.shifts, encoding="utf-8"))["solved"]
    # Precedence: the shifts frozen and verified before this sweep, then the
    # shifts the official program derives, then the engine-validated search.
    if args.derived:
        derived = json.load(open(args.derived, encoding="utf-8"))
        for ordinal in derived["accepted"]:
            shifts[str(ordinal)] = derived["shifts"][str(ordinal)]
        print(f"oracle-derived shifts accepted by the engine: {len(derived['accepted'])}")
    if args.legacy:
        shifts.update(json.load(open(args.legacy, encoding="utf-8")))
        print(f"legacy verified shifts pinned: {len(json.load(open(args.legacy, encoding='utf-8')))}")
    empty = [
        json.loads(line) for line in open(args.empty, encoding="utf-8")
    ]

    entries = {}
    for record in census:
        if record.get("u") is None:
            continue
        numerator, denominator = rationalise(record["u"])
        shift = reduce_shift(shifts.get(str(record["ordinal"]), [0, 0, 0, 1]))
        entries[record["ordinal"]] = (
            record["parent"], record["child"], numerator, denominator, shift,
        )
    added = 0
    for record in empty:
        numerator, denominator = rationalise(record["u"])
        shift = reduce_shift(shifts.get(str(record["ordinal"]), [0, 0, 0, 1]))
        entries[record["ordinal"]] = (
            record["parent"], record["child"], numerator, denominator, shift,
        )
        added += 1
    print(f"entries: {len(entries)} (census {len(entries) - added}, compact-label {added})")

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
    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("".join(lines))
    missing = [r["ordinal"] for r in census
               if r.get("u") is None and r["ordinal"] not in entries]
    print(f"records still without a frozen convention: {len(missing)}")
    if missing:
        kinds = {}
        for record in census:
            if record["ordinal"] in missing:
                kinds[record["status"]] = kinds.get(record["status"], 0) + 1
        print("  by census status:", kinds)


if __name__ == "__main__":
    main()
