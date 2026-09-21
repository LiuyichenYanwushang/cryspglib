#!/usr/bin/env python3
"""Emit the frozen embedding table from the OR-1 census (candidate records only).

This is the fast experimental path: every census ``candidate`` already carries the
official U, so the table can be emitted without re-deriving geometry.  Records the
census could not classify keep the engine's own candidate search.
"""
import argparse
import json

HEADER = '''//! Generated per-ordinal embedding settings for the full-subduction engine.
//!
//! `DO NOT EDIT`: regenerate with the task-9 settings generator.
//!
//! EXPERIMENTAL: emitted from the pinned OR-1 setting census.  Every entry is the
//! exact change of basis `U = W . P_parent . (P_sub . B_oracle)^-1` derived from
//! the official `iso` 9.6.1 table recorded in `isotropy_subgroup/iso.zip`
//! (SHA-256 `568667bfc8027095537d642297b319c872d00016b868143c666f90d5931d9f7b`).
//! `child_shift` is not yet derived for the newly added records.

use crate::mathfunc::Mat3I;

/// Identity setting transform.
const IDENTITY: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// No child origin shift: the isotropy record already uses the Hall frame.
const NO_SHIFT: [i32; 4] = [0, 0, 0, 1];

/// `(ordinal, parent SG, subgroup SG, U numerator, U denominator, child shift)`.
///
/// The convention is the exact rational matrix `U numerator / U denominator`;
/// the denominator is one for every integral record.
pub type FrozenEmbeddingSetting = (usize, u8, u8, Mat3I, i32, [i32; 4]);

/// One entry per covered isotropy record ordinal.
pub static FROZEN_EMBEDDING_SETTINGS: &[FrozenEmbeddingSetting] = &[
'''


def rationalise(entries):
    """`(numerator matrix, denominator)` for a matrix of fraction strings."""
    from fractions import Fraction

    values = [[Fraction(v) for v in row] for row in entries]
    denominator = 1
    for row in values:
        for value in row:
            denominator = denominator * value.denominator // gcd(denominator, value.denominator)
    numerator = [[int(value * denominator) for value in row] for row in values]
    return numerator, denominator


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def rust_matrix(u):
    if u == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
        return "IDENTITY"
    return "[" + ", ".join(
        "[" + ", ".join(str(v) for v in row) + "]" for row in u
    ) + "]"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("census")
    parser.add_argument("--out", required=True)
    parser.add_argument("--all", action="store_true",
                        help="also emit non-candidate records with an identity U")
    args = parser.parse_args()

    records = [
        json.loads(line)
        for line in open(args.census, encoding="utf-8")
        if json.loads(line).get("type") == "record"
    ]
    covered = [r for r in records if r.get("u") is not None]
    skipped = [r for r in records if r["status"] != "candidate"]
    covered.sort(key=lambda r: r["ordinal"])
    lines = [HEADER]
    seen = set()
    for record in covered:
        if record["ordinal"] in seen:
            raise SystemExit(f"duplicate ordinal {record['ordinal']}")
        seen.add(record["ordinal"])
        num, den = rationalise(record["u"])
        lines.append(
            f"    // {record['parent']} {record['irrep']} {record['direction_label']}"
            f" -> #{record['child']} (size {record['size']})"
            + (f" U = N/{den}" if den != 1 else "")
            + "\n"
            f"    ({record['ordinal']}, {record['parent']}, {record['child']},"
            f" {rust_matrix(num)}, {den}, NO_SHIFT),\n"
        )
    lines.append("];\n")
    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("".join(lines))
    print(f"emitted {len(covered)} entries, {len(skipped)} records left to the search")
    kinds = {}
    for record in skipped:
        kinds[record["detail"]] = kinds.get(record["detail"], 0) + 1
    print("skipped:", kinds)


if __name__ == "__main__":
    main()
