#!/usr/bin/env python3
"""Derive the child shift from the recorded/engine frame difference.

The isotropy table records each subgroup in the ISOTROPY frame, while the engine
maps into the parent frame of ``SG_DATA_HALL``.  Where the two frames differ, the
stored origin is expressed in the other frame and the frozen ``child_shift`` has
to undo exactly that difference:

    delta = T^-1 (o_printed - o_stored),   T = B^T

with `B` the printed basis of the record.  The four sign/transpose readings are
all probed, so the convention is settled by the engine's own validation rather
than by this file's algebra.
"""
import argparse
import json
import sys
from fractions import Fraction

import probe_census as probe_mod


def transpose(matrix):
    return [[matrix[j][i] for j in range(3)] for i in range(3)]


def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def inverse(matrix):
    det = (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )
    if det == 0:
        return None
    cofactors = [
        [
            (
                matrix[(i + 1) % 3][(j + 1) % 3] * matrix[(i + 2) % 3][(j + 2) % 3]
                - matrix[(i + 1) % 3][(j + 2) % 3] * matrix[(i + 2) % 3][(j + 1) % 3]
            )
            for j in range(3)
        ]
        for i in range(3)
    ]
    adjugate = transpose(cofactors)
    return [[adjugate[i][j] / det for j in range(3)] for i in range(3)]


def matvec(matrix, vector):
    return [sum(matrix[i][j] * vector[j] for j in range(3)) for i in range(3)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", required=True)
    parser.add_argument("--rejected", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    census = {
        json.loads(line)["ordinal"]: json.loads(line)
        for line in open(args.census, encoding="utf-8")
        if json.loads(line).get("type") == "record"
    }
    records = [json.loads(line) for line in open(args.rejected, encoding="utf-8")]
    print(f"records: {len(records)}")

    variants = ("B^T", "B", "B^-1", "B^-T")
    lines = []
    for record in records:
        source = census[record["ordinal"]]
        if source.get("basis_printed") is None:
            continue
        basis = [[Fraction(v) for v in row] for row in source["basis_printed"]]
        stored = [Fraction(v) for v in source["origin_machine"]]
        printed = [Fraction(v) for v in source["origin_printed"]]
        delta_origin = [printed[i] - stored[i] for i in range(3)]
        for variant in variants:
            matrix = transpose(basis) if variant in ("B^T", "B") else inverse(basis)
            if matrix is None:
                continue
            if variant.endswith("^-T"):
                matrix = transpose(matrix)
            if variant in ("B^-1", "B^-T"):
                pass
            moved = matvec(matrix, delta_origin)
            denominator = 1
            for value in moved:
                denominator = denominator * value.denominator // _gcd(denominator, value.denominator)
            numerator = [int(value * denominator) for value in moved]
            u, ud = probe_mod.rationalise(record["u"])
            for sign in (1, -1):
                lines.append(probe_mod.line(
                    f"{record['ordinal']}|{variant}|{sign}", record["ordinal"], u, ud,
                    (sign * numerator[0], sign * numerator[1], sign * numerator[2], denominator),
                ))

    outcome = probe_mod.probe(lines)
    solved = {}
    stats = {}
    for key, (ok, _payload) in outcome.items():
        if not ok:
            continue
        ordinal, variant, sign = key.split("|")
        solved.setdefault(int(ordinal), (variant, int(sign)))
        stats[variant] = stats.get(variant, 0) + 1
    print(f"solved {len(solved)}/{len(records)}")
    print("by convention:", stats)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump({str(k): list(v) for k, v in sorted(solved.items())}, handle, indent=1)


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


if __name__ == "__main__":
    sys.exit(main())
