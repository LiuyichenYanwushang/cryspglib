#!/usr/bin/env python3
"""Derive each record's child shift from the official program's two ITA settings.

The engine maps the subgroup's shipped Hall row into the parent frame of
``SG_DATA_HALL``.  Where the isotropy table and the parent's data-Hall frame use
different ITA origin choices, the stored origin is expressed in the other frame
and the frozen ``child_shift`` has to undo exactly that difference:

    s     = o(origin choice 1) - o(origin choice 2)     [parent frame]
    delta = -(B^T)^-1 s                                 [child frame]

with `B` the basis the program prints for the record.  Both inputs come from the
official program; the engine's validation is still the gate, and a record whose
two settings print different cells (so the origin difference is not a pure frame
shift) is reported instead of guessed.
"""
import argparse
import collections
import concurrent.futures
import json
import os
import sys
import tempfile
from fractions import Fraction

REPO = "/home/liuyichen/TB_rs/cryspglib"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "scripts"))

import probe_census as probe_mod  # noqa: E402
import generate_subduction_settings as generator  # noqa: E402
import audit_subduction_settings as census_mod  # noqa: E402


def gcd(left, right):
    while right:
        left, right = right, left % right
    return left


def transpose(matrix):
    return [[matrix[j][i] for j in range(3)] for i in range(3)]


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


def as_shift(vector):
    denominator = 1
    for value in vector:
        denominator = denominator * value.denominator // gcd(denominator, value.denominator)
    numerator = [int(value * denominator) for value in vector]
    divisor = gcd(gcd(abs(numerator[0]), abs(numerator[1])),
                  gcd(abs(numerator[2]), denominator))
    if divisor == 0:
        return (0, 0, 0, 1)
    return (numerator[0] // divisor, numerator[1] // divisor,
            numerator[2] // divisor, denominator // divisor)


def derive_one(record, source, legend, oracle_dir):
    """The candidate shift for one record, or a reason it cannot be derived."""
    row2 = generator.query_record(record["parent"], record["irrep"],
                                  record["label"], legend,
                                  origin_choice=2, oracle_dir=oracle_dir)
    basis = [[Fraction(str(value)) for value in row] for row in source["basis_printed"]]
    printed = [Fraction(value) for value in source["origin_printed"]]
    # The two origin choices may re-express the subgroup cell (a different but
    # equivalent basis); only the origin difference is used here, and the engine
    # still validates the resulting convention.
    if row2["subgroup"] != source["child"] or row2["size"] != source["size"]:
        return None, "origin-2 prints a different subgroup"
    shift = [printed[i] - row2["origin"][i] for i in range(3)]
    if all(value == 0 for value in shift):
        return (0, 0, 0, 1), None
    matrix = inverse(transpose(basis))
    if matrix is None:
        return None, "singular basis"
    return as_shift([-value for value in matvec(matrix, shift)]), None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", required=True)
    parser.add_argument("--rejected", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    census = {}
    for line in open(args.census, encoding="utf-8"):
        blob = json.loads(line)
        if blob.get("type") == "record":
            census[blob["ordinal"]] = blob
    records = [json.loads(line) for line in open(args.rejected, encoding="utf-8")]
    print(f"records to derive: {len(records)}")

    derived = {}
    failures = collections.Counter()
    with tempfile.TemporaryDirectory(prefix="derive-shift-") as temporary:
        oracle_dir, _ = census_mod.extract_oracle(temporary)
        legend = generator.load_point_op_legend()

        def work(record):
            try:
                return record, derive_one(record, census[record["ordinal"]], legend,
                                          oracle_dir)
            except Exception as error:  # noqa: BLE001 - reported per record
                return record, (None, f"{type(error).__name__}: {error}"[:120])

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            for index, (record, (shift, reason)) in enumerate(pool.map(work, records)):
                if shift is None:
                    failures[reason] += 1
                else:
                    derived[record["ordinal"]] = shift
                if (index + 1) % 500 == 0:
                    print(f"  {index + 1}/{len(records)}: derived {len(derived)}, "
                          f"failed {sum(failures.values())}", flush=True)
    print(f"derived {len(derived)}/{len(records)}")
    for reason, count in failures.most_common(6):
        print(f"   {count:5d}  {reason}")

    lines = []
    for record in records:
        if record["ordinal"] not in derived:
            continue
        u, denominator = probe_mod.rationalise(record["u"])
        lines.append(probe_mod.line(record["ordinal"], record["ordinal"], u,
                                    denominator, derived[record["ordinal"]]))
    outcome = probe_mod.probe(lines)
    accepted = [o for o in derived if outcome.get(str(o), (False,))[0]]
    print(f"engine accepts {len(accepted)}/{len(derived)} of the derived shifts")
    rejects = collections.Counter(
        outcome[str(o)][1].split(" for subgroup")[0] for o in derived
        if not outcome.get(str(o), (False,))[0]
    )
    for reason, count in rejects.most_common(4):
        print(f"   {count:5d}  {reason}")
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump({"shifts": {str(k): list(v) for k, v in sorted(derived.items())},
                   "accepted": sorted(accepted),
                   "failures": dict(failures)}, handle, indent=1, sort_keys=True)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
