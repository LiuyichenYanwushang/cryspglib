#!/usr/bin/env python3
"""Derive the recorded-frame child shift for every record, oracle-free.

The machine table stores each isotropy record's origin in the frame the tables
were recorded in; the engine maps into the parent frame of `SG_DATA_HALL` and
uses that stored origin.  Where the official program prints a *different* origin
for the same record, the frozen convention has to move the subgroup by

    delta = (B^T)^-1 (o_printed - o_stored)

(`B` the printed basis, the stored origin converted through the parent primitive
basis).  Both inputs come from the pinned census, so this runs over the whole
table without a single oracle query -- it only needs the engine to accept the
result.

This matters even where `delta = 0` also validates: two placements can both be
valid subgroups, and only the printed origin is the official one.  Freezing the
stored-origin placement silently moves the trivial content onto other irreps of
the star (ordinal 27 is the minimal witness).
"""
import argparse
import json
import os
import sys
from fractions import Fraction

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import probe_census as probe_mod          # noqa: E402
import verify_isotropy_oracle as geometry  # noqa: E402
from derive_shift import inverse, matvec, transpose  # noqa: E402
from derive_shift_oracle import as_shift  # noqa: E402


def recorded_frame_shift(source):
    parent = geometry.CENTERING_LETTER[source["parent"]]
    primitive = geometry.PRIMITIVE_BASIS[parent]
    stored = [Fraction(v) for v in source["origin_machine"]]
    converted = [
        sum(stored[i] * Fraction(primitive[i][j]) for i in range(3)) for j in range(3)
    ]
    printed = [Fraction(v) for v in source["origin_printed"]]
    difference = [printed[i] - converted[i] for i in range(3)]
    if all(value == 0 for value in difference):
        return None
    basis = [[Fraction(v) for v in row] for row in source["basis_printed"]]
    matrix = inverse(transpose(basis))
    if matrix is None:
        return None
    return as_shift(matvec(matrix, difference))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    sources = {}
    for line in open(args.census, encoding="utf-8"):
        blob = json.loads(line)
        if blob.get("type") == "record":
            sources[blob["ordinal"]] = blob
    try:
        for line in open("empty_fixed.jsonl", encoding="utf-8"):
            blob = json.loads(line)
            sources.setdefault(blob["ordinal"], blob)
    except FileNotFoundError:
        pass

    candidates = {}
    for ordinal, source in sources.items():
        if not source.get("basis_printed") or not source.get("u"):
            continue
        shift = recorded_frame_shift(source)
        if shift is not None:
            candidates[ordinal] = shift
    print(f"records whose stored origin is not the printed one: {len(candidates)}")

    lines, meta = [], {}
    for ordinal, shift in sorted(candidates.items()):
        u, denominator = probe_mod.rationalise(sources[ordinal]["u"])
        tag = str(ordinal)
        lines.append(probe_mod.line(tag, ordinal, u, denominator,
                                    (shift[0], shift[1], shift[2], shift[3])))
        meta[tag] = (ordinal, shift)
    outcome = probe_mod.probe(lines)
    accepted = {meta[tag][0]: meta[tag][1] for tag in meta
                if outcome.get(tag, (False,))[0]}
    print(f"engine accepts the recorded-frame shift for {len(accepted)}/{len(meta)}")
    rejected = [meta[tag][0] for tag in meta if not outcome.get(tag, (False,))[0]]
    print(f"rejected ordinals: {len(rejected)} {rejected[:20]}")
    json.dump({"shifts": {str(k): list(v) for k, v in sorted(accepted.items())},
               "accepted": sorted(accepted), "rejected": sorted(rejected)},
              open(args.out, "w", encoding="utf-8"), indent=1, sort_keys=True)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
