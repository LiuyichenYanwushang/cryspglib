#!/usr/bin/env python3
"""Per-record child-shift search for the conventions the engine rejects.

The child shift is a per-ordinal convention (the isotropy table, the child Hall
row and the parent's data-Hall frame can differ by an ITA origin choice in a way
that depends on the record's own cell), so it is searched per record and the
result is accepted only when the engine's own validation accepts it.
"""
import argparse
import json
from fractions import Fraction
from math import gcd

import probe_census as probe_mod


def grid(wide=False):
    """Small shifts first, then the full ITA origin grid.

    The origin denominators of the pinned tables live on the 1/12 and 1/8 grids
    (and their 1/16, 1/24 refinements); the narrow pass covers every correction
    whose components stay within half a cell, and the wide pass is only needed
    where a full-cell origin choice is involved.
    """
    bounds = ((12, 12), (8, 8), (16, 8), (24, 12)) if wide else ((12, 6), (8, 4))
    shifts = []
    for denominator, bound in bounds:
        for x in range(-bound, bound + 1):
            for y in range(-bound, bound + 1):
                for z in range(-bound, bound + 1):
                    shifts.append((x, y, z, denominator))
    return shifts


GRID = grid(wide=__import__("os").environ.get("WIDE_GRID") == "1")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("rejected")
    parser.add_argument("--out", required=True)
    parser.add_argument("--per-chunk", type=int, default=60)
    args = parser.parse_args()

    records = [json.loads(text) for text in open(args.rejected, encoding="utf-8")]
    print(f"records to solve: {len(records)}; grid {len(GRID)} shifts each")
    solved = {}
    unresolved = []
    for start in range(0, len(records), args.per_chunk):
        chunk = records[start:start + args.per_chunk]
        lines = []
        for position, record in enumerate(chunk):
            u, denominator = probe_mod.rationalise(record["u"])
            for slot, shift in enumerate(GRID):
                lines.append(probe_mod.line(
                    f"{position}_{slot}", record["ordinal"], u,
                    denominator, shift,
                ))
        outcome = probe_mod.probe(lines)
        for position, record in enumerate(chunk):
            hit = None
            for slot, shift in enumerate(GRID):
                if outcome.get(f"{position}_{slot}", (False,))[0]:
                    hit = shift
                    break
            if hit is None:
                unresolved.append(record["ordinal"])
            else:
                solved[record["ordinal"]] = list(hit)
        print(f"  {start + len(chunk)}/{len(records)}: solved {len(solved)}, "
              f"unresolved {len(unresolved)}", flush=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump({"solved": {str(k): v for k, v in sorted(solved.items())},
                   "unresolved": unresolved}, handle, indent=1, sort_keys=True)
    print(f"solved {len(solved)}, unresolved {len(unresolved)}")
    print("unresolved ordinals:", unresolved[:40])


if __name__ == "__main__":
    main()
