#!/usr/bin/env python3
"""Pick the convention whose condensing irrep contains the trivial rep.

The embedding check accepts any valid placement, and two placements can both be
valid; what picks the official one is the defining property of an isotropy
subgroup -- the condensing irrep must subduce the trivial representation of the
subgroup.  `examples/probe_subduction_settings.rs` reports that multiplicity
(`trivial=<n>`) for every candidate it is handed, so the candidate list can be
searched for the conventions the table's own definition endorses.

Candidates: the recorded-frame shift, the ITA origin-choice shift, `delta = 0`,
and (only when none of those works) the full origin grid.
"""
import argparse
import collections
import concurrent.futures
import json
import os
import sys
import tempfile
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "scripts"))

import probe_census as probe_mod          # noqa: E402
import generate_subduction_settings as generator  # noqa: E402
import audit_subduction_settings as census_mod  # noqa: E402
import verify_isotropy_oracle as geometry  # noqa: E402
from derive_shift import inverse, matvec, transpose  # noqa: E402
from derive_shift_oracle import as_shift  # noqa: E402


def recorded_frame_shift(source):
    if not source.get("origin_machine") or not source.get("basis_printed"):
        return None
    primitive = geometry.PRIMITIVE_BASIS[geometry.CENTERING_LETTER[source["parent"]]]
    stored = [Fraction(v) for v in source["origin_machine"]]
    converted = [sum(stored[i] * Fraction(primitive[i][j]) for i in range(3))
                 for j in range(3)]
    printed = [Fraction(v) for v in source["origin_printed"]]
    difference = [printed[i] - converted[i] for i in range(3)]
    basis = [[Fraction(v) for v in row] for row in source["basis_printed"]]
    matrix = inverse(transpose(basis))
    if matrix is None:
        return None
    return as_shift(matvec(matrix, difference))


def origin_choice_shift(record, source, legend, oracle_dir):
    if not source.get("origin_printed") or not source.get("basis_printed"):
        return None
    row2 = generator.query_record(record["parent"], record["ml"], record["label"],
                                  legend, origin_choice=2, oracle_dir=oracle_dir)
    printed = [Fraction(v) for v in source["origin_printed"]]
    difference = [printed[i] - row2["origin"][i] for i in range(3)]
    basis = [[Fraction(v) for v in row] for row in source["basis_printed"]]
    matrix = inverse(transpose(basis))
    if matrix is None:
        return None
    return as_shift([-v for v in matvec(matrix, difference)])


def grid():
    shifts = []
    for denominator, bound in ((12, 6), (8, 4), (24, 12), (16, 8)):
        for x in range(-bound, bound + 1):
            for y in range(-bound, bound + 1):
                for z in range(-bound, bound + 1):
                    shifts.append((x, y, z, denominator))
    return shifts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ordinals", required=True)
    parser.add_argument("--census", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    wanted = json.load(open(args.ordinals, encoding="utf-8"))
    records = {r["ordinal"]: r for r in generator.load_machine_records()}
    sources = {}
    for line in open(args.census, encoding="utf-8"):
        blob = json.loads(line)
        if blob.get("type") == "record":
            sources[blob["ordinal"]] = blob
    try:
        for line in open("empty_fixed.jsonl", encoding="utf-8"):
            blob = json.loads(line)
            # The census leaves the legacy compound spellings without a U; the
            # compact-label pass supplies it, so it has to win there.
            if not sources.get(blob["ordinal"], {}).get("u"):
                sources[blob["ordinal"]] = blob
    except FileNotFoundError:
        pass

    resolved = {}
    with tempfile.TemporaryDirectory(prefix="trivial-census-") as temporary:
        oracle_dir, _ = census_mod.extract_oracle(temporary)
        legend = generator.load_point_op_legend()
        wide = grid()

        def candidate_shifts(ordinal):
            source = sources[ordinal]
            out = [(0, 0, 0, 1)]
            try:
                choice = origin_choice_shift(records[ordinal], source, legend, oracle_dir)
                if choice is not None:
                    out.append(tuple(choice))
            except Exception:  # noqa: BLE001 - a missing route is not fatal
                pass
            frame = recorded_frame_shift(source)
            if frame is not None:
                out.append(tuple(frame))
            return out

        def attempt(ordinal):
            source = sources[ordinal]
            u, denominator = probe_mod.rationalise(source["u"])
            for stage, shifts in (("principled", candidate_shifts(ordinal)),
                                  ("grid", wide)):
                lines = [probe_mod.line(f"{ordinal}|{index}", ordinal, u, denominator, shift)
                         for index, shift in enumerate(shifts)]
                outcome = probe_mod.probe(lines)
                hits = []
                for index, shift in enumerate(shifts):
                    result = outcome.get(f"{ordinal}|{index}")
                    if not result or not result[0]:
                        continue
                    trivial = result[1].split("trivial=")[-1]
                    if trivial.isdigit() and int(trivial) > 0:
                        hits.append((list(shift), int(trivial), stage))
                if hits:
                    return ordinal, hits
            return ordinal, []

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            for index, (ordinal, hits) in enumerate(pool.map(attempt, wanted)):
                resolved[ordinal] = hits
                print(f"  [{index + 1}/{len(wanted)}] ordinal {ordinal}: "
                      f"{len(hits)} conventions with trivial > 0"
                      + (f" (first {hits[0][0]} from the {hits[0][2]} stage)" if hits else ""),
                      flush=True)
    if args.out.endswith(".jsonl"):
        with open(args.out, "w", encoding="utf-8") as handle:
            for ordinal, hits in sorted(resolved.items()):
                handle.write(json.dumps({"ordinal": ordinal, "hits": hits},
                                        sort_keys=True) + "\n")
    else:
        json.dump({str(k): v for k, v in sorted(resolved.items())},
                  open(args.out, "w", encoding="utf-8"), indent=1, sort_keys=True)
    solved = sum(1 for hits in resolved.values() if hits)
    print(f"records with an endorsed convention: {solved}/{len(resolved)}")


if __name__ == "__main__":
    main()
