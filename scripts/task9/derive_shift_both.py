#!/usr/bin/env python3
"""Two oracle-derived routes to the child shift, both offered to the engine.

The engine maps the subgroup's shipped Hall row into the parent frame of
`SG_DATA_HALL`, so the frozen `child_shift` has to undo whichever of two frame
differences the record carries:

* *origin choice*: the two ITA origin choices print the same cell at different
  origins, `s = o(OR1) - o(OR2)`;
* *recorded frame*: the pinned table's stored origin is expressed in the frame
  the tables were recorded in, which for part of the table is not the frame the
  engine maps into, `s = o_printed - o_stored` (the stored origin converted
  through the parent primitive basis).

Both are mapped into the child frame the same way, `delta = -(B^T)^-1 s`, and
each is probed separately, so the engine decides which difference this record
actually has instead of this file guessing.
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
import verify_isotropy_oracle as geometry  # noqa: E402
from derive_shift import inverse, matvec, transpose  # noqa: E402
from derive_shift_oracle import as_shift  # noqa: E402


def origin_choice_shift(record, source, legend, oracle_dir):
    row2 = generator.query_record(record["parent"], record["ml"], record["label"],
                                  legend, origin_choice=2, oracle_dir=oracle_dir)
    if row2["subgroup"] != source["child"] or row2["size"] != source["size"]:
        raise generator.GenerationError("origin-2 prints a different subgroup")
    printed = [Fraction(v) for v in source["origin_printed"]]
    difference = [printed[i] - row2["origin"][i] for i in range(3)]
    basis = [[Fraction(v) for v in row] for row in source["basis_printed"]]
    matrix = inverse(transpose(basis))
    if matrix is None:
        raise generator.GenerationError("singular basis")
    return as_shift([-v for v in matvec(matrix, difference)])


def recorded_frame_shift(record, source):
    parent = geometry.CENTERING_LETTER[record["parent"]]
    primitive = geometry.PRIMITIVE_BASIS[parent]
    stored = [Fraction(v) for v in source["origin_machine"]]
    converted = [
        sum(stored[i] * Fraction(primitive[i][j]) for i in range(3)) for j in range(3)
    ]
    printed = [Fraction(v) for v in source["origin_printed"]]
    difference = [printed[i] - converted[i] for i in range(3)]
    basis = [[Fraction(v) for v in row] for row in source["basis_printed"]]
    matrix = inverse(transpose(basis))
    if matrix is None:
        raise generator.GenerationError("singular basis")
    return as_shift(matvec(matrix, difference))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rejected", required=True, help="JSON list of ordinals")
    parser.add_argument("--census", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    wanted = set(json.load(open(args.rejected, encoding="utf-8")))
    records = {r["ordinal"]: r for r in generator.load_machine_records()}
    source = {}
    for line in open(args.census, encoding="utf-8"):
        blob = json.loads(line)
        if blob.get("type") == "record":
            source[blob["ordinal"]] = blob
    try:
        for line in open("empty_fixed.jsonl", encoding="utf-8"):
            blob = json.loads(line)
            source.setdefault(blob["ordinal"], blob)
    except FileNotFoundError:
        pass
    selected = [o for o in sorted(wanted)
                if o in records and source.get(o, {}).get("basis_printed")]
    print(f"records to derive: {len(selected)}")

    candidates = collections.defaultdict(dict)
    errors = collections.Counter()
    with tempfile.TemporaryDirectory(prefix="derive-both-") as temporary:
        oracle_dir, _ = census_mod.extract_oracle(temporary)
        legend = generator.load_point_op_legend()

        def work(ordinal):
            record = records[ordinal]
            found = {}
            try:
                found["origin-choice"] = origin_choice_shift(
                    record, source[ordinal], legend, oracle_dir)
            except Exception as error:  # noqa: BLE001 - reported per record
                found["origin-choice"] = f"{type(error).__name__}: {error}"[:100]
            try:
                found["recorded-frame"] = recorded_frame_shift(record, source[ordinal])
            except Exception as error:  # noqa: BLE001 - reported per record
                found["recorded-frame"] = f"{type(error).__name__}: {error}"[:100]
            return ordinal, found

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            for index, (ordinal, found) in enumerate(pool.map(work, selected)):
                candidates[ordinal] = found
                if (index + 1) % 250 == 0:
                    print(f"  {index + 1}/{len(selected)} derived", flush=True)

    lines, meta = [], {}
    for ordinal, found in candidates.items():
        blob = source[ordinal]
        u, denominator = probe_mod.rationalise(blob["u"])
        for route, shift in found.items():
            if isinstance(shift, str):
                errors[f"{route}: {shift}"] += 1
                continue
            tag = f"{ordinal}|{route}"
            lines.append(probe_mod.line(tag, ordinal, u, denominator,
                                        (shift[0], shift[1], shift[2], shift[3])))
            meta[tag] = (ordinal, route, shift)
    outcome = probe_mod.probe(lines)
    accepted = {}
    for tag, (ordinal, route, shift) in meta.items():
        if outcome.get(tag, (False,))[0]:
            accepted.setdefault(ordinal, (route, shift))
    print(f"engine accepts a derived shift for {len(accepted)}/{len(selected)}")
    print("by route:", collections.Counter(r for r, _ in accepted.values()))
    unresolved = [o for o in selected if o not in accepted]
    print(f"unresolved: {len(unresolved)} {unresolved[:20]}")
    for reason, count in errors.most_common(4):
        print(f"   {count:5d}  {reason}")
    json.dump({"shifts": {str(o): {"route": r, "shift": list(s)}
                          for o, (r, s) in sorted(accepted.items())},
               "accepted": sorted(accepted),
               "unresolved": unresolved},
              open(args.out, "w", encoding="utf-8"), indent=1, sort_keys=True)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
