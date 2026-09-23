#!/usr/bin/env python3
"""Derive a child shift from `SHOW ELEMENTS` for records whose origin differs.

The printed-origin route (`derive_shift_oracle.py`) needs the two ITA origin
choices to print the *same* cell, which fails for the monoclinic records whose
stored basis/origin do not match the program's default frame.  This route works
from the operations instead:

1. query the record at both ITA origin choices, keeping `SHOW ELEMENTS`,
2. expand those operations into the subgroup's own conventional cell,
3. select the printed setting whose full expanded operation set equals the
   subgroup's shipped Hall row (whole-set equality, never one operation),
4. `delta = T^-1 (o_setting - o_recorded)`.

A record whose Hall row matches no printed setting, or whose matching settings
disagree on the shift, is reported instead of guessed.
"""
import argparse
import collections
import concurrent.futures
import json
import os
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "scripts"))

import probe_census as probe_mod  # noqa: E402
import generate_subduction_settings as generator  # noqa: E402
import audit_subduction_settings as census_mod  # noqa: E402


def derive_one(record, legend, data_hall, hall_table, oracle_dir):
    parent = record["parent"]
    ml = record["ml"]
    label = record["label"]
    row1 = generator.query_record(parent, ml, label, legend, oracle_dir=oracle_dir)
    row2 = generator.query_record(parent, ml, label, legend, origin_choice=2,
                                  oracle_dir=oracle_dir)
    if row2["subgroup"] != row1["subgroup"] or row2["size"] != row1["size"]:
        raise generator.GenerationError("the two settings print different subgroups")
    expansions = {}
    for choice, chosen in ((1, row1), (2, row2)):
        expansions[choice] = generator.child_frame_operations(
            record, chosen["basis"], chosen["origin"], chosen["elements"])
    origins = {1: row1["origin"], 2: row2["origin"]}
    shift = generator.origin_choice_shift(
        record, row1["basis"], expansions, origins,
        hall_table[data_hall[record["child"]]]["operations"])
    return generator.encode_shift(shift)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rejected", required=True, help="JSON list of ordinals")
    parser.add_argument("--census", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    wanted = set(json.load(open(args.rejected, encoding="utf-8")))
    records = {r["ordinal"]: r for r in generator.load_machine_records()}
    selected = [records[o] for o in sorted(wanted) if o in records]
    print(f"records to derive: {len(selected)}")
    data_hall = generator.load_data_hall()
    hall_table = generator.load_hall_operations()
    table = {}
    for line in open(args.census, encoding="utf-8"):
        blob = json.loads(line)
        if blob.get("type") == "record":
            table[blob["ordinal"]] = blob

    derived = {}
    failures = collections.Counter()
    with tempfile.TemporaryDirectory(prefix="derive-elements-") as temporary:
        oracle_dir, _ = census_mod.extract_oracle(temporary)
        legend = generator.load_point_op_legend()
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(derive_one, r, legend, data_hall, hall_table, oracle_dir): r
                for r in selected
            }
            for index, future in enumerate(concurrent.futures.as_completed(futures)):
                record = futures[future]
                try:
                    derived[record["ordinal"]] = list(future.result())
                except Exception as error:  # noqa: BLE001 - reported per record
                    failures[f"{type(error).__name__}: {error}"[:110]] += 1
                if (index + 1) % 25 == 0:
                    print(f"  {index + 1}/{len(selected)}: derived {len(derived)}, "
                          f"failed {sum(failures.values())}", flush=True)
    print(f"derived {len(derived)}/{len(selected)}")
    for reason, count in failures.most_common(6):
        print(f"   {count:5d}  {reason}")

    lines = []
    for ordinal, shift in derived.items():
        source = table[ordinal]
        u, denominator = probe_mod.rationalise(source["u"])
        lines.append(probe_mod.line(str(ordinal), ordinal, u, denominator,
                                    (shift[0], shift[1], shift[2], shift[3])))
    outcome = probe_mod.probe(lines)
    accepted = [o for o in derived if outcome.get(str(o), (False,))[0]]
    print(f"engine accepts {len(accepted)}/{len(derived)}")
    for reason, count in collections.Counter(
        outcome[str(o)][1].split(" for subgroup")[0] for o in derived
        if not outcome.get(str(o), (False,))[0]
    ).most_common(4):
        print(f"   {count:5d}  {reason}")
    json.dump({"shifts": {str(k): v for k, v in sorted(derived.items())},
               "accepted": sorted(accepted),
               "failures": dict(failures)},
              open(args.out, "w", encoding="utf-8"), indent=1, sort_keys=True)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
