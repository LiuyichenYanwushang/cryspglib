#!/usr/bin/env python3
"""Read the committed frozen table and ask the engine which entries it accepts."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import generate_subduction_settings as settings
import probe_census as pm


def parse(path):
    with open(path, encoding="utf-8") as handle:
        rows = settings.parse_committed(handle.read())
    return [
        (row["ordinal"], row["parent"], row["child"], row["setting"],
         row["setting_denominator"], row["child_shift"])
        for row in rows
    ]


def main(path, out):
    entries = parse(path)
    print("entries:", len(entries))
    lines = [pm.line(str(o), o, u, d, s) for o, _p, _c, u, d, s in entries]
    outcome = pm.probe(lines)
    rejected = [e for e in entries if not outcome.get(str(e[0]), (False,))[0]]
    print("accepted:", len(entries) - len(rejected), "rejected:", len(rejected))
    import collections
    reasons = collections.Counter(
        outcome[str(e[0])][1].split(" for subgroup")[0] for e in rejected)
    for reason, count in reasons.most_common(5):
        print(f"   {count:5d}  {reason}")
    by_parent = collections.Counter(e[1] for e in rejected)
    print("rejected by parent:", by_parent.most_common(12))
    json.dump([e[0] for e in rejected], open(out, "w"))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
