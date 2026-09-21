#!/usr/bin/env python3
"""Analyse the task-9 census: which parents still fail the OR-1 exact gate."""
import collections
import json
import sys
from fractions import Fraction


def load(path):
    records = []
    summary = manifest = None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            blob = json.loads(line)
            if blob.get("type") == "record":
                records.append(blob)
            elif blob.get("type") == "summary":
                summary = blob
            else:
                manifest = blob
    return manifest, records, summary


def main(path):
    manifest, records, summary = load(path)
    print("manifest:", json.dumps(manifest, sort_keys=True)[:200])
    print("records:", len(records))
    if summary:
        print("summary:", json.dumps(summary["classification"]))
        print("detail_counts:", json.dumps(summary["detail_counts"], indent=1))
        print("queries:", summary["queries_total"], "nonempty", summary["queries_nonempty"],
              "empty", summary["queries_empty"], "error", summary["queries_error"])
    status = collections.Counter(r["status"] for r in records)
    print("status:", dict(status))
    by_parent = collections.Counter(r["parent"] for r in records if r["status"] != "candidate")
    print("\nparents with non-candidates:", len(by_parent))
    for parent, count in sorted(by_parent.items()):
        rows = [r for r in records if r["parent"] == parent]
        bad = [r for r in rows if r["status"] != "candidate"]
        kinds = collections.Counter(r["detail"] for r in bad)
        dets = collections.Counter(r["u_det"] for r in bad)
        deltas = collections.Counter(
            tuple(r["origin_delta"]) if r["origin_delta"] else None for r in bad
        )
        print(f"  SG {parent:3d}: {len(bad):4d}/{len(rows):4d} bad  {dict(kinds)}  det={dict(dets)}")
        for delta, n in deltas.most_common(3):
            print(f"        delta {delta} x{n}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
