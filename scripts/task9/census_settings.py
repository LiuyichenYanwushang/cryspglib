#!/usr/bin/env python3
"""Setting-aware census: choose, per parent, the oracle setting the pinned
isotropy table was recorded in, then classify every record under it.

The rule is *evidence driven*: the ITA setting of ``SG_DATA_HALL[parent]`` fixes
the origin/cell/axis choice the machine tables were generated in, and a record
may only be accepted when the official basis/origin reproduce the stored ones
exactly (``U`` integral unimodular and the origin equal after conversion through
the parent primitive basis).

Output records use the same schema as ``audit_subduction_settings.py`` so the
existing emitters and analysers can consume them.
"""
import argparse
import collections
import concurrent.futures
import json
import os
import re
import sys
import tempfile
import zipfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import audit_subduction_settings as audit  # noqa: E402

VARIANT_COMMANDS = [
    ("OR1", ["SET I ALL OR 1"]),
    ("OR2", ["SET I ALL OR 2"]),
    ("OR1 AXIS b", ["SET I ALL OR 1", "SET I {sg} AXIS b"]),
    ("OR1 AXIS c", ["SET I ALL OR 1", "SET I {sg} AXIS c"]),
    ("OR2 AXIS b", ["SET I ALL OR 2", "SET I {sg} AXIS b"]),
    ("OR2 AXIS c", ["SET I ALL OR 2", "SET I {sg} AXIS c"]),
]
for _cell in ("1", "2", "3"):
    VARIANT_COMMANDS.append(
        (f"OR1 CELL {_cell}", ["SET I ALL OR 1", "SET I {sg} CELL " + _cell])
    )
for _cell in ("1", "2", "3"):
    VARIANT_COMMANDS.append(
        (f"OR2 CELL {_cell}", ["SET I ALL OR 2", "SET I {sg} CELL " + _cell])
    )
for _cell in ("1", "2", "3"):
    for _axis in ("b", "c"):
        VARIANT_COMMANDS.append((
            f"OR1 AXIS {_axis} CELL {_cell}",
            ["SET I ALL OR 1", f"SET I {{sg}} AXIS {_axis}", f"SET I {{sg}} CELL {_cell}"],
        ))

BASE = list(audit.ORACLE_COMMANDS)


def load_data_hall():
    text = open(os.path.join(REPO, "src", "spg_database.rs"), encoding="utf-8").read()
    body = text.split("pub static SPACEGROUP_TYPES", 1)[1]
    rows = []
    for chunk in body.split("RawSpacegroupType {")[1:]:
        number = re.search(r"number: (\d+)", chunk)
        choice = re.search(r'choice: "([^"]*)"', chunk)
        rows.append((int(number.group(1)), choice.group(1).strip()))
    provenance = json.load(
        open(os.path.join(REPO, "scripts", "data", "iso_irrep_data_hall_v1.json"))
    )
    return {r["spacegroup"]: rows[r["data_hall"]][1] for r in provenance["spacegroups"]}


def setting_head(sg, name):
    for candidate, commands in VARIANT_COMMANDS:
        if candidate == name:
            return [c.format(sg=sg) for c in commands]
    raise KeyError(name)


def run_variant(sg, irrep, head, work_root, data_dir, binary):
    commands = [
        "PAGE 1000", "SC 250", *head,
        "VALUE PARENT {parent}", "VALUE IRREP {irrep}",
        "SHOW SUBGROUP", "SHOW BASIS", "SHOW ORIGIN", "SHOW SIZE",
        "SHOW DIRECTION", "DISPLAY ISOTROPY", "QUIT",
    ]
    audit.ORACLE_COMMANDS = tuple(
        c.format(parent=sg, irrep=irrep) if "{" in c else c for c in commands
    )
    try:
        return audit.run_query(binary, data_dir, work_root, {"parent": sg, "irrep": irrep})
    finally:
        audit.ORACLE_COMMANDS = tuple(BASE)


def classify_parent(sg, by_irrep, head, work_root, data_dir, binary):
    """Return (evidence_by_ordinal, counts)."""
    evidence = {}
    counts = collections.Counter()
    for irrep, rows in by_irrep.items():
        outcome = run_variant(sg, irrep, head, work_root, data_dir, binary)
        if outcome["error"] is not None:
            counts["query_error"] += len(rows)
            for row in rows:
                blank = dict(audit.blank_evidence())
                blank.update(status="query_error", detail=outcome["error"]["kind"])
                evidence[row["ordinal"]] = blank
            continue
        status, per_row, error = audit.classify_query(rows, outcome["rows"], sg)
        if error is not None:
            counts["query_error"] += len(rows)
            for row in rows:
                blank = dict(audit.blank_evidence())
                blank.update(status="query_error", detail=error["kind"])
                evidence[row["ordinal"]] = blank
            continue
        for row, item in zip(rows, per_row):
            evidence[row["ordinal"]] = item
            counts[item["status"]] += 1
    return evidence, counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--parents", default="problem")
    parser.add_argument("--baseline", help="census to take the problem-parent list from")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--report", action="store_true")
    args = parser.parse_args()

    inventory = audit.load_inventory()
    records = audit.load_machine_records(inventory)
    by_parent = collections.defaultdict(list)
    for record in records:
        by_parent[record["parent"]].append(record)
    if args.parents == "problem":
        baseline = [
            json.loads(line)
            for line in open(args.baseline, encoding="utf-8")
            if json.loads(line).get("type") == "record"
        ]
        bad = {r["parent"] for r in baseline if r["status"] != "candidate"}
        parents = sorted(bad)
    else:
        parents = [int(x) for x in args.parents.split(",")]
    print(f"parents in scope: {len(parents)}", flush=True)

    data_hall_choice = load_data_hall()
    per_parent = collections.defaultdict(lambda: collections.defaultdict(list))
    for parent in parents:
        for record in by_parent[parent]:
            per_parent[parent][record["irrep"]].append(record)

    chosen = {}
    with tempfile.TemporaryDirectory() as root:
        data_dir = os.path.join(root, "data")
        os.makedirs(data_dir)
        with zipfile.ZipFile(os.path.join(REPO, "isotropy_subgroup", "iso.zip")) as z:
            z.extractall(data_dir)
        binary = os.path.join(data_dir, "iso")
        os.chmod(binary, 0o755)
        work_root = os.path.join(root, "work")
        os.makedirs(work_root)
        for index, parent in enumerate(parents):
            preferred = f"OR{data_hall_choice[parent]}" if data_hall_choice[parent] in ("1", "2") else "OR1"
            order = sorted(
                VARIANT_COMMANDS,
                key=lambda item: (item[0] != preferred, item[0] != "OR1"),
            )
            results = {}
            for name, _ in order:
                head = setting_head(parent, name)
                evidence, counts = classify_parent(
                    parent, per_parent[parent], head, work_root, data_dir, binary
                )
                results[name] = (evidence, counts)
            best = max(
                results,
                key=lambda name: (results[name][1]["candidate"], -len(results[name][1])),
            )
            chosen[parent] = best
            if args.report or True:
                summary = {
                    name: results[name][1]["candidate"] for name, _ in order
                }
                print(
                    f"  [{index + 1}/{len(parents)}] SG {parent:3d} "
                    f"data-Hall choice {data_hall_choice[parent]!r} "
                    f"best={best} candidates={summary}",
                    flush=True,
                )
            per_parent[parent] = results[best][0]

    out_records = []
    for record in records:
        parent = record["parent"]
        if parent not in parents:
            continue
        evidence = per_parent[parent].get(record["ordinal"])
        item = dict(record)
        if evidence is not None:
            item.update(evidence)
        out_records.append(item)
    out_records.sort(key=lambda r: r["ordinal"])
    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(json.dumps({"type": "manifest", "kind": "setting-aware-census"},
                                sort_keys=True) + "\n")
        for item in out_records:
            handle.write(json.dumps(item, sort_keys=True) + "\n")
        handle.write(json.dumps({
            "type": "summary",
            "chosen": {str(k): v for k, v in sorted(chosen.items())},
        }, sort_keys=True) + "\n")
    print("wrote", args.output, len(out_records), "records")
    print("chosen settings:", json.dumps(chosen, sort_keys=True))


if __name__ == "__main__":
    main()
