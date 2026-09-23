#!/usr/bin/env python3
"""Resolve empty official tables by trying the compact little-irrep spellings.

112 parent irreps of the pinned machine table carry the legacy Miller-Love
compound spelling (``W1W1``), which the official program does not accept: it
answers with an empty table.  The same irrep is addressable through the compact
``little_irr_full_label`` spelling (``W1WA1``).  Which compact label belongs to
which parent is not recorded in the machine table, so the mapping is *proved*
per record: a candidate spelling is accepted only when the printed table matches
the stored rows exactly (row count, direction labels, subgroup number and size).
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

BASE = list(audit.ORACLE_COMMANDS)


def compact_labels():
    """Every label of ``little_irr_full_label``, in file order, deduplicated."""
    with zipfile.ZipFile(os.path.join(REPO, "isotropy_subgroup", "iso.zip")) as archive:
        text = archive.read("data_little.txt").decode("latin-1")
    start = text.index("little_irr_full_label")
    tail = text[start + len("little_irr_full_label"):]
    stop = re.search(r"\n[a-z_]+[a-z_0-9]*\n", tail)
    segment = tail[: stop.start()] if stop else tail
    labels = re.findall(r'"([^"]*)"', segment)
    seen = []
    for label in labels:
        stripped = label.strip()
        if stripped and stripped not in seen:
            seen.append(stripped)
    return seen


def run_query(sg, irrep, root, data_dir, binary):
    commands = [
        "PAGE 1000", "SC 250", "SET I ALL OR 1",
        "VALUE PARENT {parent}", "VALUE IRREP {irrep}",
        "SHOW SUBGROUP", "SHOW BASIS", "SHOW ORIGIN", "SHOW SIZE",
        "SHOW DIRECTION", "DISPLAY ISOTROPY", "QUIT",
    ]
    audit.ORACLE_COMMANDS = tuple(
        c.format(parent=sg, irrep=irrep) if "{" in c else c for c in commands
    )
    try:
        return audit.run_query(binary, data_dir, root, {"parent": sg, "irrep": irrep})
    finally:
        audit.ORACLE_COMMANDS = tuple(BASE)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    inventory = audit.load_inventory()
    records = audit.load_machine_records(inventory)
    by_irrep = collections.defaultdict(list)
    for record in records:
        by_irrep[(record["parent"], record["irrep"])].append(record)
    census = [
        json.loads(line) for line in open(args.census, encoding="utf-8")
        if json.loads(line).get("type") == "record"
    ]
    empty = sorted({
        (row["parent"], row["irrep"]) for row in census
        if row["status"] == "unreturned_empty"
    })
    print(f"empty-table parent irreps: {len(empty)}")

    labels = compact_labels()
    print(f"compact labels available: {len(labels)}")
    by_prefix = collections.defaultdict(list)
    for label in labels:
        by_prefix[label[:2]].append(label)

    work = []
    for parent, irrep in empty:
        prefix = irrep[:2]
        # A compound spelling repeats its factor, so the compact spelling keeps
        # the same k-point prefix; the record match below is the real gate.
        for candidate in by_prefix.get(prefix, []):
            work.append((parent, irrep, candidate))
    counts = collections.Counter(len(by_prefix.get(i[:2], [])) for _, i in empty)
    print(f"candidate queries: {len(work)}  (per irrep: {dict(counts)})")

    results = {}
    with tempfile.TemporaryDirectory() as root:
        data_dir = os.path.join(root, "data")
        os.makedirs(data_dir)
        with zipfile.ZipFile(os.path.join(REPO, "isotropy_subgroup", "iso.zip")) as z:
            z.extractall(data_dir)
        binary = os.path.join(data_dir, "iso")
        os.chmod(binary, 0o755)
        work_root = os.path.join(root, "work")
        os.makedirs(work_root)

        def attempt(item):
            parent, irrep, candidate = item
            outcome = run_query(parent, candidate, work_root, data_dir, binary)
            if outcome["error"] is not None:
                return item, None
            rows = by_irrep[(parent, irrep)]
            status, evidence, error = audit.classify_query(rows, outcome["rows"], parent)
            # An empty answer is not a match: the candidate spelling has to
            # print this record's rows, with the same labels, subgroup and size.
            if error is not None or status != "nonempty":
                return item, None
            if any(item_evidence.get("basis_printed") is None for item_evidence in evidence):
                return item, None
            return item, evidence

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            for index, (item, evidence) in enumerate(pool.map(attempt, work)):
                if evidence is None:
                    continue
                key = (item[0], item[1])
                results.setdefault(key, []).append((item[2], evidence))
                if (index + 1) % 500 == 0:
                    print(f"  {index + 1}/{len(work)} attempted, "
                          f"{len(results)} irreps matched", flush=True)

    # Several compact spellings can print the same table (the official program
    # accepts both a compound label and its factors); what has to be unique is
    # the *geometry*, because that is what the frozen convention carries.
    resolved = {}
    ambiguous = []
    for key, hits in sorted(results.items()):
        # The stored basis/origin pin the recorded frame: prefer the spellings
        # whose printed geometry is exactly the one the machine table recorded.
        exact = [hit for hit in hits
                 if all(item["status"] == "candidate" for item in hit[1])]
        pool = exact if exact else hits
        # One signature per candidate spelling: the *whole* set of printed rows
        # this irrep's records get.  Comparing per record would confuse two
        # directions of the same table with two competing conventions.
        signatures = {}
        for name, evidence in pool:
            signature = tuple(
                (
                    tuple(tuple(str(v) for v in row) for row in item["basis_printed"]),
                    tuple(str(v) for v in item["origin_printed"]),
                )
                for item in evidence
            )
            signatures.setdefault(signature, name)
        if len(signatures) > 1:
            ambiguous.append((key, [h[0] for h in pool][:8]))
            continue
        resolved[key] = pool[0]
    unresolved = [key for key in empty if key not in resolved]
    print(f"resolved {len(resolved)}/{len(empty)} empty-table irreps "
          f"(geometry-unique)")
    print(f"geometry-ambiguous: {len(ambiguous)}")
    for key, names in ambiguous[:10]:
        print(f"   {key} -> {names}")
    print(f"unresolved: {unresolved[:20]}")

    payload = []
    for (parent, irrep), (candidate, evidence) in sorted(resolved.items()):
        for record, item in zip(by_irrep[(parent, irrep)], evidence):
            entry = dict(record)
            entry.update({k: v for k, v in item.items() if k != "basis_printed"
                          and k != "origin_printed"})
            entry["basis_printed"] = [[str(v) for v in row] for row in item["basis_printed"]]
            entry["origin_printed"] = [str(v) for v in item["origin_printed"]]
            entry["oracle_label"] = candidate
            payload.append(entry)
    with open(args.out, "w", encoding="utf-8") as handle:
        for entry in payload:
            handle.write(json.dumps(entry, sort_keys=True, default=str) + "\n")
    print(f"wrote {len(payload)} record evidence rows to {args.out}")


if __name__ == "__main__":
    main()
