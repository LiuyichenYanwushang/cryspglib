#!/usr/bin/env python3
"""For one parent, derive U under a setting variant and ask the engine."""
import collections, json, os, sys, tempfile, zipfile
sys.path.insert(0, '.')
sys.path.insert(0, '/home/liuyichen/TB_rs/cryspglib/scripts')
import audit_subduction_settings as audit
import probe_census as pm

BASE = list(audit.ORACLE_COMMANDS)

def run(sg, irrep, head, root, dd, binary):
    cmds = ["PAGE 1000", "SC 250", *head, "VALUE PARENT {parent}", "VALUE IRREP {irrep}",
            "SHOW SUBGROUP", "SHOW BASIS", "SHOW ORIGIN", "SHOW SIZE", "SHOW DIRECTION",
            "DISPLAY ISOTROPY", "QUIT"]
    audit.ORACLE_COMMANDS = tuple(c.format(parent=sg, irrep=irrep) if "{" in c else c for c in cmds)
    try:
        return audit.run_query(binary, dd, root, {"parent": sg, "irrep": irrep})
    finally:
        audit.ORACLE_COMMANDS = tuple(BASE)

def main(sg, heads):
    inventory = audit.load_inventory()
    records = [r for r in audit.load_machine_records(inventory) if r["parent"] == sg]
    by_irrep = collections.defaultdict(list)
    for r in records:
        by_irrep[r["irrep"]].append(r)
    with tempfile.TemporaryDirectory() as root:
        dd = os.path.join(root, "data"); os.makedirs(dd)
        with zipfile.ZipFile('/home/liuyichen/TB_rs/cryspglib/isotropy_subgroup/iso.zip') as z:
            z.extractall(dd)
        b = os.path.join(dd, "iso"); os.chmod(b, 0o755)
        wr = os.path.join(root, "w"); os.makedirs(wr)
        for head in heads:
            lines = []
            meta = {}
            counts = collections.Counter()
            for irrep, rows in by_irrep.items():
                out = run(sg, irrep, head, wr, dd, b)
                if out["error"] is not None:
                    counts["query_error"] += len(rows); continue
                status, ev, err = audit.classify_query(rows, out["rows"], sg)
                if err is not None:
                    counts["classify_error"] += len(rows); continue
                for row, item in zip(rows, ev):
                    counts[item["status"]] += 1
                    if item.get("u") is None:
                        continue
                    u, ud = pm.rationalise(item["u"])
                    tag = f"{row['ordinal']}"
                    lines.append(pm.line(tag, row["ordinal"], u, ud, (0, 0, 0, 1)))
                    meta[tag] = (row["ordinal"], item["status"])
            outcome = pm.probe(lines)
            ok = sum(1 for t in meta if outcome.get(t, (False,))[0])
            print(f"SG {sg} {head}: census {dict(counts)} -> engine accepts {ok}/{len(meta)}")
            for t, (o, st) in sorted(meta.items(), key=lambda kv: kv[1][0]):
                if outcome.get(t, (False,))[0] and st != "candidate":
                    print(f"    non-candidate accepted: ordinal {o} ({st})")

if __name__ == "__main__":
    main(int(sys.argv[1]), [h.split(";") for h in sys.argv[2:]])
