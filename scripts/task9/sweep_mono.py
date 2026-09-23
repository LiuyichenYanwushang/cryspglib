#!/usr/bin/env python3
"""Per-parent setting sweep for the monoclinic parents (SG 3-15)."""
import collections, json, os, sys, tempfile, zipfile
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, '.')
sys.path.insert(0, os.path.join(REPO, "scripts"))
import audit_subduction_settings as audit
import probe_census as pm

BASE = list(audit.ORACLE_COMMANDS)
VARIANTS = [("OR1", ["SET I ALL OR 1"]),
            ("OR2", ["SET I ALL OR 2"])]
for axis in ("b", "c"):
    VARIANTS.append((f"OR1 AXIS {axis}", ["SET I ALL OR 1", "SET I {sg} AXIS " + axis]))
    VARIANTS.append((f"OR2 AXIS {axis}", ["SET I ALL OR 2", "SET I {sg} AXIS " + axis]))
for cell in ("1", "2", "3"):
    VARIANTS.append((f"OR1 CELL {cell}", ["SET I ALL OR 1", "SET I {sg} CELL " + cell]))
    for axis in ("b", "c"):
        VARIANTS.append((f"OR1 AXIS {axis} CELL {cell}",
                         ["SET I ALL OR 1", f"SET I {{sg}} AXIS {axis}", f"SET I {{sg}} CELL {cell}"]))

def run(sg, irrep, head, root, dd, binary):
    cmds = ["PAGE 1000", "SC 250", *head, "VALUE PARENT {parent}", "VALUE IRREP {irrep}",
            "SHOW SUBGROUP", "SHOW BASIS", "SHOW ORIGIN", "SHOW SIZE", "SHOW DIRECTION",
            "DISPLAY ISOTROPY", "QUIT"]
    audit.ORACLE_COMMANDS = tuple(c.format(parent=sg, irrep=irrep) if "{" in c else c for c in cmds)
    try:
        return audit.run_query(binary, dd, root, {"parent": sg, "irrep": irrep})
    finally:
        audit.ORACLE_COMMANDS = tuple(BASE)

def main(parents, out_path):
    inventory = audit.load_inventory()
    records = audit.load_machine_records(inventory)
    by_parent = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in records:
        if r["parent"] in parents:
            by_parent[r["parent"]][r["irrep"]].append(r)
    chosen = {}
    with tempfile.TemporaryDirectory() as root:
        dd = os.path.join(root, "data"); os.makedirs(dd)
        with zipfile.ZipFile(os.path.join(REPO, "isotropy_subgroup", "iso.zip")) as z:
            z.extractall(dd)
        b = os.path.join(dd, "iso"); os.chmod(b, 0o755)
        wr = os.path.join(root, "w"); os.makedirs(wr)
        for sg in parents:
            best = None
            for name, template in VARIANTS:
                head = [c.format(sg=sg) for c in template]
                lines, meta = [], {}
                for irrep, rows in by_parent[sg].items():
                    out = run(sg, irrep, head, wr, dd, b)
                    if out["error"] is not None:
                        continue
                    status, ev, err = audit.classify_query(rows, out["rows"], sg)
                    if err is not None:
                        continue
                    for row, item in zip(rows, ev):
                        if item.get("u") is None:
                            continue
                        u, ud = pm.rationalise(item["u"])
                        tag = f"{row['ordinal']}"
                        lines.append(pm.line(tag, row["ordinal"], u, ud, (0, 0, 0, 1)))
                        meta[tag] = (row["ordinal"], row["child"], item["status"])
                outcome = pm.probe(lines)
                ok = [t for t in meta if outcome.get(t, (False,))[0]]
                print(f"  SG {sg:3d} {name:22s}: accepts {len(ok):3d}/{len(meta):3d}", flush=True)
                if best is None or len(ok) > best[1]:
                    best = (name, len(ok), len(meta), [t for t in ok],
                            {t: meta[t] for t in ok})
            chosen[sg] = {
                "variant": best[0], "accepted": best[1], "total": best[2],
                "accepted_ordinals": sorted(v[0] for v in best[4].values()),
                "accepted_detail": {v[0]: {"child": v[1], "status": v[2]}
                                    for v in best[4].values()},
            }
            print(f"SG {sg:3d} best {best[0]} {best[1]}/{best[2]}", flush=True)
    json.dump(chosen, open(out_path, "w"), indent=1, sort_keys=True)

if __name__ == "__main__":
    parents = [int(x) for x in sys.argv[1].split(",")]
    main(parents, sys.argv[2])
