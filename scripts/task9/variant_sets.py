#!/usr/bin/env python3
"""Print the accepted-ordinal set of each setting variant for one parent."""
import collections, os, sys, tempfile, zipfile
sys.path.insert(0, '.')
sys.path.insert(0, '/home/liuyichen/TB_rs/cryspglib/scripts')
import audit_subduction_settings as audit
import probe_census as pm
from sweep_mono import VARIANTS, run

def main(sg):
    inventory = audit.load_inventory()
    records = [r for r in audit.load_machine_records(inventory) if r["parent"] == sg]
    by_irrep = collections.defaultdict(list)
    for r in records:
        by_irrep[r["irrep"]].append(r)
    sets = {}
    with tempfile.TemporaryDirectory() as root:
        dd = os.path.join(root, "data"); os.makedirs(dd)
        with zipfile.ZipFile('/home/liuyichen/TB_rs/cryspglib/isotropy_subgroup/iso.zip') as z:
            z.extractall(dd)
        b = os.path.join(dd, "iso"); os.chmod(b, 0o755)
        wr = os.path.join(root, "w"); os.makedirs(wr)
        for name, template in VARIANTS:
            head = [c.format(sg=sg) for c in template]
            lines, meta = [], {}
            for irrep, rows in by_irrep.items():
                out = run(sg, irrep, head, wr, dd, b)
                if out["error"] is not None:
                    continue
                _status, ev, err = audit.classify_query(rows, out["rows"], sg)
                if err is not None:
                    continue
                for row, item in zip(rows, ev):
                    if item.get("u") is None:
                        continue
                    u, ud = pm.rationalise(item["u"])
                    lines.append(pm.line(str(row["ordinal"]), row["ordinal"], u, ud, (0, 0, 0, 1)))
                    meta[str(row["ordinal"])] = row["ordinal"]
            outcome = pm.probe(lines)
            ok = {meta[t] for t in meta if outcome.get(t, (False,))[0]}
            if ok:
                sets[name] = ok
                print(f"{name:22s} {len(ok):3d}  {sorted(ok)}")
    union = set().union(*sets.values()) if sets else set()
    print("union:", len(union), sorted(union))
    all_ord = sorted(r["ordinal"] for r in records)
    print("never accepted:", sorted(set(all_ord) - union))

if __name__ == "__main__":
    main(int(sys.argv[1]))
