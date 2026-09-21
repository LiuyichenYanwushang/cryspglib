#!/usr/bin/env python3
"""Search a validating convention for the records task 9 still cannot freeze.

These are the few records whose stored basis/origin does not reproduce the
official cell in the parent frame the engine maps into; no oracle-derived shift
exists for them.  The search space is the signed permutations plus the
census-derived rational U, crossed with the ITA origin grid, and only conventions
the engine accepts are reported.
"""
import itertools, json, sys
sys.path.insert(0, '.')
import probe_census as pm

def permutations():
    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            matrix = [[0] * 3 for _ in range(3)]
            for row in range(3):
                matrix[row][perm[row]] = signs[row]
            out.append(matrix)
    return out

def shift_grid(bound, denominator):
    return [(x, y, z, denominator)
            for x in range(-bound, bound + 1)
            for y in range(-bound, bound + 1)
            for z in range(-bound, bound + 1)]

def main(rejected, census_path, out):
    rejected = json.load(open(rejected))
    census = {}
    for line in open(census_path):
        blob = json.loads(line)
        if blob.get("type") == "record":
            census[blob["ordinal"]] = blob
    extra = {}
    try:
        for line in open("empty_fixed.jsonl", encoding="utf-8"):
            blob = json.loads(line)
            extra[blob["ordinal"]] = blob
    except FileNotFoundError:
        pass
    perms = permutations()
    grids = [ (0, 0, 0, 1) ] + shift_grid(6, 12) + shift_grid(4, 8)
    solved = {}
    for ordinal in rejected:
        source = census.get(ordinal)
        if source is None or source.get("u") is None:
            source = extra.get(ordinal)
        if source is None or source.get("u") is None:
            continue
        census_u, census_d = pm.rationalise(source["u"])
        candidates = [(census_u, census_d)] + [(m, 1) for m in perms]
        lines, meta = [], {}
        for index, (u, d) in enumerate(candidates):
            for slot, shift in enumerate(grids):
                tag = f"{index}_{slot}"
                lines.append(pm.line(tag, ordinal, u, d, shift))
                meta[tag] = (u, d, shift)
        outcome = pm.probe(lines)
        hit = next((meta[t] for t in meta if outcome.get(t, (False,))[0]), None)
        if hit:
            solved[ordinal] = {"u": hit[0], "denominator": hit[1], "shift": list(hit[2]),
                               "source": "census" if hit[0] is census_u else "permutation"}
            print(f"  ordinal {ordinal}: solved ({solved[ordinal]['source']})", flush=True)
        else:
            print(f"  ordinal {ordinal}: no convention in the search space", flush=True)
    json.dump({str(k): v for k, v in sorted(solved.items())}, open(out, "w"), indent=1)
    print(f"solved {len(solved)}/{len(rejected)}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
