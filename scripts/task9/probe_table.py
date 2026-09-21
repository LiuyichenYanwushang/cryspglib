#!/usr/bin/env python3
"""Read the committed frozen table and ask the engine which entries it accepts."""
import json, re, sys
sys.path.insert(0, '.')
import probe_census as pm

CONST = {
    "IDENTITY": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    "CYCLIC": [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
}

def split_top_level(text):
    """Split a Rust tuple body on commas that are not inside brackets."""
    parts, depth, current = [], 0, []
    for char in text:
        if char in "[(":
            depth += 1
        elif char in "])":
            depth -= 1
        if char == "," and depth == 0:
            parts.append("".join(current).strip())
            current = []
        else:
            current.append(char)
    if current:
        parts.append("".join(current).strip())
    return parts


def parse(path):
    entries = []
    for line in open(path, encoding="utf-8"):
        stripped = line.strip()
        if not stripped.startswith("(") or not stripped.endswith("),"):
            continue
        fields = split_top_level(stripped[1:-2])
        if len(fields) != 6 or not fields[0].isdigit():
            continue
        ordinal, parent, child = int(fields[0]), int(fields[1]), int(fields[2])
        token, denominator, shift_token = fields[3], int(fields[4]), fields[5]
        if token in CONST:
            u = CONST[token]
        else:
            nums = [int(v) for v in re.findall(r"-?\d+", token)]
            assert len(nums) == 9, line
            u = [nums[0:3], nums[3:6], nums[6:9]]
        shift = [0, 0, 0, 1] if shift_token == "NO_SHIFT" else [
            int(v) for v in re.findall(r"-?\d+", shift_token)]
        entries.append((ordinal, parent, child, u, denominator, shift))
    return entries

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
