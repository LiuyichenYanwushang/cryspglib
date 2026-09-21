#!/usr/bin/env python3
"""Probe every census-derived (U, delta) convention with the engine's validator."""
import argparse
import collections
import json
import subprocess
from fractions import Fraction
from math import gcd

REPO = "/home/liuyichen/TB_rs/cryspglib"
PROBE = f"{REPO}/target/release/examples/probe_subduction_settings"


def rationalise(entries):
    values = [[Fraction(v) for v in row] for row in entries]
    denominator = 1
    for row in values:
        for value in row:
            denominator = denominator * value.denominator // gcd(denominator, value.denominator)
    return [[int(value * denominator) for value in row] for row in values], denominator


def line(tag, ordinal, u, denominator, shift):
    flat = [str(v) for row in u for v in row]
    return " ".join([str(tag), str(ordinal), *flat, str(denominator),
                     *(str(v) for v in shift)])


def probe(lines, chunk=200000):
    out = {}
    for start in range(0, len(lines), chunk):
        result = subprocess.run(
            [PROBE], input="\n".join(lines[start:start + chunk]) + "\n",
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise SystemExit(f"probe failed: {result.stderr[:400]}")
        for text in result.stdout.splitlines():
            parts = text.split(None, 2)
            out[parts[0]] = (parts[1] == "OK", parts[2] if len(parts) > 2 else "")
    return out


def load(census):
    records = [
        json.loads(text)
        for text in open(census, encoding="utf-8")
        if json.loads(text).get("type") == "record"
    ]
    usable = []
    for record in records:
        if record.get("u") is None:
            continue
        u, denominator = rationalise(record["u"])
        usable.append((record, u, denominator))
    return records, usable


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("census")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    records, usable = load(args.census)
    print(f"records: {len(records)}  with an oracle U: {len(usable)}")
    lines = [line(f"{r['ordinal']}", r["ordinal"], u, d, (0, 0, 0, 1))
             for r, u, d in usable]
    outcome = probe(lines)
    ok = [r for r, _, _ in usable if outcome[str(r["ordinal"])][0]]
    bad = [(r, u, d) for r, u, d in usable if not outcome[str(r["ordinal"])][0]]
    print(f"  accepted with delta = 0: {len(ok)}")
    print(f"  rejected:                {len(bad)}")
    reasons = collections.Counter(
        outcome[str(r["ordinal"])][1].split(" for subgroup")[0] for r, _, _ in bad
    )
    for reason, count in reasons.most_common(6):
        print(f"    {count:5d}  {reason}")
    with open(args.out, "w", encoding="utf-8") as handle:
        for record, u, denominator in bad:
            handle.write(json.dumps({
                "ordinal": record["ordinal"], "parent": record["parent"],
                "child": record["child"], "irrep": record["irrep"],
                "label": record["direction_label"],
                "u": [[str(v) for v in row] for row in record["u"]],
                "denominator": denominator,
            }, sort_keys=True) + "\n")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
