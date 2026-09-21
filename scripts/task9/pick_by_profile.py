#!/usr/bin/env python3
"""Pick the convention whose whole trivial-content profile matches the table.

`pick_by_trivial.py` only asked that the record's own condensing irrep contain the
trivial representation, which still leaves conventions that flip the parity of
*other* irreps of the star (ordinal 2200: the engine puts the content on
`Y1-..Y4-` where the table has `Y1+..Y4+`).  This driver compares the engine's
full profile -- every parent irrep whose subduction contains the trivial
representation -- with the profile the pinned table records, and stops at the
first candidate that reproduces it exactly.
"""
import argparse
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/liuyichen/TB_rs/cryspglib/scripts')
REPO = '/home/liuyichen/TB_rs/cryspglib'

import probe_census as probe_mod          # noqa: E402
import generate_subduction_settings as generator  # noqa: E402
import audit_subduction_settings as census_mod  # noqa: E402

PROBE = f'{REPO}/target/release/examples/probe_subduction_settings'


def stored_profile(ordinal):
    """`ml=frequency` for the record's stored identity subduction, or None."""
    subgroup = SUBGROUPS.get(ordinal)
    if subgroup is None:
        return None
    parts = {}
    for entry in subgroup.identity_subduction():
        parts[entry.parent_ml] = entry.frequency
    return ",".join(sorted(f"{ml}={freq}" for ml, freq in parts.items()))


SUBGROUPS = {}


def load_subgroups():
    import cryspglib_stub  # noqa: F401  (placeholder, never imported)
    return {}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ordinals", required=True)
    parser.add_argument("--census", required=True)
    parser.add_argument("--profiles", required=True,
                        help="JSON {ordinal: 'ml=freq,...'} from the audit")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    wanted = json.load(open(args.ordinals, encoding="utf-8"))
    profiles = json.load(open(args.profiles, encoding="utf-8"))
    sources = {}
    for line in open(args.census, encoding="utf-8"):
        blob = json.loads(line)
        if blob.get("type") == "record":
            sources[blob["ordinal"]] = blob
    try:
        for line in open("empty_fixed.jsonl", encoding="utf-8"):
            blob = json.loads(line)
            if not sources.get(blob["ordinal"], {}).get("u"):
                sources[blob["ordinal"]] = blob
    except FileNotFoundError:
        pass

    import subprocess
    import itertools

    def shifts():
        grid = [(0, 0, 0, 1)]
        for denominator, bound in ((4, 4), (8, 6), (12, 6), (24, 12)):
            for x in range(-bound, bound + 1):
                for y in range(-bound, bound + 1):
                    for z in range(-bound, bound + 1):
                        grid.append((x, y, z, denominator))
        seen, unique = set(), []
        for shift in grid:
            if shift in seen:
                continue
            seen.add(shift)
            unique.append(shift)
        return unique

    GRID = shifts()
    print(f"candidate shifts: {len(GRID)}")

    solved = {}
    for ordinal in wanted:
        source = sources.get(ordinal)
        if not source or not source.get("u"):
            print(f"  ordinal {ordinal}: no U available")
            continue
        target = profiles.get(str(ordinal))
        u, denominator = probe_mod.rationalise(source["u"])
        lines = [probe_mod.line(str(index), ordinal, u, denominator, shift)
                 for index, shift in enumerate(GRID)]
        result = subprocess.run([PROBE, "--profile"], input="\n".join(lines) + "\n",
                                capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ordinal {ordinal}: probe failed: {result.stderr[:120]}")
            continue
        match = None
        for line in result.stdout.splitlines():
            parts = line.split()
            if len(parts) < 2 or parts[1] != "OK":
                continue
            token = next((part for part in parts if part.startswith("profile=")), None)
            if token is None:
                continue
            # Compare as sets of `label=mult` terms: the engine emits them in
            # the parent's irrep order, the stored profile is sorted by label.
            profile = ",".join(sorted(token.split("profile=")[-1].split(",")))
            if profile == target:
                match = GRID[int(parts[0])]
                break
        if match is None:
            print(f"  ordinal {ordinal}: no candidate reproduces {target}")
        else:
            solved[ordinal] = list(match)
            print(f"  ordinal {ordinal}: {match} reproduces {target}")
    json.dump({str(k): v for k, v in solved.items()}, open(args.out, "w"), indent=1)
    print(f"solved {len(solved)}/{len(wanted)}")


if __name__ == "__main__":
    main()
