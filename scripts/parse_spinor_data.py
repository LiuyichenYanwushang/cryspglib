#!/usr/bin/env python3
"""Parse double-valued (spinor) irrep data from irrepTables package.

The irrepTables package (pip install irreptables) contains data files
sourced from the Bilbao Crystallographic Server. Each file covers one
space group:

  ~/.local/lib/python*/site-packages/irreptables/tables/irreps-SG=*-spin.dat

Output format:
- SG#, k-point label (GM, X, ...), coords, irrep label, dim, character table
"""
import os, re, sys, glob
from collections import defaultdict


def find_tables_dir():
    """Locate the irrepTables data directory."""
    candidates = []
    for p in sys.path:
        d = os.path.join(p, "irreptables", "tables")
        if os.path.isdir(d):
            candidates.append(d)
    # Also try user site-packages
    import site
    for sp in site.getusersitepackages(), site.getsitepackages()[0] if site.getsitepackages() else []:
        d = os.path.join(sp, "irreptables", "tables")
        if os.path.isdir(d) and d not in candidates:
            candidates.append(d)
    if candidates:
        return candidates[0]
    raise FileNotFoundError("irreptables/tables directory not found. pip install irreptables")


def _round_char(x, eps=1e-8):
    """Round character value to clean float."""
    if abs(x) < eps:
        return 0.0
    r = round(x)
    if abs(x - r) < eps:
        return float(r)
    for n in range(-12, 13):
        for d in (2, 3, 4, 6, 8):
            if abs(x - n / d) < eps:
                return n / d
    return round(x, 10)


def parse_spinor_file(filepath):
    """Parse one spin.dat file.

    Returns:
        sg: int, space group number
        irreps: list of dicts with keys:
            k_label, kx/ky/kz/kd, ml_label, dim, characters
    """
    with open(filepath) as f:
        lines = f.readlines()

    sg = None
    irreps = []

    # Parse header
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("SG="):
            sg = int(line.split("=")[1])
        elif line.startswith("name="):
            pass
        elif line.startswith("nsym="):
            pass
        elif line.startswith("spinor="):
            pass
        elif line.startswith("symmetries="):
            i += 1
            break
        i += 1

    # Skip symmetry operations (nsym lines)
    # Each line: R(3x3) + t(3) + SU2(4) + ... (variable format)
    # Count until we hit a blank line or kpoint
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("kpoint"):
            break
        i += 1

    # Parse kpoints and irreps
    current_k = None
    current_kvec = None
    current_op_indices = None

    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("kpoint"):
            # kpoint  GM : 0.0 0.0 0.0  : 1 2 3 4 ...
            parts = line.split(":")
            k_name = parts[0].replace("kpoint", "").strip()
            coords = [float(x) for x in parts[1].strip().split()]
            op_indices = [int(x) - 1 for x in parts[2].strip().split()]  # 0-indexed
            current_k = k_name
            current_kvec = coords
            current_op_indices = op_indices
        elif line.startswith("-"):
            # -GM6 2    2.0  0.0  ...
            parts = line.split()
            if len(parts) < 3:
                i += 1
                continue
            label = parts[0][1:]  # strip leading '-'
            dim = int(parts[1])
            chars_raw = [float(x) for x in parts[2:]]
            # Characters are in order of the little group operations
            # We need to store them in that order
            chars = [_round_char(c) for c in chars_raw]

            # Compute k-vector denominator from coords
            if current_kvec:
                kx, ky, kz = current_kvec
                # Find common denominator (up to 6)
                kd = 1
                for d in [6, 4, 3, 2, 1]:
                    ok = True
                    for v in [kx, ky, kz]:
                        v_scaled = v * d
                        if abs(v_scaled - round(v_scaled)) > 1e-6:
                            ok = False
                            break
                    if ok:
                        kd = d
                        break
                kx_i = int(round(kx * kd))
                ky_i = int(round(ky * kd))
                kz_i = int(round(kz * kd))
            else:
                kx_i = ky_i = kz_i = 0
                kd = 1

            irreps.append({
                "sg": sg,
                "k_label": current_k,
                "kx": kx_i, "ky": ky_i, "kz": kz_i, "kd": kd,
                "ml_label": label,
                "dim": dim,
                "characters": chars,
            })

        i += 1

    return sg, irreps


def parse_all_spinor():
    """Parse all 230 spin.dat files."""
    tables_dir = find_tables_dir()
    files = sorted(glob.glob(os.path.join(tables_dir, "irreps-SG=*-spin.dat")))

    all_irreps = []
    for f in files:
        sg, irreps = parse_spinor_file(f)
        all_irreps.extend(irreps)

    # Sort by SG then by k-label for contiguity
    all_irreps.sort(key=lambda x: (x["sg"], x["k_label"], x["ml_label"]))

    return all_irreps


if __name__ == "__main__":
    irreps = parse_all_spinor()
    print(f"Parsed {len(irreps)} spinor irreps from 230 SGs")
    by_sg = defaultdict(int)
    for ir in irreps:
        by_sg[ir["sg"]] += 1
    print(f"  SG range: {min(by_sg)}-{max(by_sg)}")
    print(f"  Max per SG: {max(by_sg.values())} (SG {max(by_sg, key=by_sg.get)})")
    # Show sample
    for ir in irreps[:5]:
        print(f"  SG{ir['sg']} {ir['ml_label']} dim={ir['dim']} "
              f"k=({ir['kx']}/{ir['kd']},{ir['ky']}/{ir['kd']},{ir['kz']}/{ir['kd']}) "
              f"chars={ir['characters'][:5]}...")
