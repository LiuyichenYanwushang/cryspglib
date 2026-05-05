#!/usr/bin/env python3
"""Parse double-valued (spinor) irrep data from irrepTables package.

The irrepTables package (pip install irreptables) contains data files
sourced from the Bilbao Crystallographic Server. Each file covers one
space group:

  ~/.local/lib/python*/site-packages/irreptables/tables/irreps-SG=*-spin.dat

Output format:
- SG#, k-point label (GM, X, ...), coords, irrep label, dim, character table

## Bilbao SU(2) storage convention (verified by scripts/test_su2_closure.py)

Each symmetry operation line in spin.dat stores 20 values:

    rot[9] trans[3] amp[4] phase[4]

where:
  amp[4]   = |U_ij| — amplitudes of the 4 complex 2×2 matrix elements
  phase[4] = arg(U_ij)/π — phases in units of π

The complex SU(2) matrix is reconstructed as:

    U_ij = amp[ij] · exp(iπ · phase[ij])

Then decomposed into real Pauli coefficients (u₀, u₁, u₂, u₃):

    U = u₀·I + i(u₁·σx + u₂·σy + u₃·σz)

      = [[u₀ + iu₃,    u₂ + iu₁],
         [-u₂ + iu₁,    u₀ - iu₃]]

For crystallographic point-group rotations, the Pauli coefficients
take values only from the set {0, ±½, ±1/√2, ±√3/2, ±1}.
These are stored as f64 rounded to the nearest exact algebraic value.
"""
import math, os, re, sys, glob
from collections import defaultdict

# ── Exact Pauli coefficient rounding ────────────────────────────────────
# For crystallographic double-group operations, the Pauli coefficients
# u₀,u₁,u₂,u₃ ∈ {0, ±½, ±1/√2, ±√3/2, ±1}.
# We round floating-point values to the nearest exact algebraic number
# to eliminate numerical noise from sin/cos computations.

_SQRT2 = math.sqrt(2.0)
_SQRT3 = math.sqrt(3.0)

_EXACT_PAULI_VALUES = [
    0.0,
    0.5, -0.5,
    1.0 / _SQRT2, -1.0 / _SQRT2,
    _SQRT3 / 2.0, -_SQRT3 / 2.0,
    1.0, -1.0,
]


def _round_to_exact_pauli(val, tol=1e-10):
    """Round a floating-point Pauli coefficient to the nearest exact value."""
    for exact in _EXACT_PAULI_VALUES:
        if abs(val - exact) < tol:
            return exact
    return val


def _round_amplitude(val):
    """Round amplitude to exact value {0, 1/√2, 1}.

    Uses a relaxed tolerance because spin.dat files store amplitudes with
    only ~5 significant digits (e.g. 0.70711 for 1/√2 ≈ 0.70710678).
    """
    for exact in [0.0, 1.0 / _SQRT2, 1.0]:
        if abs(val - exact) < 1e-4:
            return exact
    return val


def _amp_phase_to_pauli(amp, phase):
    """Convert Bilbao polar encoding (amp[4] + phase[4]/π) to Pauli coefficients.

    Returns [u₀, u₁, u₂, u₃] as exact f64 values.

    The amp values in spin.dat files are rounded to ~5 significant digits
    (e.g. 0.70711 for 1/√2).  We round them to the exact algebraic value
    {0, 1/√2, 1} before multiplying, so that uᵢ = amp × cos/sin(π·phase)
    comes out exact (e.g. 0.5 instead of 0.500002276).
    """
    # Round amplitudes to exact {0, 1/√2, 1} before multiplying.
    # This eliminates the ~5-digit precision loss in spin.dat amplitude values.
    exact_amp = [_round_amplitude(a) for a in amp]

    u0 = _round_to_exact_pauli(exact_amp[0] * math.cos(math.pi * phase[0]))
    u3 = _round_to_exact_pauli(exact_amp[0] * math.sin(math.pi * phase[0]))
    u2 = _round_to_exact_pauli(exact_amp[1] * math.cos(math.pi * phase[1]))
    u1 = _round_to_exact_pauli(exact_amp[1] * math.sin(math.pi * phase[1]))

    return [u0, u1, u2, u3]


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
        spin_ops: list of dicts with keys: rot[9], trans[3], su2[4]
                  su2[4] = Pauli coefficients (u₀, u₁, u₂, u₃), exact f64
        irreps: list of dicts with keys:
            k_label, kx/ky/kz/kd, ml_label, dim, characters, op_indices
    """
    with open(filepath) as f:
        lines = f.readlines()

    sg = None
    spin_ops = []  # global symmetry operations with SU(2) lifts
    irreps = []

    # Parse header
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("SG="):
            sg = int(line.split("=")[1])
        elif line.startswith("nsym="):
            pass
        elif line.startswith("spinor="):
            pass
        elif line.startswith("symmetries="):
            i += 1
            break
        i += 1

    # Parse symmetry operations (nsym lines)
    # Format: R(3x3) 9ints | t(3) 3floats | amp(4) 4floats | phase(4)/π 4floats
    # Converted to Pauli coefficients: SU(2) = u₀I + i(u₁σx + u₂σy + u₃σz)
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("kpoint"):
            break
        parts = line.split()
        if len(parts) >= 20:
            rot = [int(x) for x in parts[0:9]]
            trans = [float(x) for x in parts[9:12]]
            amp = [float(x) for x in parts[12:16]]
            phase = [float(x) for x in parts[16:20]]
            su2 = _amp_phase_to_pauli(amp, phase)
            spin_ops.append({"rot": rot, "trans": trans, "su2": su2})
        elif len(parts) >= 16:
            # Fallback for files without the extra 4 columns
            rot = [int(x) for x in parts[0:9]]
            trans = [float(x) for x in parts[9:12]]
            amp = [float(x) for x in parts[12:16]]
            phase = [0.0, 0.0, 0.0, 0.0]
            su2 = _amp_phase_to_pauli(amp, phase)
            spin_ops.append({"rot": rot, "trans": trans, "su2": su2})
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

            # Compute k-vector denominator from coords.
            # Use the SMALLEST common denominator so that spinor and scalar
            # irreps at the same k-point get the same (kx,ky,kz,kd) tuple
            # and are correctly grouped by kpoints_of().
            if current_kvec:
                kx, ky, kz = current_kvec
                from math import gcd
                kd = 1
                for d in [1, 2, 3, 4, 6]:
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
                # Reduce to simplest form
                g = gcd(abs(kx_i), abs(ky_i))
                g = gcd(g, abs(kz_i))
                g = gcd(g, kd)
                if g > 1:
                    kx_i //= g
                    ky_i //= g
                    kz_i //= g
                    kd //= g
            else:
                kx_i = ky_i = kz_i = 0
                kd = 1

            irreps.append({
                "sg": sg,
                "k_label": current_k,
                "kx": kx_i, "ky": ky_i, "kz": kz_i, "kd": kd,
                "ml_label": label,
                "dim": dim,
                "characters": chars,  # full: standard + extra (if any)
                "op_indices": current_op_indices,
            })

        i += 1

    return sg, spin_ops, irreps


def parse_all_spinor():
    """Parse all 230 spin.dat files.

    Returns:
        all_irreps: list of spinor irrep dicts
        all_spin_ops: dict SG# -> list of spin op dicts (rot, trans, su2)
                      su2[4] = Pauli coefficients (u₀, u₁, u₂, u₃)
    """
    tables_dir = find_tables_dir()
    files = sorted(glob.glob(os.path.join(tables_dir, "irreps-SG=*-spin.dat")))

    all_irreps = []
    all_spin_ops = {}  # SG# -> list of spin ops
    for f in files:
        sg, spin_ops, irreps = parse_spinor_file(f)
        all_spin_ops[sg] = spin_ops
        all_irreps.extend(irreps)

    # Sort by SG then by k-label for contiguity
    all_irreps.sort(key=lambda x: (x["sg"], x["k_label"], x["ml_label"]))

    return all_irreps, all_spin_ops


if __name__ == "__main__":
    irreps, spin_ops = parse_all_spinor()
    print(f"Parsed {len(irreps)} spinor irreps from {len(spin_ops)} SGs")
    # Show sample ops
    for sg in sorted(spin_ops.keys())[:3]:
        ops = spin_ops[sg]
        print(f"  SG{sg}: {len(ops)} spin ops")
        if ops:
            print(f"    op[0]: rot={ops[0]['rot'][:3]}... trans={ops[0]['trans']} su2={ops[0]['su2']}")
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
