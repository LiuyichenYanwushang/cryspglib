#!/usr/bin/env python3
"""Diagnose the correct antiunitary SU(2) composition formula.

For spinor irreps with Bilbao extra chars, the extra values are the
pre-computed Wigner contributions χ((a₀h)²).  We test multiple
candidate formulas for computing (a₀h)² from spinor SU(2) Pauli
coefficients, and compare the resulting Wigner classification with
the extra-char reference.

Formulas tested:
  F0  = (A·H)²           — current Rust, pure unitary
  -F0 = Θ²·(A·H)²        — with Θ²=-1 central flip
  F1  = A·H*·A*·H        — antiunitary: conjugate H
  F2  = J·(A·H)*·J·(A·H) — explicit J=iσy time reversal
  F3  = (H·A)²            — reversed order
  -F3 = Θ²·(H·A)²
  F4  = H·A*·H*·A         — reversed + antiunitary
  F5  = J·(H·A)*·J·(H·A) — reversed + J

Usage:
  python3 scripts/test_spinor_wigner_formula.py [--sg SG] [--max-cases N]
"""
import math, sys, os, argparse

# ── Pauli coefficient helpers ──────────────────────────────────────────

def su2_compose(a, b):
    """Quaternion multiply: (u₀,u)·(v₀,v) = (u₀v₀-u·v, u₀v+v₀u+u×v)"""
    u0, u1, u2, u3 = a
    v0, v1, v2, v3 = b
    return (
        u0*v0 - u1*v1 - u2*v2 - u3*v3,
        u0*v1 + u1*v0 + u2*v3 - u3*v2,
        u0*v2 - u1*v3 + u2*v0 + u3*v1,
        u0*v3 + u1*v2 - u2*v1 + u3*v0,
    )

def su2_conj(u):
    """Complex conjugate: U* = u₀I + i((-u₁)σx + u₂σy + (-u₃)σz)"""
    return (u[0], -u[1], u[2], -u[3])

def su2_neg(u):
    """Central element Ē = -I."""
    return (-u[0], -u[1], -u[2], -u[3])

# J = iσy in Pauli coefficients
J = (0.0, 0.0, 1.0, 0.0)

def su2_same_up_to_sign(a, b):
    """Check if a ≈ ±b. Returns (matched: bool, sign_flip: bool or None)."""
    na = math.sqrt(sum(x*x for x in a))
    nb = math.sqrt(sum(x*x for x in b))
    if na < 1e-12 or nb < 1e-12:
        return (False, None)
    dot = sum(a[i]*b[i] for i in range(4)) / (na * nb)
    if abs(dot - 1.0) < 1e-6:
        return (True, False)   # a ≈ +b
    if abs(dot + 1.0) < 1e-6:
        return (True, True)    # a ≈ -b (central Ē)
    return (False, None)


# ── Formula definitions ────────────────────────────────────────────────

def make_formulas():
    """Return list of (name, fn, uses_explicit_Theta2) tuples.

    fn(A, H) -> U_sq (the spinor lift of (a₀h)² in Pauli coefficients).
    uses_explicit_Theta2=True means the formula ALREADY includes Θ²=-1.
    """
    formulas = []

    # F0: (A·H)²  — current Rust
    formulas.append(("F0_(AH)^2", lambda A, H: su2_compose(su2_compose(A, H), su2_compose(A, H)), False))

    # -F0: Θ² · (A·H)²
    formulas.append(("-F0_Theta2*(AH)^2", lambda A, H: su2_neg(su2_compose(su2_compose(A, H), su2_compose(A, H))), True))

    # F1: A·H*·A*·H
    formulas.append(("F1_AH*A*H", lambda A, H: su2_compose(su2_compose(su2_compose(A, su2_conj(H)), su2_conj(A)), H), True))

    # F2: J·(A·H)*·J·(A·H)
    formulas.append(("F2_J(AH)*J(AH)", lambda A, H: su2_compose(su2_compose(su2_compose(J, su2_conj(su2_compose(A, H))), J), su2_compose(A, H)), True))

    # F3: (H·A)²
    formulas.append(("F3_(HA)^2", lambda A, H: su2_compose(su2_compose(H, A), su2_compose(H, A)), False))

    # -F3
    formulas.append(("-F3_Theta2*(HA)^2", lambda A, H: su2_neg(su2_compose(su2_compose(H, A), su2_compose(H, A))), True))

    # F4: H·A*·H*·A
    formulas.append(("F4_HA*H*A", lambda A, H: su2_compose(su2_compose(su2_compose(H, su2_conj(A)), su2_conj(H)), A), True))

    # F5: J·(H·A)*·J·(H·A)
    formulas.append(("F5_J(HA)*J(HA)", lambda A, H: su2_compose(su2_compose(su2_compose(J, su2_conj(su2_compose(H, A))), J), su2_compose(H, A)), True))

    return formulas


# ── Load spinor data ───────────────────────────────────────────────────

def load_spinor_data():
    """Load from parse_spinor_data.py, preserving all 8 SU(2) fields."""
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
    from parse_spinor_data import parse_all_spinor
    irreps, spin_ops = parse_all_spinor()
    return irreps, spin_ops


# ── Get magnetic group ops from Rust (via a small dump) ──
# We'll use the Rust test infrastructure to dump the operations we need.
# For now, we hardcode a few well-known cases.

# SG 118, UNI 1066 (BNS 128.406): unitary subgroup of the magnetic group
# The antiunitary coset representative a₀ = T·g₀

# We'll get the magnetic ops by running the Rust diagnostic test.
# But for a quick start, let's use the specific test cases we already know.

def load_test_cases():
    """Return list of (name, h_sg, uni, a0_seitz, h_seitz_list, k) for known test cases.

    These need to be filled in from the Rust diagnostics or from known data.
    """
    # We'll populate this from the Rust BCS tests
    return []


# ── Main diagnostic ────────────────────────────────────────────────────

def diagnose_one_irrep(ir, spin_ops_by_sg, mag_ops_data):
    """Test all formulas on a single spinor irrep with extra chars.

    Returns dict of formula_name → classification result.
    """
    extra = ir.get('spin_extra', [])  # need to add this to parse_spinor_data
    # For now, use characters beyond n_lg as extra
    # Actually, extra chars are stored separately in the Rust data.
    # In the parsed Python data, ir['characters'] includes BOTH lg and extra.
    pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sg', type=int, default=0, help='Space group to test (0=all)')
    parser.add_argument('--max-cases', type=int, default=50, help='Max test cases')
    args = parser.parse_args()

    irreps, spin_ops = load_spinor_data()

    # Count spinor irreps with extra chars
    with_extra = []
    for ir in irreps:
        ops = ir.get('op_indices', [])
        chars = ir.get('characters', [])
        if len(chars) > len(ops):  # has extra chars
            with_extra.append(ir)

    print(f"Total spinor irreps: {len(irreps)}")
    print(f"With extra chars: {len(with_extra)}")
    print()

    # Show formula definitions
    formulas = make_formulas()
    print(f"Testing {len(formulas)} formulas:")
    for name, _, explicit in formulas:
        print(f"  {name:>25s}  {'explicit Θ²' if explicit else 'needs manual Θ²'}")
    print()

    # For now, we need magnetic group data. Print what's needed.
    print("To complete this test, we need magnetic group operations for")
    print("specific UNI numbers.  Running Rust diagnostic to extract key cases...")
    print()
    print("Key test cases from Rust BCS validation:")
    print("  128.406 Z-point (UNI 1066, SG 118): Z6, Z7 spinor irreps")
    print("  165.95  L-point (UNI 1272, SG 165): spinor irreps at L")
    print()

    # Show spinor SU(2) data for Z6 (SG 118) as an example
    ops_118 = spin_ops.get(118, [])
    print(f"SG 118: {len(ops_118)} spin ops")
    for ir in irreps:
        if ir['sg'] == 118 and ir['ml_label'].startswith('Z') and ir.get('spinor'):
            extra_start = len(ir.get('op_indices', []))
            chars = ir.get('characters', [])
            lg_chars = chars[:extra_start]
            extra_chars = chars[extra_start:]
            print(f"  {ir['ml_label']}: dim={ir['dim']} k=({ir['kx']}/{ir['kd']},{ir['ky']}/{ir['kd']},{ir['kz']}/{ir['kd']})")
            print(f"    op_indices={ir.get('op_indices', [])[:8]}...")
            print(f"    lg_chars={[round(c,2) for c in lg_chars[:8]]}...")
            print(f"    extra_chars={[round(c,2) for c in extra_chars[:8]]}")
            print(f"    extra_sum={sum(extra_chars):.4f}")

    print()
    print("=" * 60)
    print("NEXT STEP: need to run Rust code to dump magnetic ops for UNI 1066")
    print("then test formulas on Z6/Z7.  This script provides the formula")
    print("framework; the magnetic group data needs to come from spglib.")
    print("=" * 60)


if __name__ == '__main__':
    main()
