#!/usr/bin/env python3
"""Validate Bilbao SU(2) convention by testing group closure.

For every pair of spin operations (i,j) in a spin.dat file,
compute the spatial product {R_i|t_i}{R_j|t_j} and the SU(2)
product U_i·U_j.  Verify that U_i·U_j equals ±U_k for the
spin operation k with the matching spatial product.

Tests both real 4-number and complex 8-number interpretations
of the SU(2) data in columns 12-19.
"""
import math, sys, glob, os

EPS = 1e-7

TABLES_DIR = None
for d in [
    "/home/liuyichen/.local/lib/python3.7/site-packages/irreptables/tables",
]:
    if os.path.isdir(d):
        TABLES_DIR = d
        break


def find_tables_dir():
    candidates = []
    for p in sys.path:
        d = os.path.join(p, "irreptables", "tables")
        if os.path.isdir(d):
            candidates.append(d)
    import site
    for sp in [site.getusersitepackages()] + site.getsitepackages():
        d = os.path.join(sp, "irreptables", "tables")
        if os.path.isdir(d) and d not in candidates:
            candidates.append(d)
    if TABLES_DIR and os.path.isdir(TABLES_DIR) and TABLES_DIR not in candidates:
        candidates.append(TABLES_DIR)
    return candidates[0] if candidates else None


def parse_spin_dat(path):
    """Parse spin.dat, returning ops with all 8 SU(2) numbers (4+4)."""
    ops = []
    with open(path) as f:
        lines = f.readlines()
    in_sym = False
    for line in lines:
        s = line.strip()
        if not s:
            continue
        if s.startswith("symmetries="):
            in_sym = True
            continue
        if in_sym and s.startswith("kpoint"):
            break
        if in_sym:
            parts = s.split()
            if len(parts) >= 20:
                rot = tuple(int(x) for x in parts[0:9])
                trans = tuple(float(x) for x in parts[9:12])
                su2_real = tuple(float(x) for x in parts[12:16])
                su2_imag = tuple(float(x) for x in parts[16:20])
                ops.append((rot, trans, su2_real, su2_imag))
            elif len(parts) >= 16:
                rot = tuple(int(x) for x in parts[0:9])
                trans = tuple(float(x) for x in parts[9:12])
                su2_real = tuple(float(x) for x in parts[12:16])
                su2_imag = (0.0, 0.0, 0.0, 0.0)
                ops.append((rot, trans, su2_real, su2_imag))
    return ops


def matmul3(a, b):
    out = []
    for i in range(3):
        for j in range(3):
            out.append(sum(a[3*i+k] * b[3*k+j] for k in range(3)))
    return tuple(out)


def matvec3(r, v):
    return tuple(
        r[3*i+0] * v[0] + r[3*i+1] * v[1] + r[3*i+2] * v[2]
        for i in range(3)
    )


def mod1(x):
    y = x - math.floor(x)
    if abs(y - 1.0) < EPS:
        y = 0.0
    if abs(y) < EPS:
        y = 0.0
    return y


def same_trans(a, b):
    return all(
        abs(mod1(a[i] - b[i])) < EPS or abs(mod1(a[i] - b[i]) - 1.0) < EPS
        for i in range(3)
    )


def compose_seitz(op1, op2):
    r1, t1, _, _ = op1
    r2, t2, _, _ = op2
    r = matmul3(r1, r2)
    rt2 = matvec3(r1, t2)
    t = tuple(mod1(t1[i] + rt2[i]) for i in range(3))
    return r, t


def find_spatial_product(ops, r, t):
    return [idx for idx, (rr, tt, _, _) in enumerate(ops)
            if rr == r and same_trans(tt, t)]


# ── Rotation matrix utilities ──

def trace_rot(r):
    return r[0] + r[4] + r[8]


def det3(r):
    return (
        r[0]*(r[4]*r[8] - r[5]*r[7])
        - r[1]*(r[3]*r[8] - r[5]*r[6])
        + r[2]*(r[3]*r[7] - r[4]*r[6])
    )


def mat_identity():
    return (1,0,0, 0,1,0, 0,0,1)


def rotation_order(r, max_order=12):
    cur = mat_identity()
    for n in range(1, max_order + 1):
        cur = matmul3(cur, r)
        if cur == mat_identity():
            return n
    return None


# ── Real 4-number conventions ──

def su2_mul_flat2x2(a, b):
    return (
        a[0]*b[0] + a[1]*b[2],
        a[0]*b[1] + a[1]*b[3],
        a[2]*b[0] + a[3]*b[2],
        a[2]*b[1] + a[3]*b[3],
    )


def su2_mul_quat_wxyz(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return (
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
    )


def su2_mul_pauli(a, b):
    ar, ai, br, bi = a
    cr, ci, dr, di = b
    return (
        ar*cr - ai*ci - br*dr - bi*di,
        ar*ci + ai*cr - br*di + bi*dr,
        ar*dr - ai*di + br*cr + bi*ci,
        ar*di + ai*dr + br*ci - bi*cr,
    )


# ── Complex 2×2 conventions (8 numbers → complex 2×2) ──

def to_complex2x2_real_first(u_real, u_imag):
    """[r00,r01,r10,r11] + [i00,i01,i10,i11] → 4 complex numbers."""
    return (
        complex(u_real[0], u_imag[0]),
        complex(u_real[1], u_imag[1]),
        complex(u_real[2], u_imag[2]),
        complex(u_real[3], u_imag[3]),
    )


def to_complex2x2_interleaved(u_real, u_imag):
    """[r00,i00,r01,i01] + [r10,i10,r11,i11] → 4 complex numbers."""
    return (
        complex(u_real[0], u_real[1]),
        complex(u_imag[0], u_imag[1]),
        complex(u_real[2], u_real[3]),
        complex(u_imag[2], u_imag[3]),
    )


def cmul2(a, b):
    """Complex 2×2 matrix multiply, flattened."""
    return (
        a[0]*b[0] + a[1]*b[2],
        a[0]*b[1] + a[1]*b[3],
        a[2]*b[0] + a[3]*b[2],
        a[2]*b[1] + a[3]*b[3],
    )


def cvec_norm(v):
    n = math.sqrt(sum(abs(x)**2 for x in v))
    return tuple(x / n for x in v) if n > 1e-15 else v


def complex_same_up_to_sign(a, b):
    a = cvec_norm(a)
    b = cvec_norm(b)
    dot = sum(a[i].conjugate() * b[i] for i in range(4))
    if abs(dot - 1.0) < 1e-6:
        return +1
    if abs(dot + 1.0) < 1e-6:
        return -1
    return 0


# Complex convention factory: (real_first | interleaved) × (forward | reverse)

def make_complex_mul(convert, reverse):
    def fn(su2_a, su2_b):
        A = convert(su2_a[0], su2_a[1])
        B = convert(su2_b[0], su2_b[1])
        if reverse:
            return cmul2(B, A)
        else:
            return cmul2(A, B)
    return fn


# ── Real vector comparison ──

def vec_norm(v):
    n = math.sqrt(sum(x*x for x in v))
    return tuple(x / n for x in v) if n > 1e-15 else v


def same_up_to_sign(a, b):
    a = vec_norm(a)
    b = vec_norm(b)
    dot = sum(a[i] * b[i] for i in range(4))
    if abs(dot - 1.0) < 1e-6:
        return +1
    if abs(dot + 1.0) < 1e-6:
        return -1
    return 0


# ── Convention registry ──

# Real conventions: compare_fn gets (su2_real_tuple_of_4, su2_real_tuple_of_4)
# Complex conventions: compare_fn gets ((su2_real_4, su2_imag_4), (su2_real_4, su2_imag_4))

# Real conventions: get_su2 returns 4-tuple, mul_fn(4-tuple, 4-tuple) → 4-tuple
# Complex conventions: get_su2 returns (4-tuple, 4-tuple), mul_fn(pair, pair) → complex result

REAL_CONVENTIONS = [
    ("flat2x2", su2_mul_flat2x2, same_up_to_sign),
    ("quat-wxyz", su2_mul_quat_wxyz, same_up_to_sign),
    ("pauli", su2_mul_pauli, same_up_to_sign),
]

def to_complex2x2_polar(amp, phase):
    """Bilbao spin.dat amplitude + phase/pi encoding.

    amp[4]   = |U_ij|
    phase[4] = arg(U_ij) / π
    → U_ij = amp * exp(iπ * phase)
    """
    return tuple(
        amp[i] * complex(
            math.cos(math.pi * phase[i]),
            math.sin(math.pi * phase[i]),
        )
        for i in range(4)
    )


def _make_cplx_mul_cmp(convert, reverse):
    """Return (mul_fn, cmp_fn) for a complex convention.

    mul_fn: ((real4, imag4), (real4, imag4)) → 4-complex result
    cmp_fn: (4-complex result, (real4, imag4) target) → -1/0/+1
    """
    def mul_fn(a, b):
        A = convert(a[0], a[1])
        B = convert(b[0], b[1])
        if reverse:
            return cmul2(B, A)
        else:
            return cmul2(A, B)

    def cmp_fn(result_4cplx, target_pair):
        target_4cplx = convert(target_pair[0], target_pair[1])
        return complex_same_up_to_sign(result_4cplx, target_4cplx)

    return mul_fn, cmp_fn

def su2_mul_polar(a, b):
    """Amplitude+phase/pi encoding → complex 2×2 multiply."""
    A = to_complex2x2_polar(a[0], a[1])
    B = to_complex2x2_polar(b[0], b[1])
    return cmul2(A, B)


def su2_mul_polar_rev(a, b):
    """Same but reversed multiplication order."""
    A = to_complex2x2_polar(a[0], a[1])
    B = to_complex2x2_polar(b[0], b[1])
    return cmul2(B, A)


def cmp_polar(result_4cplx, target_pair):
    target_4cplx = to_complex2x2_polar(target_pair[0], target_pair[1])
    return complex_same_up_to_sign(result_4cplx, target_4cplx)


CMPLX_CONVENTIONS = [
    ("cplx-real_first", *_make_cplx_mul_cmp(to_complex2x2_real_first, False)),
    ("cplx-real_first_rev", *_make_cplx_mul_cmp(to_complex2x2_real_first, True)),
    ("cplx-interleaved", *_make_cplx_mul_cmp(to_complex2x2_interleaved, False)),
    ("cplx-interleaved_rev", *_make_cplx_mul_cmp(to_complex2x2_interleaved, True)),
    ("polar_phase_pi", su2_mul_polar, cmp_polar),
    ("polar_phase_pi_rev", su2_mul_polar_rev, cmp_polar),
]


# ── Core test logic ──

def test_file(ops, mul_fn, cmp_fn, get_su2):
    """Generic closure test.

    get_su2(op) returns the SU(2) tuple (real-only or real+imag)
    to pass to mul_fn and cmp_fn.
    """
    if len(ops) < 2:
        return None
    spatial_missing = 0
    spatial_matched = 0
    su2_mismatch = 0
    su2_matched = 0
    sign_p = 0
    sign_m = 0
    for i, op_i in enumerate(ops):
        for j, op_j in enumerate(ops):
            r, t = compose_seitz(op_i, op_j)
            ks = find_spatial_product(ops, r, t)
            if not ks:
                spatial_missing += 1
                continue
            spatial_matched += 1
            qij = mul_fn(get_su2(op_i), get_su2(op_j))
            matched = False
            for k in ks:
                s = cmp_fn(qij, get_su2(ops[k]))
                if s != 0:
                    matched = True
                    if s > 0: sign_p += 1
                    else: sign_m += 1
                    break
            if matched:
                su2_matched += 1
            else:
                su2_mismatch += 1
    total = spatial_matched + spatial_missing
    return {
        "total": total,
        "spatial_matched": spatial_matched,
        "spatial_missing": spatial_missing,
        "su2_matched": su2_matched,
        "su2_mismatch": su2_mismatch,
        "+": sign_p, "-": sign_m,
    }


# ── SG 221 detailed diagnosis ──

def diagnose_sg221(path, mul_fn, cmp_fn, get_su2, max_print=30):
    """Print detailed mismatch info for SG 221."""
    ops = parse_spin_dat(path)
    print(f"\n  Total spin operations: {len(ops)}")

    # Show a few sample ops
    print(f"\n  Sample ops (first 3):")
    for idx in range(min(3, len(ops))):
        r, t, sr, si = ops[idx]
        print(f"    [{idx}] R={r} t={t} su2_re={sr} su2_im={si} "
              f"tr={trace_rot(r)} ord={rotation_order(r)} det={det3(r)}")

    buckets = {}
    printed = 0
    for i, op_i in enumerate(ops):
        for j, op_j in enumerate(ops):
            r, t = compose_seitz(op_i, op_j)
            ks = find_spatial_product(ops, r, t)
            if not ks:
                continue
            qij = mul_fn(get_su2(op_i), get_su2(op_j))
            matched = False
            for k in ks:
                if cmp_fn(qij, get_su2(ops[k])) != 0:
                    matched = True
                    break
            if matched:
                continue

            k = ks[0]
            ri, rj, rk = op_i[0], op_j[0], ops[k][0]
            key = (
                trace_rot(ri), rotation_order(ri), det3(ri),
                trace_rot(rj), rotation_order(rj), det3(rj),
                trace_rot(rk), rotation_order(rk), det3(rk),
            )
            buckets[key] = buckets.get(key, 0) + 1

            if printed < max_print:
                print(f"\n  mismatch #{printed}: i={i} j={j} -> k={ks}")
                print(f"    Ri: tr={trace_rot(ri)} ord={rotation_order(ri)} det={det3(ri):+d}  "
                      f"Rj: tr={trace_rot(rj)} ord={rotation_order(rj)} det={det3(rj):+d}  "
                      f"Rk: tr={trace_rot(rk)} ord={rotation_order(rk)} det={det3(rk):+d}")
                su2_i = get_su2(op_i)
                su2_j = get_su2(op_j)
                print(f"    Ui={su2_i}")
                print(f"    Uj={su2_j}")
                print(f"    UiUj={qij}")
                for kk in ks:
                    print(f"    Uk[{kk}]={get_su2(ops[kk])}")
                printed += 1

    print(f"\n  Mismatch buckets by (tri,ordi,deti, trj,ordj,detj, trk,ordk,detk):")
    for key, count in sorted(buckets.items(), key=lambda x: -x[1])[:20]:
        tri, ordi, deti, trj, ordj, detj, trk, ordk, detk = key
        print(f"    {count:>5d}  "
              f"Ri=(tr={tri:>2d} ord={ordi} det={deti:+d})  "
              f"Rj=(tr={trj:>2d} ord={ordj} det={detj:+d})  "
              f"Rk=(tr={trk:>2d} ord={ordk} det={detk:+d})")

    # Summary: which (ord_i, ord_j) combos dominate?
    combo = {}
    for (tri, ordi, deti, trj, ordj, detj, trk, ordk, detk), count in buckets.items():
        ck = (ordi, ordj)
        combo[ck] = combo.get(ck, 0) + count
    print(f"\n  Mismatch by (order_i, order_j):")
    for ck, count in sorted(combo.items(), key=lambda x: -x[1]):
        print(f"    ({ck[0]:>2d} × {ck[1]:>2d}): {count:>5d}")


# ── Main ──

def main():
    tables = find_tables_dir()
    if not tables:
        print("irreptables not found")
        return
    files = sorted(glob.glob(os.path.join(tables, "irreps-SG=*-spin.dat")))

    # ── Part 1: Full scan with real conventions ──
    print("=" * 70)
    print("PART 1: Real 4-number SU(2) conventions (parts[12:16] only)")
    print("=" * 70)

    def get_real4(op):
        return op[2]  # su2_real tuple of 4

    for name, mul_fn, cmp_fn in REAL_CONVENTIONS:
        print(f"\n--- {name} ---")
        sgs_ok = 0
        sgs_total = 0
        total_matched = 0
        total_mismatch = 0
        total_missing = 0
        bad_sgs = []
        for f in files:
            ops = parse_spin_dat(f)
            if len(ops) < 2:
                continue
            r = test_file(ops, mul_fn, cmp_fn, get_real4)
            if r is None:
                continue
            sgs_total += 1
            if r["su2_mismatch"] == 0:
                sgs_ok += 1
            else:
                sg = os.path.basename(f).replace("irreps-SG=", "").replace("-spin.dat", "")
                r["sg"] = sg
                bad_sgs.append(r)
            total_matched += r["spatial_matched"]
            total_mismatch += r["su2_mismatch"]
            total_missing += r["spatial_missing"]
        print(f"  OK: {sgs_ok}/{sgs_total} ({sgs_ok/sgs_total*100:.1f}%)" if sgs_total else "  no data")
        print(f"  matched={total_matched} mismatched={total_mismatch} missing={total_missing}")
        if bad_sgs:
            bad_sgs.sort(key=lambda r: -r["su2_mismatch"])
            print(f"  Bad SGs ({len(bad_sgs)}):")
            for r in bad_sgs[:10]:
                pct = r["su2_mismatch"] / max(r["spatial_matched"], 1) * 100
                print(f"    SG{r['sg']:>4s} mism={r['su2_mismatch']:>5d}/{r['spatial_matched']:>5d} ({pct:.1f}%)")

    # ── Part 2: SG 221 diagnosis with flat2x2 ──
    print("\n" + "=" * 70)
    print("PART 2: SG 221 mismatch diagnosis (flat2x2)")
    print("=" * 70)

    sg221_path = os.path.join(tables, "irreps-SG=221-spin.dat")
    if os.path.exists(sg221_path):
        diagnose_sg221(sg221_path, su2_mul_flat2x2, same_up_to_sign, get_real4, max_print=20)
    else:
        print("  SG 221 spin.dat not found")

    # ── Part 3: Complex 2×2 conventions (all 8 numbers) ──
    print("\n" + "=" * 70)
    print("PART 3: Complex 2×2 conventions (parts[12:20] = 8 numbers)")
    print("=" * 70)

    def get_complex8(op):
        return (op[2], op[3])  # (su2_real_4, su2_imag_4)

    for name, mul_fn, cmp_fn in CMPLX_CONVENTIONS:
        print(f"\n--- {name} ---")
        sgs_ok = 0
        sgs_total = 0
        total_matched = 0
        total_mismatch = 0
        total_missing = 0
        bad_sgs = []
        for f in files:
            ops = parse_spin_dat(f)
            if len(ops) < 2:
                continue
            r = test_file(ops, mul_fn, cmp_fn, get_complex8)
            if r is None:
                continue
            sgs_total += 1
            if r["su2_mismatch"] == 0:
                sgs_ok += 1
            else:
                sg = os.path.basename(f).replace("irreps-SG=", "").replace("-spin.dat", "")
                r["sg"] = sg
                bad_sgs.append(r)
            total_matched += r["spatial_matched"]
            total_mismatch += r["su2_mismatch"]
            total_missing += r["spatial_missing"]
        print(f"  OK: {sgs_ok}/{sgs_total} ({sgs_ok/sgs_total*100:.1f}%)" if sgs_total else "  no data")
        print(f"  matched={total_matched} mismatched={total_mismatch} missing={total_missing}")
        if bad_sgs:
            bad_sgs.sort(key=lambda r: -r["su2_mismatch"])
            print(f"  Bad SGs ({len(bad_sgs)}):")
            for r in bad_sgs[:15]:
                pct = r["su2_mismatch"] / max(r["spatial_matched"], 1) * 100
                print(f"    SG{r['sg']:>4s} mism={r['su2_mismatch']:>5d}/{r['spatial_matched']:>5d} ({pct:.1f}%)")


if __name__ == "__main__":
    main()
