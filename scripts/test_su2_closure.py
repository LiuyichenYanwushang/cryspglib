#!/usr/bin/env python3
"""Validate Bilbao SU(2) convention by testing group closure.

For every pair of spin operations (i,j) in a spin.dat file,
compute the spatial product {R_i|t_i}{R_j|t_j} and the SU(2)
product U_i·U_j.  Verify that U_i·U_j equals ±U_k for the
spin operation k with the matching spatial product.

This is a closure test: it does not assume any specific SU(2)
convention — it tries multiple conventions and reports which
ones pass.
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
            if len(parts) >= 16:
                rot = tuple(int(x) for x in parts[0:9])
                trans = tuple(float(x) for x in parts[9:12])
                su2 = tuple(float(x) for x in parts[12:16])
                ops.append((rot, trans, su2))
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
    r1, t1, _ = op1
    r2, t2, _ = op2
    r = matmul3(r1, r2)
    rt2 = matvec3(r1, t2)
    t = tuple(mod1(t1[i] + rt2[i]) for i in range(3))
    return r, t


def find_spatial_product(ops, r, t):
    return [idx for idx, (rr, tt, _) in enumerate(ops)
            if rr == r and same_trans(tt, t)]


# ── SU(2) composition conventions to test ──

def su2_mul_flat2x2(a, b):
    """2×2 matrix multiply, flattened [a,b,c,d] → [[a,b],[c,d]]."""
    return (
        a[0]*b[0] + a[1]*b[2],
        a[0]*b[1] + a[1]*b[3],
        a[2]*b[0] + a[3]*b[2],
        a[2]*b[1] + a[3]*b[3],
    )


def su2_mul_quat_wxyz(a, b):
    """Quaternion multiply (w,x,y,z)."""
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return (
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
    )


def su2_mul_pauli(a, b):
    """SU(2) as [Re(α), Im(α), Re(β), Im(β)] for [[α,β],[-β*,α*]]."""
    ar, ai, br, bi = a
    cr, ci, dr, di = b
    # [[ar+i·ai, br+i·bi], [-br+i·bi, ar-i·ai]] · [[cr+i·ci, dr+i·di], [-dr+i·di, cr-i·ci]]
    return (
        ar*cr - ai*ci - br*dr - bi*di,   # Re(α_result)
        ar*ci + ai*cr - br*di + bi*dr,   # Im(α_result)
        ar*dr - ai*di + br*cr + bi*ci,   # Re(β_result)
        ar*di + ai*dr + br*ci - bi*cr,   # Im(β_result)
    )


CONVENTIONS = {
    "flat2x2": su2_mul_flat2x2,
    "quat-wxyz": su2_mul_quat_wxyz,
    "pauli": su2_mul_pauli,
}


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


def test_file(path, mul_fn, max_fail=10):
    ops = parse_spin_dat(path)
    if len(ops) < 2:
        return None
    failures = []
    total = 0
    sign_p = 0
    sign_m = 0
    for i, op_i in enumerate(ops):
        for j, op_j in enumerate(ops):
            total += 1
            r, t = compose_seitz(op_i, op_j)
            ks = find_spatial_product(ops, r, t)
            if not ks:
                failures.append((i, j, "no spatial match"))
                if len(failures) >= max_fail:
                    break
                continue
            qij = mul_fn(op_i[2], op_j[2])
            ok = False
            for k in ks:
                s = same_up_to_sign(qij, ops[k][2])
                if s != 0:
                    ok = True
                    if s > 0: sign_p += 1
                    else: sign_m += 1
                    break
            if not ok:
                failures.append((i, j, f"SU2 mismatch; spatial→{ks}"))
                if len(failures) >= max_fail:
                    break
        if len(failures) >= max_fail:
            break
    rate = (total - len(failures)) / total * 100 if total > 0 else 0
    return {"total": total, "+": sign_p, "-": sign_m,
            "failures": len(failures), "rate": rate,
            "examples": failures[:3]}


def main():
    tables = find_tables_dir()
    if not tables:
        print("irreptables not found")
        return
    files = sorted(glob.glob(os.path.join(tables, "irreps-SG=*-spin.dat")))
    for name, mul_fn in CONVENTIONS.items():
        print(f"\n=== Testing convention: {name} ===")
        ok = 0
        total_sg = 0
        worst = None
        for f in files:
            r = test_file(f, mul_fn, max_fail=3)
            if r is None:
                continue
            total_sg += 1
            if r["failures"] == 0:
                ok += 1
            else:
                if worst is None or r["rate"] < worst["rate"]:
                    sg = os.path.basename(f).replace("irreps-SG=", "").replace("-spin.dat", "")
                    r["sg"] = sg
                    worst = r
        print(f"  Passed: {ok}/{total_sg} ({ok/total_sg*100:.1f}%)")
        if worst:
            print(f"  Worst: SG{worst.get('sg','?')} rate={worst['rate']:.1f}% fails={worst['failures']}/{worst['total']} +{worst['+']} -{worst['-']}")
            for ex in worst["examples"]:
                print(f"    {ex}")


if __name__ == "__main__":
    main()
