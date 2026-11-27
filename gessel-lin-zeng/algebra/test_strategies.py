import math
import random
import numpy as np
import mpmath as mp

mp.mp.dps = 120

def poly_trim(c):
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c

def apply_D(j, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return poly_trim(out)

def apply_S(L, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (L - m) * cm
    return poly_trim(out)

def add_poly(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return poly_trim(out)

def sub_poly(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] -= b[i]
    return poly_trim(out)

def scale_poly(c, s):
    return poly_trim([v * s for v in c])

def mul_poly(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return poly_trim(out)

def to_x_poly(h_coeffs):
    """Convert h(t) of degree d to A(x) = (1-x)^d h(x/(1-x)) = sum c_m x^m (1-x)^{d-m}."""
    d = len(h_coeffs) - 1
    # Build (1-x)^r table
    pows = [[1]]
    for r in range(1, d + 1):
        pows.append(mul_poly(pows[-1], [1, -1]))
    out = [0] * (d + 1)
    for m, cm in enumerate(h_coeffs):
        # cm * x^m * (1-x)^{d-m}
        pw = pows[d - m]
        for j, pj in enumerate(pw):
            out[m + j] += cm * pj
    return poly_trim(out)

def eval_mp(coeffs, z):
    acc = mp.mpc(0)
    for c in reversed(coeffs):
        acc = acc * z + c
    return acc

# === TEST 1: High-precision (40-digit mpmath!) check of Q_k(t, y) and Q_k^A(x, y) at q=1/2, 1/4, 1/8! ===
print("=== TEST 1: 40-digit mpmath UHP stability of Q_k^A(x, y) and Q_k(t, y) ===")
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 16):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

A = {(0, 0): [1]}
for k in range(1, 16):
    for i in range(k + 1):
        t1 = apply_S(3*k - 3 - i, A[(k-1, i)]) if i <= k - 1 else [0]
        t2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S(3*k - 2 - i, add_poly(t1, t2))

# Check k=15 specifically at the point z = -0.95 + 1e-4j that triggered np.roots!
for name, pdict, test_zs in [
    ("htilde(t)", ht, [mp.mpc(-0.95, 1e-4), mp.mpc(-0.5, 1e-3), mp.mpc(0, 1e-4), mp.mpc(1, 1e-4), mp.mpc(-2, 1e-3), mp.mpc(10, 0.1)]),
    ("A(x)", A, [mp.mpc(-20, 1e-3), mp.mpc(-1, 1e-4), mp.mpc(-0.1, 1e-4), mp.mpc(0, 1e-4), mp.mpc(1, 1e-4), mp.mpc(10, 0.1)]),
]:
    for q_val in [mp.mpf('0.5'), mp.mpf('0.25'), mp.mpf('0.125')]:
        max_im = mp.mpf('-1e9')
        worst = None
        for k in [5, 10, 15]:
            w = [q_val ** (i * (i - 1) // 2) for i in range(k + 1)]
            for zv in test_zs:
                yc = [w[i] * eval_mp(pdict[(k, i)], zv) for i in range(k, -1, -1)]
                r = mp.polyroots(yc, maxsteps=500, extraprec=60)
                mi = max(z.imag for z in r)
                if mi > max_im:
                    max_im = mi
                    worst = (k, complex(zv))
        print(f"  {name:10s} q={float(q_val):.3f}: max Im(y) = {float(max_im):+.4e} at {worst}")

# === TEST 2: Check M_{b,q}(F) = D_b F(t,y) + y F(t, q y) on Wild prod (D_{b_j} + y)(1) with q-weights! ===
print("\n=== TEST 2: Does Q_k(t, y) = sum q^{i(i-1)/2} q_{k,i}(t) y^i stay real-stable for WILD b_j? ===")
for name, bs in [
    ("Const [2]*6", [2]*6),
    ("Wild [2,100,2,100,2,100]", [2,100,2,100,2,100]),
    ("Wild [100,2,100,2,100,2]", [100,2,100,2,100,2]),
    ("Drop [100,50,25,12,6,3]", [100,50,25,12,6,3]),
]:
    row = [[1]]
    for bj in bs:
        k = len(row)
        new_row = []
        for i in range(k + 1):
            t1 = apply_D(bj, row[i]) if i <= k - 1 else [0]
            t2 = row[i - 1] if i >= 1 else [0]
            new_row.append(add_poly(t1, t2))
        row = new_row
    # Test UHP stability of sum q^{i(i-1)/2} row[i](t) y^i
    k = len(bs)
    test_zs = [mp.mpc(re, im) for re in [-5, -1.2, -0.8, -0.5, -0.1, 0, 0.1, 1, 5] for im in [1e-4, 0.05, 1.0]]
    res_str = []
    for q_val in [mp.mpf('0.5'), mp.mpf('0.25'), mp.mpf('0.1'), mp.mpf('0.01')]:
        w = [q_val ** (i * (i - 1) // 2) for i in range(k + 1)]
        max_im = mp.mpf('-1e9')
        for zv in test_zs:
            yc = [w[i] * eval_mp(row[i], zv) for i in range(k, -1, -1)]
            r = mp.polyroots(yc, maxsteps=500, extraprec=60)
            mi = max(z.imag for z in r)
            if mi > max_im: max_im = mi
        res_str.append(f"q={float(q_val)}:{'OK' if max_im < 1e-12 else f'FAIL({float(max_im):.2e})'}")
    print(f"  {name:25s}: " + ", ".join(res_str))

# === TEST 3: Check Turan (p_i^2 - p_{i-1} p_{i+1}) and M^2 - XY in x-coeffs and on R for WILD b_j! ===
print("\n=== TEST 3: Check coefficientwise positivity in x of p_i^2 - p_{i-1}p_{i+1} and (D p_{i-1})^2 - p_{i-2}(D^2 p_i) ===")
for name, bs in [
    ("GLZ inner [2,3,4,5,6,7]", [2,3,4,5,6,7]),
    ("Const [2]*6", [2]*6),
    ("Wild [2,100,2,100,2,100]", [2,100,2,100,2,100]),
    ("Wild [100,2,100,2,100,2]", [100,2,100,2,100,2]),
    ("Drop [100,50,25,12,6,3]", [100,50,25,12,6,3]),
]:
    row = [[1]]
    turan_x_pos = True
    m2xy_x_pos = True
    m2xy_R_pos = True
    for step, bj in enumerate(bs, 1):
        # Check M^2 - XY before updating! M = D_{bj} row[i-1], X = row[i-2], Y = D_{bj}^2 row[i]
        if step >= 3:
            for i in range(2, step):
                M = apply_D(bj, row[i-1])
                X = row[i-2]
                Y = apply_D(bj, apply_D(bj, row[i]))
                diff_t = sub_poly(mul_poly(M, M), mul_poly(X, Y))
                diff_x = to_x_poly(diff_t)
                if min(diff_x) < 0:
                    m2xy_x_pos = False
                # Check on (-1, 0) numerically
                for tv in [-0.99, -0.9, -0.7, -0.5, -0.3, -0.1, -0.01]:
                    val = sum(c * (tv**m) for m, c in enumerate(diff_t))
                    if val < 0:
                        m2xy_R_pos = False
        new_row = []
        for i in range(step + 1):
            t1 = apply_D(bj, row[i]) if i <= step - 1 else [0]
            t2 = row[i - 1] if i >= 1 else [0]
            new_row.append(add_poly(t1, t2))
        row = new_row
        for i in range(1, step):
            tur_t = sub_poly(mul_poly(row[i], row[i]), mul_poly(row[i-1], row[i+1]))
            tur_x = to_x_poly(tur_t)
            if min(tur_x) < 0:
                turan_x_pos = False
    print(f"  {name:25s}: turan_x>=0: {str(turan_x_pos):5s} | M^2-XY_x>=0: {str(m2xy_x_pos):5s} | M^2-XY on (-1,0)>=0: {str(m2xy_R_pos):5s}")

# === TEST 4: Check coefficientwise positivity in x of 2x2 minors p_i p_j - p_{i-1} p_{j+1} for WILD! ===
print("\n=== TEST 4: Check ALL 2x2 minors p_i(x) p_j(x) - p_{i-1}(x) p_{j+1}(x) >= 0 (coeffwise in x) ===")
for name, bs in [
    ("GLZ full D_j(D_j+y) k=7", "GLZ"),
    ("Wild [2,100,2,100,2,100]", [2,100,2,100,2,100]),
    ("Wild [100,2,100,2,100,2]", [100,2,100,2,100,2]),
    ("Drop [100,50,25,12,6,3]", [100,50,25,12,6,3]),
]:
    if bs == "GLZ":
        row = [ht[(7, i)] for i in range(8)]
    else:
        row = [[1]]
        for bj in bs:
            k = len(row)
            row = [add_poly(apply_D(bj, row[i]) if i <= k-1 else [0], row[i-1] if i >= 1 else [0]) for i in range(k+1)]
    xrow = [to_x_poly(p) for p in row]
    ok_2x2 = True
    worst = None
    for i in range(1, len(xrow)):
        for j in range(i, len(xrow) - 1):
            det = sub_poly(mul_poly(xrow[i], xrow[j]), mul_poly(xrow[i-1], xrow[j+1]))
            if min(det) < 0:
                ok_2x2 = False
                if worst is None: worst = (i, j, min(det))
    print(f"  {name:25s}: all 2x2 minors coeffwise >= 0 in x: {ok_2x2} (worst={worst})")
