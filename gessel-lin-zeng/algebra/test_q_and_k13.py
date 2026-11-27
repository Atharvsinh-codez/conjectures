import math
import random
import numpy as np

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
    for i in range(len(a)):
        out[i] += a[i]
    for i in range(len(b)):
        out[i] += b[i]
    return poly_trim(out)

def scale_poly(c, s):
    return poly_trim([v * s for v in c])

def exact_strict_interlaces(p, q):
    dp, dq = len(p) - 1, len(q) - 1
    if p[-1] <= 0 or q[-1] <= 0:
        return False, "nonpos_lc"
    if dq == 0:
        return dp == 0, "deg0"
    if dq == dp:
        lcp, lcq = p[-1], q[-1]
        r1 = [lcq * p[m] - lcp * q[m] for m in range(dq)]
        while len(r1) > 1 and r1[-1] == 0:
            r1.pop()
        if len(r1) - 1 != dq - 1 or r1[-1] <= 0:
            return False, "eq_lc_fail"
        r0 = q[:]
    elif dq == dp + 1:
        r0 = q[:]
        r1 = p[:]
    else:
        return False, f"deg_jump({dp}->{dq})"
    cur0, cur1 = r0[:], r1[:]
    while len(cur1) > 1:
        d0, d1 = len(cur0) - 1, len(cur1) - 1
        if d0 != d1 + 1:
            return False, f"deg_drop({d0},{d1})"
        a0, a1 = cur0[-1], cur1[-1]
        tmp = [a1 * cur0[m] for m in range(d0 + 1)]
        for m in range(d1 + 1):
            tmp[m + 1] -= a0 * cur1[m]
        b = tmp[d1]
        nxt = [-(a1 * tmp[m] - b * cur1[m]) for m in range(d1)]
        while len(nxt) > 1 and nxt[-1] == 0:
            nxt.pop()
        if len(nxt) - 1 != d1 - 1 or nxt[-1] <= 0:
            return False, f"sturm_fail_deg{d1-1}"
        g = math.gcd(*[abs(v) for v in nxt])
        if g > 1:
            nxt = [v // g for v in nxt]
        cur0, cur1 = cur1, nxt
    return True, "ok"

def is_real_rooted(p):
    d = len(p) - 1
    if d <= 1:
        return True
    if p[-1] < 0:
        p = [-v for v in p]
    dp = [m * p[m] for m in range(1, d + 1)]
    return exact_strict_interlaces(dp, p)[0]

# Compute A_{k,i} and htilde_{k,i} up to k=15
A = {(0, 0): [1]}
for k in range(1, 16):
    for i in range(k + 1):
        t1 = apply_S(3*k - 3 - i, A[(k-1, i)]) if i <= k - 1 else [0]
        t2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S(3*k - 2 - i, add_poly(t1, t2))

ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 16):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

def eval_poly_c(coeffs, z):
    acc = 0j
    for c in reversed(coeffs):
        acc = acc * z + c
    return acc

# === TEST 1: Adversarial UHP stability of Q_k(t, y) = sum q^{i(i-1)/2} htilde_{k,i}(t) y^i AND in x ===
print("=== TEST 1: Adversarial UHP stability of q-weighted Q_k(t, y) and Q_k^A(x, y) ===")
# Dense grid of t/x near real axis (including negative, positive, near -1, near 0, large)
pts = []
for re in [-50, -10, -3, -1.5, -1.05, -1.0, -0.95, -0.7, -0.5, -0.3, -0.1, -0.01, 0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 50.0]:
    for im in [1e-6, 1e-4, 1e-2, 0.1, 1.0, 10.0]:
        pts.append(complex(re, im))

for name, pdict in [("A_{k,i}(x)", A), ("htilde_{k,i}(t)", ht)]:
    for q_val in [0.5, 0.25, 0.1, 0.01]:
        first_fail = None
        for k in range(2, 16):
            w = [q_val ** (i * (i - 1) // 2) for i in range(k + 1)]
            max_im = -1e9
            worst_z = None
            for zv in pts:
                yc = [w[i] * eval_poly_c(pdict[(k, i)], zv) for i in range(k, -1, -1)]
                if abs(yc[0]) == 0: continue
                # Normalize by max abs to prevent overflow
                mx = max(abs(c) for c in yc)
                yc = [c / mx for c in yc]
                r = np.roots(yc)
                mi = max(z.imag for z in r)
                if mi > max_im:
                    max_im = mi
                    worst_z = zv
            if max_im > 1e-7:
                first_fail = (k, max_im, worst_z)
                break
        print(f"  {name:16s} q={q_val:4.2f}: {'ALL PASS k=2..15!' if first_fail is None else f'FAIL at k={first_fail[0]} max_Im={first_fail[1]:.3e} at z={first_fail[2]}'}")

# === TEST 2: Check EXACT Sturm for q-weighted Q_k(t, y) at t = fixed real values (including t=1, t=10, t=-2, t=-0.5) ===
print("\n=== TEST 2: Exact integer Sturm check of Q_k(t, y) in y for fixed real t (q=1/4, 1/2) ===")
for q_den in [2, 4, 8]:
    # Multiply by q_den^{k(k-1)/2} to make exact ints!
    bad = None
    for k in range(2, 14):
        E = k * (k - 1) // 2
        for tnum, tden in [(-3,1), (-2,1), (-1,2), (-1,10), (0,1), (1,10), (1,1), (10,1)]:
            dmax = 2*k - 1
            y_poly = []
            for i in range(k + 1):
                val = sum(c * (tnum**m) * (tden**(dmax - m)) for m, c in enumerate(ht[(k, i)]))
                y_poly.append(val * (q_den ** (E - i*(i-1)//2)))
            y_poly = poly_trim(y_poly)
            if not is_real_rooted(y_poly):
                bad = (k, tnum/tden)
                break
        if bad: break
    print(f"  htilde q=1/{q_den}: {'ALL PASS k=2..13' if bad is None else f'FAIL at k={bad[0]}, t={bad[1]}'}")

# === TEST 3: Check 1-step operator F(t,y) -> D_b F(t,y) + y F(t, q*y) on REAL-STABLE F(t, y)! ===
print("\n=== TEST 3: Test if M_{b,q}(F) = D_b F(t, y) + y F(t, q y) preserves real stability! ===")
# Take a simple real-stable F(t, y) with t-roots in (-1, 0) when y=0, e.g. F(t, y) = (1 + 2t) + y (root at t = -(1+y)/2)
# Check M_{b,q}(F) = D_b(1 + 2t + y) + y(1 + 2t + q*y) = (1 + (b+2)t + (2b+2)t^2) + (2 + (b+2)t)y + q*y^2?
# Wait: D_b(y) = (1 + b*t)*y!
# So M_{b,q}(1 + 2t + y) = (1 + (b+4)t + 2(b+1)t^2) + (2 + (b+2)t)y + q*y^2!
# Compute discriminant in y of M_{b,q}(1+2t+y)!
import sympy as sp
tsym, qsym, bsym = sp.symbols('t q b')
F0 = 1 + 2*tsym + sp.Symbol('y')
y = sp.Symbol('y')
DbF0 = tsym*(1+tsym)*sp.diff(F0, tsym) + (1 + bsym*tsym)*F0
M1 = sp.expand(DbF0 + y * F0.subs(y, qsym*y))
pM1 = sp.Poly(M1, y)
disc_y = sp.factor(pM1.discriminant())
print("  M_{b,q}(1 + 2t + y) =", M1)
print("  disc_y =", disc_y)
for bv in [2, 3, 5, 10]:
    for qv in [sp.Rational(1,2), sp.Rational(1,4), sp.Rational(1,10)]:
        d_sub = sp.expand(disc_y.subs({bsym: bv, qsym: qv}))
        # Check if d_sub >= 0 for ALL real t!
        r = np.roots([float(c) for c in sp.Poly(d_sub, tsym).all_coeffs()])
        real_r = [z.real for z in r if abs(z.imag) < 1e-9]
        # evaluate at t=100, t=-100, t=0, t=-0.5
        min_val = min(float(d_sub.subs(tsym, tv)) for tv in [-10, -2, -0.5, 0, 0.5, 2, 10])
        print(f"    b={bv:2d}, q={str(qv):4s}: disc_y = {d_sub} | real roots={len(real_r)} | sample_min={min_val:.3f}")

# === TEST 4: Test extreme weights w_j in prod_{j=1}^{13} D_j(D_j + w_j y)(1) (testing per-subset barrier at k=13!) ===
print("\n=== TEST 4: Stress test prod_{j=1}^{13} D_{j+1}(D_{j+1} + w_j y)(1) with extreme weights w_j ===")
for name, ws in [
    ("w=[1]*6 + [100]*7", [1]*6 + [100]*7),
    ("w=[100]*6 + [1]*7", [100]*6 + [1]*7),
    ("w=[1]*6 + [10000]*7", [1]*6 + [10000]*7),
    ("w=[1,1000,1,1000,...]", [1 if j%2==0 else 1000 for j in range(13)]),
]:
    row = [[1]]
    ok_rr = True
    ok_int = True
    first_fail = None
    for step in range(1, 14):
        cj = step + 1
        wj = ws[step - 1]
        new_row = []
        for i in range(step + 1):
            t1 = apply_D(cj, row[i]) if i <= step - 1 else [0]
            t2 = scale_poly(row[i - 1], wj) if i >= 1 else [0]
            new_row.append(apply_D(cj, add_poly(t1, t2)))
        if step == 13:
            for i in range(1, 14):
                ok, msg = exact_strict_interlaces(new_row[i], new_row[i-1])
                if not ok:
                    ok_int = False
                    if first_fail is None: first_fail = (step, i, msg)
        row = new_row
    print(f"  {name:25s} (k=13): h_i<<h_{{i-1}} = {ok_int} (fail={first_fail})")
