import math
import random
import numpy as np
import sympy as sp

# === TEST 1: Exact sympy check of whether d/dt and T ACTUALLY failed in test_exact_turan_and_firstletter.py ===
print("=== TEST 1: Exact sympy verification of d/dt and T on Turan ===")
t = sp.Symbol('t')
# Take the exact worst T(a=1, b=10) example from before:
# r0=[-0.969, -0.964, -0.894, -0.809, -0.407], r1=[-0.968, -0.942, -0.837, -0.734], r2=[-0.946, -0.874, -0.829]
def check_exact_example(name, r0, r1, r2, mode, a=1, b=10):
    p0 = sp.Integer(1)
    for r in r0: p0 *= (t - sp.Rational(r))
    p1 = sp.Integer(1)
    for r in r1: p1 *= (t - sp.Rational(r))
    p2 = sp.Integer(1)
    for r in r2: p2 *= (t - sp.Rational(r))
    p0, p1, p2 = sp.expand(p0), sp.expand(p1), sp.expand(p2)
    # Find max of (p0*p2)/(p1^2) over all real t using sympy real_roots of critical poly!
    N = sp.Poly(sp.expand(sp.diff(p0*p2, t)*p1 - 2*p0*p2*sp.diff(p1, t)), t)
    crit = [r.evalf(30) for r in N.real_roots()]
    max_in = sp.Float(1, 30)
    for tc in crit:
        v1 = p1.subs(t, tc)
        if abs(v1) > 1e-25:
            val = (p0*p2).subs(t, tc) / (v1**2)
            if val > max_in: max_in = val
    # Scale p2 by 1/max_in
    p2_s = p2 / max_in
    if mode == "deriv":
        q0, q1, q2 = sp.diff(p0, t), sp.diff(p1, t), sp.diff(p2_s, t)
    else:
        D = lambda f: sp.expand(t*(1+t)*sp.diff(f, t) + (a + b*t)*f)
        q0, q1, q2 = D(p0), D(p1), D(p2_s)
    # Evaluate q1^2 - q0*q2 at its critical points
    N2 = sp.Poly(sp.expand(sp.diff(q0*q2, t)*q1 - 2*q0*q2*sp.diff(q1, t)), t)
    crit2 = [r.evalf(30) for r in N2.real_roots()]
    max_out = (sp.Poly(q0, t).LC() * sp.Poly(q2, t).LC()) / (sp.Poly(q1, t).LC()**2)
    for tc in crit2:
        v1 = q1.subs(t, tc)
        if abs(v1) > 1e-25:
            val = (q0*q2).subs(t, tc) / (v1**2)
            if val > max_out: max_out = val
    print(f"  {name}: exact max_in={float(max_in):.8f}, exact max_out={float(max_out):.8f} ({'GENUINE FAIL!' if max_out > 1 + 1e-9 else 'WAS FLOAT NOISE -> PASS!'})")

check_exact_example("T(1,10) ex", ['-0.969', '-0.964', '-0.894', '-0.809', '-0.407'],
                    ['-0.968', '-0.942', '-0.837', '-0.734'], ['-0.946', '-0.874', '-0.829'], "T", 1, 10)
check_exact_example("T(1,3) ex", ['-0.589', '-0.425', '-0.346', '-0.342', '-0.249'],
                    ['-0.427', '-0.365', '-0.344', '-0.341'], ['-0.412', '-0.35', '-0.344'], "T", 1, 3)
check_exact_example("deriv ex", ['-0.969', '-0.964', '-0.894', '-0.809', '-0.407'],
                    ['-0.968', '-0.942', '-0.837', '-0.734'], ['-0.946', '-0.874', '-0.829'], "deriv")

# === TEST 2: Check F_y(t) = ht_{k,i-1}(t) + 2*y*ht_{k,i}(t) + y^2*ht_{k,i+1}(t) real-rooted in t for ALL real y! ===
print("\n=== TEST 2: Check if F_y(t) = p_{i-1}(t) + 2*y*p_i(t) + y^2*p_{i+1}(t) is REAL-ROOTED in t for ALL real y! ===")
def apply_D_int(j, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return out
def add_int(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return out

ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 13):
    for i in range(k + 1):
        t1 = apply_D_int(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D_int(k + 1, add_int(t1, t2))

def is_rr_np(p):
    r = np.roots(p[::-1])
    return np.max(np.abs(r.imag)) < 1e-6

ys = [-100.0, -10.0, -3.0, -2.0, -1.5, -1.0, -0.8, -0.5, -0.2, -0.05, 0.05, 0.2, 0.5, 1.0, 2.0, 10.0, 100.0]
for k in range(2, 13):
    bad = None
    for i in range(1, k):
        p0 = np.array([float(c) for c in ht[(k, i-1)]])
        p1 = np.array([float(c) for c in ht[(k, i)]])
        p2 = np.array([float(c) for c in ht[(k, i+1)]])
        for yv in ys:
            fy = p0.copy()
            fy[:len(p1)] += 2.0 * yv * p1
            fy[:len(p2)] += (yv * yv) * p2
            if not is_rr_np(fy):
                bad = (i, yv)
                break
        if bad: break
    print(f"  GLZ k={k:2d}: F_y(t) real-rooted in t for all y in test set: {'PASS ALL!' if bad is None else f'FAIL at i={bad[0]}, y={bad[1]}'}")

# Also test for Wild [2, 100, 2, 100, 2, 100]!
row = [[1]]
for bj in [2, 100, 2, 100, 2, 100]:
    k = len(row)
    row = [add_int(apply_D_int(bj, row[i]) if i <= k-1 else [0], row[i-1] if i >= 1 else [0]) for i in range(k+1)]
bad_w = None
for i in range(1, len(row) - 1):
    p0 = np.array([float(c) for c in row[i-1]])
    p1 = np.array([float(c) for c in row[i]])
    p2 = np.array([float(c) for c in row[i+1]])
    for yv in ys:
        fy = p0.copy()
        fy[:len(p1)] += 2.0 * yv * p1
        fy[:len(p2)] += (yv * yv) * p2
        if not is_rr_np(fy):
            bad_w = (i, yv)
            break
    if bad_w: break
print(f"  Wild [2,100,2,100,2,100]: F_y(t) real-rooted in t for all y: {'PASS ALL!' if bad_w is None else f'FAIL at i={bad_w[0]}, y={bad_w[1]}'}")

# Also check for GLZ: is ALL-DEGREE H_k(t, y) = sum ht_{k,i}(t) y^i real-rooted in t for ALL y >= 0?!
# AND what about (1+y)^k-weighted or plain sum ht_{k,i}(t) y^i for y >= 0?!
print("\n=== TEST 3: Is plain H_k(t, y) = sum ht_{k,i}(t) y^i real-rooted in t for all y > 0? And for Wild? ===")
for k in range(2, 12):
    bad = None
    for yv in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0, 1000.0]:
        f = np.zeros(len(ht[(k, 0)]), dtype=float)
        for i in range(k + 1):
            pi = np.array([float(c) for c in ht[(k, i)]])
            f[:len(pi)] += (yv ** i) * pi
        if not is_rr_np(f):
            bad = yv
            break
    print(f"  GLZ k={k:2d}: H_k(t, y) real-rooted in t for y > 0: {'PASS ALL!' if bad is None else f'FAIL at y={bad}'}")

# Check roots of H_k(t, y) for y > 0: do they stay in (-inf, 0) (not just (-1, 0))?
for yv in [0.1, 1.0, 10.0, 100.0]:
    f = np.zeros(len(ht[(4, 0)]), dtype=float)
    for i in range(5):
        pi = np.array([float(c) for c in ht[(4, i)]])
        f[:len(pi)] += (yv ** i) * pi
    r = np.sort(np.roots(f[::-1]).real)
    print(f"  H_4(t, y={yv:5.1f}) roots: {[round(float(z), 4) for z in r]}")
