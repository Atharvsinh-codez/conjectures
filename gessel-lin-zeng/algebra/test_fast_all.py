import math
import random
import numpy as np
from itertools import combinations

def poly_trim(c):
    while len(c) > 1 and c[-1] == 0: c.pop()
    return c

def apply_D_int(j, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return poly_trim(out)

def add_int(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return poly_trim(out)

def sub_int(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] -= b[i]
    return poly_trim(out)

def mul_int(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return poly_trim(out)

def eval_grid(p, grid):
    acc = np.zeros_like(grid, dtype=float)
    for c in reversed(p):
        acc = acc * grid + float(c)
    return acc

def is_rr_np(p):
    r = np.roots(p[::-1])
    return np.max(np.abs(r.imag)) < 1e-6

# === TEST 1: Fast check of F_y(t) = ht_{k,i-1}(t) + 2*y*ht_{k,i}(t) + y^2*ht_{k,i+1}(t) real-rooted in t for ALL real y! ===
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 13):
    for i in range(k + 1):
        t1 = apply_D_int(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D_int(k + 1, add_int(t1, t2))

print("=== TEST 1: Is F_y(t) = p_{i-1} + 2*y*p_i + y^2*p_{i+1} real-rooted in t for ALL real y? ===")
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
    print(f"  GLZ k={k:2d}: {'PASS ALL!' if bad is None else f'FAIL at i={bad[0]}, y={bad[1]}'}")

# === TEST 2: Is plain H_k(t, y) = sum ht_{k,i}(t) y^i real-rooted in t for ALL y > 0? And for Wild? ===
print("\n=== TEST 2: Is plain H_k(t, y) = sum ht_{k,i}(t) y^i real-rooted in t for ALL y > 0? ===")
for k in range(2, 13):
    bad = None
    for yv in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0, 1000.0]:
        f = np.zeros(len(ht[(k, 0)]), dtype=float)
        for i in range(k + 1):
            pi = np.array([float(c) for c in ht[(k, i)]])
            f[:len(pi)] += (yv ** i) * pi
        if not is_rr_np(f):
            bad = yv
            break
    print(f"  GLZ k={k:2d}: H_k(t, y) real-rooted in t for y > 0: {'PASS ALL!' if bad is None else f'FAIL at y={bad}'}")

# Check H_k(t, y) for Wild [2, 100, 2, 100, 2, 100] and y > 0!
row = [[1]]
for bj in [2, 100, 2, 100, 2, 100]:
    k = len(row)
    row = [add_int(apply_D_int(bj, row[i]) if i <= k-1 else [0], row[i-1] if i >= 1 else [0]) for i in range(k+1)]
bad_w = None
for yv in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0, 1000.0]:
    f = np.zeros(len(row[0]), dtype=float)
    for i in range(len(row)):
        pi = np.array([float(c) for c in row[i]])
        f[:len(pi)] += (yv ** i) * pi
    if not is_rr_np(f):
        bad_w = yv
        break
print(f"  Wild [2,100,2,100,2,100]: P_6(t, y) real-rooted in t for y > 0: {'PASS ALL!' if bad_w is None else f'FAIL at y={bad_w}'}")

# === TEST 3: Check pairwise Rayleigh R = p_{I U a} p_{I U b} - p_I p_{I U ab} >= 0 on (-1, 0) for Full GLZ k=6, 7 ===
print("\n=== TEST 3: Check pairwise Rayleigh R >= 0 on (-1, 0) for Full GLZ k=6, 7 ===")
t_int = np.linspace(-1.0, 0.0, 601)
for k in [6, 7]:
    pS = {frozenset(): [1]}
    for j in range(1, k + 1):
        nxt = {}
        for S, poly in pS.items():
            nxt[S] = apply_D_int(j + 1, apply_D_int(j + 1, poly))
            nxt[S | {j}] = apply_D_int(j + 1, poly)
        pS = nxt
    fails = 0
    tot = 0
    worst = 999.0
    for a, b in combinations(range(1, k + 1), 2):
        rem = [x for x in range(1, k + 1) if x not in (a, b)]
        for r in range(len(rem) + 1):
            for I_tup in combinations(rem, r):
                I = frozenset(I_tup)
                tot += 1
                R = sub_int(mul_int(pS[I | {a}], pS[I | {b}]), mul_int(pS[I], pS[I | {a, b}]))
                m = float(np.min(eval_grid(R, t_int)))
                if m < worst: worst = m
                if m < -1e-6: fails += 1
    print(f"  Full GLZ k={k}: Rayleigh R >= 0 on (-1,0): fails = {fails}/{tot} | worst = {worst:+.3e}")

# === TEST 4: Check if q_i^2 >= q_{i-1} q_{i+1} on (-1, 0) EVER FAILS for weighted Wild prod (D_{b_j} + w_j y)(1)! ===
print("\n=== TEST 4: Adversarial test of q_i^2 >= q_{i-1} q_{i+1} on (-1, 0) for weighted Wild prod (D_{b_j} + w_j y)(1) ===")
random.seed(999)
fails_turan_w = 0
for trial in range(100):
    bs = [random.choice([2, 3, 10, 100, 1000]) for _ in range(6)]
    ws = [random.choice([1, 10, 1000, 100000]) for _ in range(6)]
    row = [[1]]
    for bj, wj in zip(bs, ws):
        k = len(row)
        new_row = []
        for i in range(k + 1):
            t1 = apply_D_int(bj, row[i]) if i <= k - 1 else [0]
            t2 = [v * wj for v in row[i - 1]] if i >= 1 else [0]
            new_row.append(add_int(t1, t2))
        row = new_row
    for i in range(1, 6):
        diff = sub_int(mul_int(row[i], row[i]), mul_int(row[i-1], row[i+1]))
        if np.min(eval_grid(diff, t_int)) < -1e-6:
            fails_turan_w += 1
            break
print(f"  100 wild (b_j, w_j) trials (k=6): fails of q_i^2 >= q_{{i-1}} q_{{i+1}} on (-1, 0) = {fails_turan_w}/100")

# === TEST 5: Exact mpmath check of whether d/dt ACTUALLY failed Turan on the test_exact_turan_and_firstletter example! ===
print("\n=== TEST 5: Exact 80-digit mpmath check of T(1,10) and d/dt on that tight example ===")
import mpmath as mp
mp.mp.dps = 80
def poly_from_roots_mp(rs):
    p = [mp.mpf(1)]
    for r in rs:
        # multiply by (t - r)
        nxt = [mp.mpf(0)] * (len(p) + 1)
        for m, c in enumerate(p):
            nxt[m] -= c * r
            nxt[m + 1] += c
        p = nxt
    return p
def eval_mp(p, x):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc
def deriv_mp(p):
    return [m * p[m] for m in range(1, len(p))]
def apply_T_mp(a, b, p):
    out = [mp.mpf(0)] * (len(p) + 1)
    for m, c in enumerate(p):
        out[m] += (m + a) * c
        out[m + 1] += (m + b) * c
    return out
def mul_mp(a, b):
    out = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b): out[i + j] += ai * bj
    return out
def sub_mp(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] -= b[i]
    return out
def max_ratio_mp(p0, p1, p2):
    p02 = mul_mp(p0, p2)
    N = sub_mp(mul_mp(deriv_mp(p02), p1), [2*c for c in mul_mp(p02, deriv_mp(p1))])
    roots = mp.polyroots(list(reversed(N)), maxsteps=200, extraprec=30)
    mr = (p0[-1]*p2[-1]) / (p1[-1]**2)
    for z in roots:
        if abs(z.imag) < mp.mpf('1e-40'):
            v1 = eval_mp(p1, z.real)
            if abs(v1) > mp.mpf('1e-50'):
                v02 = eval_mp(p02, z.real)
                if v02 > 0:
                    mr = max(mr, v02 / (v1*v1))
    return mr

r0 = [mp.mpf(s) for s in ['-0.969', '-0.964', '-0.894', '-0.809', '-0.407']]
r1 = [mp.mpf(s) for s in ['-0.968', '-0.942', '-0.837', '-0.734']]
r2 = [mp.mpf(s) for s in ['-0.946', '-0.874', '-0.829']]
p0, p1, p2 = poly_from_roots_mp(r0), poly_from_roots_mp(r1), poly_from_roots_mp(r2)
mr_in = max_ratio_mp(p0, p1, p2)
p2s = [c / mr_in for c in p2]
mr_d = max_ratio_mp(deriv_mp(p0), deriv_mp(p1), deriv_mp(p2s))
mr_T = max_ratio_mp(apply_T_mp(1, 10, p0), apply_T_mp(1, 10, p1), apply_T_mp(1, 10, p2s))
print(f"  80-digit exact: mr_in={float(mr_in):.8f} -> mr_deriv={float(mr_d):.8f}, mr_T(1,10)={float(mr_T):.8f}")
