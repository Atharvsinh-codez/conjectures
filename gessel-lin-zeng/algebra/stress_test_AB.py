import math
import random
import numpy as np
import mpmath as mp

mp.mp.dps = 50

def apply_D_f(j, p):
    d = len(p) - 1
    out = np.zeros(d + 2, dtype=float)
    for m, cm in enumerate(p):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return out

def poly_from_roots_f(rs):
    p = np.array([1.0])
    for r in rs: p = np.polymul(p, np.array([1.0, -r]))
    return p[::-1]

def wronsk_f(f, g):
    fh, gh = f[::-1], g[::-1]
    return np.polysub(np.polymul(fh, np.polyder(gh)), np.polymul(np.polyder(fh), gh))

# Safe grid in [-5, 5] + asymptotic check (no overflow!)
t_grid = np.concatenate([np.linspace(-5.0, -1.0, 100), np.linspace(-1.0, 0.0, 401), np.linspace(0.0, 5.0, 100)])

def quad_rayleigh_min(X, M, Y):
    """Returns (min W(X,M), min W(M,Y), min (4 W(X,M)W(M,Y) - W(X,Y)^2), lc_disc)."""
    w1, w2, w3 = wronsk_f(X, M), wronsk_f(M, Y), wronsk_f(X, Y)
    disc = np.polysub(4.0 * np.polymul(w1, w2), np.polymul(w3, w3))
    # Normalize disc by max abs coeff to avoid overflow!
    m1 = float(np.min(np.polyval(w1 / max(1e-300, np.max(np.abs(w1))), t_grid)))
    m2 = float(np.min(np.polyval(w2 / max(1e-300, np.max(np.abs(w2))), t_grid)))
    md = float(np.min(np.polyval(disc / max(1e-300, np.max(np.abs(disc))), t_grid)))
    # Also check real roots of disc!
    r = np.roots(disc)
    has_real = any(abs(z.imag) < 1e-6 for z in r)
    return m1, m2, md, float(disc[0]), has_real

# === TEST 1: Heavy adversarial test of Sub-step B: (p0, D_j p1, D_j^2 p2) stable => (p0, D_{j+1} p1, D_{j+1}^2 p2) stable ===
print("=== TEST 1: Heavy stress test of Sub-step B (D_j -> D_{j+1} on random stable Quads) ===")
random.seed(111)
tested_B = 0
fails_B = 0
ex_B = None
# Generate stable triples smartly: start from a random p2 (deg d), pick p1 slightly to left of D_j p2 (and p2 << p1),
# pick p0 slightly to left of D_j p1 (and p1 << p0)!
for trial in range(5000):
    if tested_B >= 400: break
    d = random.choice([0, 1, 2, 3])
    r2 = sorted([-random.uniform(0.05, 0.95) for _ in range(d)])
    p2 = poly_from_roots_f(r2)
    j = random.choice([2, 3, 5, 10, 50])
    Dp2 = apply_D_f(j, p2)
    rDp2 = np.sort(np.roots(Dp2[::-1]).real)
    # Pick r1[m] in (alpha_m, rDp2[m]) so p2 << p1 << D_j p2!
    edges = [-1.0] + r2 + [0.0]
    r1 = [random.uniform(edges[m] + 1e-3*(rDp2[m]-edges[m]), rDp2[m] - 1e-3*(rDp2[m]-edges[m])) for m in range(d + 1)]
    p1 = poly_from_roots_f(r1)
    Dp1 = apply_D_f(j, p1)
    rDp1 = np.sort(np.roots(Dp1[::-1]).real)
    edges1 = [-1.0] + r1 + [0.0]
    r0 = [random.uniform(edges1[m] + 1e-3*(rDp1[m]-edges1[m]), rDp1[m] - 1e-3*(rDp1[m]-edges1[m])) for m in range(d + 2)]
    p0 = poly_from_roots_f(r0)
    # Random positive scales
    p1 *= 10 ** random.uniform(-1, 1)
    p2 *= 10 ** random.uniform(-1, 1)

    X, M, Y = p0, apply_D_f(j, p1), apply_D_f(j, apply_D_f(j, p2))
    m1, m2, md, lc_d, has_r = quad_rayleigh_min(X, M, Y)
    if m1 <= 1e-7 or m2 <= 1e-7 or md <= 1e-7 or lc_d <= 0 or has_r:
        continue
    tested_B += 1
    # Shift j -> j + 1 (and also test j -> j + 10!)
    M2, Y2 = apply_D_f(j + 1, p1), apply_D_f(j + 1, apply_D_f(j + 1, p2))
    m1b, m2b, mdb, lc_db, has_rb = quad_rayleigh_min(X, M2, Y2)
    if m1b <= -1e-6 or m2b <= -1e-6 or mdb <= -1e-6 or lc_db <= 0 or has_rb:
        fails_B += 1
        if ex_B is None:
            ex_B = (d, j, r0, r1, r2, m1b, m2b, mdb, lc_db, has_rb)
print(f"  Sub-step B (D_j -> D_{{j+1}}): fails = {fails_B}/{tested_B}")
if ex_B:
    print(f"    COUNTEREXAMPLE B: d={ex_B[0]}, j={ex_B[1]}, m1b={ex_B[5]:.3e}, m2b={ex_B[6]:.3e}, mdb={ex_B[7]:.3e}, lc={ex_B[8]:.3e}, has_r={ex_B[9]}")

# === TEST 2: Heavy adversarial test of Sub-step A: Q1 stable & Q2 stable => Q1 + w * D_j Q2 stable ===
print("\n=== TEST 2: Heavy stress test of Sub-step A (Q1, Q2 stable => Q1 + w*D_j Q2 stable) ===")
random.seed(222)
tested_A = 0
fails_A = 0
fails_A_w1 = 0
ex_A = None
for trial in range(15000):
    if tested_A >= 300: break
    d = random.choice([0, 1, 2])
    r3 = sorted([-random.uniform(0.05, 0.95) for _ in range(d)])
    p3 = poly_from_roots_f(r3)
    j = random.choice([2, 3, 5, 10])
    # Build p2 << D_j p3, p1 << D_j p2, p0 << D_j p1
    def step_left(p, rs):
        Dp = apply_D_f(j, p)
        rDp = np.sort(np.roots(Dp[::-1]).real)
        edges = [-1.0] + list(rs) + [0.0]
        rn = [random.uniform(edges[m] + 0.05*(rDp[m]-edges[m]), rDp[m] - 0.05*(rDp[m]-edges[m])) for m in range(len(rDp))]
        return poly_from_roots_f(rn), rn
    p2, r2 = step_left(p3, r3)
    p1, r1 = step_left(p2, r2)
    p0, r0 = step_left(p1, r1)
    # Random positive scales (including extreme!)
    p1 *= 10 ** random.uniform(-2, 2)
    p2 *= 10 ** random.uniform(-2, 2)
    p3 *= 10 ** random.uniform(-2, 2)
    X1, M1, Y1 = p0, apply_D_f(j, p1), apply_D_f(j, apply_D_f(j, p2))
    X2, M2, Y2 = p1, apply_D_f(j, p2), apply_D_f(j, apply_D_f(j, p3))
    a1, a2, ad, alc, ar = quad_rayleigh_min(X1, M1, Y1)
    if a1 <= 1e-7 or a2 <= 1e-7 or ad <= 1e-7 or alc <= 0 or ar: continue
    b1, b2, bd, blc, br = quad_rayleigh_min(X2, M2, Y2)
    if b1 <= 1e-7 or b2 <= 1e-7 or bd <= 1e-7 or blc <= 0 or br: continue
    tested_A += 1
    # Test w = 1 and random w > 0!
    for w in [1.0, 10 ** random.uniform(-2, 2)]:
        Xn = X1 + w * apply_D_f(j, X2)
        Mn = M1 + w * apply_D_f(j, M2)
        Yn = Y1 + w * apply_D_f(j, Y2)
        c1, c2, cd, clc, cr = quad_rayleigh_min(Xn, Mn, Yn)
        if c1 <= -1e-6 or c2 <= -1e-6 or cd <= -1e-6 or clc <= 0 or cr:
            if w == 1.0:
                fails_A_w1 += 1
            fails_A += 1
            if ex_A is None:
                ex_A = (d, j, w, c1, c2, cd, clc, cr)
            break
print(f"  Sub-step A (Q1, Q2 stable => Q1 + w*D Q2 stable): fails_any_w={fails_A}/{tested_A}, fails_w1={fails_A_w1}/{tested_A}")
if ex_A:
    print(f"    COUNTEREXAMPLE A: d={ex_A[0]}, j={ex_A[1]}, w={ex_A[2]:.3g}, c1={ex_A[3]:.3e}, c2={ex_A[4]:.3e}, cd={ex_A[5]:.3e}, clc={ex_A[6]:.3e}, has_r={ex_A[7]}")

# === TEST 3: 50-digit mpmath check of binom-E(t,-y) at k=15, t=-0.99+1e-7j AND B_10(x,-xy) at x=-0.0005+1e-9j! ===
print("\n=== TEST 3: 50-digit mpmath check of binom-E at k=15, t=-0.99+1e-7j & B_10(x,-xy) at x=-0.0005+1e-9j ===")
def poly_trim(c):
    while len(c) > 1 and c[-1] == 0: c.pop()
    return c
def apply_D_i(j, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return poly_trim(out)
def add_i(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return poly_trim(out)
def eval_mp(p, z):
    acc = mp.mpc(0)
    for c in reversed(p): acc = acc * z + c
    return acc

ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 16):
    for i in range(k + 1):
        t1 = apply_D_i(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D_i(k + 1, add_i(t1, t2))

E15 = []
for i in range(15):
    p = ht[(14, i)][:]
    for _ in range(i): p = apply_D_i(16, p)
    E15.append(p)
for tv in [mp.mpc('-0.99', '1e-7'), mp.mpc('-0.999', '1e-7'), mp.mpc('-0.001', '1e-7')]:
    yc = [math.comb(14, i) * eval_mp(E15[i], tv) for i in range(14, -1, -1)]
    r = mp.polyroots(yc, maxsteps=400, extraprec=40)
    print(f"  50-digit binom-E k=15 at t={complex(tv)}: minIm(y) = {float(min(w.imag for w in r)):+.4e}")

# Check B_10(x, -x*y) at x = -0.0005 + 1e-9j, -0.0001 + 1e-9j, +0.0001 + 1e-9j
def apply_S_i(L, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (L - m) * cm
    return poly_trim(out)
A = {(0, 0): [1]}
for k in range(1, 11):
    for i in range(k + 1):
        t1 = apply_S_i(3*k - 3 - i, A[(k-1, i)]) if i <= k - 1 else [0]
        t2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S_i(3*k - 2 - i, add_i(t1, t2))
for xv in [mp.mpc('-0.0005', '1e-9'), mp.mpc('-0.0002', '1e-9'), mp.mpc('-0.0001', '1e-9'), mp.mpc('0.0001', '1e-9')]:
    yc = [math.comb(10, i) * eval_mp(A[(10, i)], xv) * ((-xv)**i) for i in range(10, -1, -1)]
    r = mp.polyroots(yc, maxsteps=400, extraprec=40)
    print(f"  50-digit B_10(x,-xy) at x={complex(xv)}: maxIm(y) = {float(max(w.imag for w in r)):+.4e}")
