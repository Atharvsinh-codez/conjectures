import math
import random
import numpy as np

def poly_trim(c):
    while len(c) > 1 and c[-1] == 0: c.pop()
    return c

def apply_D(j, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
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

def mul_poly(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b): out[i + j] += ai * bj
    return poly_trim(out)

def deriv(p):
    return poly_trim([m * p[m] for m in range(1, len(p))]) if len(p) > 1 else [0]

def wronskian(f, g):
    return sub_poly(mul_poly(f, deriv(g)), mul_poly(deriv(f), g))

def eval_c(p, z):
    acc = 0j
    for c in reversed(p): acc = acc * z + c
    return acc

def eval_np(p, x):
    acc = np.zeros_like(x, dtype=float)
    for c in reversed(p): acc = acc * x + float(c)
    return acc

ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 16):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

# === TEST 1: Fast check of 4*W(X,M)*W(M,Y) - W(X,Y)^2 > 0 on all R for k=3..15 (and coeffs in t & x!) ===
print("=== TEST 1: Fast check of 4*W(X,M)*W(M,Y) - W(X,Y)^2 for Quad X - 2My + Yy^2 (k=3..15) ===")
t_R = np.concatenate([-10.0**np.linspace(3, -5, 100), [0.0], 10.0**np.linspace(-5, 3, 100), np.linspace(-1.5, 0.5, 400)])
for k in range(3, 16):
    min_val = 1e300
    min_t_c = 10**99
    bad = None
    for i in range(2, k):
        X = ht[(k-1, i-2)]
        M = apply_D(k+1, ht[(k-1, i-1)])
        Y = apply_D(k+1, apply_D(k+1, ht[(k-1, i)]))
        wxm = wronskian(X, M)
        wmy = wronskian(M, Y)
        wxy = wronskian(X, Y)
        disc = sub_poly([4 * v for v in mul_poly(wxm, wmy)], mul_poly(wxy, wxy))
        min_t_c = min(min_t_c, min(disc))
        vals = eval_np(disc, t_R)
        mv = float(np.min(vals))
        if mv <= 0:
            bad = (i, mv)
        if mv < min_val: min_val = mv
    print(f"  k={k:2d}: 4W(X,M)W(M,Y)-W(X,Y)^2 > 0 on R: {'PASS!' if bad is None else f'FAIL {bad}'} | all t-coeffs >= 0: {min_t_c >= 0}")

# === TEST 2: Fine stress test of binom-E(t, -y) at k=10..15 near t=0, t=-1, (-1,0) ===
print("\n=== TEST 2: Fine stress test of binom-E(t, -y) at k=10..15 ===")
fine_ts = [complex(re, im) for re in [-1.01, -1.001, -0.999, -0.99, -0.5, -0.1, -0.01, -0.001, -1e-4, 0.0, 1e-4, 0.001, 0.01, 0.1] for im in [1e-7, 1e-4]]
for k in range(10, 16):
    E = []
    for i in range(k):
        p = ht[(k-1, i)][:]
        for _ in range(i): p = apply_D(k + 1, p)
        E.append(p)
    min_im = 1e9
    worst = None
    for tv in fine_ts:
        yc = [math.comb(k-1, i) * eval_c(E[i], tv) for i in range(k-1, -1, -1)]
        yc = [c / max(abs(z) for z in yc) for c in yc]
        mi = min(r.imag for r in np.roots(yc))
        if mi < min_im: min_im = mi; worst = tv
    print(f"  k={k:2d}: binom-E(t,-y) fine minIm(y) = {min_im:+.4e} at t={worst}")

# === TEST 3: Test Sub-step B (shifting D_j -> D_{j+1} on ANY real-stable Quad p_0 - 2(D_j p_1)y + (D_j^2 p_2)y^2)! ===
print("\n=== TEST 3: Does shifting D_j -> D_{j+1} ALWAYS preserve real stability of p_0 - 2(D_j p_1)y + (D_j^2 p_2)y^2?! ===")
def poly_from_roots_f(rs):
    p = np.array([1.0])
    for r in rs: p = np.polymul(p, np.array([1.0, -r]))
    return p[::-1]
def apply_D_f(j, p):
    d = len(p) - 1
    out = np.zeros(d + 2, dtype=float)
    for m, cm in enumerate(p):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return out
def wronsk_f(f, g):
    # f, g low->high -> convert to high->low
    fh, gh = f[::-1], g[::-1]
    return np.polysub(np.polymul(fh, np.polyder(gh)), np.polymul(np.polyder(fh), gh))

random.seed(404)
tested_B = 0
fails_B = 0
for trial in range(500):
    d = random.choice([0, 1, 2])
    r0 = sorted([-random.uniform(0.02, 0.98) for _ in range(d + 2)])
    r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 1)]
    r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d)]
    p0, p1, p2 = poly_from_roots_f(r0), poly_from_roots_f(r1), poly_from_roots_f(r2)
    j = random.choice([2, 3, 5, 10])
    X = p0
    M = apply_D_f(j, p1)
    Y = apply_D_f(j, apply_D_f(j, p2))
    # Check if X << M << Y and 4 W(X,M) W(M,Y) - W(X,Y)^2 > 0 on t_R!
    wxm, wmy, wxy = wronsk_f(X, M), wronsk_f(M, Y), wronsk_f(X, Y)
    if np.min(np.polyval(wxm, t_R)) <= 0 or np.min(np.polyval(wmy, t_R)) <= 0:
        continue
    disc = np.polysub(4.0 * np.polymul(wxm, wmy), np.polymul(wxy, wxy))
    if np.min(np.polyval(disc, t_R)) <= 1e-6:
        continue
    tested_B += 1
    # Now shift j -> j + 1!
    M2 = apply_D_f(j + 1, p1)
    Y2 = apply_D_f(j + 1, apply_D_f(j + 1, p2))
    wxm2, wmy2, wxy2 = wronsk_f(X, M2), wronsk_f(M2, Y2), wronsk_f(X, Y2)
    disc2 = np.polysub(4.0 * np.polymul(wxm2, wmy2), np.polymul(wxy2, wxy2))
    if np.min(np.polyval(wxm2, t_R)) <= -1e-6 or np.min(np.polyval(wmy2, t_R)) <= -1e-6 or np.min(np.polyval(disc2, t_R)) <= -1e-6:
        fails_B += 1
print(f"  Sub-step B (j -> j+1 on random stable Quads): fails = {fails_B}/{tested_B}")

# === TEST 4: Test Sub-step A (if Q1 = p0 - 2(D p1)y + (D^2 p2)y^2 AND Q2 = p1 - 2(D p2)y + (D^2 p3)y^2 are BOTH real-stable,
#             is Q1 + D Q2 ALWAYS real-stable?!) ===
print("\n=== TEST 4: Does Q1 stable & Q2 stable => Q1 + D Q2 ALWAYS real-stable?! ===")
tested_A = 0
fails_A = 0
for trial in range(1500):
    d = random.choice([0, 1])
    r0 = sorted([-random.uniform(0.02, 0.98) for _ in range(d + 3)])
    r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 2)]
    r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d + 1)]
    r3 = [random.uniform(r2[m]+1e-4, r2[m+1]-1e-4) for m in range(d)]
    p0, p1, p2, p3 = poly_from_roots_f(r0), poly_from_roots_f(r1), poly_from_roots_f(r2), poly_from_roots_f(r3)
    # Random positive scales
    p1 *= random.uniform(0.5, 2.0)
    p2 *= random.uniform(0.5, 2.0)
    p3 *= random.uniform(0.5, 2.0)
    j = random.choice([2, 3, 5, 10])
    # Q1: (X1, M1, Y1) = (p0, D p1, D^2 p2)
    X1, M1, Y1 = p0, apply_D_f(j, p1), apply_D_f(j, apply_D_f(j, p2))
    # Q2: (X2, M2, Y2) = (p1, D p2, D^2 p3)
    X2, M2, Y2 = p1, apply_D_f(j, p2), apply_D_f(j, apply_D_f(j, p3))
    def is_quad_stable(X, M, Y):
        w1, w2, w3 = wronsk_f(X, M), wronsk_f(M, Y), wronsk_f(X, Y)
        if np.min(np.polyval(w1, t_R)) <= 1e-7 or np.min(np.polyval(w2, t_R)) <= 1e-7: return False
        disc = np.polysub(4.0 * np.polymul(w1, w2), np.polymul(w3, w3))
        return np.min(np.polyval(disc, t_R)) > 1e-7
    if not (is_quad_stable(X1, M1, Y1) and is_quad_stable(X2, M2, Y2)):
        continue
    tested_A += 1
    # Check Q1 + D Q2: (X1 + D X2, M1 + D M2, Y1 + D Y2) = (p0 + D p1, D p1 + D^2 p2, D^2 p2 + D^3 p3)!
    Xn = X1 + apply_D_f(j, X2)
    Mn = M1 + apply_D_f(j, M2)
    Yn = Y1 + apply_D_f(j, Y2)
    if not is_quad_stable(Xn, Mn, Yn):
        fails_A += 1
print(f"  Sub-step A (Q1, Q2 stable => Q1 + D Q2 stable): fails = {fails_A}/{tested_A}")
