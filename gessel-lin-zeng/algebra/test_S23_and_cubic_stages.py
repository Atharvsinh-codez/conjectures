"""
Targeted 15-second algebraic experiments:
1. Adversarial test of Stage 2 + Stage 3 COMBINED on barely-stable Quads:
   (U0, D_j U1, D_j^2 U2) real-stable => (D_j U0, D_{j+1} D_j U1, D_{j+1}^2 D_j U2) real-stable?
   Compare failure rate & worst drop in min(4 W_XM W_MY / W_XY^2) vs j = 4, 8, 16, 32!
2. Test Cubic Rayleigh addition in R^2 when u, v, w >= 0.85 (and check why the 38/951 failed).
3. Test whether Quad (X, M, Y) survives wild / step-jump b_j = [2]^6 + [10^5]^6 (where A_k broke!).
4. Check first_base combinatorial refinement of JSP_{k,i} at k=5 (since it passed k=3 and k=4!).
"""
import math, random
from collections import defaultdict
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

t_grid = np.concatenate([np.linspace(-10.0, -1.01, 100), np.linspace(-0.999, -0.001, 300), np.linspace(0.001, 10.0, 100)])

def quad_min_ratio(X, M, Y):
    w1, w2, w3 = wronsk_f(X, M), wronsk_f(M, Y), wronsk_f(X, Y)
    v1 = np.polyval(w1, t_grid)
    v2 = np.polyval(w2, t_grid)
    v3 = np.polyval(w3, t_grid)
    if np.min(v1) <= 0 or np.min(v2) <= 0:
        return -1.0, True
    disc = np.polysub(4.0 * np.polymul(w1, w2), np.polymul(w3, w3))
    r = np.roots(disc)
    has_real = any(abs(z.imag) < 1e-6 for z in r)
    # Also leading coeff
    if disc[0] <= 0: has_real = True
    ratio = np.min((4.0 * v1 * v2) / (v3 * v3 + 1e-300))
    return float(ratio), has_real

print("=== 1. ADVERSARIAL TEST OF S3-ONLY vs S2+S3 (WITH D_j APPLIED FIRST!) ===", flush=True)
random.seed(111)
for j in [3, 6, 12, 24]:
    tested = 0
    fails_S3 = 0
    fails_S23 = 0
    max_drop_S3 = 0.0
    max_drop_S23 = 0.0
    for trial in range(4000):
        if tested >= 150: break
        d = random.choice([0, 1, 2])
        r2 = sorted([-random.uniform(0.05, 0.95) for _ in range(d)])
        p2 = poly_from_roots_f(r2)
        Dp2 = apply_D_f(j, p2)
        rDp2 = np.sort(np.roots(Dp2[::-1]).real)
        edges = [-1.0] + r2 + [0.0]
        r1 = [random.uniform(edges[m] + 1e-3*(rDp2[m]-edges[m]), rDp2[m] - 1e-3*(rDp2[m]-edges[m])) for m in range(d + 1)]
        p1 = poly_from_roots_f(r1)
        Dp1 = apply_D_f(j, p1)
        rDp1 = np.sort(np.roots(Dp1[::-1]).real)
        edges1 = [-1.0] + r1 + [0.0]
        r0 = [random.uniform(edges1[m] + 1e-3*(rDp1[m]-edges1[m]), rDp1[m] - 1e-3*(rDp1[m]-edges1[m])) for m in range(d + 2)]
        p0 = poly_from_roots_f(r0)
        p1 *= 10 ** random.uniform(-0.8, 0.8)
        p2 *= 10 ** random.uniform(-0.8, 0.8)
        X, M, Y = p0, apply_D_f(j, p1), apply_D_f(j, apply_D_f(j, p2))
        r_in, hr_in = quad_min_ratio(X, M, Y)
        if hr_in or r_in <= 1.001: continue
        tested += 1
        # S3 only
        r_s3, hr_s3 = quad_min_ratio(p0, apply_D_f(j+1, p1), apply_D_f(j+1, apply_D_f(j+1, p2)))
        if hr_s3 or r_s3 < 1.0: fails_S3 += 1
        if r_in - r_s3 > max_drop_S3: max_drop_S3 = r_in - r_s3
        # S2 + S3
        h0, h1, h2 = apply_D_f(j, p0), apply_D_f(j, p1), apply_D_f(j, p2)
        r_s23, hr_s23 = quad_min_ratio(h0, apply_D_f(j+1, h1), apply_D_f(j+1, apply_D_f(j+1, h2)))
        if hr_s23 or r_s23 < 1.0: fails_S23 += 1
        if r_in - r_s23 > max_drop_S23: max_drop_S23 = r_in - r_s23
    print(f"  j={j:2d} ({tested:3d} tested): S3_fails={fails_S3:2d} (max_drop={max_drop_S3:.4f}), S2+S3_fails={fails_S23:2d} (max_drop={max_drop_S23:.4f})", flush=True)

print("\n=== 2. INSPECT WHY 38/951 FAILED CUBIC RAYLEIGH ADDITION IN R^2 ===", flush=True)
random.seed(123)
def W2(u, v): return u[0]*v[1] - u[1]*v[0]
def check_rayleigh_cubic(v0, v1, v2, v3):
    c0, c1, c2, c3, c4 = W2(v0,v1), -2*W2(v0,v2), 3*W2(v1,v2)+W2(v0,v3), -2*W2(v1,v3), W2(v2,v3)
    if c0 <= 0 or c4 <= 0: return False
    crts = np.roots([4*c4, 3*c3, 2*c2, c1])
    for z in crts:
        if abs(z.imag) < 1e-7:
            y = z.real
            if (((c4*y + c3)*y + c2)*y + c1)*y + c0 < -1e-9 * max(c0, c4): return False
    return True

for label, u_min in [("all", 0.0), ("u,v,w >= 0.7", 0.7), ("u,v,w >= 0.9", 0.9), ("x1,x2 <= 1", -1.0)]:
    tested = 0; fails = 0
    for _ in range(30000):
        if tested >= 1000: break
        ths = sorted([random.uniform(0.01, math.pi - 0.01) for _ in range(5)])
        rads = [10 ** random.uniform(-1, 1) for _ in range(5)]
        vecs = [(r*math.cos(th), r*math.sin(th)) for r, th in zip(rads, ths)]
        if not (check_rayleigh_cubic(*vecs[:4]) and check_rayleigh_cubic(*vecs[1:])): continue
        u_v = W2(vecs[0], vecs[2]) / (2*math.sqrt(W2(vecs[0],vecs[1])*W2(vecs[1],vecs[2])))
        v_v = W2(vecs[1], vecs[3]) / (2*math.sqrt(W2(vecs[1],vecs[2])*W2(vecs[2],vecs[3])))
        w_v = W2(vecs[2], vecs[4]) / (2*math.sqrt(W2(vecs[2],vecs[3])*W2(vecs[3],vecs[4])))
        x1 = math.sqrt(W2(vecs[0],vecs[1])*W2(vecs[2],vecs[3])) / W2(vecs[1],vecs[2])
        x2 = math.sqrt(W2(vecs[1],vecs[2])*W2(vecs[3],vecs[4])) / W2(vecs[2],vecs[3])
        if u_min > 0 and min(u_v, v_v, w_v) < u_min: continue
        if u_min < 0 and max(x1, x2) > 1.0: continue
        tested += 1
        u = [(vecs[m][0]+vecs[m+1][0], vecs[m][1]+vecs[m+1][1]) for m in range(4)]
        if not check_rayleigh_cubic(*u): fails += 1
    print(f"  Cubic R^2 addition ({label:12s}): fails = {fails}/{tested}", flush=True)

print("\n=== 3. DOES QUAD SURVIVE STEP-JUMP b_j = [2]^5 + [10^4]^5 (WHERE A_k BROKE)? ===", flush=True)
def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p
def poly_add(a, b):
    n = max(len(a), len(b))
    return poly_trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def poly_deriv(a):
    return [i * a[i] for i in range(1, len(a))] if len(a) > 1 else [0]
def apply_Dj(p, j):
    dp = poly_deriv(p)
    res = [0] * (len(p) + 2)
    for i, c in enumerate(p):
        res[i] += c; res[i+1] += j * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c; res[i+2] += c
    return poly_trim(res)
def poly_eval_mp(p, x):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

for name, bs in [("step-jump [2]*5 + [10000]*5", [2]*5 + [10000]*5),
                 ("step-jump [10000]*5 + [2]*5", [10000]*5 + [2]*5),
                 ("alternating [2, 10000]*5", [2, 10000]*5)]:
    polys = {0: [1, 2], 1: [1]} # k=1
    for step, b in enumerate(bs, start=2):
        new_p = {}
        for i in range(step + 1):
            p1 = apply_Dj(polys[i], b) if i <= step - 1 else [0]
            p2 = polys[i-1] if i >= 1 else [0]
            new_p[i] = apply_Dj(poly_add(p1, p2), b)
        polys = new_p
    # Check Quad (polys[i-2], D_b polys[i-1], D_b^2 polys[i]) with b_next = bs[-1] and b_next = 10000!
    k = len(bs) + 1
    b_next = bs[-1]
    min_r = mp.mpf('1e100')
    w_neg = 0
    grid = [mp.mpf(v) for v in [-10, -2, -1.1, -0.9, -0.5, -0.1, -0.01, -1e-4, 0, 1e-4, 0.01, 0.1, 1, 10]]
    for i in range(2, k + 1):
        X = polys[i-2]
        M = apply_Dj(polys[i-1], b_next)
        Y = apply_Dj(apply_Dj(polys[i], b_next), b_next)
        dX, dM, dY = poly_deriv(X), poly_deriv(M), poly_deriv(Y)
        for tv in grid:
            x, dx = poly_eval_mp(X, tv), poly_eval_mp(dX, tv)
            m, dm = poly_eval_mp(M, tv), poly_eval_mp(dM, tv)
            y, dy = poly_eval_mp(Y, tv), poly_eval_mp(dY, tv)
            wXM = x*dm - dx*m; wMY = m*dy - dm*y; wXY = x*dy - dx*y
            if wXM <= 0 or wMY <= 0: w_neg += 1
            if wXY != 0:
                r = (4*wXM*wMY)/(wXY*wXY)
                if r < min_r: min_r = r
    print(f"  {name:30s}: w_neg={w_neg}, min 4*W_XM*W_MY/W_XY^2 = {mp.nstr(min_r, 6)}", flush=True)
