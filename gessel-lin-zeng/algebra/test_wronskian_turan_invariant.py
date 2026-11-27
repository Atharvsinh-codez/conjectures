"""
TEST THE EXACT ALGEBRAIC INVARIANT DISCOVERED FROM EXPANDING (T p_i)^2 - (T p_{i-1})(T p_{i+1}) IN Phi:
  (T p_i)^2 - (T p_{i-1})(T p_{i+1}) >= 0 for ALL Phi = a + b*t in R (at fixed t in (-1,0))
  <==> p_i^2 >= p_{i-1} p_{i+1}  AND  4 W(p_{i-1}, p_i) W(p_i, p_{i+1}) >= W(p_{i-1}, p_{i+1})^2  at t!

Tests:
1. Check J_i(t) := 4 W(p_{i-1}, p_i) W(p_i, p_{i+1}) - W(p_{i-1}, p_{i+1})^2 on (-1, 0) for:
   (a) GLZ htilde_{k,i} for k=3..20 and k=35!
   (b) Universal prod_{j=1}^k (T_j + w_j y)(1) (50 wild trials, k=6)!
   (c) Universal Laguerre prod_{j=1}^k (L_j + w_j y)(1) on (-inf, 0) (50 wild trials, k=6)!
2. Check C_0, C_1, C_2 and the 3x3 matrix of q_i^2 - q_{i-1} q_{i+1} in (1, Phi, w) on (-1,0)!
3. Check if (Turan + J_i >= 0 on (-1,0)) PROPAGATES in 1 step:
   (a) under pure T: p_i -> T p_i! (WAIT: does J_i(T p) >= 0 on (-1,0) follow from J_i(p) >= 0 on (-1,0)?!)
   (b) under q_i = T p_i + w p_{i-1} on random 4-tuples satisfying Turan + J >= 0 on (-1,0)!
   (c) ALSO check the Plucker ratio W(p_{i-2}, p_{i-1}) W(p_i, p_{i+1}) / W(p_{i-1}, p_i)^2 on (-1,0) for universal T_j & Laguerre!
"""
import math, random
import numpy as np
import mpmath as mp
mp.mp.dps = 50

def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p
def poly_add(a, b):
    n = max(len(a), len(b))
    return poly_trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def poly_sub(a, b):
    n = max(len(a), len(b))
    return poly_trim([(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(n)])
def poly_mul(a, b):
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j, bj in enumerate(b):
            res[i+j] += ai * bj
    return poly_trim(res)
def poly_deriv(a):
    return [i * a[i] for i in range(1, len(a))] if len(a) > 1 else [0]
def wronskian(f, g):
    return poly_sub(poly_mul(f, poly_deriv(g)), poly_mul(poly_deriv(f), g))
def apply_T(p, a, b):
    # t(1+t)p' + (a + b t) p
    dp = poly_deriv(p)
    res = [0] * (len(p) + 2)
    for i, c in enumerate(p):
        res[i] += a * c
        res[i+1] += b * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c
        res[i+2] += c
    return poly_trim(res)
def apply_L(p, a, b):
    res = [0] * (len(p) + 1)
    for i, c in enumerate(p):
        res[i] += (i + a) * c
        res[i+1] += b * c
    return poly_trim(res)
def poly_eval_mp(p, x):
    acc = mp.mpf(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc

def compute_htilde(K):
    h = {1: {0: [1, 2], 1: [1]}}
    for k in range(2, K + 1):
        h[k] = {}
        for i in range(k + 1):
            p1 = apply_T(h[k-1][i], 1, k + 1) if i <= k - 1 else [0]
            p2 = h[k-1][i-1] if i >= 1 else [0]
            h[k][i] = apply_T(poly_add(p1, p2), 1, k + 1)
    return h

print("=== TEST 1a: J_i(t) = 4 W(p_{i-1},p_i) W(p_i,p_{i+1}) - W(p_{i-1},p_{i+1})^2 on (-1,0) for GLZ ===", flush=True)
ht = compute_htilde(35)
tgrid = [-mp.mpf(m)/50 for m in range(1, 50)] + [-mp.mpf(10)**(-p) for p in [3, 6, 10, 16]] + [mp.mpf(-1) + mp.mpf(10)**(-p) for p in [3, 6, 10]]
for k in [5, 10, 15, 20, 35]:
    # Precompute derivatives
    dh = {i: poly_deriv(ht[k][i]) for i in range(k+1)}
    min_ratio = mp.mpf('1e100')
    neg_cnt = 0
    worst = None
    for i in range(1, k):
        for tv in tgrid:
            f0, df0 = poly_eval_mp(ht[k][i-1], tv), poly_eval_mp(dh[i-1], tv)
            f1, df1 = poly_eval_mp(ht[k][i], tv),   poly_eval_mp(dh[i], tv)
            f2, df2 = poly_eval_mp(ht[k][i+1], tv), poly_eval_mp(dh[i+1], tv)
            w01 = f0*df1 - df0*f1
            w12 = f1*df2 - df1*f2
            w02 = f0*df2 - df0*f2
            # Wait: check sign of w01*w12! Since h_{i+1} << h_i << h_{i-1}, W(h_{i-1}, h_i) < 0 and W(h_i, h_{i+1}) < 0, so product > 0!
            r = (4 * w01 * w12) / (w02 * w02) if w02 != 0 else mp.mpf('1e100')
            if 4*w01*w12 - w02*w02 < 0:
                neg_cnt += 1
            if r < min_ratio:
                min_ratio = r
                worst = (i, float(tv))
    print(f"  GLZ k={k:2d}: J_i < 0 count on (-1,0) = {neg_cnt}, min 4*W01*W12/W02^2 = {mp.nstr(min_ratio, 8)} at {worst}", flush=True)

print("\n=== TEST 1b: J_i(t) on (-1,0) for Universal prod_{j=1}^k (T_j + w_j y)(1) (50 wild trials, k=6) ===", flush=True)
random.seed(7)
fail_J_T = 0
min_r_T = mp.mpf('1e100')
worst_T = None
tgrid_fast = [-mp.mpf(m)/30 for m in range(1, 30)] + [-mp.mpf(10)**(-p) for p in [3, 6]] + [mp.mpf(-1) + mp.mpf(10)**(-p) for p in [3, 6]]
for tr in range(50):
    k = 6
    polys = {0: [mp.mpf(1)]}
    for j in range(1, k+1):
        a_j = mp.mpf(10 ** random.uniform(-2, 2))
        b_j = a_j + mp.mpf(10 ** random.uniform(-2, 2))
        w_j = mp.mpf(10 ** random.uniform(-2, 2))
        new_p = {}
        for i in range(j+1):
            p1 = apply_T(polys[i], a_j, b_j) if i in polys else [mp.mpf(0)]
            p2 = [w_j * c for c in polys[i-1]] if (i-1) in polys else [mp.mpf(0)]
            new_p[i] = poly_add(p1, p2)
        polys = new_p
    dp = {i: poly_deriv(polys[i]) for i in range(k+1)}
    bad = False
    for i in range(1, k):
        for tv in tgrid_fast:
            f0, df0 = poly_eval_mp(polys[i-1], tv), poly_eval_mp(dp[i-1], tv)
            f1, df1 = poly_eval_mp(polys[i], tv),   poly_eval_mp(dp[i], tv)
            f2, df2 = poly_eval_mp(polys[i+1], tv), poly_eval_mp(dp[i+1], tv)
            w01 = f0*df1 - df0*f1
            w12 = f1*df2 - df1*f2
            w02 = f0*df2 - df0*f2
            if 4*w01*w12 - w02*w02 < 0:
                bad = True
            if w02 != 0:
                r = (4*w01*w12)/(w02*w02)
                if r < min_r_T:
                    min_r_T = r
                    worst_T = (tr, i, float(tv))
    if bad:
        fail_J_T += 1
print(f"  Universal T_j (50 wild trials, k=6): J_i < 0 fails = {fail_J_T}/50, min 4*W01*W12/W02^2 = {mp.nstr(min_r_T, 8)} at {worst_T}", flush=True)

print("\n=== TEST 1c: J_i(x) on (-inf,0) for Universal Laguerre prod (L_j + w_j y)(1) (50 wild trials, k=6) ===", flush=True)
fail_J_L = 0
min_r_L = mp.mpf('1e100')
worst_L = None
xgrid_fast = [-mp.mpf(10)**(mp.mpf(p)/4) for p in range(-16, 17)]
for tr in range(50):
    k = 6
    polys = {0: [mp.mpf(1)]}
    for j in range(1, k+1):
        a_j = mp.mpf(10 ** random.uniform(-2, 2))
        b_j = mp.mpf(10 ** random.uniform(-2, 2))
        w_j = mp.mpf(10 ** random.uniform(-2, 2))
        new_p = {}
        for i in range(j+1):
            p1 = apply_L(polys[i], a_j, b_j) if i in polys else [mp.mpf(0)]
            p2 = [w_j * c for c in polys[i-1]] if (i-1) in polys else [mp.mpf(0)]
            new_p[i] = poly_add(p1, p2)
        polys = new_p
    dp = {i: poly_deriv(polys[i]) for i in range(k+1)}
    bad = False
    for i in range(1, k):
        for xv in xgrid_fast:
            f0, df0 = poly_eval_mp(polys[i-1], xv), poly_eval_mp(dp[i-1], xv)
            f1, df1 = poly_eval_mp(polys[i], xv),   poly_eval_mp(dp[i], xv)
            f2, df2 = poly_eval_mp(polys[i+1], xv), poly_eval_mp(dp[i+1], xv)
            w01 = f0*df1 - df0*f1
            w12 = f1*df2 - df1*f2
            w02 = f0*df2 - df0*f2
            if 4*w01*w12 - w02*w02 < 0:
                bad = True
            if w02 != 0:
                r = (4*w01*w12)/(w02*w02)
                if r < min_r_L:
                    min_r_L = r
                    worst_L = (tr, i, float(xv))
    if bad:
        fail_J_L += 1
print(f"  Universal Laguerre L_j (50 wild trials, k=6): J_i < 0 fails = {fail_J_L}/50, min 4*W01*W12/W02^2 = {mp.nstr(min_r_L, 8)} at {worst_L}", flush=True)

print("\n=== TEST 2: Does J_i >= 0 on (-1,0) propagate under PURE T (p_i -> T p_i)? ===", flush=True)
# Test 200 random triples p2 << p1 << p0 in (-1,0) (deg 2,3,4) satisfying Turan + J >= 0 on (-1,0),
# check if T p2, T p1, T p0 ALSO satisfy J >= 0 on (-1,0) for random T(a,b)!
def poly_from_roots(rts):
    p = [1.0]
    for r in rts:
        p = poly_mul(p, [-r, 1.0])
    return p

random.seed(99)
valid_triples = 0
fail_pureT_J = 0
fail_pureT_turan = 0
ts = np.linspace(-0.999, -0.001, 400)
for trial in range(2000):
    if valid_triples >= 200: break
    # Generate p2 << p1 << p0 in (-1,0) of degrees 2, 3, 4
    # Roots r0_1 < r1_1 < r0_2 < r1_2 < r0_3 < r1_3 < r0_4, and r2_1 in (r1_1, r1_2), r2_2 in (r1_2, r1_3)
    pts = sorted([random.uniform(-0.98, -0.02) for _ in range(7)])
    r0 = [pts[0], pts[2], pts[4], pts[6]]
    r1 = [pts[1], pts[3], pts[5]]
    r2 = [random.uniform(r1[0], r1[1]), random.uniform(r1[1], r1[2])]
    c0, c1, c2 = 10**random.uniform(-1,1), 10**random.uniform(-1,1), 10**random.uniform(-1,1)
    p0 = [c0*x for x in poly_from_roots(r0)]
    p1 = [c1*x for x in poly_from_roots(r1)]
    p2 = [c2*x for x in poly_from_roots(r2)]
    # Check Turan + J >= 0 on (-1, 0)
    v0 = np.polyval(list(reversed(p0)), ts)
    v1 = np.polyval(list(reversed(p1)), ts)
    v2 = np.polyval(list(reversed(p2)), ts)
    if np.min(v1*v1 - v0*v2) < 0: continue
    dv0 = np.polyval(list(reversed(poly_deriv(p0))), ts)
    dv1 = np.polyval(list(reversed(poly_deriv(p1))), ts)
    dv2 = np.polyval(list(reversed(poly_deriv(p2))), ts)
    w01 = v0*dv1 - dv0*v1
    w12 = v1*dv2 - dv1*v2
    w02 = v0*dv2 - dv0*v2
    if np.min(4*w01*w12 - w02*w02) < 0: continue
    valid_triples += 1
    # Apply random T(a,b) with 0 < a < b
    a = 10**random.uniform(-1, 1)
    b = a + 10**random.uniform(-1, 1)
    q0, q1, q2 = apply_T(p0, a, b), apply_T(p1, a, b), apply_T(p2, a, b)
    qv0 = np.polyval(list(reversed(q0)), ts)
    qv1 = np.polyval(list(reversed(q1)), ts)
    qv2 = np.polyval(list(reversed(q2)), ts)
    if np.min(qv1*qv1 - qv0*qv2) < -1e-9 * np.max(qv1*qv1):
        fail_pureT_turan += 1
    dq0 = np.polyval(list(reversed(poly_deriv(q0))), ts)
    dq1 = np.polyval(list(reversed(poly_deriv(q1))), ts)
    dq2 = np.polyval(list(reversed(poly_deriv(q2))), ts)
    qw01 = qv0*dq1 - dq0*qv1
    qw12 = qv1*dq2 - dq1*qv2
    qw02 = qv0*dq2 - dq0*qv2
    if np.min(4*qw01*qw12 - qw02*qw02) < -1e-9 * np.max(np.abs(4*qw01*qw12)):
        fail_pureT_J += 1
print(f"  Random valid triples ({valid_triples} tested): pure-T Turan fails={fail_pureT_turan}, pure-T J_i fails={fail_pureT_J}", flush=True)
