import numpy as np
import mpmath as mp
import sympy as sp
import math

# 1. TEST MODIFIED TURAN M_i(t) := (-t)(1+t) W(P_{i+1}, P_i)(t) + P_i(t)^2 - P_{i-1}(t) P_{i+1}(t)
#    (where W(f, g) = f g' - f' g, so W(P_{i+1}, P_i) > 0 on R when P_{i+1} << P_i!)
#    Recall: M_i(t) > 0 on (-1, 0) ==> D_b P_{i+1} + P_i << D_b P_i + P_{i-1} FOR ALL b > 1 (2-line proof!).
#    Test:
#      (a) Does D_b preserve M_i(t) > 0 on (-1, 0) for arbitrary triples (P_{i-1}, P_i, P_{i+1})?
#      (b) Does D_b + y (or D_b(D_b + y)) preserve M_i(t) > 0 on (-1, 0) for arbitrary (P_{i-2}, P_{i-1}, P_i, P_{i+1})?
#      (c) What if we ALSO assume J_raw >= gamma (e.g. gamma = 0.88 or 0.5) AND M_i > 0?

def apply_D_np(b, p):
    res = np.zeros(len(p)+1)
    for n, c in enumerate(p):
        res[n] += c*(n+1); res[n+1] += c*(n+b)
    return res

def add_np(a, b):
    n = max(len(a), len(b))
    r = np.zeros(n); r[:len(a)] += a; r[:len(b)] += b
    return r

def wronskian_np(f, g):
    fp = np.arange(1, len(f))*f[1:] if len(f) > 1 else np.array([0.0])
    gp = np.arange(1, len(g))*g[1:] if len(g) > 1 else np.array([0.0])
    a = np.convolve(f, gp); b = np.convolve(fp, g)
    n = max(len(a), len(b))
    res = np.zeros(n); res[:len(a)] += a; res[:len(b)] -= b
    return res

def eval_np(p, t):
    v = np.zeros_like(t)
    for c in p[::-1]: v = v*t + c
    return v

def from_roots(rs, scale=1.0):
    p = np.array([scale], dtype=np.float64)
    for r in rs: p = np.convolve(p, [-r, 1.0])
    return p

ts = np.concatenate([
    -1.0 + np.geomspace(1e-5, 0.49, 80),
    -np.geomspace(1e-5, 0.5, 80),
])

def eval_M(p0, p1, p2, ts):
    # deg p0 = d+1, p1 = d, p2 = d-1 (p2 << p1 << p0)
    w21 = eval_np(wronskian_np(p2, p1), ts)
    v0, v1, v2 = eval_np(p0, ts), eval_np(p1, ts), eval_np(p2, ts)
    return (-ts)*(1.0 + ts)*w21 + v1*v1 - v0*v2

def min_J_np(p2, p1, p0, ts):
    w1 = eval_np(wronskian_np(p2, p1), ts)
    w2 = eval_np(wronskian_np(p1, p0), ts)
    w12 = eval_np(wronskian_np(p2, p0), ts)
    if np.min(w1) <= 0 or np.min(w2) <= 0: return -1.0
    return float(np.min(4.0*w1*w2 / (w12*w12 + 1e-300)))

print("=== 1. DOES D_b PRESERVE M_i(t) > 0 ON (-1, 0)? ===")
rng = np.random.default_rng(42)
for b in [2.0, 10.0, 100.0, 10000.0]:
    valid = 0; fail_D = 0; valid_J88 = 0; fail_D_J88 = 0
    for _ in range(3000):
        r0 = np.sort(-rng.uniform(0.01, 0.99, size=3))
        r1 = np.array([rng.uniform(r0[0], r0[1]), rng.uniform(r0[1], r0[2])])
        r2 = np.array([rng.uniform(r1[0], r1[1])])
        s1 = np.exp(rng.uniform(-3, 3)); s2 = np.exp(rng.uniform(-3, 3))
        p0, p1, p2 = from_roots(r0, 1.0), from_roots(r1, s1), from_roots(r2, s2)
        if np.min(eval_M(p0, p1, p2, ts)) <= 0: continue
        valid += 1
        q0, q1, q2 = apply_D_np(b, p0), apply_D_np(b, p1), apply_D_np(b, p2)
        ok_q = np.min(eval_M(q0, q1, q2, ts)) > 0
        if not ok_q: fail_D += 1
        if min_J_np(p2, p1, p0, ts) >= 0.88:
            valid_J88 += 1
            if not ok_q: fail_D_J88 += 1
    print(f"  D_b (b={b:7.1f}): valid={valid:4d}, fail_M={fail_D:4d} | with J>=0.88: valid={valid_J88:4d}, fail_M={fail_D_J88:4d}")

print("\n=== 2. DOES (D_b + y) OR D_b(D_b + y) PRESERVE M_i(t) > 0 ON (-1, 0)? ===")
for b in [2.0, 10.0, 100.0]:
    valid = 0; fail_sgl = 0; fail_dbl = 0
    valid_J88 = 0; fail_sgl_J88 = 0; fail_dbl_J88 = 0
    for _ in range(4000):
        r0 = np.sort(-rng.uniform(0.01, 0.99, size=4))
        r1 = np.array([rng.uniform(r0[m], r0[m+1]) for m in range(3)])
        r2 = np.array([rng.uniform(r1[m], r1[m+1]) for m in range(2)])
        r3 = np.array([rng.uniform(r2[0], r2[1])])
        s1, s2, s3 = np.exp(rng.uniform(-2, 2, size=3))
        h0, h1, h2, h3 = from_roots(r0, 1.0), from_roots(r1, s1), from_roots(r2, s2), from_roots(r3, s3)
        if np.min(eval_M(h0, h1, h2, ts)) <= 0 or np.min(eval_M(h1, h2, h3, ts)) <= 0: continue
        valid += 1
        g0 = add_np(apply_D_np(b, h1), h0)
        g1 = add_np(apply_D_np(b, h2), h1)
        g2 = add_np(apply_D_np(b, h3), h2)
        ok_s = np.min(eval_M(g0, g1, g2, ts)) > 0
        if not ok_s: fail_sgl += 1
        d0, d1, d2 = apply_D_np(b, g0), apply_D_np(b, g1), apply_D_np(b, g2)
        ok_d = np.min(eval_M(d0, d1, d2, ts)) > 0
        if not ok_d: fail_dbl += 1
        if min_J_np(h2, h1, h0, ts) >= 0.88 and min_J_np(h3, h2, h1, ts) >= 0.88:
            valid_J88 += 1
            if not ok_s: fail_sgl_J88 += 1
            if not ok_d: fail_dbl_J88 += 1
    print(f"  b={b:5.1f}: valid={valid:4d} (sgl_fail={fail_sgl:4d}, dbl_fail={fail_dbl:4d}) | with J>=0.88: valid={valid_J88:4d} (sgl_fail={fail_sgl_J88:4d}, dbl_fail={fail_dbl_J88:4d})")

print("\n=== 3. MULTI-STEP ADVERSARIAL ATTACK ON b_j >= 2 STARTING FROM [2]*60 (10 MORE WILD STEPS!) ===")
mp.mp.dps = 120
def apply_D_mp(b, p):
    res = [mp.mpf(0)]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+b)
    return res

def add_mp(a, b):
    res = [mp.mpf(0)]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v

def deriv_mp(p):
    return [r*p[r] for r in range(1, len(p))] if len(p) > 1 else [mp.mpf(0)]

def eval_W_mp(f, g, t):
    return eval_mp(f, t)*eval_mp(deriv_mp(g), t) - eval_mp(deriv_mp(f), t)*eval_mp(g, t)

# Start at K=60 with b = [2]*60, single-operator, keep indices 15..45, apply 6 targeted huge steps!
K = 60
D2 = [[mp.mpf(1)]]
for m in range(1, K + 2):
    D2.append(apply_D_mp(mp.mpf(2), D2[-1]))
P = {r: [mp.mpf(math.comb(K, r)) * c for c in D2[K - r]] for r in range(18, 43)}

# At each of 6 steps, find the slope scale at i=30 and hit with b = 4 * slope_scale (or 2, or 100 * slope_scale)
# Check BOTH min rel Turan AND min M_i AND min J_raw!
for step in range(1, 7):
    i0 = 30
    s_sc = eval_mp(deriv_mp(P[i0]), mp.mpf(0)) / eval_mp(P[i0], mp.mpf(0))
    b_cur = mp.mpf(4) * s_sc
    newP = {}
    for r in range(18 + step, 43):
        newP[r] = add_mp(apply_D_mp(b_cur, P[r]), P[r-1])
    P = newP
    # Check at t near -3/b_cur, -1/b_cur, -0.5/b_cur, 0, -0.5, -0.9999
    t_list = [mp.mpf(0), -mp.mpf('0.2')/b_cur, -mp.mpf('1')/b_cur, -mp.mpf('3')/b_cur, -mp.mpf('10')/b_cur, mp.mpf('-0.5'), mp.mpf('-0.9999')]
    mn_tur = mp.mpf('1e100'); mn_M = mp.mpf('1e100'); mn_J = mp.mpf('1e100')
    for tv in t_list:
        v0, v1, v2 = eval_mp(P[i0-1], tv), eval_mp(P[i0], tv), eval_mp(P[i0+1], tv)
        w21 = eval_W_mp(P[i0+1], P[i0], tv)
        w10 = eval_W_mp(P[i0], P[i0-1], tv)
        w20 = eval_W_mp(P[i0+1], P[i0-1], tv)
        tur = (v1*v1 - v0*v2) / (v1*v1 + abs(v0*v2) + mp.mpf('1e-300'))
        M = ((-tv)*(1+tv)*w21 + v1*v1 - v0*v2) / (abs((-tv)*(1+tv)*w21) + v1*v1 + abs(v0*v2) + mp.mpf('1e-300'))
        J = 4*w21*w10/(w20*w20) if (w21 > 0 and w10 > 0 and w20 != 0) else mp.mpf(-1)
        if tur < mn_tur: mn_tur = tur
        if M < mn_M: mn_M = M
        if J < mn_J: mn_J = J
    print(f"  Step K+{step} (b={float(b_cur):.2e}): min_rel_Tur={float(mn_tur):+.5f}, min_rel_M={float(mn_M):+.5f}, min_Jraw={float(mn_J):+.5f}")
