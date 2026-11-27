import numpy as np
import mpmath as mp

# High-precision (100 dps) check of whether double_same_T or double_Dm_w1 REALLY breaks
# Turan or interlacing or RZ (to separate genuine counterexamples from float64 noise near t=-1)!

def apply_T_mp(p, a, b):
    d = len(p) - 1
    res = [mp.mpf(0)] * (d + 2)
    for n in range(d + 1):
        res[n] += (a + n) * p[n]
        res[n + 1] += (b + n) * p[n]
    return res

def add_mp(p, q):
    n = max(len(p), len(q))
    r = [mp.mpf(0)] * n
    for i, v in enumerate(p): r[i] += v
    for i, v in enumerate(q): r[i] += v
    return r

def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p):
        v = v * t + c
    return v

def step_double_mp(P, a, b, w):
    k = len(P)
    newP = []
    for i in range(k + 1):
        t2 = apply_T_mp(P[i], a, b) if i < k else [mp.mpf(0)]
        wy = [w * c for c in P[i - 1]] if i >= 1 else [mp.mpf(0)]
        newP.append(apply_T_mp(add_mp(t2, wy), a, b))
    return newP

def check_rz_exact_sturm_mp(p, dps=100):
    mp.mp.dps = dps
    while len(p) > 1 and p[-1] == 0:
        p = p[:-1]
    if len(p) <= 2:
        return True, 0.0
    roots = mp.polyroots(p[::-1], maxsteps=200, extraprec=50)
    mx = float(max(abs(mp.im(z)) for z in roots))
    return mx < 1e-20, mx

mp.mp.dps = 100
rng = np.random.default_rng(2026)

ts_mp = [mp.mpf(x) for x in [
    "-0.999999", "-0.9999", "-0.99", "-0.9", "-0.5", "-0.1", "-0.01", "-0.001", "-1e-5", "-1e-7"
]]

for mode in ["double_same_T", "double_Dm_w1"]:
    neg_turan_mp = 0
    min_turan_mp = 1e9
    rz_fails_mp = 0
    for trial in range(40):
        P = [[mp.mpf(1)]]
        params = []
        for j in range(9):
            if mode == "double_Dm_w1":
                a = mp.mpf(1)
                b = mp.mpf(1) + mp.power(10, rng.uniform(-3, 5))
                w = mp.mpf(1)
            else:
                a = mp.power(10, rng.uniform(-4, 2))
                b = a + mp.power(10, rng.uniform(-3, 5))
                w = mp.power(10, rng.uniform(-4, 5))
            params.append((float(a), float(b), float(w)))
            P = step_double_mp(P, a, b, w)
        # Check Turan in 100-dps!
        for i in range(1, len(P) - 1):
            for tv in ts_mp:
                v0, v1, v2 = eval_mp(P[i-1], tv), eval_mp(P[i], tv), eval_mp(P[i+1], tv)
                denom = v1*v1 + abs(v0*v2)
                if denom > 0:
                    rel = float((v1*v1 - v0*v2) / denom)
                    if rel < min_turan_mp:
                        min_turan_mp = rel
                    if rel < -1e-12:
                        neg_turan_mp += 1
        # Check RZ of middle polynomials in 100-dps!
        for i in [2, 3, 4, 5, 6]:
            ok, im = check_rz_exact_sturm_mp(P[i], dps=100)
            if not ok:
                rz_fails_mp += 1
                print(f"  *** GENUINE 100-DPS RZ FAIL in {mode}! trial={trial}, i={i}, max_im={im:.4e}")
                print(f"      params={params}")
                break
    print(f"[{mode:14s}] k=9, 40 trials (100 dps): rz_fails={rz_fails_mp}/40, neg_turan_hits={neg_turan_mp}, min_rel_turan={min_turan_mp:+.6f}")
