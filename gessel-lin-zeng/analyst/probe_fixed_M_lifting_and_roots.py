import mpmath as mp
import numpy as np

mp.mp.dps = 80

def apply_D_int(m, p):
    res = [0]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv_int(p): return [r*p[r] for r in range(1, len(p))] if len(p) > 1 else [0]
def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v
def eval_W(f, g, t):
    return eval_mp(f, t)*eval_mp(deriv_int(g), t) - eval_mp(deriv_int(f), t)*eval_mp(g, t)

KMAX = 14
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    H[k][0] = apply_D_int(m, apply_D_int(m, H[k-1][0]))
    for i in range(1, k):
        H[k][i] = apply_D_int(m, poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1]))
    H[k][k] = apply_D_int(m, H[k-1][k-1])

ts = [mp.mpf(s) for s in ["-10", "-2", "-1.1", "-0.99", "-0.9", "-0.5", "-0.1", "-0.01", "-1e-5", "0", "0.01", "0.1", "1", "10"]]

print("=== 1. FIXED-M LIFTING: J(tilde_h_{k,i-1}, D_M tilde_h_{k,i}, D_M^2 tilde_h_{k,i+1}) FOR FIXED M = 1, 2, k+1, 20, 50 ===")
for M in [1, 2, 3, 15, 30, 100]:
    print(f"\n  --- Fixed M = {M} ---")
    for k in [4, 7, 10, 14]:
        mn_J = mp.mpf('1e100'); loc = None
        for i in range(1, k):
            X = H[k][i-1]
            Mid = apply_D_int(M, H[k][i])
            Y = apply_D_int(M, apply_D_int(M, H[k][i+1]))
            for tv in ts:
                w1 = eval_W(X, Mid, tv); w2 = eval_W(Mid, Y, tv); w12 = eval_W(X, Y, tv)
                if w1 <= 0 or w2 <= 0:
                    mn_J = mp.mpf(-1); loc = (i, float(tv), "W<=0"); break
                J = 4*w1*w2/(w12*w12) if w12 != 0 else mp.mpf('1e10')
                if J < mn_J: mn_J = J; loc = (i, float(tv))
            if mn_J < 0: break
        print(f"    k={k:2d}: min_J(h_{{i-1}}, D_{M} h_i, D_{M}^2 h_{{i+1}}) = {float(mn_J):+9.5f} at {loc}")

print("\n=== 2. ALSO TEST FULL D_M^i LIFTING: (D_M^{i-1} h_{k,i-1}, D_M^i h_{k,i}, D_M^{i+1} h_{k,i+1}) AT M = 1, 2, 15 ===")
for M in [1, 2, 15]:
    for k in [6, 10, 14]:
        # build F[i] = D_M^i H[k][i]
        F = []
        for i in range(k + 1):
            p = H[k][i]
            for _ in range(i): p = apply_D_int(M, p)
            F.append(p)
        mn_J = mp.mpf('1e100'); loc = None
        for i in range(1, k):
            for tv in ts:
                w1 = eval_W(F[i-1], F[i], tv); w2 = eval_W(F[i], F[i+1], tv); w12 = eval_W(F[i-1], F[i+1], tv)
                if w1 <= 0 or w2 <= 0:
                    mn_J = mp.mpf(-1); loc = (i, float(tv), "W<=0"); break
                J = 4*w1*w2/(w12*w12) if w12 != 0 else mp.mpf('1e10')
                if J < mn_J: mn_J = J; loc = (i, float(tv))
            if mn_J < 0: break
        print(f"  Full D_{M}^i at k={k:2d}: min_J = {float(mn_J):+9.5f} at {loc}")

print("\n=== 3. INSPECT s-ROOTS OF A_{k, i}(-e^s) AT k = 12 FOR i = 0, 4, 6, 9, 12 ===")
for i in [0, 4, 6, 9, 12]:
    # roots of tilde_h_{12, i}(t) in (-1, 0), mapped to s = log(-x) = log(-t/(1+t))!
    p = H[12][i]
    rts = sorted(float(mp.re(z)) for z in mp.polyroots(p[::-1], maxsteps=200, extraprec=40))
    s_rts = [np.log(-t / (1 + t)) for t in rts]
    diffs = np.diff(s_rts)
    print(f"  k=12, i={i:2d} (deg={len(s_rts):2d}): s in [{s_rts[0]:+6.2f}, {s_rts[-1]:+6.2f}], min_gap={diffs.min():.3f}, max_gap={diffs.max():.3f}, median_gap={np.median(diffs):.3f}")
