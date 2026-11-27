"""Diagnostic of:
1. Normalized D-lifted family ghat_{k, i}(t) = (1/binom(k,i)) D_{k+1}^i htilde_{k,i}(t) for k = 2..16:
   - min_{t in R} 4 W(ghat_{i-1}, ghat_i) W(ghat_i, ghat_{i+1}) / W(ghat_{i-1}, ghat_{i+1})^2  (vs 1 - C/k^2!)
   - min/max_{t in R} W(ghat_{i-1}, ghat_i) / W(ghat_i, ghat_{i+1})  (how close to 1?!)
   - max_{t in R} W(ghat_{i-1}, ghat_i) W(ghat_{i+1}, ghat_{i+2}) / W(ghat_i, ghat_{i+1})^2
2. Unlifted inf_{t in R} W(h_{i-1}, h_i) W(h_i, h_{i+1}) / W(h_{i-1}, h_{i+1})^2 for k = 2..16
   (gives exact upper bound on q for Q_{k, q}(t, y) stability!).
3. Test if GLZ-shaped product prod_{j=1}^k D_{m_j}(D_{m_j} + y)(1) (with outer D_{m_j}, c=1, w=1)
   can be broken for adversarial wild m_j (e.g. step jumps, huge spikes, alternating)!
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 60

def apply_D_int(m, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv_int(p): return [r*p[r] for r in range(1, len(p))]

def eval_mp(p, tv):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * tv + c
    return acc

KMAX = 16
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D_int(m, H[k-1][0])
    for i in range(1, k): U[i] = poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1): H[k][i] = apply_D_int(m, U[i])

ts_chk = [mp.mpf(s) for s in ["-10", "-2", "-1.1", "-1", "-0.99", "-0.9", "-0.7", "-0.5", "-0.3", "-0.1", "-0.02", "-0.002", "-1e-6", "0", "1e-6", "0.01", "0.2", "1", "10"]]

print("=== 1. NORMALIZED ghat_{k,i} = (1/C(k,i)) D_{k+1}^i htilde_{k,i} ON R (k=2..16) ===")
for k in range(2, KMAX+1):
    m = k + 1
    G = []
    for i in range(k+1):
        p = list(H[k][i])
        for _ in range(i): p = apply_D_int(m, p)
        G.append(p)
    Gp = [deriv_int(p) for p in G]
    min_Jhat = 1e9; worst_J_info = None
    min_wrat = 1e9; max_wrat = 0.0
    max_Rhat = 0.0
    for tv in ts_chk:
        vs = [eval_mp(G[i], tv) / mp.mpf(math.comb(k, i)) for i in range(k+1)]
        vps = [eval_mp(Gp[i], tv) / mp.mpf(math.comb(k, i)) for i in range(k+1)]
        w1 = [vps[i]*vs[i+1] - vs[i]*vps[i+1] for i in range(k)]
        w2 = [vps[i]*vs[i+2] - vs[i]*vps[i+2] for i in range(k-1)]
        for i in range(k-1):
            if w2[i] != 0:
                jr = float(4 * w1[i] * w1[i+1] / (w2[i] * w2[i]))
                if jr < min_Jhat:
                    min_Jhat = jr; worst_J_info = (i+1, float(tv))
            if w1[i+1] != 0:
                wr = float(w1[i] / w1[i+1])
                min_wrat = min(min_wrat, wr); max_wrat = max(max_wrat, wr)
        for i in range(k-2):
            if w1[i+1] != 0:
                rp = float(w1[i] * w1[i+2] / (w1[i+1] * w1[i+1]))
                max_Rhat = max(max_Rhat, rp)
    defect = (1.0 - min_Jhat) * (k * k)
    print(f"  k={k:2d}: min_Jhat={min_Jhat:.6f} (k^2*(1-Jhat)={defect:+.3f} at {worst_J_info}), "
          f"W_i/W_{{i+1}} in [{min_wrat:.4f}, {max_wrat:.4f}], max_Rhat={max_Rhat:.5f}", flush=True)

print("\n=== 2. UNLIFTED htilde_{k,i}: inf_{t in R} W_{i-1,i} W_{i,i+1} / W_{i-1,i+1}^2 (k=2..16) ===")
for k in range(2, KMAX+1):
    Hp = [deriv_int(H[k][i]) for i in range(k+1)]
    min_q = 1e9; worst_info = None
    for tv in ts_chk:
        vs = [eval_mp(H[k][i], tv) for i in range(k+1)]
        vps = [eval_mp(Hp[i], tv) for i in range(k+1)]
        w1 = [vps[i]*vs[i+1] - vs[i]*vps[i+1] for i in range(k)]
        w2 = [vps[i]*vs[i+2] - vs[i]*vps[i+2] for i in range(k-1)]
        for i in range(k-1):
            if w2[i] != 0:
                val = float(w1[i] * w1[i+1] / (w2[i] * w2[i]))
                if val < min_q: min_q = val; worst_info = (i+1, float(tv))
    print(f"  k={k:2d}: min W01*W12/W02^2 = {min_q:.6f} (4*min = {4*min_q:.5f}) at {worst_info}", flush=True)

print("\n=== 3. ADVERSARIAL WILD m_j TEST ON GLZ-SHAPED prod D_{m_j}(D_{m_j} + y)(1) ===")
# Try to break real-rootedness or Turan on (-1,0) of prod_{j=1}^k D_{m_j}(D_{m_j} + y)(1) with wild m_j > 1!
def test_wild_double_D(label, ms):
    polys = [[1]]
    for m in ms:
        U = [[0]] * (len(polys) + 1)
        for i, p in enumerate(polys):
            U[i] = poly_add_int(U[i], apply_D_int(m, p))
            U[i+1] = poly_add_int(U[i+1], p)
        polys = [apply_D_int(m, u) for u in U]
    k = len(ms)
    # Check real-rootedness via scaled mpmath? Faster: check Turan on (-1,0) AND roots for each i!
    mp.mp.dps = 50
    nonreal = 0; turan_neg = 0
    for i in range(1, k):
        # Check Turan f_i^2 - f_{i-1} f_{i+1} on (-1, 0)
        for tv in [mp.mpf(s) for s in ["-0.999", "-0.99", "-0.9", "-0.5", "-0.1", "-0.01", "-1e-4", "-1e-7", "-1e-10"]]:
            vi = eval_mp(polys[i], tv); vm = eval_mp(polys[i-1], tv); vp = eval_mp(polys[i+1], tv)
            if vi*vi - vm*vp < 0:
                turan_neg += 1; break
    # Check roots of middle polynomial i = k//2
    p = polys[k // 2]; d = len(p) - 1
    g = (mp.mpf(p[0]) / mp.mpf(p[-1])) ** (mp.mpf(1) / d)
    cs = [mp.mpf(p[r]) * (g ** r) for r in range(d + 1)]
    rts = mp.polyroots(cs[::-1], maxsteps=300, extraprec=40)
    max_im = float(max(abs(mp.im(r)) for r in rts))
    print(f"  {label:38s}: turan_neg_i={turan_neg}, mid max|Im(t)|={max_im:.2e}", flush=True)

test_wild_double_D("step [2]*5 + [10000]*5 (k=10)", [2]*5 + [10000]*5)
test_wild_double_D("step [10000]*5 + [2]*5 (k=10)", [10000]*5 + [2]*5)
test_wild_double_D("alt [2, 10000]*5 (k=10)", [2, 10000]*5)
test_wild_double_D("exp 2^j (k=10)", [2**j for j in range(1, 11)])
test_wild_double_D("rev exp 2^(11-j) (k=10)", [2**(11-j) for j in range(1, 11)])
test_wild_double_D("spike [2]*4 + [10^6] + [2]*5", [2]*4 + [10**6] + [2]*5)
