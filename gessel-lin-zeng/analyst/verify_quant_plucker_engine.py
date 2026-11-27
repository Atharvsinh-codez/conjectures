"""Verify the Quantitative Plucker Engine:
1. Check GLOBAL minimum of J_new = 4(1+a+u)(1+b+v) / ((1+u)(1+v) - ab)^2
   over ALL a, b > 0 with a*b = R in [0, 4] and u in [0, 2*sqrt(a/J_in)], v in [0, 2*sqrt(b/J_in)]!
   Does J_new >= 1 + eps - eps^2 / 4 (or similar exact bound) hold GLOBALLY?!
2. Measure the EXACT shift loss Delta_shift(k) = min_t J(D_{k+1}) - min_t J(D_{k+2}) on a FINE grid of t in R
   for k = 3..15, AND find the EXACT point (i, t) where shift loss is largest!
3. Check the ratio: (shift perturbation W(X, t*f_1)) / W(X, D_{k+1} f_1) on ALL t in R!
   WAIT A SECOND!! STOP!! Look at that ratio:
     W(f_0, t f_1) / W(f_0, D_{k+1} f_1) !!
   If that ratio is uniformly O(1/k) (or O(1/k^2)?!) on ALL t in R, we can BOUND the shift loss in 5 lines!!
"""
import math
import numpy as np
import mpmath as mp

print("=== 1. GLOBAL MIN OF J_new OVER ALL a, b > 0 WITH ab <= 4, u <= 2*sqrt(a/J_in), v <= 2*sqrt(b/J_in) ===")
for eps in [0.4, 0.2, 0.1, 0.04, 0.02, 0.01]:
    J_in = 1.0 + eps
    min_Jn = 1e9; worst = None
    for R in [0.0, 0.1, 0.5, 0.8, max(0.0, 1.0 - 2*eps), 1.0 - eps, 1.0, 1.5, 2.0, 4.0]:
        a_list = [1.0] if R == 0 else [math.sqrt(R) * math.exp(lr) for lr in np.linspace(-4, 4, 41)]
        for a in a_list:
            b = 0.0 if R == 0 else R / a
            um = 2.0 * math.sqrt(a / J_in); vm = 2.0 * math.sqrt(b / J_in)
            for fu in [0.0, 0.5, 0.9, 1.0]:
                for fv in [0.0, 0.5, 0.9, 1.0]:
                    u = fu * um; v = fv * vm
                    den = (1 + u)*(1 + v) - R
                    if den <= 0: continue
                    jn = 4*(1 + a + u)*(1 + b + v) / (den * den)
                    if jn < min_Jn:
                        min_Jn = jn; worst = (R, a, b, fu, fv)
    drop = J_in - min_Jn
    print(f"  eps={eps:.3f}: J_in={J_in:.6f}, min J_new={min_Jn:.6f}, drop/eps^2 = {drop/(eps*eps):.5f} (compare 0.25), worst R={worst[0]:.3f}, a/b={worst[1]/(worst[2]+1e-9):.2f}")

print("\n=== 2. EXACT SHIFT PERTURBATION RATIO W(f_0, t*f_1) / W(f_0, D_{k+1} f_1) AND SHIFT LOSS ON GLZ (k=3..14) ===")
mp.mp.dps = 50
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

KMAX = 14
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D_int(m, H[k-1][0])
    for i in range(1, k): U[i] = poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1): H[k][i] = apply_D_int(m, U[i])

ts_dense = [mp.mpf(s) for s in [
    "-100", "-10", "-2", "-1.2", "-1.02", "-1", "-0.99", "-0.95", "-0.85", "-0.7", "-0.5",
    "-0.3", "-0.15", "-0.05", "-0.01", "-0.002", "-1e-5", "0", "1e-5", "0.002", "0.01", "0.03",
    "0.07", "0.15", "0.3", "0.7", "2", "10", "100"
]]

for k in range(3, KMAX + 1):
    max_w_ratio = 0.0; worst_wr_info = None
    max_J_drop = -1e9; worst_drop_info = None
    min_J_before = 1e9; min_J_after = 1e9
    for i in range(2, k + 1):
        f0 = H[k][i-2]; f1 = H[k][i-1]; f2 = H[k][i]
        tf1 = [0] + list(f1)
        M1 = apply_D_int(k+1, f1); Y1 = apply_D_int(k+1, apply_D_int(k+1, f2))
        M2 = apply_D_int(k+2, f1); Y2 = apply_D_int(k+2, apply_D_int(k+2, f2))
        f0p, tf1p = deriv_int(f0), deriv_int(tf1)
        M1p, Y1p, M2p, Y2p = deriv_int(M1), deriv_int(Y1), deriv_int(M2), deriv_int(Y2)
        for tv in ts_dense:
            v0, v0p = eval_mp(f0, tv), eval_mp(f0p, tv)
            vt, vtp = eval_mp(tf1, tv), eval_mp(tf1p, tv)
            m1, m1p = eval_mp(M1, tv), eval_mp(M1p, tv)
            y1, y1p = eval_mp(Y1, tv), eval_mp(Y1p, tv)
            m2, m2p = eval_mp(M2, tv), eval_mp(M2p, tv)
            y2, y2p = eval_mp(Y2, tv), eval_mp(Y2p, tv)
            w_shift = v0p*vt - v0*vtp
            w_xm1 = v0p*m1 - v0*m1p; w_my1 = m1p*y1 - m1*y1p; w_xy1 = v0p*y1 - v0*y1p
            w_xm2 = v0p*m2 - v0*m2p; w_my2 = m2p*y2 - m2*y2p; w_xy2 = v0p*y2 - v0*y2p
            if w_xm1 != 0:
                wr = float(w_shift / w_xm1)
                if wr > max_w_ratio: max_w_ratio = wr; worst_wr_info = (i, float(tv))
            if w_xy1 != 0 and w_xy2 != 0:
                j1 = float(4 * w_xm1 * w_my1 / (w_xy1 * w_xy1))
                j2 = float(4 * w_xm2 * w_my2 / (w_xy2 * w_xy2))
                min_J_before = min(min_J_before, j1)
                min_J_after = min(min_J_after, j2)
                if j1 - j2 > max_J_drop:
                    max_J_drop = j1 - j2; worst_drop_info = (i, float(tv), j1, j2)
    print(f"  k={k:2d}: max W(f0,tf1)/W(f0,D f1) = {max_w_ratio:.5f} (k*ratio={k*max_w_ratio:.3f} at {worst_wr_info}) | "
          f"minJ: {min_J_before:.4f}->{min_J_after:.4f} (k^2*drop_min={k*k*(min_J_before-min_J_after):.3f}) | "
          f"max pointwise drop={max_J_drop:.5f} at {worst_drop_info[:2]}", flush=True)
