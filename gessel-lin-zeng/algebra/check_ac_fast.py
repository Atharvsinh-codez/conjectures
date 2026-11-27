"""Fast mpmath check for (a) Plucker ratio at k=8..11 and (c) Turan/ULC on [-1,0] at k=35,36."""
import math
import mpmath as mp
mp.mp.dps = 50

def poly_add(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]
def poly_sub(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(n)]
def poly_mul(a, b):
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j, bj in enumerate(b):
            res[i+j] += ai * bj
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res
def poly_deriv(a):
    return [i * a[i] for i in range(1, len(a))] if len(a) > 1 else [0]
def apply_Dj(p, j):
    dp = poly_deriv(p)
    return poly_add(poly_add([0]+dp, [0,0]+dp), poly_add(p, [0]+[j*c for c in p]))
def wronskian(f, g):
    return poly_sub(poly_mul(f, poly_deriv(g)), poly_mul(poly_deriv(f), g))
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
            p1 = apply_Dj(h[k-1][i], k + 1) if i <= k - 1 else [0]
            p2 = h[k-1][i-1] if i >= 1 else [0]
            h[k][i] = apply_Dj(poly_add(p1, p2), k + 1)
    return h

print("=== CHECK (a): Plucker ratio W_AB*W_CE/W_BC^2 and sign of W_AC at k=8..11 (50-digit mpmath) ===", flush=True)
ht = compute_htilde(11)
grid = [mp.mpf(x) for x in [-1000, -100, -20, -5, -2, -1.5, -1.1, -1.01, -1, -0.99, -0.9, -0.7, -0.5, -0.3, -0.1, -0.01, 0, 0.1, 1, 10, 100, 1000]]
grid += [-mp.mpf(m)/40 for m in range(1, 40)]
for k in [8, 9, 10, 11]:
    max_ratio = mp.mpf(0)
    worst = None
    wac_neg = 0
    for i in range(3, k):
        A_p = ht[k-1][i-3]
        B_p = apply_Dj(ht[k-1][i-2], k+1)
        C_p = apply_Dj(apply_Dj(ht[k-1][i-1], k+1), k+1)
        E_p = apply_Dj(apply_Dj(apply_Dj(ht[k-1][i], k+1), k+1), k+1)
        wAB = wronskian(A_p, B_p)
        wBC = wronskian(B_p, C_p)
        wCE = wronskian(C_p, E_p)
        wAC = wronskian(A_p, C_p)
        s0 = 1 if wAC[0] > 0 else -1
        for tv in grid:
            vAB = poly_eval_mp(wAB, tv)
            vBC = poly_eval_mp(wBC, tv)
            vCE = poly_eval_mp(wCE, tv)
            vAC = poly_eval_mp(wAC, tv)
            if vAC * s0 <= 0:
                wac_neg += 1
            r = (vAB * vCE) / (vBC * vBC)
            if r > max_ratio:
                max_ratio = r
                worst = (i, float(tv))
    print(f"  k={k}: max W_AB*W_CE/W_BC^2 = {mp.nstr(max_ratio, 10)} at {worst}; W_AC sign flips = {wac_neg}", flush=True)

print("\n=== CHECK (c): Turan & ULC of htilde_{k,i}(t) on [-1, 0] at k=35, 36 ===", flush=True)
ht36 = compute_htilde(36)
tgrid = [-mp.mpf(10)**(-p) for p in [1, 2, 3, 4, 6, 8, 12, 16, 20]] + [-mp.mpf(m)/50 for m in range(1, 50)] + [mp.mpf(-1) + mp.mpf(10)**(-p) for p in [3, 6, 10]]
for k in [35, 36]:
    min_turan_ratio = mp.mpf('1e100') # h_i^2 / (h_{i-1} h_{i+1}) when h_{i-1}h_{i+1} > 0
    min_ulc_ratio = mp.mpf('1e100')
    turan_neg = 0
    ulc_neg = 0
    worst_u = None
    for i in range(1, k):
        c_ratio = mp.mpf(math.comb(k, i)**2) / (math.comb(k, i-1) * math.comb(k, i+1))
        for tv in tgrid:
            v0 = poly_eval_mp(ht36[k][i], tv)
            vm = poly_eval_mp(ht36[k][i-1], tv)
            vp = poly_eval_mp(ht36[k][i+1], tv)
            if v0*v0 - vm*vp < 0:
                turan_neg += 1
            if v0*v0 - c_ratio * vm*vp < 0:
                ulc_neg += 1
            if vm*vp > 0:
                r = (v0*v0) / (vm*vp)
                if r < min_turan_ratio:
                    min_turan_ratio = r
                ru = r / c_ratio
                if ru < min_ulc_ratio:
                    min_ulc_ratio = ru
                    worst_u = (i, float(tv))
    print(f"  k={k}: turan_neg={turan_neg}, ulc_neg={ulc_neg}, min(h^2/(h- h+))={mp.nstr(min_turan_ratio, 8)}, min(ULC ratio)={mp.nstr(min_ulc_ratio, 8)} at {worst_u}", flush=True)
