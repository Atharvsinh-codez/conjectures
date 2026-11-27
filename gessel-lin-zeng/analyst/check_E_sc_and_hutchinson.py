import mpmath as mp
import numpy as np

mp.mp.dps = 120

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

# 1. Check E_{k, i} (M = k+2) at t = +-alpha / s_{k,i} AND across all R up to k = 40!
KMAX = 40
row = [[1, 2], [1]]
print("=== 1. SANITY CHECK OF E_{k+1, i} BASE QUAD (h_{i-1}, D_{k+2} h_i, D_{k+2}^2 h_{i+1}) AT t = +-alpha/s_{k,i} UP TO k = 40 ===")
for k in range(2, KMAX + 1):
    m = k + 1
    new_row = [None]*(k + 1)
    new_row[0] = apply_D_int(m, apply_D_int(m, row[0]))
    for i in range(1, k):
        new_row[i] = apply_D_int(m, poly_add_int(apply_D_int(m, row[i]), row[i-1]))
    new_row[k] = apply_D_int(m, row[k-1])
    row = new_row

    if k in [10, 20, 30, 35, 40]:
        j = k + 2
        mn_J = mp.mpf('1e100'); loc = None
        for i in [1, 2, k//4, k//2, 3*k//4, k-1]:
            X = row[i-1]
            M = apply_D_int(j, row[i])
            Y = apply_D_int(j, apply_D_int(j, row[i+1]))
            s1 = mp.mpf(row[i][1]) / mp.mpf(row[i][0])
            t_pts = [mp.mpf(s) for s in ["-10", "-1.1", "-0.9", "-0.5", "-0.1", "0", "0.1", "1", "10"]]
            for mult in ["0.01", "0.1", "0.5", "1", "2", "10", "100"]:
                t_pts.append(mp.mpf(mult) / s1)
                t_pts.append(-mp.mpf(mult) / s1)
            for tv in t_pts:
                w1 = eval_W(X, M, tv); w2 = eval_W(M, Y, tv); w12 = eval_W(X, Y, tv)
                if w1 <= 0 or w2 <= 0:
                    mn_J = mp.mpf(-1); loc = (i, float(tv), "W<=0"); break
                J = 4*w1*w2/(w12*w12) if w12 != 0 else mp.mpf('1e10')
                if J < mn_J: mn_J = J; loc = (i, float(tv))
        print(f"  k={k:2d} (j={j:2d}): min_J = {float(mn_J):+.6f} (k*(J-1) = {k*(float(mn_J)-1):+.4f}) at {loc}")

# 2. Check Hutchinson ratio a_r^2 / (a_{r-1} a_{r+1}) for A_{k,i}(x) at k = 10, 20!
def apply_S_int(p, L):
    d = len(p) - 1
    res = [0]*(d + 2)
    for n in range(d + 2):
        v = 0
        if n <= d: v += (n + 1) * p[n]
        if 1 <= n <= d + 1: v += (L - n + 1) * p[n - 1]
        res[n] = v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

A = {(0, 0): [1]}
for k in range(1, 21):
    for i in range(k + 1):
        u1 = apply_S_int(A[(k-1, i)], 3*k - 3 - i) if i <= k - 1 else [0]
        u2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S_int(poly_add_int(u1, u2), 3*k - 2 - i)

print("\n=== 2. HUTCHINSON RATIOS a_r^2 / (a_{r-1} a_{r+1}) OF A_{k,i}(x) AT k = 20 ===")
for i in [0, 4, 10, 17, 20]:
    a = A[(20, i)]
    rats = [float(mp.mpf(a[r])**2 / (mp.mpf(a[r-1])*mp.mpf(a[r+1]))) for r in range(1, len(a)-1)]
    n_gt4 = sum(1 for x in rats if x >= 4.0)
    print(f"  k=20, i={i:2d} (deg={len(a)-1:2d}): min_rat={min(rats):.3f}, count(>=4)={n_gt4}/{len(rats)}, left5={[round(x,2) for x in rats[:5]]}, right5={[round(x,2) for x in rats[-5:]]}")
