"""1. FIXED Stage 3 pole-gap test at d = 2 (proper b1, b2 sampling so valid > 1000!):
   Given f2=1, f1=1+a*t, f0=(1+b1*t)(1+b2*t) with 1 < b2 < a < b1,
   r_{XM}, r_{MY} <= 0.75 * j, AND J(j) >= 1:
   What is max (J(j) - J(j+1)) as a function of j = 20, 40, 80, 160?!
   Is it O(1/j) or O(1/j^2)?!
2. Also test with STRONGER hypothesis: J(r) >= 1 at BOTH r = 0.75 * j AND r = j!
   Does J(0.75*j) >= 1 AND J(j) >= 1 force J(j) - J(j+1) = O(1/j^2)?!
3. Also check on GLZ (k = 3..14): what is the exact smallest r_J(k) such that
   J(htilde_{k,i-2}, D_r htilde_{k,i-1}, D_r^2 htilde_{k,i}) >= 1 on ALL t in R?!
"""
import numpy as np
import mpmath as mp

def apply_D_flt(m, p):
    res = np.zeros(len(p) + 1)
    for r, c in enumerate(p): res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def eval_W_flt(f, g, ts):
    fp = np.arange(1, len(f))*f[1:]; gp = np.arange(1, len(g))*g[1:]
    return np.polyval(fp[::-1], ts)*np.polyval(g[::-1], ts) - np.polyval(f[::-1], ts)*np.polyval(gp[::-1], ts)

ts_r = np.concatenate([-10**np.linspace(2, -4, 30), np.linspace(-0.99, -0.01, 60), [0.0], 10**np.linspace(-4, 2, 30)])
rng = np.random.default_rng(2026)

print("=== 1. FIXED d=2 STAGE 3 POLE-GAP TEST (r* <= 0.75*j, AND J(0.75*j) >= 1) ===")
for j in [20, 40, 80, 160]:
    r_cut = 0.75 * j
    valid_A = 0; max_drop_A = -1e9; cfg_A = None
    valid_B = 0; max_drop_B = -1e9; cfg_B = None
    for _ in range(15000):
        # Need 1 < b2 < a < b1 with b1 + b2 < 2a + r_cut and b1*b2 < a*(1 + r_cut) and a < r_cut!
        a = rng.uniform(1.05, r_cut * 0.95)
        b2 = rng.uniform(1.01, a * 0.99)
        b1_max = min(2*a + r_cut - b2, a*(1 + r_cut)/b2) * 0.98
        if b1_max <= a * 1.01: continue
        b1 = rng.uniform(a * 1.01, b1_max)
        f2 = np.array([1.0]); f1 = np.array([1.0, a]); f0 = np.array([1.0, b1+b2, b1*b2])
        Mr = apply_D_flt(r_cut, f1); Yr = apply_D_flt(r_cut, apply_D_flt(r_cut, f2))
        w1r = eval_W_flt(f0, Mr, ts_r); w2r = eval_W_flt(Mr, Yr, ts_r); w12r = eval_W_flt(f0, Yr, ts_r)
        if np.max(w1r) >= 0 or np.max(w2r) >= 0: continue
        Jr = float(np.min(4*w1r*w2r / (w12r*w12r + 1e-300)))
        Mj = apply_D_flt(j, f1); Yj = apply_D_flt(j, apply_D_flt(j, f2))
        w1 = eval_W_flt(f0, Mj, ts_r); w2 = eval_W_flt(Mj, Yj, ts_r); w12 = eval_W_flt(f0, Yj, ts_r)
        if np.max(w1) >= 0 or np.max(w2) >= 0: continue
        Jj = float(np.min(4*w1*w2 / (w12*w12 + 1e-300)))
        if Jj < 1.0: continue
        Mj1 = apply_D_flt(j+1, f1); Yj1 = apply_D_flt(j+1, apply_D_flt(j+1, f2))
        w1p = eval_W_flt(f0, Mj1, ts_r); w2p = eval_W_flt(Mj1, Yj1, ts_r); w12p = eval_W_flt(f0, Yj1, ts_r)
        Jj1 = float(np.min(4*w1p*w2p / (w12p*w12p + 1e-300)))
        drop = Jj - Jj1
        valid_A += 1
        if drop > max_drop_A: max_drop_A = drop; cfg_A = (round(b2,2), round(a,2), round(b1,2), round(Jr,4), round(Jj,5), round(Jj1,5))
        if Jr >= 1.0:
            valid_B += 1
            if drop > max_drop_B: max_drop_B = drop; cfg_B = (round(b2,2), round(a,2), round(b1,2), round(Jr,4), round(Jj,5), round(Jj1,5))
    print(f"  j={j:3d} [A: r*<=0.75j, valid={valid_A:5d}]: max drop={max_drop_A:+.6f} (j*drop={j*max_drop_A:+.4f}, j^2*drop={j*j*max_drop_A:+.2f}) cfg={cfg_A}")
    print(f"  j={j:3d} [B: J(0.75j)>=1, val={valid_B:5d}]: max drop={max_drop_B:+.6f} (j*drop={j*max_drop_B:+.4f}, j^2*drop={j*j*max_drop_B:+.2f}) cfg={cfg_B}", flush=True)

print("\n=== 2. EXACT GLZ THRESHOLD r_J(k) WHERE J(htilde_{k,i-2}, D_r htilde_{k,i-1}, D_r^2 htilde_{k,i}) >= 1 ON R ===")
mp.mp.dps = 50
def apply_D_mp(m, poly):
    res = [mp.mpf(0)]*(len(poly)+1)
    for r, c in enumerate(poly): res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def apply_D_int(m, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly): res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv(p): return [r*p[r] for r in range(1, len(p))]
def eval_mp(p, tv):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * tv + c
    return acc

H = [[None]*(k+1) for k in range(13)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, 13):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D_int(m, H[k-1][0])
    for i in range(1, k): U[i] = poly_add(apply_D_int(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1): H[k][i] = apply_D_int(m, U[i])

ts_mp = [mp.mpf(s) for s in ["-5", "-1.5", "-1", "-0.9", "-0.5", "-0.1", "-0.01", "0", "0.01", "0.1", "0.5", "2"]]
for k in [4, 6, 8, 10, 12]:
    i = k // 2 + 1
    f0 = [mp.mpf(c) for c in H[k][i-2]]
    f1 = [mp.mpf(c) for c in H[k][i-1]]
    f2 = [mp.mpf(c) for c in H[k][i]]
    # Binary search smallest r in [0, k+1] where min_t J(r) >= 1 and W1, W2 < 0!
    def ok_r(rv):
        rv = mp.mpf(rv)
        M = apply_D_mp(rv, f1); Y = apply_D_mp(rv, apply_D_mp(rv, f2))
        f0p, Mp, Yp = deriv(f0), deriv(M), deriv(Y)
        for tv in ts_mp:
            v0, v0p = eval_mp(f0, tv), eval_mp(f0p, tv)
            mv, mpv = eval_mp(M, tv), eval_mp(Mp, tv)
            yv, ypv = eval_mp(Y, tv), eval_mp(Yp, tv)
            w1 = v0p*mv - v0*mpv; w2 = mpv*yv - mv*ypv; w12 = v0p*yv - v0*ypv
            if w1 >= 0 or w2 >= 0 or 4*w1*w2 < w12*w12:
                return False
        return True
    lo, hi = 0.0, float(k + 1)
    for _ in range(18):
        mid = 0.5*(lo + hi)
        if ok_r(mid): hi = mid
        else: lo = mid
    print(f"  k={k:2d}, i={i:2d}: r_J(k) = {hi:.4f} (r_J/k = {hi/k:.4f}, gap k+1-r_J = {k+1-hi:.4f})", flush=True)
