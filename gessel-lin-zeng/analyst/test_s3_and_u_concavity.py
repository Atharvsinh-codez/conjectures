"""1. FIXED d=2 Stage 3 pole-gap test (sample f0 roots directly inside M_r's roots so valid ~ 5000!):
   Test THREE hypotheses on (f0, f1, f2):
   (A) Pole gap only: r_{XM}, r_{MY} <= 0.75 * j and J(j) >= 1.
   (B) Interval J >= 1: J(0.75 * j) >= 1 AND J(j) >= 1.
   (C) Stronger interval J >= 1: J(0.6 * j) >= 1 AND J(j) >= 1.
   Measure max (J(j) - J(j+1)) for j = 20, 40, 80, 160 -> is it O(1/j) or O(1/j^2)?!
2. Also test on GLZ (k = 4..12): is u_i(t) = W(v_{i-1}, v_{i+1}) / (2 sqrt(W_{i-1,i} W_{i,i+1})) CONCAVE in i
   (u_{i-1} + u_{i+1} <= 2 u_i) for ALL t in R?!
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

print("=== 1. FIXED d=2 STAGE 3 POLE-GAP TEST (PROPER ROOT-INTERLACING SAMPLING!) ===")
for j in [20, 40, 80, 160]:
    r_cut = 0.75 * j
    r_cut2 = 0.60 * j
    valA = valB = valC = 0
    dropA = dropB = dropC = -1e9
    cfgA = cfgB = cfgC = None
    for _ in range(15000):
        # f2 = 1, f1 = 1 + a*t (root r1 = -1/a in (-1, 0), a > 1).
        # D_{r_cut} f1 has 2 roots m1 < -1/a < m2 < 0!
        # Pick a in (1.05, 0.9 * r_cut):
        a = 1.0 + rng.uniform(0.05, r_cut * 0.9 - 1.0)
        f2 = np.array([1.0]); f1 = np.array([1.0, a])
        Mr = apply_D_flt(r_cut, f1)
        m_rts = np.sort(np.roots(Mr[::-1]).real)  # -1 < m1 < -1/a < m2 < 0
        m1, m2 = m_rts[0], m_rts[1]
        # For f0 << M_r (roots x1 < m1 < x2 < m2) AND f1 << f0 (x1 < -1/a < x2):
        # We need: -1 < x1 < m1, AND max(-1/a, m1) = -1/a < x2 < m2!
        #Wait: CHECK! -1/a is ALREADY > m1 and < m2!
        # So picking x1 in (-1, m1) and x2 in (-1/a, m2) GUARANTEES BOTH f1 << f0 AND f0 < M_r!!
        x1 = rng.uniform(-0.9999, m1 - 1e-6 * (m1 + 1.0))
        x2 = rng.uniform(-1.0/a + 1e-6*(m2 + 1.0/a), m2 - 1e-6*(m2 + 1.0/a))
        b2 = -1.0 / x1; b1 = -1.0 / x2
        f0 = np.array([1.0, b1 + b2, b1 * b2])

        # Verify r_cut interlacing (W_XM < 0 and W_MY < 0)
        Yr = apply_D_flt(r_cut, apply_D_flt(r_cut, f2))
        w1r = eval_W_flt(f0, Mr, ts_r); w2r = eval_W_flt(Mr, Yr, ts_r); w12r = eval_W_flt(f0, Yr, ts_r)
        if np.max(w1r) >= -1e-12 or np.max(w2r) >= -1e-12: continue
        Jr = float(np.min(4*w1r*w2r / (w12r*w12r + 1e-300)))

        Mj = apply_D_flt(j, f1); Yj = apply_D_flt(j, apply_D_flt(j, f2))
        w1 = eval_W_flt(f0, Mj, ts_r); w2 = eval_W_flt(Mj, Yj, ts_r); w12 = eval_W_flt(f0, Yj, ts_r)
        if np.max(w1) >= -1e-12 or np.max(w2) >= -1e-12: continue
        Jj = float(np.min(4*w1*w2 / (w12*w12 + 1e-300)))
        if Jj < 1.0: continue

        Mj1 = apply_D_flt(j+1, f1); Yj1 = apply_D_flt(j+1, apply_D_flt(j+1, f2))
        w1p = eval_W_flt(f0, Mj1, ts_r); w2p = eval_W_flt(Mj1, Yj1, ts_r); w12p = eval_W_flt(f0, Yj1, ts_r)
        Jj1 = float(np.min(4*w1p*w2p / (w12p*w12p + 1e-300)))
        drop = Jj - Jj1

        valA += 1
        if drop > dropA: dropA = drop; cfgA = (round(a,2), round(x1,4), round(x2,4), round(Jr,4), round(Jj,5), round(Jj1,5))
        if Jr >= 1.0:
            valB += 1
            if drop > dropB: dropB = drop; cfgB = (round(a,2), round(x1,4), round(x2,4), round(Jr,4), round(Jj,5), round(Jj1,5))
            # Check C: J(0.6 * j) >= 1
            Mr2 = apply_D_flt(r_cut2, f1); Yr2 = apply_D_flt(r_cut2, apply_D_flt(r_cut2, f2))
            w1c = eval_W_flt(f0, Mr2, ts_r); w2c = eval_W_flt(Mr2, Yr2, ts_r); w12c = eval_W_flt(f0, Yr2, ts_r)
            if np.max(w1c) < -1e-12 and np.max(w2c) < -1e-12:
                Jc = float(np.min(4*w1c*w2c / (w12c*w12c + 1e-300)))
                if Jc >= 1.0:
                    valC += 1
                    if drop > dropC: dropC = drop; cfgC = (round(a,2), round(Jc,4), round(Jj,5), round(Jj1,5))
    print(f"  j={j:3d} [A r*<=0.75j, n={valA:5d}]: drop={dropA:+.6f} (j*d={j*dropA:+.3f}, j^2*d={j*j*dropA:+.1f}) {cfgA}")
    print(f"  j={j:3d} [B J(0.75j)>=1, n={valB:5d}]: drop={dropB:+.6f} (j*d={j*dropB:+.3f}, j^2*d={j*j*dropB:+.1f}) {cfgB}")
    print(f"  j={j:3d} [C J(0.60j)>=1, n={valC:5d}]: drop={dropC:+.6f} (j*d={j*dropC:+.3f}, j^2*d={j*j*dropC:+.1f}) {cfgC}", flush=True)

print("\n=== 2. CHECK CONCAVITY OF u_i(t) = 1/sqrt(J_i(t)) IN i FOR GLZ (k = 4..12) ===")
mp.mp.dps = 50
def apply_D_int(m, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly): res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv_int(p): return [r*p[r] for r in range(1, len(p))]
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

ts_mp = [mp.mpf(s) for s in ["-5", "-1.1", "-1", "-0.8", "-0.5", "-0.1", "-0.01", "0", "0.01", "0.2", "2"]]
for k in range(4, 12):
    m = k + 1
    G = []
    for i in range(k+1):
        p = list(H[k][i])
        for _ in range(i): p = apply_D_int(m, p)
        G.append(p)
    Gp = [deriv_int(p) for p in G]
    min_conc = 1e9; worst_info = None
    for tv in ts_mp:
        vs = [eval_mp(p, tv) for p in G]; vps = [eval_mp(dp, tv) for dp in Gp]
        w1 = [vps[i]*vs[i+1] - vs[i]*vps[i+1] for i in range(k)]
        w2 = [vps[i]*vs[i+2] - vs[i]*vps[i+2] for i in range(k-1)]
        us = [float(abs(w2[i]) / (2 * mp.sqrt(abs(w1[i] * w1[i+1])))) for i in range(k-1)]
        for i in range(1, len(us) - 1):
            c_val = 2*us[i] - us[i-1] - us[i+1]
            if c_val < min_conc: min_conc = c_val; worst_info = (i+1, float(tv))
    print(f"  k={k:2d}: min(2 u_i - u_{{i-1}} - u_{{i+1}}) = {min_conc:+.6f} at {worst_info}")
