"""
Three razor-sharp 1-step closure tests (< 15s total):
TEST 1 (Route III - Derivative Sandwich):
  Given h (deg d+1), g (deg d), f (deg d-1) with ALL negative roots satisfying:
    h' << g << P_{d+1} h   AND   g' << f << P_d g,
  does S_{L-1} f << h ALWAYS hold (for L-1 = d + k - 1, like GLZ)?
  Also test if (S_{L-1} f + g)' << ... (i.e. can Derivative Sandwich itself close in 1 step?!)

TEST 2 (Route I - Link 2 via R <= 1 instead of unbounded W_{0,r} > 0!):
  Remember @verifier (#23): under R = W01*W23/W12^2 <= 1, Quad Theorem B has ZERO defects
  and propagates W(v0+v1, v2+v3) > 0 automatically!
  Test on 5 vectors v0..v4 in R^2 with W(v_r, v_{r+1}) > 0 and W(v_r, v_{r+2}) > 0:
  Suppose BOTH 4-tuples (v0..v3) and (v1..v4) satisfy J >= 1 + eps AND R <= 1 (or R <= R_max).
  What is max R_out = R(v0+v1, v1+v2, v2+v3, v3+v4)? Does R <= 1 (or R <= 4) CLOSE under Stage 1?!

TEST 3 (Route I - Stage 3 via delta_i(t) = W(h_{i-1}, t h_i) / W(h_{i-1}, D_j h_i)):
  Test whether j^2 |delta_{i+1}(t) - delta_i(t)| <= C CLOSES under Stage 1 (addition) and Stage 2 (D_j)!
  Specifically, compute how delta_i(t) transforms under:
    (a) Stage 1: h_i -> U_i = D_j h_i + h_{i-1}
    (b) Stage 2: U_i -> D_j U_i
    (c) Stage 3: j -> j + 1
  on GLZ (k = 5..15) AND on random d=3 triples!
"""
import numpy as np
import mpmath as mp
from scipy.optimize import minimize

mp.mp.dps = 50

# ============================================================
# TEST 1: DERIVATIVE SANDWICH -> 4TH CROSS-PAIR S_{L-1} f << h
# ============================================================
print("=== TEST 1: Does Derivative Sandwich (h'<<g<<P_{d+1}h, g'<<f<<P_d g) imply S_{L-1} f << h? ===")
def poly_deriv(p):
    return [m * p[m] for m in range(1, len(p))]

def apply_Pc(p, c):
    out = [(c - m) * p[m] for m in range(len(p))]
    while len(out) > 1 and abs(out[-1]) < 1e-13 * max(1.0, max(abs(x) for x in out)):
        out.pop()
    return out

def apply_SL(p, L):
    out = [0.0] * (len(p) + 1)
    for m, c in enumerate(p):
        out[m] += (1 + m) * c
        out[m + 1] += (L - m) * c
    return out

def roots_sorted(p):
    r = np.roots(p[::-1])
    if np.max(np.abs(r.imag)) > 1e-6:
        return None
    return np.sort(r.real)

def interlaces_deg1(rp, rq):
    if rp is None or rq is None or len(rq) != len(rp) + 1:
        return False
    return all(rq[m] < rp[m] < rq[m + 1] for m in range(len(rp)))

rng = np.random.default_rng(42)
for d, L_minus_1 in [(4, 6), (5, 8), (7, 11), (9, 14)]:
    fails_4th = 0
    fails_B = 0
    fails_sand_B = 0
    trials = 500
    for _ in range(trials):
        # random negative roots for h (deg d+1)
        rh = -np.sort(np.exp(rng.normal(0, 1.5, d + 1)))[::-1]
        h = list(np.poly(rh)[::-1])
        rhp = roots_sorted(poly_deriv(h))
        rPh = roots_sorted(apply_Pc(h, d + 1))
        # pick g inside (rhp[m], rPh[m])!
        rg = [rng.uniform(rhp[m], rPh[m]) for m in range(d)]
        g = list(np.poly(rg)[::-1])
        rgp = roots_sorted(poly_deriv(g))
        rPg = roots_sorted(apply_Pc(g, d))
        # pick f inside (rgp[m], rPg[m])!
        rf = [rng.uniform(rgp[m], rPg[m]) for m in range(d - 1)]
        f = list(np.poly(rf)[::-1])
        # Check 4th cross-pair: S_{L-1} f << h!
        Sf = apply_SL(f, L_minus_1)
        if not interlaces_deg1(roots_sorted(Sf), rh):
            fails_4th += 1
        # Also check B_1 = S_{L-1} f + g vs B_0 = S_L g + h (with random positive scale factor on f, g, h)!
        B1 = [Sf[m] + g[m] for m in range(d + 1)]
        Sg = apply_SL(g, L_minus_1 + 1)
        B0 = [Sg[m] + h[m] for m in range(d + 2)]
        rB1 = roots_sorted(B1)
        rB0 = roots_sorted(B0)
        if not interlaces_deg1(rB1, rB0):
            fails_B += 1
        # Check if B0' << B1 << P_{d+1} B0!
        if rB1 is not None and rB0 is not None:
            rB0p = roots_sorted(poly_deriv(B0))
            rPB0 = roots_sorted(apply_Pc(B0, d + 1))
            if not (all(rB0p[m] < rB1[m] < rPB0[m] for m in range(d))):
                fails_sand_B += 1
        else:
            fails_sand_B += 1
    print(f"  d={d}, L-1={L_minus_1:2d}: fails(S_{{L-1}} f << h) = {fails_4th}/{trials}, fails(B1 << B0) = {fails_B}/{trials}, fails(B0'<<B1<<P B0) = {fails_sand_B}/{trials}")

# ============================================================
# TEST 2: DOES R <= 1 (or R <= 4) CLOSE UNDER PASCAL ADDITION OF 5 VECTORS?!
# ============================================================
print("\n=== TEST 2: Does R <= 1 (with J >= 1) CLOSE under Stage 1 Pascal addition? ===")
def W2(a, b): return a[0]*b[1] - a[1]*b[0]

def get_JR(vs):
    w01 = W2(vs[0], vs[1]); w12 = W2(vs[1], vs[2]); w23 = W2(vs[2], vs[3])
    w02 = W2(vs[0], vs[2]); w13 = W2(vs[1], vs[3])
    if w01 <= 0 or w12 <= 0 or w23 <= 0 or w02 <= 0 or w13 <= 0:
        return None
    J1 = 4*w01*w12 / (w02**2)
    J2 = 4*w12*w23 / (w13**2)
    R = (w01 * w23) / (w12**2)
    return min(J1, J2), R

for R_cap in [1.0, 2.0, 4.0]:
    max_Rout = 0.0
    valid = 0
    exceed = 0
    for _ in range(50000):
        ang = np.sort(rng.uniform(0.01, np.pi - 0.01, 5))
        lens = np.exp(rng.normal(0, 1.5, 5))
        vs = [(lens[i]*np.cos(ang[i]), lens[i]*np.sin(ang[i])) for i in range(5)]
        res1 = get_JR(vs[:4])
        res2 = get_JR(vs[1:])
        if res1 is None or res2 is None:
            continue
        if res1[0] < 1.0 or res2[0] < 1.0 or res1[1] > R_cap or res2[1] > R_cap:
            continue
        valid += 1
        ws = [(vs[i][0]+vs[i+1][0], vs[i][1]+vs[i+1][1]) for i in range(4)]
        reso = get_JR(ws)
        if reso is None:
            exceed += 1
        else:
            if reso[1] > max_Rout:
                max_Rout = reso[1]
            if reso[1] > R_cap + 1e-12:
                exceed += 1
    print(f"  R_in <= {R_cap:.1f} (with J>=1, open half-plane): valid={valid:4d}, max R_out = {max_Rout:.6f}, exceed={exceed}")

# Also run adversarial Nelder-Mead trying to push R_out > 1 when J1, J2 >= 1 and R1, R2 <= 1!
best_Rout = 0.0
best_info = None
for seed in range(25):
    rng2 = np.random.default_rng(500 + seed)
    found = False
    for _ in range(1000):
        ang = np.sort(rng2.uniform(0.05, np.pi - 0.05, 5))
        lens = np.exp(rng2.normal(0, 1.5, 5))
        vs = [(lens[i]*np.cos(ang[i]), lens[i]*np.sin(ang[i])) for i in range(5)]
        r1 = get_JR(vs[:4]); r2 = get_JR(vs[1:])
        if r1 and r2 and r1[0] >= 1.0 and r2[0] >= 1.0 and r1[1] <= 1.0 and r2[1] <= 1.0:
            found = True
            break
    if not found: continue
    x0 = np.concatenate([np.log(np.diff(ang)), np.log(lens[[0,1,3,4]] / lens[2])])
    def objR(z):
        g = np.exp(np.clip(z[:4], -12, 4))
        tot = np.sum(g)
        if tot >= np.pi - 1e-4: g = g * (np.pi - 0.01) / (tot + 0.01)
        th = np.zeros(5); th[0] = 0.005; th[1:] = th[0] + np.cumsum(g)
        L = np.exp(np.clip([z[4], z[5], 0.0, z[6], z[7]], -10, 10))
        vs = [(L[i]*np.cos(th[i]), L[i]*np.sin(th[i])) for i in range(5)]
        r1 = get_JR(vs[:4]); r2 = get_JR(vs[1:])
        if r1 is None or r2 is None: return 10.0
        pen = 0.0
        if r1[0] < 1.0: pen += 100*(1.0 - r1[0])
        if r2[0] < 1.0: pen += 100*(1.0 - r2[0])
        if r1[1] > 1.0: pen += 100*(r1[1] - 1.0)
        if r2[1] > 1.0: pen += 100*(r2[1] - 1.0)
        if pen > 0: return 10.0 + pen
        ws = [(vs[i][0]+vs[i+1][0], vs[i][1]+vs[i+1][1]) for i in range(4)]
        ro = get_JR(ws)
        return -ro[1] if ro else 10.0
    res = minimize(objR, x0, method='Nelder-Mead', options={'maxiter': 1500})
    if -res.fun > best_Rout:
        best_Rout = -res.fun
print(f"  Adversarial R_in <= 1, J_in >= 1: max R_out = {best_Rout:.8f}!")

# ============================================================
# TEST 3: TRACK j^2 |delta_{i+1}(t) - delta_i(t)| ACROSS STAGES 1, 2, 3 ON GLZ!
# ============================================================
print("\n=== TEST 3: GLZ j^2 |delta_{i+1}(t) - delta_i(t)| ACROSS STAGE 1 (h->U), STAGE 2 (U->D_j U), STAGE 3 (j->j+1) ===")
def apply_Dj_mp(poly, j):
    out = [mp.mpf(0)] * (len(poly) + 1)
    for m, c in enumerate(poly):
        if c != 0:
            out[m] += (1 + m) * c
            out[m + 1] += (j + m) * c
    return out

def padd_mp(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for idx in range(len(a)): out[idx] += a[idx]
    for idx in range(len(b)): out[idx] += b[idx]
    return out

def mul_t_mp(poly):
    return [mp.mpf(0)] + list(poly)

def eval_pd_mp(poly, tv):
    v = mp.mpf(0); dv = mp.mpf(0)
    for c in reversed(poly):
        dv = dv * tv + v
        v = v * tv + c
    return v, dv

def max_delta_diff(polys, j, ts):
    # delta_i(t) = W(polys[i-1], t*polys[i]) / W(polys[i-1], D_j polys[i])
    K = len(polys) - 1
    deltas = {}
    max_d = 0.0
    for i in range(1, K + 1):
        f0, f1 = polys[i - 1], polys[i]
        tf1 = mul_t_mp(f1)
        Df1 = apply_Dj_mp(f1, j)
        deltas[i] = []
        for t in ts:
            v0 = eval_pd_mp(f0, t); vt = eval_pd_mp(tf1, t); vD = eval_pd_mp(Df1, t)
            wn = v0[0]*vt[1] - v0[1]*vt[0]
            wd = v0[0]*vD[1] - v0[1]*vD[0]
            dval = float(wn / wd)
            deltas[i].append(dval)
            if j * dval > max_d: max_d = j * dval
    max_diff = 0.0
    for i in range(1, K):
        for ti in range(len(ts)):
            df = abs(deltas[i + 1][ti] - deltas[i][ti]) * (j ** 2)
            if df > max_diff: max_diff = df
    return max_d, max_diff

ts_d = [mp.mpf(x) for x in [-10.0, -1.2, -0.95, -0.8, -0.5, -0.2, -0.05, -0.005, 0.0, 0.05, 0.2, 1.0, 10.0]]
h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
for k in range(2, 13):
    j = k + 1
    U = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = padd_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = padd_mp(u, h[i - 1])
        U.append(u)
    hn = [apply_Dj_mp(u, j) for u in U]
    if k in [4, 6, 8, 10, 12]:
        d_in = max_delta_diff(h, j, ts_d)
        d_S1 = max_delta_diff(U, j, ts_d)
        d_S2 = max_delta_diff(hn, j, ts_d)
        d_S3 = max_delta_diff(hn, j + 1, ts_d)
        print(f"  k={k:2d} (j={j:2d}): [j*max_d, j^2*max|d_{{i+1}}-d_i|] | in={d_in[0]:.3f},{d_in[1]:.3f} -> S1(U)={d_S1[0]:.3f},{d_S1[1]:.3f} -> S2(h_k,j)={d_S2[0]:.3f},{d_S2[1]:.3f} -> S3(h_k,j+1)={d_S3[0]:.3f},{d_S3[1]:.3f}")
    h = hn
