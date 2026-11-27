"""Fast decisive check for collaborating with @algebraist:
1. Test Q_{k, q}(t, y) = sum q^{i(i-1)/2} htilde_{k,i}(t) y^i for q = 1/2, 1/4, 1/8, 1/16:
   (a) Full Rayleigh R(t, y) on R^2 (including TINY t > 0 like 1e-12, 1e-8, 1e-5, 1e-3, AND negative t!) for k = 2..16!
   (b) Single-step Herglotz residue R(t0, y0) for q = 1/4, 1/8, 1/16!
2. Test D-lifted q-deformed G_{k, q}(t, y) = sum q^{i(i-1)/2} (D_{k+1}^i htilde_{k,i}(t)) y^i for q = 1/2, 1/4, 1/8 (k=2..14)!
3. Test how D_m transforms the Plucker ratio R = W_{AB} W_{CE} / W_{BC}^2 on REAL-STABLE 4-tuples (A, B, C, E)
   (constructed properly from real-stable polynomials P(t, y) = sum c_i(t) y^i evaluated at geometric y_0 < y_1 < y_2 < y_3!).
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

def deriv_int(p):
    return [r*p[r] for r in range(1, len(p))]

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

Hp = [[deriv_int(H[k][i]) if H[k][i] is not None else None for i in range(k+1)] for k in range(KMAX+1)]

print("=== 1. RAYLEIGH CHECK OF Q_{k, q}(t, y) = sum q^{i(i-1)/2} htilde_{k,i}(t) y^i ON R^2 (INCL. TINY t > 0) ===")
# Note: at t0 > 0 small, check y0 near the roots of Q_{k, q}(t0, y) = 0!!
# At a root y0 of Q(t0, y) = 0, Rayleigh is R(t0, y0) = Q_t(t0, y0) * Q_y(t0, y0) (NO cancellation!!).
# Even better: check Q_t * Q_y - Q * Q_{ty} at dense (t, y) AND specifically at the y-roots of Q(t0, y) = 0!
ts_probe = [mp.mpf(s) for s in ["-2", "-0.99", "-0.5", "-0.05", "-1e-5", "0", "1e-12", "1e-9", "1e-7", "1e-5", "1e-4", "1e-3", "0.01", "0.1", "1"]]
for q_str in ["0.5", "0.25", "0.125", "0.0625"]:
    q = mp.mpf(q_str)
    fails = []
    for k in range(2, 17):
        mp.mp.dps = 100
        # For each t0 in ts_probe, find y-roots of Q(t0, y) via Newton from geometric guesses and check sign of Q_t * Q_y!
        bad_for_k = None
        for tv in ts_probe:
            c0 = [q**(i*(i-1)//2) * eval_mp(H[k][i], tv) for i in range(k+1)]
            c1 = [q**(i*(i-1)//2) * eval_mp(Hp[k][i], tv) for i in range(k+1)]
            # Also directly check Rayleigh at geometric y points y = -alpha * q^{-j} * c0[j]/c0[j+1]
            # First try mp.polyroots for k <= 12 (or Newton for all k when tv >= 0)
            if tv >= 0:
                # For tv >= 0, all c0[i] > 0! Approximate roots: y_j ~ -c0[j-1] / c0[j] (for small q!)
                # Check signs of Q_t * Q_y at Newton-refined roots AND check if polynomial has k real roots
                # Even simpler: check sign changes of Q(tv, y) at midpoints y_mid = -sqrt(c0[j-1]*c0[j+1])/c0[j] (nope, geometric mean of consecutive ratios)
                rts = []
                for j in range(1, k + 1):
                    yj = -c0[j-1] / c0[j]
                    for _ in range(40):
                        v = mp.mpf(0); dv = mp.mpf(0)
                        for idx in range(k, -1, -1):
                            if idx < k: dv = dv * yj + v
                            v = v * yj + c0[idx]
                        if dv == 0: break
                        step = v / dv; yj -= step
                        if abs(step) < mp.mpf('1e-40') * (abs(yj) + 1): break
                    rts.append(yj)
                # Deduplicate distinct roots found
                rts.sort()
                uniq = []
                for r in rts:
                    if not uniq or abs(r - uniq[-1]) > mp.mpf('1e-20') * (abs(r) + 1):
                        uniq.append(r)
                if len(uniq) < k and q <= mp.mpf("0.25"):
                    # Maybe complex roots! Check with polyroots
                    pass
                signs = []
                for yj in uniq:
                    v0 = mp.mpf(0); dv0 = mp.mpf(0); v1 = mp.mpf(0)
                    for idx in range(k, -1, -1):
                        if idx < k: dv0 = dv0 * yj + v0
                        v0 = v0 * yj + c0[idx]; v1 = v1 * yj + c1[idx]
                    prod = dv0 * v1
                    signs.append(1 if prod > 0 else (-1 if prod < 0 else 0))
                if len(set(signs)) > 1 or 0 in signs:
                    bad_for_k = (k, mp.nstr(tv, 3), len(uniq), signs.count(1), signs.count(-1))
                    break
        if bad_for_k:
            fails.append(bad_for_k)
            if len(fails) >= 3: break
    print(f"  Q_{{k, q}} (q={q_str:6s}) Rayleigh fails on t>=0 (k=2..16): {fails}", flush=True)

print("\n=== 2. TEST D-LIFTED q-DEFORMED G_{k, q}(t, y) = sum q^{i(i-1)/2} (D_{k+1}^i htilde_{k,i}) y^i ===")
mp.mp.dps = 80
for q_str in ["0.5", "0.25"]:
    q = mp.mpf(q_str)
    fails = []
    for k in range(2, 13):
        m = k + 1
        Gs = []; Gps = []
        for i in range(k+1):
            p = list(H[k][i])
            for _ in range(i): p = apply_D_int(m, p)
            Gs.append(p); Gps.append(deriv_int(p))
        for tv in [mp.mpf("0"), mp.mpf("1e-7"), mp.mpf("1e-3"), mp.mpf("0.1"), mp.mpf("1"), mp.mpf("-0.5"), mp.mpf("-1")]:
            c0 = [q**(i*(i-1)//2) * eval_mp(Gs[i], tv) for i in range(k+1)]
            c1 = [q**(i*(i-1)//2) * eval_mp(Gps[i], tv) for i in range(k+1)]
            # Sample y near -c0[j-1]/c0[j] if all same sign, plus standard grid
            ys = [-abs(c0[j-1]/c0[j]) for j in range(1, k+1) if c0[j] != 0]
            ys += [mp.mpf(s) for s in ["-100", "-10", "-1", "-0.1", "-0.01", "0.1", "1"]]
            pos = 0; neg = 0
            for yv in ys:
                A = mp.mpf(0); At = mp.mpf(0); Ay = mp.mpf(0); Aty = mp.mpf(0)
                yp = mp.mpf(1)
                for idx in range(k+1):
                    A += c0[idx]*yp; At += c1[idx]*yp
                    if idx > 0:
                        ypp = yp / yv if yv != 0 else (mp.mpf(1) if idx == 1 else mp.mpf(0))
                        Ay += idx * c0[idx] * ypp; Aty += idx * c1[idx] * ypp
                    yp *= yv
                ray = At * Ay - A * Aty
                if ray > 0: pos += 1
                elif ray < 0: neg += 1
            if pos > 0 and neg > 0:
                fails.append((k, mp.nstr(tv, 3), pos, neg)); break
    print(f"  G_{{k, q}} (q={q_str:4s}) Rayleigh fails (k=2..12): {fails[:4]}", flush=True)

print("\n=== 3. HOW D_m TRANSFORMS PLUCKER RATIO R = W_{AB} W_{CE} / W_{BC}^2 ON GENUINE REAL-STABLE QUADS ===")
# Construct genuine real-stable 4-tuples (A, B, C, E) by taking P(t, y) = prod_{r=1}^d (t + a_r + b_r y) (a_r, b_r > 0)
# and evaluating at y = y0 < y1 < y2 < y3 (WAIT: by real stability of P, ANY 3 evaluations at y0 < y1 < y2 give
# a real-stable triple? Let's check! Or take P(t, y) = A + 2By + Cy^2? No, for a 4-tuple, take P(t, y) = prod (t + a_r + b_r y)
# and take [y^0], [y^1], [y^2], [y^3] coefficients, OR evaluate at 4 points!)
def apply_D_flt(m, p):
    res = np.zeros(len(p) + 1)
    for r, c in enumerate(p):
        res[r] += c * (r + 1); res[r + 1] += c * (r + m)
    return res

def eval_W(f, g, ts):
    fp = np.arange(1, len(f)) * f[1:]; gp = np.arange(1, len(g)) * g[1:]
    return np.polyval(fp[::-1], ts)*np.polyval(g[::-1], ts) - np.polyval(f[::-1], ts)*np.polyval(gp[::-1], ts)

rng = np.random.default_rng(123)
ts_grid = np.concatenate([-10**np.linspace(2, -4, 40), np.linspace(-0.99, -0.01, 100), [0.0], 10**np.linspace(-4, 2, 40)])
valid = 0; inc_cnt = 0; max_amp = 0.0; max_out_when_in_le1 = 0.0
for _ in range(200):
    d = 4
    # Build real-stable P(t, y) = prod_{r=1}^d (1 + u_r t + v_r y) with roots in (-1, 0) at y=0 => u_r > 1, v_r > 0
    us = rng.uniform(1.1, 10.0, d); vs = rng.uniform(0.1, 5.0, d)
    # Wait: evaluate P(t, y) at 4 increasing y values y = 0, s, 2s, 3s? Or geometric?
    # First check [y^j] coeffs (they have decreasing degree unless we multiply by (1+t)^j? No, use evaluations at 4 y's!)
    ys_eval = np.sort(rng.uniform(0.0, 2.0, 4))
    polys = []
    for yv in ys_eval:
        p = np.array([1.0])
        for r in range(d):
            p = np.convolve(p, [1.0 + vs[r]*yv, us[r]])
        polys.append(p)
    A, B, C, E = polys
    wab = eval_W(A, B, ts_grid); wbc = eval_W(B, C, ts_grid); wce = eval_W(C, E, ts_grid)
    wac = eval_W(A, C, ts_grid); wbe = eval_W(B, E, ts_grid)
    if np.min(4*wab*wbc - wac**2) < -1e-9 or np.min(4*wbc*wce - wbe**2) < -1e-9:
        continue
    if np.min(wbc) <= 0 or np.min(wab) <= 0 or np.min(wce) <= 0:
        continue
    valid += 1
    rho_in = float(np.max((wab * wce) / (wbc**2)))
    m = float(10 ** rng.uniform(0.1, 1.5)) + 1.0
    DA, DB, DC, DE = apply_D_flt(m, A), apply_D_flt(m, B), apply_D_flt(m, C), apply_D_flt(m, E)
    dwab = eval_W(DA, DB, ts_grid); dwbc = eval_W(DB, DC, ts_grid); dwce = eval_W(DC, DE, ts_grid)
    rho_out = float(np.max((dwab * dwce) / (dwbc**2)))
    if rho_out > rho_in + 1e-9: inc_cnt += 1
    if rho_out / rho_in > max_amp: max_amp = rho_out / rho_in
    if rho_in <= 1.0 and rho_out > max_out_when_in_le1: max_out_when_in_le1 = rho_out
print(f"  Valid real-stable quads = {valid}: D_m increased max(rho) in {inc_cnt}/{valid}, max amp={max_amp:.4f}, max rho_out(when rho_in<=1)={max_out_when_in_le1:.4f}")
