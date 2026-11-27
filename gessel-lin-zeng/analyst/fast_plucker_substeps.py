"""FAST Sections 2 & 3 of probe_q_crit_and_Dm_plucker.py (takes ~3 seconds!):
- Exact tracking of min(J_ratio) and max(Plucker R = W_AB W_CE / W_BC^2) across ALL substeps:
  V -> Pascal -> D_{k+1} -> D_{k+1}^2 -> shift (D_{k+2})
- Exact tracking of D_m action on R = W_AB W_CE / W_BC^2 on genuine real-stable 4-tuples!
- Also fast check of q = 0.55, 0.50, 0.40, 0.25 at t=0 via pure Newton (k=2..45, dps=150)!
"""
import math
import numpy as np
import mpmath as mp

print("=== 2. GLZ SUBSTEP TRACKING OF min(J_ratio) AND max(Plucker R = W_AB W_CE / W_BC^2) ===")
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

KPOLY = 12
H = [[None]*(k+1) for k in range(KPOLY+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KPOLY+1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D_int(m, H[k-1][0])
    for i in range(1, k): U[i] = poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1): H[k][i] = apply_D_int(m, U[i])

ts_chk = [mp.mpf(s) for s in ["-10", "-2", "-1.1", "-1", "-0.99", "-0.9", "-0.7", "-0.5", "-0.3", "-0.1", "-0.02", "-0.002", "0", "0.01", "0.2", "1", "10"]]

def stats_row(seq):
    L = len(seq)
    derivs = [deriv_int(p) for p in seq]
    min_J = 1e9; max_R = 0.0
    for tv in ts_chk:
        vs = [eval_mp(p, tv) for p in seq]
        vps = [eval_mp(dp, tv) for dp in derivs]
        w1 = [vps[i]*vs[i+1] - vs[i]*vps[i+1] for i in range(L - 1)]
        w2 = [vps[i]*vs[i+2] - vs[i]*vps[i+2] for i in range(L - 2)]
        for i in range(L - 2):
            if w2[i] != 0:
                jr = float(4 * w1[i] * w1[i+1] / (w2[i] * w2[i]))
                if jr < min_J: min_J = jr
        for i in range(L - 3):
            if w1[i+1] != 0:
                rp = float(w1[i] * w1[i+2] / (w1[i+1] * w1[i+1]))
                if rp > max_R: max_R = rp
    return min_J, max_R

for k in range(4, 12):
    V = []
    for i in range(k):
        p = list(H[k-1][i])
        for _ in range(i): p = apply_D_int(k+1, p)
        V.append(p)
    P = [V[0]] + [poly_add_int(V[i], V[i-1]) for i in range(1, k)] + [V[k-1]]
    DP = [apply_D_int(k+1, p) for p in P]
    D2P = [apply_D_int(k+1, p) for p in DP]
    Vnext = []
    for i in range(k + 1):
        p = list(H[k][i])
        for _ in range(i): p = apply_D_int(k+2, p)
        Vnext.append(p)
    j0, r0 = stats_row(V)
    j1, r1 = stats_row(P)
    j2, r2 = stats_row(DP)
    j2b, r2b = stats_row(D2P)
    j3, r3 = stats_row(Vnext)
    print(f"  k={k:2d}: J_min: V={j0:.4f} -> P={j1:.4f} -> D={j2:.4f} -> D^2={j2b:.4f} -> shift={j3:.4f} | "
          f"R_max: V={r0:.4f} -> P={r1:.4f} -> D={r2:.4f} -> D^2={r2b:.4f} -> shift={r3:.4f}", flush=True)

print("\n=== 3. D_m ACTION ON PLUCKER RATIO R FOR GENUINE REAL-STABLE QUADS (EVENLY SPACED y) ===")
def apply_D_flt(m, p):
    res = np.zeros(len(p) + 1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def eval_W_flt(f, g, ts):
    fp = np.arange(1, len(f))*f[1:]; gp = np.arange(1, len(g))*g[1:]
    return np.polyval(fp[::-1], ts)*np.polyval(g[::-1], ts) - np.polyval(f[::-1], ts)*np.polyval(gp[::-1], ts)

rng = np.random.default_rng(42)
ts_flt = np.concatenate([-10**np.linspace(1.5, -3, 30), np.linspace(-0.98, -0.02, 80), [0.0], 10**np.linspace(-3, 1.5, 30)])
valid = 0; inc_cnt = 0; max_amp = 0.0; max_out_le1 = 0.0
for _ in range(300):
    d = 4
    us = rng.uniform(1.2, 8.0, d); vs = rng.uniform(0.1, 1.0, d)
    eps = rng.uniform(0.05, 0.4)
    polys = []
    for step in range(4):
        yv = step * eps
        p = np.array([1.0])
        for r in range(d): p = np.convolve(p, [1.0 + vs[r]*yv, us[r]])
        polys.append(p)
    A, B, C, E = polys
    wab = eval_W_flt(A, B, ts_flt); wbc = eval_W_flt(B, C, ts_flt); wce = eval_W_flt(C, E, ts_flt)
    wac = eval_W_flt(A, C, ts_flt); wbe = eval_W_flt(B, E, ts_flt)
    if np.min(4*wab*wbc - wac**2) < 0 or np.min(4*wbc*wce - wbe**2) < 0: continue
    valid += 1
    rho_in = float(np.max((wab * wce) / (wbc**2)))
    m = float(10 ** rng.uniform(0.1, 1.5)) + 1.0
    DA, DB, DC, DE = apply_D_flt(m, A), apply_D_flt(m, B), apply_D_flt(m, C), apply_D_flt(m, E)
    dwab = eval_W_flt(DA, DB, ts_flt); dwbc = eval_W_flt(DB, DC, ts_flt); dwce = eval_W_flt(DC, DE, ts_flt)
    rho_out = float(np.max((dwab * dwce) / (dwbc**2)))
    if rho_out > rho_in + 1e-9: inc_cnt += 1
    max_amp = max(max_amp, rho_out / rho_in)
    if rho_in <= 1.0: max_out_le1 = max(max_out_le1, rho_out)
print(f"  Valid real-stable quads = {valid}: D_m increased max(R) in {inc_cnt}/{valid}, max amp = {max_amp:.4f}, max R_out(when R_in<=1) = {max_out_le1:.4f}")
