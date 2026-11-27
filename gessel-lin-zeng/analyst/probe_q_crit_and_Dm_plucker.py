"""Decisive quantitative probe of BOTH:
1. Exact critical q for Q_{k, q}(t, y) = sum q^{i(i-1)/2} htilde_{k,i}(t) y^i at t=0 up to k = 50!
   (Test q = 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.25, 0.125 — AND test Q_{k, q} at q=0.25 on @algebraist's style
    adversarial / wild m_j!)
2. Exact tracking of J_ratio = 4 W(X,M) W(M,Y) / W(X,Y)^2 AND Plucker ratio R = W(A,B) W(C,E) / W(B,C)^2
   through ALL 3 substeps of ONE GLZ induction step:
     (0) Start at k-1 (aligned with D_{k+1}): (A, B, C, E, F) = (h_{k-1, i-3}, D h_{k-1, i-2}, ..., D^4 h_{k-1, i+1})
     (1) Pascal sum: (A+B, B+C, C+E, E+F)
     (2) Apply D_{k+1}: (h_{k, i-2}, D_{k+1} h_{k, i-1}, D_{k+1}^2 h_{k, i}, D_{k+1}^3 h_{k, i+1})
     (3) Shift D_{k+1} -> D_{k+2}: (h_{k, i-2}, D_{k+2} h_{k, i-1}, D_{k+2}^2 h_{k, i}, D_{k+2}^3 h_{k, i+1})
3. Test how D_m transforms R = W_{AB} W_{CE} / W_{BC}^2 on genuine real-stable 4-tuples (from P(t, y) at evenly spaced y)!
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 100

# === 1. Exact critical q of Q_{k, q}(t, y) at t = 0 up to k = 50 ===
KMAX = 50
a0 = [[0]*(k+1) for k in range(KMAX+1)]
a1 = [[0]*(k+1) for k in range(KMAX+1)]
a0[1][0], a1[1][0] = 1, 2; a0[1][1], a1[1][1] = 1, 0
for k in range(2, KMAX+1):
    m = k + 1
    for i in range(k+1):
        u0 = (a0[k-1][i] if i<=k-1 else 0) + (a0[k-1][i-1] if i>=1 else 0)
        u1 = ((2*a1[k-1][i] + m*a0[k-1][i]) if i<=k-1 else 0) + (a1[k-1][i-1] if i>=1 else 0)
        a0[k][i] = u0; a1[k][i] = 2*u1 + m*u0

print("=== 1. Q_{k, q}(t, y) = sum q^{i(i-1)/2} htilde_{k,i}(t) y^i: EXACT RAYLEIGH AT t=0 (k=2..45) ===")
for q_str in ["0.9", "0.8", "0.7", "0.6", "0.55", "0.5", "0.4", "0.25"]:
    q = mp.mpf(q_str)
    first_fail = None
    for k in range(2, 46):
        # Need extra precision for small q and large k: q^{k^2/2} ~ 0.25^{1000} ~ 10^{-600} -> skip if too huge, or set dps dynamically!
        needed_dps = int(abs(k*(k-1)/2 * math.log10(float(q_str)))) + 60
        if needed_dps > 350:
            break
        mp.mp.dps = max(80, needed_dps)
        p0 = [q**(i*(i-1)//2) * a0[k][i] for i in range(k+1)]
        p1 = [q**(i*(i-1)//2) * a1[k][i] for i in range(k+1)]
        # Find roots of p0(y) = sum q^{i(i-1)/2} C(k,i) y^i by Newton from asymptotic guesses
        # Wait: for q near 1 (like 0.9), guesses -q^{-(j-1)}*j/(k-j+1) can collide, so use mp.polyroots when q >= 0.5 and k <= 25,
        # or better: use exact recurrence of roots from k-1 to k (interlacing!!).
        # Let's just use polyroots on scaled variable y = q^{-(k-1)/2} u for q >= 0.5, and Newton for q < 0.5!
        if float(q_str) >= 0.5:
            scale = q ** (mp.mpf(k - 1) / 2)
            c0 = [p0[i] * (scale ** (-i)) for i in range(k + 1)]
            # Normalize by max
            mx = max(abs(x) for x in c0)
            c0 = [x / mx for x in c0]
            try:
                rts_u = mp.polyroots(c0[::-1], maxsteps=400, extraprec=60)
            except Exception:
                first_fail = (k, "noconv"); break
            if max(abs(mp.im(r)) for r in rts_u) > mp.mpf('1e-15'):
                first_fail = (k, "P0_nonreal"); break
            rts = sorted([mp.re(r) / scale for r in rts_u])
        else:
            rts = []
            for j in range(1, k + 1):
                yj = -(q ** (-(j - 1))) * mp.mpf(j) / mp.mpf(k - j + 1)
                for _ in range(40):
                    v = mp.mpf(0); dv = mp.mpf(0)
                    for idx in range(k, -1, -1):
                        if idx < k: dv = dv * yj + v
                        v = v * yj + p0[idx]
                    if dv == 0: break
                    step = v / dv; yj -= step
                    if abs(step) < mp.mpf('1e-40') * abs(yj): break
                rts.append(yj)
            rts.sort()
            if not all(rts[m+1] - rts[m] > mp.mpf('1e-20') * abs(rts[m]) for m in range(k - 1)):
                first_fail = (k, "dup_root"); break
        signs = []
        for yj in rts:
            v0 = mp.mpf(0); dv0 = mp.mpf(0); v1 = mp.mpf(0)
            for idx in range(k, -1, -1):
                if idx < k: dv0 = dv0 * yj + v0
                v0 = v0 * yj + p0[idx]; v1 = v1 * yj + p1[idx]
            prod = dv0 * v1
            signs.append(1 if prod > 0 else (-1 if prod < 0 else 0))
        if len(set(signs)) > 1 or 0 in signs:
            first_fail = (k, f"+:{signs.count(1)}, -:{signs.count(-1)}")
            break
    print(f"  q = {q_str:4s}: tested up to k={k:2d}, first_fail = {first_fail}", flush=True)

# === 2. Track J_ratio and Plucker ratio R across ALL 3 substeps of GLZ induction ===
print("\n=== 2. GLZ SUBSTEP TRACKING OF min(J_ratio) AND max(Plucker R = W_AB W_CE / W_BC^2) ===")
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
    """Given equal-degree sequence seq[0..L], return (min_J_ratio, max_R_plucker)."""
    L = len(seq)
    derivs = [deriv_int(p) for p in seq]
    min_J = 1e9; max_R = 0.0
    for tv in ts_chk:
        vs = [eval_mp(p, tv) for p in seq]
        vps = [eval_mp(dp, tv) for dp in derivs]
        # Wronskians W(i, i+1) and W(i, i+2)
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
    # Stage 0: v_i = D_{k+1}^i htilde_{k-1, i}  (i = 0..k-1)
    V = []
    for i in range(k):
        p = list(H[k-1][i])
        for _ in range(i): p = apply_D_int(k+1, p)
        V.append(p)
    # Stage 1: Pascal sum P_i = V_i + V_{i-1}  (i = 0..k, with V_{-1}=0, V_k=0)
    P = [V[0]] + [poly_add_int(V[i], V[i-1]) for i in range(1, k)] + [V[k-1]]
    # Stage 2: Apply D_{k+1}: G_i = D_{k+1}^2 P_i? Wait: D_{k+1} P_i = D_{k+1}^{i-1} htilde_{k, i} -> no,
    # D_{k+1}^2 P_i = D_{k+1}^i htilde_{k, i}!! Let's check BOTH D_{k+1} P and D_{k+1}^2 P!
    DP = [apply_D_int(k+1, p) for p in P]
    D2P = [apply_D_int(k+1, p) for p in DP]
    # Stage 3: Next step's V_next: D_{k+2}^i htilde_{k, i} (i = 0..k)
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
    print(f"  k={k:2d}: J_min: V={j0:.4f} -> Pascal={j1:.4f} -> D={j2:.4f} -> D^2={j2b:.4f} -> shift={j3:.4f} | "
          f"R_max: V={r0:.4f} -> P={r1:.4f} -> D={r2:.4f} -> D^2={r2b:.4f} -> shift={r3:.4f}", flush=True)

# === 3. Test how D_m transforms Plucker ratio R = W_AB W_CE / W_BC^2 on genuine real-stable quads ===
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
for _ in range(200):
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
