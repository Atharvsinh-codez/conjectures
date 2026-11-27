"""EXACT mpmath verification of:
1. Did [1.0001]*8 + [1e6]*2 ACTUALLY break Turan/interlacing/real-rootedness of prod_{j=1}^k D_{m_j}(D_{m_j}+y)(1),
   or was float64 losing precision at t = -0.999?!
2. Check if W^{(3)}_{k, i}(t, y) = g_{k,i-3} + 3 g_{k,i-2} y + 3 g_{k,i-1} y^2 + g_{k,i} y^3 is real-stable
   for GLZ k = 3..12 (and Rayleigh at t=0 for k = 3..60)!
3. Check the universal cubic Wronskian conjecture W(p0,p1)*W(p2,p3) <= W(p1,p2)^2 for ALL real-stable cubics
   (including degenerate / extreme boundary cases)!
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 80

def apply_D_mp(m, poly):
    res = [mp.mpf(0)]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add_mp(a, b):
    res = [mp.mpf(0)]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def step_double_D_mp(polys, m):
    k = len(polys)
    U = [[mp.mpf(0)] for _ in range(k + 1)]
    for i, p in enumerate(polys):
        U[i] = poly_add_mp(U[i], apply_D_mp(m, p))
        U[i+1] = poly_add_mp(U[i+1], p)
    return [apply_D_mp(m, u) for u in U]

def eval_mp(p, tv):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * tv + c
    return acc

print("=== 1. EXACT MPMATH (80 dps) CHECK OF [1.0001]*8 + [1e6]*2 AND [1e6]*2 + [1.0001]*8 ===")
for label, ms in [
    ("[1.0001]*8 + [1e6]*2", [mp.mpf("1.0001")]*8 + [mp.mpf("1e6")]*2),
    ("[1e6]*2 + [1.0001]*8", [mp.mpf("1e6")]*2 + [mp.mpf("1.0001")]*8),
    ("[1.00001]*6 + [1e8]*4", [mp.mpf("1.00001")]*6 + [mp.mpf("1e8")]*4),
    ("[1e8]*4 + [1.00001]*6", [mp.mpf("1e8")]*4 + [mp.mpf("1.00001")]*6),
]:
    polys = [[mp.mpf(1)]]
    for mv in ms: polys = step_double_D_mp(polys, mv)
    k = len(ms)
    ts_neg = [mp.mpf(s) for s in ["-0.99999", "-0.9999", "-0.999", "-0.99", "-0.9", "-0.5", "-0.1", "-0.01", "-1e-5", "-1e-9"]]
    mn = mp.mpf('1e100'); worst = None
    for tv in ts_neg:
        vals = [eval_mp(p, tv) for p in polys]
        for i in range(1, k):
            num = vals[i]**2 - vals[i-1]*vals[i+1]
            den = vals[i]**2 + abs(vals[i-1]*vals[i+1]) + mp.mpf('1e-200')
            r = num / den
            if r < mn: mn = r; worst = (i, float(tv))
    print(f"  {label:26s}: EXACT min_rel_Turan = {mp.nstr(mn, 6)} at {worst}", flush=True)

print("\n=== 2. CUBIC WINDOW g_{k,i-3} + 3 g_{k,i-2} y + 3 g_{k,i-1} y^2 + g_{k,i} y^3 RAYLEIGH AT t=0 (k=3..60) ===")
KMAX = 60
a0 = [[0]*(k+1) for k in range(KMAX+1)]; a1 = [[0]*(k+1) for k in range(KMAX+1)]
a0[1][0], a1[1][0] = 1, 2; a0[1][1], a1[1][1] = 1, 0
for k in range(2, KMAX+1):
    m = k + 1
    for i in range(k+1):
        u0 = (a0[k-1][i] if i<=k-1 else 0) + (a0[k-1][i-1] if i>=1 else 0)
        u1 = ((2*a1[k-1][i] + m*a0[k-1][i]) if i<=k-1 else 0) + (a1[k-1][i-1] if i>=1 else 0)
        a0[k][i] = u0; a1[k][i] = 2*u1 + m*u0

cubic_fails = []
for k in range(3, KMAX+1):
    m = k + 1
    for i in range(3, k + 1):
        # g_j at t=0: (g0, g1) for j = i-3, i-2, i-1, i
        gs0 = []; gs1 = []
        for j in range(i-3, i+1):
            g0 = a0[k][j]
            g1 = (1<<j)*a1[k][j] + ((1<<j)-1)*m*a0[k][j]
            gs0.append(mp.mpf(g0)); gs1.append(mp.mpf(g1))
        # P0(y) = g0[0] + 3 g0[1] y + 3 g0[2] y^2 + g0[3] y^3
        # P1(y) = g1[0] + 3 g1[1] y + 3 g1[2] y^2 + g1[3] y^3
        p0 = [gs0[0], 3*gs0[1], 3*gs0[2], gs0[3]]
        p1 = [gs1[0], 3*gs1[1], 3*gs1[2], gs1[3]]
        s0 = max(p0); s1 = max(p1)
        c0 = [x/s0 for x in p0]; c1 = [x/s1 for x in p1]
        rts = mp.polyroots(c0[::-1], maxsteps=200, extraprec=30)
        if max(abs(mp.im(r)) for r in rts) > mp.mpf('1e-20'):
            cubic_fails.append((k, i, "P0_nonreal")); break
        signs = []
        for yj in rts:
            yj = mp.re(yj)
            v0 = ((c0[3]*yj + c0[2])*yj + c0[1])*yj + c0[0]
            dv0 = (3*c0[3]*yj + 2*c0[2])*yj + c0[1]
            v1 = ((c1[3]*yj + c1[2])*yj + c1[1])*yj + c1[0]
            signs.append(1 if dv0*v1 > 0 else (-1 if dv0*v1 < 0 else 0))
        if len(set(signs)) > 1 or 0 in signs:
            cubic_fails.append((k, i, signs)); break
    if len(cubic_fails) >= 5: break
print(f"  Aligned cubic window Rayleigh at t=0 fails (k=3..60): {cubic_fails}", flush=True)

print("\n=== 3. STRESS-TEST CONJECTURE: W(p0,p1)*W(p2,p3) <= W(p1,p2)^2 FOR ALL REAL-STABLE CUBICS ===")
rng = np.random.default_rng(999)
ts_r = np.linspace(-3, 3, 150)
max_R = 0.0
for _ in range(1000):
    # Extreme PSD matrices (rank-1 + tiny eps, wild eigenvalues 10^{-4}..10^{4})
    U1 = np.linalg.qr(rng.normal(0, 1, (3, 3)))[0]
    U2 = np.linalg.qr(rng.normal(0, 1, (3, 3)))[0]
    e1 = 10 ** rng.uniform(-3, 3, 3)
    e2 = 10 ** rng.uniform(-3, 3, 3)
    M1 = U1 @ np.diag(e1) @ U1.T
    M2 = U2 @ np.diag(e2) @ U2.T
    t_pts = np.array([-2.0, -1.0, 0.0, 1.0]); y_pts = np.array([0.0, 1.0, -1.0, 2.0])
    V_y = np.vander(y_pts, 4, increasing=True); V_t = np.vander(t_pts, 4, increasing=True)
    grid = np.array([[np.linalg.det(np.eye(3) + tv*M1 + yv*M2) for yv in y_pts] for tv in t_pts])
    C_mat = np.linalg.solve(V_t, np.linalg.solve(V_y, grid.T).T)
    p = [C_mat[:, i] for i in range(4)]
    def eval_W_p(f, g):
        fp = np.arange(1, 4)*f[1:]; gp = np.arange(1, 4)*g[1:]
        return np.polyval(fp[::-1], ts_r)*np.polyval(g[::-1], ts_r) - np.polyval(f[::-1], ts_r)*np.polyval(gp[::-1], ts_r)
    w01 = eval_W_p(p[0], p[1]); w12 = eval_W_p(p[1], p[2]); w23 = eval_W_p(p[2], p[3])
    mask = np.abs(w12) > 1e-12 * (np.max(np.abs(w12)) + 1e-30)
    if not np.any(mask): continue
    r = float(np.max((w01[mask] * w23[mask]) / (w12[mask] ** 2)))
    if r > max_R: max_R = r
print(f"  1000 extreme PSD 3x3 cubics: max W(p0,p1)*W(p2,p3) / W(p1,p2)^2 = {max_R:.6f}")
