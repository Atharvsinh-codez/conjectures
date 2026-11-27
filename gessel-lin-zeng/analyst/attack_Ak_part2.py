"""Continuation of attack_Ak_stability.py: Sections 3, 4, 5 (with mpmath ratio to avoid int->float overflow)."""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 60

def apply_D_at_m1(m, b0, b1):
    return (1 - m) * b0, m * b0 - m * b1

def apply_D_at_inf(m, d, cd, cd1):
    return d + 1, (d + m) * cd, (d - 1 + m) * cd1 + (d + 1) * cd

KMAX = 65
b0 = [[0]*(k+1) for k in range(KMAX+1)]
b1 = [[0]*(k+1) for k in range(KMAX+1)]
b0[1][0], b1[1][0] = -1, 2
b0[1][1], b1[1][1] = 1, 0

cd = [[0]*(k+1) for k in range(KMAX+1)]
cd1 = [[0]*(k+1) for k in range(KMAX+1)]
cd[1][0], cd1[1][0] = 2, 1
cd[1][1], cd1[1][1] = 1, 0

a0 = [[0]*(k+1) for k in range(KMAX+1)]
a1 = [[0]*(k+1) for k in range(KMAX+1)]
a0[1][0], a1[1][0] = 1, 2
a0[1][1], a1[1][1] = 1, 0

for k in range(2, KMAX+1):
    m = k + 1
    for i in range(k+1):
        u0, u1 = 0, 0
        if i <= k - 1:
            t0, t1 = apply_D_at_m1(m, b0[k-1][i], b1[k-1][i])
            u0 += t0; u1 += t1
        if i >= 1:
            u0 += b0[k-1][i-1]; u1 += b1[k-1][i-1]
        b0[k][i], b1[k][i] = apply_D_at_m1(m, u0, u1)

        du = 2*k - i - 2
        uc, uc1 = 0, 0
        if i <= k - 1:
            _, tc, tc1 = apply_D_at_inf(m, 2*(k-1)-i-1, cd[k-1][i], cd1[k-1][i])
            uc += tc; uc1 += tc1
        if i >= 1:
            uc += cd[k-1][i-1]; uc1 += cd1[k-1][i-1]
        _, cd[k][i], cd1[k][i] = apply_D_at_inf(m, du, uc, uc1)

        au0 = (a0[k-1][i] if i <= k-1 else 0) + (a0[k-1][i-1] if i >= 1 else 0)
        au1 = ((2*a1[k-1][i] + m*a0[k-1][i]) if i <= k-1 else 0) + (a1[k-1][i-1] if i >= 1 else 0)
        a0[k][i] = au0; a1[k][i] = 2*au1 + m*au0

def safe_ratio(num, den):
    return float(mp.mpf(num) / mp.mpf(den))

print("=== 3. TRIPLE QUADRATIC Q_{k,i} RAYLEIGH AT t = -1 AND t = INF FOR k = 3..60 ===")
q_m1_fails = []; q_inf_fails = []
max_r_m1 = (0.0, None); max_r_inf = (0.0, None)
for k in range(3, 61):
    m = k + 1
    for i in range(2, k):
        x0, x1 = b0[k-1][i-2], b1[k-1][i-2]
        m0, m1 = apply_D_at_m1(m, b0[k-1][i-1], b1[k-1][i-1])
        y0, y1 = apply_D_at_m1(m, *apply_D_at_m1(m, b0[k-1][i], b1[k-1][i]))
        c0 = x1*m0 - x0*m1; c1 = x1*y0 - x0*y1; c2 = m1*y0 - m0*y1
        if c0*c2 <= 0 or c1*c1 > 4*c0*c2:
            q_m1_fails.append((k, i, safe_ratio(c1*c1, c0*c2) if c0*c2 else 999))
        else:
            r = safe_ratio(c1*c1, c0*c2)
            if r > max_r_m1[0]: max_r_m1 = (r, (k, i))

        _, xc, xc1 = 2*(k-1)-(i-2)-1, cd[k-1][i-2], cd1[k-1][i-2]
        dm, mc, mc1 = apply_D_at_inf(m, 2*(k-1)-(i-1)-1, cd[k-1][i-1], cd1[k-1][i-1])
        dy, yc, yc1 = apply_D_at_inf(m, *apply_D_at_inf(m, 2*(k-1)-i-1, cd[k-1][i], cd1[k-1][i]))
        c0 = xc1*mc - xc*mc1; c1 = xc1*yc - xc*yc1; c2 = mc1*yc - mc*yc1
        if c0*c2 <= 0 or c1*c1 > 4*c0*c2:
            q_inf_fails.append((k, i, safe_ratio(c1*c1, c0*c2) if c0*c2 else 999))
        else:
            r = safe_ratio(c1*c1, c0*c2)
            if r > max_r_inf[0]: max_r_inf = (r, (k, i))
print(f"  Triple Q at t=-1 fails (k=3..60): {len(q_m1_fails)} {q_m1_fails[:5]}, max ratio={max_r_m1}")
print(f"  Triple Q at t=inf fails (k=3..60): {len(q_inf_fails)} {q_inf_fails[:5]}, max ratio={max_r_inf}", flush=True)

def check_rayleigh_1var(p0_list, p1_list):
    k = len(p0_list) - 1
    s0 = max(abs(x) for x in p0_list); s1 = max(abs(x) for x in p1_list)
    c0 = [mp.mpf(x)/mp.mpf(s0) for x in p0_list]
    c1 = [mp.mpf(x)/mp.mpf(s1) for x in p1_list]
    rts = mp.polyroots(c0[::-1], maxsteps=300, extraprec=40)
    if max(abs(mp.im(r)) for r in rts) > mp.mpf('1e-20'):
        return False
    signs = []
    for yj in rts:
        yj = mp.re(yj)
        v0 = mp.mpf(0); dv0 = mp.mpf(0); v1 = mp.mpf(0)
        for idx in range(k, -1, -1):
            if idx < k: dv0 = dv0 * yj + v0
            v0 = v0 * yj + c0[idx]
            v1 = v1 * yj + c1[idx]
        prod = dv0 * v1
        signs.append(1 if prod > 0 else (-1 if prod < 0 else 0))
    return (len(set(signs)) == 1 and signs[0] != 0)

print("\n=== 4. GENERAL m-ALIGNED A_{k,m}(t,y) = sum C(k,i) D_m^i htilde_{k,i} y^i AT t=0,-1,inf (k=10, 20) ===")
for k in [10, 20]:
    for m_test in [2, k//2, k, k+1, k+2, k+5, 2*k, 5*k]:
        p0_0, p1_0, p0_m1, p1_m1, p0_inf, p1_inf = [], [], [], [], [], []
        for i in range(k+1):
            w = math.comb(k, i)
            g0, g1 = a0[k][i], (1<<i)*a1[k][i] + ((1<<i)-1)*m_test*a0[k][i]
            p0_0.append(w*g0); p1_0.append(w*g1)
            v0, v1 = b0[k][i], b1[k][i]
            for _ in range(i): v0, v1 = apply_D_at_m1(m_test, v0, v1)
            p0_m1.append(w*v0); p1_m1.append(w*v1)
            d, vc, vc1 = 2*k - i - 1, cd[k][i], cd1[k][i]
            for _ in range(i): d, vc, vc1 = apply_D_at_inf(m_test, d, vc, vc1)
            p0_inf.append(w*vc); p1_inf.append(w*vc1)
        print(f"  k={k:2d}, m={m_test:3d}: t=0:{check_rayleigh_1var(p0_0,p1_0)}, t=-1:{check_rayleigh_1var(p0_m1,p1_m1)}, t=inf:{check_rayleigh_1var(p0_inf,p1_inf)}", flush=True)

print("\n=== 5. TEST A_k(t, y) FOR WILD MONOTONE / DECREASING / RANDOM m_j ===")
def apply_D_int(m, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + 1)
        res[r + 1] += coef * (r + m)
    return res

def poly_add_int(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def test_A_uhp(ms, align_m):
    k = len(ms)
    polys = [[1]]
    for m in ms:
        U = [[0]] * (len(polys) + 1)
        for i, p in enumerate(polys):
            U[i] = poly_add_int(U[i], apply_D_int(m, p))
            U[i+1] = poly_add_int(U[i+1], p)
        polys = [apply_D_int(m, u) for u in U]
    gs = []
    for i, p in enumerate(polys):
        q = list(p)
        for _ in range(i): q = apply_D_int(align_m, q)
        gs.append(q)
    d = len(gs[0]) - 1
    rng = np.random.default_rng(99)
    ys = [complex(rx, ix) for rx in [-10, -2, -1, -0.5, 0, 1, 5] for ix in [1e-3, 0.1, 1.0, 5.0]]
    for _ in range(40):
        ys.append(complex(rng.normal(0, 2), 10**rng.uniform(-3, 1)))
    min_im, max_im = 1e9, -1e9
    for yv in ys:
        c = np.zeros(d+1, dtype=np.complex128)
        for i in range(k+1):
            w = math.comb(k, i) * (yv**i)
            for r, coef in enumerate(gs[i]): c[r] += w * float(coef)
        if abs(c[-1]) < 1e-30: continue
        rts = np.roots(c[::-1])
        min_im = min(min_im, float(np.min(rts.imag)))
        max_im = max(max_im, float(np.max(rts.imag)))
    return min_im, max_im

for label, ms, am in [
    ("const m=3", [3,3,3,3,3], 3),
    ("GLZ k=5 (m=2..6)", [2,3,4,5,6], 6),
    ("squares m=2,5,10,17,26", [2,5,10,17,26], 26),
    ("exp m=2,4,8,16,32", [2,4,8,16,32], 32),
    ("wild jump m=2,3,100,101", [2,3,100,101], 101),
    ("decreasing m=10,8,6,4,2", [10,8,6,4,2], 2),
    ("decreasing align 10", [10,8,6,4,2], 10),
    ("wild m=2,50,3,100", [2,50,3,100], 100),
]:
    mn, mx = test_A_uhp(ms, am)
    stable = (mn >= -1e-9 or mx <= 1e-9)
    print(f"  {label:28s}: Im(t) in [{mn:+.2e}, {mx:+.2e}]  stable={stable}", flush=True)
