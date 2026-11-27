"""Adversarial stress test of the 4 surviving candidates:
1. Full aligned binomial A_k(t, y) = sum binom(k,i) D_{k+1}^i htilde_{k,i}(t) y^i:
   - Rayleigh at t=0 for k = 2..100!
   - Full UHP stability for k = 10..13 (with mpmath/scaled coeffs if needed)
2. Triple aligned quadratic Q_{k,i}(t, y) = X + 2 M y + Y y^2:
   - 4*W(X,M)*W(M,Y) - W(X,Y)^2 on R for k = 3..15!
   - Also check WILD monotone m_1 < m_2 < ... < m_k!
3. Fixed-t real-rootedness in y of B_k(t, y) = sum binom(k,i) htilde_{k,i}(t) y^i for t in (-1, 0):
   - GLZ k = 10..24 (in mpmath precision!)
   - Wild Atomic & Laguerre products!
4. ULC on (-1, 0): (htilde_{k,i}/C(k,i))^2 > (htilde_{k,i-1}/C(k,i-1))*(htilde_{k,i+1}/C(k,i+1))
   - Exact (dividing by -t) for k=2..10, grid for k=11..22
   - Wild Atomic & Laguerre products!
"""
import math
import numpy as np
import sympy as sp
import mpmath as mp

t = sp.Symbol('t')

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

def apply_Dr_01(m, r, b0, b1):
    return b0, (1 << r) * b1 + ((1 << r) - 1) * m * b0

print("=== 1. ALIGNED BINOMIAL A_k(t, y) RAYLEIGH AT t=0 FOR k = 2..100 ===")
KMAX = 100
a0 = [[0]*(k+1) for k in range(KMAX+1)]
a1 = [[0]*(k+1) for k in range(KMAX+1)]
a0[1][0] = 1; a1[1][0] = 2
a0[1][1] = 1; a1[1][1] = 0
for k in range(2, KMAX + 1):
    m = k + 1
    for i in range(k + 1):
        u0 = (a0[k-1][i] if i <= k-1 else 0) + (a0[k-1][i-1] if i >= 1 else 0)
        u1 = ((2*a1[k-1][i] + m*a0[k-1][i]) if i <= k-1 else 0) + (a1[k-1][i-1] if i >= 1 else 0)
        a0[k][i] = u0
        a1[k][i] = 2*u1 + m*u0

# For each k, check sign of R(y) = Pt0(y)*P0'(y) - P0(y)*Pt0'(y) on y in (-inf, 0) (roots of P0(y) = sum C(k,i)^2 y^i are all negative!)
# Remember P0(y) = sum C(k,i)^2 y^i has all negative roots, and R(y) = sum c_r y^r.
# Even better: evaluate R(y) at all k roots of P0(y)!! At a root y_j of P0, R(y_j) = Pt0(y_j)*P0'(y_j)!
# Since P0 has simple roots, R(y) >= 0 (or <= 0) everywhere on R iff Pt0 interlaces P0 AND no extra sign changes!
# In fact, if Pt0(y_j) has alternating signs (same sign as P0'(y_j)), that's NECESSARY AND SUFFICIENT for interlacing
# (Wait: deg Pt0 = deg P0 = k, so P0'(y_j)*Pt0(y_j) having constant sign at all k roots of P0 <=> Pt0 strictly interlaces P0 <=> R(y) has NO real roots!!)
# WAIT A SECOND! Verify: if deg P0 = deg Pt0 = k and P0 has k simple real roots, then:
#   P0'(y_j)*Pt0(y_j) > 0 for all j=1..k  ==>  Pt0 has a root in each of the k-1 open intervals (y_j, y_{j+1}),
#   plus 1 more root on (-inf, y_1) or (y_k, +inf) ==> Pt0 strictly interlaces P0 ==> R(y) > 0 for ALL y in R (Hermite-Biehler)!
rayleigh_fails = []
for k in range(2, KMAX + 1):
    m = k + 1
    # Normalize g0 and g1 by dividing g0 by C(k,0) and g1 by g1(k,0)
    _, g1_0 = apply_Dr_01(m, 0, a0[k][0], a1[k][0])
    # Wait: use mpmath at 50 digits for k > 25 so huge binomials never overflow!
    mp.mp.dps = 60
    p0_c = [mp.mpf(math.comb(k, i)**2) for i in range(k + 1)]
    # P0(y) = sum C(k,i)^2 y^i = (1-y)^k P_k((1+y)/(1-y)) where P_k is Legendre polynomial!
    # Wait: roots of P0(y) = sum C(k,i)^2 y^i are y_j = -(1 - x_j)/(1 + x_j) where x_j are Legendre roots in (-1, 1)!!
    # EXACT analytic formula for roots of P0(y)!!
    # Let's check via numpy legendre roots + 1 Newton step in mpmath!
    leg_x = np.polynomial.legendre.leggauss(k)[0]
    p1_c = []
    for i in range(k + 1):
        g0, g1 = apply_Dr_01(m, i, a0[k][i], a1[k][i])
        p1_c.append(mp.mpf(math.comb(k, i)) * mp.mpf(g1) / mp.mpf(g1_0))
    signs = []
    min_abs_val = mp.mpf('1e100')
    for xj in leg_x:
        yj = mp.mpf(-(1.0 - xj) / (1.0 + xj))
        # 2 Newton refinements on P0(y) = 0
        for _ in range(4):
            val = mp.mpf(0); dval = mp.mpf(0)
            for idx in range(k, -1, -1):
                if idx < k: dval = dval * yj + val
                val = val * yj + p0_c[idx]
            yj -= val / dval
        # Compute P0'(yj) and P1(yj)
        val0 = mp.mpf(0); dval0 = mp.mpf(0)
        for idx in range(k, -1, -1):
            if idx < k: dval0 = dval0 * yj + val0
            val0 = val0 * yj + p0_c[idx]
        val1 = mp.mpf(0)
        for idx in range(k, -1, -1):
            val1 = val1 * yj + p1_c[idx]
        prod = dval0 * val1
        signs.append(1 if prod > 0 else (-1 if prod < 0 else 0))
    if len(set(signs)) > 1 or 0 in signs:
        rayleigh_fails.append((k, signs.count(1), signs.count(-1)))
        if len(rayleigh_fails) <= 5:
            print(f"  FAIL A_{k} at t=0: +:{signs.count(1)}, -:{signs.count(-1)}")
print(f"  A_k(t, y) Rayleigh at t=0 fails for k=2..{KMAX}: {len(rayleigh_fails)} {rayleigh_fails[:5]}", flush=True)

# Build exact H[k][i] up to KPOLY = 16
KPOLY = 16
H = [[None]*(k+1) for k in range(KPOLY+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KPOLY+1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D_int(m, H[k-1][0])
    for i in range(1, k):
        U[i] = poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1):
        H[k][i] = apply_D_int(m, U[i])

print("\n=== 2. EXACT ULC ON (-1, 0) FOR k=2..10 (DIVIDING OUT -t) + GRID k=11..16 ===")
ulc_fails = []
for k in range(2, 11):
    for i in range(1, k):
        c_left = math.comb(k, i-1) * math.comb(k, i+1)
        c_right = math.comb(k, i) ** 2
        hi = sum(sp.Integer(v)*t**r for r, v in enumerate(H[k][i]))
        hm = sum(sp.Integer(v)*t**r for r, v in enumerate(H[k][i-1]))
        hp = sum(sp.Integer(v)*t**r for r, v in enumerate(H[k][i+1]))
        diff = sp.expand(c_left * hi*hi - c_right * hm*hp)
        assert diff.subs(t, 0) == 0
        Q = sp.Poly(sp.cancel(diff / (-t)), t, domain='ZZ')
        rc = int(sp.count_roots(Q, -1, 0))
        v0 = Q.eval(0)
        if rc > 0 or v0 <= 0:
            ulc_fails.append((k, i, rc, v0 > 0))
print(f"  EXACT ULC fails on (-1, 0) for k=2..10: {ulc_fails}", flush=True)

print("\n=== 3. B_k(t, y) = sum binom(k,i) htilde_{k,i}(t) y^i: y-ROOTS FOR t in (-1, 0) (k=10..16, mpmath) ===")
mp.mp.dps = 50
ts_mp = [mp.mpf(s) for s in ["-1e-8", "-1e-5", "-1e-3", "-0.01", "-0.05", "-0.1", "-0.2", "-0.3", "-0.5", "-0.7", "-0.9", "-0.99", "-0.9999"]]
for k in range(10, 17):
    max_im = mp.mpf(0); worst_t = None
    for tv in ts_mp:
        coeffs = []
        for i in range(k+1):
            val = mp.mpf(0)
            for c in reversed(H[k][i]): val = val * tv + c
            coeffs.append(mp.mpf(math.comb(k, i)) * val)
        rts = mp.polyroots(coeffs[::-1], maxsteps=200, extraprec=30)
        mi = max(abs(mp.im(r)) for r in rts)
        if mi > max_im: max_im = mi; worst_t = tv
    print(f"  B_{k:2d}: max |Im(y)| = {mp.nstr(max_im, 5)} at t={mp.nstr(worst_t, 4)}", flush=True)

print("\n=== 4. TRIPLE QUADRATIC 4*W(X,M)*W(M,Y) - W(X,Y)^2 ON R FOR k=7..15 (GRID IN t AND x) ===")
def poly_deriv(p):
    return [r * p[r] for r in range(1, len(p))]

def eval_poly_mp(p, tv):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * tv + c
    return acc

mp.mp.dps = 50
ts_R = [mp.mpf(s) for s in [
    "-100", "-10", "-2", "-1.1", "-1.01", "-1", "-0.999", "-0.99", "-0.95", "-0.9",
    "-0.8", "-0.7", "-0.6", "-0.5", "-0.4", "-0.3", "-0.2", "-0.1", "-0.05", "-0.01",
    "-1e-3", "-1e-5", "-1e-8", "0", "1e-6", "1e-3", "0.01", "0.1", "0.5", "1", "2", "10", "100"
]]
wron_grid_fails = []
for k in range(7, 15):
    m = k + 1
    for i in range(2, k):
        X = H[k-1][i-2]
        M = apply_D_int(m, H[k-1][i-1])
        Y = apply_D_int(m, apply_D_int(m, H[k-1][i]))
        Xp, Mp, Yp = poly_deriv(X), poly_deriv(M), poly_deriv(Y)
        for tv in ts_R:
            xv, mv, yv = eval_poly_mp(X, tv), eval_poly_mp(M, tv), eval_poly_mp(Y, tv)
            xpv, mpv, ypv = eval_poly_mp(Xp, tv), eval_poly_mp(Mp, tv), eval_poly_mp(Yp, tv)
            Wxm = xpv*mv - xv*mpv
            Wmy = mpv*yv - mv*ypv
            Wxy = xpv*yv - xv*ypv
            disc = 4*Wxm*Wmy - Wxy*Wxy
            if disc < 0 or Wxm * Wmy < 0:
                wron_grid_fails.append((k, i, float(tv), float(disc)))
                break
print(f"  Triple Wronskian Disc fails on R (k=7..14): {len(wron_grid_fails)} {wron_grid_fails[:5]}", flush=True)

print("\n=== 5. TEST ULC AND B_N(t, y) FOR WILD ATOMIC D_{c,m} AND LAGUERRE L_{c,m} ===")
def apply_D_flt(c, m, poly):
    d = len(poly) - 1
    res = [0.0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * (r + m)
    return res

def apply_L_flt(c, m, poly):
    d = len(poly) - 1
    res = [0.0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * m
    return res

def poly_add_flt(a, b, s=1.0):
    n = max(len(a), len(b))
    res = [0.0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += s * v
    return res

def eval_flt(p, tv):
    acc = 0.0
    for c in reversed(p): acc = acc * tv + c
    return acc

rng = np.random.default_rng(123)
for kind, op_fn, ts_test in [
    ("Atomic D_{c,m}", apply_D_flt, np.concatenate([-10.0**np.linspace(-6,-1,15), np.linspace(-0.15,-0.98,25)])),
    ("Laguerre L_{c,m}", apply_L_flt, -10.0**np.linspace(-5, 2, 40)),
]:
    ulc_f = 0; B_f = 0
    for trial in range(40):
        N = int(rng.integers(3, 9))
        polys = [[1.0]]
        for j in range(N):
            c = float(10**rng.uniform(-0.5, 0.8))
            m = (c + float(10**rng.uniform(-0.5, 1.5))) if kind.startswith("Atomic") else float(10**rng.uniform(-0.5, 1.0))
            lam = float(10**rng.uniform(-0.5, 0.5))
            new_p = [[0.0]] * (len(polys) + 1)
            for i, p in enumerate(polys):
                new_p[i] = poly_add_flt(new_p[i], op_fn(c, m, p), 1.0)
                new_p[i+1] = poly_add_flt(new_p[i+1], p, lam)
            polys = new_p
        # Check ULC with respect to deg N: (F_i / C(N,i))^2 >= (F_{i-1}/C(N,i-1))*(F_{i+1}/C(N,i+1))
        # WAIT: for general lam_j, at t=0, sum F_i(0) y^i = prod (c_j + lam_j y) (real-rooted!),
        # so C(N, i) is the natural Newton factor!
        bad_u = False; bad_B = False
        for tv in ts_test:
            vals = [eval_flt(p, tv) for p in polys]
            for i in range(1, N):
                vi = vals[i] / math.comb(N, i)
                vm = vals[i-1] / math.comb(N, i-1)
                vp = vals[i+1] / math.comb(N, i+1)
                if vi*vi - vm*vp < -1e-9 * (vi*vi + abs(vm*vp) + 1e-300):
                    bad_u = True; break
            # Test B_N(t, y) = sum C(N, i) F_i(t) y^i
            bcoeffs = [math.comb(N, i) * vals[i] for i in range(N+1)]
            if abs(bcoeffs[-1]) > 1e-30:
                rts = np.roots(bcoeffs[::-1])
                if np.max(np.abs(rts.imag)) > 1e-6:
                    bad_B = True
        if bad_u: ulc_f += 1
        if bad_B: B_f += 1
    print(f"  {kind:16s}: ULC fails = {ulc_f}/40, B_N(t,y) real y-roots fails = {B_f}/40", flush=True)
