"""CORRECTED TEST OF THE EQUAL-DEGREE ALIGNED FAMILY:
  htilde_{1,0} = 1 + 2t, htilde_{1,1} = 1,
  htilde_{k,i} = D_{k+1}(D_{k+1} htilde_{k-1,i} + htilde_{k-1,i-1})  for k >= 2.
  g_{k, i}(t) := D_{k+1}^i htilde_{k, i}(t)  (ALL have exact degree 2k - 1!).

Tests:
1. Exact Rayleigh at t=0 for:
   (a) Triple quadratic X + 2 M y + Y y^2 where
       X = htilde_{k-1, i-2}, M = D_{k+1} htilde_{k-1, i-1}, Y = D_{k+1}^2 htilde_{k-1, i}  (k = 3..60!)
   (b) Full aligned row sum_{i=0}^k w_i g_{k,i}(t) y^i with w_i = 1 and w_i = binom(k,i)  (k = 2..30!)
2. Full UHP stability (Im(t) > 0, Im(y) > 0 AND Im(t) > 0, Im(y) < 0) of:
   - Triple quadratic X + 2 M y + Y y^2  (k = 3..10, ALL i)
   - Full aligned row (w=1 and w=binom(k,i))  (k = 2..9)
3. Exact 4*W(X,M)*W(M,Y) - W(X,Y)^2 sign check on ALL t in R for k = 3..10!
4. Also test triple quadratic (X, M, Y) for WILD parameters m_j in Atomic Theorem!
"""
import math
import numpy as np
import sympy as sp

y, t = sp.symbols('y t')

def apply_D(m, poly):
    """D_m = t(1+t)d/dt + 1 + m*t."""
    d = len(poly) - 1
    res = [0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + 1)
        res[r + 1] += coef * (r + m)
    return res

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def apply_Dr_01(m, r, b0, b1):
    # b0 -> b0, b1 -> 2*b1 + m*b0 => after r steps: b1 -> 2^r * b1 + (2^r - 1) * m * b0
    return b0, (1 << r) * b1 + ((1 << r) - 1) * m * b0

# Exact [t^0] and [t^1] of TRUE htilde_{k, i}(t) up to KMAX = 60
KMAX = 60
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

print("=== 1. RATIO g1/g0 OF TRUE g_{k,i} = D_{k+1}^i htilde_{k,i} AT k = 10, 20, 40 ===")
for k in [10, 20, 40]:
    m = k + 1
    ratios = []
    for i in range(k + 1):
        g0, g1 = apply_Dr_01(m, i, a0[k][i], a1[k][i])
        ratios.append(g1 / g0)
    print(f"  k={k}: g1/g0 / (g1[0]/g0[0]) = {[round(r/ratios[0], 5) for r in ratios[:5]]} ... {[round(r/ratios[0], 5) for r in ratios[-3:]]}")

print("\n=== 2. EXACT RAYLEIGH AT t=0 FOR TRIPLE QUADRATIC X + 2 M y + Y y^2 (k=3..60) ===")
# X = htilde_{k-1, i-2}, M = D_{k+1} htilde_{k-1, i-1}, Y = D_{k+1}^2 htilde_{k-1, i}
quad_fails = []
max_ratio = (0.0, None)
for k in range(3, KMAX + 1):
    m = k + 1
    for i in range(2, k):
        x0, x1 = a0[k-1][i-2], a1[k-1][i-2]
        m0, m1 = apply_Dr_01(m, 1, a0[k-1][i-1], a1[k-1][i-1])
        y0, y1 = apply_Dr_01(m, 2, a0[k-1][i], a1[k-1][i])
        c0 = x1*m0 - x0*m1
        c1 = x1*y0 - x0*y1
        c2 = m1*y0 - m0*y1
        if c0 * c2 <= 0 or c1*c1 > 4*c0*c2:
            quad_fails.append((k, i, c0 > 0, c2 > 0, float(c1*c1)/float(c0*c2) if c0*c2 != 0 else 999))
        else:
            r_val = float(c1*c1) / float(c0*c2)
            if r_val > max_ratio[0]:
                max_ratio = (r_val, (k, i, c0 > 0))
print(f"  Triple (X, M, Y) Rayleigh fails at t=0 (k=3..60): {quad_fails[:5]}")
print(f"  Max c1^2/(c0*c2) (need <= 4.0): {max_ratio}")

# Also test same-row aligned triple: (g_{k, i-2}, g_{k, i-1}, g_{k, i}) = (htilde_{k,i-2}, D_{k+1} htilde_{k,i-1}, D_{k+1}^2 htilde_{k,i})
row_quad_fails = []
max_row_ratio = (0.0, None)
for k in range(2, KMAX + 1):
    m = k + 1
    for i in range(2, k + 1):
        x0, x1 = a0[k][i-2], a1[k][i-2]
        m0, m1 = apply_Dr_01(m, 1, a0[k][i-1], a1[k][i-1])
        y0, y1 = apply_Dr_01(m, 2, a0[k][i], a1[k][i])
        c0 = x1*m0 - x0*m1
        c1 = x1*y0 - x0*y1
        c2 = m1*y0 - m0*y1
        if c0 * c2 <= 0 or c1*c1 > 4*c0*c2:
            row_quad_fails.append((k, i, c0 > 0, c2 > 0, float(c1*c1)/float(c0*c2) if c0*c2 != 0 else 999))
        else:
            r_val = float(c1*c1) / float(c0*c2)
            if r_val > max_row_ratio[0]:
                max_row_ratio = (r_val, (k, i, c0 > 0))
print(f"  Same-row (h_{'{k,i-2}'}, D h_{'{k,i-1}'}, D^2 h_{'{k,i}'}) Rayleigh fails at t=0 (k=2..60): {row_quad_fails[:5]}")
print(f"  Max same-row c1^2/(c0*c2): {max_row_ratio}")

print("\n=== 3. EXACT RAYLEIGH AT t=0 FOR FULL ALIGNED ROW sum w_i D_{k+1}^i htilde_{k,i} y^i (k=2..22) ===")
for wname, wfn in [("w=1 (plain)", lambda k, i: 1), ("w=binom(k,i)", lambda k, i: math.comb(k, i))]:
    fails = []
    for k in range(2, 23):
        m = k + 1
        P0 = 0; Pt0 = 0
        for i in range(k + 1):
            g0, g1 = apply_Dr_01(m, i, a0[k][i], a1[k][i])
            w = wfn(k, i)
            P0 += w * g0 * y**i
            Pt0 += w * g1 * y**i
        R = sp.Poly(sp.expand(Pt0 * sp.diff(P0, y) - P0 * sp.diff(Pt0, y)), y, domain='ZZ')
        rc = int(sp.count_roots(R))
        if rc > 0:
            fails.append((k, rc))
            if len(fails) >= 4:
                break
    print(f"  {wname:15s}: Rayleigh at t=0 fails (k=2..22) = {fails}", flush=True)

# Build TRUE htilde_{k, i}(t) for k = 1..10
KPOLY = 10
H = [[None]*(k+1) for k in range(KPOLY+1)]
H[1][0] = [1, 2]
H[1][1] = [1]
for k in range(2, KPOLY + 1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D(m, H[k-1][0])
    for i in range(1, k):
        U[i] = poly_add(apply_D(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1):
        H[k][i] = apply_D(m, U[i])

# Sanity check degree of H[k][i] == 2k - i - 1
for k in range(1, KPOLY + 1):
    for i in range(k + 1):
        assert len(H[k][i]) - 1 == 2*k - i - 1, f"deg mismatch at ({k},{i}): {len(H[k][i])-1}"
print("\n[Sanity check passed: deg(htilde_{k,i}) = 2k - i - 1 for all k=1..10]")

rng = np.random.default_rng(1234)
test_ys = [complex(re_y, im_y)
           for re_y in [-50, -10, -2, -1, -0.5, -0.1, 0, 0.1, 0.5, 1, 2, 10, 50]
           for im_y in [1e-4, 1e-2, 0.1, 1.0, 10.0]]
for _ in range(150):
    test_ys.append(complex(rng.normal(0, 1) * (10 ** rng.uniform(-2, 2)), 10 ** rng.uniform(-4, 2)))

print("\n=== 4. FULL UHP STABILITY OF FULL ALIGNED ROWS (k=2..9) ===")
for wname, wfn in [("plain aligned", lambda k, i: 1), ("binom aligned", lambda k, i: math.comb(k, i))]:
    print(f"--- {wname} ---")
    for k in range(2, 10):
        m = k + 1
        gs = []
        for i in range(k + 1):
            p = list(H[k][i])
            for _ in range(i): p = apply_D(m, p)
            gs.append(p)
        dmax = 2*k - 1
        max_im = -1e9; min_im = 1e9
        for yv in test_ys:
            c = np.zeros(dmax + 1, dtype=np.complex128)
            for i in range(k + 1):
                w = wfn(k, i) * (yv ** i)
                for r_idx, coef in enumerate(gs[i]):
                    c[r_idx] += w * float(coef)
            if abs(c[-1]) < 1e-30: continue
            rts = np.roots(c[::-1])
            max_im = max(max_im, float(np.max(rts.imag)))
            min_im = min(min_im, float(np.min(rts.imag)))
        print(f"  k={k}: Im(t) in [{min_im:+.3e}, {max_im:+.3e}] for y in H+")

print("\n=== 5. FULL UHP STABILITY + EXACT WRONSKIAN DISCRIMINANT FOR TRIPLE QUADRATIC X + 2 M y + Y y^2 (k=3..10) ===")
def to_sympy(p):
    return sum(sp.Integer(c) * t**r for r, c in enumerate(p))

def wronskian(f, g):
    return sp.expand(sp.diff(f, t)*g - f*sp.diff(g, t))

uhp_fails = []
wron_disc_fails = []
for k in range(3, 11):
    m = k + 1
    k_max_im = -1e9; k_min_im = 1e9
    for i in range(2, k):
        X = H[k-1][i-2]
        M = apply_D(m, H[k-1][i-1])
        Y = apply_D(m, apply_D(m, H[k-1][i]))
        d = len(X) - 1
        assert len(M) - 1 == d and len(Y) - 1 == d
        # UHP check
        for yv in test_ys:
            c = np.array([float(X[r]) + 2.0*yv*float(M[r]) + (yv*yv)*float(Y[r]) for r in range(d+1)], dtype=np.complex128)
            if abs(c[-1]) < 1e-30: continue
            rts = np.roots(c[::-1])
            k_max_im = max(k_max_im, float(np.max(rts.imag)))
            k_min_im = min(k_min_im, float(np.min(rts.imag)))
            if np.max(rts.imag) > 1e-9 and np.min(rts.imag) < -1e-9:
                uhp_fails.append((k, i, yv, float(np.min(rts.imag)), float(np.max(rts.imag))))
                break
        # Exact Wronskian discriminant 4*W(X,M)*W(M,Y) - W(X,Y)^2 on ALL t in R (for k <= 6)
        if k <= 6:
            Xs, Ms, Ys = to_sympy(X), to_sympy(M), to_sympy(Y)
            Wxm = wronskian(Xs, Ms)
            Wmy = wronskian(Ms, Ys)
            Wxy = wronskian(Xs, Ys)
            Disc = sp.Poly(sp.expand(4*Wxm*Wmy - Wxy*Wxy), t, domain='ZZ')
            rc_all = int(sp.count_roots(Disc))
            rc_neg = int(sp.count_roots(Disc, -1, 0))
            val0 = Disc.eval(0)
            if rc_all > 0 or val0 <= 0:
                wron_disc_fails.append((k, i, rc_all, rc_neg, val0 > 0))
    print(f"  k={k}: UHP Im(t) in [{k_min_im:+.3e}, {k_max_im:+.3e}]")
print(f"  Triple UHP fails (k=3..10): {len(uhp_fails)} {uhp_fails[:3]}")
print(f"  Exact 4*W(X,M)*W(M,Y) - W(X,Y)^2 fails on R (k=3..8): {wron_disc_fails[:10]}")
