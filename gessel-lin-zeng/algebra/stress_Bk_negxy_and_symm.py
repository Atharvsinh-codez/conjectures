"""
Ruthless stress test of surviving bivariate stability candidates:
1. B_k(x, -x*y) = sum_{i=0}^k C(k,i) A_{k,i}(x) (-x)^i y^i for GLZ (k=10..24)
   - Test real-rootedness in z of B_k(x, z) for dense real x in (-100, 100)
   - Test Rayleigh det of B_k(x, -x*y) on real (x, y) grid (especially near x=0, x>0, x<0!)
   - Test UHP x -> max Im(y) and UHP y -> max Im(x)
2. E_k^{binom}(t, y) = sum_{i=0}^{k-1} C(k-1,i) D_{k+1}^i htilde_{k-1,i}(t) (-y)^i for k=16..24
3. Q_k(t, y) = sum_{i=0}^k (1/2)^{C(i,2)} htilde_{k,i}(t) y^i for k=16..24
4. Check roots of Jacobi-Stirling Ehrhart poly p_{k,i}(n) for k=2,3.
"""
import math
import mpmath as mp
mp.mp.dps = 60

def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p
def poly_add(a, b):
    n = max(len(a), len(b))
    return poly_trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def poly_deriv(a):
    return [i * a[i] for i in range(1, len(a))] if len(a) > 1 else [0]
def apply_SL_x(p, L):
    dp = poly_deriv(p)
    return poly_trim(poly_add(poly_add(p, [0]+[L*c for c in p]), poly_add([0]+dp, [0,0]+[-c for c in dp])))
def apply_Dj(p, j):
    dp = poly_deriv(p)
    return poly_trim(poly_add(poly_add([0]+dp, [0,0]+dp), poly_add(p, [0]+[j*c for c in p])))
def poly_eval(p, x):
    acc = mp.mpc(0) if isinstance(x, mp.mpc) else mp.mpf(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc

def compute_A(K):
    A = {0: {0: [1]}}
    for k in range(1, K + 1):
        A[k] = {}
        for i in range(k + 1):
            p1 = apply_SL_x(A[k-1][i], 3*k - 3 - i) if i <= k - 1 else [0]
            p2 = A[k-1][i-1] if i >= 1 else [0]
            A[k][i] = apply_SL_x(poly_add(p1, p2), 3*k - 2 - i)
    return A

def compute_htilde(K):
    h = {1: {0: [1, 2], 1: [1]}}
    for k in range(2, K + 1):
        h[k] = {}
        for i in range(k + 1):
            p1 = apply_Dj(h[k-1][i], k + 1) if i <= k - 1 else [0]
            p2 = h[k-1][i-1] if i >= 1 else [0]
            h[k][i] = apply_Dj(poly_add(p1, p2), k + 1)
    return h

A = compute_A(22)

print("=== 1. RUTHLESS STRESS TEST OF B_k(x, -x*y) ===", flush=True)
# First: for fixed REAL x != 0, B_k(x, -x*y) has all real y-roots <=> B_k(x, z) has all real z-roots!
# Also for real stability of F(x, y) = B_k(x, -x*y):
# For EVERY real x != 0, EVERY real root y_r(x) MUST satisfy dy_r/dx <= 0!
# Equivalently, if F(x, y) = 0, then F_x * F_y - F * F_xy >= 0 (Wait: on F=0, F_x * F_y >= 0 <=> dy/dx = -F_x/F_y <= 0!)
# Let's check BOTH: (a) nonreal z-roots of B_k(x, z) at real x, and (b) sign of dy_r/dx = -F_x/F_y at ALL real roots y_r(x)!
x_grid = [mp.mpf(v) for v in [-100, -10, -3, -1, -0.5, -0.2, -0.1, -0.05, -0.01, -0.001,
                               0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 50, 100]]
for k in [6, 8, 10, 12, 14, 16, 18, 20, 22]:
    max_im_z = mp.mpf(0)
    max_dydx = mp.mpf('-1e100')
    worst_im = None
    worst_deriv = None
    # Precompute derivatives of A_{k,i}(x)
    dA = {i: poly_deriv(A[k][i]) for i in range(k+1)}
    for xv in x_grid:
        # Build F(x, y) = sum_i C(k,i) A_{k,i}(x) (-x)^i y^i
        # and F_x(x, y) = sum_i C(k,i) [A'_{k,i}(x)(-x)^i - i A_{k,i}(x)(-x)^{i-1}] y^i
        cy = []
        cy_x = []
        for i in range(k+1):
            c = math.comb(k, i)
            av = poly_eval(A[k][i], xv)
            dav = poly_eval(dA[i], xv)
            negx_i = (-xv)**i
            cy.append(c * av * negx_i)
            term_x = dav * negx_i + (c * 0 if i == 0 else -i * av * ((-xv)**(i-1)))
            cy_x.append(c * term_x)
        cy_y = [i * cy[i] for i in range(1, k+1)]
        rts = mp.polyroots(list(reversed(cy)), maxsteps=200, extraprec=30)
        for rt in rts:
            if abs(rt.imag) > max_im_z:
                max_im_z = abs(rt.imag)
                worst_im = (float(xv), complex(rt))
            if abs(rt.imag) < mp.mpf('1e-20'):
                yr = rt.real
                fx = poly_eval(cy_x, yr)
                fy = poly_eval(cy_y, yr)
                if fy != 0:
                    dydx = -fx / fy
                    if dydx > max_dydx:
                        max_dydx = dydx
                        worst_deriv = (float(xv), float(yr), float(dydx))
    print(f"  k={k:2d}: max|Im(y)|={mp.nstr(max_im_z, 5)}, max(dy/dx)={mp.nstr(max_dydx, 6)} at {worst_deriv}", flush=True)
    if max_im_z > mp.mpf('1e-15'):
        print(f"     NONREAL ROOT at {worst_im}!", flush=True)

print("\n=== 2. CHECK B_k(x, -x*y) RAYLEIGH & UHP AT MORE x NEAR ROOTS OF A_{k,0} ===", flush=True)
# Check k=10, 14, 18 at x = roots of A_{k,0} and midpoints, plus complex UHP points
for k in [10, 14, 18]:
    rts_A0 = mp.polyroots(list(reversed(A[k][0])), maxsteps=200, extraprec=30)
    rts_A0 = sorted([r.real for r in rts_A0])
    # Test midpoints between consecutive roots of A_{k,0}
    mids = [(rts_A0[m] + rts_A0[m+1])/2 for m in range(len(rts_A0)-1)]
    max_im = mp.mpf(0)
    max_dydx = mp.mpf('-1e100')
    dA = {i: poly_deriv(A[k][i]) for i in range(k+1)}
    for xv in mids:
        cy = [math.comb(k, i) * poly_eval(A[k][i], xv) * ((-xv)**i) for i in range(k+1)]
        cy_x = [math.comb(k, i) * (poly_eval(dA[i], xv)*((-xv)**i) + (0 if i==0 else -i*poly_eval(A[k][i], xv)*((-xv)**(i-1)))) for i in range(k+1)]
        cy_y = [i * cy[i] for i in range(1, k+1)]
        rts = mp.polyroots(list(reversed(cy)), maxsteps=200, extraprec=30)
        for rt in rts:
            if abs(rt.imag) > max_im: max_im = abs(rt.imag)
            if abs(rt.imag) < mp.mpf('1e-20'):
                dydx = -poly_eval(cy_x, rt.real) / poly_eval(cy_y, rt.real)
                if dydx > max_dydx: max_dydx = dydx
    print(f"  k={k} (at all {len(mids)} root midpoints of A_{k,0}): max|Im(y)|={mp.nstr(max_im, 5)}, max(dy/dx)={mp.nstr(max_dydx, 6)}", flush=True)

print("\n=== 3. CHECK B_k(x, -x*y) AT x -> +infty AND x = +0.01..1000 FOR k=24, 30 ===", flush=True)
A30 = compute_A(30)
for k in [24, 28, 30]:
    max_im = mp.mpf(0)
    max_dydx = mp.mpf('-1e100')
    worst = None
    dA = {i: poly_deriv(A30[k][i]) for i in range(k+1)}
    for xv in [mp.mpf(v) for v in [-10, -1, -0.1, -0.01, 0.001, 0.01, 0.1, 1, 10, 100]]:
        cy = [math.comb(k, i) * poly_eval(A30[k][i], xv) * ((-xv)**i) for i in range(k+1)]
        cy_x = [math.comb(k, i) * (poly_eval(dA[i], xv)*((-xv)**i) + (0 if i==0 else -i*poly_eval(A30[k][i], xv)*((-xv)**(i-1)))) for i in range(k+1)]
        cy_y = [i * cy[i] for i in range(1, k+1)]
        rts = mp.polyroots(list(reversed(cy)), maxsteps=300, extraprec=50)
        for rt in rts:
            if abs(rt.imag) > max_im: max_im = abs(rt.imag)
            if abs(rt.imag) < mp.mpf('1e-20'):
                dydx = -poly_eval(cy_x, rt.real) / poly_eval(cy_y, rt.real)
                if dydx > max_dydx:
                    max_dydx = dydx
                    worst = (float(xv), float(rt.real), float(dydx))
    print(f"  k={k}: max|Im(y)|={mp.nstr(max_im, 5)}, max(dy/dx)={mp.nstr(max_dydx, 6)} at {worst}", flush=True)

print("\n=== 4. CHECK binom-E_k(t, y) FOR k=18, 22, 26 (dy/dt <= 0 & real roots) ===", flush=True)
ht26 = compute_htilde(26)
for k in [18, 22, 26]:
    # E_{k,i}(t) = D_{k+1}^i htilde_{k-1, i}(t), i=0..k-1
    E = {}
    dE = {}
    for i in range(k):
        p = ht26[k-1][i]
        for _ in range(i):
            p = apply_Dj(p, k+1)
        E[i] = p
        dE[i] = poly_deriv(p)
    max_im = mp.mpf(0)
    max_dydt = mp.mpf('-1e100')
    worst = None
    for tv in [mp.mpf(v) for v in [-10, -2, -1.01, -0.99, -0.9, -0.5, -0.1, -0.01, 0, 0.01, 0.5, 5]]:
        cy = [math.comb(k-1, i) * poly_eval(E[i], tv) * ((-1)**i) for i in range(k)]
        cy_t = [math.comb(k-1, i) * poly_eval(dE[i], tv) * ((-1)**i) for i in range(k)]
        cy_y = [i * cy[i] for i in range(1, k)]
        rts = mp.polyroots(list(reversed(cy)), maxsteps=300, extraprec=50)
        for rt in rts:
            if abs(rt.imag) > max_im: max_im = abs(rt.imag)
            if abs(rt.imag) < mp.mpf('1e-20'):
                dydt = -poly_eval(cy_t, rt.real) / poly_eval(cy_y, rt.real)
                if dydt > max_dydt:
                    max_dydt = dydt
                    worst = (float(tv), float(rt.real), float(dydt))
    print(f"  binom-E k={k}: max|Im(y)|={mp.nstr(max_im, 5)}, max(dy/dt)={mp.nstr(max_dydt, 6)} at {worst}", flush=True)

print("\n=== 5. ROOTS OF JACOBI-STIRLING EHRHART POLY p_{k,i}(n) ===", flush=True)
import sympy as sp
n_sym = sp.Symbol('n')
def falling_binom(expr, m):
    res = sp.Integer(1)
    for r in range(m):
        res *= (expr - r)
    return sp.expand(res / sp.factorial(m))
for k in [2, 3]:
    for i in range(k+1):
        d = 3*k - i
        poly_n = sp.Integer(0)
        for j, c in enumerate(A[k][i]):
            poly_n += c * falling_binom(n_sym + d - j - 1, d)
        p_obj = sp.Poly(sp.expand(poly_n), n_sym)
        rts = [complex(r.evalf()) for r in p_obj.nroots()]
        nonreal = [r for r in rts if abs(r.imag) > 1e-7]
        print(f"  p_{{{k},{i}}}(n) (deg {d}): #nonreal={len(nonreal)}, sample={nonreal[:2]}", flush=True)
