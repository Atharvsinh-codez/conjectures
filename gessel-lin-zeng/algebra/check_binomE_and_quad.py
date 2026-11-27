"""
Fast check of:
1. Newton disproof of B_35(x, -x*y) at x = 7e-17 (confirming complex y-roots!)
2. binom-E_k(t, y) stability & Newton inequality E_{k,i}(t)^2 >= E_{k,i-1}(t) E_{k,i+1}(t) on R for k=10..40!
3. Roots of Jacobi-Stirling Ehrhart poly p_{k,i}(n) for k=2,3.
4. Check M^2 - X*Y and Quad stability for Laguerre L_j = x d/dx + a_j + b_j x (equal & wild).
"""
import math, random
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

print("=== 1. CONFIRM NEWTON DISPROOF OF B_k(x, -x*y) AT k=35, x=7e-17 ===", flush=True)
A35 = compute_A(35)
xv = mp.mpf('7e-17')
v14 = poly_eval(A35[35][14], xv)
v15 = poly_eval(A35[35][15], xv)
v16 = poly_eval(A35[35][16], xv)
print(f"  k=35, x=7e-17: A_15^2 - A_14*A_16 = {mp.nstr(v15*v15 - v14*v16, 8)} (ratio = {mp.nstr(v15*v15/(v14*v16), 10)})", flush=True)
# Also find exact first k where A_{k,i}(x)^2 < A_{k,i-1}(x) A_{k,i+1}(x) for some x > 0!
# Note: at x = 0+, [x^1] of (A_i^2 - A_{i-1} A_{i+1}) can be negative starting k=19, wait:
# at x=0, A_i(0)^2 - A_{i-1}(0) A_{i+1}(0) = C(k,i)^2 - C(k,i-1)C(k,i+1) > 0!
# Wait: for Newton inequality on B_k(x, -x*y) = sum C(k,i) A_{k,i}(x)(-x)^i y^i,
# the coefficient of y^i is c_i = C(k,i) A_{k,i}(x)(-x)^i, so normalized by C(k,i) it is a_i = A_{k,i}(x)(-x)^i,
# and a_i^2 - a_{i-1} a_{i+1} = x^{2i} (A_i^2 - A_{i-1} A_{i+1})!
# Wait: check at x = 7e-17 if B_35(x, -x*y) actually has complex y-roots!
cy = [math.comb(35, i) * poly_eval(A35[35][i], xv) * ((-xv)**i) for i in range(36)]
# Rescale y = z / xv so coeffs are C(35,i) A_{35,i}(xv) (-1)^i (O(1) magnitude!)
cz = [math.comb(35, i) * poly_eval(A35[35][i], xv) * ((-1)**i) for i in range(36)]
rts = mp.polyroots(list(reversed(cz)), maxsteps=400, extraprec=60)
nonreal = [r for r in rts if abs(r.imag) > mp.mpf('1e-15')]
print(f"  B_35(7e-17, -7e-17 * y): #nonreal y-roots = {len(nonreal)}, max|Im(z)| = {mp.nstr(max(abs(r.imag) for r in rts), 6)}", flush=True)

# Check earlier k for complex y-roots of B_k(x, -x*y) at small x > 0!
for k in [20, 24, 28, 30, 32, 34]:
    found = False
    for p in range(1, 20):
        xv2 = mp.mpf(10)**(-p)
        cz2 = [math.comb(k, i) * poly_eval(A35[k][i], xv2) * ((-1)**i) for i in range(k+1)]
        rts2 = mp.polyroots(list(reversed(cz2)), maxsteps=300, extraprec=50)
        mx = max(abs(r.imag) for r in rts2)
        if mx > mp.mpf('1e-15'):
            print(f"  B_{k}(1e-{p}, -x*y): COMPLEX y-roots! max|Im(z)|={mp.nstr(mx, 6)}", flush=True)
            found = True
            break
    if not found:
        print(f"  B_{k}(x, -x*y): all real on x=10^(-1..-19)", flush=True)

print("\n=== 2. STRESS TEST binom-E_k(t, y) FOR k=18, 22, 26, 30, 35 ===", flush=True)
ht35 = compute_htilde(35)
for k in [18, 22, 26, 30, 35]:
    E = {}
    dE = {}
    for i in range(k):
        p = ht35[k-1][i]
        for _ in range(i):
            p = apply_Dj(p, k+1)
        E[i] = p
        dE[i] = poly_deriv(p)
    # First check Newton inequality E_{k,i}(t)^2 >= E_{k,i-1}(t) E_{k,i+1}(t) on grid including tiny t>0, t<0, large t
    tgrid = [mp.mpf(v) for v in [-100, -10, -2, -1.01, -0.99, -0.5, -0.1, -0.01, -1e-6, 0, 1e-16, 1e-10, 1e-6, 1e-3, 0.01, 0.1, 1, 10, 100]]
    min_newton_ratio = mp.mpf('1e100')
    worst_n = None
    for i in range(1, k-1):
        for tv in tgrid:
            v0 = poly_eval(E[i], tv)
            vm = poly_eval(E[i-1], tv)
            vp = poly_eval(E[i+1], tv)
            if vm*vp > 0:
                r = (v0*v0) / (vm*vp)
                if r < min_newton_ratio:
                    min_newton_ratio = r
                    worst_n = (i, float(tv))
    # Also check y-roots & dy/dt at selected t
    max_im = mp.mpf(0)
    max_dydt = mp.mpf('-1e100')
    worst_d = None
    for tv in [mp.mpf(v) for v in [-10, -1.5, -0.99, -0.5, -0.1, 0, 1e-6, 0.1, 1, 10]]:
        cy = [math.comb(k-1, i) * poly_eval(E[i], tv) * ((-1)**i) for i in range(k)]
        cy_t = [math.comb(k-1, i) * poly_eval(dE[i], tv) * ((-1)**i) for i in range(k)]
        cy_y = [i * cy[i] for i in range(1, k)]
        rts = mp.polyroots(list(reversed(cy)), maxsteps=400, extraprec=60)
        for rt in rts:
            if abs(rt.imag) > max_im: max_im = abs(rt.imag)
            if abs(rt.imag) < mp.mpf('1e-15'):
                dydt = -poly_eval(cy_t, rt.real) / poly_eval(cy_y, rt.real)
                if dydt > max_dydt:
                    max_dydt = dydt
                    worst_d = (float(tv), float(rt.real), float(dydt))
    print(f"  binom-E k={k:2d}: min E_i^2/(E- E+)={mp.nstr(min_newton_ratio, 6)} at {worst_n}; max|Im(y)|={mp.nstr(max_im, 5)}; max(dy/dt)={mp.nstr(max_dydt, 6)} at {worst_d}", flush=True)

print("\n=== 3. ROOTS OF JACOBI-STIRLING EHRHART POLY p_{k,i}(n) ===", flush=True)
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
        for j, c in enumerate(A35[k][i]):
            poly_n += c * falling_binom(n_sym + d - j - 1, d)
        p_obj = sp.Poly(sp.expand(poly_n), n_sym)
        rts = [complex(r.evalf()) for r in p_obj.nroots()]
        nonreal = [r for r in rts if abs(r.imag) > 1e-7]
        print(f"  p_{{{k},{i}}}(n) (deg {d}): #nonreal={len(nonreal)}, sample={nonreal[:2]}", flush=True)
