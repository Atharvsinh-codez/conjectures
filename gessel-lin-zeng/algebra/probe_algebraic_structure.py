"""
Probe algebraic & operator structure of:
1. Universal operator product prod_{j=1}^k (T_j + w_j y)(1) for generalized operators:
   - Laguerre L_j = x d/dx + a_j + b_j x (roots in (-inf, 0))
   - Hermite H_j = d/dx + a_j + b_j x (roots in R)
   - Different branch points T_j = (t+u_j)(t+v_j) d/dt + a_j + b_j t
   - Pure differential in u: W_j^{-1} d/du W_j for general log-concave W_j!
2. Single-operator Turan: (W^{(n)}(u))^2 >= W^{(n-1)}(u) W^{(n+1)}(u) for W(u) = e^{au}/(1+e^u)^b,
   and for general PF_inf functions (e.g. e^{-u^2}, e^{au} e^{-e^u})!
3. Check (b) B_k(x, -xy) fast via mpmath (without slow sympy count_roots) and Check (c) k=35,36 Turan/ULC on [-1,0].
4. Check roots of Jacobi-Stirling Ehrhart polynomial p_{k,i}(n).
"""
import math, random
import mpmath as mp
import numpy as np

mp.mp.dps = 50

def poly_add(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]

def poly_sub(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(n)]

def poly_mul(a, b):
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j, bj in enumerate(b):
            res[i+j] += ai * bj
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def poly_deriv(a):
    if len(a) <= 1: return [0]
    return [i * a[i] for i in range(1, len(a))]

def poly_eval_mp(p, x):
    acc = mp.mpc(0) if isinstance(x, mp.mpc) else mp.mpf(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc

def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p

def apply_SL_x(p, L):
    dp = poly_deriv(p)
    x_dp = [0] + dp
    x2_dp = [0, 0] + [-c for c in dp]
    Lx_p = [0] + [L * c for c in p]
    return poly_trim(poly_add(poly_add(p, Lx_p), poly_add(x_dp, x2_dp)))

def apply_Dj(p, j):
    dp = poly_deriv(p)
    t_dp = [0] + dp
    t2_dp = [0, 0] + dp
    jt_p = [0] + [j * c for c in p]
    return poly_add(poly_add(t_dp, t2_dp), poly_add(p, jt_p))

def compute_htilde(K):
    h = {1: {0: [1, 2], 1: [1]}}
    for k in range(2, K + 1):
        h[k] = {}
        for i in range(k + 1):
            p1 = apply_Dj(h[k-1][i], k + 1) if i <= k - 1 else [0]
            p2 = h[k-1][i-1] if i >= 1 else [0]
            h[k][i] = apply_Dj(poly_add(p1, p2), k + 1)
    return h

def compute_A(K):
    A = {0: {0: [1]}}
    for k in range(1, K + 1):
        A[k] = {}
        for i in range(k + 1):
            p1 = apply_SL_x(A[k-1][i], 3*k - 3 - i) if i <= k - 1 else [0]
            p2 = A[k-1][i-1] if i >= 1 else [0]
            A[k][i] = apply_SL_x(poly_add(p1, p2), 3*k - 2 - i)
    return A

# ==========================================================
# TEST 1: Single-operator T^n(1) Turan on (-1,0) and Laguerre L^n(1) on (-inf, 0)
# ==========================================================
print("=== TEST 1: Single-operator Turan (T^n(1))^2 >= T^{n-1}(1) T^{n+1}(1) ===")
# 1a. T = t(1+t)d/dt + a + b*t (on (-1,0))
for (a, b) in [(1, 2), (1, 5), (3, 4), (0.5, 10.0), (5.0, 5.5)]:
    polys = [[mp.mpf(1)]]
    for n in range(1, 13):
        p = polys[-1]
        dp = [i * p[i] for i in range(1, len(p))]
        res = [mp.mpf(0)] * (len(p) + 1)
        for i, c in enumerate(p):
            res[i] += a * c
            res[i+1] += b * c
        for i, c in enumerate(dp):
            res[i+1] += c
            res[i+2] += c
        polys.append(res)
    # Check Turan on (-1, 0)
    min_val = mp.mpf('1e100')
    worst = None
    grid = [-mp.mpf(m)/200 for m in range(1, 200)] + [-mp.mpf(10)**(-p) for p in [3, 5, 8]] + [mp.mpf(-1) + mp.mpf(10)**(-p) for p in [3, 5, 8]]
    for n in range(1, 12):
        for tv in grid:
            v0 = poly_eval_mp(polys[n], tv)
            vm = poly_eval_mp(polys[n-1], tv)
            vp = poly_eval_mp(polys[n+1], tv)
            val = v0*v0 - vm*vp
            if val < min_val:
                min_val = val
                worst = (n, float(tv))
    print(f"  T(a={a}, b={b}): min(p_n^2 - p_{{n-1}} p_{{n+1}}) on (-1,0) = {mp.nstr(min_val, 6)} at {worst}")

# 1b. Laguerre L = x d/dx + a + b*x on (-inf, 0)
for (a, b) in [(1, 1), (0.5, 1), (3, 1), (0.1, 2)]:
    polys = [[mp.mpf(1)]]
    for n in range(1, 13):
        p = polys[-1]
        res = [mp.mpf(0)] * (len(p) + 1)
        for i, c in enumerate(p):
            res[i] += (a + i) * c
            res[i+1] += b * c
        polys.append(res)
    min_val = mp.mpf('1e100')
    worst = None
    grid = [-mp.mpf(x) for x in [1e-4, 1e-2, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 5, 8, 12, 20, 50, 100]] + [-mp.mpf(m)/20 for m in range(1, 200)]
    for n in range(1, 12):
        for xv in grid:
            v0 = poly_eval_mp(polys[n], xv)
            vm = poly_eval_mp(polys[n-1], xv)
            vp = poly_eval_mp(polys[n+1], xv)
            val = v0*v0 - vm*vp
            if val < min_val:
                min_val = val
                worst = (n, float(xv))
    print(f"  Laguerre L(a={a}, b={b}): min(p_n^2 - p_{{n-1}} p_{{n+1}}) on (-inf,0) = {mp.nstr(min_val, 6)} at {worst}")

# 1c. Hermite H = d/dx + 2*x (wait: H_n(x) = (-1)^n e^{x^2} d^n/dx^n e^{-x^2}, so H = d/dx - 2x gives (-1)^n H_n!)
# Wait: for W(u) = e^{-u^2}, (W^{(n)})^2 - W^{(n-1)} W^{(n+1)} = e^{-2u^2} (H_n^2 - H_{n-1} H_{n+1}) > 0 EVERYWHERE by Hermite Turan!
# WAIT A SECOND!! WHAT ABOUT GENERAL prod_{j=1}^k (L_j + w_j y)(1) AND prod_{j=1}^k (H_j + w_j y)(1)?!

print("\n=== TEST 2: Universal product prod (L_j + w_j y)(1) for Laguerre L_j = x d/dx + a_j + b_j x ===")
random.seed(42)
fail_L_interlace = 0
fail_L_turan = 0
for trial in range(50):
    k = 6
    seq = [1.0] # dict i -> poly
    polys = {0: [1.0]}
    for j in range(1, k + 1):
        a_j = 10 ** random.uniform(-1, 2)
        b_j = 10 ** random.uniform(-1, 2)
        w_j = 10 ** random.uniform(-1, 2)
        new_p = {}
        for i in range(j + 1):
            # L_j(polys[i]) + w_j * polys[i-1]
            res = [0.0] * (j - i + 1)
            if i in polys:
                p = polys[i]
                for m, c in enumerate(p):
                    res[m] += (m + a_j) * c
                    res[m+1] += b_j * c
            if (i - 1) in polys:
                p = polys[i-1]
                for m, c in enumerate(p):
                    res[m] += w_j * c
            new_p[i] = res
        polys = new_p
    # Check real-rootedness & interlacing q_{k,i} << q_{k,i-1} in (-inf, 0)
    roots = {}
    ok = True
    for i in range(k + 1):
        p = polys[i]
        if len(p) <= 1:
            roots[i] = []
            continue
        rts = np.roots(list(reversed(p)))
        if np.max(np.abs(np.imag(rts))) > 1e-6:
            ok = False; break
        rts = np.sort(np.real(rts))
        if rts[-1] >= 0:
            ok = False; break
        roots[i] = rts
    if ok:
        for i in range(1, k):
            r_lo = roots[i]
            r_hi = roots[i-1]
            # Expect r_hi[0] < r_lo[0] < r_hi[1] < ... < r_lo[-1] < r_hi[-1]
            for m in range(len(r_lo)):
                if not (r_hi[m] < r_lo[m] < r_hi[m+1]):
                    ok = False; break
    if not ok:
        fail_L_interlace += 1
    # Check Turan on (-inf, 0)
    tok = True
    xs = -np.logspace(-4, 3, 200)
    for i in range(1, k):
        v0 = np.polyval(list(reversed(polys[i])), xs)
        vm = np.polyval(list(reversed(polys[i-1])), xs)
        vp = np.polyval(list(reversed(polys[i+1])), xs)
        if np.min(v0*v0 - vm*vp) < -1e-9 * np.max(np.abs(v0*v0)):
            tok = False; break
    if not tok:
        fail_L_turan += 1
print(f"  Laguerre L_j = x d/dx + a_j + b_j x (50 wild trials, k=6): interlace fails={fail_L_interlace}, Turan(-inf,0) fails={fail_L_turan}")

print("\n=== TEST 3: Universal product prod (H_j + w_j y)(1) for Hermite H_j = d/dx + a_j + b_j x (b_j > 0) ===")
# Wait: H_j f = e^{-b_j x^2/2 - a_j x} d/dx (e^{b_j x^2/2 + a_j x} f) -- NO! For real roots we need e^{-b_j x^2/2} inside d/dx,
# i.e. H_j = d/dx - b_j x + a_j (or -d/dx + b_j x + a_j)!
# Check: (-d/dx + b x)(1) = b x; (-d/dx + b x)(b x) = -b + b^2 x^2 (2 real roots +/- 1/sqrt(b)!)
fail_H_interlace = 0
fail_H_turan = 0
for trial in range(50):
    k = 6
    polys = {0: [1.0]}
    for j in range(1, k + 1):
        a_j = random.uniform(-5, 5)
        b_j = 10 ** random.uniform(-1, 1)
        w_j = 10 ** random.uniform(-1, 1)
        new_p = {}
        for i in range(j + 1):
            res = [0.0] * (j - i + 1)
            if i in polys:
                p = polys[i]
                for m, c in enumerate(p):
                    if m >= 1:
                        res[m-1] += -m * c
                    res[m] += a_j * c
                    res[m+1] += b_j * c
            if (i - 1) in polys:
                p = polys[i-1]
                for m, c in enumerate(p):
                    res[m] += w_j * c
            new_p[i] = res
        polys = new_p
    roots = {}
    ok = True
    for i in range(k + 1):
        p = polys[i]
        if len(p) <= 1:
            roots[i] = []
            continue
        rts = np.roots(list(reversed(p)))
        if np.max(np.abs(np.imag(rts))) > 1e-6:
            ok = False; break
        roots[i] = np.sort(np.real(rts))
    if ok:
        for i in range(1, k):
            r_lo = roots[i]
            r_hi = roots[i-1]
            for m in range(len(r_lo)):
                if not (r_hi[m] < r_lo[m] < r_hi[m+1]):
                    ok = False; break
    if not ok:
        fail_H_interlace += 1
    tok = True
    xs = np.linspace(-20, 20, 400)
    for i in range(1, k):
        v0 = np.polyval(list(reversed(polys[i])), xs)
        vm = np.polyval(list(reversed(polys[i-1])), xs)
        vp = np.polyval(list(reversed(polys[i+1])), xs)
        if np.min(v0*v0 - vm*vp) < -1e-9 * np.max(np.abs(v0*v0)):
            tok = False; break
    if not tok:
        fail_H_turan += 1
print(f"  Hermite H_j = -d/dx + b_j x + a_j (50 wild trials, k=6): interlace fails={fail_H_interlace}, Turan(R) fails={fail_H_turan}")

print("\n=== TEST 4: Check B_k(x, -x*y) stability for k=10..18 (fast mpmath) ===")
A = compute_A(18)
for k in [8, 10, 12, 14, 16, 18]:
    # Check leading coeff poly in y
    lc_poly = [mp.mpf(math.comb(k, i) * A[k][i][-1] * ((-1)**i)) for i in range(k+1)]
    rts_lc = mp.polyroots(list(reversed(lc_poly)), maxsteps=200, extraprec=30)
    max_im_lc = max(abs(r.imag) for r in rts_lc)
    max_im_y = mp.mpf('-1e100')
    worst_x = None
    for xval in [mp.mpc(-10, '1e-6'), mp.mpc(-1, '1e-6'), mp.mpc(-0.1, '1e-6'), mp.mpc(-0.01, '1e-6'), mp.mpc(-0.0001, '1e-8'),
                 mp.mpc(0.01, '1e-6'), mp.mpc(0.1, '1e-6'), mp.mpc(1, '1e-6'), mp.mpc(10, '1e-6'), mp.mpc(0, 1), mp.mpc(-0.5, 0.5), mp.mpc(0.5, 0.5)]:
        coeffs_y = [math.comb(k, i) * poly_eval_mp(A[k][i], xval) * ((-xval)**i) for i in range(k+1)]
        roots = mp.polyroots(list(reversed(coeffs_y)), maxsteps=200, extraprec=30)
        for rt in roots:
            if rt.imag > max_im_y:
                max_im_y = rt.imag
                worst_x = (xval, rt)
    print(f"  k={k}: lc_poly max|Im|={mp.nstr(max_im_lc, 5)}; UHP x -> max Im(y)={mp.nstr(max_im_y, 6)}")
    if max_im_y > mp.mpf('1e-12'):
        print(f"    FAIL at x={worst_x[0]} -> y={worst_x[1]}")

print("\n=== TEST 5: Roots of Jacobi-Stirling Ehrhart polynomial p_{k,i}(n) for k=1..4 ===")
# sum_{n>=0} p_{k,i}(n) x^n = x A_{k,i}(x) / (1-x)^{3k-i+1}
# Recover p_{k,i}(n) by binomial interpolation: p_{k,i}(n) = sum_j a_j C(n + 3k-i - 1 - j, 3k-i)
import sympy as sp
n_sym = sp.Symbol('n')
for k in [2, 3]:
    for i in range(k + 1):
        d_den = 3*k - i
        poly_n = sp.Integer(0)
        for j, c in enumerate(A[k][i]):
            # [x^n] x^{j+1} / (1-x)^{d_den+1} = C(n - j - 1 + d_den, d_den)
            poly_n += c * sp.binomial(n_sym + d_den - j - 1, d_den)
        poly_n = sp.expand(poly_n)
        # Numeric roots
        p_obj = sp.Poly(poly_n, n_sym)
        rts = [complex(r.evalf()) for r in p_obj.nroots()]
        nonreal = [r for r in rts if abs(r.imag) > 1e-7]
        print(f"  p_{{{k},{i}}}(n) (deg {d_den}): #nonreal roots = {len(nonreal)}, sample={nonreal[:2]}")
