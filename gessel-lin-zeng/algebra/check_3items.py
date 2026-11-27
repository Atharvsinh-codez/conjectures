"""
Fast exact checks for 3 remaining items:
(a) Plucker ratio W(A,B)*W(C,E) / W(B,C)^2 and positivity of W(A,C) at k=10,11 using exact integer / mpmath
(b) B_k(x, -x*y) = sum_{i=0}^k C(k,i) A_{k,i}(x) (-x)^i y^i stability & Sturm tests k=2..18
(c) Turan h_{k,i}(t)^2 >= h_{k,i-1}(t)*h_{k,i+1}(t) on (-1,0) at k=35,36 (and ULC on [-1,0])
"""
import math
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 60

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

def poly_scale(a, c):
    return [c * x for x in a]

def apply_Dj(p, j):
    # D_j p = t(1+t)p' + (1 + j t) p
    dp = poly_deriv(p)
    # t(1+t)dp = t dp + t^2 dp
    t_dp = [0] + dp
    t2_dp = [0, 0] + dp
    p1 = p
    jt_p = [0] + [j * c for c in p]
    return poly_add(poly_add(t_dp, t2_dp), poly_add(p1, jt_p))

def apply_SL_x(p, L):
    # S_L p(x) = (1 + L x) p(x) + x(1-x) p'(x)
    dp = poly_deriv(p)
    x_dp = [0] + dp
    x2_dp = [0, 0] + [-c for c in dp]
    Lx_p = [0] + [L * c for c in p]
    return poly_add(poly_add(p, Lx_p), poly_add(x_dp, x2_dp))

def compute_htilde(K):
    # htilde_{1,0} = 1 + 2t, htilde_{1,1} = 1
    h = {1: {0: [1, 2], 1: [1]}}
    for k in range(2, K + 1):
        h[k] = {}
        for i in range(k + 1):
            p1 = apply_Dj(h[k-1][i], k + 1) if i <= k - 1 else [0]
            p2 = h[k-1][i-1] if i >= 1 else [0]
            U = poly_add(p1, p2)
            h[k][i] = apply_Dj(U, k + 1)
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

def wronskian(f, g):
    # W(f, g) = f g' - f' g
    return poly_sub(poly_mul(f, poly_deriv(g)), poly_mul(poly_deriv(f), g))

def poly_eval_mp(p, x):
    acc = mp.mpc(0) if isinstance(x, mp.mpc) else mp.mpf(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc

def exact_sturm_count_real_roots(p):
    # Exact Sturm using Fractions on all of R (after removing common gcd)
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    if len(p) <= 1:
        return 0
    # Primitive integer poly
    g = 0
    for c in p:
        g = math.gcd(g, abs(c))
    if g > 1:
        p = [c // g for c in p]
    # Use sympy if fast, or rational Sturm
    import sympy as sp
    t = sp.Symbol('t')
    poly = sp.Poly(list(reversed(p)), t)
    return sp.count_roots(poly)

print("=== CHECK (a): Plucker ratio W(A,B)W(C,E)/W(B,C)^2 and W(A,C) > 0 at k=8..11 ===")
ht = compute_htilde(11)
for k in [8, 9, 10, 11]:
    # E_{k,i} = D_{k+1}^i htilde_{k-1,i}
    E = {}
    for i in range(k):
        p = ht[k-1][i]
        for _ in range(i):
            p = apply_Dj(p, k+1)
        E[i] = p
    # Check quad quadruples A=E[i-3], B=E[i-2], C=E[i-1], E_poly=E[i]
    # Actually let's check both (E_{k,i-3}, E_{k,i-2}, E_{k,i-1}, E_{k,i})
    # Wait: in Q_{k,i} = X - 2 M y + Y y^2, X = h_{k-1,i-2}, M = D h_{k-1,i-1}, Y = D^2 h_{k-1,i}.
    # In Sub-step A: Q_{k,i-1} + D Q_{k,i} uses the 4 polynomials:
    # A = h_{k-1,i-3}, B = D h_{k-1,i-2}, C = D^2 h_{k-1,i-1}, E_p = D^3 h_{k-1,i}
    # (which are D^{-(i-3)} of E_{k,i-3..i}).
    max_ratio = mp.mpf(0)
    min_wac = mp.mpf('1e100')
    worst_info = None
    for i in range(3, k):
        A_p = ht[k-1][i-3]
        B_p = apply_Dj(ht[k-1][i-2], k+1)
        C_p = apply_Dj(apply_Dj(ht[k-1][i-1], k+1), k+1)
        E_p = apply_Dj(apply_Dj(apply_Dj(ht[k-1][i], k+1), k+1), k+1)
        # Recall interlacing direction: h_{k-1,i} << h_{k-1,i-1}, so D^3 h_{k-1,i} << D^2 h_{k-1,i-1} << D h_{k-1,i-2} << h_{k-1,i-3} ??
        # Wait: let's check sign of W(A_p, B_p)!
        wAB = wronskian(A_p, B_p)
        wBC = wronskian(B_p, C_p)
        wCE = wronskian(C_p, E_p)
        wAC = wronskian(A_p, C_p)
        # Evaluate on grid including near roots in (-1,0) and outside
        grid = [mp.mpf(x) for x in [-100, -10, -5, -2, -1.5, -1.1, -1.01, -1, -0.99, -0.9, -0.7, -0.5, -0.3, -0.1, -0.01, 0, 0.1, 1, 10, 100]]
        # Also finer grid in (-1, 0)
        grid += [mp.mpf(-m)/40 for m in range(1, 40)]
        for tv in grid:
            vAB = poly_eval_mp(wAB, tv)
            vBC = poly_eval_mp(wBC, tv)
            vCE = poly_eval_mp(wCE, tv)
            vAC = poly_eval_mp(wAC, tv)
            # Note sign: if A << B then W(A,B) = A B' - A' B < 0 or > 0?
            # If A << B (root of A < root of B): A(t) = t+2, B(t) = t+1 => A B' - A' B = (t+2)-(t+1) = +1 > 0.
            # Here h_{k-1,i} << h_{k-1,i-1}, wait: D shifts right! So check sign at t=0.
            r = (vAB * vCE) / (vBC * vBC)
            if r > max_ratio:
                max_ratio = r
                worst_info = (i, float(tv), float(vAB), float(vBC), float(vCE), float(vAC))
            if abs(vAC) < abs(min_wac) or (vAC * poly_eval_mp(wAC, mp.mpf(0)) < 0):
                if vAC * poly_eval_mp(wAC, mp.mpf(0)) < 0:
                    min_wac = -abs(vAC)
                elif min_wac > 0 and abs(vAC) < min_wac:
                    min_wac = abs(vAC)
    # Exact Sturm check on wAC and (9*wBC^2 - wAB*wCE) for worst i at k=10
    print(f"k={k}: max W_AB*W_CE/W_BC^2 = {mp.nstr(max_ratio, 10)}, worst={worst_info}, min_wac_signed={mp.nstr(min_wac, 6)}")

# Do exact root count on wAC and 9*wBC^2 - wAB*wCE for k=10, i=3..9
print("\nExact Sturm check for k=10 (A,B,C,E = D^{0,1,2,3} h_{9,i-3..i}):")
k = 10
for i in range(3, k):
    A_p = ht[k-1][i-3]
    B_p = apply_Dj(ht[k-1][i-2], k+1)
    C_p = apply_Dj(apply_Dj(ht[k-1][i-1], k+1), k+1)
    E_p = apply_Dj(apply_Dj(apply_Dj(ht[k-1][i], k+1), k+1), k+1)
    wAB = wronskian(A_p, B_p)
    wBC = wronskian(B_p, C_p)
    wCE = wronskian(C_p, E_p)
    wAC = wronskian(A_p, C_p)
    # Check 1*wBC^2 - wAB*wCE and 9*wBC^2 - wAB*wCE
    diff1 = poly_sub(poly_mul(wBC, wBC), poly_mul(wAB, wCE))
    diff9 = poly_sub(poly_scale(poly_mul(wBC, wBC), 9), poly_mul(wAB, wCE))
    r_ac = exact_sturm_count_real_roots(wAC)
    r_d1 = exact_sturm_count_real_roots(diff1)
    r_d9 = exact_sturm_count_real_roots(diff9)
    print(f"  k=10, i={i}: roots(W_AC)={r_ac} (val0={wAC[0]}), roots(W_BC^2 - W_AB*W_CE)={r_d1} (val0={diff1[0]}), roots(9*W_BC^2 - W_AB*W_CE)={r_d9}")

print("\n=== CHECK (b): B_k(x, -x*y) = sum C(k,i) A_{k,i}(x) (-x)^i y^i ===")
A = compute_A(20)
# Check at x = 0+? No, at x->0, (-x)^i vanishes! Wait: substitute y -> y/(-x) recovers B_k(x, y)!
# HOLD ON! B_k(x, -x y) is a nonlinear substitution of B_k(x, y):
# Wait: for x in UHP (Im(x) > 0) and y in UHP (Im(y) > 0), what is Im(-x y)?
# If x = r e^{i theta}, y = s e^{i phi} with theta, phi in (0, pi), then -x y = r s e^{i(pi + theta + phi)}!
# Its argument can be in (pi, 3pi), so Im(-x y) CAN be > 0 OR < 0!
# Let's check Rayleigh determinant and UHP roots for B_k(x, -x y) for k=10..20!
for k in [8, 10, 12, 14, 16, 18, 20]:
    # Check at x -> -infty (leading terms of A_{k,i}(x) (-x)^i):
    # deg A_{k,i} = 2k-i-1, so deg(A_{k,i}(x) (-x)^i) = 2k-1 for ALL i!
    # Wait: leading coeff of A_{k,i}(x) is lc_i, and (-x)^i gives (-1)^i lc_i!
    # So as x -> infty, B_k(x, -xy) / x^{2k-1} -> sum_{i=0}^k C(k,i) lc(A_{k,i}) (-y)^i!
    # Let's check if THAT polynomial in y is real-rooted!
    lc_poly = [math.comb(k, i) * A[k][i][-1] * ((-1)**i) for i in range(k+1)]
    r_lc = exact_sturm_count_real_roots(lc_poly)
    # Also check subleading term ratio (x^{2k-2})!
    sub_poly = [math.comb(k, i) * A[k][i][-2] * ((-1)**i) for i in range(k+1)]
    # Check UHP stability of B_k(x, -x y) at several x in UHP:
    # Real-stable <=> for Im(x) > 0, ALL y-roots satisfy Im(y) <= 0 (so max_im_y <= 0)
    # and for Im(y) > 0, ALL x-roots satisfy Im(x) <= 0 (so max_im_x <= 0).
    max_im_y = mp.mpf('-1e100')
    worst_x = None
    for xval in [mp.mpc(-10, '1e-6'), mp.mpc(-1, '1e-6'), mp.mpc(-0.1, '1e-6'), mp.mpc(-0.01, '1e-6'), mp.mpc(-0.001, '1e-7'),
                 mp.mpc(0.01, '1e-6'), mp.mpc(0.1, '1e-6'), mp.mpc(1, '1e-6'), mp.mpc(10, '1e-6'), mp.mpc(0, 1), mp.mpc(-0.5, 0.5), mp.mpc(0.5, 0.5)]:
        coeffs_y = [math.comb(k, i) * poly_eval_mp(A[k][i], xval) * ((-xval)**i) for i in range(k+1)]
        roots = mp.polyroots(list(reversed(coeffs_y)), maxsteps=200, extraprec=40)
        for rt in roots:
            if rt.imag > max_im_y:
                max_im_y = rt.imag
                worst_x = (xval, rt)
    max_im_x = mp.mpf('-1e100')
    worst_y = None
    for yval in [mp.mpc(-1, '1e-6'), mp.mpc(-0.1, '1e-6'), mp.mpc(0.1, '1e-6'), mp.mpc(1, '1e-6'), mp.mpc(10, '1e-6'), mp.mpc(0, 1)]:
        px = [mp.mpc(0)] * (2*k)
        for i in range(k+1):
            w = math.comb(k, i) * ((-1)**i) * (yval**i)
            for d, c in enumerate(A[k][i]):
                px[d + i] += w * c
        roots_x = mp.polyroots(list(reversed(px)), maxsteps=200, extraprec=40)
        for rt in roots_x:
            if rt.imag > max_im_x:
                max_im_x = rt.imag
                worst_y = (yval, rt)
    print(f"k={k}: lc_poly real roots={r_lc}/{k}; UHP x -> max Im(y)={mp.nstr(max_im_y, 6)}; UHP y -> max Im(x)={mp.nstr(max_im_x, 6)}")
    if max_im_y > mp.mpf('1e-12'):
        print(f"   FAIL in y at x={worst_x[0]}: root y={worst_x[1]}")
    if max_im_x > mp.mpf('1e-12'):
        print(f"   FAIL in x at y={worst_y[0]}: root x={worst_y[1]}")

print("\n=== CHECK (c): Turan & ULC of htilde_{k,i}(t) on [-1, 0] at k=35, 36 ===")
ht36 = compute_htilde(36)
for k in [35, 36]:
    min_ulc = mp.mpf('1e100')
    min_turan = mp.mpf('1e100')
    worst_u = None
    worst_t = None
    # Check both near t=0 (x=0-) and dense grid in [-1, 0]
    tgrid = [-mp.mpf(10)**(-p) for p in [1, 2, 3, 4, 6, 8, 12, 16, 20]] + [-mp.mpf(m)/100 for m in range(1, 100)] + [mp.mpf(-1) + mp.mpf(10)**(-p) for p in [3, 6, 10]]
    for i in range(1, k):
        c_ratio = mp.mpf(math.comb(k, i)**2) / (math.comb(k, i-1) * math.comb(k, i+1))
        for tv in tgrid:
            v0 = poly_eval_mp(ht36[k][i], tv)
            vm = poly_eval_mp(ht36[k][i-1], tv)
            vp = poly_eval_mp(ht36[k][i+1], tv)
            turan_val = v0*v0 - vm*vp
            ulc_val = v0*v0 - c_ratio * vm*vp
            if turan_val < 0:
                if turan_val < min_turan:
                    min_turan = turan_val
                    worst_t = (i, float(tv), float(turan_val))
            elif min_turan > 0 and turan_val < min_turan:
                min_turan = turan_val
                worst_t = (i, float(tv), float(turan_val))
            if ulc_val < 0:
                if ulc_val < min_ulc:
                    min_ulc = ulc_val
                    worst_u = (i, float(tv), float(ulc_val))
            elif min_ulc > 0 and ulc_val < min_ulc:
                min_ulc = ulc_val
                worst_u = (i, float(tv), float(ulc_val))
    print(f"k={k}: min(h^2 - h- h+) = {mp.nstr(min_turan, 6)} at {worst_t}; min(ULC) = {mp.nstr(min_ulc, 6)} at {worst_u}")

