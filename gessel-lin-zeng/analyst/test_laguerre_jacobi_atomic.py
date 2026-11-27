"""CRITICAL DIAGNOSTIC:
1. Test Laguerre D_j = t*d/dt + c_j + m_j*t (c_j > 0, m_j > 0, lam_j >= 0) for WILD sequences (c_j, m_j, lam_j)!
2. Test Jacobi (-1, 1): D_j = (1-t^2)*d/dt + c_j - m_j*t (m_j > |c_j|) for WILD sequences!
3. Test 1000 FAST float64 adversarial trials on GLZ Atomic Theorem (N = 5, 6, 7, 8) with WILD (c_j, m_j, lam_j)
   — and ANY suspect immediately checked in EXACT sympy arithmetic!
4. ALSO check: WHY did Cauchy fail at m=[1, 1000, 1, 1000]? Look at the failure we got: (3, 0) (i=0!!)!
   Wait: at i=0, Cauchy ALREADY fails for a SINGLE product C_1000^2 C_1^2(1)!!
   (Because C_m = (1+s^2)d/ds - 2*m*s does NOT preserve real-rootedness when m drops from 1000 to 1, since (1+s^2)^{-m} has poles at +-i!)
"""
import math
import numpy as np
import sympy as sp

t = sp.Symbol('t')

def apply_Lag(c, m, poly):
    # t*d/dt + c + m*t
    d = len(poly) - 1
    res = [sp.Rational(0)] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * m
    return res

def apply_Jac(c, m, poly):
    # (1-t^2)*d/dt + c - m*t  (roots in (-1, 1) when m > |c|!)
    d = len(poly) - 1
    res = [sp.Rational(0)] * (d + 2)
    for r, coef in enumerate(poly):
        if r >= 1:
            res[r - 1] += coef * r
        res[r] += coef * c
        res[r + 1] += coef * (-r - m)
    return res

def apply_GLZ(c, m, poly):
    d = len(poly) - 1
    res = [sp.Rational(0)] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * (r + m)
    return res

def poly_add_q(a, b):
    n = max(len(a), len(b))
    res = [sp.Rational(0)] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def poly_sub_q(a, b):
    return poly_add_q(a, [-v for v in b])

def poly_mul_q(a, b):
    res = [sp.Rational(0)] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            res[i + j] += u * v
    return res

def test_family_exact(name, op_fn, steps, lo, hi):
    N = len(steps)
    F = [[None]*(k+1) for k in range(N+1)]
    F[0][0] = [sp.Rational(1)]
    for k in range(1, N + 1):
        c, m, lam = steps[k - 1]
        F[k][0] = op_fn(c, m, F[k-1][0])
        for i in range(1, k):
            F[k][i] = poly_add_q(op_fn(c, m, F[k-1][i]), [lam * v for v in F[k-1][i-1]])
        F[k][k] = [lam * v for v in F[k-1][k-1]]
    rr_f = []; lc_f = []
    for i in range(N + 1):
        P = sp.Poly(list(reversed(F[N][i])), t, domain='QQ')
        if P.degree() > 0:
            rc = int(sp.count_roots(P)) if lo is None else int(sp.count_roots(P, lo, hi))
            if rc < P.degree():
                rr_f.append((i, rc, P.degree()))
    for i in range(1, N):
        D = poly_sub_q(poly_mul_q(F[N][i], F[N][i]), poly_mul_q(F[N][i-1], F[N][i+1]))
        Pd = sp.Poly(list(reversed(D)), t, domain='QQ')
        mid = sp.Rational(-1, 2) if lo == -1 else (sp.Rational(0) if lo is None else sp.Rational(-1))
        rc = int(sp.count_roots(Pd)) if lo is None else int(sp.count_roots(Pd, lo, hi))
        if Pd.eval(mid) <= 0 or rc > 0:
            lc_f.append((i, rc))
    print(f"  {name:48s}: RR={rr_f[:3]}, LC={lc_f[:3]}")

print("=== 1. LAGUERRE D_j = t*d/dt + c_j + m_j*t (c_j, m_j > 0, roots in (-inf, 0)) ===")
test_family_exact("Lag m=[1,2,3,4,5], c=1, lam=1", apply_Lag, [(1, k, 1) for k in range(1, 6)], None, None)
test_family_exact("Lag m=[1,1000,1,1000,1], c=1, lam=1", apply_Lag, [(1, 1 if k%2 else 1000, 1) for k in range(1, 6)], None, None)
test_family_exact("Lag m=[1000,1,1000,1,1000], c=1, lam=1", apply_Lag, [(1, 1000 if k%2 else 1, 1) for k in range(1, 6)], None, None)
test_family_exact("Lag c=[1,1000,1,1000], m=[1000,1,1000,1], lam=1", apply_Lag, [(1,1000,1),(1000,1,1),(1,1000,1),(1000,1,1)], None, None)
test_family_exact("Lag c=[1000,1,1000,1], m=[1,1000,1,1000], lam=1", apply_Lag, [(1000,1,1),(1,1000,1),(1000,1,1),(1,1000,1)], None, None)
test_family_exact("Lag m=[1,1000,1,1000], lam=[1,10^6,1,10^6]", apply_Lag, [(1,1,1),(1,1000,10**6),(1,1,1),(1,1000,10**6)], None, None)
test_family_exact("Lag m=[1,1000,1,1000], lam=[10^6,1,10^6,1]", apply_Lag, [(1,1,10**6),(1,1000,1),(1,1,10**6),(1,1000,1)], None, None)

print("\n=== 2. JACOBI (-1,1): D_j = (1-t^2)*d/dt + c_j - m_j*t (|c_j| < m_j) ===")
test_family_exact("Jac m=[2,1000,2,1000], c=[1,500,-1,-500]", apply_Jac, [(1,2,1),(500,1000,1),(-1,2,1),(-500,1000,1)], -1, 1)
test_family_exact("Jac m=[1000,2,1000,2], c=[999,-1,-999,1]", apply_Jac, [(999,1000,1),(-1,2,1),(-999,1000,1),(1,2,1)], -1, 1)

print("\n=== 3. TRY EVEN HARDER TO BREAK LAGUERRE & GLZ AT N = 4, 5, 6 WITH ASYMMETRIC EXTREME PARAMS ===")
# Try separating two subsets by HUGE factor (like the (13,7) barrier!)
# Look! What if we pick N = 8 with 4 steps of m=1 and 4 steps of m=10^4 (or alternating)?!
test_family_exact("Lag 4x(m=1) then 4x(m=10^4)", apply_Lag, [(1,1,1)]*4 + [(1,10**4,1)]*4, None, None)
test_family_exact("Lag 4x(m=10^4) then 4x(m=1)", apply_Lag, [(1,10**4,1)]*4 + [(1,1,1)]*4, None, None)
test_family_exact("GLZ 4x(m=2) then 4x(m=10^4)", apply_GLZ, [(1,2,1)]*4 + [(1,10**4,1)]*4, -1, 0)
test_family_exact("GLZ 4x(m=10^4) then 4x(m=2)", apply_GLZ, [(1,10**4,1)]*4 + [(1,2,1)]*4, -1, 0)
# What about NEGATIVE c? No, c>0 is required so D_j(1) has root in (-1,0).
# What about c_j VERY close to m_j vs VERY close to 0?
test_family_exact("GLZ c=[1,999,1,999], m=[1000,1000,1000,1000]", apply_GLZ, [(1,1000,1),(999,1000,1),(1,1000,1),(999,1000,1)], -1, 0)
test_family_exact("GLZ c=[1,1,999,999], m=[1000,1000,1000,1000]", apply_GLZ, [(1,1000,1),(1,1000,1),(999,1000,1),(999,1000,1)], -1, 0)
test_family_exact("GLZ c=[999,999,1,1], m=[1000,1000,1000,1000]", apply_GLZ, [(999,1000,1),(999,1000,1),(1,1000,1),(1,1000,1)], -1, 0)
