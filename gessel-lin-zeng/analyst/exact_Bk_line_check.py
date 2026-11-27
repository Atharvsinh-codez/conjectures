"""Exact integer Sturm check of B_7, B_8, B_10 on the suspect positive lines & exact Rayleigh check!"""
import math
import sympy as sp

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_S(L, p):
    d = len(p) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(p):
        res[m] += c * (m + 1)
        res[m + 1] += c * (L - m)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def build_A(kmax):
    A = [[None] * (k + 1) for k in range(kmax + 1)]
    A[0][0] = [1]
    for k in range(1, kmax + 1):
        A[k][0] = apply_S(3*k - 2, apply_S(3*k - 3, A[k-1][0]))
        for i in range(1, k):
            L = 3*k - 2 - i
            A[k][i] = apply_S(L, poly_add(apply_S(L - 1, A[k-1][i]), A[k-1][i-1]))
        A[k][k] = apply_S(2*k - 2, A[k-1][k-1])
    return A

A = build_A(10)
s = sp.Symbol('s')
for k, (a, b, c) in [
    (7, (sp.Rational(-4), sp.Rational(-50), sp.Rational(13))),
    (8, (sp.Rational(4), sp.Rational(50), sp.Rational(12))),
    (10, (sp.Rational(9, 2), sp.Rational(-2), sp.Rational(1, 50))),
]:
    expr = sum(math.comb(k, i) * sum(coeff * (a + s)**m for m, coeff in enumerate(A[k][i])) * (b + c * s)**i for i in range(k + 1))
    P = sp.Poly(sp.expand(expr), s, domain='QQ')
    rc = sp.count_roots(P)
    print(f"EXACT B_{k}({a}+s, {b}+{c}*s): deg={P.degree()}, real_roots={rc} -> {'PASS' if rc == P.degree() else 'FAIL!!'}")
