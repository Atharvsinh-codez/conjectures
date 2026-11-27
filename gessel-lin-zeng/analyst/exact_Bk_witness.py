"""Exact Sturm witness killing real stability of B_10(x, y) = sum binom(10,i) A_{10,i}(x) y^i!"""
import math
import numpy as np
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
k = 10
y = sp.Symbol('y')
P0 = sum(math.comb(k, i) * A[k][i][0] * y**i for i in range(k + 1))
Px0 = sum(math.comb(k, i) * A[k][i][1] * y**i for i in range(k + 1))
Py0 = sp.diff(P0, y)
Pxy0 = sp.diff(Px0, y)
R0 = sp.expand(Px0 * Py0 - P0 * Pxy0)
print(f"Exact R(0, -3/4) = {R0.subs(y, sp.Rational(-3, 4))}")

# Compute y_r'(0) = -Px0(y_r) / Py0(y_r) at each root y_r of P0(y) = sum binom(10,i)^2 y^i!
roots_y = sorted(float(r.evalf()) for r in sp.real_roots(sp.Poly(P0, y)))
print("Roots y_r of B_10(0, y) and y_r'(0) = -Px(0, y_r) / Py(0, y_r):")
for r in roots_y:
    dy = -float(Px0.subs(y, r)) / float(Py0.subs(y, r))
    print(f"  y_r = {r:+.8f}, y_r'(0) = {dy:+.8e}")

# Also test exact positive line near the root with y_r'(0) > 0!
# Wait: if y_r'(0) > 0 at some root y_r, then for y = y_r + i*eps (wait: for x = i*eps, y(i*eps) = y_r + i*eps*y_r'(0) in H+!)
# Even simpler: check exact Sturm of B_10(x, y_0) for y_0 slightly to the right/left, OR along line (s, y_r + c*s)!
# Wait: if y_r'(0) > 0, then on the line y = y_r + c*x with c = y_r'(0) > 0, the curve y = y_r(x) is TANGENT to the line at x=0!
# Slightly shifting y_0 up or down creates a pair of COMPLEX CONJUGATE roots of B_10(s, y_0 + c*s)!!
for r in roots_y:
    dy = -float(Px0.subs(y, r)) / float(Py0.subs(y, r))
    if dy > 0:
        # Try exact rational b near r and c near dy
        s = sp.Symbol('s')
        c_rat = sp.Rational(int(round(dy * 100)), 100)
        if c_rat <= 0: c_rat = sp.Rational(1, 100)
        for delta in [-1, 0, 1, -2, 2, -5, 5, -10, 10]:
            b_rat = sp.Rational(int(round(r * 1000)) + delta, 1000)
            expr = sum(math.comb(k, i) * sum(coeff * s**m for m, coeff in enumerate(A[k][i])) * (b_rat + c_rat * s)**i for i in range(k + 1))
            P = sp.Poly(sp.expand(expr), s, domain='QQ')
            rc = int(sp.count_roots(P))
            if rc < P.degree():
                print(f"EXACT LINE KILLER FOR B_10: B_10(s, {b_rat} + {c_rat}*s): deg={P.degree()}, real_roots={rc} < {P.degree()}!!")
                break
        break
