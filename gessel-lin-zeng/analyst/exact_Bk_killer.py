"""Exact rational check of Rayleigh determinant Px*Py - P*Pxy for B_k(x, y) at x = 0!
Wait: at x = 0:
  A_{k,i}(0) = binom(k, i),
  A'_{k,i}(0) = [x^1] A_{k,i}(x) = ???
Let's compute EXACT B_k(0, y) and B_{k,x}(0, y), and check EXACTLY when W_y(B_k(0,y), B_{k,x}(0,y)) < 0!!
"""
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

A = build_A(12)
y = sp.Symbol('y')
for k in range(2, 13):
    P = sum(math.comb(k, i) * A[k][i][0] * y**i for i in range(k + 1))
    Px = sum(math.comb(k, i) * (A[k][i][1] if len(A[k][i]) > 1 else 0) * y**i for i in range(k + 1))
    Py = sp.diff(P, y)
    Pxy = sp.diff(Px, y)
    R = sp.Poly(sp.expand(Px * Py - P * Pxy), y, domain='ZZ')
    # Check if R has any real roots (sign changes!)
    rc = sp.count_roots(R)
    # Check sign at y = -1, -3/4, -1/2, 1
    vals = {v: int(R.eval(sp.Rational(v[0], v[1]))) for v in [(-1, 1), (-3, 4), (-1, 2), (-4, 1), (1, 1)]}
    neg_at = [v for v, val in vals.items() if val < 0]
    print(f"k={k:2d}: deg(R)={R.degree():2d}, real_roots_of_R={int(rc):2d}, neg_samples={neg_at}")
