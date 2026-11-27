"""Inspect:
1. Roots of the Hilbert polynomial P_{k,i}(n) = [y^i] JS(n+k, n; y) in n (for k=1..6)!
2. Check if P_{k,i}(n) has real roots, or what its roots in C look like!
3. Check the EXACT operator action of D_j on 2x2 determinants / Turan inequalities!
4. Check central-slice A_{k, k//2}(x) real-rootedness & Pair(2,1) at k = 25, 30!
"""
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

A = build_A(6)
n_sym = sp.Symbol('n')

def falling_binom(expr, D):
    res = sp.Integer(1)
    for j in range(D):
        res *= (expr - j)
    return res / math.factorial(D)

print("=== 1. ROOTS OF HILBERT POLYNOMIAL P_{k,i}(n) IN n (k=1..5) ===")
for k in range(1, 6):
    for i in range(k + 1):
        D = 3*k - i
        P = sum(c * falling_binom(n_sym - 1 - m + D, D) for m, c in enumerate(A[k][i]))
        P_poly = sp.Poly(sp.expand(P), n_sym, domain='QQ')
        # Factor out n*(n+1)*...*(n+k)
        triv = 1
        for r in range(0, k + 1):
            if P_poly.eval(-r) == 0:
                triv *= (n_sym + r)
        # Actually divide out all integer roots
        int_roots = [r for r in range(-D, 2) if P_poly.eval(r) == 0]
        # Remaining polynomial
        Q = P_poly
        for r in int_roots:
            while Q.eval(r) == 0:
                Q = Q.exquo(sp.Poly(n_sym - r, n_sym, domain='QQ'))
        # Roots of Q
        q_roots = [complex(z.evalf()) for z in Q.nroots()] if Q.degree() > 0 else []
        n_real = sum(1 for z in q_roots if abs(z.imag) < 1e-8)
        print(f"  k={k}, i={i}: deg={P_poly.degree():2d}, int_roots={int_roots}, rem_deg={Q.degree():2d}, rem_real={n_real:2d}, "
              f"rem_Re_range=({min(z.real for z in q_roots):.2f}, {max(z.real for z in q_roots):.2f})" if q_roots else
              f"  k={k}, i={i}: deg={P_poly.degree():2d}, int_roots={int_roots}, rem_deg=0")

print("\n=== 2. CHECK SYMMETRY / PALINDROMICITY OF P_{k,i}(n) ===")
for k in [2, 3, 4]:
    for i in range(k + 1):
        D = 3*k - i
        P = sp.expand(sum(c * falling_binom(n_sym - 1 - m + D, D) for m, c in enumerate(A[k][i])))
        # Check P(-n - c) == +- P(n)
        for c_shift in [0, 1, 2, k, k+1]:
            if sp.expand(P.subs(n_sym, -n_sym - c_shift) - P) == 0:
                print(f"  k={k}, i={i}: P(-n-{c_shift}) = +P(n)")
            elif sp.expand(P.subs(n_sym, -n_sym - c_shift) + P) == 0:
                print(f"  k={k}, i={i}: P(-n-{c_shift}) = -P(n)")
