"""Explore THREE big analytic/structural angles:
1. Check roots of the Hilbert polynomial P_{k,i}(n) = [y^i] JS(n+k, n; y) in n!
2. Check real stability of B_k(x, y) = sum_i binom(k,i) A_{k,i}(x) y^i AND other weightings for k up to 10!
3. Check the EXACT Wronskian / Cauchy-index / root margins of A_{k,i} << A_{k,i-1} and Pair (2,1) up to k = 25!
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

# 1. Test B_k(x, y) = sum binom(k,i) A_{k,i}(x) y^i for real stability!
# Necessary condition for real stability: for EVERY real x, y -> B_k(x, y) is real-rooted!
# AND for EVERY real y, x -> B_k(x, y) is real-rooted (actually only for positive directions: x=a+t, y=b+c*t with c>0)!
# Wait: for a 2-variable real polynomial to be real-stable, for EVERY real x_0, B_k(x_0, y) MUST be real-rooted in y,
# AND for EVERY real y_0, B_k(x, y_0) MUST be real-rooted in x!
A = build_A(10)
print("=== 1. TEST REAL STABILITY OF B_k(x,y) = sum binom(k,i) A_{k,i}(x) y^i ===")
for k in range(2, 10):
    # Test B_k(x0, y) at various real x0!
    found_bad = None
    for x0 in [1.0, 2.0, 5.0, 10.0, 100.0, -0.5, -1.0, -2.0, -5.0, 0.5, 0.1]:
        # poly in y: coeff of y^i is binom(k,i) * A_{k,i}(x0)
        ycoeffs = [math.comb(k, i) * sum(c * (x0**m) for m, c in enumerate(A[k][i])) for i in range(k + 1)]
        r = np.roots(list(reversed(ycoeffs)))
        if np.max(np.abs(r.imag)) > 1e-7:
            found_bad = (f"x0={x0}", float(np.max(np.abs(r.imag))))
            break
    if not found_bad:
        for y0 in [-100.0, -10.0, -2.0, -1.0, -0.5, -0.1, 0.1, 0.5, 1.0, 5.0, 50.0]:
            # poly in x
            dmax = 2*k - 1
            xcoeffs = [0.0] * (dmax + 1)
            for i in range(k + 1):
                w = math.comb(k, i) * (y0 ** i)
                for m, c in enumerate(A[k][i]):
                    xcoeffs[m] += w * c
            while len(xcoeffs) > 1 and abs(xcoeffs[-1]) < 1e-12:
                xcoeffs.pop()
            r = np.roots(list(reversed(xcoeffs)))
            if np.max(np.abs(r.imag)) > 1e-6:
                found_bad = (f"y0={y0}", float(np.max(np.abs(r.imag))))
                break
    print(f"  k={k}: {'FAIL at ' + str(found_bad) if found_bad else 'PASSED all axis tests!'}")

# Also test general weighted sum_i w_{k,i} A_{k,i}(x) y^i for real stability in y at x=1!
# Note A_{k,i}(1) = number of JSP_{k,i}! What is A_{k,i}(1)?
print("\n=== 2. CHECK A_{k,i}(1) AND B_k(x, y) WITH BINOM(k,i) AT x=0 AND x->inf ===")
for k in range(1, 7):
    vals1 = [sum(A[k][i]) for i in range(k + 1)]
    vals0 = [A[k][i][0] for i in range(k + 1)]
    print(f"  k={k}: A(0)={vals0}, A(1)={vals1}")

# 3. Check roots of Hilbert polynomial P_{k,i}(n) in n for k=1..4!
print("\n=== 3. ROOTS OF HILBERT POLYNOMIAL P_{k,i}(n) ===")
n_sym = sp.Symbol('n')
for k in range(1, 5):
    for i in range(k + 1):
        # sum_{n>=0} P(n) x^n = x A_{k,i}(x) / (1-x)^{3k-i+1}
        # => P(n) = sum_{m=0}^{deg} a_m * binom(n - 1 - m + (3k-i), 3k-i)
        D = 3*k - i
        P = sum(c * sp.binomial(n_sym - 1 - m + D, D) for m, c in enumerate(A[k][i]))
        P_poly = sp.Poly(sp.expand(P), n_sym)
        rr = sp.count_roots(P_poly)
        print(f"  k={k}, i={i}: deg(P)={P_poly.degree()}, real_roots={rr}")
