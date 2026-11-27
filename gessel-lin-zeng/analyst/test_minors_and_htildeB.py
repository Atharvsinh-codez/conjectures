"""Test:
1. UHP stability of B_k(x,y) = sum binom(k,i) A_{k,i}(x) y^i for k=2..14
   AND Htilde_B_k(t,y) = sum binom(k,i) htilde_{k,i}(t) y^i for k=2..10!
2. Check exact sign of 2x2 minors M_{a,b} = f_{a,b}^2 - f_{a-1,b-1} f_{a+1,b+1} and
   general 2x2 determinants in the grid (D_{k+1}^r htilde_{k-1, i})!
3. Check whether B_k(x,y) satisfies a low-order PDE relating B_k and B_{k-1}!
4. Also check whether the 3-term or 2-term Turan / log-concavity in x has ALL POSITIVE COEFFICIENTS!
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

def poly_sub(a, b):
    return poly_add(a, [-v for v in b])

def poly_mul(a, b):
    res = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            res[i + j] += u * v
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

def apply_Dtilde(j, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(poly):
        res[m] += c * (m + 1)
        res[m + 1] += c * (m + j + 1)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def build_htilde(kmax):
    p = [[None] * (k + 1) for k in range(kmax + 1)]
    p[1][0] = [1, 2]
    p[1][1] = [1]
    for k in range(2, kmax + 1):
        p[k][0] = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][0]))
        for i in range(1, k):
            p[k][i] = apply_Dtilde(k, poly_add(apply_Dtilde(k, p[k - 1][i]), p[k - 1][i - 1]))
        p[k][k] = apply_Dtilde(k, p[k - 1][k - 1])
    return p

A = build_A(14)
H = build_htilde(10)

print("=== 1. UHP STABILITY TEST OF Htilde_B_k(t, y) = sum binom(k,i) htilde_{k,i}(t) y^i ===")
# Note: in t, deg_t(htilde_{k,0}) = 2k-1 is CONSTANT for all y!
# Check Im(t) <= 0 for all y in H+!
rng = np.random.default_rng(123)
test_ys = [complex(re_y, im_y)
           for re_y in [-50, -10, -2, -1, -0.5, -0.1, 0, 0.1, 0.5, 1, 2, 10, 50]
           for im_y in [1e-4, 1e-2, 0.1, 1.0, 10.0]]
for _ in range(200):
    test_ys.append(complex(rng.normal(0, 1) * (10 ** rng.uniform(-2, 2)), 10 ** rng.uniform(-4, 2)))

for k in range(1, 9):
    max_im_t = -1e9
    worst_y = None
    dmax = 2*k - 1
    for y in test_ys:
        coeffs = np.zeros(dmax + 1, dtype=np.complex128)
        for i in range(k + 1):
            w = math.comb(k, i) * (y ** i)
            for m, c in enumerate(H[k][i]):
                coeffs[m] += w * c
        r = np.roots(coeffs[::-1])
        mi = float(np.max(r.imag))
        if mi > max_im_t:
            max_im_t = mi
            worst_y = y
    print(f"  Htilde_B_{k}(t,y): max Im(t) = {max_im_t:+.3e} (at y={worst_y})")

print("\n=== 2. CHECK COEFF POSITIVITY OF A_{k,i}^2 - A_{k,i-1} A_{k,i+1} IN x (k=2..10) ===")
for k in range(2, 11):
    neg_cases = []
    for i in range(1, k):
        D = poly_sub(poly_mul(A[k][i], A[k][i]), poly_mul(A[k][i-1], A[k][i+1]))
        if any(c < 0 for c in D):
            neg_cases.append((i, min(D)))
    print(f"  k={k:2d}: A_{{k,i}}^2 - A_{{k,i-1}} A_{{k,i+1}} neg coeffs = {neg_cases}")

print("\n=== 3. CHECK COEFF POSITIVITY OF (S_{L-1} A_{k-1,i-1})^2 - A_{k-1,i-2}*(S_L S_{L-1} A_{k-1,i}) IN x ===")
# Careful: what is the x-variable version of M^2 - X*Y?
# Recall in t: X = htilde_{k-1,i-2}, M = D_{k+1} htilde_{k-1,i-1}, Y = D_{k+1}^2 htilde_{k-1,i}.
# In x: D_{k+1} on htilde_{k-1,j} (deg 2k-3-j) corresponds to S_{(2k-3-j) + k} = S_{3k-3-j}!
# Check! For j = i-1: L = 3k-2-i. So M <-> S_{3k-2-i} A_{k-1,i-1}!
# For j = i: first D_{k+1} is S_{3k-3-i}, second D_{k+1} (on deg 2k-2-i) is S_{3k-2-i}!
# So Y <-> S_{3k-2-i} S_{3k-3-i} A_{k-1,i}!!
# And X <-> A_{k-1, i-2}!
for k in range(3, 11):
    neg_cases = []
    for i in range(2, k):
        L = 3*k - 2 - i
        X = A[k-1][i-2]
        M = apply_S(L, A[k-1][i-1])
        Y = apply_S(L, apply_S(L - 1, A[k-1][i]))
        D = poly_sub(poly_mul(M, M), poly_mul(X, Y))
        if any(c < 0 for c in D):
            neg_cases.append((i, min(D)))
    print(f"  k={k:2d}: M^2 - X*Y (in x) neg coeffs = {neg_cases}")

print("\n=== 4. CHECK 2x2 MINORS A_{k,i} A_{k,j} - A_{k,i-1} A_{k,j+1} (i <= j) IN x FOR k=2..8 ===")
for k in range(2, 9):
    neg_minors = []
    for i in range(1, k):
        for j in range(i, k):
            D = poly_sub(poly_mul(A[k][i], A[k][j]), poly_mul(A[k][i-1], A[k][j+1]))
            if any(c < 0 for c in D):
                neg_minors.append((i, j, min(D)))
    print(f"  k={k}: neg 2x2 Hankel/Toeplitz minors in x = {neg_minors}")

print("\n=== 5. CHECK EVEN BETTER: WHAT ABOUT THE NORMALIZED / BINOMIAL RATIO MINOR? ===")
# Check (A_{k,i}/binom(k,i))^2 vs (A_{k,i-1}/binom(k,i-1))*(A_{k,i+1}/binom(k,i+1))!
# We already know at x=1 that's NEGATIVE (log-convex)!
# What about (A_{k,i}/A_{k,i}(0))? Same as A_{k,i}/binom(k,i).
# What about the Wronskian of A_{k,i-1} and A_{k,i}: A_{k,i-1}' A_{k,i} - A_{k,i-1} A_{k,i}'!
# Wait: since A_{k,i} << A_{k,i-1}, (A_{k,i-1}/A_{k,i})' > 0 on all R => W = A_{k,i-1}' A_{k,i} - A_{k,i-1} A_{k,i}' > 0 on R!
# Does W have ALL POSITIVE COEFFICIENTS in x?!
def poly_deriv(p):
    return [m * c for m, c in enumerate(p)][1:] if len(p) > 1 else [0]

for k in range(2, 9):
    neg_W = []
    for i in range(1, k + 1):
        f, g = A[k][i-1], A[k][i]
        W = poly_sub(poly_mul(poly_deriv(f), g), poly_mul(f, poly_deriv(g)))
        if any(c < 0 for c in W):
            neg_W.append((i, min(W)))
    print(f"  k={k}: Wronskian W(A_{{k,i-1}}, A_{{k,i}}) neg x-coeffs = {neg_W}")
