"""Test:
1. q-deformed bivariate GF R_k(x, y) = sum_i q^{i(i-1)/2} A_{k,i}(x) y^i
   and Q_k(t, y) = sum_i q^{i(i-1)/2} htilde_{k,i}(t) y^i
   for real stability (Rayleigh at x=0 and on R^2, plus UHP y test) for q in [0.5, 0.25, 0.2, 0.1]!
2. Also test w_{k,i} = (k-i)! and w_{k,i} = 1/(k-i)! just to finish exhausting simple weights.
3. Also test whether the operator T_{j,q}(F) = D_j F(t,y) + y F(t, q*y) preserves real stability on Q_{k-1}
   and check exact Rayleigh of Q_k(t, y) at t = 0 and on (-1, 0)!
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

A = build_A(12)

print("=== 1. EXACT RAYLEIGH AT x=0 FOR q^{i(i-1)/2} A_{k,i}(x) y^i (k=2..12) ===")
y = sp.Symbol('y')
for q_num, q_den in [(1, 2), (1, 3), (1, 4), (1, 5), (1, 10)]:
    q = sp.Rational(q_num, q_den)
    fails = []
    for k in range(2, 13):
        P = sum((q ** (i*(i-1)//2)) * A[k][i][0] * y**i for i in range(k + 1))
        Px = sum((q ** (i*(i-1)//2)) * A[k][i][1] * y**i for i in range(k + 1))
        R = sp.Poly(sp.expand(Px * sp.diff(P, y) - P * sp.diff(Px, y)), y, domain='QQ')
        rc = int(sp.count_roots(R))
        if rc > 0:
            fails.append((k, rc))
    print(f"  q={q_num}/{q_den}: Rayleigh at x=0 fails for k in {fails}")

print("\n=== 2. UHP STABILITY TEST FOR R_k(x, y) = sum q^{i(i-1)/2} A_{k,i}(x) y^i (k=2..8) ===")
rng = np.random.default_rng(99)
test_ys = [complex(re_y, im_y)
           for re_y in [-1000, -100, -20, -5, -1, -0.2, 0, 0.2, 1, 5, 20, 100]
           for im_y in [1e-4, 1e-2, 0.1, 1.0, 10.0]]
for _ in range(200):
    test_ys.append(complex(rng.normal(0, 1) * (10 ** rng.uniform(-2, 4)), 10 ** rng.uniform(-4, 3)))

for q in [0.5, 0.25, 0.1]:
    print(f"--- q = {q} ---")
    for k in range(2, 9):
        max_im = -1e9
        worst_y = None
        dmax = 2*k - 1
        for yv in test_ys:
            coeffs = np.zeros(dmax + 1, dtype=np.complex128)
            for i in range(k + 1):
                w = (q ** (i*(i-1)//2)) * (yv ** i)
                for m, c in enumerate(A[k][i]):
                    coeffs[m] += w * float(c)
            r = np.roots(coeffs[::-1])
            mi = float(np.max(r.imag))
            if mi > max_im:
                max_im = mi
                worst_y = yv
        print(f"  k={k}: max Im(x) = {max_im:+.3e} (at y={worst_y})")

print("\n=== 3. ALSO CHECK w=(k-i)! AND w=1/(k-i)! AT k=3,4 ===")
for k in [3, 4]:
    for name, wfn in [("(k-i)!", lambda i: math.factorial(k-i)), ("1/(k-i)!", lambda i: 1.0/math.factorial(k-i))]:
        max_im = -1e9
        dmax = 2*k - 1
        for yv in test_ys[:150]:
            coeffs = np.zeros(dmax + 1, dtype=np.complex128)
            for i in range(k + 1):
                w = wfn(i) * (yv ** i)
                for m, c in enumerate(A[k][i]):
                    coeffs[m] += w * float(c)
            r = np.roots(coeffs[::-1])
            mi = float(np.max(r.imag))
            if mi > max_im: max_im = mi
        print(f"  k={k} [{name}]: max Im(x) = {max_im:+.3e}")
