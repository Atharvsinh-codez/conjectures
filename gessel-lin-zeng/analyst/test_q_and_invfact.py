"""Test:
1. Exact Rayleigh at x=0 and full R^2 for w = 1/(k-i)! (k=2..14) — confirm whether it fails at larger k!
2. Adversarial test of R_{k,q}(x, y) = sum q^{i(i-1)/2} A_{k,i}(x) y^i (and Q_{k,q}(t, y)) on ALL R^2
   (including x=0, x<0, x>0, Rayleigh Px*Py - P*Pxy, and positive lines) for k=2..12 and q = 1/2, 1/3, 1/4, 1/10!
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

A = build_A(14)
y = sp.Symbol('y')

print("=== 1. EXACT RAYLEIGH AT x=0 FOR w_i = 1/(k-i)! (k=2..14) ===")
for k in range(2, 15):
    P = sum(sp.Rational(A[k][i][0], math.factorial(k - i)) * y**i for i in range(k + 1))
    Px = sum(sp.Rational(A[k][i][1], math.factorial(k - i)) * y**i for i in range(k + 1))
    R = sp.Poly(sp.expand(Px * sp.diff(P, y) - P * sp.diff(Px, y)), y, domain='QQ')
    rc = int(sp.count_roots(R))
    if rc > 0 or k in [2, 5, 8, 10, 14]:
        print(f"  k={k:2d}: Rayleigh real roots at x=0 = {rc}")

print("\n=== 2. ADVERSARIAL RAYLEIGH TEST ON R^2 FOR R_{k,q}(x, y) = sum q^{i(i-1)/2} A_{k,i}(x) y^i ===")
def eval_poly_and_deriv(coeffs, x0):
    val = 0.0; dval = 0.0
    for c in reversed(coeffs):
        dval = dval * x0 + val
        val = val * x0 + float(c)
    return val, dval

for q in [0.5, 0.25, 0.1]:
    print(f"--- q = {q} ---")
    for k in range(2, 13):
        w = [q ** (i*(i-1)//2) for i in range(k + 1)]
        min_rel = (1e300, None)
        # Test x across (-10, 10) and y near roots (-q^{-j} * scale)!
        xs = [-5.0, -2.0, -1.0, -0.5, -0.2, -0.05, 0.0, 0.1, 0.5, 1.0, 2.0, 5.0]
        # y grid spanning exponential range up to q^{-k}
        ys = []
        for e in range(-2, k + 2):
            base = (1.0 / q) ** e
            for mult in [-3.0, -1.5, -1.0, -0.7, -0.5, -0.3, 0.5, 1.0]:
                ys.append(mult * base)
        for x0 in xs:
            vals = [eval_poly_and_deriv(A[k][i], x0) for i in range(k + 1)]
            for y0 in ys:
                # Compute in log-safe scaling
                terms = [w[i] * (y0 ** i) for i in range(k + 1)]
                scale = max(abs(terms[i] * vals[i][0]) for i in range(k + 1)) + 1e-300
                P = sum(terms[i] * vals[i][0] / scale for i in range(k + 1))
                Px = sum(terms[i] * vals[i][1] / scale for i in range(k + 1))
                Py = sum((i * w[i] * (y0 ** (i-1)) if i >= 1 else 0.0) * vals[i][0] / scale for i in range(k + 1))
                Pxy = sum((i * w[i] * (y0 ** (i-1)) if i >= 1 else 0.0) * vals[i][1] / scale for i in range(k + 1))
                ray = Px * Py - P * Pxy
                rel = ray / (abs(Px * Py) + abs(P * Pxy) + 1e-300)
                if rel < min_rel[0]:
                    min_rel = (rel, (x0, y0))
        print(f"  k={k:2d}: min norm Rayleigh = {min_rel[0]:+.4e} at x={min_rel[1][0]}, y={min_rel[1][1]:.3e}")
