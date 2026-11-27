"""CRITICAL DIAGNOSTIC:
1. Check if prod_{j=1}^N (D_j + lam_j*y)(1) has ALL REAL ROOTS IN y for EVERY t in (-1, 0)!!
   (Remember: earlier we saw sum A_{k,i}(x) y^i FAILED in y at x=-0.01, k=3!
    Wait: THAT had D_j(D_j + y) [missing factor (1-x)^i? NO! In h_{k,i}(t), sum h_{k,i}(t) y^i = sum (1+t)^{2k-i} A_{k,i}(x) y^i
    = (1+t)^{2k} sum A_{k,i}(x) (y/(1+t))^i — so rescaling y by 1/(1+t) DOES NOT change reality of y-roots!!]
    SO sum h_{k,i}(t) y^i ALREADY has complex y-roots at t = -0.01, k = 3!!)
2. THEREFORE, test multiplying H_3(t, y) (which has COMPLEX y-roots at t = -0.01!) by (r + t + alpha*y)^b
   for b = 1, 2, 3, 4, 5, 6 and check if:
   (a) Row log-concavity F_i^2 - F_{i-1} F_{i+1} >= 0 on (-1, 0) BREAKS!
   (b) Real-rootedness of F_i(t) in (-1, 0) BREAKS!!
"""
import math
import numpy as np
import sympy as sp

t = sp.Symbol('t')
y = sp.Symbol('y')

def apply_GLZ(c, m, poly):
    d = len(poly) - 1
    res = [sp.Rational(0)] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * (r + m)
    return res

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [sp.Rational(0)] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def poly_mul(a, b):
    res = [sp.Rational(0)] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            res[i + j] += u * v
    return res

# 1. Check y-roots of (D_2 + y)^a(1) for t in (-1, 0)!
print("=== 1. CHECK y-ROOTS OF (D_{1,2} + y)^a(1) FOR t in (-1, 0) ===")
D2 = [[sp.Rational(1)]]
for _ in range(8): D2.append(apply_GLZ(1, 2, D2[-1]))
for a in range(2, 8):
    max_im = 0.0; worst_t = None
    for t0 in np.linspace(-0.99, -0.01, 200):
        ycoeffs = [math.comb(a, j) * sum(float(v) * (t0**m) for m, v in enumerate(D2[a - j])) for j in range(a + 1)]
        r = np.roots(ycoeffs[::-1])
        mi = float(np.max(np.abs(r.imag)))
        if mi > max_im: max_im = mi; worst_t = t0
    print(f"  (D_{{1,2}} + y)^{a}(1): max |Im(y)| on (-1,0) = {max_im:.4e} (at t={worst_t})")

# 2. Build GLZ htilde_{3, i}(t) (which HAS complex y-roots at t = -0.01!)
P = [[None]*(k+1) for k in range(6)]
P[0][0] = [sp.Rational(1)]
for k in range(1, 6):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_GLZ(1, m, P[k-1][0])
    for i in range(1, k):
        U[i] = poly_add(apply_GLZ(1, m, P[k-1][i]), P[k-1][i-1])
    U[k] = P[k-1][k-1]
    for i in range(k+1):
        P[k][i] = apply_GLZ(1, m, U[i])

# Check y-roots of sum P[3][i](t) y^i at t = -0.01!
t0 = -0.01
yc3 = [sum(float(v)*(t0**m) for m, v in enumerate(P[3][i])) for i in range(4)]
print(f"\n=== 2. y-ROOTS OF GLZ Htilde_3(-0.01, y): {np.roots(yc3[::-1])} ===")

# Now multiply Htilde_3(t, y) by (r + t + alpha*y)^b for various r in (0,1), alpha > 0, b = 1..6!
# WAIT: if Htilde_3(-0.01, y) has complex root y_0 = u + i*v, can multiplying by (r + t + alpha*y)^b
# make some coefficient F_i(t) have a COMPLEX root in t (or fail LC)?!
print("\n=== 3. ATTACK Htilde_3(t, y) * (r + t + alpha*y)^b FOR b = 1..8 ===")
rng = np.random.default_rng(123)
found_rr = None; found_lc = None
for b in [1, 2, 3, 4, 6, 8]:
    for trial in range(150):
        r = rng.uniform(0.001, 0.999)
        alpha = 10 ** rng.uniform(-3, 3)
        # Multiply polynomial list F[i](t) by (r + t + alpha*y) b times (in float64 first!)
        F = [np.array([float(v) for v in P[3][i]], dtype=np.float64) for i in range(4)]
        for _ in range(b):
            newF = [None] * (len(F) + 1)
            newF[0] = np.polymul(F[0][::-1], [1.0, r])[::-1]
            for i in range(1, len(F)):
                a1 = np.polymul(F[i][::-1], [1.0, r])
                a2 = alpha * F[i-1][::-1]
                newF[i] = np.polyadd(a1, a2)[::-1]
            newF[-1] = alpha * F[-1]
            F = newF
        # Check RR and LC on (-1, 0)
        for i in range(len(F)):
            rts = np.roots(F[i][::-1])
            if np.max(np.abs(rts.imag)) > 1e-4:
                found_rr = (b, r, alpha, i, float(np.max(np.abs(rts.imag))))
                break
        if found_rr: break
        ts = np.linspace(-0.99, -0.01, 200)
        for i in range(1, len(F) - 1):
            v0 = np.polyval(F[i-1][::-1], ts)
            v1 = np.polyval(F[i][::-1], ts)
            v2 = np.polyval(F[i+1][::-1], ts)
            if np.min(v1*v1 - v0*v2) < -1e-7 * np.max(v1*v1):
                found_lc = (b, r, alpha, i)
                break
    if found_rr: break
print("  found_rr fail:", found_rr)
print("  found_lc fail:", found_lc)
