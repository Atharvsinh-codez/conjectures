"""Diagnostic script:
1. Check log-space root margins of Pair (2,1) and summands (confirming whether 0.9999 was just geometric root clustering near 0!).
2. Check EXACTly at what r does D^r htilde_{k,i} << htilde_{k, i-r-1} fail!
3. Check if ANY power binom(k,i)^p makes sum binom(k,i)^p A_{k,i}(x) y^i pass Rayleigh at x=0 for k=10..15!
4. Check positivity of Term 1, Term 2, Term 3 of U_i^2 - U_{i-1} U_{i+1} on (-1, 0) AND check if M^2 - X*Y is ALWAYS > 0 on (-1, 0) up to k=16!
"""
import math
import numpy as np
import mpmath as mp
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

A = build_A(18)

# 1. Check log-space margins of Pair (2,1) at k=10, 14, 18!
def get_roots_x_mp(poly, dps=50):
    mp.mp.dps = dps
    d = len(poly) - 1
    if d == 0: return []
    g = (mp.mpf(poly[0]) / mp.mpf(poly[-1])) ** (mp.mpf(1) / d)
    scaled = [mp.mpf(c) * (g ** m) for m, c in enumerate(poly)]
    r = mp.polyroots(list(reversed(scaled)), maxsteps=200, extraprec=40)
    return sorted(mp.re(z) * g for z in r)

print("=== 1. LOG-SPACE (-log(-x)) FRACTIONS FOR PAIR (2,1) ===")
for k in [10, 14, 18]:
    for i in [2, k//2, k-1]:
        L = 3*k - 2 - i
        rB = [float(mp.log(-x)) for x in get_roots_x_mp(A[k-1][i-2], dps=45)] # decreasing!
        rS = [float(mp.log(-x)) for x in get_roots_x_mp(apply_S(L-1, A[k-1][i]), dps=45)]
        # rB[0] > rS[0] > rB[1] > ... (since x < 0, larger |x| = more negative x = larger log(-x))
        # wait: get_roots_x_mp returns ascending x (most negative first), so log(-x) is DESCENDING!
        fracs = [(rB[m] - rS[m]) / (rB[m] - rB[m+1]) for m in range(len(rS))]
        print(f"  k={k:2d}, i={i:2d}: log(-x) fracs min={min(fracs):.4f}, max={max(fracs):.4f}, "
              f"ends=[{fracs[0]:.3f}, {fracs[1]:.3f}, ..., {fracs[-2]:.3f}, {fracs[-1]:.3f}]")

# 2. Check Rayleigh at x=0 for w_i = binom(k,i)^p, p = 1, 2, 3, 4!
print("\n=== 2. RAYLEIGH AT x=0 FOR w_i = binom(k,i)^p (k=10..15) ===")
y = sp.Symbol('y')
for p_exp in [1, 2, 3, 4]:
    bad_ks = []
    for k in range(2, 16):
        P = sum((math.comb(k, i)**p_exp) * A[k][i][0] * y**i for i in range(k + 1))
        Px = sum((math.comb(k, i)**p_exp) * A[k][i][1] * y**i for i in range(k + 1))
        # Check if y_r'(0) = -Px(y_r)/P'(y_r) <= 0 at all roots of P!
        # Wait: P(y) = sum binom(k,i)^{p_exp+1} y^i. Let's check exact roots of R = Px*Py - P*Pxy!
        R = sp.Poly(sp.expand(Px * sp.diff(P, y) - P * sp.diff(Px, y)), y, domain='ZZ')
        rc = int(sp.count_roots(R))
        if rc > 0:
            bad_ks.append((k, rc))
    print(f"  p={p_exp}: Rayleigh fails at x=0 for k in {bad_ks[:8]}")

# 3. Check exact coefficient [x^1] A_{k,i}(x) / [x^0] A_{k,i}(x)!!
print("\n=== 3. RATIO [x^1] A_{k,i}(x) / [x^0] A_{k,i}(x) = sum_r (-1/r) ===")
for k in [4, 6, 10]:
    ratios = [sp.Rational(A[k][i][1], A[k][i][0]) for i in range(k)]
    print(f"  k={k}: a1/a0 = {[float(r) for r in ratios]}")
    print(f"       exact  = {ratios}")
