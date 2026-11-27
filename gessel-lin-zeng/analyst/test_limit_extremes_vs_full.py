"""CRITICAL DIAGNOSTIC:
1. For (D_2 + alpha*y)^2 * D_M^6 * (D_2 + y)^2 (1), [y^2] = f + 4*alpha*g + alpha^2*h
   where f = D_2^2 D_M^6(1), g = D_2 D_M^6 D_2(1), h = D_M^6 D_2^2(1).
   Check:
   (a) Do f and h even interlace / have a common interlacer (does f + beta*h have complex roots)?!
   (b) What is the smallest c such that f + c*alpha*g + alpha^2*h is real-rooted for ALL alpha > 0?!
2. NOW TEST EVEN BIGGER SEPARATION:
   (D_2 + alpha*y)^r * D_M^s * (D_2 + y)^r (1) for r = 3, 4, 5, 6 and s = 6, 10, M = 10^4 (and M = 10^8)!
   Check [y^r] in EXACT / high-precision arithmetic!!
   WAIT A SECOND!! Look at the M -> inf limit of (D_2 + alpha*y)^r * D_M^s * (D_2 + y)^r (1)!!
   As M -> inf: D_M^s / M^s -> t^s (multiplication by t^s!!)!
   Wait: check! D_M(f) = t(1+t) f' + (1 + M t) f, so D_M^s(f) / M^s -> t^s * f!!
   THEREFORE, AS M -> inf:
     M^{-s} * (D_2 + alpha*y)^r * D_M^s * (D_2 + y)^r (1)
     ---> (D_2 + alpha*y)^r [ t^s * (D_2 + y)^r (1) ]!!
   HOLD EVERYTHING!! WE CAN COMPUTE THAT LIMIT *EXACTLY* IN 0.1 SECONDS FOR r = 2..8, s = 2..20!!
   TEST THAT LIMIT RIGHT NOW!!
"""
import math
import numpy as np
import sympy as sp

t = sp.Symbol('t')

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

def count_real_with_mult(P):
    total = 0
    for fac, mult in P.sqf_list()[1]:
        total += mult * int(sp.count_roots(fac))
    return total

print("=== 1. TEST M->inf LIMIT: [y^r] (D_2 + alpha*y)^r [ t^s * (D_2 + y)^r (1) ] ===")
# Precompute D_2^j(1) and D_2^j(t^s * D_2^l(1))
# [y^r] = sum_{j=0}^r binom(r, j)^2 * alpha^j * D_2^{r-j} ( t^s * D_2^j (1) )!!
D2_1 = [[sp.Rational(1)]]
for _ in range(10): D2_1.append(apply_GLZ(1, 2, D2_1[-1]))

for r in [2, 3, 4, 5, 6]:
    for s_pow in [4, 10, 20]:
        # Build terms T_j = D_2^{r-j}( t^s * D_2^j(1) )
        terms = []
        for j in range(r + 1):
            p = [sp.Rational(0)] * s_pow + list(D2_1[j])
            for _ in range(r - j):
                p = apply_GLZ(1, 2, p)
            terms.append(p)
        # First check if extremes T_0 + beta * T_r fail real-rootedness!
        ext_fail = None
        for exp_b in range(-12, 13, 2):
            beta = sp.Rational(10)**exp_b
            comb = poly_add(terms[0], [beta * v for v in terms[r]])
            # strip t factors
            while len(comb) > 1 and comb[0] == 0: comb = comb[1:]
            P = sp.Poly(list(reversed(comb)), t, domain='QQ')
            rc = count_real_with_mult(P)
            if rc < P.degree():
                ext_fail = (exp_b, rc, P.degree())
                break
        # Now check FULL binomial^2 sum: sum_{j=0}^r binom(r, j)^2 * alpha^j * terms[j]!!
        full_fail = None
        for exp_a in range(-6, 7):
            alpha = sp.Rational(10)**exp_a
            comb = [sp.Rational(0)]
            for j in range(r + 1):
                w = (math.comb(r, j) ** 2) * (alpha ** j)
                comb = poly_add(comb, [w * v for v in terms[j]])
            while len(comb) > 1 and comb[0] == 0: comb = comb[1:]
            P = sp.Poly(list(reversed(comb)), t, domain='QQ')
            rc = count_real_with_mult(P)
            if rc < P.degree():
                full_fail = (exp_a, rc, P.degree())
                break
        print(f"  r={r}, s={s_pow:2d}: extremes T_0+beta*T_r FAIL={ext_fail} | FULL sum FAIL={full_fail}")
