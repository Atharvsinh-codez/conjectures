"""M -> infinity limit test of the Atomic Theorem for BOTH Laguerre AND GLZ!!
1. Laguerre limit: P(t) = [y^i] ( (t + alpha*y)^b * (D_1 + y)^a (1) ) where D_1 = t*d/dt + 1 + t.
   Check if P(t) is ALWAYS real-rooted for ALL a, b <= 12 and ALL alpha > 0 (e.g. alpha = 0.001..1000)!
2. GLZ limit: For D_M = t(1+t)d/dt + c + M*t with M -> inf, c/M -> r in (0, 1), lam/M -> alpha > 0:
   (D_M + lam*y) / M -> (r + t + alpha*y)!!
   So the M -> inf limit of the GLZ Atomic Theorem is:
     P(t) = [y^i] ( prod_{s=1}^b (r_s + t + alpha_s * y) * prod_{j=1}^a (D_j + y)(1) )!!
   TEST THIS FOR r_s in (0, 1), alpha_s > 0, a, b <= 10!
   WAIT A SECOND!! Look at what happens if you set r_s = 1/2 (so r_s + t has root -1/2 in (-1, 0)):
     Does [y^i] ( (1/2 + t + alpha*y)^b * (D_1 + y)^a (1) ) ALWAYS have all real roots in (-1, 0)?!
   WAIT A SECOND!! WHAT ABOUT AT b = 1:
     [y^i] ( (r + t + alpha*y) * sum F_j y^j ) = (r + t) F_i + alpha F_{i-1}!!
   CHECK IF (r + t) F_i + alpha F_{i-1} HAS ALL REAL ROOTS FOR EVERY r in (0, 1) AND alpha > 0!!
   WAIT A SECOND!! STOP EVERYTHING!!
   If F_i << F_{i-1} (deg F_{i-1} = deg F_i + 1), DOES (r + t) F_i + alpha F_{i-1} ALWAYS HAVE ALL REAL ROOTS?!
   NO!! Because (r + t) adds an ARBITRARY root -r in (-1, 0), which does NOT have to interlace F_{i-1}!!
   TEST THAT RIGHT THIS SECOND!!
"""
import math
import numpy as np
import sympy as sp

t = sp.Symbol('t')

def apply_Lag(c, m, poly):
    d = len(poly) - 1
    res = [sp.Rational(0)] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * m
    return res

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

# 1. Test (r + t) F_{a, i}(t) + alpha * F_{a, i-1}(t) for GLZ F_{a, j} = binom(a, j) D_2^{a-j}(1)!!
# THIS IS LITERALLY THE M -> inf LIMIT OF ONE EXTRA FACTOR (D_{r*M, M} + alpha*M*y) ON TOP OF (D_{1, 2} + y)^a (1)!!
print("=== 1. TEST GLZ M->inf LIMIT: (r + t) F_{a,i} + alpha * F_{a, i-1} ===")
D2_pows = [[sp.Rational(1)]]
for j in range(1, 10):
    D2_pows.append(apply_GLZ(1, 2, D2_pows[-1]))

# Try a = 4, i = 2: F_{4,2} = 6 D_2^2(1), F_{4,1} = 4 D_2^3(1)
# Check roots of D_2^2(1) and D_2^3(1), then place -r in (-1, 0) OUTSIDE the interlacing gap!!
for a in [2, 3, 4, 5]:
    for i in range(1, a):
        Fi = [math.comb(a, i) * v for v in D2_pows[a - i]]
        Fim1 = [math.comb(a, i - 1) * v for v in D2_pows[a - i + 1]]
        # Test r in {1/100, 1/10, 1/2, 9/10, 99/100} and alpha in {1/100, 1/10, 1, 10, 100}
        found = None
        for r in [sp.Rational(1, 100), sp.Rational(1, 10), sp.Rational(1, 2), sp.Rational(9, 10), sp.Rational(99, 100)]:
            for alpha in [sp.Rational(1, 100), sp.Rational(1, 10), sp.Rational(1), sp.Rational(10), sp.Rational(100)]:
                comb = poly_add(poly_mul([r, sp.Rational(1)], Fi), [alpha * v for v in Fim1])
                P = sp.Poly(list(reversed(comb)), t, domain='QQ')
                rc = int(sp.count_roots(P, -1, 0))
                if rc < P.degree():
                    found = (r, alpha, rc, P.degree())
                    break
            if found: break
        if found:
            print(f"  GLZ M->inf COUNTEREXAMPLE at a={a}, i={i}: r={found[0]}, alpha={found[1]}, rc={found[2]} < deg={found[3]}!!")
            # ALSO check FINITE M = 10^6!!
            M = 10**6
            r_v, al_v = found[0], found[1]
            c_big = r_v * M
            lam_big = al_v * M
            comb_fin = poly_add(apply_GLZ(c_big, M, Fi), [lam_big * v for v in Fim1])
            Pf = sp.Poly(list(reversed(comb_fin)), t, domain='QQ')
            rcf = int(sp.count_roots(Pf, -1, 0))
            print(f"    FINITE M=10^6 CHECK: rc={rcf} < deg={Pf.degree()}!!")
            break
    if found: break

# 2. ALSO TEST LAGUERRE M->inf LIMIT: (t + alpha*y)^b * (D_1 + y)^a (1)!!
# Note: in Laguerre, (r + t) has r = 0 (root at 0, NOT arbitrary -r in (-1,0)!), so t * F_i ALWAYS interlaces F_{i-1}!!
print("\n=== 2. TEST LAGUERRE M->inf LIMIT: [y^i] (t + alpha*y)^b (D_1 + y)^a (1) ===")
L_pows = [[sp.Rational(1)]]
for j in range(1, 12):
    L_pows.append(apply_Lag(1, 1, L_pows[-1]))

lag_fails = []
for a in [2, 4, 6, 8]:
    for b in [1, 2, 4, 6, 8]:
        for i in range(1, a + b):
            for alpha in [sp.Rational(1, 100), sp.Rational(1, 5), sp.Rational(1), sp.Rational(5), sp.Rational(100)]:
                # sum_j binom(b, i-j) alpha^{i-j} t^{b-(i-j)} * binom(a, j) L_pows[a-j]
                acc = [sp.Rational(0)]
                for j in range(max(0, i - b), min(a, i) + 1):
                    coeff = math.comb(b, i - j) * (alpha ** (i - j)) * math.comb(a, j)
                    tpow = b - (i - j)
                    term = [sp.Rational(0)] * tpow + [coeff * v for v in L_pows[a - j]]
                    acc = poly_add(acc, term)
                P = sp.Poly(list(reversed(acc)), t, domain='QQ')
                if P.degree() > 0 and int(sp.count_roots(P)) < P.degree():
                    lag_fails.append((a, b, i, alpha, int(sp.count_roots(P)), P.degree()))
                    break
            if lag_fails: break
        if lag_fails: break
    if lag_fails: break
print(f"  Laguerre M->inf limit fails: {lag_fails}")

# 3. ALSO TEST GLZ WITH FIXED c_j = 1 (so r = c/M -> 0, SAME AS LAGUERRE!!):
# In GLZ, c_j = 1 is FIXED! So as M -> inf, (D_{1, M} + alpha*M*y)/M -> (t + alpha*y), root at 0 (NOT arbitrary -r)!
print("\n=== 3. TEST GLZ WITH FIXED c_j = 1, M->inf LIMIT: [y^i] (t + alpha*y)^b (D_{1,2} + y)^a (1) ===")
glz_c1_fails = []
for a in [2, 4, 6, 8]:
    for b in [1, 2, 4, 6]:
        for i in range(1, a + b):
            for alpha in [sp.Rational(1, 100), sp.Rational(1, 5), sp.Rational(1), sp.Rational(5), sp.Rational(100)]:
                acc = [sp.Rational(0)]
                for j in range(max(0, i - b), min(a, i) + 1):
                    coeff = math.comb(b, i - j) * (alpha ** (i - j)) * math.comb(a, j)
                    tpow = b - (i - j)
                    term = [sp.Rational(0)] * tpow + [coeff * v for v in D2_pows[a - j]]
                    acc = poly_add(acc, term)
                P = sp.Poly(list(reversed(acc)), t, domain='QQ')
                if P.degree() > 0 and int(sp.count_roots(P, -1, 0)) < P.degree():
                    # Note: t^tpow has root at 0, check total real roots!
                    if int(sp.count_roots(P)) < P.degree():
                        glz_c1_fails.append((a, b, i, alpha, int(sp.count_roots(P)), P.degree()))
                        break
            if glz_c1_fails: break
        if glz_c1_fails: break
    if glz_c1_fails: break
print(f"  GLZ (c=1) M->inf limit fails: {glz_c1_fails}")
