"""CRITICAL DIAGNOSTIC:
1. For r = 2, 3, 4, 5, 6, 7, 8, 9, 10, test Rayleigh at t=0 of
   Q^{(r)}_{k,i}(t, y) = sum_{j=0}^r binom(r, j) P_{k, i+j}(t) y^j
   for GLZ k up to 25! Find the EXACT smallest r that ever fails!
2. Test r = 2 (Q_{k,i} = P_{k,i-1} + 2 P_{k,i} y + P_{k,i+1} y^2) ADVERSARIALLY:
   - Exact Rayleigh P_t P_y - P P_{ty} >= 0 on ALL R^2 (not just t=0!) for GLZ k=2..12 AND wild sequences!
   - Exact Sturm check of t -> Q_{k,i}(t, y0) at negative y0 for GLZ k=2..10 AND wild sequences!
3. Test the 1-step closure for r = 2:
   Given an ARBITRARY 4-tuple (f0, f1, f2, f3) on (-1, 0) where BOTH Q_1 = f0 + 2 f1 y + f2 y^2
   and Q_2 = f1 + 2 f2 y + f3 y^2 are real-stable (and Q_2 << Q_1? test with & without Q_2 << Q_1!),
   is D Q_2 + Q_1 ALWAYS real-stable?!
"""
import math
import numpy as np
import sympy as sp

t = sp.Symbol('t')
y = sp.Symbol('y')

def apply_D(c_val, m_val, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c_val)
        res[r + 1] += coef * (r + m_val)
    return res

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def build_general(pairs):
    kmax = len(pairs)
    P = [[None]*(k+1) for k in range(kmax+1)]
    P[0][0] = [1]
    for k in range(1, kmax + 1):
        co, mo, ci, mi = pairs[k - 1]
        U = [None]*(k+1)
        U[0] = apply_D(ci, mi, P[k-1][0])
        for i in range(1, k):
            U[i] = poly_add(apply_D(ci, mi, P[k-1][i]), P[k-1][i-1])
        U[k] = P[k-1][k-1]
        for i in range(k+1):
            P[k][i] = apply_D(co, mo, U[i])
    return P

P_glz = build_general([(1, k+1, 1, k+1) for k in range(1, 21)])

print("=== 1. EXACT RAYLEIGH AT t=0 FOR DEGREE-r WINDOW Q^{(r)}_{k,i}(0, y) = sum_{j=0}^r binom(r,j) P_{k,i+j}(0) y^j ===")
for r_deg in range(2, 11):
    first_fail = None
    for k in range(r_deg, 21):
        for i in range(0, k - r_deg + 1):
            P0 = sum(math.comb(r_deg, j) * P_glz[k][i+j][0] * y**j for j in range(r_deg + 1))
            Pt0 = sum(math.comb(r_deg, j) * P_glz[k][i+j][1] * y**j for j in range(r_deg + 1))
            R = sp.Poly(sp.expand(Pt0 * sp.diff(P0, y) - P0 * sp.diff(Pt0, y)), y, domain='ZZ')
            if int(sp.count_roots(R)) > 0:
                first_fail = (k, i)
                break
        if first_fail: break
    print(f"  r = {r_deg:2d}: first Rayleigh fail at t=0 (k<=20) = {first_fail}")

print("\n=== 2. EXACT RAYLEIGH ON ALL R^2 FOR r=2 QUADRATIC Q_{k,i}(t,y) = P_{k,i-1} + 2 P_{k,i} y + P_{k,i+1} y^2 ===")
# For Q(t, y) = f0(t) + 2 f1(t) y + f2(t) y^2:
# Q_t Q_y - Q Q_{ty} = 2 [ W(f0, f1) + W(f0, f2) y + W(f1, f2) y^2 ]!
# WHERE W(f, g) = f' g - f g'!!
# WAIT A SECOND!! STOP EVERYTHING!! LOOK AT THAT IDENTITY!!
# Q_t Q_y - Q Q_{ty} >= 0 FOR ALL (t, y) IN R^2  <==>
#   W(f0, f1)(t) >= 0,  W(f1, f2)(t) >= 0,  AND  4 W(f0, f1)(t) W(f1, f2)(t) - W(f0, f2)(t)^2 >= 0  FOR ALL t IN R!!
# OH!! MY!! GOD!! COMPUTE 4 W(f0, f1) W(f1, f2) - W(f0, f2)^2 EXACTLY!!
def poly_deriv(p):
    return [m * c for m, c in enumerate(p)][1:] if len(p) > 1 else [0]
def poly_mul(a, b):
    res = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            res[i + j] += u * v
    return res
def poly_sub(a, b):
    return poly_add(a, [-v for v in b])
def Wronsk(f, g):
    return poly_sub(poly_mul(poly_deriv(f), g), poly_mul(f, poly_deriv(g)))

for k in range(2, 11):
    fails_m10 = []; fails_R = []
    for i in range(1, k):
        f0, f1, f2 = P_glz[k][i-1], P_glz[k][i], P_glz[k][i+1]
        W01 = Wronsk(f0, f1)
        W12 = Wronsk(f1, f2)
        W02 = Wronsk(f0, f2)
        D = poly_sub([4*c for c in poly_mul(W01, W12)], poly_mul(W02, W02))
        P_d = sp.Poly(list(reversed(D)), t, domain='ZZ')
        val_m = P_d.eval(sp.Rational(-1, 2))
        rc_m10 = int(sp.count_roots(P_d, -1, 0))
        rc_R = int(sp.count_roots(P_d))
        if val_m <= 0 or rc_m10 > 0: fails_m10.append((i, rc_m10, val_m > 0))
        if val_m <= 0 or rc_R > 0: fails_R.append((i, rc_R, val_m > 0))
    print(f"  GLZ k={k:2d}: 4*W01*W12 - W02^2 fails on (-1,0)={fails_m10}, on ALL R={fails_R}")

print("\n=== 3. SAME EXACT CHECK OF 4*W01*W12 - W02^2 ON WILD SEQUENCES!! ===")
wilds = [
    ("m=[2, 10^4, 2, 10^4, 2]", [(1, m, 1, m) for m in [2, 10**4, 2, 10**4, 2]]),
    ("m=[10^4, 2, 10^4, 2, 10^4]", [(1, m, 1, m) for m in [10**4, 2, 10**4, 2, 10**4]]),
    ("m_out=2, m_in=1000 (k=1..5)", [(1, 2, 1, 1000)]*5),
    ("m_out=1000, m_in=2 (k=1..5)", [(1, 1000, 1, 2)]*5),
    ("c=[1,100,1,100], m=[2,200,2,200]", [(1,2,1,2), (100,200,100,200), (1,2,1,2), (100,200,100,200)]),
]
for name, pairs in wilds:
    Pw = build_general(pairs)
    fails_m10 = []; fails_R = []
    for kk in range(2, len(pairs) + 1):
        for i in range(1, kk):
            f0, f1, f2 = Pw[kk][i-1], Pw[kk][i], Pw[kk][i+1]
            W01 = Wronsk(f0, f1); W12 = Wronsk(f1, f2); W02 = Wronsk(f0, f2)
            D = poly_sub([4*c for c in poly_mul(W01, W12)], poly_mul(W02, W02))
            P_d = sp.Poly(list(reversed(D)), t, domain='ZZ')
            val_m = P_d.eval(sp.Rational(-1, 2))
            rc_m10 = int(sp.count_roots(P_d, -1, 0))
            rc_R = int(sp.count_roots(P_d))
            if val_m <= 0 or rc_m10 > 0: fails_m10.append((kk, i, rc_m10))
            if val_m <= 0 or rc_R > 0: fails_R.append((kk, i, rc_R))
    print(f"  {name:35s}: fails(-1,0)={fails_m10}, fails(R)={fails_R}")
