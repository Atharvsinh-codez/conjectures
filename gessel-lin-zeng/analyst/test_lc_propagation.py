"""CRITICAL TEST OF OUR MIRACLE INTERLACING LEMMA + LOG-CONCAVITY PROPAGATION:
We ALREADY PROVED BY HAND (5 lines!):
  If f_2 << f_1 << f_0 on (-1, 0) and f_1^2 >= f_0 * f_2 on (-1, 0),
  then for ANY D = t(1+t)d/dt + c + m*t (0 < c < m):
    (D f_2 + f_1) << (D f_1 + f_0)!!

NOW TEST:
1. Does row log-concavity P_{k,i}^2 >= P_{k,i-1} P_{k,i+1} on (-1, 0) hold for ALL wild (c_k, m_k, c'_k, m'_k)?!
2. Does U_{k,i}^2 >= U_{k,i-1} U_{k,i+1} on (-1, 0) hold for ALL wild (c_k, m_k, c'_k, m'_k)?!
3. Does the SINGLE-STEP operator D preserve row log-concavity on (-1, 0):
   For arbitrary f_2 << f_1 << f_0 on (-1,0) with f_1^2 >= f_0*f_2 on (-1,0), is (D f_1)^2 >= (D f_0)(D f_2) on (-1,0)?!
4. Does the QUADRATIC Q_{k,i}(t, y) = P_{k,i-1}(t) + 2 P_{k,i}(t) y + P_{k,i+1}(t) y^2 have ONLY REAL ROOTS in t
   for EVERY real y (and in y for EVERY real t, and Rayleigh >= 0) for GLZ k = 2..18 AND for wild m_k?!
   (Wait: check for lambda = 2 AND for general lambda! Note: if Q(t,y) = f_0 + 2 f_1 y + f_2 y^2 has all real t-roots
   for every real y in R, then for EVERY real y, D(Q(t, y)) = D f_0 + 2 (D f_1) y + (D f_2) y^2 ALSO has all real t-roots
   for every real y!Wait: is real-rooted in t for all real y SUFFICIENT to get (D f_1)^2 >= (D f_0)(D f_2) on (-1,0)?!)
"""
import math
import numpy as np
import sympy as sp

t = sp.Symbol('t')

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

def poly_sub(a, b):
    return poly_add(a, [-v for v in b])

def poly_mul(a, b):
    res = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            res[i + j] += u * v
    return res

def check_nonneg_m10(poly):
    P = sp.Poly(list(reversed(poly)), t, domain='QQ')
    # Check value at -1/2 and root count in (-1, 0)
    val = P.eval(sp.Rational(-1, 2))
    rc = int(sp.count_roots(P, -1, 0))
    rc_all = int(sp.count_roots(P))
    return (val > 0 and rc == 0), (val > 0 and rc_all == 0), rc, rc_all

print("=== 1. TEST ROW LOG-CONCAVITY OF P AND U FOR WILD SEQUENCES ===")
wild_tests = [
    ("GLZ k=1..6", [(1, k+1, 1, k+1) for k in range(1, 7)]),
    ("m=[2, 10^5, 2, 10^5, 2]", [(1, m, 1, m) for m in [2, 10**5, 2, 10**5, 2]]),
    ("m=[10^5, 2, 10^5, 2, 10^5]", [(1, m, 1, m) for m in [10**5, 2, 10**5, 2, 10**5]]),
    ("m_out=2, m_in=1000 (k=1..5)", [(1, 2, 1, 1000)]*5),
    ("m_out=1000, m_in=2 (k=1..5)", [(1, 1000, 1, 2)]*5),
    ("c=[1,100,1,100], m=[2,200,2,200]", [(1,2,1,2), (100,200,100,200), (1,2,1,2), (100,200,100,200)]),
]
for name, pairs in wild_tests:
    kmax = len(pairs)
    P = [[None]*(k+1) for k in range(kmax+1)]
    P[0][0] = [sp.Rational(1)]
    p_fail_m10 = []; p_fail_R = []; u_fail_m10 = []; u_fail_R = []
    for k in range(1, kmax + 1):
        co, mo, ci, mi = pairs[k - 1]
        U = [None] * (k + 1)
        U[0] = apply_D(ci, mi, P[k-1][0])
        for i in range(1, k):
            U[i] = poly_add(apply_D(ci, mi, P[k-1][i]), P[k-1][i-1])
        U[k] = P[k-1][k-1]
        for i in range(k + 1):
            P[k][i] = apply_D(co, mo, U[i])
        for i in range(1, k):
            dp = poly_sub(poly_mul(P[k][i], P[k][i]), poly_mul(P[k][i-1], P[k][i+1]))
            ok1, ok2, rc1, rc2 = check_nonneg_m10(dp)
            if not ok1: p_fail_m10.append((k, i, rc1))
            if not ok2: p_fail_R.append((k, i, rc2))
            du = poly_sub(poly_mul(U[i], U[i]), poly_mul(U[i-1], U[i+1]))
            ok1, ok2, rc1, rc2 = check_nonneg_m10(du)
            if not ok1: u_fail_m10.append((k, i, rc1))
            if not ok2: u_fail_R.append((k, i, rc2))
    print(f"  {name:35s}: P_m10={p_fail_m10}, P_R={p_fail_R} | U_m10={u_fail_m10}, U_R={u_fail_R}")

print("\n=== 2. TEST QUADRATIC STABILITY Q_{k,i}(t,y) = P_{k,i-1} + 2 P_{k,i} y + P_{k,i+1} y^2 FOR GLZ ===")
# For fixed real y, is t -> P_{k,i-1}(t) + 2 P_{k,i}(t) y + P_{k,i+1}(t) y^2 REAL-ROOTED in t?!
# (Even stronger: check Rayleigh Q_t Q_y - Q Q_{ty} >= 0 at t=0 and on R^2!)
pairs_glz = [(1, k+1, 1, k+1) for k in range(1, 15)]
P = [[None]*(k+1) for k in range(15)]
P[0][0] = [1]
for k in range(1, 15):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D(1, m, P[k-1][0])
    for i in range(1, k):
        U[i] = poly_add(apply_D(1, m, P[k-1][i]), P[k-1][i-1])
    U[k] = P[k-1][k-1]
    for i in range(k+1):
        P[k][i] = apply_D(1, m, U[i])

# Check Rayleigh of Q_{k,i}(t, y) = f0 + 2 f1 y + f2 y^2 at t = 0!!
# At t = 0: Q(0, y) = a0 + 2 a1 y + a2 y^2, Q_t(0, y) = b0 + 2 b1 y + b2 y^2.
# Wait: Rayleigh Q_t Q_y - Q Q_{ty} at t = 0 is a QUADRATIC in y:
#   2 [ (b0 a1 - a0 b1) + (b0 a2 - a0 b2) y + (b1 a2 - a1 b2) y^2 ]!
# It is >= 0 for ALL real y IFF (b0 a2 - a0 b2)^2 <= 4 (b0 a1 - a0 b1)(b1 a2 - a1 b2)!!
y = sp.Symbol('y')
for k in range(2, 15):
    fails_ray0 = []
    for i in range(1, k):
        a0, a1, a2 = P[k][i-1][0], P[k][i][0], P[k][i+1][0]
        b0, b1, b2 = P[k][i-1][1], P[k][i][1], P[k][i+1][1]
        c0 = b0*a1 - a0*b1
        c1 = b0*a2 - a0*b2
        c2 = b1*a2 - a1*b2
        disc = c1*c1 - 4*c0*c2
        if c0 <= 0 or c2 <= 0 or disc > 0:
            fails_ray0.append((i, c0 > 0, c2 > 0, float(disc) / (float(c1*c1) + 1)))
    print(f"  GLZ k={k:2d}: quadratic Q_{{k,i}}(0,y) Rayleigh fails = {fails_ray0[:5]}")

print("\n=== 3. TEST RANDOM ARBITRARY f_2 << f_1 << f_0 WITH f_1^2 >= f_0 f_2 ON (-1,0): DOES (D f_1)^2 >= (D f_0)(D f_2)? ===")
rng = np.random.default_rng(42)
fails_step2 = 0
for trial in range(500):
    # Pick deg f_2 = 2, f_1 = 3, f_0 = 4: roots in (-1, 0) interlacing!
    r0 = sorted(rng.uniform(-0.98, -0.02, 4))
    r1 = [rng.uniform(r0[m], r0[m+1]) for m in range(3)]
    r2 = [rng.uniform(r1[m], r1[m+1]) for m in range(2)]
    def from_roots(rs):
        p = np.array([1.0])
        for r in rs:
            # factor (1 - t/r) so p(0) = 1
            p = np.polymul(p, [-1.0/r, 1.0])
        return p[::-1] # low-to-high
    f0 = from_roots(r0)
    f1 = from_roots(r1)
    f2 = from_roots(r2)
    # Check if f1^2 >= f0*f2 on (-1, 0); if not, scale f1 up slightly so min(f1^2 - f0*f2) >= 0!
    ts = np.linspace(-0.999, -0.001, 400)
    v0 = np.polyval(f0[::-1], ts)
    v1 = np.polyval(f1[::-1], ts)
    v2 = np.polyval(f2[::-1], ts)
    # Find min scale s such that s^2 v1^2 >= v0 v2 everywhere on (-1,0)
    mask = (v0 * v2 > 0)
    if np.any(mask):
        need = np.sqrt(np.max((v0[mask] * v2[mask]) / (v1[mask]**2 + 1e-30)))
        if need > 1.0:
            f1 = f1 * (need * 1.01)
    # Now f1^2 > f0*f2 on (-1, 0)! Apply D = t(1+t)d/dt + 1 + m*t for m = 5!
    m = rng.uniform(1.5, 50.0)
    Df0 = apply_D(1.0, m, f0)
    Df1 = apply_D(1.0, m, f1)
    Df2 = apply_D(1.0, m, f2)
    dv0 = np.polyval(Df0[::-1], ts)
    dv1 = np.polyval(Df1[::-1], ts)
    dv2 = np.polyval(Df2[::-1], ts)
    if np.min(dv1*dv1 - dv0*dv2) < -1e-9 * np.max(dv1*dv1):
        fails_step2 += 1
print(f"  Random f_2 << f_1 << f_0 with f_1^2 >= f_0 f_2: Step 2 (Df_1)^2 >= (Df_0)(Df_2) fails = {fails_step2}/500")
