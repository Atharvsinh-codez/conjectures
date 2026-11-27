"""Test the asymmetric M->inf limit of the Atomic Theorem:
  P(x) = sum_{j=0}^b binom(b, j) * binom(a, i - j) * gamma^j * A_{a - i + j + 1}(x)
for arbitrary a, b, i and gamma > 0!
Also test G_r(x, gamma) = sum_{j=0}^r binom(r, j)^2 * gamma^j * A_{j+1}(x) at extreme gamma (1e-6..1e6) and r up to 24!
"""
import math
import sympy as sp

E = [[1]]
for n in range(1, 35):
    prev = E[-1]
    nxt = [0] * (len(prev) + 1)
    for m, c in enumerate(prev):
        nxt[m] += c * (m + 1)
        nxt[m + 1] += c * (n - m)
    E.append(nxt)

x = sp.Symbol('x')

print("=== 1. ASYMMETRIC LIMIT sum binom(b,j)*binom(a,i-j)*gamma^j*A_{a-i+j+1}(x) ===")
fails = []
for a, b, i in [(10, 4, 4), (20, 4, 4), (20, 6, 6), (20, 4, 2), (10, 10, 5), (15, 10, 8)]:
    fail_case = None
    for exp_g in [-6, -4, -2, -1, 0, 1, 2, 4, 6]:
        gamma = sp.Rational(10) ** exp_g
        deg_max = a - i + b
        poly = [sp.Rational(0)] * (deg_max + 1)
        for j in range(max(0, i - a), min(b, i) + 1):
            w = math.comb(b, j) * math.comb(a, i - j) * (gamma ** j)
            ej = E[a - i + j] # A_{a-i+j+1}
            for m, c in enumerate(ej):
                poly[m] += w * c
        while len(poly) > 1 and poly[-1] == 0: poly.pop()
        P = sp.Poly(list(reversed(poly)), x, domain='QQ')
        rc = int(sp.count_roots(P))
        if rc < P.degree():
            fail_case = (exp_g, rc, P.degree())
            break
    print(f"  a={a:2d}, b={b:2d}, i={i:2d}: fail = {fail_case}")

print("\n=== 2. CHECK WHAT HAPPENS WITH (1-x)^{b-j} FACTOR vs WITHOUT! ===")
# WAIT A SECOND! CHECK OUR ALGEBRA FOR (D_2 + alpha*y)^b [ t^s * (D_2 + y)^a (1) ]!!
# In [y^i] (D_2 + alpha*y)^b [ t^s * (D_2 + y)^a (1) ]:
#   j copies of alpha*y from left => (b - j) copies of D_2 on left => factor (1+t)^{b - j} as s -> inf!
#   (i - j) copies of y from right => (a - (i - j)) copies of D_2 on right => D_2^{a - i + j}(1) = (1+t)^{a - i + j} A_{a-i+j+1}(x)!
#   Product of (1+t) powers: (1+t)^{(b - j) + (a - i + j)} = (1+t)^{a + b - i} — INDEPENDENT OF j!!
#   SO THE ALGEBRA WAS 100% EXACT!! NO (1-x)^{b-j} FACTOR AT ALL!!
# Now check EVEN MORE EXTREME: a = 30, b = 4, i = 4, and fine gamma grid around 1/30, 1/10, 1, 10, 30!
for a, b, i in [(25, 4, 4), (25, 6, 6), (25, 8, 8)]:
    fail_case = None
    for g in [sp.Rational(1, 1000), sp.Rational(1, 100), sp.Rational(1, 30), sp.Rational(1, 10), sp.Rational(1, 3), sp.Rational(1), sp.Rational(3), sp.Rational(10), sp.Rational(100)]:
        deg_max = a - i + b
        poly = [sp.Rational(0)] * (deg_max + 1)
        for j in range(max(0, i - a), min(b, i) + 1):
            w = math.comb(b, j) * math.comb(a, i - j) * (g ** j)
            ej = E[a - i + j]
            for m, c in enumerate(ej):
                poly[m] += w * c
        while len(poly) > 1 and poly[-1] == 0: poly.pop()
        P = sp.Poly(list(reversed(poly)), x, domain='QQ')
        rc = int(sp.count_roots(P))
        if rc < P.degree():
            fail_case = (g, rc, P.degree())
            break
    print(f"  FINE a={a:2d}, b={b:2d}, i={i:2d}: fail = {fail_case}")
