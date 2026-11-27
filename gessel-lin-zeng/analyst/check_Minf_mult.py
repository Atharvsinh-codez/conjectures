"""Check whether test_Minf_limit.py had COMPLEX roots or only multiple root at t=0,
AND test FINITE M = 10^6 with a=2, b=4, i=1, alpha=1/100 AND also test b >= 2 with complex-root check (counting multiplicity)!"""
import math
import sympy as sp

t = sp.Symbol('t')

def apply_GLZ(c, m, poly):
    d = len(poly) - 1
    res = [sp.Rational(0)] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * (r + m)
    return res

def apply_Lag(c, m, poly):
    d = len(poly) - 1
    res = [sp.Rational(0)] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * m
    return res

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [sp.Rational(0)] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def count_real_with_mult(P):
    """Exact count of real roots WITH multiplicity!"""
    total = 0
    for fac, mult in P.sqf_list()[1]:
        total += mult * int(sp.count_roots(fac))
    return total

print("=== CHECK sp.count_roots(t^3) BEHAVIOR ===")
print("  sp.count_roots(Poly(t^3)) =", sp.count_roots(sp.Poly(t**3, t)))
print("  sp.count_roots(Poly((t+1)*t^3)) =", sp.count_roots(sp.Poly((t+1)*t**3, t)))

# Now properly re-test M -> inf limit COUNTING MULTIPLICITY (and check finite M=10^6)!
D2_pows = [[sp.Rational(1)]]
L_pows = [[sp.Rational(1)]]
for j in range(1, 14):
    D2_pows.append(apply_GLZ(1, 2, D2_pows[-1]))
    L_pows.append(apply_Lag(1, 1, L_pows[-1]))

for label, pows in [("Laguerre (D_1 + y)^a", L_pows), ("GLZ (D_{1,2} + y)^a", D2_pows)]:
    fails = []
    for a in [2, 3, 4, 6, 8, 10]:
        for b in [2, 3, 4, 6, 8, 10]:
            for i in range(1, a + b):
                for alpha in [sp.Rational(1, 100), sp.Rational(1, 10), sp.Rational(1), sp.Rational(10), sp.Rational(100)]:
                    acc = [sp.Rational(0)]
                    for j in range(max(0, i - b), min(a, i) + 1):
                        coeff = math.comb(b, i - j) * (alpha ** (i - j)) * math.comb(a, j)
                        tpow = b - (i - j)
                        term = [sp.Rational(0)] * tpow + [coeff * v for v in pows[a - j]]
                        acc = poly_add(acc, term)
                    P = sp.Poly(list(reversed(acc)), t, domain='QQ')
                    # strip t factor!
                    while len(acc) > 1 and acc[0] == 0:
                        acc = acc[1:]
                    Q = sp.Poly(list(reversed(acc)), t, domain='QQ')
                    if Q.degree() > 0:
                        rc = count_real_with_mult(Q)
                        if rc < Q.degree():
                            fails.append((a, b, i, alpha, rc, Q.degree()))
                            break
                if fails: break
            if fails: break
        if fails: break
    print(f"  {label}: TRUE complex-root fails = {fails}")
