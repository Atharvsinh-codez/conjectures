"""2-second check: Is sum_{j=0}^r binom(r, j)^2 * gamma^j * A_{j+1}(x) ALWAYS real-rooted for r up to 25 and all gamma > 0?!
(Where A_{j+1}(x) is the Eulerian polynomial: A_1=1, A_2=1+x, A_3=1+4x+x^2, ...)
Also check WITHOUT the square: sum_{j=0}^r binom(r, j) * gamma^j * A_{j+1}(x)!
"""
import math
import numpy as np
import sympy as sp

# Build Eulerian polynomials A_1, A_2, ..., A_{26} (low-to-high)
# A_1 = [1], A_{n+1}(x) = (1 + n x) A_n(x) + x(1-x) A_n'(x)
E = [[1]] # E[0] = A_1 = 1
for n in range(1, 26):
    prev = E[-1]
    nxt = [0] * (len(prev) + 1)
    for m, c in enumerate(prev):
        nxt[m] += c * (m + 1)
        nxt[m + 1] += c * (n - m)
    E.append(nxt)

print("Check E[0..3]:", E[:4])
x = sp.Symbol('x')

for mode, power in [("binom(r,j)^2", 2), ("binom(r,j)^1", 1)]:
    print(f"\n=== TEST {mode} * gamma^j * A_{{j+1}}(x) ===")
    for r in [4, 8, 12, 16, 20]:
        fail = None
        for gamma in [sp.Rational(1, 100), sp.Rational(1, 10), sp.Rational(1, 3), sp.Rational(1), sp.Rational(3), sp.Rational(10), sp.Rational(100)]:
            poly = [sp.Rational(0)] * (r + 1)
            for j in range(r + 1):
                w = (math.comb(r, j) ** power) * (gamma ** j)
                for m, c in enumerate(E[j]):
                    poly[m] += w * c
            P = sp.Poly(list(reversed(poly)), x, domain='QQ')
            rc = int(sp.count_roots(P))
            if rc < P.degree():
                fail = (gamma, rc, P.degree())
                break
        print(f"  r={r:2d}: fail = {fail}")
