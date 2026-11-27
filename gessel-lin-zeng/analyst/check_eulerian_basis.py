"""3-second miracle check:
1. Does the linear operator E_1(z^n) = A_{n+1}(x) / n! (and E_2(z^n) = A_{n+1}(x) / (n+1)!) preserve RZ(-inf, 0]?
   Test on (1 + z)^d and random real-negative-rooted polynomials!
2. Expand GLZ A_{k, i}(x) in the basis { A_{n+1}(x) / n! } and { A_{n+1}(x) / (n+1)! } and { A_{n+1, n+1}(x) / n! }!
   CHECK IF THE COEFFICIENT POLYNOMIAL sum alpha_n z^n IS REAL-ROOTED!!
3. ALSO check the expansion of htilde_{k, i}(t) in the basis { D_2^n(1) / n! } and { D_1^n(1) / n! }!!
"""
import math
import sympy as sp

x = sp.Symbol('x')
z = sp.Symbol('z')

# Build Eulerian polynomials E[n] = A_{n+1}(x) (deg n), n = 0..20
E = [[1]]
for n in range(1, 22):
    prev = E[-1]
    nxt = [0] * (len(prev) + 1)
    for m, c in enumerate(prev):
        nxt[m] += c * (m + 1)
        nxt[m + 1] += c * (n - m)
    E.append(nxt)

# Build GLZ A_{k,i}(x)
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

A = [[None]*(k+1) for k in range(10)]
A[0][0] = [1]
for k in range(1, 10):
    A[k][0] = apply_S(3*k - 2, apply_S(3*k - 3, A[k-1][0]))
    for i in range(1, k):
        L = 3*k - 2 - i
        A[k][i] = apply_S(L, poly_add(apply_S(L - 1, A[k-1][i]), A[k-1][i-1]))
    A[k][k] = apply_S(2*k - 2, A[k-1][k-1])

print("=== 1. TEST IF E_1(z^n) = A_{n+1}(x)/n! AND E_2(z^n) = A_{n+1}(x)/(n+1)! PRESERVE RZ(-inf,0] ON (1+z)^d ===")
for label, denom_fn in [("A_{n+1}(x)/n!", lambda n: math.factorial(n)),
                        ("A_{n+1}(x)/(n+1)!", lambda n: math.factorial(n + 1))]:
    fails = []
    for d in range(2, 12):
        for gamma in [sp.Rational(1, 10), sp.Rational(1), sp.Rational(10)]:
            poly = [sp.Rational(0)] * (d + 1)
            for n in range(d + 1):
                w = sp.Rational(math.comb(d, n), denom_fn(n)) * (gamma ** n)
                for m, c in enumerate(E[n]):
                    poly[m] += w * c
            P = sp.Poly(list(reversed(poly)), x, domain='QQ')
            if int(sp.count_roots(P)) < P.degree():
                fails.append((d, gamma, int(sp.count_roots(P)), P.degree()))
                break
    print(f"  {label:20s}: fails on (1+gamma*z)^d (d<=11) = {fails[:4]}")

def expand_in_basis(target, basis):
    """Given target poly of deg d and basis[0..d] where deg(basis[n]) = n, solve triangular system."""
    d = len(target) - 1
    rem = [sp.Rational(c) for c in target]
    coeffs = [sp.Rational(0)] * (d + 1)
    for n in range(d, -1, -1):
        c = rem[n] / basis[n][n]
        coeffs[n] = c
        for m in range(n + 1):
            rem[m] -= c * basis[n][m]
    return coeffs

print("\n=== 2. EXPAND A_{k, i}(x) IN EULERIAN BASIS { A_{n+1}(x) / n! } AND CHECK COEFF POLYNOMIAL ROOTS ===")
basis_E1 = [[sp.Rational(c, math.factorial(n)) for c in E[n]] for n in range(20)]
basis_E2 = [[sp.Rational(c, math.factorial(n + 1)) for c in E[n]] for n in range(20)]
# Also second-order Eulerian basis: A[n+1][n+1] has deg n!
basis_S2 = [[sp.Rational(c, math.factorial(n)) for c in A[n + 1][n + 1]] for n in range(9)]

for bname, bas in [("E1: A_{n+1}/n!", basis_E1), ("E2: A_{n+1}/(n+1)!", basis_E2), ("S2: A_{n+1,n+1}/n!", basis_S2)]:
    print(f"\n--- Basis {bname} ---")
    for k in range(2, 6):
        for i in range(k + 1):
            if len(A[k][i]) > len(bas): continue
            cs = expand_in_basis(A[k][i], bas)
            all_pos = all(c >= 0 for c in cs)
            P = sp.Poly(list(reversed(cs)), z, domain='QQ')
            rc = int(sp.count_roots(P)) if P.degree() > 0 else 0
            if k <= 4 or rc < P.degree() or not all_pos:
                print(f"  k={k}, i={i}: all_nonneg={str(all_pos):5s}, real_roots={rc}/{P.degree()}, coeffs={[float(c) for c in cs[:5]]}")
