"""
Initial exact & analytic exploration of Bardakov-Kozlovskaya Conjecture 3.3(2):
For even n = 2m >= 2:
  p_n(x,y,z) = prod_{k=0}^{n-1} (z - (alpha + zeta^k beta)^n)
             = sum_{a+2b+3c=n} c_{a,b,c}^{(n)} e_1^a e_2^b e_3^c
  where alpha^n = x, beta^n = y, zeta = exp(2*pi*i/n).

1. Verify exact Toeplitz determinant identity:
   p_n(x,y,z) = det M_n(x,y,z), where
     M_{ij} = (z - x - y) if i==j, -binom(n, j-i)*x if j>i, -binom(n, n+j-i)*y if j<i.
2. Compute exact c_{a,b,c}^{(2m)} for 2m = 2, 4, 6, 8, 10, 12 using fast bivariate evaluation!
   Wait: to get ALL c_{a,b,c}^{(n)} from det M_n(x, y, z):
   Set x + y + z = 1 (homogeneity!), OR even better:
   Set x = u, y = v, z = 1 - u - v!
   Wait: in terms of (u, v) with x=u, y=v, z=1-u-v:
   e_1 = 1, e_2 = u*v + (u+v)*(1-u-v), e_3 = u*v*(1-u-v)!
   EVEN BETTER: set x = 1 (no, set x, y, z to roots of T^3 - e_1 T^2 + e_2 T - e_3 = 0,
   OR set z = 1 and let s = x+y, p = xy! Then:
     e_1 = s + 1, e_2 = s + p, e_3 = p!
   LOOK HOW CLEAN THAT IS!!
     p = e_3, s = e_2 - e_3, e_1 = 1 + e_2 - e_3!!
   AND in M_n(x, y, z) with z = 1, x+y = s, xy = p:
   WAIT A SECOND!! Can we compute p_n(x, y, 1) as a polynomial in (s, p) directly?!
   YES! Since p_n(x, y, 1) is symmetric in (x, y) of degree n:
   Set s = e_2 - e_3, p = e_3? NO, even simpler:
   Set x = t, y = p/t, z = 1 - t - p/t? No, polynomial interpolation in (e_2, e_3) with e_1 = 1!
   Wait: with z = 1, x+y = s, xy = p, we just solve x = (s + sqrt(s^2-4p))/2, y = (s - sqrt(s^2-4p))/2
   — EVEN CLEANER: set x = u, y = v, z = 1, compute p_n(u, v, 1) in SymPy/exact ints,
   then rewrite in e_1 = u+v+1, e_2 = uv+u+v, e_3 = uv!
   WAIT: EVEN SIMPLER!! Set u = r, v = 1 -> no, set x = 1, y = u, z = v:
   Then e_1 = 1 + u + v, e_2 = u + v + uv, e_3 = uv!
   Actually, with e_1 = 1 + s, e_2 = s + p, e_3 = p (where s = u+v, p = uv, z = 1):
   Notice s = e_2 - e_3, and 1 = e_1 - e_2 + e_3!
   Let's extract exact c_{a,b,c}^{(2m)} for 2m = 2, 4, 6, 8, 10, 12 and inspect patterns!
"""
import sympy as sp
from math import comb

def toeplitz_pn(n, x, y, z):
    d = z - x - y
    M = sp.Matrix(n, n, lambda i, j: d if i == j else (-comb(n, j - i) * x if j > i else -comb(n, n + j - i) * y))
    return sp.expand(M.det())

# 1. Verify Toeplitz formula against direct definition for n=2, 4, 6
a_sym, b_sym, x, y, z = sp.symbols('a b x y z')
for n in [2, 4, 6]:
    zeta = sp.exp(2 * sp.pi * sp.I / n)
    prod_expr = sp.Integer(1)
    for k in range(n):
        prod_expr *= (z - (a_sym + zeta**k * b_sym)**n)
    prod_exp = sp.expand(prod_expr)
    # Substitute a^n -> x, b^n -> y
    p_def = sp.expand(prod_exp.subs({a_sym**n: x, b_sym**n: y}))
    # Better: collect powers of a, b
    poly_ab = sp.Poly(prod_exp, a_sym, b_sym)
    p_def2 = sp.Integer(0)
    for (pa, pb), coeff in poly_ab.terms():
        c_simp = sp.simplify(coeff)
        if c_simp == 0:
            continue
        assert pa % n == 0 and pb % n == 0, f"Non-multiple power {pa}, {pb}"
        p_def2 += c_simp * x**(pa // n) * y**(pb // n)
    p_toep = toeplitz_pn(n, x, y, z)
    assert sp.expand(p_def2 - p_toep) == 0, f"Mismatch at n={n}!"
    print(f"n={n}: Toeplitz determinant EXACTLY matches definition!")

# 2. Extract exact c_{a,b,c}^{(n)} for n = 2, 4, 6, 8, 10!
# Fast method: set z = 1, x = u, y = v, convert symmetric poly in (u, v) to (s, p) = (u+v, uv),
# then substitute s = e2 - e3, and homogenize with e1 = 1 + e2 - e3 -> wait:
# P_n(1 + s, s + p, p) = Q_n(s, p)!
# Wait: since e1 = 1 + s, e2 = s + p, e3 = p, put s = e2 - e3, then e1 = 1 + e2 - e3!
# So to recover P_n(e1, e2, e3) from Q_n(s, p) = p_n(x, y, 1) (with s=x+y, p=xy):
# Notice P_n(e1, e2, e3) is homogeneous of weighted degree n in (e1:1, e2:2, e3:3),
# so P_n(e1, e2, e3) = (e1 - e2 + e3)... wait: at z = 1, z = e1 - s - ... NO!
# e1 = x+y+z, e2 = xy + z(x+y), e3 = xyz -> z satisfies z^3 - e1*z^2 + e2*z - e3 = 0.
# Putting z = 1 gives 1 - e1 + e2 - e3 = 0, i.e. e1 = 1 + e2 - e3!
# That gives P_n(1 + e2 - e3, e2, e3) = Q_n(e2 - e3, e3)!
# Wait: does P_n(1 + e2 - e3, e2, e3) UNIQUELY determine P_n(e1, e2, e3)?!
# YES!! Put e1 = 1 + t, e2 = t*E2... NO! Put e2 -> e2/e1^2, e3 -> e3/e1^3?
# Even simpler: P_n(1, E2, E3) = (1 / (1+s)^n) * Q_n(s, p) with E2 = (s+p)/(1+s)^2, E3 = p/(1+s)^3!
# No, sympy.symmetrize(p_n(x,y,z), [x,y,z], formal=True) does it directly in 1 line for n<=10!

e1, e2, e3 = sp.symbols('e1 e2 e3')
for n in [2, 4, 6, 8, 10]:
    p = toeplitz_pn(n, x, y, z)
    # Express in elementary symmetric polynomials via sympy!
    sym_expr, rem, mapping = sp.symmetrize(p, [x, y, z], formal=True)
    assert rem == 0
    # mapping maps formal symbols to e1, e2, e3
    subs_dict = {sym: target for (sym, expr), target in zip(mapping, [e1, e2, e3])}
    P = sp.Poly(sp.expand(sym_expr.subs(subs_dict)), e1, e2, e3)
    print(f"\n=== n = {n} (total terms = {len(P.terms())}, expected = {sum(1 for b in range(n//2+1) for c in range(n//3+1) if 2*b+3*c<=n)}) ===")
    sign_ok = True
    for c in range(n // 3 + 1):
        row = []
        for b in range((n - 3*c) // 2 + 1):
            a = n - 2*b - 3*c
            coeff = int(P.coeff_monomial(e1**a * e2**b * e3**c))
            exp_sign = (-1)**(b*(c+1) + c)
            if coeff == 0 or (1 if coeff > 0 else -1) != exp_sign:
                sign_ok = False
            row.append((b, coeff))
        print(f"  c={c}: " + ", ".join(f"b={b}:{val}" for b, val in row))
    print(f"  Sign law (-1)^(b(c+1)+c) & nonvanishing: {'PASS' if sign_ok else 'FAIL'}")
