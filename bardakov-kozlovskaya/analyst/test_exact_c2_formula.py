"""
Explore the EXACT c = 2 formula using our exponential H_{2m, r}(u) identity:
  p_{2m}(1, u, z) = (1 - u)^{2m} * exp( -2m * sum_{r>=1} (H_{2m,r}(u) / r) * (z / (1-u)^{2m})^r )
where H_{n, r}(u) = (1-u)^{nr} * sum_{j>=0} comb(n(j+r)-1, nr-1) u^j  (palindromic, deg nr - r).

Expand p_{2m}(1, u, z) to O(z^2):
  [z^0]: (1 - u)^{2m}
  [z^1]: -2m * H_{2m, 1}(u)
  [z^2]: (1 - u)^{-2m} * [ (2m)^2/2 * H_{2m, 1}(u)^2 - m * H_{2m, 2}(u) ]
       = (m / (1 - u)^{2m}) * [ 2m * H_{2m, 1}(u)^2 - H_{2m, 2}(u) ] !!
WAIT A SECOND!! CHECK THAT FORMULA FOR [z^2] p_{2m}(1, u, z) AND EXTRACT Q_{2m, 2}(1+u, u)!!
Remember:
  p_{2m}(1, u, z) = Q_0(e1, e2) + e3 * Q_1(e1, e2) + e3^2 * Q_2(e1, e2) + O(z^3)
  with e1 = 1 + u + z, e2 = u + (1+u)z, e3 = u*z!
So [z^2] p_{2m}(1, u, z) = [z^2] Q_0(1+u+z, u+(1+u)z) + u * [z^1] Q_1(1+u+z, u+(1+u)z) + u^2 * Q_2(1+u, u)!!
Compute BOTH subtraction terms explicitly in SymPy for n = 4, 6, 8, 10 and inspect the structure of Q_2(1+u, u)!!
"""
import sympy as sp
from math import comb
import sys
sys.path.insert(0, r'd:\solving-conjectures\work\bardakov-kozlovskaya\verify')
import engine

u, z, e1, e2 = sp.symbols('u z e1 e2')

def H_nr(n, r, u_sym):
    # deg = n*r - r
    deg = n * r - r
    coeffs = [sum((-1)**j * comb(n*r, j) * comb(n*(k - j + r) - 1, n*r - 1) for j in range(k + 1)) for k in range(deg + 1)]
    return sum(coeffs[k] * u_sym**k for k in range(deg + 1))

for n in [4, 6, 8, 10]:
    m = n // 2
    H1 = H_nr(n, 1, u)
    H2 = H_nr(n, 2, u)
    # Check [z^2] p_n(1, u, z) = m * (2m * H1^2 - H2) / (1 - u)^{2m}
    num = sp.expand(m * (2*m * H1**2 - H2))
    q, rem = sp.div(num, (1 - u)**n, u)
    assert rem == 0, f"Divisibility failed at n={n}!"
    # Compare with engine E_table: [z^2] p_n(1, u, z) = (-1)^{n-2} E_{n-2}(1, u) = E_{n-2}(1, u)
    E = engine.compute_E_table(n)
    E_z2 = sum(E[n - 2][j] * u**j for j in range(n - 1))
    assert sp.expand(q - E_z2) == 0, f"Mismatch [z^2] at n={n}!"
    # Now check Q_2(1+u, u)!
    # [z^2] Q_0(1+u+z, u+(1+u)z) where Q_0 = (e1^2 - 4*e2)^m:
    # (1+u+z)^2 - 4*(u+(1+u)z) = (1-u)^2 - 2*(1+u)*z + z^2!
    # So [z^2] ((1-u)^2 - 2(1+u)z + z^2)^m = m*(1-u)^{2m-2} + comb(m, 2)*4*(1+u)^2*(1-u)^{2m-4}
    Q0_z2 = sp.expand(m*(1-u)**(n-2) + 2*m*(m-1)*(1+u)**2*(1-u)**(n-4))
    # What is [z^1] Q_1(1+u+z, u+(1+u)z)?
    # Notice (d/de1 + (1+u)*d/de2) Q_1(e1, e2) at (1+u, u)!
    # For Q_1(e1, e2) = sum c_b e1^{d-2b} e2^b, f_1(u) = Q_1(1+u, u):
    # Wait: d/du Q_1(1+u, u) = (d/de1 + d/de2) Q_1! Not (d/de1 + (1+u)d/de2)!
    # Wait: Euler homogeneity: (e1 d/de1 + 2 e2 d/de2) Q_1 = (n-3) Q_1!
    # At e1=1+u, e2=u: (1+u) Q1_1 + 2u Q1_2 = (n-3) f_1(u)
    # and              Q1_1 + Q1_2 = f_1'(u)!
    # Solve this 2x2 linear system for Q1_1 + (1+u) Q1_2!!
    # Wait: det = (1+u) - 2u = 1 - u!!
    # Q1_1 = ((n-3) f_1 - 2u f_1') / (1 - u)
    # Q1_2 = ((1+u) f_1' - (n-3) f_1) / (1 - u)
    # => Q1_1 + (1+u) Q1_2 = [ (n-3)(1 - (1+u)) f_1 + ((1+u)^2 - 2u) f_1' ] / (1 - u)
    #                      = [ -(n-3) u f_1(u) + (1 + u^2) f_1'(u) ] / (1 - u) !!
    # WAIT A SECOND!! CHECK THAT EXACT IDENTITY FOR Q_2(1+u, u)!!
    f1 = sp.cancel(-n * (H1 - (1+u)*(1-u)**(n-2)) / u)
    Q1_z1 = sp.cancel((-(n - 3)*u*f1 + (1 + u**2)*sp.diff(f1, u)) / (1 - u))
    f2_pred = sp.cancel((q - Q0_z2 - u * Q1_z1) / u**2)
    C = engine.compute_elementary_coeffs(n)
    f2_exact = sum(C[(n - 6 - 2*b, b, 2)] * u**b * (1 + u)**(n - 6 - 2*b) for b in range(max(0, (n-6)//2 + 1))) if n >= 6 else sp.Integer(0)
    assert sp.expand(f2_pred - f2_exact) == 0, f"Mismatch f2 at n={n}!"
    print(f"n={n:2d}: EXACT closed form for f_{{n,2}}(u) in terms of H1, H2 VERIFIED! deg={n-6}")
    # Also check roots of E_{n-2}(1, u) = [z^2] p_n(1, u, z)!
    import numpy as np
    if n >= 6:
        r_E = np.roots([float(E[n-2][j]) for j in range(n-1)][::-1])
        r_H2 = np.roots([float(sp.Poly(H2, u).nth(j)) for j in range(2*n-2+1)][::-1])
        print(f"   E_{{n-2}}([z^2]) max_im={np.max(np.abs(r_E.imag)):.2e} (Re:[{r_E.real.min():.3f},{r_E.real.max():.3f}]) | H2 max_im={np.max(np.abs(r_H2.imag)):.2e}")
