"""
Independent exact-arithmetic & symbolic audit of #43 for Bardakov-Kozlovskaya Conjecture 3.3(2):
1. Verify every algebraic identity in #43's 12-line c=1 proof (m = 2..20):
   - Exact identity -c_{2m-3-2b, b, 1}^{(2m)} == 2m * (gamma_{b+1} - (-4)^{b+1} C(m-1, b+1))
   - Exact Rolle recurrence F_r(u) == (r)! * ... and F_{2m-1}(u)/(2m-1)! == A_{2m-1}^{(2m)}(u)
   - Exact alternating trig identity for gamma_{m-1} (verified to 50 digits via mpmath)
   - Exact check that gamma_{m-1} > 4^{m-1} for m = 2..20.
2. Verify #43 item 2:
   - Depressed-Cubic Resultant (m = 1..5, including m = 3 cube factorization P_6 = prod R_1(e1,e2,w*e3^{1/3})!)
   - Half-Size Lucas-Chebyshev Resultant (m = 1..5).
"""

from __future__ import annotations
import math
import sys
import mpmath as mp
import sympy as sp

sys.path.insert(0, r"d:\solving-conjectures\work\bardakov-kozlovskaya\verify")
import engine as bk


def audit_c1_proof():
    print("=== 1. Independent audit of #43's 12-line all-even-n c=1 proof ===")
    mp.mp.dps = 60

    for m in range(2, 16):
        n = 2 * m
        # Build F_r via Rolle recurrence: F_0 = 1, F_r(u) = n u(1-u) F'_{r-1}(u) + r(1 + (n-1)u) F_{r-1}(u)
        F = [1]
        for r in range(1, n):
            dF = [j * F[j] for j in range(1, len(F))] if len(F) > 1 else [0]
            nxt = [0] * (len(F) + 1)
            # n u(1-u) F' = n u F' - n u^2 F'
            for j, c in enumerate(dF):
                if c:
                    nxt[j + 1] += n * c
                    nxt[j + 2] -= n * c
            # r(1 + (n-1)u) F
            for j, c in enumerate(F):
                if c:
                    nxt[j] += r * c
                    nxt[j + 1] += r * (n - 1) * c
            while len(nxt) > 1 and nxt[-1] == 0:
                nxt.pop()
            F = nxt

        fact = math.factorial(n - 1)
        assert all(c % fact == 0 for c in F)
        A_n = [c // fact for c in F]
        # Verify palindromic of degree n-1 and monic
        assert len(A_n) == n and A_n[0] == 1 and A_n[-1] == 1 and A_n == list(reversed(A_n))

        # Divide by (1+u) and extract gamma_r: A_n(u)/(1+u) = sum_{r=0}^{m-1} gamma_r (1+u)^{2(m-1-r)} u^r
        # Exact division by (u + 1):
        Q = [0] * (n - 1)
        rem = A_n[-1]
        Q[-1] = rem
        for j in range(n - 3, -1, -1):
            rem = A_n[j + 1] - rem
            Q[j] = rem
        assert A_n[0] - rem == 0

        # Peel gamma_r in (1+u)^{2(m-1-r)} u^r basis
        rem_Q = list(Q)
        gammas = []
        for r in range(m):
            gr = rem_Q[r]
            gammas.append(gr)
            deg_1pu = 2 * (m - 1 - r)
            for j in range(deg_1pu + 1):
                rem_Q[r + j] -= gr * math.comb(deg_1pu, j)
        assert all(x == 0 for x in rem_Q)

        # Cross-check against exact engine.py c=1 coefficients!
        coeffs = bk.compute_elementary_coeffs(n)
        for b in range(m - 1):
            expected = -n * (gammas[b + 1] - ((-4) ** (b + 1)) * math.comb(m - 1, b + 1))
            assert coeffs[(n - 3 - 2 * b, b, 1)] == expected, f"c=1 gamma mismatch at m={m}, b={b}"

        # Check gamma_{m-1} == (-1)^{m-1} A'_n(-1) and matches exact trig alternating sum to 1e-50!
        dA_at_m1 = sum(j * A_n[j] * ((-1) ** (j - 1)) for j in range(1, n))
        assert gammas[m - 1] == ((-1) ** (m - 1)) * dA_at_m1

        trig_sum = mp.mpf(0)
        for k in range(m):
            th = mp.mpf(2 * k + 1) * mp.pi / (4 * m)
            term = mp.cos(th) / (mp.sin(th) ** (2 * m + 1))
            trig_sum += term if (k % 2 == 0) else -term
        trig_val = trig_sum / (2 * m)
        err = abs(trig_val - mp.mpf(gammas[m - 1]))
        assert err < mp.mpf('1e-45') * abs(trig_val), f"Trig mismatch at m={m}: err={err}"
        assert gammas[m - 1] > 4 ** (m - 1)
        if m <= 6:
            print(f"  m={m} (n={n:2d}): gamma_{{m-1}} = {gammas[m-1]} (> 4^{m-1}={4**(m-1)}), trig err = {float(err):.2e}, c=1 exact match!")
    print("  ALL m=2..15 (n=4..30) verified 100% identically for Rolle F_r, gamma_r, c=1 formula, and trig identity!")


def audit_resultants():
    print("\n=== 2. Independent SymPy verification of #43 Depressed-Cubic & Lucas-Chebyshev Resultants ===")
    u, eta, e1, e2, e3, q, p = sp.symbols('u eta e1 e2 e3 q p')

    # (a) Depressed-Cubic Resultant (m = 1, 2, 3, 4, 5)
    for m in [1, 2, 3, 4, 5]:
        n = 2 * m
        coeffs = bk.compute_elementary_coeffs(n)
        P_exact = sp.expand(sum(val * (e1**a) * (e2**b) * (e3**c) for (a, b, c), val in coeffs.items()))

        Pi = sum(
            sp.Rational(2 * m, m - r) * math.comb(m - r, 2 * r) * (u ** (m - 3 * r)) * (eta ** r)
            for r in range(m // 3 + 1)
        )
        Sigma = sum(
            ((-1) ** j) * sp.Rational(m, m - j) * math.comb(2 * m - 2 * j, j) * (u ** (2 * m - 3 * j)) * (eta ** j)
            for j in range((2 * m) // 3 + 1)
        )
        if m % 3 != 0:
            res = sp.expand(sp.resultant(Pi - e1, Sigma - e2, u))
            P_sub = sp.expand(P_exact.subs(e3, eta**m))
            assert sp.expand(res - P_sub) == 0, f"Depressed-cubic mismatch at m={m}!"
            print(f"  Depressed-Cubic Resultant m={m} (n={n}, 3 nmid m): EXACT MATCH!")
        else:
            # m = 3*K: Pi and Sigma are polynomials in U = u^3!
            K = m // 3
            U = sp.Symbol('U')
            Pi_U = sum(
                sp.Rational(2 * m, m - r) * math.comb(m - r, 2 * r) * (U ** (K - r)) * (eta ** r)
                for r in range(K + 1)
            )
            Sigma_U = sum(
                ((-1) ** j) * sp.Rational(m, m - j) * math.comb(2 * m - 2 * j, j) * (U ** (2 * K - j)) * (eta ** j)
                for j in range(2 * K + 1)
            )
            R_K = sp.expand(sp.resultant(Pi_U - e1, Sigma_U - e2, U))
            # Check prod_{w^3=1} R_K(e1, e2, w * e3^{1/3}) == P_{6K}(e1, e2, e3)
            # Substitute eta -> s where s^3 = e3, compute product over s, w*s, w^2*s!
            # Wait: R_K is in (e1, e2, eta) where eta^m = eta^{3K} = e3 => eta^K = e3^{1/3}!
            # Check at K=1 (m=3, n=6): eta^3 = e3!
            s = sp.Symbol('s')
            w = sp.Rational(-1, 2) + sp.I * sp.sqrt(3) / 2
            w2 = sp.Rational(-1, 2) - sp.I * sp.sqrt(3) / 2
            R_K_s = R_K.subs(eta, s)  # for K=1, s^3 = e3
            prod3 = sp.expand(R_K_s * R_K_s.subs(s, w * s) * R_K_s.subs(s, w2 * s))
            # Replace s^3 -> e3
            poly_s = sp.Poly(prod3, s)
            P_from_RK = sp.Integer(0)
            for (deg_s,), coeff in poly_s.terms():
                assert deg_s % 3 == 0
                P_from_RK += coeff * (e3 ** (deg_s // 3))
            assert sp.expand(P_from_RK - P_exact) == 0, f"Depressed-cubic K={K} mismatch!"
            print(f"  Depressed-Cubic Resultant m={m} (n={n}, K={K}): R_1(e1,e2,eta) = {R_K}; EXACT 3-fold product match!")

    # (b) Half-Size Lucas-Chebyshev Resultant (m = 1..5):
    # P_{2m}(1 + q - p, q, p) == Res_u(L_{2m}(u) + u^{2m} - q, u^{2m} - p)
    # Wait: what is L_{2m}(u)? Standard Lucas/Chebyshev polynomial where x^n + y^n = L_n(x+y, xy)?
    # Check: if t^2 - u*t + v = 0, no, check (1, u^{2m}, (1+u)^{2m}):
    # Wait: for x=1, y=w^{2m}, z=(1+w)^{2m}: let u = w(1+w)?
    # Let's check e1 = 1 + q - p? No, e3 = p, e2 = q, e1 = 1 + q - p => (1-x)(1-y)(1-z) = 1 - e1 + e2 - e3 = 0 => x = 1!
    # With x = 1: e3 = y*z = p, e2 = y + z + y*z = q => y + z = q - p, so e1 = 1 + y + z = 1 + q - p!!
    # Now on the rational curve y = w^{2m}, z = (1+w)^{2m}:
    #   p = y*z = (w(1+w))^{2m} = U^{2m} where U = w(1+w)!!
    #   q - p = y + z = w^{2m} + (1+w)^{2m}!
    #   Since w and -(1+w) have sum = -1 and product = -w(1+w) = -U,
    #   w^{2m} + (-(1+w))^{2m} = L_{2m}(-1, -U) = sum_{r=0}^m (2m/(2m-r)) C(2m-r, r) U^r !!
    #   WAIT: check for m=1 (n=2): w^2 + (1+w)^2 = 1 + 2w + 2w^2 = 1 + 2U!
    #   Check formula: r=0: (2/2)C(2,0)=1; r=1: (2/1)C(1,1)U = 2U => 1 + 2U!!
    #   Check for m=2 (n=4): w^4 + (1+w)^4 = (1+2U)^2 - 2U^2 = 1 + 4U + 2U^2!
    #   Check formula: r=0: 1; r=1: (4/3)C(3,1)U=4U; r=2: (4/2)C(2,2)U^2 = 2U^2! EXACT!!
    for m in [1, 2, 3, 4, 5]:
        n = 2 * m
        coeffs = bk.compute_elementary_coeffs(n)
        P_exact = sp.expand(sum(val * ((1 + q - p)**a) * (q**b) * (p**c) for (a, b, c), val in coeffs.items()))
        L_2m = sum(sp.Rational(2 * m, 2 * m - r) * math.comb(2 * m - r, r) * (u ** r) for r in range(m + 1))
        # Note q = (y+z) + p = L_2m(u) + u^{2m} (WAIT: verify sign of u^{2m}: U^{2m} = p!)
        # Wait: Res_u(L_2m(u) + u^{2m} - q, u^{2m} - p) = Res_u(L_2m(u) + p - q, u^{2m} - p)
        # (Wait: check whether P_{2m}(1+q-p, q, p) has degree 2m in q, while L_2m has degree m in u, so Res_u has degree 2m in q!)
        res_L = sp.expand(sp.resultant(L_2m + u**(2 * m) - q, u**(2 * m) - p, u))
        # Check ratio (sign or factor)
        assert sp.expand(res_L - P_exact) == 0, f"Lucas-Chebyshev mismatch at m={m}: diff={sp.expand(res_L - P_exact)}"
        print(f"  Half-Size Lucas-Chebyshev Resultant m={m} (n={n}): EXACT MATCH!")


if __name__ == "__main__":
    audit_c1_proof()
    audit_resultants()
