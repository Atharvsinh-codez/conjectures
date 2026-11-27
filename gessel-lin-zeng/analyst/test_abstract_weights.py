"""GENIUS DIAGNOSTIC EXPERIMENT:
Does the recurrence
  phi_{0,0}(s) = w(s),   phi_{k,i}(s) = w(s) * (phi''_{k-1,i}(s) + phi'_{k-1,i-1}(s))
produce real-rooted interlacing functions for OTHER log-concave / PF_inf weights w(s)?
Test:
  1. Gaussian: w(s) = exp(-s^2/2)  => phi_{k,i}(s) = exp(-(k+1)s^2/2) * P_{k,i}(s) (exact polynomials in s!)
     WAIT: Let's check exact recurrence for P_{k,i}(s) when w(s) = exp(-s^2/2)!
     phi_{k,i} = e^{-(k+1)s^2/2} P_{k,i}(s).
     (e^{-m s^2/2} P)' = e^{-m s^2/2} (P' - m s P).
     So define L_m(P) = P' - m s P!
     Then P_{k,i} = L_k( L_k P_{k-1,i} + P_{k-1,i-1} )!!
     WHOA!! LOOK HOW IDENTICAL THAT IS TO h_{k,i} = D_k( D_k h_{k-1,i} + h_{k-1,i-1} )!!
  2. Exponential-of-exp (Gumbel): w(s) = e^s * exp(-e^s), phi_{k,i}(s) = e^s exp(-(k+1)e^s) Q_{k,i}(e^s)!
     Let u = e^s, d/ds = u d/du:
     (u e^{-m u} Q(u))' = u e^{-m u} [ u Q'(u) + (1 - m u) Q(u) ]!
     So M_m(Q) = u Q' + (1 - m u) Q (LAGUERRE-TYPE OPERATOR!!)
     and Q_{k,i} = M_k( M_k Q_{k-1,i} + Q_{k-1,i-1} )!!
  3. Rational/Cauchy (NOT log-concave!): w(s) = 1/(1+s^2)!
Let's test ALL THREE for real-rootedness, row interlacing, Pair (2,1), AND Turan M^2 - XY > 0!!
"""
import numpy as np
import sympy as sp

s = sp.Symbol('s')

def test_family(name, op_fn, kmax=10, domain_check="R"):
    """Build P_{k,i} via P_{0,0} = 1, P_{k,i} = op(k, op(k, P_{k-1,i}) + P_{k-1,i-1}) and test!"""
    P = [[None] * (k + 1) for k in range(kmax + 1)]
    P[0][0] = sp.Integer(1)
    for k in range(1, kmax + 1):
        P[k][0] = sp.expand(op_fn(k, op_fn(k, P[k-1][0])))
        for i in range(1, k):
            P[k][i] = sp.expand(op_fn(k, op_fn(k, P[k-1][i]) + P[k-1][i-1]))
        P[k][k] = sp.expand(op_fn(k, P[k-1][k-1]))

    print(f"\n=== FAMILY: {name} (k=1..{kmax}) ===")
    rr_fails = []
    for k in range(1, kmax + 1):
        for i in range(k + 1):
            poly = sp.Poly(P[k][i], s, domain='ZZ')
            d = poly.degree()
            if d == 0: continue
            rc = int(sp.count_roots(poly))
            if rc < d:
                rr_fails.append((k, i, rc, d))
    print(f"  Real-rootedness fails: {rr_fails[:10]} (total {len(rr_fails)})")

    # Check row outer-interlacing P_{k,i} << P_{k,i-1} and Pair(2,1) op(k, P_{k-1,i}) << P_{k-1,i-2} via numpy roots
    row_fails = []
    p21_fails = []
    for k in range(2, kmax + 1):
        for i in range(1, k + 1):
            r_big = sorted(np.roots([float(c) for c in sp.Poly(P[k][i-1], s).all_coeffs()]).real)
            r_sm = sorted(np.roots([float(c) for c in sp.Poly(P[k][i], s).all_coeffs()]).real)
            if len(r_sm) > 0 and len(r_big) == len(r_sm) + 1:
                ok = all(r_big[m] < r_sm[m] < r_big[m+1] for m in range(len(r_sm)))
                if not ok: row_fails.append((k, i))
        for i in range(2, k):
            A_p = sp.Poly(sp.expand(op_fn(k, P[k-1][i])), s)
            B_p = sp.Poly(P[k-1][i-2], s)
            r_sm = sorted(np.roots([float(c) for c in A_p.all_coeffs()]).real)
            r_big = sorted(np.roots([float(c) for c in B_p.all_coeffs()]).real)
            if len(r_sm) > 0 and len(r_big) == len(r_sm) + 1:
                ok = all(r_big[m] < r_sm[m] < r_big[m+1] for m in range(len(r_sm)))
                if not ok: p21_fails.append((k, i))
    print(f"  Row interlacing fails: {row_fails[:10]}")
    print(f"  Pair(2,1) fails:       {p21_fails[:10]}")

# 1. Gaussian: L_m(P) = P' - m * s * P
test_family("Gaussian w(s)=exp(-s^2/2), L_m = d/ds - m*s", lambda m, P: sp.diff(P, s) - m * s * P, kmax=8)

# 2. Gumbel / Laguerre: M_m(Q) = s*Q' + (1 - m*s)*Q (where s = u = e^x > 0)
test_family("Gumbel w=u*exp(-u), M_m = u*d/du + (1 - m*u)", lambda m, P: s * sp.diff(P, s) + (1 - m * s) * P, kmax=8)

# 3. Shifted Gaussian with CONSTANT m=1 (time-independent! Wait: what if w(s)=exp(-s^2/2) is NOT multiplied, i.e. L_1 = d/ds - s?):
# No, multiplying by w(s) each step increases m by 1! What if we multiply by w(s) = 1/(1+s^2) (NOT log-concave!)?
# ( (1+s^2)^{-m} P )' = (1+s^2)^{-m-1} [ (1+s^2) P' - 2 m s P ]
test_family("Cauchy w(s)=1/(1+s^2) [NOT log-concave], C_m = (1+s^2)d/ds - 2*m*s",
            lambda m, P: (1 + s**2) * sp.diff(P, s) - 2 * m * s * P, kmax=7)

# 4. Generalized Logistic: w(s) = e^{a*s} / (1 + e^s)^b  => M_m = s(1+s)d/ds + a - b*m*s?
# Wait: for GLZ, psi_{0,0} = e^s/(1+e^s), each step multiplies by 1/(1+e^s), so weight at step k is e^s/(1+e^s)^{k+1}.
# What about w(s) = 1/cosh(s) = 2 e^s/(1+e^{2s})? That's same as logistic!
# What about Gamma/Laguerre with alpha=0: M_m(Q) = s*Q' - m*s*Q? No, need constant term for zeros > 0!
# What about Hermite-shifted: L_m = d/ds - (m + 0.5)*s? Let's test m=2*k (faster weight growth!)
test_family("Fast Gaussian w(s)=exp(-s^2), L_m = d/ds - 2*m*s", lambda m, P: sp.diff(P, s) - 2 * m * s * P, kmax=7)
