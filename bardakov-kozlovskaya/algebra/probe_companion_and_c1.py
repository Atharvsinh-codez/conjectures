"""
1. Inspect the exact m x m (and K x K for m = 3K) matrix M(e1) = Sigma_{2m}(C_Pi(e1))
   whose characteristic polynomial det(e2 I - M(e1)) IS P_{2m}(e1, e2, 1) (up to 2^{2m})!
2. Inspect exact gamma-coefficients of A_{2m-1}^{(2m)}(u) vs 4^b C(m-1, b) for c = 1 proof,
   AND check the exact c=1 alternating-sign margin at r_j = sin(j pi/2m)/sin((j+1) pi/2m) for m = 2..30!
3. Check the EXACT algebraic form of the e2 = 1 (and e3 = 1) Vandermonde-Cauchy integral
   to support @analyst's 2-chart nonnegative-square engine!
"""
import math
import sympy as sp
import mpmath as mp

mp.mp.dps = 50

print("=== 1. COMPANION MATRIX M(e1) = Sigma_{2m}(C_Pi(e1)) FOR m = 2, 3, 4, 5, 6 ===")
e1, e2, eta = sp.symbols('e1 e2 eta')
for m in [2, 3, 4, 5, 6]:
    if m % 3 == 0:
        K = m // 3
        # Pi_{6K}(U) = sum_{r=0}^K (6K/(3K-r)) C(3K-r, 2r) U^{K-r} eta^r, lc = 2
        # Monic: U^K + sum_{r=1}^K c_r U^{K-r} - e1/2
        cPi = [sp.Rational(m, m - r) * sp.binomial(m - r, 2*r) * eta**r for r in range(K + 1)]
        # cPi[0] = 1
        # Companion matrix C of U^K + c_1 U^{K-1} + ... + (c_K - e1/2)
        C = sp.zeros(K, K)
        for i in range(K - 1): C[i + 1, i] = 1
        for r in range(1, K + 1):
            val = cPi[r] - (e1 / 2 if r == K else 0)
            C[K - r, K - 1] = -val
        # Evaluate Sigma_{6K}(C) = sum_{j=0}^{2K} (-1)^j (m/(m-j)) C(2m-2j, j) eta^j C^{2K-j}
        M = sp.zeros(K, K)
        Cp = sp.eye(K)
        pows = [Cp]
        for _ in range(2 * K):
            Cp = sp.expand(Cp * C)
            pows.append(Cp)
        for j in range(2 * K + 1):
            coeff = (-1)**j * sp.Rational(m, m - j) * sp.binomial(2*m - 2*j, j) * eta**j
            M += coeff * pows[2*K - j]
        M = M.applyfunc(sp.expand)
        print(f"\n  m={m} (K={K}, {K}x{K} matrix M(e1) for R_{K}(e1,e2,eta) = (-2)^{{2K}} det(e2 I - M)):")
        for row in range(K):
            print("   ", [M[row, col] for col in range(K)])
    else:
        # m x m companion matrix of (Pi_{2m}(u) - e1)/2 = u^m + sum_{r=1}^{m//3} c_r u^{m-3r} - e1/2
        C = sp.zeros(m, m)
        for i in range(m - 1): C[i + 1, i] = 1
        for r in range(1, m // 3 + 1):
            c_r = sp.Rational(m, m - r) * sp.binomial(m - r, 2*r) * eta**r
            C[m - 3*r, m - 1] = -c_r
        C[0, m - 1] += e1 / 2
        M = sp.zeros(m, m)
        Cp = sp.eye(m)
        pows = [Cp]
        for _ in range(2 * m):
            Cp = sp.expand(Cp * C)
            pows.append(Cp)
        for j in range((2*m) // 3 + 1):
            coeff = (-1)**j * sp.Rational(m, m - j) * sp.binomial(2*m - 2*j, j) * eta**j
            M += coeff * pows[2*m - 3*j]
        M = M.applyfunc(sp.expand)
        print(f"\n  m={m} ({m}x{m} matrix M(e1) for P_{2*m}(e1,e2,eta^{m})):")
        for row in range(m):
            print("   ", [M[row, col] for col in range(m)])

print("\n=== 2. EXACT c=1 GAMMA-COEFFS vs 4^b C(m-1, b) AND EXACT SIGN-CHANGE MARGIN AT r_j ===")
for m in range(2, 13):
    n = 2 * m
    # Compute A_{n-1}^{(n)}(u)
    A = [sum(((-1)**r) * math.comb(n, r) * math.comb(n*(j - r + 1) - 1, n - 1) for r in range(j + 1)) for j in range(n)]
    # Divide by (1 + u) to get palindromic deg 2m-2 poly B(u)
    B = [0] * (n - 1)
    rem = list(A)
    for j in range(n - 1):
        B[j] = rem[j]
        rem[j + 1] -= rem[j]
    assert rem[-1] == 0
    # Convert B(u) to gamma-basis: B(u) = sum_{b=0}^{m-1} gamma_b (1+u)^{2m-2-2b} u^b
    rem = list(B)
    gam = [0] * m
    for b in range(m):
        gam[b] = rem[b]
        deg_s = 2*m - 2 - 2*b
        for r in range(deg_s + 1):
            rem[b + r] -= gam[b] * math.comb(deg_s, r)
    # Compare gam[b] with 4^b C(m-1, b) for EVEN b >= 2!
    even_ratios = [gam[b] / ((4**b) * math.comb(m - 1, b)) for b in range(2, m, 2)]
    min_er = min(even_ratios) if even_ratios else float('inf')
    # Also check B(u) coeff-wise? Wait: check min_{b>=1} gam[b] / (4^b C(m-1, b))
    all_ratios = [gam[b] / ((4**b) * math.comb(m - 1, b)) for b in range(1, m)]
    print(f"  m={m:2d}: min_even_b>=2 (gam_b / (4^b C(m-1,b))) = {min_er:.2f}, min_all_b>=1 = {min(all_ratios):.2f}")

print("\n=== 3. EXACT CHECK OF c=1 SIGN-CHANGE AT r_j (ACTUAL SUM, INCLUDING REINFORCEMENT!) FOR m=2..25 ===")
for m in [2, 3, 4, 5, 10, 15, 20, 25]:
    th0 = mp.pi / (2 * m)
    min_ratio_actual = mp.mpf('1e99')
    worst_j = -1
    for j in range(1, m):
        psi = j * mp.pi / (2 * m)
        s = mp.sin(psi) / mp.sin(psi + th0)
        t = s**(2 * m)
        # Actual real sum over k=0..m-1, times (-1)^j:
        k0 = (mp.sin(psi + th0) / mp.sin(th0))**(2 * m) / m
        actual_sum = sum(mp.re((1 - s * mp.exp(1j * (2*k+1)*mp.pi/(2*m)))**(-2*m)) for k in range(m)) / m
        sub = (1 - t) / ((1 + t)**2)
        signed_val = ((-1)**j) * (actual_sum - sub)
        ratio = signed_val / k0
        if ratio < min_ratio_actual:
            min_ratio_actual = ratio
            worst_j = j
    print(f"  m={m:2d}: min_j [(-1)^j (actual_sum - sub) / |k=0 term|] = {float(min_ratio_actual):.6f} (at j={worst_j})")
