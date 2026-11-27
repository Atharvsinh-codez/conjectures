"""
Verify and explore our BRAND-NEW Depressed-Cubic Resultant Identity for P_{2m}(e1, e2, e3):
  Let Pi_{2m}(u, eta)    = sum_{r=0}^{floor(m/3)}  (2m/(m-r)) C(m-r, 2r) u^{m-3r} eta^r
      Sigma_{2m}(u, eta) = sum_{j=0}^{floor(2m/3)} (-1)^j (m/(m-j)) C(2m-2j, j) u^{2m-3j} eta^j
  (where eta = v^2, e3 = v^{2m} = eta^m).

1. Verify in exact SymPy for m = 1, 2, 3, 4, 5, 6:
   - When 3 does not divide m: Res_u(Pi_{2m}(u, eta) - e1, Sigma_{2m}(u, eta) - e2) == P_{2m}(e1, e2, eta^m)!
   - When m = 3K: let R_K(e1, e2, eta) = Res_U(Pi_{6K}(U, eta) - e1, Sigma_{6K}(U, eta) - e2) (U = u^3);
     then prod_{w^3 = 1} R_K(e1, e2, w * eta) == P_{6K}(e1, e2, eta^{3K})!
2. Inspect R_K(e1, e2, eta) for K = 1, 2, 3 (n = 6, 12, 18) — what are its signs?!
3. Check the roots of Pi_{2m}(u, 1) = 0 and Sigma_{2m}(u, 1) = 0 (real-rootedness in u^3?!)!
4. Check our corrected peeling-layer sign law sign(F^{(c)}[k][b]) == (-1)^{(b+c)(k+c+1) + (k+c)} for n = 4..24!
"""
import math
import numpy as np
import sympy as sp

def poly_mul_xy(A, B):
    out = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if a:
            for j, b in enumerate(B):
                if b: out[i + j] += a * b
    return out

def compute_p_n_fast(n):
    P = [None]
    for r in range(1, n + 1):
        P.append([n * math.comb(n * r, n * j) for j in range(r + 1)])
    E = [[1]]
    for d in range(1, n + 1):
        acc = [0] * (d + 1)
        for r in range(1, d + 1):
            prod = poly_mul_xy(P[r], E[d - r])
            if (r - 1) % 2 == 0:
                for idx, v in enumerate(prod): acc[idx] += v
            else:
                for idx, v in enumerate(prod): acc[idx] -= v
        E.append([v // d for v in acc])
    def xy_to_sp(poly_xy):
        d = len(poly_xy) - 1; rem = list(poly_xy); res = [0] * (d // 2 + 1)
        for b in range(d // 2 + 1):
            c = rem[b]; res[b] = c
            if c:
                for r in range(d - 2*b + 1): rem[b + r] -= c * math.comb(d - 2*b, r)
        return res
    F = [[(-1 if d%2 else 1)*v for v in xy_to_sp(E[d])] for d in range(n, -1, -1)]
    coeffs = {}; m = n; c = 0
    bad_layers = 0
    while m >= 0:
        for k in range(m + 1):
            for b, val in enumerate(F[k]):
                want_pos = ((((b + c) * (k + c + 1) + (k + c)) % 2) == 0)
                if val == 0 or (val > 0) != want_pos:
                    bad_layers += 1
        row0 = list(F[0])
        for b, val in enumerate(row0): coeffs[(m - 2*b, b, c)] = val
        for b, val in enumerate(row0):
            if not val: continue
            a = m - 2*b
            for j in range(b + 1):
                bp = b - j; vj = val * math.comb(b, j)
                for i in range(a + 1): F[i + j][bp] -= vj * math.comb(a, i)
        if m < 3: break
        m_new = m - 3
        F = [[F[k+1][bp+1] for bp in range((m_new-k)//2+1)] for k in range(m_new+1)]
        m = m_new; c += 1
    return coeffs, bad_layers

print("=== 1. CORRECTED INTERMEDIATE PEELING LAYER SIGN LAW (n = 4..24) ===")
saved = {}
for n in range(4, 26, 2):
    coeffs, bad = compute_p_n_fast(n)
    saved[n] = coeffs
    if n <= 16 or n == 24 or bad > 0:
        print(f"  n={n:2d}: bad_layers = {bad}")

print("\n=== 2. EXACT SYMPY CHECK OF DEPRESSED-CUBIC RESULTANT (m = 1..6) ===")
u, U, e1, e2, eta = sp.symbols('u U e1 e2 eta')
for m in range(1, 7):
    n = 2 * m
    coeffs = saved[n] if n in saved else {(2,0,0): 1, (0,1,0): -4}
    P_target = sp.expand(sum(v * e1**a * e2**b * (eta**m)**c for (a, b, c), v in coeffs.items()))
    # Note: for j = m (only possible when m=0, skipped), 2m - 2j = 0. When j = 2m/3, m-j = m/3 > 0!
    Pi_u = sum(sp.Rational(2*m, m - r) * sp.binomial(m - r, 2*r) * u**(m - 3*r) * eta**r for r in range(m // 3 + 1))
    Sig_u = sum((-1)**j * sp.Rational(m, m - j) * sp.binomial(2*m - 2*j, j) * u**(2*m - 3*j) * eta**j for j in range((2*m) // 3 + 1))
    if m % 3 != 0:
        res = sp.expand(sp.resultant(Pi_u - e1, Sig_u - e2, u))
        ok = (sp.simplify(res - P_target) == 0)
        print(f"  m={m} (n={n}, 3 nmid m): Res_u(Pi - e1, Sig - e2) == P_{n}(e1, e2, eta^{m})? {ok}")
    else:
        K = m // 3
        Pi_U = sum(sp.Rational(2*m, m - r) * sp.binomial(m - r, 2*r) * U**(K - r) * eta**r for r in range(K + 1))
        Sig_U = sum((-1)**j * sp.Rational(m, m - j) * sp.binomial(2*m - 2*j, j) * U**(2*K - j) * eta**j for j in range(2*K + 1))
        R_K = sp.expand(sp.resultant(Pi_U - e1, Sig_U - e2, U))
        # Product over w^3 = 1: write R_K = A + B*eta + C*eta^2 (mod eta^3 -> use exact cube roots w = (-1+sqrt(-3))/2)
        w = sp.exp(2 * sp.pi * sp.I / 3)
        # Better: group R_K by eta^r mod 3: R_K = A0 + A1 + A2 where A_t has eta^{t mod 3}
        R_poly = sp.Poly(R_K, eta)
        A0 = sum(c_ * eta**d_ for (d_,), c_ in R_poly.terms() if d_ % 3 == 0)
        A1 = sum(c_ * eta**d_ for (d_,), c_ in R_poly.terms() if d_ % 3 == 1)
        A2 = sum(c_ * eta**d_ for (d_,), c_ in R_poly.terms() if d_ % 3 == 2)
        prod3 = sp.expand(A0**3 + A1**3 + A2**3 - 3 * A0 * A1 * A2)
        ok = (sp.simplify(prod3 - P_target) == 0)
        print(f"  m={m} (n={n}, K={K}): prod_{{w^3=1}} R_K(e1, e2, w*eta) == P_{n}(e1, e2, eta^{m})? {ok}")
        R_e = sp.Poly(R_K, e1, e2, eta)
        print(f"    R_{K}(e1, e2, eta) terms ({len(R_e.terms())}):")
        for (a, b, c_), val in sorted(R_e.terms(), key=lambda t: (t[0][2], t[0][1])):
            print(f"      e1^{a} e2^{b} eta^{c_}: {val}")

print("\n=== 3. ROOTS OF Pi_{2m}(u, 1) = 0 AND Sigma_{2m}(u, 1) = 0 IN U = u^3 (m = 3..15) ===")
for m in range(3, 16):
    # Pi_{2m}(u, 1) = u^{m%3} * sum_{r=0}^{m//3} coeff * U^{m//3 - r}
    K1 = m // 3
    cPi = [float(sp.Rational(2*m, m - r) * sp.binomial(m - r, 2*r)) for r in range(K1 + 1)]
    rPi = np.roots(cPi)
    imPi = np.max(np.abs(rPi.imag)) if len(rPi) else 0.0
    maxPi = np.max(rPi.real) if len(rPi) else 0.0
    K2 = (2*m) // 3
    cSig = [float((-1)**j * sp.Rational(m, m - j) * sp.binomial(2*m - 2*j, j)) for j in range(K2 + 1)]
    rSig = np.roots(cSig)
    imSig = np.max(np.abs(rSig.imag)) if len(rSig) else 0.0
    minSig = np.min(rSig.real) if len(rSig) else 0.0
    print(f"  m={m:2d}: Pi(U) deg={K1} max_im={imPi:.1e} max_real={maxPi:.4f} | Sig(U) deg={K2} max_im={imSig:.1e} min_real={minSig:.4f}")
