"""
Probe our new Lucas-Chebyshev resultant identity:
  P_{2m}(1 + q - p, q, p) = Res_u(L_{2m}(u) + u^{2m} - q, u^{2m} - p)
1. Verify exact identity for m = 1..5 in SymPy!
2. Verify our 10-line proofs of BOTH a=0 (e1=0) and b=0 (e2=0) slices for m = 1..20:
   - Count real roots of 1 + L_{2m}(u) = 0 in (-inf, -1) vs floor(m/3)
   - Count real roots of L_{2m}(u) + u^{2m} = 0 in (-1, -1/4) vs floor(2m/3)
3. Check real-rootedness of ALL a-slices S_a(e2, e3) in W = e3^2 / e2^3 (a=0,1,2,...) for n = 6..30!
4. Check our c=1 gamma-positivity proof: check min c_j = (1 - r_j)^2 / r_j for roots -r_j of A_{2m-1}^{(2m)}(u)!
5. Check complete homogeneous h_lambda and forgotten f_lambda basis signs for n = 2..10!
6. Check if P_{2m}(e1, e2, e3) has a clean recursion/relation for fixed c = 2, 3!
"""
import math
import numpy as np
import sympy as sp

def L_2m_poly(m, u):
    return sum((2 * m // (2 * m - k) if (2*m) % (2*m-k) == 0 else sp.Rational(2*m, 2*m-k)) * sp.binomial(2*m - k, k) * u**k for k in range(m + 1))

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
    while m >= 0:
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
    return coeffs

print("=== 1. EXACT SYMPY CHECK OF P_{2m}(1+q-p, q, p) = Res_u(L_{2m}(u) + u^{2m} - q, u^{2m} - p) ===")
u, q, p = sp.symbols('u q p')
for m in range(1, 6):
    n = 2 * m
    coeffs = compute_p_n_fast(n)
    lhs = sp.expand(sum(v * (1 + q - p)**a * q**b * p**c for (a, b, c), v in coeffs.items()))
    F2m = sp.expand(L_2m_poly(m, u) + u**(2*m) - q)
    rhs = sp.expand(sp.resultant(F2m, u**(2*m) - p, u))
    print(f"  m={m} (n={n}): diff == 0? {sp.simplify(lhs - rhs) == 0}")

print("\n=== 2. CHECK ROOT COUNTS FOR a=0 AND b=0 PROOFS (m = 1..15) ===")
for m in range(1, 16):
    # (a) 1 + L_{2m}(u) = 0 real roots in (-inf, -1)
    cL = [int(sp.Rational(2*m, 2*m - k) * sp.binomial(2*m - k, k)) for k in range(m + 1)]
    cL[0] += 1
    rL = np.roots(cL[::-1])
    real_L = [r.real for r in rL if abs(r.imag) < 1e-7 and r.real < -1 - 1e-7]
    # (b) L_{2m}(u) + u^{2m} = 0 real roots in (-1, -1/4)
    cF = [0] * (2*m + 1)
    for k in range(m + 1):
        cF[k] = int(sp.Rational(2*m, 2*m - k) * sp.binomial(2*m - k, k))
    cF[2*m] += 1
    rF = np.roots(cF[::-1])
    real_F = [r.real for r in rF if abs(r.imag) < 1e-7 and -1 + 1e-7 < r.real < -0.25 + 1e-7]
    ok_a = (len(real_L) == m // 3)
    ok_b = (len(real_F) == (2*m) // 3)
    if not (ok_a and ok_b) or m <= 10 or m == 15:
        print(f"  m={m:2d}: a=0 roots in (-inf,-1): {len(real_L)} == {m//3} ({ok_a}) | b=0 roots in (-1,-1/4): {len(real_F)} == {(2*m)//3} ({ok_b})")

print("\n=== 3. CHECK REAL-ROOTEDNESS OF ALL a-SLICES S_a(e2, e3) IN W = e3^2 / e2^3 (n = 6..30) ===")
def get_roots_scaled(row):
    d = len(row) - 1
    if d <= 0: return np.array([])
    s = (abs(row[0]) / abs(row[-1])) ** (1.0 / d)
    sc = [float(row[k] * (s ** k)) for k in range(d + 1)]
    mx = max(abs(x) for x in sc)
    sc = [x / mx for x in sc]
    rts = np.roots(sc[::-1])
    if np.max(np.abs(rts.imag)) > 1e-5:
        return None
    return np.sort(rts.real * s)

for n in range(6, 32, 2):
    coeffs = compute_p_n_fast(n)
    bad_a = []
    for a in range(0, n - 5):
        # For fixed a, c == a (mod 2), 2b + 3c = n - a -> c = c0 + 2k!
        c0 = a % 2
        row = []
        c = c0
        while 3 * c <= n - a:
            b = (n - a - 3 * c) // 2
            row.append(coeffs[(a, b, c)])
            c += 2
        if len(row) >= 3:
            r = get_roots_scaled(row)
            if r is None:
                bad_a.append(a)
    print(f"  n={n:2d}: non-real-rooted a-slices = {bad_a}")

print("\n=== 4. CHECK c=1 PROOF: MIN c_j = (1-r_j)^2 / r_j FOR ROOTS -r_j OF A_{2m-1}^{(2m)}(u) ===")
for m in range(2, 13):
    n = 2 * m
    # A_{n-1}^{(n)}(u) = sum_{j=0}^{n-1} [sum_{r=0}^j (-1)^r C(n, r) C(n*(j-r+1)-1, n-1)] u^j
    A = []
    for j in range(n):
        val = sum(((-1)**r) * math.comb(n, r) * math.comb(n*(j - r + 1) - 1, n - 1) for r in range(j + 1))
        A.append(val)
    # Roots of A(u)
    rts = np.roots(A[::-1])
    max_im = np.max(np.abs(rts.imag))
    # filter negative roots r < -1 (since roots come in pairs -r_j, -1/r_j plus -1)
    neg_r = sorted([-r.real for r in rts if r.real < -1.0001])
    cj = [(1 - r)**2 / r for r in neg_r]
    # Also compute exact gamma_b of A(u)/(1+u) vs 4^b C(m-1, b)
    print(f"  m={m:2d} (n={n:2d}): max_im={max_im:.1e}, count(r>1)={len(neg_r)}, min c_j={min(cj):.4f}, max c_j={max(cj):.2f}")
