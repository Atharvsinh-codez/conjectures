"""
Structural algebra probe of P_{2m}(e1, e2, e3) for Bardakov-Kozlovskaya Conjecture 3.3(2):
1. Check real-rootedness & interlacing of:
   - c-slices Q_c(1, u) = sum_b c_{2m-3c-2b, b, c} u^b (through n = 40!)
   - b-slices R_b(1, v) = sum_c c_{2m-2b-3c, b, c} v^c (Wait: v = e3 / e1^3!)
   - a-slices S_a(1, w) = sum_{2b+3c=2m-a} c_{a, b, c} w^c (even vs odd c!)
2. Check if P_{2m}(1, u, v), P_{2m}^even(1, -u, v), -P_{2m}^odd(1, u, v) are REAL-STABLE in (u, v)!
3. Check the EXACT algebraic sign pattern of p_{2m}(x, y, z) in OTHER symmetric bases:
   - monomial m_lambda(x, y, z)
   - Schur s_lambda(x, y, z)
   - power-sum p_lambda(x, y, z)
   - complete homogeneous h_lambda(x, y, z)
4. Check our EXACT Lucas resultant formula P_{2m}(1+s, s+p, p) = Res_u(L_{2m}(u) - s, u^{2m} - p)
   and inspect the exact c = 1, 2, 3 slices!
5. Check whether P_{2m}(e1, e2, e3) satisfies any clean recurrence/relation in m or differential equation in (e1, e2, e3)!
"""
import math
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 50

def poly_mul_xy(A, B):
    out = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if a:
            for j, b in enumerate(B):
                if b:
                    out[i + j] += a * b
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
        d = len(poly_xy) - 1
        rem = list(poly_xy)
        res = [0] * (d // 2 + 1)
        for b in range(d // 2 + 1):
            c = rem[b]
            res[b] = c
            if c:
                deg_s = d - 2 * b
                for r in range(deg_s + 1):
                    rem[b + r] -= c * math.comb(deg_s, r)
        return res

    F = [[0] * ((n - k) // 2 + 1) for k in range(n + 1)]
    for d in range(n + 1):
        k = n - d
        sgn = -1 if (d % 2 == 1) else 1
        sp_row = xy_to_sp(E[d])
        for b, v in enumerate(sp_row):
            F[k][b] = sgn * v

    coeffs = {}
    m = n
    c = 0
    while m >= 0:
        row0 = list(F[0])
        for b, val in enumerate(row0):
            coeffs[(m - 2 * b, b, c)] = val
        for b, val in enumerate(row0):
            if not val: continue
            a = m - 2 * b
            for j in range(b + 1):
                bp = b - j
                vj = val * math.comb(b, j)
                for i in range(a + 1):
                    F[i + j][bp] -= vj * math.comb(a, i)
        if m < 3:
            break
        m_new = m - 3
        F = [[F[k + 1][bp + 1] for bp in range((m_new - k) // 2 + 1)] for k in range(m_new + 1)]
        m = m_new
        c += 1
    return coeffs, E

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

# ============================================================
# 1. REAL-ROOTEDNESS & INTERLACING OF c-, b-, a-SLICES (n = 6..36)
# ============================================================
print("=== 1. REAL-ROOTEDNESS OF c-SLICES, b-SLICES, AND EVEN/ODD b-SLICES ===")
saved = {}
for n in range(2, 38, 2):
    coeffs, E = compute_p_n_fast(n)
    saved[n] = (coeffs, E)
    if n < 6: continue
    # c-slices Q_c(1, u)
    bad_c = []
    rts_c = {}
    for c in range(0, n // 3 + 1):
        deg_u = (n - 3 * c) // 2
        row = [coeffs[(n - 3*c - 2*b, b, c)] for b in range(deg_u + 1)]
        r = get_roots_scaled(row)
        rts_c[c] = r
        if r is None: bad_c.append(c)
    # Check interlacing of Q_c and Q_{c+2} (same sign family in u!)
    inter_c2 = 0; inter_c2_tot = 0
    for c in range(0, n // 3 - 1):
        r1, r2 = rts_c.get(c), rts_c.get(c + 2)
        if r1 is None or r2 is None or len(r2) == 0: continue
        # deg(Q_c) - deg(Q_{c+2}) = 3! Wait: too big to interlace directly!
        pass

    # b-slices R_b(1, v) = sum_c c_{n-2b-3c, b, c} v^c (and separately even c / odd c!)
    bad_b = []
    bad_b_even = []
    bad_b_odd = []
    for b in range(0, n // 2 + 1):
        max_c = (n - 2 * b) // 3
        row = [coeffs[(n - 2*b - 3*c, b, c)] for c in range(max_c + 1)]
        if len(row) >= 3 and get_roots_scaled(row) is None:
            bad_b.append(b)
        row_e = [coeffs[(n - 2*b - 3*c, b, c)] for c in range(0, max_c + 1, 2)]
        if len(row_e) >= 3 and get_roots_scaled(row_e) is None:
            bad_b_even.append(b)
        row_o = [coeffs[(n - 2*b - 3*c, b, c)] for c in range(1, max_c + 1, 2)]
        if len(row_o) >= 3 and get_roots_scaled(row_o) is None:
            bad_b_odd.append(b)
    if n <= 18 or n % 6 == 0 or bad_c or bad_b_even or bad_b_odd:
        print(f"  n={n:2d}: bad_c={bad_c}, bad_b(all)={bad_b}, bad_b(even_c)={bad_b_even}, bad_b(odd_c)={bad_b_odd}")

# ============================================================
# 2. TEST REAL-STABILITY OF P_even(1, -u, v) AND -P_odd(1, u, v) IN (u, v)
# ============================================================
print("\n=== 2. REAL-STABILITY TEST OF P_even(1, -u, w=v^2) AND -P_odd(1, u, w=v^2)/v ===")
rng = np.random.default_rng(42)
for n in [6, 8, 10, 12, 14, 16, 18]:
    coeffs, _ = saved[n]
    # P_even(1, -u, v) = sum_{c even} |c_{a,b,c}| u^b w^{c/2} where w = v^2!
    # Wait: test BOTH as poly in (u, w=v^2) AND in (u, v)!
    # Note: a polynomial with all positive coeffs in (u, w) is real-stable iff for all Im(u) > 0, roots in w have Im(w) <= 0!
    fail_even_uw = 0; fail_odd_uw = 0
    for _ in range(200):
        u_val = rng.uniform(-5, 5) + 1j * rng.uniform(0.05, 5)
        # evaluate coeff of w^k for k = 0..(n//6)
        max_ke = n // 6
        if max_ke >= 1:
            cw = [0j] * (max_ke + 1)
            for (a, b, c), val in coeffs.items():
                if c % 2 == 0:
                    cw[c // 2] += abs(val) * (u_val ** b)
            rts = np.roots(cw[::-1])
            if np.any(rts.imag > 1e-8):
                fail_even_uw += 1
        max_ko = (n - 3) // 6
        if max_ko >= 1:
            cw = [0j] * (max_ko + 1)
            for (a, b, c), val in coeffs.items():
                if c % 2 == 1:
                    cw[(c - 1) // 2] += abs(val) * (u_val ** b)
            rts = np.roots(cw[::-1])
            if np.any(rts.imag > 1e-8):
                fail_odd_uw += 1
    # Also test P_2m(1, u, v) directly for real-stability in (u, v) under all 4 sign flips (s1*u, s2*v)!
    rs_full = {}
    for s1 in (1, -1):
        for s2 in (1, -1):
            fl = 0
            for _ in range(100):
                u_val = rng.uniform(-3, 3) + 1j * rng.uniform(0.05, 3)
                cv = [0j] * (n // 3 + 1)
                for (a, b, c), val in coeffs.items():
                    cv[c] += val * ((s1 * u_val) ** b) * (s2 ** c)
                rts = np.roots(cv[::-1])
                if np.any(rts.imag > 1e-8):
                    fl += 1
            rs_full[(s1, s2)] = fl
    print(f"  n={n:2d}: P_even(1,-u,w) fails={fail_even_uw}/200, -P_odd(1,u,w)/v fails={fail_odd_uw}/200, P_full fails={rs_full}")

# ============================================================
# 3. CHECK SYMMETRIC BASES: m_lambda, s_lambda, p_lambda, h_lambda FOR n = 2, 4, 6, 8!
# ============================================================
print("\n=== 3. SYMMETRIC BASIS SIGN PATTERNS FOR n = 2, 4, 6, 8 ===")
x, y, z = sp.symbols('x y z')
for n in [2, 4, 6, 8]:
    coeffs, E = saved[n]
    # Build p_n(x, y, z) from E
    poly = 0
    for d in range(n + 1):
        sgn = -1 if d % 2 == 1 else 1
        poly += sgn * z**(n - d) * sum(E[d][j] * x**(d - j) * y**j for j in range(d + 1))
    poly = sp.Poly(sp.expand(poly), x, y, z)
    # m_lambda: read off x^p y^q z^r for p >= q >= r
    m_coeffs = {}
    for p in range(n, -1, -1):
        for q in range(min(p, n - p), -1, -1):
            r = n - p - q
            if 0 <= r <= q:
                m_coeffs[(p, q, r)] = int(poly.coeff_monomial(x**p * y**q * z**r))
    m_pos = sum(1 for v in m_coeffs.values() if v > 0)
    m_neg = sum(1 for v in m_coeffs.values() if v < 0)
    m_zer = sum(1 for v in m_coeffs.values() if v == 0)
    print(f"  n={n}: monomial m_lambda: +{m_pos}, -{m_neg}, 0:{m_zer}")
    if n <= 6:
        print(f"    m_lambda = {m_coeffs}")
    # Schur s_lambda: multiply p_n(x,y,z) by Vandermonde (x-y)(x-z)(y-z) and read coeff of x^{p+2} y^{q+1} z^r!
    V = (x - y) * (x - z) * (y - z)
    pV = sp.Poly(sp.expand(poly.as_expr() * V), x, y, z)
    s_coeffs = {}
    for (p, q, r) in m_coeffs:
        s_coeffs[(p, q, r)] = int(pV.coeff_monomial(x**(p+2) * y**(q+1) * z**r))
    s_pos = sum(1 for v in s_coeffs.values() if v > 0)
    s_neg = sum(1 for v in s_coeffs.values() if v < 0)
    s_zer = sum(1 for v in s_coeffs.values() if v == 0)
    print(f"  n={n}: Schur s_lambda:    +{s_pos}, -{s_neg}, 0:{s_zer}")
    if n <= 6:
        print(f"    s_lambda = {s_coeffs}")
    # Check signs of s_lambda vs (-1)^... or parity!
    if n == 8:
        print(f"    n=8 s_lambda = {s_coeffs}")
