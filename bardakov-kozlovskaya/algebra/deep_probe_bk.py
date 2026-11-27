"""
Deep algebraic probe of Bardakov-Kozlovskaya Conjecture 3.3(2):
1. Push exact sign-law check from n = 26 up to n = 50 (pure Python int, ~2-3s!).
2. Check exact values at special slices:
   - e1 = 0 (a = 0): P_{2m}(0, e2, e3)
   - e2 = 0 (b = 0): P_{2m}(e1, 0, e3)
   - e1 = e2 = 0 (a = b = 0, 3 | 2m): c_{0,0,2m/3}^{(2m)} (WAIT: at n=6 it was 57^3 = 185193! Is it ALWAYS an m-th power?!)
   - c = 0, 1, 2, 3 slices Q_c(1, u): check real-rootedness in u (positive roots for even c, negative for odd c)!
3. Check prime factorizations & small-m patterns of c_{a,b,c}^{(2m)}!
4. Check symmetric-function expansions of p_{2m}(x, y, z) in m_lambda, s_lambda, h_lambda, p_lambda!
5. Check P_{2m}(e1, -e2, e3) vs P_{2m}(e1, e2, e3) and check odd n = 3, 5, 7, 9 for structural contrast!
"""
import math
import time
import numpy as np
import sympy as sp

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
    return coeffs

print("=== 1. EXACT SIGN LAW & NONVANISHING CHECK n = 2..50 ===")
t0 = time.time()
saved = {}
for n in range(2, 52, 2):
    coeffs = compute_p_n_fast(n)
    if n <= 24:
        saved[n] = coeffs
    zeros = [k for k, v in coeffs.items() if v == 0]
    sign_bad = [(k, v) for k, v in coeffs.items() if (v > 0) != (((k[1] * (k[2] + 1) + k[2]) % 2) == 0)]
    if n <= 16 or n % 6 == 0 or n >= 40 or zeros or sign_bad:
        print(f"  n={n:2d}: terms={len(coeffs):3d}, zeros={len(zeros)}, sign_bad={len(sign_bad)}, t={time.time()-t0:.2f}s")
    if zeros or sign_bad:
        print(f"    !!! ANOMALY at n={n}: zeros={zeros[:5]}, sign_bad={[(k, 1 if v>0 else -1) for k,v in sign_bad[:5]]}")

print("\n=== 2. REAL-ROOTEDNESS OF c-SLICES Q_c(1, u) FOR n = 6..24 ===")
for n in [6, 8, 10, 12, 14, 16, 18, 20, 22, 24]:
    coeffs = saved[n]
    max_c = n // 3
    non_rr = []
    for c in range(1, max_c + 1):
        # Q_c(1, u) = sum_b c_{n-3c-2b, b, c} u^b
        deg_u = (n - 3 * c) // 2
        if deg_u < 2: continue
        # rescale to float safely
        row = [coeffs[(n - 3*c - 2*b, b, c)] for b in range(deg_u + 1)]
        # normalize geometric scale
        s = (abs(row[0]) / abs(row[-1])) ** (1.0 / deg_u)
        sc = [float(row[b] * (s**b)) for b in range(deg_u + 1)]
        mx = max(abs(x) for x in sc)
        sc = [x / mx for x in sc]
        rts = np.roots(sc[::-1])
        max_im = np.max(np.abs(rts.imag))
        if max_im > 1e-6:
            non_rr.append((c, deg_u, float(max_im)))
    print(f"  n={n:2d} (max_c={max_c}): non-real-rooted c-slices = {non_rr}")

print("\n=== 3. FACTORIZATIONS OF SPECIAL COEFFICIENTS ===")
print("  (a) c_{0, 0, 2m/3}^{(2m)} for 6 | n (pure e3^{n/3} term):")
for n in [6, 12, 18, 24]:
    val = saved[n][(0, 0, n // 3)]
    # check if it is an (n/2)-th power? Wait: at n=6, 57^3 (m=3)! At n=12, check 6th root!
    m = n // 2
    root_est = round(val ** (1.0 / m))
    is_mth = (root_est ** m == val)
    # Wait: check 3rd root or (n/3)-th root!
    k3 = n // 3
    r3 = round(val ** (1.0 / k3))
    print(f"    n={n:2d}: val={val}, is_{m}th_pow={is_mth} ({root_est if is_mth else 'no'}), is_{k3}th_pow={r3**k3==val}")
    if n <= 18:
        print(f"      factors: {sp.factorint(val)}")

print("  (b) Leading c=1 coeff c_{n-3, 0, 1}^{(n)}:")
for n in range(4, 20, 2):
    v = -saved[n][(n - 3, 0, 1)]
    print(f"    n={n:2d}: -c_{{n-3,0,1}} = {v} = {sp.factorint(v)}")

print("  (c) Leading c=2 coeff c_{n-6, 0, 2}^{(n)}:")
for n in range(6, 20, 2):
    v = saved[n][(n - 6, 0, 2)]
    print(f"    n={n:2d}: c_{{n-6,0,2}} = {v} = {sp.factorint(v)}")

print("\n=== 4. ODD n = 3, 5, 7, 9 FOR STRUCTURAL CONTRAST ===")
for n in [3, 5, 7, 9]:
    coeffs = compute_p_n_fast(n)
    zeros = [k for k, v in coeffs.items() if v == 0]
    print(f"  odd n={n}: terms={len(coeffs)}, zeros={zeros}")
    for c in range(n // 3 + 1):
        row = sorted([(a, b, coeffs[(a, b, c)]) for (a, b, cc) in coeffs if cc == c], key=lambda t: t[1])
        print(f"    c={c}: " + ", ".join(f"({a},{b}):{v}" for a, b, v in row))
