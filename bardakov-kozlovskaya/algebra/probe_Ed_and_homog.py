"""
Check:
1. Real-rootedness & gamma-signs of E_d(x, y) (the z^{2m-d} slice of p_{2m}(x, y, z)) for n = 2..30!
2. Signs of complete homogeneous h_lambda and forgotten f_lambda bases for n = 2..10!
3. Signs of -log P_{2m}(1, u, v) as a power series in (u, v) for n = 2..12!
4. Check if P_{2m}(e1, e2, e3) satisfies a clean identity in the Wendt matrix or Lucas resultant or differential operator!
5. Check the exact sign of E_d(1, u) at u = 1 (i.e. p_{2m}(1, 1, z)!) and u = -1 (i.e. p_{2m}(1, -1, z)!)!
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

def compute_E(n):
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
    return E

def xy_to_sp(poly_xy):
    d = len(poly_xy) - 1; rem = list(poly_xy); res = [0] * (d // 2 + 1)
    for b in range(d // 2 + 1):
        c = rem[b]; res[b] = c
        if c:
            for r in range(d - 2*b + 1): rem[b + r] -= c * math.comb(d - 2*b, r)
    return res

def get_roots_scaled(row):
    d = len(row) - 1
    if d <= 0: return np.array([])
    if row[0] == 0 or row[-1] == 0: return None
    s = (abs(row[0]) / abs(row[-1])) ** (1.0 / d)
    sc = [float(row[k] * (s ** k)) for k in range(d + 1)]
    mx = max(abs(x) for x in sc)
    sc = [x / mx for x in sc]
    rts = np.roots(sc[::-1])
    if np.max(np.abs(rts.imag)) > 1e-5:
        return None
    return np.sort(rts.real * s)

print("=== 1. CHECK REAL-ROOTEDNESS & GAMMA-SIGNS OF E_d(x, y) FOR n = 4..26 ===")
for n in range(4, 28, 2):
    E = compute_E(n)
    bad_rr = []
    bad_gamma = []
    for d in range(1, n):  # skip d=n which is (x-y)^n (multiple root at 1)
        r = get_roots_scaled(E[d])
        if r is None:
            bad_rr.append(d)
        sp_row = xy_to_sp(E[d])
        for b, val in enumerate(sp_row):
            want_pos = (True if d % 2 == 1 else (b % 2 == 0))
            if val == 0 or (val > 0) != want_pos:
                bad_gamma.append((d, b, val))
    print(f"  n={n:2d}: E_d non-real-rooted={bad_rr}, gamma-sign fails={len(bad_gamma)} {bad_gamma[:3]}")

print("\n=== 2. CHECK -log P_{2m}(1, u, v) COEFFICIENTS IN (u, v) FOR n = 6, 8, 10 ===")
u, v = sp.symbols('u v')
for n in [6, 8, 10]:
    E = compute_E(n)
    F = [[(-1 if d%2 else 1)*v_ for v_ in xy_to_sp(E[d])] for d in range(n, -1, -1)]
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
    # Compute -log(P(1, u, v)) truncated to u^6 v^4 using formal power series!
    # P(1, u, v) = 1 + R(u, v) where R(0, 0) = 0
    # Wait: compute [u^b v^c] of -log(1 + R) up to b <= 5, c <= 4!
    max_b, max_c = 5, 4
    R = [[sp.Integer(0)]*(max_c+1) for _ in range(max_b+1)]
    for (a, b, c), val in coeffs.items():
        if (b, c) != (0, 0) and b <= max_b and c <= max_c:
            R[b][c] = sp.Integer(val)
    def mul_trunc(A, B):
        C = [[sp.Integer(0)]*(max_c+1) for _ in range(max_b+1)]
        for i1 in range(max_b+1):
            for j1 in range(max_c+1):
                if A[i1][j1] == 0: continue
                for i2 in range(max_b+1 - i1):
                    for j2 in range(max_c+1 - j1):
                        if B[i2][j2] != 0:
                            C[i1+i2][j1+j2] += A[i1][j1] * B[i2][j2]
        return C
    # -log(1 + R) = sum_{k=1}^{max_b+max_c} (-1)^k / k * R^k
    L = [[sp.Integer(0)]*(max_c+1) for _ in range(max_b+1)]
    Rk = R
    for k in range(1, max_b + max_c + 1):
        sgn = sp.Rational((-1)**k, k)
        for i in range(max_b+1):
            for j in range(max_c+1):
                L[i][j] += sgn * Rk[i][j]
        if k < max_b + max_c:
            Rk = mul_trunc(Rk, R)
    neg = [(b, c, float(L[b][c])) for b in range(max_b+1) for c in range(max_c+1) if (b,c)!=(0,0) and L[b][c] <= 0]
    print(f"  n={n:2d}: -log P(1,u,v) non-positive coeffs (b<=5, c<=4): {len(neg)} {neg[:5]}")

print("\n=== 3. CHECK p_{2m}(1, -1, z) AND p_{2m}(1, 1, z) ROOTS & COEFFS (n = 2..12) ===")
for n in [2, 4, 6, 8, 10, 12]:
    E = compute_E(n)
    # p_{2m}(1, 1, z) = sum (-1)^d z^{n-d} E_d(1, 1)
    # p_{2m}(1, -1, z) = sum (-1)^d z^{n-d} E_d(1, -1)
    vals_p1 = [(-1 if d%2 else 1)*sum(E[d]) for d in range(n + 1)]
    vals_m1 = [(-1 if d%2 else 1)*sum(((-1)**j)*v for j, v in enumerate(E[d])) for d in range(n + 1)]
    print(f"  n={n:2d}: signs of p(1,1,z)={''.join('+' if v>0 else ('0' if v==0 else '-') for v in vals_p1)} | p(1,-1,z)={''.join('+' if v>0 else ('0' if v==0 else '-') for v in vals_m1)}")
