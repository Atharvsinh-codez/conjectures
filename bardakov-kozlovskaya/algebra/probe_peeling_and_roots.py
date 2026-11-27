"""
Quick 4-second structural probe:
1. Check whether EVERY intermediate peeling layer G_c(z, s, p) in
     e3 * G_{c+1} = G_c - Q_c(s+z, p+s*z)
   obeys the exact sign law sign(F^{(c)}[k][b]) = (-1)^{b*(k+c+1) + (k+c)} for n = 4..24!
2. Track the roots v_j(u) of P_{2m}(1, u, v) = 0 for n = 6, 8, 10, 12 at u in (-2, 0.3)!
   Specifically: for u in (0, 1/4), are ALL roots v_j(u) REAL?!
3. Check our E_d dominant-product margin at arg(w_r(s)) = j*pi for n = 6, 8, 10, 12!
4. Check if P_{2m}(e1, e2, e3) has a clean positive expression in ANY other polynomial basis
   (e.g. h_lambda, (e1^2-4e2, e2, e3), (e1, e1^2-3e2, e3), cubic discriminant Delta_3)!
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

print("=== 1. SIGN LAW OF EVERY INTERMEDIATE PEELING LAYER G_c(z, s, p) (n = 4..24) ===")
saved_c = {}
for n in range(4, 26, 2):
    E = compute_E(n)
    F = [[(-1 if d%2 else 1)*v for v in xy_to_sp(E[d])] for d in range(n, -1, -1)]
    coeffs = {}; m = n; c = 0
    bad_layers = []
    while m >= 0:
        # Check EVERY entry F[k][b] in layer c!
        # Wait: in G_c(x,y,z), monomial z^k s^{m-k-2b} p^b came from z^{k+c} s^{m-k-2b} p^{b+c} in p_{2m}!
        # So effective (k_eff, b_eff) = (k+c, b+c)? Wait: at k=0, it's c_{m-2b, b, c} with sign (-1)^{b(c+1)+c}!
        # Check BOTH formulas: sign == (-1)^{b*(k+c+1) + (k+c)}!
        for k in range(m + 1):
            for b, val in enumerate(F[k]):
                want_pos = (((b * (k + c + 1) + (k + c)) % 2) == 0)
                if val == 0 or (val > 0) != want_pos:
                    bad_layers.append((c, k, b, val))
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
    saved_c[n] = coeffs
    print(f"  n={n:2d}: intermediate layer sign fails = {len(bad_layers)} {bad_layers[:3]}")

print("\n=== 2. ROOTS v_j(u) OF P_{2m}(1, u, v) = 0 FOR u IN (-0.5, 0.35) (n = 6, 8, 10, 12) ===")
for n in [6, 8, 10, 12]:
    coeffs = saved_c[n]
    deg_v = n // 3
    print(f"  n={n} (deg_v={deg_v}):")
    for u_val in [-0.5, -0.1, -0.01, 0.0, 0.05, 0.15, 0.24, 0.25, 0.30, 0.33]:
        cv = [0.0] * (deg_v + 1)
        for (a, b, c), val in coeffs.items():
            cv[c] += float(val) * (u_val ** b)
        # scale roots
        rts = np.roots(cv[::-1])
        n_real = sum(1 for r in rts if abs(r.imag) < 1e-7 * max(1.0, abs(r)))
        n_pos = sum(1 for r in rts if abs(r.imag) < 1e-7 * max(1.0, abs(r)) and r.real > 0)
        n_neg = sum(1 for r in rts if abs(r.imag) < 1e-7 * max(1.0, abs(r)) and r.real < 0)
        print(f"    u={u_val:6.2f}: real={n_real}/{deg_v} (pos={n_pos}, neg={n_neg})")

print("\n=== 3. CHECK COMPLETE HOMOGENEOUS h_lambda BASIS SIGNS (n = 2, 4, 6, 8) ===")
# Recall h_1 = e_1, h_2 = e_1^2 - e_2, h_3 = e_1^3 - 2 e_1 e_2 + e_3!
# Inverse: e_1 = h_1, e_2 = h_1^2 - h_2, e_3 = h_1^3 - 2 h_1 h_2 + h_3!
# Wait: for 3 variables, h_1, h_2, h_3 ALSO generate Lambda_3! Let's check P_{2m} in (h1, h2, h3)!
h1, h2, h3 = sp.symbols('h1 h2 h3')
e1, e2, e3 = sp.symbols('e1 e2 e3')
for n in [2, 4, 6, 8, 10]:
    coeffs = saved_c[n] if n in saved_c else {(2,0,0):1, (0,1,0):-4}
    P_e = sum(v * e1**a * e2**b * e3**c for (a, b, c), v in coeffs.items())
    P_h = sp.Poly(sp.expand(P_e.subs({e1: h1, e2: h1**2 - h2, e3: h1**3 - 2*h1*h2 + h3})), h1, h2, h3)
    pos = sum(1 for _, v in P_h.terms() if v > 0)
    neg = sum(1 for _, v in P_h.terms() if v < 0)
    print(f"  n={n:2d}: in (h1, h2, h3): terms={len(P_h.terms())}, +{pos}, -{neg}")
    if n <= 6:
        print(f"    P_h = {dict(P_h.terms())}")
