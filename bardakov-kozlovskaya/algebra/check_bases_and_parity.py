"""
Check ALL symmetric bases (e_lambda, m_lambda, s_lambda, h_lambda, p_lambda) for n = 2..30!
Does the EXACT parity sign law:
    sign = +1 if (k1, k2, k3) all EVEN,
    sign = -1 if TWO of (k1, k2, k3) are ODD
hold for m_lambda, s_lambda, e_lambda, h_lambda, p_lambda?!
Also check:
- Does P_{2m}(e1, e2, e3) in (D = e1^2 - 4 e2, e1, e3) have a clean sign law?
- What is the sign of c_{a,b,c}^{(2m)} modulo powers of 2 / prime factors?
- Check the exact c = 2 formula from Newton/Veronese power sums!
"""
import math
import sympy as sp

def poly_mul_xy(A, B):
    out = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if a:
            for j, b in enumerate(B):
                if b:
                    out[i + j] += a * b
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

print("=== 1. CHECK MONOMIAL m_lambda & SCHUR s_lambda PARITY SIGN LAW FOR n = 2..30 ===")
for n in range(2, 32, 2):
    E = compute_E(n)
    # Monomial coeff M[(p, q, r)] = [x^p y^q z^r] p_n(x, y, z) = (-1)^{n-r} E[n-r][q] (since z^r means d = n-r, x^{d-q} y^q)
    m_bad = []
    m_zero = []
    for p in range(n, -1, -1):
        for q in range(min(p, n - p), -1, -1):
            r = n - p - q
            if 0 <= r <= q:
                d = n - r
                val = (-1 if d % 2 else 1) * E[d][q]
                want_pos = (p % 2 == 0 and q % 2 == 0 and r % 2 == 0)
                if val == 0:
                    m_zero.append((p, q, r))
                elif (val > 0) != want_pos:
                    m_bad.append(((p, q, r), 1 if val > 0 else -1))
    # Schur coeff S[(p, q, r)] = [x^{p+2} y^{q+1} z^r] p_n(x, y, z) * (x^2 y - x^2 z - x y^2 + x z^2 + y^2 z - y z^2)
    def get_M(a, b, c):
        if a < 0 or b < 0 or c < 0 or a + b + c != n: return 0
        d = n - c
        return (-1 if d % 2 else 1) * E[d][b]
    s_bad = []
    s_zero = []
    for p in range(n, -1, -1):
        for q in range(min(p, n - p), -1, -1):
            r = n - p - q
            if 0 <= r <= q:
                # Vandermonde V = x^2 y - x^2 z - x y^2 + x z^2 + y^2 z - y z^2
                # coeff of x^{p+2} y^{q+1} z^r in p_n * V:
                val = (get_M(p, q, r)
                       - get_M(p, q+1, r-1)
                       - get_M(p+1, q-1, r)
                       + get_M(p+1, q+1, r-2)
                       + get_M(p+2, q-1, r-1)
                       - get_M(p+2, q, r-2))
                want_pos = (p % 2 == 0 and q % 2 == 0 and r % 2 == 0)
                if val == 0:
                    s_zero.append((p, q, r))
                elif (val > 0) != want_pos:
                    s_bad.append(((p, q, r), 1 if val > 0 else -1))
    print(f"  n={n:2d}: m_zero={len(m_zero)}, m_bad={len(m_bad)} {m_bad[:3]} | s_zero={len(s_zero)}, s_bad={len(s_bad)} {s_bad[:3]}")

print("\n=== 2. CHECK P_{2m} IN (D = e1^2 - 4 e2, e1, e3) AND IN POWER SUMS p_lambda (n = 2..12) ===")
e1, e2, e3, D = sp.symbols('e1 e2 e3 D')
for n in [2, 4, 6, 8, 10, 12]:
    E = compute_E(n)
    # Quick exact e-coeffs
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
    P_expr = sum(v * e1**a * e2**b * e3**c for (a, b, c), v in coeffs.items())
    # Substitute e2 = (e1^2 - D)/4
    P_D = sp.Poly(sp.expand(P_expr.subs(e2, (e1**2 - D) / 4)), e1, D, e3)
    # Check signs of (a', b', c') in e1^{a'} D^{b'} e3^{c'}!
    # Notice: at c=0, only D^{n/2}!
    d_signs = []
    for (aa, bb, cc), val in sorted(P_D.terms(), key=lambda t: (t[0][2], t[0][1])):
        d_signs.append(((aa, bb, cc), '+' if val > 0 else '-'))
    # Check if D-coeffs have a clean sign pattern:
    d_bad = [((aa,bb,cc), val) for (aa,bb,cc), val in P_D.terms() if cc > 0 and ((val > 0) != ((bb + cc) % 2 == 0))]
    print(f"  n={n:2d}: in (e1, D=e1^2-4e2, e3): terms={len(P_D.terms())}, (sign==(-1)^(b'+c') for c>0) fails={len(d_bad)} {[(k, '+' if v>0 else '-') for k,v in d_bad[:4]]}")
    if n <= 8:
        print(f"    P_D = {dict(P_D.terms())}")
    # Also check power-sum basis (p1, p2, p3): e1 = p1, e2 = (p1^2 - p2)/2, e3 = (p1^3 - 3 p1 p2 + 2 p3)/6!
    p1, p2, p3 = sp.symbols('p1 p2 p3')
    P_p = sp.Poly(sp.expand(P_expr.subs({e1: p1, e2: (p1**2 - p2)/2, e3: (p1**3 - 3*p1*p2 + 2*p3)/6})), p1, p2, p3)
    p_pos = sum(1 for _, v in P_p.terms() if v > 0)
    p_neg = sum(1 for _, v in P_p.terms() if v < 0)
    p_bad = [((aa,bb,cc), v) for (aa,bb,cc), v in P_p.terms() if (v > 0) != (((aa + cc) % 2 == 0) and (bb % 2 == 0))]
    print(f"  n={n:2d}: in (p1, p2, p3): terms={len(P_p.terms())}, +{p_pos}, -{p_neg}")
    if n <= 6:
        print(f"    P_p = {dict(P_p.terms())}")
