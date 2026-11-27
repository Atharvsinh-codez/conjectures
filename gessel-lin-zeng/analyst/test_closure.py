"""Test the exact closure conditions:
1. X < Y: p_{k-1,i-2} < Dtilde_k^2 p_{k-1,i}
2. Turan / log-concavity: M^2 - X*Y > 0 on (-1, 0), where
   X = p_{k-1,i-2}, M = Dtilde_k p_{k-1,i-1}, Y = Dtilde_k^2 p_{k-1,i}
   AND check its coefficients in x = t/(1+t) (i.e. A_{k,i}!) and in t!
3. Also check row log-concavity: p_{k,i}^2 - p_{k,i-1} p_{k,i+1} on (-1,0) and in x!
4. Also check if p_{k,i-1} < Dtilde_k p_{k,i} holds for ALL k up to 22!
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 50

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for idx, v in enumerate(a): res[idx] += v
    for idx, v in enumerate(b): res[idx] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def poly_sub(a, b):
    return poly_add(a, [-v for v in b])

def poly_mul(a, b):
    res = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            res[i + j] += u * v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_Dtilde(j, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(poly):
        res[m] += c * (m + 1)
        res[m + 1] += c * (m + j + 1)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def build_p(kmax):
    p = [[None] * (k + 1) for k in range(kmax + 1)]
    p[1][0] = [1, 2]
    p[1][1] = [1]
    for k in range(2, kmax + 1):
        p[k][0] = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][0]))
        for i in range(1, k):
            p[k][i] = apply_Dtilde(k, poly_add(apply_Dtilde(k, p[k - 1][i]), p[k - 1][i - 1]))
        p[k][k] = apply_Dtilde(k, p[k - 1][k - 1])
    return p

def get_roots(poly):
    if len(poly) <= 1: return []
    coeffs = [mp.mpf(c) for c in reversed(poly)]
    r = mp.polyroots(coeffs, maxsteps=200, extraprec=30)
    assert max(abs(mp.im(z)) for z in r) < mp.mpf('1e-25')
    return sorted(mp.re(z) for z in r)

def poly_eval(poly, tval):
    acc = mp.mpf(0)
    for c in reversed(poly):
        acc = acc * tval + c
    return acc

def to_x_coeffs(poly):
    """Convert p(t) of deg d to A(x) = (1-x)^d p(x/(1-x)) = sum c_m x^m (1-x)^{d-m}."""
    d = len(poly) - 1
    # Expand in x
    x = sp.Symbol('x')
    expr = sum(c * x**m * (1 - x)**(d - m) for m, c in enumerate(poly))
    P = sp.Poly(sp.expand(expr), x)
    # return low-to-high list
    coeffs = [0] * (d + 1)
    for (deg,), val in P.terms():
        coeffs[deg] = int(val)
    return coeffs

KMAX = 12
p = build_p(KMAX)

print("=== TEST 1: DOES X < Y (i.e. p_{k-1,i-2} < Dtilde_k^2 p_{k-1,i}) HOLD? ===")
xy_fails = []
for k in range(3, KMAX + 1):
    for i in range(2, k):
        X = p[k - 1][i - 2]
        Y = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][i]))
        rX = get_roots(X)
        rY = get_roots(Y)
        d = len(rX)
        ok = all(rX[m] < rY[m] for m in range(d)) and all(rY[m] < rX[m+1] for m in range(d-1))
        if not ok:
            xy_fails.append((k, i))
print(f"X < Y fails for k=3..{KMAX}: {xy_fails[:15]} (total {len(xy_fails)})")

print("\n=== TEST 2: DOES M^2 - X*Y > 0 ON (-1, 0) HOLD? AND COEFF SIGN IN t AND x? ===")
turan_neg_on_m1_0 = []
turan_t_neg_coeff = []
turan_x_neg_coeff = []
min_ratio_turan = (1e9, None)
for k in range(3, KMAX + 1):
    for i in range(2, k):
        X = p[k - 1][i - 2]
        M = apply_Dtilde(k, p[k - 1][i - 1])
        Y = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][i]))
        Dpoly = poly_sub(poly_mul(M, M), poly_mul(X, Y))
        # Check exact Sturm roots of Dpoly on (-1, 0)
        t = sp.Symbol('t')
        P = sp.Poly(list(reversed(Dpoly)), t, domain='ZZ')
        nroots = sp.count_roots(P, -1, 0)
        # Also check value at t = -0.5
        val_mid = poly_eval(Dpoly, mp.mpf('-0.5'))
        if nroots > 0 or val_mid <= 0:
            turan_neg_on_m1_0.append((k, i, nroots, float(val_mid)))
        if any(c < 0 for c in Dpoly):
            turan_t_neg_coeff.append((k, i))
        if k <= 9:
            xc = to_x_coeffs(Dpoly)
            if any(c < 0 for c in xc):
                turan_x_neg_coeff.append((k, i, [c for c in xc if c < 0][:3]))
print(f"M^2 - X*Y <= 0 somewhere on (-1,0) for k=3..{KMAX}: {turan_neg_on_m1_0[:10]}")
print(f"M^2 - X*Y has neg t-coeff: {len(turan_t_neg_coeff)} pairs, sample={turan_t_neg_coeff[:5]}")
print(f"M^2 - X*Y has neg x-coeff (k<=9): {len(turan_x_neg_coeff)} pairs, sample={turan_x_neg_coeff[:5]}")

print("\n=== TEST 3: WHAT ABOUT ROW LOG-CONCAVITY p_{k,i}^2 - p_{k,i-1}*p_{k,i+1} ON (-1,0) & ALL R? ===")
row_lc_fail = []
for k in range(2, KMAX + 1):
    for i in range(1, k):
        Dpoly = poly_sub(poly_mul(p[k][i], p[k][i]), poly_mul(p[k][i-1], p[k][i+1]))
        t = sp.Symbol('t')
        P = sp.Poly(list(reversed(Dpoly)), t, domain='ZZ')
        nroots_m1_0 = sp.count_roots(P, -1, 0)
        nroots_all = sp.count_roots(P)
        val_mid = poly_eval(Dpoly, mp.mpf('-0.5'))
        if nroots_m1_0 > 0 or val_mid <= 0:
            row_lc_fail.append((k, i, nroots_m1_0, nroots_all, float(val_mid)))
print(f"p_{{k,i}}^2 - p_{{k,i-1}}*p_{{k,i+1}} fails on (-1,0): {row_lc_fail[:10]}")

print("\n=== TEST 4: CHECK IF Dtilde_k p_{k-1,i} < Dtilde_k^2 p_{k-1,i} AND HOW MUCH MARGIN IN b_{r-1} < a_r ===")
# Specifically check min over r of (a_r - b_{r-1}) / (m_r - m_{r-1})
min_margin = (1e9, None)
for k in range(3, KMAX + 1):
    for i in range(2, k):
        X = p[k - 1][i - 2]
        M = apply_Dtilde(k, p[k - 1][i - 1])
        Y = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][i]))
        rM = get_roots(M)
        rXM = get_roots(poly_add(X, M))
        rMY = get_roots(poly_add(M, Y))
        for r in range(1, len(rM)):
            margin = float((rXM[r] - rMY[r - 1]) / (rM[r] - rM[r - 1]))
            if margin < min_margin[0]:
                min_margin = (margin, (k, i, r, len(rM)))
print(f"Min relative margin (a_r - b_{{r-1}})/(m_r - m_{{r-1}}) for k=3..{KMAX}: {min_margin}")
