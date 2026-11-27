"""CRITICAL EXPERIMENT:
1. Test if Gauss H_k(s, -z) = prod_{j=1}^k L_{m_j}(L_{m_j} - z)(1) is REAL-STABLE for arbitrary m_j > 0!
2. Test GLZ D_{m_k} with EXTREME sequences m_k (e.g. [1, 10^6, 1, 10^6], [10^6, 1, 10^6, 1], [10^6, 10^6, 1, 1],
   [1, 1, 10^6, 10^6], random m_k, fractional m_k > 0) for:
   - Real-rootedness of P_{k,i}
   - Row interlacing P_{k,i} << P_{k,i-1}
   - Pair (2,1) D_{m_k} P_{k-1,i} << P_{k-1,i-2}
   - Shifted equal-degree P_{k-1,i-1} < D_{m_k} P_{k-1,i}
   - Turan (D_{m_k} P_{k-1,i-1})^2 - P_{k-1,i-2} D_{m_k}^2 P_{k-1,i} > 0 on (-1, 0)!
3. ALSO test the SIMPLEST possible 1st-order operator: E_{m_k} = d/ds + m_k (constant coeff!) and
   F_{m_k} = s*d/ds + m_k*s (Laguerre) and general affine-drift operator!
"""
import numpy as np
import sympy as sp

s = sp.Symbol('s')
t = sp.Symbol('t')
z = sp.Symbol('z')

# 1. Test Gauss H_k(s, -z) real stability!
print("=== 1. TEST GAUSS H_k(s, -z) REAL STABILITY ===")
for name, ms in [("m=[1,2,3,4,5]", [1,2,3,4,5]), ("m=[1,100,1,100]", [1,100,1,100]), ("m=[100,1,100,1]", [100,1,100,1])]:
    H = sp.Integer(1)
    for m in ms:
        L = lambda f, mm=m: sp.diff(f, s) - mm * s * f
        H = sp.expand(L(L(H) - z * H))
    # Test UHP stability: fix z in H+, check all roots s have Im(s) <= 0!
    poly_s = sp.Poly(H, s)
    # Build coeff functions of z
    deg_s = poly_s.degree()
    coeffs_z = [sp.lambdify(z, poly_s.nth(d)) for d in range(deg_s, -1, -1)]
    rng = np.random.default_rng(7)
    max_im = -1e9; worst = None
    for _ in range(400):
        zv = complex(rng.normal(0, 5), 10 ** rng.uniform(-3, 2))
        c = [fn(zv) for fn in coeffs_z]
        r = np.roots(c)
        mi = float(np.max(r.imag))
        if mi > max_im: max_im = mi; worst = zv
    print(f"  Gauss {name}: max Im(s) for z in H+ = {max_im:+.3e} (at z={worst})")

# 2. Test EXTREME m_k sequences for GLZ D_m = t(1+t)d/dt + 1 + m*t!
print("\n=== 2. TEST EXTREME m_k SEQUENCES FOR GLZ D_m = t(1+t)d/dt + 1 + m*t ===")
def D_op(m, f):
    return t * (1 + t) * sp.diff(f, t) + (1 + m * t) * f

def exact_cauchy_interlace_t(p0, p1):
    p0 = sp.Poly(p0, t, domain='QQ')
    p1 = sp.Poly(p1, t, domain='QQ')
    da, db = p0.degree(), p1.degree()
    if da == 0 or db == 0: return True, "triv"
    if p0.LC() < 0: p0 = -p0
    if p1.LC() < 0: p1 = -p1
    chain = [p0, p1]
    while chain[-1].degree() > 0:
        _, r = sp.div(chain[-2], chain[-1], domain='QQ')
        if r.is_zero: return False, "common"
        r = -r
        r = r.exquo(sp.Poly(abs(r.LC()), t, domain='QQ'))
        chain.append(r)
    def sc(seq):
        nz = [u for u in seq if u != 0]
        return sum(1 for a, b in zip(nz, nz[1:]) if a * b < 0)
    vp = sc([1 if q.LC() > 0 else -1 for q in chain])
    vn = sc([(1 if q.LC() > 0 else -1) * ((-1)**q.degree()) for q in chain])
    ci = vn - vp
    if da == db:
        if ci == da: return True, "B<A"
        if ci == -da: return True, "A<B"
        return False, f"FAIL({ci}/{da})"
    if da == db + 1:
        if ci == da: return True, "A_outer"
        return False, f"FAIL({ci}/{da})"
    return False, f"deg({da},{db})"

extreme_seqs = [
    ("m=[2, 10000, 2, 10000, 2]", [2, 10000, 2, 10000, 2]),
    ("m=[10000, 2, 10000, 2, 10000]", [10000, 2, 10000, 2, 10000]),
    ("m=[10000, 10000, 2, 2, 2]", [10000, 10000, 2, 2, 2]),
    ("m=[2, 2, 2, 10000, 10000]", [2, 2, 2, 10000, 10000]),
    ("m=[1, 100, 10000, 1, 1]", [1, 100, 10000, 1, 1]),
    ("m=[10000, 100, 1, 1, 1]", [10000, 100, 1, 1, 1]),
    ("m=[1, 1, 1, 1, 1] (const)", [1, 1, 1, 1, 1]),
    ("m=[2, 2, 10^6, 2, 2, 2]", [2, 2, 10**6, 2, 2, 2]),
    ("m=[10^6, 10^6, 10^6, 2, 2, 2]", [10**6, 10**6, 10**6, 2, 2, 2]),
]
for name, ms in extreme_seqs:
    kmax = len(ms)
    P = [[None]*(k+1) for k in range(kmax+1)]
    P[0][0] = sp.Integer(1)
    for k in range(1, kmax + 1):
        m = ms[k - 1]
        P[k][0] = sp.expand(D_op(m, D_op(m, P[k-1][0])))
        for i in range(1, k):
            P[k][i] = sp.expand(D_op(m, D_op(m, P[k-1][i]) + P[k-1][i-1]))
        P[k][k] = sp.expand(D_op(m, P[k-1][k-1]))
    # Note: remove common factor (1+t) ONLY if m_0 == 1! Instead just count real roots of P[k][i]!
    rr_fails = []; row_fails = []; p21_fails = []; eq_fails = []
    for k in range(1, kmax + 1):
        for i in range(k + 1):
            poly = sp.Poly(P[k][i], t, domain='QQ')
            rc = int(sp.count_roots(poly))
            if rc < poly.degree():
                rr_fails.append((k, i, rc, poly.degree()))
        for i in range(1, k + 1):
            ok, d = exact_cauchy_interlace_t(P[k][i-1], P[k][i])
            if not ok: row_fails.append((k, i, d))
        if k >= 2:
            m = ms[k - 1]
            for i in range(1, k):
                ok, d = exact_cauchy_interlace_t(D_op(m, P[k-1][i]), P[k-1][i-1])
                if not ok: eq_fails.append((k, i, d))
            for i in range(2, k):
                ok, d = exact_cauchy_interlace_t(P[k-1][i-2], D_op(m, P[k-1][i]))
                if not ok: p21_fails.append((k, i, d))
    print(f"  {name:32s}: RR={rr_fails[:2]}, row={row_fails[:2]}, eq={eq_fails[:2]}, p21={p21_fails[:2]}")
