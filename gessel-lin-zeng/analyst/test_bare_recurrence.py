"""CRITICAL TEST: Drop the outer D'_k completely!
Test the bare recurrence:
  F_{0,0} = 1,   F_{k,i} = D_{c_k, m_k} F_{k-1, i} + lam_k * F_{k-1, i-1}
where D_{c,m} = t(1+t)d/dt + c + m*t (0 < c < m) and lam_k > 0!
Equivalently: sum_i F_{k,i}(t) y^i = prod_{j=1}^k (D_{c_j, m_j} + lam_j * y)(1)!!

Check for WILD sequences (c_k, m_k, lam_k):
1. Are ALL F_{k,i}(t) real-rooted in (-1, 0)?
2. Does row interlacing F_{k,i} << F_{k,i-1} hold?
3. Does row log-concavity F_{k,i}^2 >= F_{k,i-1} F_{k,i+1} hold on (-1, 0)? On ALL R?!
4. ALSO test: WHAT IF lam_k IS HUGE OR TINY (e.g. lam = [1, 10^6, 1, 10^6] or [10^6, 1, 1, 1])?!
5. ALSO test: WHAT IF WE REPLACE D_{c,m} BY A *PURE MULTIPLICATION + DERIVATIVE* OPERATOR
   WITH DIFFERENT SUPPORT, e.g.,
   - Laguerre: t d/dt + c + m t  (roots in (-inf, 0))
   - Pure derivative on arbitrary polynomial: d/dt (Wait: d/dt lowers degree; instead (a+b*t) + c*d/dt!)
   - Even simpler: WHAT IF F_{0,0} IS AN *ARBITRARY* POLYNOMIAL WITH SIMPLE ROOTS IN (-1, 0)?!
"""
import sympy as sp

t = sp.Symbol('t')

def apply_D(c_val, m_val, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c_val)
        res[r + 1] += coef * (r + m_val)
    return res

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def poly_sub(a, b):
    return poly_add(a, [-v for v in b])

def poly_mul(a, b):
    res = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            res[i + j] += u * v
    return res

def exact_interlace(A, B):
    p0 = sp.Poly(list(reversed(A)), t, domain='QQ')
    p1 = sp.Poly(list(reversed(B)), t, domain='QQ')
    da, db = p0.degree(), p1.degree()
    if da == 0 or db == 0: return True
    if p0.LC() < 0: p0 = -p0
    if p1.LC() < 0: p1 = -p1
    chain = [p0, p1]
    while chain[-1].degree() > 0:
        _, r = sp.div(chain[-2], chain[-1], domain='QQ')
        if r.is_zero: return False
        r = -r
        r = r.exquo(sp.Poly(abs(r.LC()), t, domain='QQ'))
        chain.append(r)
    def sc(seq):
        nz = [u for u in seq if u != 0]
        return sum(1 for a, b in zip(nz, nz[1:]) if a * b < 0)
    vp = sc([1 if q.LC() > 0 else -1 for q in chain])
    vn = sc([(1 if q.LC() > 0 else -1) * ((-1)**q.degree()) for q in chain])
    return (vn - vp) == da

def test_bare(name, steps, init_poly=None):
    """steps = list of (c_k, m_k, lam_k)."""
    kmax = len(steps)
    F = [[None]*(k+1) for k in range(kmax+1)]
    F[0][0] = [sp.Rational(1)] if init_poly is None else init_poly
    rr_f = []; row_f = []; lc_m10 = []; lc_R = []
    for k in range(1, kmax + 1):
        c, m, lam = steps[k - 1]
        F[k][0] = apply_D(c, m, F[k-1][0])
        for i in range(1, k):
            F[k][i] = poly_add(apply_D(c, m, F[k-1][i]), [lam * v for v in F[k-1][i-1]])
        F[k][k] = [lam * v for v in F[k-1][k-1]]
        for i in range(k + 1):
            P = sp.Poly(list(reversed(F[k][i])), t, domain='QQ')
            if P.degree() > 0 and int(sp.count_roots(P, -1, 0)) < P.degree():
                rr_f.append((k, i))
        for i in range(1, k + 1):
            if not exact_interlace(F[k][i-1], F[k][i]):
                row_f.append((k, i))
        for i in range(1, k):
            D = poly_sub(poly_mul(F[k][i], F[k][i]), poly_mul(F[k][i-1], F[k][i+1]))
            Pd = sp.Poly(list(reversed(D)), t, domain='QQ')
            val = Pd.eval(sp.Rational(-1, 2))
            if val <= 0 or int(sp.count_roots(Pd, -1, 0)) > 0:
                lc_m10.append((k, i))
            if val <= 0 or int(sp.count_roots(Pd)) > 0:
                lc_R.append((k, i))
    print(f"  {name:46s}: RR={rr_f[:2]}, row={row_f[:2]}, lc(-1,0)={lc_m10[:2]}, lc(R)={lc_R[:2]}")

print("=== BARE RECURRENCE F_{k,i} = D_{c_k, m_k} F_{k-1,i} + lam_k F_{k-1,i-1} ===")
test_bare("c=1, m=k+1, lam=1 (k=1..7)", [(1, k+1, 1) for k in range(1, 8)])
test_bare("c=1, m=[2,10^5,2,10^5,2,10^5], lam=1", [(1, m, 1) for m in [2, 10**5, 2, 10**5, 2, 10**5]])
test_bare("c=1, m=[10^5,2,10^5,2,10^5,2], lam=1", [(1, m, 1) for m in [10**5, 2, 10**5, 2, 10**5, 2]])
test_bare("c=1, m=[2,3,4,5,6], lam=[1,10^6,1,10^6,1]", [(1, k+1, 1 if k%2 else 10**6) for k in range(1, 6)])
test_bare("c=1, m=[2,3,4,5,6], lam=[10^6,1,10^6,1,10^6]", [(1, k+1, 10**6 if k%2 else 1) for k in range(1, 6)])
test_bare("c=[1,100,1,100,1], m=[2,200,2,200,2], lam=1", [(1 if k%2 else 100, 2 if k%2 else 200, 1) for k in range(1, 6)])
test_bare("c=1, m=[2,1000,2,1000], lam=[1,10^6,1,10^6]", [(1, 2, 1), (1, 1000, 10**6), (1, 2, 1), (1, 1000, 10**6)])
test_bare("c=1, m=[1000,2,1000,2], lam=[1,10^6,1,10^6]", [(1, 1000, 1), (1, 2, 10**6), (1, 1000, 1), (1, 2, 10**6)])
test_bare("c=1, m=[2,1000,2,1000], lam=[10^6,1,10^6,1]", [(1, 2, 10**6), (1, 1000, 1), (1, 2, 10**6), (1, 1000, 1)])
# Also test with NONTRIVIAL initial polynomial F_{0,0} = (1+2t)(1+5t)!
test_bare("init=(1+2t)(1+5t), m=[2,10^4,2,10^4], lam=1", [(1, m, 1) for m in [2, 10**4, 2, 10**4]], init_poly=[1, 7, 10])
# Test EXTREME: what if m_k = 1.001 vs 10^6?
test_bare("m=[1001/1000, 10^6, 1001/1000, 10^6, 1001/1000]", [(1, sp.Rational(1001,1000) if k%2 else 10**6, 1) for k in range(1, 6)])
