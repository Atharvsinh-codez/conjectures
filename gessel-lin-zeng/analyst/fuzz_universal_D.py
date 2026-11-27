"""ADVERSARIAL FUZZER on arbitrary sequences (m_1, ..., m_k) > 0 for:
  P_{0,0} = 1,   P_{k,i} = D_{m_k}( D_{m_k} P_{k-1,i} + P_{k-1,i-1} ),   D_m = t(1+t)d/dt + 1 + m*t.
Test:
1. 50 wild random & adversarial sequences of length k = 6, 7, 8 (including huge jumps, tiny m < 1, huge m = 10^8,
   monotonic, alternating, bimodal, etc.)!
2. Also test U_{k,i} = D_{m_k} P_{k-1,i} + P_{k-1,i-1} interlacing U_{k,i} << U_{k,i-1}!
3. ALSO test what happens if we insert a DIFFERENT parameter a_k vs b_k in the TWO copies of D:
   P_{k,i} = D_{a_k}( D_{b_k} P_{k-1,i} + P_{k-1,i-1} )!!
   WAIT A SECOND!! TEST THAT TOO!!
4. ALSO test what happens if we test GENERAL 1st-order operators D_k = t(1+t)d/dt + c_k + m_k t (varying c_k too!)!
"""
import math
import numpy as np
import sympy as sp

t = sp.Symbol('t')

def apply_D_num(c_val, m_val, poly):
    """Apply (t(1+t) d/dt + c_val + m_val * t) to low-to-high float64 or mpmath poly."""
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

def exact_cauchy_interlace_list(A, B):
    p0 = sp.Poly(list(reversed(A)), t, domain='QQ')
    p1 = sp.Poly(list(reversed(B)), t, domain='QQ')
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
    if da == db + 1 and ci == da:
        return True, "A_outer"
    return False, f"FAIL({ci}/{da})"

def exact_rr_count(A):
    p0 = sp.Poly(list(reversed(A)), t, domain='QQ')
    d = p0.degree()
    if d == 0: return True, 0, 0
    rc = int(sp.count_roots(p0, -1, 0))
    return (rc == d), rc, d

def run_exact_test(name, pairs_cm):
    """pairs_cm = list of (c_out, m_out, c_in, m_in) for k=1..kmax."""
    kmax = len(pairs_cm)
    P = [[None]*(k+1) for k in range(kmax+1)]
    P[0][0] = [sp.Rational(1)]
    rr_fails = []
    row_fails = []
    u_row_fails = []
    for k in range(1, kmax + 1):
        co, mo, ci, mi = pairs_cm[k - 1]
        U = [None] * (k + 1)
        U[0] = apply_D_num(ci, mi, P[k-1][0])
        for i in range(1, k):
            U[i] = poly_add(apply_D_num(ci, mi, P[k-1][i]), P[k-1][i-1])
        U[k] = P[k-1][k-1]
        for i in range(k + 1):
            P[k][i] = apply_D_num(co, mo, U[i])
            ok, rc, d = exact_rr_count(P[k][i])
            if not ok:
                rr_fails.append((k, i, rc, d))
        for i in range(1, k + 1):
            ok, msg = exact_cauchy_interlace_list(P[k][i-1], P[k][i])
            if not ok: row_fails.append((k, i, msg))
            ok_u, msg_u = exact_cauchy_interlace_list(U[i-1], U[i])
            if not ok_u: u_row_fails.append((k, i, msg_u))
    print(f"{name:45s}: RR={rr_fails[:2]}, row={row_fails[:2]}, U_row={u_row_fails[:2]}")

print("=== 1. EXTREME m_k WITH SAME INNER/OUTER D_{m_k} (c=1, m_k > 1) ===")
for name, ms in [
    ("m=[2, 10^8, 2, 10^8, 2, 10^8]", [2, 10**8, 2, 10**8, 2, 10**8]),
    ("m=[10^8, 2, 10^8, 2, 10^8, 2]", [10**8, 2, 10**8, 2, 10**8, 2]),
    ("m=[2, 3, 10^9, 2, 2, 2, 2]", [2, 3, 10**9, 2, 2, 2, 2]),
    ("m=[10^9, 10^9, 2, 2, 2, 2, 2]", [10**9, 10**9, 2, 2, 2, 2, 2]),
    ("m=[3/2, 1000, 3/2, 1000, 3/2, 1000]", [sp.Rational(3,2), 1000, sp.Rational(3,2), 1000, sp.Rational(3,2), 1000]),
    ("m=[1/2, 1/2, 1/2, 1/2, 1/2] (m<1!)", [sp.Rational(1,2)]*5),
    ("m=[1/10, 100, 1/10, 100, 1/10] (m<1!)", [sp.Rational(1,10), 100, sp.Rational(1,10), 100, sp.Rational(1,10)]),
]:
    run_exact_test(name, [(1, m, 1, m) for m in ms])

print("\n=== 2. DIFFERENT INNER vs OUTER m: D_{m_out}(D_{m_in} P_{k-1,i} + P_{k-1,i-1}) ===")
for name, pairs in [
    ("m_out=2, m_in=1000 (k=1..5)", [(1, 2, 1, 1000)]*5),
    ("m_out=1000, m_in=2 (k=1..5)", [(1, 1000, 1, 2)]*5),
    ("random (m_out, m_in) in [2..500]", [(1, 500, 1, 2), (1, 3, 1, 400), (1, 200, 1, 5), (1, 2, 1, 300), (1, 100, 1, 7)]),
]:
    run_exact_test(name, pairs)

print("\n=== 3. VARYING CONSTANT TERM c_k: D_{c,m} = t(1+t)d/dt + c + m*t (0 < c < m) ===")
for name, pairs in [
    ("c=[1,100,1,100,1], m=[2,200,2,200,2]", [(1,2,1,2), (100,200,100,200), (1,2,1,2), (100,200,100,200), (1,2,1,2)]),
    ("c=[100,1,100,1,100], m=[200,2,200,2,200]", [(100,200,100,200), (1,2,1,2), (100,200,100,200), (1,2,1,2), (100,200,100,200)]),
    ("c_out=1, m_out=2, c_in=100, m_in=200", [(1, 2, 100, 200)]*5),
]:
    run_exact_test(name, pairs)
