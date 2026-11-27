"""BRUTAL ADVERSARIAL ATTACK ON THE ATOMIC THEOREM:
Try to break real-rootedness of [y^i] prod_{j=1}^N (D_{c_j, m_j} + lam_j * y)(1)
using:
1. Bimodal lambda_j = [1]*a + [1e-12]*b + [1]*a with HUGE m_j jump in the middle!
   (This isolates early-only + late-only + split subsets, mimicking the (13,7) barrier!)
2. GLZ with double-D: prod_{j=1}^k D_j (D_j + lam_j * y)(1) for k = 14 with bimodal lam_j!
3. Laguerre with bimodal lambda_j and huge m_j jump!
"""
import numpy as np
import sympy as sp

t = sp.Symbol('t')

def apply_D_f(c, m, poly):
    d = len(poly) - 1
    res = np.zeros(d + 2, dtype=np.float64)
    for r in range(d + 1):
        res[r] += poly[r] * (r + c)
        res[r + 1] += poly[r] * (r + m)
    return res

def apply_Lag_f(c, m, poly):
    d = len(poly) - 1
    res = np.zeros(d + 2, dtype=np.float64)
    for r in range(d + 1):
        res[r] += poly[r] * (r + c)
        res[r + 1] += poly[r] * m
    return res

def pad_add(a, b):
    n = max(len(a), len(b))
    out = np.zeros(n, dtype=np.float64)
    out[:len(a)] += a
    out[:len(b)] += b
    return out

def test_exact_double_glz(name, k_steps):
    """k_steps = list of (m_k, lam_k) for prod D_{m_k}(D_{m_k} + lam_k * y)(1)."""
    # Use sympy exact rationals!
    kmax = len(k_steps)
    F = [[None]*(k+1) for k in range(kmax+1)]
    F[0][0] = [sp.Rational(1)]
    def apply_D_q(c, m, poly):
        d = len(poly) - 1
        res = [sp.Rational(0)] * (d + 2)
        for r, coef in enumerate(poly):
            res[r] += coef * (r + c)
            res[r + 1] += coef * (r + m)
        return res
    def add_q(a, b):
        n = max(len(a), len(b))
        res = [sp.Rational(0)] * n
        for i, v in enumerate(a): res[i] += v
        for i, v in enumerate(b): res[i] += v
        return res
    for k in range(1, kmax + 1):
        m, lam = k_steps[k - 1]
        U = [None]*(k+1)
        U[0] = apply_D_q(1, m, F[k-1][0])
        for i in range(1, k):
            U[i] = add_q(apply_D_q(1, m, F[k-1][i]), [lam * v for v in F[k-1][i-1]])
        U[k] = [lam * v for v in F[k-1][k-1]]
        for i in range(k+1):
            F[k][i] = apply_D_q(1, m, U[i])
    # Check RR at level kmax for all i!
    fails = []
    for i in range(1, kmax):
        P = sp.Poly(list(reversed(F[kmax][i])), t, domain='QQ')
        rc = int(sp.count_roots(P, -1, 0))
        if rc < P.degree():
            fails.append((i, rc, P.degree()))
    print(f"  {name:52s}: RR fails = {fails}")

print("=== 1. ATTACK DOUBLE-GLZ WITH BIMODAL lam_k (mimicking (13,7) barrier!) ===")
# At k=14, take lam = [1, 0, 0, ..., 0, 1] -> wait: if middle lam = 0, then ONLY i = 0, 1, 2 exist!
# To get ONLY early-deletion S={1..7} and late-deletion S={8..14}, we can't with rank-1 weights unless
# we have a=1..2 on each end and huge m jump in middle!
# Wait: with lam = [1, 0, 0, 0, 0, 0, 0, 1], i=1 is literally A_{8, {1}} + A_{8, {8}}!!
# CHECK A_{k, {1}} + alpha * A_{k, {k}} FOR k = 14, OR A_{k, {1,2}} + ...!
# Wait: for i=1, Ma-Wang ALREADY proved EVERY positive combination of {A_{k,{j}}} is real-rooted!
# For i=4,5,6,7, to put huge weight on {1..i} and {k-i+1..k}, set lam = [1]*a + [0]*b + [1]*a with a=4, b=4 (k=12), i=4!
# Even simpler: try lam = [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1] (k=12) — ONLY 4+4=8 slots have lam=1, separated by 4 forced-keep steps!
test_exact_double_glz("GLZ k=10: lam=[1,1,1,0,0,0,0,1,1,1], m=k+1", [(k+1, 1 if (k<=3 or k>=8) else 0) for k in range(1, 11)])
test_exact_double_glz("GLZ k=12: lam=[1]*4+[0]*4+[1]*4, m=k+1", [(k+1, 1 if (k<=4 or k>=9) else 0) for k in range(1, 13)])
test_exact_double_glz("GLZ k=8: m=[2,2,10^4,10^4,10^4,10^4,2,2], lam=[1,1,0,0,0,0,1,1]",
                      [(2,1),(2,1),(10**4,0),(10**4,0),(10**4,0),(10**4,0),(2,1),(2,1)])
test_exact_double_glz("GLZ k=8: m=[10^4,10^4,2,2,2,2,10^4,10^4], lam=[1,1,0,0,0,0,1,1]",
                      [(10**4,1),(10**4,1),(2,0),(2,0),(2,0),(2,0),(10**4,1),(10**4,1)])

print("\n=== 2. ATTACK BARE ATOMIC (SINGLE D_j) WITH EXTREME BIMODAL lam_j AND m_j ===")
def test_exact_bare(name, steps, op_type="GLZ"):
    N = len(steps)
    F = [[None]*(k+1) for k in range(N+1)]
    F[0][0] = [sp.Rational(1)]
    def apply_op(c, m, poly):
        d = len(poly) - 1
        res = [sp.Rational(0)] * (d + 2)
        for r, coef in enumerate(poly):
            res[r] += coef * (r + c)
            res[r + 1] += coef * (r + m if op_type == "GLZ" else m)
        return res
    def add_q(a, b):
        n = max(len(a), len(b))
        res = [sp.Rational(0)] * n
        for i, v in enumerate(a): res[i] += v
        for i, v in enumerate(b): res[i] += v
        return res
    for k in range(1, N + 1):
        c, m, lam = steps[k - 1]
        F[k][0] = apply_op(c, m, F[k-1][0])
        for i in range(1, k):
            F[k][i] = add_q(apply_op(c, m, F[k-1][i]), [lam * v for v in F[k-1][i-1]])
        F[k][k] = [lam * v for v in F[k-1][k-1]]
    fails = []
    for i in range(1, N):
        if not any(v != 0 for v in F[N][i]): continue
        P = sp.Poly(list(reversed(F[N][i])), t, domain='QQ')
        rc = int(sp.count_roots(P, -1, 0)) if op_type == "GLZ" else int(sp.count_roots(P))
        if rc < P.degree():
            fails.append((i, rc, P.degree()))
    print(f"  {name:55s}: RR fails = {fails}")

# Test bare GLZ & Laguerre with 2 deletions separated by 8 forced-keeps with huge m jump!
# Wait: with lam=[1, 0, 0, 0, 0, 0, 0, 0, 0, 1], i=1 is D_10...D_2(1) + D_10...D_1(1)? No:
# j=1 has lam=1, j=10 has lam=1, middle j=2..9 have lam=0:
# [y^1] = D_10 D_9...D_2 (1) + D_9...D_1 (1)!!
# WAIT A SECOND!! CHECK THAT!!
# If D_1 = D_2 = ... = D_9 = D_2 and D_10 = D_{10^6}:
# [y^1] = D_{10^6} D_2^8(1) + D_2^9(1)!! (They share D_2^8(1), so always real-rooted!)
# To avoid sharing, put TWO at start and TWO at end: lam = [1, 1, 0, 0, 0, 0, 0, 0, 1, 1] and check i = 2!!
# ALSO weight left end by 1 and right end by alpha = 1000 (or 1/1000) to match magnitudes!!
for alpha in [sp.Rational(1), sp.Rational(1000), sp.Rational(1, 1000), sp.Rational(10**6), sp.Rational(1, 10**6)]:
    test_exact_bare(f"Bare GLZ N=10, 2+6(0)+2, m=2/10^4, al={alpha}",
                    [(1, 2, 1), (1, 2, 1)] + [(1, 10**4, 0)]*6 + [(1, 2, alpha), (1, 2, alpha)], "GLZ")
    test_exact_bare(f"Bare Lag N=10, 2+6(0)+2, m=1/10^4, al={alpha}",
                    [(1, 1, 1), (1, 1, 1)] + [(1, 10**4, 0)]*6 + [(1, 1, alpha), (1, 1, alpha)], "Lag")
