"""CRITICAL PROBING EXPERIMENT:
Why did Gaussian, Gumbel, AND Cauchy ALL pass real-rootedness, row interlacing, AND Pair(2,1)?!
Test:
1. Gaussian L_{m_k} = d/ds - m_k * s with:
   (a) m_k = k (up to k=12, ALSO check bivariate stability in (s, y)!!)
   (b) m_k = decreasing (e.g. 10, 9, 8, 7, 6, 5, 4)
   (c) m_k = wildly alternating (1, 100, 1, 100, 1, 100)
   (d) m_k = k^2 (super-fast growth: 1, 4, 9, 16, 25, 36)
2. GLZ operator D_{m_k} = t(1+t)d/dt + 1 + m_k*t with:
   (a) m_k = wildly alternating (2, 100, 2, 100, 2, 100)
   (b) m_k = decreasing (20, 18, 16, 14, 12, 10)
   (c) m_k = 1, 1, 1, 1, ... (constant!)
3. BIMODAL weight w(s) = (1 + 20 s^2) exp(-s^2/2) (which is NOT unimodal, NOT log-concave!):
   phi_{k,i} = w(s) (phi''_{k-1,i} + phi'_{k-1,i-1})!
4. ALSO check: for Gaussian L_m = d/ds - m*s, what IS the bivariate GF sum P_{k,i}(s) y^i?! Is IT real-stable in (s, y)??!
"""
import numpy as np
import sympy as sp

s = sp.Symbol('s')
t = sp.Symbol('t')

def test_seq(name, op_at_step, kmax, var=s):
    P = [[None] * (k + 1) for k in range(kmax + 1)]
    P[0][0] = sp.Integer(1)
    for k in range(1, kmax + 1):
        P[k][0] = sp.expand(op_at_step(k, op_at_step(k, P[k-1][0])))
        for i in range(1, k):
            P[k][i] = sp.expand(op_at_step(k, op_at_step(k, P[k-1][i]) + P[k-1][i-1]))
        P[k][k] = sp.expand(op_at_step(k, P[k-1][k-1]))
    rr_fails = []
    for k in range(1, kmax + 1):
        for i in range(k + 1):
            poly = sp.Poly(P[k][i], var, domain='QQ')
            d = poly.degree()
            if d == 0: continue
            rc = int(sp.count_roots(poly))
            if rc < d:
                rr_fails.append((k, i, rc, d))
    print(f"{name:55s}: RR fails = {rr_fails[:6]} (total {len(rr_fails)})")
    return P

# 1. Gaussian with various m_k
test_seq("Gauss m_k = [10,9,8,7,6,5,4]", lambda k, f: sp.diff(f, s) - (11 - k) * s * f, 7)
test_seq("Gauss m_k = [1,100,1,100,1,100]", lambda k, f: sp.diff(f, s) - (1 if k%2==1 else 100) * s * f, 6)
test_seq("Gauss m_k = [100,1,100,1,100,1]", lambda k, f: sp.diff(f, s) - (100 if k%2==1 else 1) * s * f, 6)
test_seq("Gauss m_k = k^2 [1,4,9,16,25,36,49]", lambda k, f: sp.diff(f, s) - (k*k) * s * f, 7)
test_seq("Gauss m_k = 2^k [2,4,8,16,32,64]", lambda k, f: sp.diff(f, s) - (2**k) * s * f, 6)
test_seq("Gauss m_k = [1, 1000, 1, 1000, 1, 1000]", lambda k, f: sp.diff(f, s) - (1 if k%2==1 else 1000) * s * f, 6)

# 2. GLZ D_{m_k} = t(1+t)d/dt + 1 + m_k*t with arbitrary m_k!
def D_op(m, f):
    return t * (1 + t) * sp.diff(f, t) + (1 + m * t) * f

test_seq("GLZ D_{m_k} with m_k = [2,100,2,100,2,100]", lambda k, f: D_op(2 if k%2==1 else 100, f), 6, var=t)
test_seq("GLZ D_{m_k} with m_k = [100,2,100,2,100,2]", lambda k, f: D_op(100 if k%2==1 else 2, f), 6, var=t)
test_seq("GLZ D_{m_k} with m_k = [10,9,8,7,6,5]", lambda k, f: D_op(11 - k, f), 6, var=t)
test_seq("GLZ D_{m_k} with m_k = k^2 [1,4,9,16,25,36]", lambda k, f: D_op(k*k, f), 6, var=t)
test_seq("GLZ D_{m_k} with m_k = 10^k [10,100,1000,...]", lambda k, f: D_op(10**k, f), 5, var=t)

# 3. Check if GAUSSIAN with m_k = k has bivariate stability of sum P_{k,i}(s) y^i!!
P_g = test_seq("Gauss m_k = k (k=1..6)", lambda k, f: sp.diff(f, s) - k * s * f, 6)
y = sp.Symbol('y')
for k in [2, 3, 4]:
    H_y = sum(P_g[k][i] * y**i for i in range(k + 1))
    # check at s = 0, s = 1, s = 2
    bads = []
    for s0 in [0, 1, 2, 5]:
        poly_y = sp.Poly(sp.expand(H_y.subs(s, s0)), y, domain='QQ')
        rc = int(sp.count_roots(poly_y))
        if rc < poly_y.degree():
            bads.append((s0, rc, poly_y.degree()))
    print(f"  Gauss plain sum P_{{k,i}}(s) y^i at k={k}: fails in y at {bads}")

# 4. BIMODAL weight w(s) = (1 + 20 s^2) exp(-s^2/2)!
# If phi_{k,i} = w^{k+1} * P_{k,i} / (1+20s^2)^{...} wait:
# phi_{0,0} = (1+20s^2) e^{-s^2/2}.
# phi' = ( (1+20s^2)P' + (40s - s(1+20s^2))P ) e^{-s^2/2}...
# Simpler: phi_{k,i}(s) = e^{-(k+1)s^2/2} * R_{k,i}(s), where
# R_{0,0} = 1 + 20 s^2, and R_{k,i} = (1 + 20 s^2) * [ L_k(L_k R_{k-1,i} + R_{k-1,i-1}) ]!
# Wait: L_k is (d/ds - k*s), which commutes with multiplying by (1+20s^2) ONLY AFTER differentiating!
# Exact: phi_{k-1,i} = e^{-k s^2/2} R_{k-1,i}, so phi'_{k-1,i} = e^{-k s^2/2} L_k(R_{k-1,i}),
# phi''_{k-1,i} = e^{-k s^2/2} L_k^2(R_{k-1,i}),
# then multiply by w(s) = (1+20s^2) e^{-s^2/2} => R_{k,i} = (1+20s^2) * L_k(L_k R_{k-1,i} + R_{k-1,i-1})!
R = [[None]*(6) for _ in range(6)]
R[0][0] = 1 + 20 * s**2
for k in range(1, 5):
    L = lambda f, kk=k: sp.diff(f, s) - kk * s * f
    R[k][0] = sp.expand((1 + 20*s**2) * L(L(R[k-1][0])))
    for i in range(1, k):
        R[k][i] = sp.expand((1 + 20*s**2) * L(L(R[k-1][i]) + R[k-1][i-1]))
    R[k][k] = sp.expand((1 + 20*s**2) * L(R[k-1][k-1]))
print("\n=== BIMODAL w(s) = (1 + 20 s^2) exp(-s^2/2): strip factor (1+20s^2) and check real roots ===")
for k in range(1, 5):
    for i in range(k + 1):
        Q = sp.Poly(R[k][i], s, domain='ZZ').exquo(sp.Poly(1 + 20*s**2, s, domain='ZZ'))
        rc = int(sp.count_roots(Q))
        if rc < Q.degree():
            print(f"  BIMODAL FAIL at k={k}, i={i}: real_roots={rc} < deg={Q.degree()}")
print("  (Done bimodal check k=1..4)")
