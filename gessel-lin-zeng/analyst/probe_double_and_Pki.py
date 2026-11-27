import numpy as np
import mpmath as mp
from fractions import Fraction

# 1. Test double-operator T_j(T_j + w_j y)(1) with @algebraist's EXACT k=9 counterexample parameters
#    AND random extreme (a_j, b_j, w_j) attacks!
#    Recall T_j f = (a_j + b_j t) f + t(1+t) f' (deg +1), and @algebraist's k=9 single-operator disproof had
#    a9 = 1/1000, b9 = 1000000001/1000, w9 = 1000 (after 8 nice steps).
#    Let's test BOTH:
#      (A) prod_{j=1}^k T_j (T_j + w_j y) (1)  with varying (a_j, b_j, w_j)
#      (B) prod_{j=1}^k D_{m_j} (D_{m_j} + w_j y) (1) with wild m_j > 1 AND wild w_j > 0!
#      (C) prod_{j=1}^k D_{m_j} (D_{n_j} + y) (1) with m_j != n_j!
#      (D) prod_{j=1}^k (D_{m_j}^2 + w_j y D_{m_j}) (1) -- same as (B)!

def apply_T(p, a, b):
    # p in ascending powers of t (mpmath mpf)
    d = len(p) - 1
    res = [mp.mpf(0)] * (d + 2)
    for n in range(d + 1):
        res[n] += (a + n) * p[n]
        res[n + 1] += (b + n) * p[n]
    return res

def add_poly(p, q):
    n = max(len(p), len(q))
    res = [mp.mpf(0)] * n
    for i, v in enumerate(p): res[i] += v
    for i, v in enumerate(q): res[i] += v
    return res

def scale_poly(p, s):
    return [v * s for v in p]

def check_real_rooted_mp(p, dps=60):
    # p ascending; strip trailing zeros
    while len(p) > 1 and p[-1] == 0:
        p = p[:-1]
    deg = len(p) - 1
    if deg <= 1:
        return True, 0.0
    mp.mp.dps = dps
    # Use numpy for fast check first, then mpmath if suspicious
    arr = np.array([float(v) for v in p[::-1]], dtype=np.float64)
    if np.all(np.isfinite(arr)) and arr[0] != 0:
        arr = arr / np.max(np.abs(arr))
        try:
            r = np.roots(arr)
            max_im = float(np.max(np.abs(r.imag)))
            if max_im < 1e-4:
                return True, max_im
        except Exception:
            pass
    # High-precision mpmath polyroots
    try:
        roots = mp.polyroots(p[::-1], maxsteps=200, extraprec=40)
        max_im = float(max(abs(mp.im(z)) for z in roots))
        return max_im < 1e-10, max_im
    except Exception as e:
        return None, -1.0

def test_double_family(name, steps_list):
    # steps_list: list of (a1, b1, a2, b2, w) for T1(T2 + w y)
    # Wait: P_k(t, y) = sum_i p_{k,i}(t) y^i
    # P_{k,i} = T1( T2 P_{k-1, i} + w P_{k-1, i-1} )
    mp.mp.dps = 60
    P = [[mp.mpf(1)]]
    for k, (a1, b1, a2, b2, w) in enumerate(steps_list, 1):
        newP = []
        for i in range(k + 1):
            t2 = apply_T(P[i], a2, b2) if i < len(P) else [mp.mpf(0)]
            wy = scale_poly(P[i - 1], mp.mpf(w)) if i >= 1 else [mp.mpf(0)]
            newP.append(apply_T(add_poly(t2, wy), a1, b1))
        P = newP
    fails = []
    for i, poly in enumerate(P):
        ok, im = check_real_rooted_mp(poly, dps=80)
        if ok is False:
            fails.append((i, len(poly)-1, im))
    print(f"  {name:50s}: k={len(steps_list)}, fails={fails}")
    return fails

print("=== 1. ATTACKING DOUBLE-OPERATOR T1(T2 + w y)(1) ===")
# Attack A: @algebraist's k=9 extreme final step!
# First 8 steps GLZ-like (a=1, b=j+1, w=1), 9th step (a=1e-3, b=1e6, w=1e3)
s_alg9 = [(1, j+1, 1, j+1, 1) for j in range(1, 9)] + [(0.001, 1e6, 0.001, 1e6, 1e3)]
test_double_family("GLZ8 + alg9_extreme(T(T+wy))", s_alg9)

# Attack B: wild a, b, w on SAME D_m(D_m + w y) with a=1, wild w!
s_wild_w = [(1, 2, 1, 2, 1)]*5 + [(1, 1000, 1, 1000, 1e6)]*4
test_double_family("D_m(D_m + w y) wild w=1e6 (k=9)", s_wild_w)

s_wild_w2 = [(1, 2, 1, 2, 1)]*6 + [(1, 1.001, 1, 1.001, 1e8)] + [(1, 1e6, 1, 1e6, 1e-6)]*2
test_double_family("D_m(D_m + w y) opposing wild w (k=9)", s_wild_w2)

# Attack C: T(T + w y) with a << 1, b >> 1, w >> 1 at varying positions
s_abw1 = [(1, 2, 1, 2, 1)]*4 + [(1e-4, 1e5, 1e-4, 1e5, 1e5)]*4
test_double_family("T(T+wy) small a=1e-4, big b,w=1e5 (k=8)", s_abw1)

s_abw2 = [(1e-3, 1e4, 1e-3, 1e4, 1e4)]*3 + [(1, 2, 1, 2, 1)]*5
test_double_family("T(T+wy) early extreme then nice (k=8)", s_abw2)

s_abw3 = [(1, 2, 1, 2, 1)]*6 + [(1e-5, 1e6, 1e-5, 1e6, 1e6)] + [(1, 2, 1, 2, 1)]*2
test_double_family("T(T+wy) extreme at step 7 then nice (k=9)", s_abw3)

# Attack D: D_{m1}(D_{m2} + y) with m1 != m2!
s_diff_m = [(1, 2, 1, 1e6, 1)]*4 + [(1, 1e6, 1, 2, 1)]*4
test_double_family("D_{m1}(D_{m2}+y) alternating 2/1e6 (k=8)", s_diff_m)

s_diff_m2 = [(1, 2, 1, 2, 1)]*5 + [(1, 1.0001, 1, 1e6, 1)]*3
test_double_family("D_{m1}(D_{m2}+y) m1=1.0001, m2=1e6 (k=8)", s_diff_m2)

# Attack E: D_m(D_m + y) with EXTREME step-jumps and larger k=12!
s_Dm12 = [(1, 1.0001, 1, 1.0001, 1)]*6 + [(1, 1e8, 1, 1e8, 1)]*6
test_double_family("D_m(D_m+y) 6x1.0001 + 6x1e8 (k=12)", s_Dm12)

s_Dm12_rev = [(1, 1e8, 1, 1e8, 1)]*6 + [(1, 1.0001, 1, 1.0001, 1)]*6
test_double_family("D_m(D_m+y) 6x1e8 + 6x1.0001 (k=12)", s_Dm12_rev)

# Attack F: randomized extreme T(T+wy) with 0 < a < b, w > 0 at k=9 (20 trials)
rng = np.random.default_rng(777)
rand_fails = 0
for tr in range(20):
    steps = []
    for j in range(9):
        a = 10**rng.uniform(-4, 2)
        b = a + 10**rng.uniform(-3, 5)
        w = 10**rng.uniform(-4, 5)
        steps.append((a, b, a, b, w))
    mp.mp.dps = 60
    P = [[mp.mpf(1)]]
    for a1, b1, a2, b2, w in steps:
        newP = []
        for i in range(len(P) + 1):
            t2 = apply_T(P[i], a2, b2) if i < len(P) else [mp.mpf(0)]
            wy = scale_poly(P[i-1], mp.mpf(w)) if i >= 1 else [mp.mpf(0)]
            newP.append(apply_T(add_poly(t2, wy), a1, b1))
        P = newP
    bad = False
    for i in range(2, 7):
        ok, im = check_real_rooted_mp(P[i], dps=80)
        if ok is False:
            bad = True
            print(f"  RANDOM T(T+wy) FAIL at trial {tr}, i={i}, max_im={im:.3e}!")
            break
    if bad:
        rand_fails += 1
print(f"  Random extreme T(T+wy) k=9: {rand_fails}/20 failed")

print("\n=== 2. ROOTS OF THE EHRHART / JACOBI-STIRLING POLYNOMIAL P_{k,i}(N) ===")
# Recall sum_{N>=0} P_{k,i}(N) x^N = x A_{k,i}(x) / (1-x)^{3k-i+1}.
# Equivalently P_{k,i}(N) = sum_{r=0}^{2k-i-1} a_r * C(N - 1 - r + 3k - i, 3k - i)
import sympy as sp
N = sp.Symbol('N')

def apply_S_int(p, L):
    d = len(p) - 1
    res = [0]*(d + 2)
    for n in range(d + 2):
        v = 0
        if n <= d: v += (n + 1) * p[n]
        if 1 <= n <= d + 1: v += (L - n + 1) * p[n - 1]
        res[n] = v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def poly_add_int(a, b):
    n = max(len(a), len(b))
    res = [0]*n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

A = {(0, 0): [1]}
for k in range(1, 6):
    for i in range(k + 1):
        u1 = apply_S_int(A[(k-1, i)], 3*k - 3 - i) if i <= k - 1 else [0]
        u2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S_int(poly_add_int(u1, u2), 3*k - 2 - i)

for k in [1, 2, 3, 4]:
    for i in range(k + 1):
        D = 3*k - i
        # P(N) = sum_r A_{k,i}[r] * binomial(N - 1 - r + D, D)
        poly = 0
        for r, coeff in enumerate(A[(k, i)]):
            term = coeff
            for m in range(D):
                term *= (N - 1 - r + D - m)
            poly += term
        poly = sp.Poly(sp.expand(poly), N)
        # divide out N(N+1)...(N+k)
        q = poly
        for z in range(k + 1):
            q, rem = sp.div(q, sp.Poly(N + z, N))
            assert rem == 0
        roots = np.roots([float(c) for c in q.all_coeffs()])
        max_im = float(np.max(np.abs(roots.imag))) if len(roots) else 0.0
        print(f"  k={k}, i={i}: D={D}, trivial roots 0..-{k}, rem deg={q.degree()}, max|Im|={max_im:.4f}, roots={np.round(roots, 3)}")
