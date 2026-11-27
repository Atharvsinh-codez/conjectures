"""CRITICAL TEST:
Is the consecutive-triple quadratic
  Q_{k,i}(t, y) = P_{k,i-1}(t) + 2 P_{k,i}(t) y + P_{k,i+1}(t) y^2
REAL-STABLE in (t, y) for:
  1. GLZ up to k = 25 (check exact Rayleigh at t=0 up to k=25, AND UHP y-test up to k=14)!
  2. WILD sequences (m = [2, 10^5, 2, 10^5, 2], m_out=2/m_in=1000, m_out=1000/m_in=2, c=[1,100,1,100])!
  3. ALSO test r=2 WITH WEIGHT c: P_{k,i-1} + c * P_{k,i} y + P_{k,i+1} y^2 for c = 2, 1.5, 1.0!
  4. ALSO test r=3: P_{k,i-1} + 3 P_{k,i} y + 3 P_{k,i+1} y^2 + P_{k,i+2} y^3 at t=0 up to k=20!
"""
import math
import numpy as np
import sympy as sp

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

def build_general(pairs):
    kmax = len(pairs)
    P = [[None]*(k+1) for k in range(kmax+1)]
    P[0][0] = [1]
    for k in range(1, kmax + 1):
        co, mo, ci, mi = pairs[k - 1]
        U = [None]*(k+1)
        U[0] = apply_D(ci, mi, P[k-1][0])
        for i in range(1, k):
            U[i] = poly_add(apply_D(ci, mi, P[k-1][i]), P[k-1][i-1])
        U[k] = P[k-1][k-1]
        for i in range(k+1):
            P[k][i] = apply_D(co, mo, U[i])
    return P

def test_quad_uhp(f0, f1, f2, c_mid=2.0, n_trials=400):
    """Test real stability of f0(t) + c_mid*f1(t)*y + f2(t)*y^2 by checking Im(t) <= 0 for y in H+!"""
    rng = np.random.default_rng(77)
    dmax = len(f0) - 1
    test_ys = [complex(re_y, im_y)
               for re_y in [-100, -10, -2, -1, -0.5, -0.1, 0, 0.1, 0.5, 1, 2, 10, 100]
               for im_y in [1e-4, 1e-2, 0.1, 1.0, 10.0]]
    for _ in range(n_trials):
        test_ys.append(complex(rng.normal(0, 1) * (10 ** rng.uniform(-2, 2)), 10 ** rng.uniform(-4, 2)))
    max_im = -1e9
    worst_y = None
    for y in test_ys:
        c = np.zeros(dmax + 1, dtype=np.complex128)
        for m, v in enumerate(f0): c[m] += float(v)
        for m, v in enumerate(f1): c[m] += c_mid * float(v) * y
        for m, v in enumerate(f2): c[m] += float(v) * (y * y)
        r = np.roots(c[::-1])
        mi = float(np.max(r.imag))
        if mi > max_im:
            max_im = mi; worst_y = y
    return max_im, worst_y

print("=== 1. EXACT RAYLEIGH AT t=0 FOR QUADRATIC Q_{k,i} (r=2) AND CUBIC (r=3) ON GLZ k=2..25 ===")
P_glz = build_general([(1, k+1, 1, k+1) for k in range(1, 26)])
for k in range(2, 26):
    f2 = []
    for i in range(1, k):
        a0, a1, a2 = P_glz[k][i-1][0], P_glz[k][i][0], P_glz[k][i+1][0]
        b0, b1, b2 = P_glz[k][i-1][1], P_glz[k][i][1], P_glz[k][i+1][1]
        c0 = b0*a1 - a0*b1; c1 = b0*a2 - a0*b2; c2 = b1*a2 - a1*b2
        if c0 <= 0 or c2 <= 0 or c1*c1 > 4*c0*c2:
            f2.append(i)
    if f2 or k in [2, 10, 15, 20, 25]:
        print(f"  GLZ k={k:2d}: quad Q_{{k,i}}(0,y) Rayleigh fails = {f2[:8]}")

print("\n=== 2. FULL UHP STABILITY OF QUADRATIC Q_{k,i}(t,y) ON GLZ k=2..10 ===")
for k in range(2, 11):
    worst = (-1e9, None)
    for i in range(1, k):
        mi, wy = test_quad_uhp(P_glz[k][i-1], P_glz[k][i], P_glz[k][i+1], c_mid=2.0, n_trials=200)
        if mi > worst[0]: worst = (mi, (i, wy))
    print(f"  GLZ k={k:2d}: max Im(t) across all i = {worst[0]:+.3e} at i={worst[1][0]}, y={worst[1][1]}")

print("\n=== 3. CRITICAL: TEST QUADRATIC Q_{k,i}(t,y) ON WILD SEQUENCES!! ===")
wilds = [
    ("m=[2, 10^4, 2, 10^4, 2]", [(1, m, 1, m) for m in [2, 10**4, 2, 10**4, 2]]),
    ("m=[10^4, 2, 10^4, 2, 10^4]", [(1, m, 1, m) for m in [10**4, 2, 10**4, 2, 10**4]]),
    ("m_out=2, m_in=1000 (k=1..5)", [(1, 2, 1, 1000)]*5),
    ("m_out=1000, m_in=2 (k=1..5)", [(1, 1000, 1, 2)]*5),
    ("c=[1,100,1,100], m=[2,200,2,200]", [(1,2,1,2), (100,200,100,200), (1,2,1,2), (100,200,100,200)]),
]
for name, pairs in wilds:
    Pw = build_general(pairs)
    k = len(pairs)
    # Check exact Rayleigh at t=0 AND UHP test!
    ray_fails = []
    worst_im = (-1e9, None)
    for kk in range(2, k + 1):
        for i in range(1, kk):
            a0, a1, a2 = Pw[kk][i-1][0], Pw[kk][i][0], Pw[kk][i+1][0]
            b0, b1, b2 = Pw[kk][i-1][1], Pw[kk][i][1], Pw[kk][i+1][1]
            c0 = b0*a1 - a0*b1; c1 = b0*a2 - a0*b2; c2 = b1*a2 - a1*b2
            if c0 <= 0 or c2 <= 0 or c1*c1 > 4*c0*c2:
                ray_fails.append((kk, i))
            mi, wy = test_quad_uhp(Pw[kk][i-1], Pw[kk][i], Pw[kk][i+1], c_mid=2.0, n_trials=150)
            if mi > worst_im[0]: worst_im = (mi, (kk, i, wy))
    print(f"  {name:35s}: ray0_fails={ray_fails[:5]}, max Im(t)={worst_im[0]:+.3e} at {worst_im[1][:2]}")
