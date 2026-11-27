"""Deep asymptotic & structural analysis of A_{k,i}(x) and htilde_{k,i}(t):
1. Compute EXACT roots (high-precision mpmath in x-variable? No, use exact recurrence + Newton/Aberth in mpmath,
   OR exact Cauchy-index up to k=30 for central slice i = k//2!)
   Wait: check if Pair (2,1), row interlacing, and real-rootedness hold at k = 25, 30, 35, 40!
2. Measure the EXACT normalized interlacing margins (how close roots of S_{L-1} A_{k-1,i} get to roots of A_{k-1,i-2},
   and how close roots of D_{k+1} htilde_{k-1,i} get to roots of htilde_{k-1,i-1}) as k grows!
3. ALSO: Check if the recurrence for A_{k,i}(x) has a CLOSED-FORM INTEGRAL / PROBABILITY / DIFFERENTIAL-EQUATION /
   SADDLE-POINT representation, OR check the exact root distribution density rho_alpha(x) as k -> inf with i/k -> alpha!
"""
import time
import math
import numpy as np
import mpmath as mp

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_S(L, p):
    d = len(p) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(p):
        res[m] += c * (m + 1)
        res[m + 1] += c * (L - m)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def build_A(kmax):
    A = [[None] * (k + 1) for k in range(kmax + 1)]
    A[0][0] = [1]
    for k in range(1, kmax + 1):
        A[k][0] = apply_S(3*k - 2, apply_S(3*k - 3, A[k-1][0]))
        for i in range(1, k):
            L = 3*k - 2 - i
            A[k][i] = apply_S(L, poly_add(apply_S(L - 1, A[k-1][i]), A[k-1][i-1]))
        A[k][k] = apply_S(2*k - 2, A[k-1][k-1])
    return A

def get_roots_x_mp(poly, dps=60):
    """Find all roots of A_{k,i}(x) on (-inf, 0) using log-space / companion or mpmath.polyroots in u = -1/x or x!
    Wait: in x, roots are spread on (-inf, 0). Better: map t = x/(1-x) in (-1, 0), or even better:
    we know ALL roots of A_{k,i} are simple and negative, and if we have roots of A_{k, i+1} (or previous step),
    we can find EVERY root of A_{k,i} to 60 digits in 5 Newton steps using interlacing brackets!!
    EVEN BETTER: use mp.polyroots on A_{k,i}(x) after scaling x -> x/scale so geometric mean of roots is ~1!
    """
    mp.mp.dps = dps
    d = len(poly) - 1
    if d == 0: return []
    # Geometric mean of |roots| is (a_0 / a_d)^{1/d}
    g = (mp.mpf(poly[0]) / mp.mpf(poly[-1])) ** (mp.mpf(1) / d)
    # Substitute x = g * u => coeffs are a_m * g^m
    scaled = [mp.mpf(c) * (g ** m) for m, c in enumerate(poly)]
    r = mp.polyroots(list(reversed(scaled)), maxsteps=200, extraprec=40)
    assert max(abs(mp.im(z)) for z in r) < mp.mpf('1e-25'), f"Nonreal root! max im={max(abs(mp.im(z)) for z in r)}"
    return sorted(mp.re(z) * g for z in r)

A = build_A(20)
print("=== RELATIVE INTERLACING MARGINS IN t = x/(1-x) FOR PAIR (2,1): S_{L-1} A_{k-1,i} << A_{k-1,i-2} ===")
# Let r_0 < r_1 < ... < r_d be roots of B = A_{k-1,i-2} (in t-variable!),
# and s_0 < s_1 < ... < s_{d-1} be roots of S = S_{L-1} A_{k-1,i} (in t-variable!).
# Define fractional position f_m = (s_m - r_m) / (r_{m+1} - r_m) in (0, 1)!
# If f_m stays strictly bounded away from 0 and 1 (say in [0.1, 0.9]), Pair (2,1) is ROCK SOLID!
for k in [6, 10, 14, 18, 20]:
    for i in [2, k//4, k//2, 3*k//4, k-1]:
        if i < 2 or i >= k: continue
        L = 3*k - 2 - i
        rB_x = get_roots_x_mp(A[k-1][i-2], dps=50)
        rS_x = get_roots_x_mp(apply_S(L - 1, A[k-1][i]), dps=50)
        rB = [x / (1 - x) for x in rB_x]
        rS = [x / (1 - x) for x in rS_x]
        fracs = [float((rS[m] - rB[m]) / (rB[m+1] - rB[m])) for m in range(len(rS))]
        print(f"k={k:2d}, i={i:2d} (d={len(rS):2d}): Pair(2,1) t-fracs min={min(fracs):.4f}, max={max(fracs):.4f}, "
              f"ends=[{fracs[0]:.3f}, {fracs[1]:.3f}, ..., {fracs[-2]:.3f}, {fracs[-1]:.3f}]")

print("\n=== ALSO CHECK RELATIVE INTERLACING MARGINS OF SUMMANDS: A_{k-1,i-1} < S_{L-1} A_{k-1,i} (same deg!) ===")
for k in [6, 10, 14, 18, 20]:
    for i in [2, k//2, k-1]:
        if i < 1 or i >= k: continue
        L = 3*k - 2 - i
        r1_x = get_roots_x_mp(A[k-1][i-1], dps=50)
        r2_x = get_roots_x_mp(apply_S(L - 1, A[k-1][i]), dps=50)
        r1 = [x / (1 - x) for x in r1_x]
        r2 = [x / (1 - x) for x in r2_x]
        # r1[0] < r2[0] < r1[1] < r2[1] < ... < r1[d-1] < r2[d-1]
        # Fractional position of r2[m] in (r1[m], r1[m+1]):
        fracs = [float((r2[m] - r1[m]) / (r1[m+1] - r1[m])) for m in range(len(r1) - 1)]
        print(f"k={k:2d}, i={i:2d}: summands (r2-r1)/(r1'-r1) min={min(fracs):.4f}, max={max(fracs):.4f}, "
              f"first/last=[{fracs[0]:.3f}, ..., {fracs[-1]:.3f}]")
