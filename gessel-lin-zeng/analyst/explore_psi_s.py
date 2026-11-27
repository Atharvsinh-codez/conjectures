"""Analyze the ultra-clean s-space recurrence:
  psi_{k,i}(s) = (1 / (1 + e^s)) * (psi''_{k-1, i}(s) + psi'_{k-1, i-1}(s)),  psi_{0,0}(s) = e^s/(1+e^s),
  where psi_{k,i}(s) = e^s / (1+e^s)^{3k-i+1} * A_{k,i}(-e^s).
1. Inspect zero locations s_m = log(-x_m) of psi_{k,i}(s) for k = 10, 15, 20, 25!
2. Check Turan (psi'_{k,i})^2 - psi_{k,i-1} psi''_{k,i+1} > 0 and row LC psi_{k,i}^2 - psi_{k,i-1} psi_{k,i+1} > 0
   across s in R for k = 11..20!
3. Check the EXACT Wronskian & Cauchy-index structure of (psi_{k,0}, ..., psi_{k,k}) in s-space!
4. Check if A_{k,i}(x) stays real-rooted and interlacing at k = 25, 30 (central slice i = k//2)!
"""
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

A = build_A(25)

def eval_psi_and_derivs(k, i, s_arr):
    """Evaluate psi_{k,i}(s), psi'_{k,i}(s), psi''_{k,i}(s) via A_{k,i}(-e^s).
    Note: psi_{k,i}(s) = e^s (1+e^s)^{-(3k-i+1)} A_{k,i}(-e^s)
    And psi'_{k,i}(s) = e^s (1+e^s)^{-(3k-i+2)} (S_{3k-i} A_{k,i})(-e^s)!!
    WAIT A SECOND! Check that!
    Let's verify: (d/ds)[e^s (1+e^s)^{-M} f(-e^s)] with u = e^s, d/ds = u d/du:
    = u [ (1+u)^{-M} f(-u) - M u (1+u)^{-M-1} f(-u) - u (1+u)^{-M} f'(-u) ]
    = u (1+u)^{-M-1} [ (1 + (1-M)u) f(-u) - u(1+u) f'(-u) ]
    Set x = -u: 1 + (M-1)x + x(1-x) f'(x) = (S_{M-1} f)(x)!!
    With M = 3k - i + 1, M - 1 = 3k - i!!
    EXACTLY!!
      psi_{k,i}(s)   <-> A_{k,i}(x)              with M = 3k-i+1
      psi'_{k,i}(s)  <-> S_{3k-i} A_{k,i}(x)     with M = 3k-i+2
      psi''_{k,i}(s) <-> S_{3k-i+1} S_{3k-i} A_{k,i}(x) with M = 3k-i+3!
    (CHECK! At k-1, M = 3(k-1)-i+1 = 3k-2-i:
      psi'_{k-1,i}  <-> S_{3k-3-i} A_{k-1,i}
      psi''_{k-1,i} <-> S_{3k-2-i} S_{3k-3-i} A_{k-1,i}
      psi'_{k-1,i-1} <-> S_{3k-2-i} A_{k-1,i-1}
      SUM = S_{3k-2-i}(S_{3k-3-i} A_{k-1,i} + A_{k-1,i-1}) = A_{k,i}!!
      100% CONSISTENT!!)
    """
    pass

# Check Turan (psi'_{k-1, i-1})^2 - psi_{k-1, i-2} * psi''_{k-1, i} > 0 on s in R for k = 11..20!
# In x = -e^s < 0, factoring out e^{2s} (1+e^s)^{-(6k-2i-2)}:
#   psi'_{k-1, i-1} <-> S_{3k-2-i} A_{k-1, i-1}      (M = 3k-i)
#   psi_{k-1, i-2}  <-> A_{k-1, i-2}                 (M = 3k-i)
#   psi''_{k-1, i}  <-> S_{3k-2-i} S_{3k-3-i} A_{k-1, i}  (M = 3k-i)
# ALL THREE HAVE THE EXACT SAME DENOMINATOR POWER M = 3k-i!!
# So (psi'_{k-1,i-1})^2 - psi_{k-1,i-2} * psi''_{k-1,i} > 0 on s in R
#  <=>  (S_L A_{k-1,i-1})^2 - A_{k-1,i-2} * (S_L S_{L-1} A_{k-1,i}) > 0 on x < 0 (where L = 3k-2-i)!

def eval_poly_mp(poly, x):
    acc = mp.mpf(0)
    for c in reversed(poly):
        acc = acc * x + c
    return acc

print("=== 1. CHECK TURAN M^2 - X*Y > 0 ON x < 0 FOR k = 12, 16, 20, 25 ===")
mp.mp.dps = 50
# Sample s across the entire zero range [-2*k*log(2), +2*k*log(2)]!
for k in [12, 16, 20, 25]:
    min_ratio = (1e300, None)
    for i in [2, k//4, k//2, 3*k//4, k-1]:
        if i < 2 or i >= k: continue
        L = 3*k - 2 - i
        X = A[k-1][i-2]
        M = apply_S(L, A[k-1][i-1])
        Y = apply_S(L, apply_S(L - 1, A[k-1][i]))
        # Test 200 values of s from -2*k to +2*k
        for s in np.linspace(-1.8*k, 1.5*k, 200):
            x = -mp.exp(mp.mpf(float(s)))
            vx = eval_poly_mp(X, x)
            vm = eval_poly_mp(M, x)
            vy = eval_poly_mp(Y, x)
            diff = vm*vm - vx*vy
            denom = vm*vm + abs(vx*vy) + mp.mpf('1e-100')
            rel = float(diff / denom)
            if rel < min_ratio[0]:
                min_ratio = (rel, (i, float(s)))
    print(f"  k={k:2d}: min (M^2 - X*Y)/(M^2 + |XY|) = {min_ratio[0]:+.5e} at i={min_ratio[1][0]}, s={min_ratio[1][1]:.2f}")

print("\n=== 2. INSPECT ZERO DISTRIBUTION IN s = log(-x) FOR k = 10, 20 (i = 0, k//2, k) ===")
for k in [10, 20]:
    for i in [0, k//2, k]:
        poly = A[k][i]
        d = len(poly) - 1
        g = (mp.mpf(poly[0]) / mp.mpf(poly[-1])) ** (mp.mpf(1) / d)
        scaled = [mp.mpf(c) * (g ** m) for m, c in enumerate(poly)]
        r = mp.polyroots(list(reversed(scaled)), maxsteps=200, extraprec=40)
        xs = sorted(mp.re(z) * g for z in r) # ascending (most neg first)
        ss = sorted(float(mp.log(-x)) for x in xs) # ascending s
        gaps = [ss[m+1] - ss[m] for m in range(len(ss) - 1)]
        print(f"  k={k:2d}, i={i:2d} (d={d:2d}): s in [{ss[0]:.3f}, {ss[-1]:.3f}] (div k: [{ss[0]/k:.3f}, {ss[-1]/k:.3f}]), "
              f"gaps min={min(gaps):.3f}, med={np.median(gaps):.3f}, max={max(gaps):.3f}")
