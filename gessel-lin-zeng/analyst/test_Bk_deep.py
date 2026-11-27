"""Test:
1. B_k(x, y) = sum binom(k,i) A_{k,i}(x) y^i AND sum binom(k,i) htilde_{k,i}(t) y^i stability up to k=14!
2. CRITICAL: For B_k(x, y), test the EXACT quadratic-form / Hermite-Biehler / UHP condition near real axis!
   Specifically: P(x, y) is real-stable iff for all real y, P(x, y) is real-rooted in x AND
   Wronskian W_x(P(x, y), P_y(x, y)) = P_x(x,y) P_y(x,y) - P(x,y) P_{xy}(x,y) >= 0 (or <= 0)!
   WAIT A SECOND!! STOP!!
   For P(x, y) real-stable, for EVERY real x_0, y -> P(x_0, y) is real-rooted, AND for EVERY real y_0,
   x -> P(x, y_0) is real-rooted, AND the differential Hermite-Biehler / Wronskian condition
   P_x P_y - P P_{xy} >= 0 for ALL real (x, y) is NECESSARY (and with real-rootedness on one slice, SUFFICIENT for degree <= ... wait:
   Let's check P_x P_y - P P_{xy} on ALL real (x, y)!)!
3. Also test what happens to sum c_i A_{k,i}(x) y^i for c_i = binom(k,i)^alpha, alpha in [0.1, 0.2, ..., 1.0]!
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

def apply_Dtilde(j, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(poly):
        res[m] += c * (m + 1)
        res[m + 1] += c * (m + j + 1)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def build_htilde(kmax):
    p = [[None] * (k + 1) for k in range(kmax + 1)]
    p[1][0] = [1, 2]
    p[1][1] = [1]
    for k in range(2, kmax + 1):
        p[k][0] = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][0]))
        for i in range(1, k):
            p[k][i] = apply_Dtilde(k, poly_add(apply_Dtilde(k, p[k - 1][i]), p[k - 1][i - 1]))
        p[k][k] = apply_Dtilde(k, p[k - 1][k - 1])
    return p

A = build_A(14)
H = build_htilde(10)

print("=== 1. TEST B_k(x, y) = sum binom(k,i) A_{k,i}(x) y^i IN y AT MANY REAL x0 (k=2..14) ===")
for k in range(2, 15):
    worst = (0.0, None)
    xs = list(np.linspace(-20, 20, 401)) + [-1000, -100, -50, 50, 100, 1000] + list(np.linspace(-3, 3, 601))
    for x0 in xs:
        ycoeffs = [math.comb(k, i) * sum(float(c) * (x0**m) for m, c in enumerate(A[k][i])) for i in range(k + 1)]
        if abs(ycoeffs[-1]) < 1e-30:
            continue
        r = np.roots(list(reversed(ycoeffs)))
        mi = float(np.max(np.abs(r.imag)))
        if mi > worst[0]:
            worst = (mi, x0)
    print(f"  k={k:2d}: max |Im(y)| across {len(xs)} real x0 = {worst[0]:.3e} (at x0={worst[1]})")

print("\n=== 2. TEST Htilde_B_k(t, y) = sum binom(k,i) htilde_{k,i}(t) y^i IN y AT REAL t0 (k=2..10) ===")
for k in range(2, 11):
    worst = (0.0, None)
    ts = list(np.linspace(-5, 5, 401)) + [-50, 50]
    for t0 in ts:
        ycoeffs = [math.comb(k, i) * sum(float(c) * (t0**m) for m, c in enumerate(H[k][i])) for i in range(k + 1)]
        if abs(ycoeffs[-1]) < 1e-30:
            continue
        r = np.roots(list(reversed(ycoeffs)))
        mi = float(np.max(np.abs(r.imag)))
        if mi > worst[0]:
            worst = (mi, t0)
    print(f"  k={k:2d}: max |Im(y)| across {len(ts)} real t0 = {worst[0]:.3e} (at t0={worst[1]})")

print("\n=== 3. CRITICAL: TEST CRITICAL THRESHOLD alpha FOR sum binom(k,i)^alpha A_{k,i}(x) y^i AT x=1, x=-1, x=0 ===")
for k in [2, 3, 4, 5, 8, 12]:
    for alpha in [0.0, 0.25, 0.5, 0.75, 1.0]:
        bad = False
        for x0 in [1.0, 2.0, 10.0, -0.5, -2.0]:
            ycoeffs = [(math.comb(k, i)**alpha) * sum(float(c) * (x0**m) for m, c in enumerate(A[k][i])) for i in range(k + 1)]
            r = np.roots(list(reversed(ycoeffs)))
            if np.max(np.abs(r.imag)) > 1e-6:
                bad = True
                break
        if alpha == 0.0 or not bad:
            pass
    # Find min alpha where x0 in [-5..50] all real-rooted
    lo, hi = 0.0, 2.0
    for _ in range(25):
        mid = 0.5 * (lo + hi)
        bad = False
        for x0 in [0.01, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0, -0.5, -1.0, -2.0, -5.0]:
            ycoeffs = [(math.comb(k, i)**mid) * sum(float(c) * (x0**m) for m, c in enumerate(A[k][i])) for i in range(k + 1)]
            r = np.roots(list(reversed(ycoeffs)))
            if np.max(np.abs(r.imag)) > 1e-5:
                bad = True
                break
        if bad: lo = mid
        else: hi = mid
    print(f"  k={k:2d}: min alpha for real roots in y ~ {hi:.4f}")
