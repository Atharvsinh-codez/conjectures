"""Probe the normalized equal-degree finite-difference hierarchy:
  ghat_{k, i}(t) = (1 / binom(k, i)) D_{k+1}^i htilde_{k, i}(t),   ghat_{k, i}(0) = 1.
  Delta_1(k, i) := ghat_{k, i} - ghat_{k, i-1}  (divisible by t!)
  Delta_2(k, i) := Delta_1(k, i+1) - Delta_1(k, i)  (divisible by t? or t^2?)

Check for k = 2..12:
1. Is Delta_1(k, i) / t real-rooted in (-1, 0)?
2. Does Delta_1(k, i) / t interlace BOTH ghat_{k, i-1} and ghat_{k, i}? (Automatically yes if ghat_{k,i-1} < ghat_{k,i}!)
3. Does Delta_1(k, i) / t INTERLACE Delta_1(k, i-1) / t?! (And in WHICH direction?!)
4. Does Delta_2(k, i) have ALL REAL ROOTS in (-1, 0), and is it divisible by t^2?!
5. CRITICAL TEST: What is the EXACT relation between Delta_1(k, i) and Wronskians?
   And ALSO: check Q_m(t, y) = htilde_{k-1, i-2} + 2 (D_m htilde_{k-1, i-1}) y + (D_m^2 htilde_{k-1, i}) y^2
   as a function of m (at k = 10, 20, 35): at what m does Q_m FIRST fail real stability?!
"""
import math
import sympy as sp
import mpmath as mp

t = sp.Symbol('t')
mp.mp.dps = 60

def apply_D(m, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def poly_lincomb(a, ca, b, cb):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += ca * v
    for i, v in enumerate(b): res[i] += cb * v
    return res

KMAX = 11
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D(m, H[k-1][0])
    for i in range(1, k): U[i] = poly_add(apply_D(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1): H[k][i] = apply_D(m, U[i])

def roots_mp(p):
    # strip leading/trailing zeros
    lo = 0
    while lo < len(p) and p[lo] == 0: lo += 1
    q = p[lo:]
    while len(q) > 1 and q[-1] == 0: q.pop()
    d = len(q) - 1
    if d <= 0: return lo, [], 0.0
    g = (abs(mp.mpf(q[0])) / abs(mp.mpf(q[-1]))) ** (mp.mpf(1) / d)
    cs = [mp.mpf(q[r]) * (g ** r) for r in range(d + 1)]
    rts = mp.polyroots(cs[::-1], maxsteps=400, extraprec=50)
    max_im = float(max(abs(mp.im(r)) for r in rts))
    reals = sorted([mp.re(r) * g for r in rts])
    return lo, reals, max_im

def check_interlace_equal(r1, r2):
    """Check if r1 < r2 (each r1[j] < r2[j] < r1[j+1]) or r2 < r1."""
    if len(r1) != len(r2): return "diff_deg"
    d = len(r1)
    if all(r1[j] < r2[j] for j in range(d)) and all(r2[j] < r1[j+1] for j in range(d-1)):
        return "r1<r2"
    if all(r2[j] < r1[j] for j in range(d)) and all(r1[j] < r2[j+1] for j in range(d-1)):
        return "r2<r1"
    return "NO"

print("=== 1. FINITE DIFFERENCES Delta_1 AND Delta_2 OF ghat_{k, i} (k = 2..10) ===")
for k in range(2, 11):
    m = k + 1
    # Exact integer multiple: scale ghat_{k,i} by k! so everything stays integer!
    kfact = math.factorial(k)
    Ghat = []
    for i in range(k+1):
        p = list(H[k][i])
        for _ in range(i): p = apply_D(m, p)
        scale = kfact // math.comb(k, i)
        Ghat.append([c * scale for c in p])
    D1 = [poly_lincomb(Ghat[i], 1, Ghat[i-1], -1) for i in range(1, k+1)]
    D2 = [poly_lincomb(D1[i], 1, D1[i-1], -1) for i in range(1, len(D1))]
    d1_nonreal = 0; d1_interlace = set(); d1_tpow = set()
    d1_roots = []
    for p in D1:
        lo, rts, mi = roots_mp(p)
        d1_tpow.add(lo)
        if mi > 1e-20 or (rts and (rts[0] <= -1 or rts[-1] >= 0)):
            d1_nonreal += 1
        d1_roots.append(rts)
    for j in range(len(d1_roots) - 1):
        d1_interlace.add(check_interlace_equal(d1_roots[j], d1_roots[j+1]))

    d2_nonreal = 0; d2_tpow = set(); d2_pos_coeffs = 0; d2_neg_coeffs = 0
    for p in D2:
        lo, rts, mi = roots_mp(p)
        d2_tpow.add(lo)
        if mi > 1e-20 or (rts and (rts[0] <= -1 or rts[-1] >= 0)):
            d2_nonreal += 1
        if all(c >= 0 for c in p): d2_pos_coeffs += 1
        else: d2_neg_coeffs += 1
    print(f"  k={k:2d}: D1 t^pow={d1_tpow}, nonreal={d1_nonreal}/{len(D1)}, D1_interlace={d1_interlace} | "
          f"D2 t^pow={d2_tpow}, nonreal={d2_nonreal}/{len(D2)}, all_pos_coeffs={d2_pos_coeffs}/{len(D2)}", flush=True)

print("\n=== 2. INSPECT EXACT D1 AND D2 COEFFICIENTS AT k = 3, 4 ===")
for k in [3, 4]:
    m = k + 1; kfact = math.factorial(k)
    Ghat = []
    for i in range(k+1):
        p = list(H[k][i])
        for _ in range(i): p = apply_D(m, p)
        scale = kfact // math.comb(k, i)
        Ghat.append([c * scale for c in p])
    D1 = [poly_lincomb(Ghat[i], 1, Ghat[i-1], -1) for i in range(1, k+1)]
    D2 = [poly_lincomb(D1[i], 1, D1[i-1], -1) for i in range(1, len(D1))]
    print(f"  k={k}: D1[0]/k! = {[c/kfact for c in D1[0]]}")
    print(f"  k={k}: D2[0]/k! = {[round(c/kfact, 4) for c in D2[0]]}")

print("\n=== 3. AT WHAT m DOES Q_m = (htilde_{k-1,i-2}, D_m htilde_{k-1,i-1}, D_m^2 htilde_{k-1,i}) FAIL J >= 0? ===")
# Check at t = 0 FIRST (exact c1^2 <= 4 c0 c2) for k = 10, 20, 35, 50!
KBIG = 50
a0 = [[0]*(k+1) for k in range(KBIG+1)]; a1 = [[0]*(k+1) for k in range(KBIG+1)]
a0[1][0], a1[1][0] = 1, 2; a0[1][1], a1[1][1] = 1, 0
for k in range(2, KBIG+1):
    m = k + 1
    for i in range(k+1):
        u0 = (a0[k-1][i] if i<=k-1 else 0) + (a0[k-1][i-1] if i>=1 else 0)
        u1 = ((2*a1[k-1][i] + m*a0[k-1][i]) if i<=k-1 else 0) + (a1[k-1][i-1] if i>=1 else 0)
        a0[k][i] = u0; a1[k][i] = 2*u1 + m*u0

for k in [10, 20, 35, 50]:
    i = k // 2
    x0, x1 = a0[k-1][i-2], a1[k-1][i-2]
    b0m, b1m = a0[k-1][i-1], a1[k-1][i-1]
    b0y, b1y = a0[k-1][i], a1[k-1][i]
    # Test m = k+1, 2k, 10k, 100k, 10^6, 10^10, 10^20
    res = []
    for mult in [1, 2, 10, 100, 10**6, 10**12]:
        mv = (k + 1) * mult
        m0, m1 = b0m, 2*b1m + mv*b0m
        y0, y1 = b0y, 4*b1y + 3*mv*b0y
        c0 = x1*m0 - x0*m1; c1 = x1*y0 - x0*y1; c2 = m1*y0 - m0*y1
        rat = float(mp.mpf(c1*c1)/mp.mpf(c0*c2)) if c0*c2 > 0 else -1.0
        res.append(f"{mult}*(k+1):{rat:.4f}")
    print(f"  k={k:2d}, i={i:2d} (t=0 c1^2/(c0*c2), need<=4): " + ", ".join(res))
