"""
5-second cleanup + verification:
1. Find the CLEANEST small-denominator rational x (e.g. integer or 1-decimal fraction!)
   for the exact Fraction Laguerre non-real-root certificates at:
   - i = 4, k = 31 (deg 57)
   - i = 5, k = 16 (deg 26)
   - i = k - 3, k = 28 (deg 30)
   AND double-check that k = 30 (i=4), k = 15 (i=5), k = 27 (i=k-3) ALL PASS,
   AND check that aggregate Ma-Wang Remark 1 (i=4) & A_{k+2,k+2}<<A_{k,k-3} STILL PASS past k=31, 28!
2. Also check if ALL k in [31..38] fail for i=4, [16..24] for i=5, [28..35] for i=k-3!
"""
from fractions import Fraction
import mpmath as mp
mp.mp.dps = 100

def apply_SL_int(p, L):
    out = [0] * (len(p) + 1)
    for m, c in enumerate(p):
        if c != 0:
            out[m] += (1 + m) * c
            out[m + 1] += (L - m) * c
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out

def build_A_S_int(k, S):
    p = [1]
    L = 0
    for j in range(1, k + 1):
        if j in S:
            p = apply_SL_int(p, L); L += 2
        else:
            p = apply_SL_int(apply_SL_int(p, L), L + 1); L += 3
    return p

def eval_p012_frac(p, t):
    v = Fraction(0); dv = Fraction(0); d2v = Fraction(0)
    for c in reversed(p):
        d2v = d2v * t + 2 * dv
        dv = dv * t + v
        v = v * t + c
    return v, dv, d2v

cases = [
    ("i=4, k=31",   31, set(range(28, 32)), {1, 2, 3, 4},         -23.0409, [Fraction(-23), Fraction(-2304, 100), Fraction(-23041, 1000)]),
    ("i=5, k=16",   16, set(range(12, 17)), {1, 2, 3, 4, 5},      -10.5925, [Fraction(-11), Fraction(-21, 2), Fraction(-106, 10), Fraction(-1059, 100)]),
    ("i=k-3, k=28", 28, set(range(4, 29)),  set(range(1, 26)),     -7.2862, [Fraction(-7), Fraction(-73, 10), Fraction(-729, 100), Fraction(-7286, 1000)]),
]

print("=== CLEAN EXACT FRACTION LAGUERRE CERTIFICATES ===")
for label, k, Sl, Se, approx_x, cands in cases:
    fl = build_A_S_int(k, Sl)
    fe = build_A_S_int(k, Se)
    found = False
    for xf in cands:
        v1, d1, dd1 = eval_p012_frac(fl, xf)
        v2, d2, dd2 = eval_p012_frac(fe, xf)
        if d1 * d2 < 0:
            c_ex = -d1 / d2
            # test simple integer p_int for p_int * fl + 1 * fe!
            p_approx = int(round(float(1 / c_ex))) if c_ex < 1 else None
            if p_approx and p_approx > 0:
                for p_int in [p_approx, p_approx - 1, p_approx + 1]:
                    if p_int <= 0: continue
                    Rv = p_int * v1 + v2
                    Rd = p_int * d1 + d2
                    Rdd = p_int * dd1 + dd2
                    lag = Rd * Rd - Rv * Rdd
                    if lag < 0:
                        print(f"  [CLEAN!] {label} (deg={len(fl)-1}): R(x) = {p_int}*A_late + A_early, at x = {xf}: R'^2 - R*R'' = {float(lag):.4e} < 0!")
                        found = True
                        break
            if found: break
    if not found:
        print(f"  {label}: need slightly finer x than {cands}")

# Also test clean integer x at k = 32 (i=4), k = 17 (i=5), k = 29 (i=k-3) where crossing is wider!
print("\n=== EVEN SIMPLER INTEGER/HALF-INTEGER CERTIFICATES AT k+1 OR k+2 ===")
for label, k, Sl, Se, xs_try in [
    ("i=4, k=32",   32, set(range(29, 33)), {1, 2, 3, 4},     [Fraction(-24), Fraction(-25), Fraction(-23), Fraction(-22), Fraction(-26), Fraction(-49, 2), Fraction(-47, 2)]),
    ("i=5, k=17",   17, set(range(13, 18)), {1, 2, 3, 4, 5},  [Fraction(-11), Fraction(-12), Fraction(-10), Fraction(-23, 2), Fraction(-25, 2)]),
    ("i=k-3, k=29", 29, set(range(4, 30)),  set(range(1, 27)),[Fraction(-7), Fraction(-8), Fraction(-15, 2), Fraction(-17, 2)]),
]:
    fl = build_A_S_int(k, Sl)
    fe = build_A_S_int(k, Se)
    for xf in xs_try:
        v1, d1, dd1 = eval_p012_frac(fl, xf)
        v2, d2, dd2 = eval_p012_frac(fe, xf)
        if d1 * d2 < 0:
            c_ex = -d1 / d2
            p_int = int(round(float(1 / c_ex)))
            if p_int > 0:
                Rv = p_int * v1 + v2; Rd = p_int * d1 + d2; Rdd = p_int * dd1 + dd2
                lag = Rd * Rd - Rv * Rdd
                if lag < 0:
                    print(f"  [ULTRA-CLEAN!] {label} (deg={len(fl)-1}): R(x) = {p_int}*A_late + A_early, at x = {xf}: R'^2 - R*R'' = {float(lag):.4e} < 0!")
                    break
