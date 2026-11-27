"""FAST high-precision + exact test of the 4-pair interlacing closure:
Check:
  (A) p_{k,i} << p_{k,i-1}   (outer-interlaced, deg diff = 1)
  (B) Dtilde_{k+1} p_{k,i} << p_{k,i-2}  (Pair (2,1)! outer-interlaced, deg diff = 1!)
  (C) Dtilde_{k+1} p_{k,i} vs p_{k,i-1}  (same degree! what direction?)
  (D) What about general Dtilde_{k+m} p_{k,i} << p_{k,i-2}?
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 50

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for idx, v in enumerate(a):
        res[idx] += v
    for idx, v in enumerate(b):
        res[idx] += v
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def apply_Dtilde(j, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(poly):
        res[m] += c * (m + 1)
        res[m + 1] += c * (m + j + 1)
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def build_p(kmax):
    p = [[None] * (k + 1) for k in range(kmax + 1)]
    p[1][0] = [1, 2]
    p[1][1] = [1]
    for k in range(2, kmax + 1):
        p[k][0] = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][0]))
        for i in range(1, k):
            p[k][i] = apply_Dtilde(k, poly_add(apply_Dtilde(k, p[k - 1][i]), p[k - 1][i - 1]))
        p[k][k] = apply_Dtilde(k, p[k - 1][k - 1])
    return p

def get_roots(poly):
    if len(poly) <= 1:
        return []
    # mpmath polyroots expects highest-degree first
    coeffs = [mp.mpf(c) for c in reversed(poly)]
    r = mp.polyroots(coeffs, maxsteps=200, extraprec=30)
    assert max(abs(mp.im(z)) for z in r) < mp.mpf('1e-25'), f"Nonreal root! {r}"
    return sorted(mp.re(z) for z in r)

def check_outer(r_small, r_big):
    """Check if r_big (len d+1) strictly outer-interlaces r_small (len d):
       r_big[0] < r_small[0] < r_big[1] < ... < r_small[d-1] < r_big[d].
    """
    d = len(r_small)
    assert len(r_big) == d + 1
    if d == 0:
        return True, mp.mpf(1)
    gaps = []
    for m in range(d):
        g1 = r_small[m] - r_big[m]
        g2 = r_big[m + 1] - r_small[m]
        gaps.append(min(g1, g2))
    mg = min(gaps)
    return (mg > 0), mg

def check_same(rA, rB):
    """Check interlacing for equal-length sorted lists."""
    d = len(rA)
    assert len(rB) == d
    if d == 0:
        return True, "trivial", mp.mpf(1)
    ok_AB = all(rA[m] < rB[m] for m in range(d)) and all(rB[m] < rA[m+1] for m in range(d-1))
    ok_BA = all(rB[m] < rA[m] for m in range(d)) and all(rA[m] < rB[m+1] for m in range(d-1))
    if ok_AB:
        gaps = [rB[m] - rA[m] for m in range(d)] + [rA[m+1] - rB[m] for m in range(d-1)]
        return True, "A<B", min(gaps)
    if ok_BA:
        gaps = [rA[m] - rB[m] for m in range(d)] + [rB[m+1] - rA[m] for m in range(d-1)]
        return True, "B<A", min(gaps)
    # Report how much it fails
    return False, "FAIL", mp.mpf(0)

KMAX = 14
p = build_p(KMAX)
roots = [[None] * (k + 1) for k in range(KMAX + 1)]
for k in range(1, KMAX + 1):
    for i in range(k + 1):
        roots[k][i] = get_roots(p[k][i])

print(f"All p_{{k,i}} real-rooted in (-1,0) for k=1..{KMAX} (50-digit mpmath)!")

# Check (A): p_{k,i} << p_{k,i-1}
fails_A = []
min_gap_A = (mp.mpf(1), None)
for k in range(2, KMAX + 1):
    for i in range(1, k + 1):
        ok, g = check_outer(roots[k][i], roots[k][i-1])
        if not ok:
            fails_A.append((k, i, float(g)))
        elif g < min_gap_A[0]:
            min_gap_A = (g, (k, i))
print(f"(A) p_{{k,i}} << p_{{k,i-1}}: fails={fails_A[:10]}, min_gap={float(min_gap_A[0]):.3e} at {min_gap_A[1]}")

# Check (B): Dtilde_{k+1} p_{k,i} << p_{k,i-2}  (THE MAGIC PAIR (2,1)!)
fails_B = []
min_gap_B = (mp.mpf(1), None)
for k in range(2, KMAX + 1):
    for i in range(2, k + 1):
        r_small = get_roots(apply_Dtilde(k + 1, p[k][i]))
        r_big = roots[k][i - 2]
        ok, g = check_outer(r_small, r_big)
        if not ok:
            fails_B.append((k, i, float(g)))
        elif g < min_gap_B[0]:
            min_gap_B = (g, (k, i))
print(f"(B) Dtilde_{{k+1}} p_{{k,i}} << p_{{k,i-2}}: fails={fails_B[:10]}, min_gap={float(min_gap_B[0]):.3e} at {min_gap_B[1]}")

# Check (C): Dtilde_{k+1} p_{k,i} vs p_{k,i-1} (same deg!)
dirs_C = set(); fails_C = []
min_gap_C = (mp.mpf(1), None)
for k in range(1, KMAX + 1):
    for i in range(1, k + 1):
        rA = get_roots(apply_Dtilde(k + 1, p[k][i]))
        rB = roots[k][i - 1]
        ok, d, g = check_same(rA, rB)
        dirs_C.add(d)
        if not ok:
            fails_C.append((k, i, d))
        elif d != "trivial" and g < min_gap_C[0]:
            min_gap_C = (g, (k, i, d))
print(f"(C) Dtilde_{{k+1}} p_{{k,i}} vs p_{{k,i-1}}: dirs={dirs_C}, fails={fails_C[:10]}, min_gap={float(min_gap_C[0]):.3e} at {min_gap_C[1]}")

# Check (D): Dtilde_{k+1}^2 p_{k,i} << p_{k,i-3} ?? And what about Dtilde_k p_{k,i} vs p_{k,i-1}?
dirs_D = set(); fails_D = []
for k in range(1, KMAX + 1):
    for i in range(1, k + 1):
        rA = get_roots(apply_Dtilde(k, p[k][i]))
        rB = roots[k][i - 1]
        ok, d, g = check_same(rA, rB)
        dirs_D.add(d)
        if not ok:
            fails_D.append((k, i, d))
print(f"(D) Dtilde_k p_{{k,i}} vs p_{{k,i-1}}: dirs={dirs_D}, fails={fails_D[:10]}")

# Check (E): Does Dtilde_{k+2} p_{k,i} << p_{k,i-2}? Or Dtilde_{m} p_{k,i} << p_{k,i-2} for m >= k+1?
fails_E = []
for k in range(2, min(KMAX + 1, 11)):
    for i in range(2, k + 1):
        for m in [k, k+1, k+2, k+3]:
            r_small = get_roots(apply_Dtilde(m, p[k][i]))
            r_big = roots[k][i - 2]
            ok, g = check_outer(r_small, r_big)
            if not ok:
                fails_E.append((k, i, m, float(g)))
print(f"(E) Dtilde_m p_{{k,i}} << p_{{k,i-2}} (m in [k..k+3]): fails={fails_E[:15]}")
