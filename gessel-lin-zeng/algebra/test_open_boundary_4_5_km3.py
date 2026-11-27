"""
Fast 5-second test of:
1. Open boundary cases i = 4, i = 5, i = k - 3 (ALL below the i>=6, k-i>=4 barrier!!):
   Remember for fixed (k, i), by proved root monotonicity (SOURCES_AND_STATE.md Sec 4):
     A_{k, {k-i+1..k}} (late) <= A_{k, S} <= A_{k, {1..i}} (early)  for ALL |S| = i!!
   So a polynomial M is a common upper interlacer for ALL {A_{k,S} : |S|=i} iff:
     A_{k, {k-i+1..k}} << M   AND   A_{k, {1..i}} << M  (actually just the two extremes!!)!
   Check through k = 25 (in 60-digit mpmath!):
   - i = 4: Does A_{k, S} << H_{k-1} := S_{3k-6} A_{k-2, 0} hold for ALL |S|=4?
            Does S_{3k-7} A_{k-1, 4} << A_{k-2, 0} (Ma-Wang Remark 1!) hold?
            Does S_{3k-7} A_{k-1, {1..4}} << A_{k-2, 0} hold?
   - i = 5: Does A_{k, S} << A_{k-2, 0} hold for ALL |S|=5?
   - i = k - 3: Does A_{k+2, k+2} << A_{k, S} hold for ALL |S|=k-3?
            (And check i = k - 4: A_{k+3, k+3} << A_{k, S} — should fail since k-i=4 is the barrier!)
2. Diagonal interlacings for ALL (k, i) (k <= 16):
   - A_{k, i} << A_{k-1, i-3} (for 3 <= i <= k)
   - A_{k, i} << A_{k-1, i-2} (equal deg: which direction?)
   - A_{k+1, i+3} << A_{k, i}
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 60

def apply_SL_mp(p, L):
    out = [mp.mpf(0)] * (len(p) + 1)
    for m, c in enumerate(p):
        if c != 0:
            out[m] += (1 + m) * c
            out[m + 1] += (L - m) * c
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out

def padd_mp(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for idx in range(len(a)): out[idx] += a[idx]
    for idx in range(len(b)): out[idx] += b[idx]
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out

def apply_Dj_mp(p, j):
    out = [mp.mpf(0)] * (len(p) + 1)
    for m, c in enumerate(p):
        if c != 0:
            out[m] += (1 + m) * c
            out[m + 1] += (j + m) * c
    return out

def build_h_S(k, S):
    """Build h~_{k, S}(t) = prod_{j=k..2} D_{j+1}^{2 - 1_{j in S}} D_2^{1 - 1_{1 in S}} (1)"""
    p = [mp.mpf(1)]
    if 1 not in S:
        p = apply_Dj_mp(p, 2)
    for j in range(2, k + 1):
        p = apply_Dj_mp(p, j + 1)
        if j not in S:
            p = apply_Dj_mp(p, j + 1)
    return p

def roots_t_mp(p):
    # Compute roots in (-1, 0) via exact Sturm/bisection or mp.polyroots on normalized poly
    d = len(p) - 1
    if d == 0: return []
    mx = max(abs(c) for c in p)
    # Try float64 first, verify separation, fallback to mp if close/complex
    pf = [float(c / mx) for c in p]
    r = np.roots(pf[::-1])
    if np.max(np.abs(r.imag)) < 1e-6:
        rr = np.sort(r.real)
        if len(rr) == d and np.all(rr > -1.0) and np.all(rr < 0.0):
            return list(rr)
    rmp = mp.polyroots(list(reversed([c / mx for c in p])), maxsteps=400, extraprec=40)
    if max(abs(z.imag) for z in rmp) > 1e-12:
        return None
    return sorted([float(z.real) for z in rmp])

def check_interlace_t(p, q):
    """Check p << q in (-1, 0): deg q in {deg p, deg p+1}, rightmost root of q > rightmost root of p."""
    # Use Cauchy index / sign check at roots OR Wronskian positivity?
    # Even simpler and exact: evaluate p at roots of q (if deg q = deg p + 1)!
    # Wait: roots_t_mp is fast for k <= 16 (deg <= 30), let's check both!
    rp = roots_t_mp(p)
    rq = roots_t_mp(q)
    if rp is None or rq is None:
        return False, "non-real"
    m, n = len(rp), len(rq)
    if n == m + 1:
        ok = all(rq[idx] < rp[idx] < rq[idx + 1] for idx in range(m))
        if not ok:
            # compute min margin
            margins = [min(rp[idx] - rq[idx], rq[idx+1] - rp[idx]) for idx in range(m)]
            return False, f"min_margin={min(margins):.3e}"
        return True, f"min_margin={min(min(rp[idx]-rq[idx], rq[idx+1]-rp[idx]) for idx in range(m)):.3e}"
    elif n == m:
        ok = all(rp[idx] < rq[idx] for idx in range(m)) and all(rq[idx] < rp[idx+1] for idx in range(m-1))
        return ok, "eq-deg"
    return False, f"deg_mismatch({m},{n})"

# Build h~_{k, 0} and h~_{k, k} up to k = 18
h0 = {1: [mp.mpf(1), mp.mpf(2)]}
hk_k = {1: [mp.mpf(1)]}
for k in range(2, 19):
    h0[k] = apply_Dj_mp(apply_Dj_mp(h0[k-1], k+1), k+1)
    hk_k[k] = apply_Dj_mp(hk_k[k-1], k+1)

print("=== 1A. OPEN CASE i = 4: Does A_{k, S} << H_{k-1} (i.e. h~_{k, S} << D_k h~_{k-2, 0}) hold for ALL |S|=4? ===")
# Wait: check t-transform of H_{k-1} = S_{3k-6} A_{k-2, 0}!
# deg A_{k-2, 0} = 2k - 5. S_L on deg d -> D_{L + 1 - d}: (3k-6) + 1 - (2k-5) = k!!
# YES!! H_{k-1} <-> D_k h~_{k-2, 0} (deg 2k - 4), while h~_{k, S} (|S|=4) has deg 2k - 5!
for k in [5, 6, 8, 10, 12, 15]:
    early = build_h_S(k, {1, 2, 3, 4})
    late  = build_h_S(k, {k-3, k-2, k-1, k})
    H_km1 = apply_Dj_mp(h0[k-2], k)
    ok_e, me = check_interlace_t(early, H_km1)
    ok_l, ml = check_interlace_t(late, H_km1)
    # Also check Ma-Wang Remark 1: S_{3k-7} A_{k-1, S} << A_{k-2, 0} (i.e. D_{k+1} h~_{k-1, S} << h~_{k-2, 0})!
    # Wait: deg(D_{k+1} h~_{k-1, S}) = (2k-7)+1 = 2k-6, deg(h~_{k-2, 0}) = 2k-5!
    e_mw = apply_Dj_mp(build_h_S(k-1, {1, 2, 3, 4}), k+1)
    l_mw = apply_Dj_mp(build_h_S(k-1, {k-4, k-3, k-2, k-1}), k+1) if k >= 5 else None
    ok_mw_e, _ = check_interlace_t(e_mw, h0[k-2])
    ok_mw_l, _ = check_interlace_t(l_mw, h0[k-2])
    print(f"  k={k:2d}: early<<H_{{k-1}}:{ok_e} ({me}), late<<H_{{k-1}}:{ok_l} ({ml}) | MW-Rem1 (D_{{k+1}} h_{{k-1,S}} << h_{{k-2,0}}): early={ok_mw_e}, late={ok_mw_l}")

print("\n=== 1B. OPEN CASE i = 5: Does A_{k, S} << A_{k-2, 0} (h~_{k, S} << h~_{k-2, 0}) hold for ALL |S|=5? ===")
for k in [6, 7, 8, 10, 12, 15]:
    early = build_h_S(k, {1, 2, 3, 4, 5})
    late  = build_h_S(k, {k-4, k-3, k-2, k-1, k})
    ok_e, me = check_interlace_t(early, h0[k-2])
    ok_l, ml = check_interlace_t(late, h0[k-2])
    print(f"  k={k:2d}: early<<h_{{k-2,0}}:{ok_e} ({me}), late<<h_{{k-2,0}}:{ok_l} ({ml})")

print("\n=== 1C. OPEN CASE i = k - 3: Does A_{k+2, k+2} << A_{k, S} (h~_{k+2, k+2} << h~_{k, S}) hold for ALL |S|=k-3? ===")
# Check degrees: deg h~_{k, k-3} = 2k - (k-3) - 1 = k + 2; deg h~_{k+2, k+2} = k + 1!
for k in [4, 5, 6, 8, 10, 12, 15]:
    i = k - 3
    early = build_h_S(k, set(range(1, i + 1)))
    late  = build_h_S(k, set(range(k - i + 1, k + 1)))
    ok_e, me = check_interlace_t(hk_k[k + 2], early)
    ok_l, ml = check_interlace_t(hk_k[k + 2], late)
    print(f"  k={k:2d} (i={i:2d}): h_{{k+2,k+2}}<<early:{ok_e} ({me}), h_{{k+2,k+2}}<<late:{ok_l} ({ml})")

print("\n=== 1D. CHECK i = k - 4 (SHOULD FAIL AT LARGE k): h~_{k+3, k+3} << h~_{k, S} ===")
for k in [5, 6, 8, 10, 12, 15]:
    i = k - 4
    early = build_h_S(k, set(range(1, i + 1)))
    late  = build_h_S(k, set(range(k - i + 1, k + 1)))
    ok_e, _ = check_interlace_t(hk_k[k + 3], early)
    ok_l, _ = check_interlace_t(hk_k[k + 3], late)
    print(f"  k={k:2d} (i={i:2d}): h_{{k+3,k+3}}<<early:{ok_e}, h_{{k+3,k+3}}<<late:{ok_l}")
