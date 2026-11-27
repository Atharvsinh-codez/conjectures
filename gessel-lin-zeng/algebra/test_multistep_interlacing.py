"""
3-SECOND TEST OF THE MULTI-STEP INTERLACING HIERARCHY I(k, r):
  Does D_{k+2}^r h~_{k, i} << h~_{k, i - r - 1} hold for ALL k, i, r?!
  More generally, for each (k, i, r), what is the MAXIMUM real m = m_max(k, i, r)
  such that D_m^r h~_{k, i} << h~_{k, i - r - 1} (AND also test D_{k+1}^r h~_{k, i} << h~_{k, i-r-1})?!
  ALSO test in the x-variable: S_{L-1} ... S_{L-r} A_{k, i} << A_{k, i - r - 1}!!
  WAIT A SECOND! What IS D_{k+2}^r h~_{k, i} in the x-variable?!
  Remember: D_j on deg d poly corresponds to S_{d + j - 1}!
  For h~_{k, i} (deg d = 2k-i-1):
    D_{k+2} h~_{k, i}  <-> S_{3k - i} A_{k, i}  (deg d+1)
    D_{k+2}^2 h~_{k, i} <-> S_{3k - i + 1} S_{3k - i} A_{k, i} (deg d+2)
    ...
    D_{k+2}^r h~_{k, i} <-> S_{3k - i + r - 1} ... S_{3k - i} A_{k, i}!!
  AND check the exact commutator property: does D_{j+1} D_j f << D_j D_{j+1} f ALWAYS hold for f in RZ(-1,0)?!
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 50

def apply_Dj_mp(p, j):
    out = [mp.mpf(0)] * (len(p) + 1)
    for m, c in enumerate(p):
        if c != 0:
            out[m] += (1 + m) * c
            out[m + 1] += (j + m) * c
    return out

def poly_add_mp(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return out

def roots_t(p):
    # p is poly in t with all roots in (-1, 0), positive coeffs
    d = len(p) - 1
    if d == 0:
        return np.array([])
    # Normalize with mpmath scaling for accurate float64 roots (or use mp.polyroots if d<=20)
    # Wait: roots are in (-1, 0), let's convert to float after normalizing max coeff!
    mx = max(abs(c) for c in p)
    pf = [float(c / mx) for c in p]
    r = np.roots(pf[::-1])
    if np.max(np.abs(r.imag)) > 1e-6:
        # fallback to mpmath to avoid float64 noise!
        try:
            rmp = mp.polyroots(list(reversed([c / mx for c in p])), maxsteps=300, extraprec=30)
            if max(abs(z.imag) for z in rmp) > 1e-12:
                return None
            return np.sort([float(z.real) for z in rmp])
        except Exception:
            return None
    return np.sort(r.real)

def interlaces_deg1(rp, rq):
    # deg q = deg p + 1; check rq[0] < rp[0] < rq[1] < ... < rp[m-1] < rq[m]
    if rp is None or rq is None:
        return False
    m, n = len(rp), len(rq)
    if n != m + 1:
        return False
    return all(rq[q] < rp[q] < rq[q + 1] for q in range(m))

# 1. Test I(k, r): D_m^r h~_{k, i} << h~_{k, i - r - 1} for m = k+1, k+2, k+3, 2k, 100!
print("=== TEST 1: D_m^r h~_{k, i} << h~_{k, i - r - 1} ON GLZ (k = 2..14) ===")
h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
h_store = {1: h}
for k in range(2, 15):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = poly_add_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = poly_add_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
    h_store[k] = h

for offset_name, m_fn in [("m=k+1", lambda k, r: k+1),
                          ("m=k+2", lambda k, r: k+2),
                          ("m=k+2+r", lambda k, r: k+2+r),
                          ("m=2k", lambda k, r: 2*k)]:
    tot = 0
    fails = 0
    first_fail = None
    for k in range(2, 13):
        hk = h_store[k]
        # precompute roots of hk
        rk = [roots_t(p) for p in hk]
        for r in range(1, k):
            m = m_fn(k, r)
            for i in range(r + 1, k + 1):
                tot += 1
                p = hk[i]
                for _ in range(r):
                    p = apply_Dj_mp(p, m)
                rp = roots_t(p)
                if not interlaces_deg1(rp, rk[i - r - 1]):
                    fails += 1
                    if first_fail is None:
                        first_fail = (k, i, r, m)
    print(f"  {offset_name:10s} (k=2..12): fails = {fails}/{tot}, first_fail = {first_fail}")

# 2. Check exact failure boundary in m for k=4..12, r=1..k-1!
print("\n=== TEST 2: Exact failure (if any) of D_{k+2}^r h~_{k,i} << h~_{k,i-r-1} at k=4..12 by r ===")
for r in range(1, 9):
    tot = 0; fails = 0; ff = None
    for k in range(r + 1, 13):
        hk = h_store[k]
        rk = [roots_t(p) for p in hk]
        for i in range(r + 1, k + 1):
            tot += 1
            p = hk[i]
            for _ in range(r):
                p = apply_Dj_mp(p, k + 2)
            if not interlaces_deg1(roots_t(p), rk[i - r - 1]):
                fails += 1
                if ff is None: ff = (k, i)
    print(f"  r={r}: fails = {fails}/{tot}, first_fail={ff}")

# 3. Check commutator interlacing: does D_{j+1} D_j f << D_j D_{j+1} f ALWAYS hold for f in RZ(-1, 0)?!
print("\n=== TEST 3: Universal commutator check D_{j+1} D_j f << D_j D_{j+1} f ===")
rng = np.random.default_rng(42)
for d in [1, 2, 4, 8]:
    fails_c = 0
    for _ in range(300):
        rts = -rng.uniform(0.001, 0.999, d)
        f = [mp.mpf(c) for c in np.poly(rts)[::-1]]  # roots in (-1, 0)
        j = int(rng.integers(2, 20))
        lhs = apply_Dj_mp(apply_Dj_mp(f, j), j + 1)
        rhs = apply_Dj_mp(apply_Dj_mp(f, j + 1), j)
        rl = roots_t(lhs)
        rr = roots_t(rhs)
        # equal degree d+2: check rl[q] < rr[q] < rl[q+1]
        ok = (rl is not None and rr is not None and
              all(rl[q] < rr[q] for q in range(len(rl))) and
              all(rr[q] < rl[q + 1] for q in range(len(rl) - 1)))
        if not ok:
            fails_c += 1
    print(f"  deg={d}: fails(D_{{j+1}} D_j f << D_j D_{{j+1}} f) = {fails_c}/300")
