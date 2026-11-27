"""
Fast, exact/high-precision surgical checks (< 15s total):
1. GLZ D-lifted Wronskian span W(E_{k,i}, E_{k,i+r}) > 0 for r=1,2,3,4 in 60-digit mpmath (k=4..18).
2. Route II Q_{k,q}:
   - Exact root ratio min(r_m / r_{m+1}) of f_k(y) = sum C(k,i) q^{C(i,2)} y^i vs 1/q!
   - Test adversarial b_j >= 2 (with w_j = 1, AND with wild w_j > 0!) at q = 1/4, 1/8, 1/16!
   - Also check the EXACT deformed exponential threshold for sum (q^{C(i,2)} / i!) y^i.
3. Route I Strengthened Cubic on GLZ:
   - Exact c_min(k) for GLZ E_{k,i} (k=4..14) AND exact Stage 3 drop of c_min on GLZ!
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 60

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

def poly_eval_and_deriv_mp(p, t):
    # Horner for p(t) and p'(t) simultaneously!
    val = mp.mpf(0)
    dval = mp.mpf(0)
    for c in reversed(p):
        dval = dval * t + val
        val = val * t + c
    return val, dval

# ============================================================
# PART 1: GLZ W(E_{k,i}, E_{k,i+r}) > 0 IN 60-DIGIT MPMATH
# ============================================================
print("=== PART 1: GLZ W(E_{k,i}, E_{k,i+r}) > 0 (60-digit mpmath) ===")
h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
ts_mp = [mp.mpf(x) for x in [-2.0, -1.1, -0.999, -0.95, -0.8, -0.5, -0.2, -0.05, -0.001, -1e-6, 0.0, 0.1, 1.0, 10.0]]

h_store = {}
for k in range(2, 17):
    j = k + 1
    E = []
    for i in range(k):
        p = h[i]
        for _ in range(i):
            p = apply_Dj_mp(p, j)
        E.append(p)
    if k in [4, 6, 8, 10, 12, 14, 16]:
        # evaluate (E_i(t), E'_i(t)) at all grid points
        vals = [[poly_eval_and_deriv_mp(E[i], t) for t in ts_mp] for i in range(k)]
        msg = []
        for r in [1, 2, 3, 4]:
            if r >= k: continue
            neg = 0
            min_norm = mp.mpf('1e100')
            for i in range(k - r):
                for ti in range(len(ts_mp)):
                    a0, a1 = vals[i][ti]
                    b0, b1 = vals[i + r][ti]
                    w = a0 * b1 - a1 * b0
                    c0, c1 = vals[i + 1][ti]
                    w1 = a0 * c1 - a1 * c0
                    rat = w / w1
                    if rat < min_norm: min_norm = rat
                    if w <= 0: neg += 1
            msg.append(f"r={r}: neg={neg} (min={float(min_norm):+.3e})")
        print(f"  k={k:2d}: " + ", ".join(msg))
    # step
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = poly_add_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = poly_add_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
    h_store[k] = h

# ============================================================
# PART 2: ROUTE II Q_{k,q} - ROOT RATIO OF f_k(y) & WILD w_j TEST
# ============================================================
print("\n=== PART 2A: Root ratio min(r_m/r_{m+1}) of f_k(y) = sum C(k,i) q^{C(i,2)} y^i vs 1/q ===")
for q_val in [0.5, 0.35, 0.30, 0.25, 0.20, 0.125]:
    qm = mp.mpf(str(q_val))
    # Check k = 10, 20, 40
    ratios = []
    for k in [8, 16, 30]:
        coeffs = [mp.binomial(k, i) * (qm ** (i*(i-1)//2)) for i in range(k + 1)]
        # rescale y -> s*u
        s = (coeffs[0] / coeffs[-1]) ** (mp.mpf(1) / k)
        sc = [coeffs[i] * (s ** i) for i in range(k + 1)]
        mx = max(sc)
        sc = [x / mx for x in sc]
        rts = mp.polyroots(list(reversed(sc)), maxsteps=400, extraprec=40)
        max_im = max(abs(r.imag) for r in rts)
        if max_im > 1e-20:
            ratios.append(f"k={k}:COMPLEX")
        else:
            reals = sorted([r.real for r in rts]) # most negative to least negative
            # ratio r_m / r_{m+1} > 1 (since r_m < r_{m+1} < 0)
            min_rat = min(reals[m] / reals[m+1] for m in range(k - 1))
            ratios.append(f"k={k}:{float(min_rat):.4f}")
    print(f"  q={q_val:.3f} (need > 1/q = {1/q_val:.3f} for b_k->inf!): " + ", ".join(ratios))

print("\n=== PART 2B: Wild w_j test on f_k(y) = prod (1 + w_j Y_q)(1) at q=1/4, 1/8 ===")
# At t=0, if b_k -> infty, we NEED f_{k-1}(y) << y f_{k-1}(q y), i.e. root separation > 1/q!
# What if w_1 = 1, w_2 = 1e-6, w_3 = 1, ...?!
for q_val in [0.25, 0.125, 0.05]:
    # Test f_2 = 1 + (w1+w2)y + q w1 w2 y^2 -> roots of (1 + w1 y) + w2 y(1 + q w1 y)
    # Wait: if w1 = 1, w2 = 1, w3 = 1e-6, ... what is min root ratio?
    # Even simpler: f_2 with w1=1, w2=1 has roots ratio (1+sqrt(1-q))^2 / q > 1/q!
    # What about w1=1, w2=w (arbitrary w > 0)?
    # f_2 = 1 + (1+w)y + q w y^2. Disc = (1+w)^2 - 4 q w >= (1+w)^2 - (1+w)^2 * q = (1-q)(1+w)^2.
    # Ratio of roots = ((1+w) + sqrt((1+w)^2 - 4qw)) / ((1+w) - sqrt((1+w)^2 - 4qw))
    #                = ((1+w) + sqrt(...))^2 / (4 q w) >= (1 + sqrt(1-q))^2 / q > 1/q (MIN AT w=1!!)!
    # WHOA! For k=2, arbitrary w1, w2 > 0 has MINIMUM root ratio at w1 = w2!
    # Check k = 3, 4, 5 with random wild w_j in [1e-4, 1e4]!
    rng = np.random.default_rng(9)
    min_rat_seen = 1e9
    fails_interlace = 0
    for _ in range(300):
        ws = np.exp(rng.normal(0, 3.0, 6))
        p = np.array([1.0])
        for w in ws:
            # p(y) + w y p(q y)
            p_qy = np.array([p[i] * (q_val ** i) for i in range(len(p))])
            p_new = np.zeros(len(p) + 1)
            p_new[:len(p)] += p
            p_new[1:] += w * p_qy
            p = p_new
        rts = np.roots(p[::-1])
        if np.max(np.abs(rts.imag)) > 1e-6:
            fails_interlace += 1
            continue
        reals = np.sort(rts.real)
        mr = np.min(reals[:-1] / reals[1:])
        if mr < min_rat_seen:
            min_rat_seen = mr
        if mr < 1.0 / q_val:
            fails_interlace += 1
    print(f"  Wild w_j (k=6), q={q_val:.3f} (1/q={1/q_val:.2f}): min root ratio = {min_rat_seen:.4f}, fails(<1/q) = {fails_interlace}/300")

# ============================================================
# PART 3: STRENGTHENED CUBIC c_min(k) & STAGE 3 DROP ON GLZ
# ============================================================
print("\n=== PART 3: GLZ STRENGTHENED CUBIC c_min(k) AND STAGE 3 DROP ===")
def bezoutiant_cubic_mp_float(v0, v1, v2, v3, c):
    # compute normalized W_{rs} in mpmath, convert to float for eigvalsh
    vs = [v0, v1, v2, v3]
    def W(i, j): return vs[i][0]*vs[j][1] - vs[i][1]*vs[j][0]
    w12 = W(1, 2)
    b = [1.0, -3.0*c, 3.0*c, -1.0]
    B = np.zeros((3, 3), dtype=float)
    for r in range(4):
        for s in range(r+1, 4):
            wrs = float(W(r, s) / w12)
            val = -b[r]*b[s]*wrs
            for m in range(s - r):
                B[r + m, r + (s - r - 1 - m)] += val
    return B

def find_cmin_glz_quad(v0, v1, v2, v3):
    lo, hi = 0.01, 1.0
    if np.linalg.eigvalsh(bezoutiant_cubic_mp_float(v0, v1, v2, v3, hi))[0] < 0:
        return 1.5
    for _ in range(45):
        mid = 0.5*(lo + hi)
        if np.linalg.eigvalsh(bezoutiant_cubic_mp_float(v0, v1, v2, v3, mid))[0] >= 0:
            hi = mid
        else:
            lo = mid
    return hi

ts_c = [mp.mpf(x) for x in [-2.0, -1.05, -0.99, -0.8, -0.5, -0.2, -0.02, 0.0, 0.2, 2.0, 20.0]]
for k in [5, 7, 9, 11, 13, 15]:
    # Compare E^{(j)} = (h_{k-1,i}, D_j h_{k-1,i+1}, D_j^2 h_{k-1,i+2}, D_j^3 h_{k-1,i+3})
    # at j = k+1 vs j = k+2 (Stage 3 shift!)
    hk1 = h_store[k - 1]
    max_c_j = 0.0
    max_c_j1 = 0.0
    max_drop = -1e9
    for i in range(k - 3):
        # Note: in our E_{k,i} = D_j^i h_{k-1,i}, shifting by D_j^i gives (h_{i}, D_j h_{i+1}, D_j^2 h_{i+2}, D_j^3 h_{i+3})
        # Wait: let's check both j = k+1 (before Stage 3) and j = k+2 (after Stage 3) on h_{k, i}!!
        # Actually Stage 3 is shifting h_k from D_{k+1} to D_{k+2}!
        hk = h_store[k]
        p0 = hk[i]
        p1_a = apply_Dj_mp(hk[i+1], k+1); p1_b = apply_Dj_mp(hk[i+1], k+2)
        p2_a = apply_Dj_mp(apply_Dj_mp(hk[i+2], k+1), k+1); p2_b = apply_Dj_mp(apply_Dj_mp(hk[i+2], k+2), k+2)
        p3_a = apply_Dj_mp(p2_a, k+1); p3_b = apply_Dj_mp(p2_b, k+2)
        # Wait: p3_a should be D_{k+1}^3 hk[i+3]! Fix!
        p3_a = apply_Dj_mp(apply_Dj_mp(apply_Dj_mp(hk[i+3], k+1), k+1), k+1)
        p3_b = apply_Dj_mp(apply_Dj_mp(apply_Dj_mp(hk[i+3], k+2), k+2), k+2)
        for t in ts_c:
            va = [poly_eval_and_deriv_mp(p, t) for p in (p0, p1_a, p2_a, p3_a)]
            vb = [poly_eval_and_deriv_mp(p, t) for p in (p0, p1_b, p2_b, p3_b)]
            ca = find_cmin_glz_quad(*va)
            cb = find_cmin_glz_quad(*vb)
            if ca > max_c_j: max_c_j = ca
            if cb > max_c_j1: max_c_j1 = cb
            if cb - ca > max_drop: max_drop = cb - ca
    print(f"  k={k:2d}: max c_min(j=k+1)={max_c_j:.6f} (k*(1-c)={k*(1-max_c_j):.3f}), max c_min(j=k+2)={max_c_j1:.6f}, max Stage3 drop={max_drop:+.6f} (k^2*drop={k*k*max_drop:+.3f})")
