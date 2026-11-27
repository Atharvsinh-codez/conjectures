"""
Two fast, decisive tests (< 15s total):
1. CAN WE BREAK real-rootedness of prod_{j=1}^k (D_{b_j} + w_j y)(1) with ALL a_j = 1, ALL b_j >= 2?
   (Test BOTH arbitrary w_j > 0 AND fixed w_j = 1!)
   Also test Q_{k, q} at q = 1/4, 1/8, 1/16 via Hermite-Biehler matrix on extreme b_j >= 2, w_j > 0!
2. Finish Part 3 of test_surgical_3.py: GLZ Strengthened Cubic c_min(k) & Stage 3 drop (k = 5..15),
   and also check how exact Cubic Stage 1 bound behaves as c -> 1^-!
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
    v = mp.mpf(0); dv = mp.mpf(0)
    for c in reversed(p):
        dv = dv * t + v
        v = v * t + c
    return v, dv

def poly_eval_c_mp(p, t):
    s = mp.mpc(0)
    for c in reversed(p):
        s = s * t + c
    return s

def check_hb_pd(coeffs_y):
    n = len(coeffs_y) - 1
    s = (abs(coeffs_y[0]) / abs(coeffs_y[-1])) ** (mp.mpf(1) / n)
    sc = [coeffs_y[j] * (s ** j) for j in range(n + 1)]
    mx = max(abs(z) for z in sc)
    sc = [z / mx for z in sc]
    A = [z.real for z in sc]
    B = [z.imag for z in sc]
    H = [[mp.mpf(0)] * n for _ in range(n)]
    for r in range(n + 1):
        for col in range(r + 1, n + 1):
            val = -(A[r]*B[col] - B[r]*A[col])
            for m in range(col - r):
                H[r + m][r + (col - r - 1 - m)] += val
    for i in range(n):
        piv = H[i][i]
        if piv <= 0:
            return False, float(piv)
        for j in range(i + 1, n):
            fac = H[j][i] / piv
            for k in range(j, n):
                H[j][k] -= fac * H[i][k]
                H[k][j] = H[j][k]
    return True, 1.0

# ============================================================
# TEST 1: TRY TO BREAK a_j = 1, b_j >= 2 (like our k=9 attack!)
# ============================================================
print("=== TEST 1A: Try to drive J_i(t) < 1 (or < 1/4) with a_j = 1, b_j >= 2 ===")
# Remember in k=9, how did we drive J_i(t) < 1?
# Let's test extreme b_j >= 2 sequences (with w_j = 1 and with wild w_j) for k = 6..12,
# measuring min_{i, t} 4 W(h_{i+1}, h_i) W(h_i, h_{i-1}) / W(h_{i+1}, h_{i-1})^2 !
seqs = [
    ("step_2_to_1e6_k8", [2,2,2,2, 1e6,1e6,1e6,1e6], [1]*8),
    ("step_1e6_to_2_k8", [1e6,1e6,1e6,1e6, 2,2,2,2], [1]*8),
    ("alt_2_1e6_k8",     [2,1e6,2,1e6,2,1e6,2,1e6], [1]*8),
    ("geom_10_k8",       [2, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7], [1]*8),
    ("rev_geom_10_k8",   [1e7, 1e6, 1e5, 1e4, 1e3, 100, 10, 2], [1]*8),
    ("wild_w_1",         [2,2,2, 1e6,1e6,1e6, 2, 1e6], [1,1,1, 1e-8, 1e-8, 1e-8, 1e8, 1e-8]),
    ("wild_w_2",         [1e6, 2, 2, 2, 2, 1e6, 1e6, 2], [1, 1e-6, 1e-6, 1e6, 1e-6, 1e6, 1e-6, 1e6]),
    ("k9_style_a1",      [2, 2, 2, 1e4, 1e4, 1e4, 1e8, 1e8, 1e8], [1, 1, 1, 1e-4, 1e-4, 1e-4, 1e-8, 1e-8, 1e-8]),
    ("k9_style_w1",      [2, 2, 2, 1e4, 1e4, 1e4, 1e8, 1e8, 1e8], [1]*9),
]
ts_probe = [mp.mpf(x) for x in [-0.999, -0.9, -0.5, -0.1, -1e-3, -1e-5, -1e-7, -1e-9, -1e-11, 0.0, 1.0, 100.0]]
ts_uhp = [mp.mpc('-1e-12', '1e-15'), mp.mpc('-1e-8', '1e-11'), mp.mpc('-1e-4', '1e-7'),
          mp.mpc('-0.5', '1e-4'), mp.mpc('-0.99', '1e-5'), mp.mpc('10', '1e-3'), mp.mpc('1e-10', '1e-13')]

for name, bs, ws in seqs:
    h = [[mp.mpf(1)]]
    for b, w in zip(bs, ws):
        bm, wm = mp.mpf(b), mp.mpf(w)
        k_cur = len(h)
        hn = []
        for i in range(k_cur + 1):
            u = [mp.mpf(0)]
            if i < k_cur: u = poly_add_mp(u, apply_Dj_mp(h[i], bm))
            if i >= 1: u = poly_add_mp(u, [wm * c for c in h[i - 1]])
            hn.append(u)
        h = hn
    # Check min J = 4 W(h_{i+1},h_i) W(h_i,h_{i-1}) / W(h_{i+1},h_{i-1})^2
    min_J = mp.mpf('1e100')
    min_info = None
    for i in range(1, len(h) - 1):
        for t in ts_probe:
            v0, d0 = poly_eval_and_deriv_mp(h[i+1], t)
            v1, d1 = poly_eval_and_deriv_mp(h[i], t)
            v2, d2 = poly_eval_and_deriv_mp(h[i-1], t)
            w01 = v0*d1 - d0*v1
            w12 = v1*d2 - d1*v2
            w02 = v0*d2 - d0*v2
            if w02 == 0: continue
            J = 4*w01*w12 / (w02**2)
            if J < min_J:
                min_J = J
                min_info = (i, float(t))
    # Check Q_{k, q} at q = 0.5, 0.25, 0.125, 0.01
    q_res = []
    for qv in [0.5, 0.25, 0.125, 0.01]:
        qm = mp.mpf(str(qv))
        ok_all = True
        for tc in ts_uhp:
            cy = [(qm**(i*(i-1)//2)) * poly_eval_c_mp(h[i], tc) for i in range(len(h))]
            ok, _ = check_hb_pd(cy)
            if not ok:
                ok_all = False
                break
        q_res.append(f"q={qv}:{'OK' if ok_all else 'FAIL'}")
    print(f"  {name:18s}: min_J={float(min_J):+.5f} at {min_info} | " + " ".join(q_res))

# ============================================================
# TEST 2: GLZ STRENGTHENED CUBIC c_min(k) & STAGE 3 DROP
# ============================================================
print("\n=== TEST 2: GLZ STRENGTHENED CUBIC c_min(k) & STAGE 3 DROP (k=5..15) ===")
def bezoutiant_cubic_mp_float(v0, v1, v2, v3, c):
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

def find_cmin(v0, v1, v2, v3):
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

h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
h_store = {1: h}
for k in range(2, 16):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = poly_add_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = poly_add_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
    h_store[k] = h

ts_c = [mp.mpf(x) for x in [-2.0, -1.05, -0.99, -0.8, -0.5, -0.2, -0.02, 0.0, 0.2, 2.0, 20.0, 1000.0]]
for k in [5, 7, 9, 11, 13, 15]:
    hk = h_store[k]
    max_ca = 0.0
    max_cb = 0.0
    max_drop = -1e9
    worst_where = None
    for i in range(k - 2):  # i+3 <= k
        p0 = hk[i]
        p1_a = apply_Dj_mp(hk[i+1], k+1); p1_b = apply_Dj_mp(hk[i+1], k+2)
        p2_a = apply_Dj_mp(apply_Dj_mp(hk[i+2], k+1), k+1); p2_b = apply_Dj_mp(apply_Dj_mp(hk[i+2], k+2), k+2)
        p3_a = apply_Dj_mp(apply_Dj_mp(apply_Dj_mp(hk[i+3], k+1), k+1), k+1)
        p3_b = apply_Dj_mp(apply_Dj_mp(apply_Dj_mp(hk[i+3], k+2), k+2), k+2)
        for t in ts_c:
            va = [poly_eval_and_deriv_mp(p, t) for p in (p0, p1_a, p2_a, p3_a)]
            vb = [poly_eval_and_deriv_mp(p, t) for p in (p0, p1_b, p2_b, p3_b)]
            ca = find_cmin(*va)
            cb = find_cmin(*vb)
            if ca > max_ca: max_ca = ca
            if cb > max_cb: max_cb = cb
            if cb - ca > max_drop:
                max_drop = cb - ca
                worst_where = (i, float(t))
    print(f"  k={k:2d}: max c(j=k+1)={max_ca:.6f} [k(1-c)={k*(1-max_ca):.3f}], max c(j=k+2)={max_cb:.6f}, max drop={max_drop:+.6f} [k^2*drop={k*k*max_drop:+.3f}] at {worst_where}")
