"""
Three crisp 5-second checks:
1. Ma-Wang (Oct 2026) Remark 1 on the AGGREGATE A_{k-1, 4}:
   Does S_{3k-7} A_{k-1, 4} << A_{k-2, 0}  (i.e. D_{k+1} h~_{k-1, 4} << h~_{k-2, 0})
   HOLD or FAIL for k = 5..35?!
   (We just proved individual A_{k-1, {1,2,3,4}} FAILS at k=10! What about the aggregate A_{k-1, 4}?)
2. Check Term I, Term II, Term III of U_{k, i}^2 - U_{k, i-1} U_{k, i+1} on (-1, 0) for GLZ (k = 4..14)!
3. Check B_k^(E)(t, y) = sum_{i=0}^k (-1)^i C(k,i) (D_{k+1}^i h~_{k,i}(t)) y^i at k = 30, 35
   (wait: also check if B_k^(E) passes at t = -1e-16 + 1e-20 j, where plain Q_k failed at k=35!).
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 80

def apply_Dj_mp(poly, j):
    out = [mp.mpf(0)] * (len(poly) + 1)
    for m, c in enumerate(poly):
        if c != 0:
            out[m] += (1 + m) * c
            out[m + 1] += (j + m) * c
    return out

def padd_mp(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for idx in range(len(a)): out[idx] += a[idx]
    for idx in range(len(b)): out[idx] += b[idx]
    return out

def eval_pd_mp(poly, tv):
    v = mp.mpf(0); dv = mp.mpf(0)
    for c in reversed(poly):
        dv = dv * tv + v
        v = v * tv + c
    return v, dv

def peval_c(poly, t):
    s = mp.mpc(0)
    for c in reversed(poly):
        s = s * t + c
    return s

def check_hb_pd_lhp(coeffs_y):
    n = len(coeffs_y) - 1
    if n <= 0: return True
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
        if piv <= 0: return False
        for j in range(i + 1, n):
            fac = H[j][i] / piv
            for k in range(j, n):
                H[j][k] -= fac * H[i][k]
                H[k][j] = H[j][k]
    return True

# Build GLZ h~_{k, i} up to k = 35
h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
h_store = {1: h}
for k in range(2, 36):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = padd_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = padd_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
    h_store[k] = h

# ============================================================
# 1. MA-WANG REMARK 1 ON AGGREGATE A_{k-1, 4}: D_{k+1} h~_{k-1, 4} << h~_{k-2, 0}
# ============================================================
print("=== 1. MA-WANG REMARK 1 (AGGREGATE): D_{k+1} h~_{k-1, 4} << h~_{k-2, 0} (k = 5..35) ===")
# Check W(D_{k+1} h~_{k-1, 4}, h~_{k-2, 0})(t) > 0 on fine grid of (-1, 0) (and exact Bezoutiant!)
# WAIT: f << g (deg g = deg f + 1, both RZ(-1,0), positive lc) <=> Hermite-Biehler matrix of (g, f) is PD!
# EVEN BETTER: check Bezoutiant of (g, f) in 80-digit mpmath — 100% EXACT NECESSARY & SUFFICIENT FOR f << g!
def check_interlace_bezoutiant(f, g):
    """Exact 80-digit check of f << g where deg g = deg f + 1, positive lc!
    Wait: for g (deg n) and f (deg n-1) with positive lc, f << g <=>
      (g(x) f(y) - g(y) f(x)) / (x - y)  has POSITIVE DEFINITE symmetric n x n matrix!
    Check on g(t) = (t+0.5)(t+0.2) = 0.1 + 0.7t + t^2, f(t) = t + 0.3 = 0.3 + t:
      g(x)f(y) - g(y)f(x) at y=x -> g'(x)f(x) - g(x)f'(x) = -W(g, f) = +W(f, g) > 0!
    Wait: check (g(x)f(y) - g(y)f(x))/(x-y) at x->y:
      (g(y)+(x-y)g'(y))f(y) - g(y)(f(y)+(x-y)f'(y)) / (x-y) = g'(y)f(y) - g(y)f'(y) = W(f, g)(y) > 0!
    YES!! Positive definite Bezoutiant of (g, f) <=> f << g (100% EXACT, NO ROOT-FINDING NEEDED!)
    """
    n = len(g) - 1
    assert len(f) == n
    # rescale t = s*u so |g[0]| ~ |g[n]|
    s = (abs(g[0]) / abs(g[-1])) ** (mp.mpf(1) / n)
    gs = [g[m] * (s ** m) for m in range(n + 1)]
    fs = [f[m] * (s ** m) for m in range(n)] + [mp.mpf(0)]
    mxg = max(abs(x) for x in gs); mxf = max(abs(x) for x in fs)
    gs = [x / mxg for x in gs]; fs = [x / mxf for x in fs]
    H = [[mp.mpf(0)] * n for _ in range(n)]
    for r in range(n + 1):
        for col in range(r + 1, n + 1):
            # (gs[r]*fs[col] - gs[col]*fs[r]) * (x^r y^col - y^r x^col)/(x - y)
            val = -(gs[r]*fs[col] - gs[col]*fs[r])
            for m in range(col - r):
                H[r + m][r + (col - r - 1 - m)] += val
    for i in range(n):
        piv = H[i][i]
        if piv <= 0: return False, i
        for j in range(i + 1, n):
            fac = H[j][i] / piv
            for k in range(j, n):
                H[j][k] -= fac * H[i][k]
                H[k][j] = H[j][k]
    return True, n

for k in [5, 6, 8, 10, 12, 15, 20, 25, 30, 35]:
    lhs = apply_Dj_mp(h_store[k - 1][4], k + 1)
    rhs = h_store[k - 2][0]
    ok, idx = check_interlace_bezoutiant(lhs, rhs)
    # Also check i=4 candidate h~_{k, 4} << D_k h~_{k-2, 0} (aggregate H_{k-1}!)
    ok_H, idx_H = check_interlace_bezoutiant(h_store[k][4], apply_Dj_mp(h_store[k - 2][0], k))
    # Also check i=5 aggregate h~_{k, 5} << h~_{k-2, 0}!
    ok_5 = check_interlace_bezoutiant(h_store[k][5], h_store[k - 2][0])[0] if k >= 6 else None
    print(f"  k={k:2d}: MW-Rem1(agg D_{{k+1}} h_{{k-1,4}} << h_{{k-2,0}}): {'PASS' if ok else f'FAIL@{idx}'} | h_{{k,4}}<<H_{{k-1}}: {'PASS' if ok_H else f'FAIL@{idx_H}'} | h_{{k,5}}<<h_{{k-2,0}}: {ok_5}")

# ============================================================
# 2. CHECK TERM I, TERM II, TERM III OF U_{k,i}^2 - U_{k,i-1} U_{k,i+1} ON (-1, 0)
# ============================================================
print("\n=== 2. SIGNS OF TERM I, TERM II, TERM III OF U-TURAN ON (-1, 0) (k = 4..12) ===")
ts_u = [mp.mpf(x) for x in np.linspace(-0.995, -0.005, 80)]
for k in [4, 6, 8, 10, 12]:
    h_prev = h_store[k - 1]
    j = k + 1
    neg_T1 = 0; neg_T3 = 0; neg_U = 0; tot = 0
    for i in range(2, k - 1):
        # U_i = D_j h_i + h_{i-1}
        a_m1, a_0, a_1, a_2 = h_prev[i-2], h_prev[i-1], h_prev[i], h_prev[i+1]
        b_0, b_1, b_2 = apply_Dj_mp(a_0, j), apply_Dj_mp(a_1, j), apply_Dj_mp(a_2, j)
        for t in ts_u:
            tot += 1
            am1 = eval_pd_mp(a_m1, t)[0]; a0 = eval_pd_mp(a_0, t)[0]; a1 = eval_pd_mp(a_1, t)[0]
            b0 = eval_pd_mp(b_0, t)[0]; b1 = eval_pd_mp(b_1, t)[0]; b2 = eval_pd_mp(b_2, t)[0]
            T1 = b1**2 - b0*b2
            T2 = a0**2 - am1*a1
            T3 = 2*a0*b1 - am1*b2 - a1*b0
            if T1 < 0: neg_T1 += 1
            if T3 < 0: neg_T3 += 1
            if T1 + T2 + T3 < 0: neg_U += 1
    print(f"  k={k:2d} (tot={tot:4d}): neg_T1={neg_T1}, neg_T3={neg_T3}, neg_U={neg_U}")

# ============================================================
# 3. STRESS-TEST B_k^(E)(t, y) AT k = 30 AND k = 35!
# ============================================================
print("\n=== 3. B_k^(E)(t, y) = sum (-1)^i C(k,i) (D_{k+1}^i h~_{k,i}(t)) y^i AT k = 30, 35 ===")
ts_uhp = [mp.mpc('-2.0', '1e-4'), mp.mpc('-0.99', '1e-5'), mp.mpc('-0.5', '1e-5'),
          mp.mpc('-1e-16', '1e-20'), mp.mpc('0.0', '1e-12'), mp.mpc('5.0', '1e-4')]
for k in [30, 35]:
    hk = h_store[k]
    E = []
    for i in range(k + 1):
        p = hk[i]
        for _ in range(i): p = apply_Dj_mp(p, k + 1)
        E.append(p)
    ok_all = True; fail_t = None
    for tc in ts_uhp:
        cy = [((-1)**i) * mp.binomial(k, i) * peval_c(E[i], tc) for i in range(k + 1)]
        if not check_hb_pd_lhp(cy):
            ok_all = False; fail_t = float(tc.real); break
    print(f"  k={k:2d}: B_k^(E)(t, y): {'PASS' if ok_all else 'FAIL at t=' + str(fail_t)}")
