"""
Fast 5-second checks:
1. Exact q_max(k) for G_{k,q}^-(t, y) = sum_{i=0}^{k-1} (-1)^i q^{C(i,2)} E_{k,i}(t) y^i (k = 3..16),
   AND compare q_max BEFORE Stage 3 (G~_{k,q}^- with D_{k+1}) vs AFTER Stage 3 (G_{k+1,q}^- with D_{k+2})!
2. Check if sum_{i=0}^k A_{k,i}(lambda^i x) y^i (or with (-1)^i) is real-stable for lambda = 2 (or other lambda)!
3. Check if G_k^-(t, y) with BINOMIAL multiplier sum_{i=0}^{k-1} (-1)^i C(k-1, i) E_{k,i}(t) y^i
   (or 1/i!, or i!) is real-stable for k = 2..16!
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 70

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
        if piv <= 0:
            return False
        for j in range(i + 1, n):
            fac = H[j][i] / piv
            for k in range(j, n):
                H[j][k] -= fac * H[i][k]
                H[k][j] = H[j][k]
    return True

ts_uhp = [mp.mpc('-2.0', '1e-4'), mp.mpc('-0.99', '1e-5'), mp.mpc('-0.8', '1e-5'),
          mp.mpc('-0.5', '1e-5'), mp.mpc('-0.2', '1e-5'), mp.mpc('-0.01', '1e-6'),
          mp.mpc('0.0', '1e-8'), mp.mpc('0.2', '1e-5'), mp.mpc('5.0', '1e-4'), mp.mpc('100.0', '1e-3')]

def find_qmax(E_list):
    # bisect q in [0.1, 1.0] for sum (-1)^i q^{C(i,2)} E_i(t) y^i
    # precompute E_i(tc)
    evals = [[peval_c(p, tc) for p in E_list] for tc in ts_uhp]
    def ok_q(qv):
        qm = mp.mpf(qv)
        n = len(E_list)
        for ev in evals:
            cy = [((-1)**i) * (qm**(i*(i-1)//2)) * ev[i] for i in range(n)]
            if not check_hb_pd_lhp(cy):
                return False
        return True
    if ok_q(1.0): return 1.0
    lo, hi = 0.05, 1.0
    if not ok_q(lo): return 0.0
    for _ in range(25):
        mid = 0.5*(lo + hi)
        if ok_q(mid): lo = mid
        else: hi = mid
    return lo

print("=== 1. q_max(k) OF G_{k,q}^- (AND STAGE 3 SHIFT D_{k+1} -> D_{k+2}) + C(k-1,i)*E_{k,i} ===")
h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
for k in range(2, 16):
    j = k + 1
    E_cur = []
    for i in range(k):
        p = h[i]
        for _ in range(i): p = apply_Dj_mp(p, j)
        E_cur.append(p)
    # step to h_k
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = padd_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = padd_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    # Build G~_k (with D_{k+1}^i h_k) and G_{k+1} (with D_{k+2}^i h_k)!
    E_tilde = []
    E_next = []
    for i in range(k + 1):
        pa = hn[i]; pb = hn[i]
        for _ in range(i):
            pa = apply_Dj_mp(pa, j)
            pb = apply_Dj_mp(pb, j + 1)
        E_tilde.append(pa)
        E_next.append(pb)
    if k in [2, 3, 4, 5, 6, 8, 10, 12, 14]:
        q_before = find_qmax(E_tilde)
        q_after  = find_qmax(E_next)
        # Also test C(k, i) * E_tilde (and i! * E_tilde) at q = 1!
        ev_t = [[peval_c(p, tc) for p in E_tilde] for tc in ts_uhp]
        ok_bin = all(check_hb_pd_lhp([((-1)**i)*mp.binomial(k, i)*ev[i] for i in range(k+1)]) for ev in ev_t)
        print(f"  k={k:2d}: q_max(D_{{k+1}}) = {q_before:.6f}, q_max(D_{{k+2}}) = {q_after:.6f}, diff={q_after-q_before:+.6f} | C(k,i)*E_i: {'PASS' if ok_bin else 'FAIL'}")
    h = hn

print("\n=== 2. CHECK sum A_{k,i}(2^i x) y^i (and (-1)^i) FOR REAL STABILITY ===")
def apply_SL(p, L):
    out = [mp.mpf(0)]*(len(p)+1)
    for m, c in enumerate(p):
        out[m] += (1+m)*c; out[m+1] += (L-m)*c
    while len(out)>1 and out[-1]==0: out.pop()
    return out

A = {(1,0): [mp.mpf(1), mp.mpf(1)], (1,1): [mp.mpf(1)]}
for k in range(2, 11):
    for i in range(k + 1):
        L = 3*k - 2 - i
        u = [mp.mpf(0)]
        if i <= k - 1: u = padd_mp(u, apply_SL(A[(k-1, i)], L - 1))
        if i >= 1: u = padd_mp(u, A[(k-1, i-1)])
        A[(k, i)] = apply_SL(u, L)

xs_uhp = [mp.mpc('-2.0', '1e-3'), mp.mpc('-0.5', '1e-4'), mp.mpc('0.0', '1e-6'), mp.mpc('1.0', '1e-4'), mp.mpc('10.0', '1e-3')]
for lam in [2.0, 1.5, 4.0]:
    lm = mp.mpf(str(lam))
    res_p = []
    for k in [2, 3, 4, 6, 8, 10]:
        ok_pos = True
        for xc in xs_uhp:
            cy = [peval_c(A[(k, i)], (lm**i) * xc) for i in range(k + 1)]
            if not check_hb_pd_lhp(cy):
                ok_pos = False; break
        res_p.append(f"k={k}:{'OK' if ok_pos else 'FAIL'}")
    print(f"  sum A_{{k,i}}({lam}^i x) y^i: " + " ".join(res_p))
