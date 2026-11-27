"""
Fast high-precision checks (< 15s total):
1. GLZ J_raw(t) < 1 zone at k = 35, 45, 55:
   - Where on (-1, 0) does J_raw(t) < 1 actually occur?
   - Whenever J_raw(t) < 1 on (-1, 0), is t ALWAYS > -1/(k+1) AND B_0(t) > 0 AND A_0(t) > 0?!
   - Does GLZ plain Turan h~_{k,i}(t)^2 - h~_{k,i-1}(t) h~_{k,i+1}(t) >= 0 hold on (-1, 0) at k=35, 45, 55?
2. Coefficientwise log-concavity of A_{k,i}(x)^2 - A_{k,i-1}(x) A_{k,i+1}(x) in x (k = 2..12),
   AND Bernstein coefficients of Delta_{k,i}(t) = h~_{k,i}(t)^2 - h~_{k,i-1}(t) h~_{k,i+1}(t) on [-1, 0]!
3. Stress-test B_k^{(E)}(t, y) = sum_{i=0}^k (-1)^i C(k, i) (D_{k+1}^i h~_{k,i}(t)) y^i up to k = 26!
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
        if piv <= 0:
            return False
        for j in range(i + 1, n):
            fac = H[j][i] / piv
            for k in range(j, n):
                H[j][k] -= fac * H[i][k]
                H[k][j] = H[j][k]
    return True

# ============================================================
# 1. BUILD GLZ h~_{k, i} UP TO k = 50 AND CHECK J_raw < 1 ZONE
# ============================================================
print("=== 1. GLZ J_raw < 1 ZONE & (A0, B0, C0) SIGNS AT k = 35, 45, 50 ===")
h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
h_store = {}
for k in range(2, 51):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = padd_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = padd_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
    if k in [18, 22, 26, 35, 45, 50]:
        h_store[k] = h

for k in [35, 45, 50]:
    hk = h_store[k]
    # Check middle i = k//2 across log-spaced negative t from -0.999 to -1e-35
    i = k // 2
    f0, f1, f2 = hk[i - 1], hk[i], hk[i + 1]
    # Find rightmost root scale ~ 1 / (f1'(0)/f1(0))
    v1_0, d1_0 = eval_pd_mp(f1, mp.mpf(0))
    scale = v1_0 / d1_0  # approx |rightmost root|
    t_grid = [-mp.mpf(x) for x in [0.999, 0.99, 0.9, 0.5, 0.1, 0.03, 0.01]]
    # Add log-spaced points from -1/(k+1) down to -1e-5 * scale
    for pwr in range(1, 32):
        tv = -mp.mpf(10)**(-pwr)
        if abs(tv) > scale * mp.mpf('1e-4'):
            t_grid.append(tv)
    for mult in [100, 30, 10, 3, 1.5, 1.0, 0.7, 0.5, 0.3, 0.1, 0.01]:
        t_grid.append(-mp.mpf(str(mult)) * scale)
    t_grid.sort()
    j_lt_1 = []
    min_C0_norm = mp.mpf('1e100')
    bad_B0_when_Jlt1 = 0
    bad_t_when_Jlt1 = 0
    for t in t_grid:
        v0, d0 = eval_pd_mp(f0, t)
        v1, d1 = eval_pd_mp(f1, t)
        v2, d2 = eval_pd_mp(f2, t)
        w21 = v2*d1 - d2*v1
        w10 = v1*d0 - d1*v0
        w20 = v2*d0 - d2*v0
        J = 4*w21*w10 / (w20**2)
        u0, u1, u2 = t*(1+t)*d0, t*(1+t)*d1, t*(1+t)*d2
        A0 = u1**2 - u0*u2
        B0 = 2*u1*v1 - (u0*v2 + u2*v0)
        C0 = v1**2 - v0*v2
        # Check D_{k+1} Turan at Phi = 1 + (k+1)t!
        Phi = 1 + (k + 1)*t
        DT = A0 + Phi*B0 + (Phi**2)*C0
        if C0 / (v1**2 + abs(v0*v2)) < min_C0_norm:
            min_C0_norm = C0 / (v1**2 + abs(v0*v2))
        if J < 1:
            j_lt_1.append((float(t), float(t/scale), float(J), A0 > 0, B0 > 0, C0 > 0, DT > 0))
            if B0 <= 0: bad_B0_when_Jlt1 += 1
            if t <= -mp.mpf(1)/(k + 1): bad_t_when_Jlt1 += 1
    print(f"  k={k}, i={i} (scale={float(scale):.2e}): min C0_norm={float(min_C0_norm):+.4e}, #pts(J<1)={len(j_lt_1)}, bad_B0={bad_B0_when_Jlt1}, bad_t(<-1/(k+1))={bad_t_when_Jlt1}")
    if j_lt_1:
        # Print leftmost (most negative t) point where J < 1!
        print(f"    leftmost J<1 pt: t={j_lt_1[0][0]:.3e} (t/scale={j_lt_1[0][1]:.3f}), J={j_lt_1[0][2]:.6f}, (A0>0,B0>0,C0>0,DT>0)={j_lt_1[0][3:]}")
        print(f"    rightmost J<1 pt: t={j_lt_1[-1][0]:.3e} (t/scale={j_lt_1[-1][1]:.3f}), J={j_lt_1[-1][2]:.6f}, (A0>0,B0>0,C0>0,DT>0)={j_lt_1[-1][3:]}")

# ============================================================
# 2. COEFFICIENTWISE LOG-CONCAVITY OF A_{k,i}(x) & BERNSTEIN OF Delta_{k,i}(t)
# ============================================================
print("\n=== 2. COEFFICIENTWISE LOG-CONCAVITY OF A_{k,i}(x) IN x AND BERNSTEIN OF Delta_{k,i}(t) ===")
def apply_SL_int(p, L):
    out = [0]*(len(p)+1)
    for m, c in enumerate(p):
        out[m] += (1+m)*c; out[m+1] += (L-m)*c
    while len(out)>1 and out[-1]==0: out.pop()
    return out

def padd_int(a, b):
    n = max(len(a), len(b))
    out = [0]*n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return out

def pmul_int(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out

def psub_int(a, b):
    n = max(len(a), len(b))
    out = [0]*n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] -= b[i]
    return out

A_int = {(1,0): [1, 1], (1,1): [1]}
for k in range(2, 13):
    for i in range(k + 1):
        L = 3*k - 2 - i
        u = [0]
        if i <= k - 1: u = padd_int(u, apply_SL_int(A_int[(k-1, i)], L - 1))
        if i >= 1: u = padd_int(u, A_int[(k-1, i-1)])
        A_int[(k, i)] = apply_SL_int(u, L)

# Check A_{k,i}(x)^2 - A_{k,i-1}(x) A_{k,i+1}(x) in x:
neg_x = 0; tot_x = 0; first_neg_x = None
# Also check Bernstein of Delta(t) = (1+t)^D * (A_i(t/(1+t))^2 - A_{i-1} A_{i+1}) on [-1, 0]!
# Wait: setting t = -u/(1) -> x = -u/(1-u):
#   Delta(t) = sum_m c_m t^m (1+t)^{D-m} = sum_m (-1)^m c_m u^m (1-u)^{D-m}!!
# WHOA!! LOOK AT THAT!!
# Since t = -u and 1+t = 1-u, the Bernstein coefficients of Delta(t) on [-1, 0] ARE LITERALLY
#   (-1)^m * [x^m](A_i(x)^2 - A_{i-1}(x) A_{i+1}(x))!!
# Wait: check! A(t/(1+t)) = sum a_m t^m (1+t)^{-m}, so (1+t)^D * [A_i^2 - A_{i-1} A_{i+1}](t/(1+t))
#   = sum_m d_m t^m (1+t)^{D - m} = sum_m (-1)^m d_m u^m (1-u)^{D - m}!!
# So for ALL Bernstein coeffs to be >= 0 on [-1,0], d_m would have to ALTERNATE in sign (-1)^m d_m >= 0 —
# whereas in x > 0, strong q-log-concavity means d_m >= 0 (all POSITIVE)! Let's check what signs d_m has!
for k in range(2, 13):
    for i in range(1, k):
        tot_x += 1
        d = psub_int(pmul_int(A_int[(k, i)], A_int[(k, i)]), pmul_int(A_int[(k, i-1)], A_int[(k, i+1)]))
        if min(d) < 0:
            neg_x += 1
            if first_neg_x is None:
                first_neg_x = (k, i, [x for x in d[:6]])
print(f"  A_{{k,i}}(x)^2 - A_{{k,i-1}}(x) A_{{k,i+1}}(x) coeffwise >= 0 in x (k=2..12): neg={neg_x}/{tot_x}, first_neg={first_neg_x}")

# ============================================================
# 3. STRESS-TEST B_k^{(E)}(t, y) = sum (-1)^i C(k,i) (D_{k+1}^i h~_{k,i}(t)) y^i AT k = 18, 22, 26!
# ============================================================
print("\n=== 3. STRESS-TEST B_k^{(E)}(t, y) = sum (-1)^i C(k,i) (D_{k+1}^i h~_{k,i}(t)) y^i ===")
ts_uhp = [mp.mpc('-2.0', '1e-4'), mp.mpc('-0.99', '1e-5'), mp.mpc('-0.5', '1e-5'),
          mp.mpc('-0.01', '1e-6'), mp.mpc('0.0', '1e-8'), mp.mpc('2.0', '1e-4'), mp.mpc('50.0', '1e-3')]
for k in [18, 22, 26]:
    hk = h_store[k]
    E = []
    for i in range(k + 1):
        p = hk[i]
        for _ in range(i): p = apply_Dj_mp(p, k + 1)
        E.append(p)
    ok_all = True
    fail_t = None
    for tc in ts_uhp:
        cy = [((-1)**i) * mp.binomial(k, i) * peval_c(E[i], tc) for i in range(k + 1)]
        if not check_hb_pd_lhp(cy):
            ok_all = False; fail_t = float(tc.real); break
    print(f"  k={k:2d}: B_k^(E)(t, y): {'PASS' if ok_all else 'FAIL at t=' + str(fail_t)}")
