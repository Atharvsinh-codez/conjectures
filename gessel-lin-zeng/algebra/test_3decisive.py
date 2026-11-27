"""
Three fast, decisive checks:
1. Strengthened Cubic Stage 1 Bezoutiant:
   - Check mixed-sign y-roots explicitly (so v_r[0] has BOTH positive and negative values!).
   - Find the exact extremal value of c_out(c_in) for c_in in {1.0, 0.95, 0.9, 0.8, 0.5}
     and inspect the extremal 5-vector configuration (v0,v1,v2,v3,v4) to guess the exact algebraic formula!
   - Also check whether B(v0,v1,v2,v3; c) >= 0 AUTOMATICALLY implies W(v0,v3) >= 0 (or can W(v0,v3) < 0?).
2. GLZ k=35 Q_{35,q} with b36 = 10^15, w36 in {1, 1e-30}, q in {1/2, 1/4, 1/8},
   using high-precision Hermite-Biehler / companion / Bezoutiant check (NO polyroots convergence failures!).
   Wait: for a complex polynomial P(y) = sum_{j=0}^k c_j y^j, all roots have Im(y) < 0 iff
   the Hermite-Biehler / Schur-Cohn / Bezoutiant matrix of (Re P, Im P) is POSITIVE DEFINITE!
   YES!! Computing the eigenvalues of a 37x37 real symmetric matrix takes 0.001 seconds and NEVER fails to converge!
3. Check the derivative sandwich A'_{k,i-1} << A_{k,i} << P_{2k-i} A_{k,i-1} on GLZ,
   AND check the EXACT direction of S_{L-1} A_{k-1,i} vs A_{k-1,i-1}!
"""
import numpy as np
import mpmath as mp
from scipy.optimize import minimize

mp.mp.dps = 80

# ============================================================
# 1. STRENGTHENED CUBIC STAGE 1: MIXED-SIGN Y-ROOTS & EXTREMIZER
# ============================================================
def bezoutiant_cubic(v0, v1, v2, v3, c):
    b = [1.0, -3.0*c, 3.0*c, -1.0]
    vs = [v0, v1, v2, v3]
    def W(i, j):
        return vs[i][0]*vs[j][1] - vs[i][1]*vs[j][0]
    B = np.zeros((3, 3), dtype=float)
    for r in range(4):
        for s in range(r+1, 4):
            val = -b[r]*b[s]*W(r, s)
            for m in range(s - r):
                B[r + m, r + (s - r - 1 - m)] += val
    return B

def min_c_psd(v0, v1, v2, v3):
    # Check if PSD at c=2.0
    if np.linalg.eigvalsh(bezoutiant_cubic(v0, v1, v2, v3, 2.0))[0] < -1e-12:
        return None
    # Note: we want the largest c <= 2 where min_eig < 0 (if PSD at c=2)
    # Wait: is PSD(c) monotone in c for W(v0,v3) >= 0? Let's verify on a grid first if needed.
    lo, hi = 1e-6, 2.0
    if np.linalg.eigvalsh(bezoutiant_cubic(v0, v1, v2, v3, lo))[0] >= 0:
        return lo
    for _ in range(60):
        mid = 0.5*(lo + hi)
        if np.linalg.eigvalsh(bezoutiant_cubic(v0, v1, v2, v3, mid))[0] >= 0:
            hi = mid
        else:
            lo = mid
    return hi

print("=== 1A. CHECK MIXED-SIGN y-ROOTS FOR CUBIC BEZOUTIANT STAGE 1 ===")
rng = np.random.default_rng(123)
mixed_valid = 0
mixed_fails = 0
max_cout_mixed = 0.0
for _ in range(30000):
    ang = np.sort(rng.uniform(0.01, np.pi - 0.01, 5))
    # Force mixed signs of cos(ang) (so v_r[0] has BOTH >0 and <0!)
    if not (ang[0] < np.pi/2 < ang[4]):
        continue
    lens = np.exp(rng.normal(0, 1.5, 5))
    vs = [(lens[i]*np.cos(ang[i]), lens[i]*np.sin(ang[i])) for i in range(5)]
    e1 = np.linalg.eigvalsh(bezoutiant_cubic(vs[0],vs[1],vs[2],vs[3],1.0))[0]
    if e1 < 0:
        continue
    e2 = np.linalg.eigvalsh(bezoutiant_cubic(vs[1],vs[2],vs[3],vs[4],1.0))[0]
    if e2 < 0:
        continue
    # Check that P1(y) actually has mixed-sign y-roots!
    p1 = [-vs[3][0], 3*vs[2][0], -3*vs[1][0], vs[0][0]]
    r1 = np.roots(p1)
    if not (np.any(r1.real > 0) and np.any(r1.real < 0)):
        continue
    mixed_valid += 1
    ws = [(vs[i][0]+vs[i+1][0], vs[i][1]+vs[i+1][1]) for i in range(4)]
    eout = np.linalg.eigvalsh(bezoutiant_cubic(ws[0],ws[1],ws[2],ws[3],1.0))[0]
    if eout < -1e-11:
        mixed_fails += 1
    co = min_c_psd(ws[0],ws[1],ws[2],ws[3])
    if co is not None and co > max_cout_mixed:
        max_cout_mixed = co
print(f"  Mixed-sign y-roots (c=1.0): valid={mixed_valid}, fails={mixed_fails}, max_cout={max_cout_mixed:.8f}")

print("\n=== 1B. EXACT EXTREMIZER OF CUBIC STAGE 1 FOR c_in IN {1.0, 0.9, 0.8, 0.6} ===")
# Note: Bezoutiant depends ONLY on the 2x2 determinants W(v_r, v_s), which are SL(2,R)-invariant!
# So mixed-sign y-roots is AUTOMATICALLY equivalent to positive y-roots by an SL(2,R) rotation
# of R^2, AS LONG AS all 5 vectors lie in an open half-plane (W(v0, v4) > 0)!
# WAIT A SECOND!! CHECK THAT!!
# If v0, v1, v2, v3, v4 lie in an open half-plane (W(v_r, v_s) > 0 for ALL 0 <= r < s <= 4),
# then there is an SL(2,R) rotation moving ALL 5 vectors into the first quadrant (v_r[0]>0, v_r[1]>0),
# where ALL y-roots of P1, P2 are POSITIVE!!
# WAIT: WHAT IF W(v0, v3) > 0 and W(v1, v4) > 0, BUT W(v0, v4) < 0 (angle span exceeds pi)?!
# TEST THAT RIGHT NOW!!
span_over_pi_valid = 0
span_over_pi_fails = 0
for _ in range(50000):
    # Pick angles in (0, 1.8*pi) with th3 - th0 < pi, th4 - th1 < pi, BUT th4 - th0 > pi!
    th0 = 0.0
    th1 = rng.uniform(0.2, 1.5)
    th2 = rng.uniform(th1 + 0.05, min(np.pi - 0.1, th1 + 1.5))
    th3 = rng.uniform(th2 + 0.05, np.pi - 0.05)
    th4 = rng.uniform(np.pi + 0.05, th1 + np.pi - 0.05)
    ang = [th0, th1, th2, th3, th4]
    lens = np.exp(rng.normal(0, 1.5, 5))
    vs = [(lens[i]*np.cos(ang[i]), lens[i]*np.sin(ang[i])) for i in range(5)]
    if np.linalg.eigvalsh(bezoutiant_cubic(vs[0],vs[1],vs[2],vs[3],1.0))[0] < 0:
        continue
    if np.linalg.eigvalsh(bezoutiant_cubic(vs[1],vs[2],vs[3],vs[4],1.0))[0] < 0:
        continue
    span_over_pi_valid += 1
    ws = [(vs[i][0]+vs[i+1][0], vs[i][1]+vs[i+1][1]) for i in range(4)]
    eout = np.linalg.eigvalsh(bezoutiant_cubic(ws[0],ws[1],ws[2],ws[3],1.0))[0]
    if eout < -1e-10:
        span_over_pi_fails += 1
print(f"  Span > pi (W(v0,v3)>0, W(v1,v4)>0, BUT W(v0,v4)<0), c=1.0: valid={span_over_pi_valid}, fails={span_over_pi_fails}")

# Now find sharp extremizer for c_in = 0.9, 0.8, 0.6 with W(v0, v4) > 0:
for c_test in [0.95, 0.90, 0.80, 0.60]:
    best_co = 0.0
    best_vs = None
    for seed in range(20):
        rng2 = np.random.default_rng(200 + seed)
        found = False
        for _ in range(1000):
            ang = np.sort(rng2.uniform(0.05, np.pi - 0.05, 5))
            lens = np.exp(rng2.normal(0, 1.5, 5))
            vs = [(lens[i]*np.cos(ang[i]), lens[i]*np.sin(ang[i])) for i in range(5)]
            if (np.linalg.eigvalsh(bezoutiant_cubic(vs[0],vs[1],vs[2],vs[3],c_test))[0] >= 1e-6 and
                np.linalg.eigvalsh(bezoutiant_cubic(vs[1],vs[2],vs[3],vs[4],c_test))[0] >= 1e-6):
                found = True
                break
        if not found:
            continue
        dth = np.diff(ang)
        logL = np.log(lens) - np.log(lens[2])
        x0 = np.concatenate([np.log(dth), [logL[0], logL[1], logL[3], logL[4]]])
        def obj(z):
            g = np.exp(np.clip(z[:4], -14, 5))
            tot = np.sum(g)
            if tot >= np.pi - 1e-4:
                g = g * (np.pi - 0.01) / (tot + 0.01)
            th = np.zeros(5)
            th[0] = 0.005
            th[1:] = th[0] + np.cumsum(g)
            l = np.clip([z[4], z[5], 0.0, z[6], z[7]], -12, 12)
            L = np.exp(l)
            vs = [(L[i]*np.cos(th[i]), L[i]*np.sin(th[i])) for i in range(5)]
            e1 = np.linalg.eigvalsh(bezoutiant_cubic(vs[0],vs[1],vs[2],vs[3],c_test))[0]
            e2 = np.linalg.eigvalsh(bezoutiant_cubic(vs[1],vs[2],vs[3],vs[4],c_test))[0]
            if e1 < 0 or e2 < 0:
                return 10.0 + 1000.0*(max(0.0, -e1) + max(0.0, -e2))
            ws = [(vs[i][0]+vs[i+1][0], vs[i][1]+vs[i+1][1]) for i in range(4)]
            co = min_c_psd(ws[0],ws[1],ws[2],ws[3])
            return -co if co is not None else -10.0
        res = minimize(obj, x0, method='Nelder-Mead', options={'maxiter': 2000})
        if -res.fun > best_co:
            best_co = -res.fun
            z = res.x
            g = np.exp(np.clip(z[:4], -14, 5))
            tot = np.sum(g)
            if tot >= np.pi - 1e-4:
                g = g * (np.pi - 0.01) / (tot + 0.01)
            th = np.zeros(5); th[0] = 0.005; th[1:] = th[0] + np.cumsum(g)
            L = np.exp(np.clip([z[4], z[5], 0.0, z[6], z[7]], -12, 12))
            best_vs = [(L[i]*np.cos(th[i]), L[i]*np.sin(th[i])) for i in range(5)]
    # Print pairwise W ratios at best_vs!
    def W(a, b): return a[0]*b[1] - a[1]*b[0]
    w01 = W(best_vs[0], best_vs[1]); w12 = W(best_vs[1], best_vs[2]); w23 = W(best_vs[2], best_vs[3]); w34 = W(best_vs[3], best_vs[4])
    w02 = W(best_vs[0], best_vs[2]); w13 = W(best_vs[1], best_vs[3]); w24 = W(best_vs[2], best_vs[4])
    d_in = 1 - c_test
    d_out = 1 - best_co
    # Compare candidate formulas:
    # Formula A: (1 + 3 c^2) / (3 + c^2)? At 0.9: 3.43/3.81 = 0.90026;
    # Formula B: (3 c + 1) / ... no
    # Formula C: 2c / (1 + c) -> at 0.9: 1.8/1.9 = 0.947
    # Formula D: (1 + 3c^2) / 4 -> at 0.9: 0.8575
    # Formula E: sqrt(c^2 + (1-c^2)/4)? No.
    # Check 4c/(3+c^2) -> at 0.9: 3.6/3.81 = 0.9448
    # Check c + (1-c)^2 / 4 -> 0.9025, 0.8100, 0.6400!
    f_cand = c_test + (1 - c_test)**2 / 4.0
    f_cand2 = 4 * c_test / (1 + c_test)**2 # wait: that's > 1
    print(f"  c_in={c_test:.2f}: best_cout={best_co:.8f}, c+(1-c)^2/4={f_cand:.8f}, 1/d_out-1/d_in={1/d_out - 1/d_in:.5f}")
    print(f"    W(01,12,23,34) norm by w12: ({w01/w12:.4f}, 1, {w23/w12:.4f}, {w34/w12:.4f}); J012={4*w01*w12/w02**2:.4f}, J123={4*w12*w23/w13**2:.4f}, J234={4*w23*w34/w24**2:.4f}")

# ============================================================
# 2. GLZ k=35 Q_{35,q} WITH HERMITE-BIEHLER BEZOUTIANT (EXACT!)
# ============================================================
print("\n=== 2. GLZ k=35 Q_{35,q} + (D_{b36} + Y_q) VIA HERMITE-BIEHLER MATRIX ===")
def apply_Dj_mp(p, j):
    n = len(p)
    out = [mp.mpf(0)] * (n + 1)
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

def poly_eval_mp(p, t):
    s = mp.mpc(0)
    for c in reversed(p):
        s = s * t + c
    return s

def check_uhp_zeros_in_lhp(coeffs_y):
    """
    coeffs_y[j] is complex coeff of y^j (j = 0..n).
    Rescale y = s*u to normalize dynamic range, then check Hermite-Biehler matrix H of (Re P, Im P)!
    Wait: all roots of P(y) have Im(y) < 0  <=>  H(P) is POSITIVE DEFINITE (or negative definite depending on sign of lc)!
    Specifically, for P(y) = A(y) + i B(y) (A, B real polys):
    Let P*(y) = conj(P(conj(y))) = A(y) - i B(y).
    Hermite matrix H_{j,k} = coeff of x^j y^k in (P(x) P*(y) - P*(x) P(y)) / (2i (x - y))
                           = (B(x) A(y) - A(x) B(y)) / (x - y).
    At x = y: B'(x) A(x) - A'(x) B(x).
    Check on P(y) = y + i (root at y = -i, in LHP!):
      A(y) = y, B(y) = 1. (B(x)A(y) - A(x)B(y))/(x-y) = (y - x)/(x - y) = -1!
      So -(B(x)A(y) - A(x)B(y))/(x-y) = (A(x)B(y) - B(x)A(y))/(x-y) = +1 > 0!
    Check on P(y) = (1 + 2t) + y at t = i*eps: A(y) = 1 + y, B(y) = 2 eps.
      (A(x)B(y) - B(x)A(y))/(x-y) = 2 eps (1+x - 1-y)/(x-y) = +2 eps > 0!
    YES! ALL roots of P(y) lie STRICTLY in Im(y) < 0  <=>  the symmetric n x n matrix
      H_{j,k} = coeff of x^j y^k in (A(x) B(y) - B(x) A(y)) / (x - y)
    is POSITIVE DEFINITE!
    AND we check positive definiteness in 80-digit mpmath via exact Cholesky / Sylvester pivots!
    """
    n = len(coeffs_y) - 1
    # geometric rescaling so |c_0| = |c_n|
    s = (abs(coeffs_y[0]) / abs(coeffs_y[-1])) ** (mp.mpf(1) / n)
    sc = [coeffs_y[j] * (s ** j) for j in range(n + 1)]
    mx = max(abs(z) for z in sc)
    sc = [z / mx for z in sc]
    A = [z.real for z in sc]
    B = [z.imag for z in sc]
    # Build n x n matrix H in mpmath
    H = [[mp.mpf(0)] * n for _ in range(n)]
    for r in range(n + 1):
        for col in range(r + 1, n + 1):
            # (A_r B_col - B_r A_col) * (x^r y^col - y^r x^col)/(x - y)
            # Note (x^r y^col - y^r x^col)/(x - y) = -(xy)^r * sum_{m=0}^{col-r-1} x^m y^{col-r-1-m}
            val = -(A[r]*B[col] - B[r]*A[col])
            for m in range(col - r):
                H[r + m][r + (col - r - 1 - m)] += val
    # Exact Cholesky in 80-digit mpmath!
    min_piv = mp.mpf('1e100')
    for i in range(n):
        piv = H[i][i]
        if piv < min_piv:
            min_piv = piv
        if piv <= 0:
            return False, float(piv), i
        for j in range(i + 1, n):
            fac = H[j][i] / piv
            for k in range(j, n):
                H[j][k] -= fac * H[i][k]
                H[k][j] = H[j][k]
    return True, float(min_piv), n

h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
for k in range(2, 36):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = poly_add_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = poly_add_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
h35 = h

for b36 in [mp.mpf(37), mp.mpf('1e15')]:
    for w36 in [mp.mpf(1), mp.mpf('1e-30')]:
        h36 = []
        for i in range(37):
            u = [mp.mpf(0)]
            if i <= 35: u = poly_add_mp(u, apply_Dj_mp(h35[i], b36))
            if i >= 1: u = poly_add_mp(u, [w36 * c for c in h35[i - 1]])
            h36.append(u)
        for q_val in [mp.mpf('0.5'), mp.mpf('0.25'), mp.mpf('0.125')]:
            ok_all = True
            worst_info = None
            for t_test in [mp.mpc('-1e-16', '1e-20'), mp.mpc('1e-16', '1e-20'), mp.mpc('-1e-6', '1e-10'), mp.mpc('10', '1e-3')]:
                coeffs = [(q_val ** (i*(i-1)//2)) * poly_eval_mp(h36[i], t_test) for i in range(37)]
                ok, piv, idx = check_uhp_zeros_in_lhp(coeffs)
                if not ok:
                    ok_all = False
                    worst_info = (t_test, piv, idx)
                    break
            print(f"  b36={float(b36):.0e}, w36={float(w36):.0e}, q={float(q_val):.3f}: {'PASS' if ok_all else 'FAIL at ' + str(worst_info)}")

# ============================================================
# 3. EXACT DIRECTION OF S_{L-1} A_{k-1,i} VS A_{k-1,i-1} ON GLZ
# ============================================================
print("\n=== 3. EXACT DIRECTION OF S_{L-1} A_{k-1,i} vs A_{k-1,i-1} ON GLZ ===")
def apply_SL_f(p, L):
    out = [0.0]*(len(p)+1)
    for m, c in enumerate(p):
        out[m] += (1+m)*c
        out[m+1] += (L-m)*c
    while len(out)>1 and abs(out[-1]) < 1e-12: out.pop()
    return out

A = {(1,0): [1.0, 2.0], (1,1): [1.0]}
for k in range(2, 8):
    for i in range(k + 1):
        L = 3*k - 2 - i
        u = [0.0]
        if i <= k - 1:
            t1 = apply_SL_f(A[(k-1,i)], L-1)
            u = [0.0]*max(len(u), len(t1))
            for m, v in enumerate(t1): u[m] += v
        if i >= 1:
            t2 = A[(k-1,i-1)]
            if len(t2) > len(u): u += [0.0]*(len(t2)-len(u))
            for m, v in enumerate(t2): u[m] += v
        A[(k,i)] = apply_SL_f(u, L)

for k in [3, 4, 5, 6]:
    for i in range(1, k):
        L = 3*k - 2 - i
        Sf = apply_SL_f(A[(k-1,i)], L-1)
        g = A[(k-1,i-1)]
        rSf = np.sort(np.roots(Sf[::-1]).real)
        rg = np.sort(np.roots(g[::-1]).real)
        # Check which is to the right!
        if k <= 4:
            print(f"  (k={k},i={i}): r(g)   = {np.round(rg, 4)}")
            print(f"           r(Sf)  = {np.round(rSf, 4)}")
        dir_g_lt_Sf = all(rg[m] < rSf[m] for m in range(len(rg)))
        dir_Sf_lt_g = all(rSf[m] < rg[m] for m in range(len(rg)))
        if k >= 5:
            print(f"  (k={k},i={i}): g << Sf? {dir_g_lt_Sf}, Sf << g? {dir_Sf_lt_g}")
