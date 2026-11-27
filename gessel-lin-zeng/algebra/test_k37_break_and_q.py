"""
5-second decisive test:
1. Using GLZ h35 (where J_17(0) = 0.9971 < 1) + b36* = 1.5862705e16:
   - With w36 = 1e-30, can we break INTERLACING at k=37 and REAL-ROOTEDNESS at k=38
     using ONLY a_j = 1, b_j >= 2 (arbitrary w_j > 0)?
   - What if w_j = 1 is FIXED? Test k=40 with b_j = 2 (constant!) and b_1 = 10^15, b_2..k = 2
     at q = 0.5, 0.4, 0.3, 0.25!
"""
import mpmath as mp
mp.mp.dps = 100

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

def poly_eval_deriv(p, t):
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
            return False, i
        for j in range(i + 1, n):
            fac = H[j][i] / piv
            for k in range(j, n):
                H[j][k] -= fac * H[i][k]
                H[k][j] = H[j][k]
    return True, n

# PART 1: Break interlacing at k=37 with a_j = 1, b_j >= 2, wild w_j > 0!
print("=== PART 1: Break interlacing at k=37 with a_j=1, b_j>=2 (wild w_j) ===")
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

# Step 36: single D_{b36*} with w36 = 1e-40 (note: GLZ used doubled D, so let's do single D_{b36*} + w36)
t_star = -mp.mpf('1e-16')
v0, d0 = poly_eval_deriv(h[16], t_star)
v1, d1 = poly_eval_deriv(h[17], t_star)
v2, d2 = poly_eval_deriv(h[18], t_star)
u0, u1, u2 = t_star*(1+t_star)*d0, t_star*(1+t_star)*d1, t_star*(1+t_star)*d2
B0 = 2*u1*v1 - (u0*v2 + u2*v0)
C0 = v1**2 - v0*v2
Phi_star = -B0 / (2*C0)
b36 = (1 - Phi_star) / (-t_star)
w36 = mp.mpf('1e-40')
h36_16 = poly_add_mp(apply_Dj_mp(h[16], b36), [w36*c for c in h[15]])
h36_17 = poly_add_mp(apply_Dj_mp(h[17], b36), [w36*c for c in h[16]])
h36_18 = poly_add_mp(apply_Dj_mp(h[18], b36), [w36*c for c in h[17]])

p0, dp0 = poly_eval_deriv(h36_16, t_star)
p1, dp1 = poly_eval_deriv(h36_17, t_star)
p2, dp2 = poly_eval_deriv(h36_18, t_star)
turan36 = p1**2 - p0*p2
W21 = p2*dp1 - dp2*p1
print(f"  h36 at t*=-1e-16: p0={float(p0):.4e}, p1={float(p1):.4e}, p2={float(p2):.4e}")
print(f"  turan36 = p1^2 - p0*p2 = {float(turan36):.4e} (< 0!), (-t*)(1+t*)W(p2,p1) = {float((-t_star)*(1+t_star)*W21):.4e}")

# For step 37: U1 = D_{b37} p1 + w37 p0, U2 = D_{b37} p2 + w37 p1.
# We want U2(t*) = 0 => w37 = -(D_{b37} p2)(t*) / p1(t*)!
# Check sign of -(D_{b37} p2)(t*) / p1(t*) for large b37 > 0:
# (D_{b37} p2)(t*) = t*(1+t*) dp2 + (1 + b37 t*) p2 ~ -b37 * 1e-16 * p2.
# So w37 ~ +b37 * 1e-16 * p2 / p1! Since p1 and p2 have the SAME sign (both < 0 or both > 0), w37 > 0!!
b37 = mp.mpf('1e30')
Dp2_val = t_star*(1+t_star)*dp2 + (1 + b37*t_star)*p2
w37 = -Dp2_val / p1
Dp1_val = t_star*(1+t_star)*dp1 + (1 + b37*t_star)*p1
U1_val = Dp1_val + w37 * p0
# Theorem A formula: U1(t*) * p1(t*) = (-t*)(1+t*) W(p2, p1) + w37 * (p1^2 - p0 p2)
prod_val = U1_val * p1
print(f"  b37=1e30 => w37 = {float(w37):+.6e} (> 0!)")
print(f"  U2(t*) = 0, U1(t*)*p1(t*) = {float(prod_val):+.6e} (< 0 MEANS INTERLACING U2 << U1 IS BROKEN!!)")

# Also check W(U2, U1)(t*) to be 100% certain!
U1_poly = poly_add_mp(apply_Dj_mp(h36_17, b37), [w37*c for c in h36_16])
U2_poly = poly_add_mp(apply_Dj_mp(h36_18, b37), [w37*c for c in h36_17])
u1v, du1v = poly_eval_deriv(U1_poly, t_star)
u2v, du2v = poly_eval_deriv(U2_poly, t_star)
W_U2_U1 = u2v*du1v - du2v*u1v
print(f"  EXACT check: U2(t*)={float(u2v):.3e}, U1(t*)={float(u1v):.3e}, W(U2, U1)(t*) = {float(W_U2_U1):+.6e} (< 0 => INTERLACING BROKEN!)")
# Wait: check sign of W(p2, p1) first to make sure W(U2, U1) should be > 0 for interlacing!
print(f"  Compare sign: W(h36_18, h36_17)(t*) = {float(W21):+.6e}")

# PART 2: Check FIXED w_j = 1 at k=25 with b_1 = 1e15, b_2..25 = 2 (our exact boundary attack for q=1/2!)
print("\n=== PART 2: Fixed w_j = 1, b_1 = 1e15, b_2..k = 2 (k = 10, 20, 30) ===")
for k_tot in [10, 20, 30]:
    h_f = [[mp.mpf(1)]]
    bs = [mp.mpf('1e15')] + [mp.mpf(2)]*(k_tot - 1)
    for b in bs:
        kc = len(h_f)
        hn = []
        for i in range(kc + 1):
            u = [mp.mpf(0)]
            if i < kc: u = poly_add_mp(u, apply_Dj_mp(h_f[i], b))
            if i >= 1: u = poly_add_mp(u, h_f[i - 1])
            hn.append(u)
        h_f = hn
    res = []
    for qv in [0.5, 0.45, 0.4, 0.35, 0.3, 0.25]:
        qm = mp.mpf(str(qv))
        ok_all = True
        for tc in [mp.mpc('0', '1e-12'), mp.mpc('-0.5', '1e-6'), mp.mpc('10', '1e-4'), mp.mpc('1000', '1e-2')]:
            cy = [(qm**(i*(i-1)//2)) * poly_eval_c_mp(h_f[i], tc) for i in range(len(h_f))]
            ok, idx = check_hb_pd(cy)
            if not ok:
                ok_all = False
                break
        res.append(f"q={qv}:{'OK' if ok_all else 'FAIL'}")
    print(f"  k={k_tot:2d} (b=[1e15, 2,...]): " + " ".join(res))
