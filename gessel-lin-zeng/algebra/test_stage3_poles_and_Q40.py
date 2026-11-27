"""
Two high-leverage investigations (< 15s total):
1. ROUTE I STAGE 3 ANATOMY:
   Compute the EXACT 6 j-poles (a1; b1, b2, b3; c1, c2) of J_j(t) = 4 W(f0, D_j f1) W(D_j f1, D_j^2 f2) / W(f0, D_j^2 f2)^2
   on GLZ (k = 6..18), check J_infty(t) = J(f0, t f1, t^2 f2) on R (even up to k=35!),
   and compute S(t) = (a1 + b1 + b2 + b3) - 2(c1 + c2) to see WHY Stage 3 drop is O(1/k^2)!
2. ROUTE II ONE-STEP CLOSURE TEST:
   - Try to BREAK prod_{j=1}^k (D_{b_j} + Y_q)(1) with FIXED w_j = 1, b_j >= 2 at q = 1/4, 1/8
     using GLZ k=40 (where J(0) < 0.96) + nasty b_41, b_42, b_43!
   - Also test ONE-STEP preservation: take a random real-stable F(t, y) = sum_{i=0}^d f_i(t) y^i
     with strong y-root separation (say all y-roots separated by >= 1/q) and check if
     D_b F(t, y) + y F(t, q y) is ALWAYS real-stable for b >= 2!
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

# Build GLZ h_k up to k = 40
h_store = {}
h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
h_store[1] = h
for k in range(2, 41):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = poly_add_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = poly_add_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
    if k in [6, 10, 14, 18, 25, 35, 40]:
        h_store[k] = h

# ============================================================
# PART 1A: CHECK J_infty(t) = J(f0, t f1, t^2 f2) ON GLZ UP TO k=40!
# ============================================================
print("=== PART 1A: J_infty(t) = J(h_{k,i-2}, t h_{k,i-1}, t^2 h_{k,i}) ON GLZ ===")
ts = [mp.mpf(x) for x in [-5.0, -1.5, -1.01, -0.99, -0.8, -0.5, -0.2, -0.05, -0.005, -1e-6, -1e-12, 0.0, 0.01, 0.2, 2.0, 20.0, 1000.0]]
for k in [6, 10, 14, 18, 25, 35, 40]:
    hk = h_store[k]
    min_Jinf = mp.mpf('1e100')
    where = None
    neg_W = 0
    for i in range(2, k + 1):
        f0, f1, f2 = hk[i - 2], hk[i - 1], hk[i]
        for t in ts:
            v0, d0 = poly_eval_deriv(f0, t)
            v1, d1 = poly_eval_deriv(f1, t)
            v2, d2 = poly_eval_deriv(f2, t)
            # W(f0, t f1) = t(v0 d1 - d0 v1) + v0 v1
            # W(t f1, t^2 f2) / t^2 = t(v1 d2 - d1 v2) + v1 v2
            # W(f0, t^2 f2) / t = t(v0 d2 - d0 v2) + 2 v0 v2
            w01 = t*(v0*d1 - d0*v1) + v0*v1
            w12 = t*(v1*d2 - d1*v2) + v1*v2
            w02 = t*(v0*d2 - d0*v2) + 2*v0*v2
            if w01 <= 0 or w12 <= 0 or w02 <= 0:
                neg_W += 1
            if w02 != 0:
                Jinf = 4*w01*w12 / (w02**2)
                if Jinf < min_Jinf:
                    min_Jinf = Jinf
                    where = (i, float(t))
    print(f"  k={k:2d}: min J_infty = {float(min_Jinf):.6f} [k*(J-1)={k*(float(min_Jinf)-1):.3f}] at {where}, neg_W={neg_W}")

# ============================================================
# PART 1B: EXACT 6 j-POLES OF J_j(t) ON GLZ (k = 6, 10, 14, 18)
# ============================================================
print("\n=== PART 1B: EXACT 6 j-POLES (a1; b1,b2,b3; c1,c2) AND S = a1+sum(b)-2*sum(c) ===")
def mul_t(p):
    return [mp.mpf(0)] + list(p)

for k in [6, 10, 14, 18]:
    hk = h_store[k]
    i = k // 2 + 1
    f0, f1, f2 = hk[i - 2], hk[i - 1], hk[i]
    # Decompose D_j f1 = D_0 f1 + j (t f1)
    # D_j^2 f2 = D_0^2 f2 + j (D_0(t f2) + t D_0 f2) + j^2 (t^2 f2)
    D0_f1 = apply_Dj_mp(f1, 0)
    t_f1 = mul_t(f1)
    D0_f2 = apply_Dj_mp(f2, 0)
    D02_f2 = apply_Dj_mp(D0_f2, 0)
    t_f2 = mul_t(f2)
    lin_f2 = poly_add_mp(apply_Dj_mp(t_f2, 0), mul_t(D0_f2))
    t2_f2 = mul_t(t_f2)

    print(f"  --- k={k}, i={i} (j=k+1={k+1}) ---")
    for t_val in [-0.8, -0.5, -0.1, 0.2, 2.0, 50.0]:
        t = mp.mpf(str(t_val))
        vX = poly_eval_deriv(f0, t)
        vM0 = poly_eval_deriv(D0_f1, t)
        vM1 = poly_eval_deriv(t_f1, t)
        vY0 = poly_eval_deriv(D02_f2, t)
        vY1 = poly_eval_deriv(lin_f2, t)
        vY2 = poly_eval_deriv(t2_f2, t)
        def W(a, b): return a[0]*b[1] - a[1]*b[0]
        # W(X, M_j) = p1 j + p0
        p1, p0 = W(vX, vM1), W(vX, vM0)
        a1 = float(-p0 / p1)
        # W(X, Y_j) = q2 j^2 + q1 j + q0
        q2, q1, q0 = W(vX, vY2), W(vX, vY1), W(vX, vY0)
        cs = np.roots([float(q2/q2), float(q1/q2), float(q0/q2)])
        # W(M_j, Y_j) = r3 j^3 + r2 j^2 + r1 j + r0
        r3 = W(vM1, vY2)
        r2 = W(vM1, vY1) + W(vM0, vY2)
        r1 = W(vM1, vY0) + W(vM0, vY1)
        r0 = W(vM0, vY0)
        bs = np.roots([1.0, float(r2/r3), float(r1/r3), float(r0/r3)])
        S_val = a1 + np.sum(bs).real - 2*np.sum(cs).real
        # Also exact d/dj ln J at j = k+1!
        j0 = k + 1
        dlnJ = (1/(j0-a1) + np.sum(1/(j0-bs)) - 2*np.sum(1/(j0-cs))).real
        bs_s = np.sort(bs.real)
        cs_s = np.sort(cs.real)
        print(f"    t={t_val:+5.1f}: a1={a1:+7.2f}, c=({cs_s[0]:+6.2f},{cs_s[1]:+6.2f}), b=({bs_s[0]:+6.2f},{bs_s[1]:+6.2f},{bs_s[2]:+6.2f}) | S={S_val:+7.3f}, k^2*dlnJ={k*k*dlnJ:+7.3f}")

# ============================================================
# PART 2: TEST FIXED w_j = 1 AFTER GLZ k=40 WITH NASTY b_{41}, b_{42}, b_{43}
# ============================================================
print("\n=== PART 2: GLZ k=40 + nasty b41, b42, b43 WITH FIXED w_j = 1 ===")
h40 = h_store[40]
# Find exact b41* where D_{b41*} h40 breaks Turan at t = -1e-18 (or -1e-20)!
# Compute slope at t=0 of h40[20]
v20, d20 = poly_eval_deriv(h40[20], mp.mpf(0))
slope = d20 / v20
t_probe = -mp.mpf('0.5') / slope
print(f"  GLZ k=40: h40[20]'(0)/h40[20](0) = {float(slope):.3e}, probing t = {float(t_probe):.3e}")
f0, f1, f2 = h40[19], h40[20], h40[21]
v0, d0 = poly_eval_deriv(f0, t_probe)
v1, d1 = poly_eval_deriv(f1, t_probe)
v2, d2 = poly_eval_deriv(f2, t_probe)
u0, u1, u2 = t_probe*(1+t_probe)*d0, t_probe*(1+t_probe)*d1, t_probe*(1+t_probe)*d2
B0 = 2*u1*v1 - (u0*v2 + u2*v0)
C0 = v1**2 - v0*v2
Phi_star = -B0 / (2*C0)
b41 = (1 - Phi_star) / (-t_probe)
print(f"  Exact critical b41* = {float(b41):.3e}")

# Now step with FIXED w = 1 using b41*, then b42 = 1e30, then b43 = 2!
h_cur = h40
for step_idx, b_val in enumerate([b41, mp.mpf('1e30'), mp.mpf(2)], start=41):
    kc = len(h_cur)
    hn = []
    for i in range(kc + 1):
        u = [mp.mpf(0)]
        if i < kc: u = poly_add_mp(u, apply_Dj_mp(h_cur[i], b_val))
        if i >= 1: u = poly_add_mp(u, h_cur[i - 1])
        hn.append(u)
    h_cur = hn
    # Check Q_{k, q} at q = 0.25, 0.125!
    for qv in [0.5, 0.25, 0.125]:
        qm = mp.mpf(str(qv))
        ok_all = True
        for tc in [mp.mpc(t_probe, abs(t_probe)*mp.mpf('1e-4')), mp.mpc(0, abs(t_probe)*mp.mpf('1e-2')), mp.mpc('-0.5', '1e-6'), mp.mpc('10', '1e-3')]:
            cy = [(qm**(i*(i-1)//2)) * poly_eval_c_mp(h_cur[i], tc) for i in range(len(h_cur))]
            ok, _ = check_hb_pd(cy)
            if not ok:
                ok_all = False
                break
        if not ok_all or step_idx == 43:
            print(f"  step {step_idx} (b={float(b_val):.2e}), w=1, q={qv}: {'OK' if ok_all else 'FAIL!'}")
