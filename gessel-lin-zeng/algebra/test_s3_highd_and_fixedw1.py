"""
Two decisive diagnostics (< 15s total):
1. STAGE 3 ANATOMY:
   - Print the exact d=2 extremizer (a, b1, b2, t) at j=80, 160 on t >= 0 and its 6 j-poles!
   - Test Stage 3 drop on HIGHER DEGREE (d = j = 20, 40) triples (f0, f1, f2) in RZ(-1,0)
     where f_r are random GLZ-like interlacing triples (or adversarial perturbations of GLZ-like roots)!
     Does degree d = j KILL the -O(1/j) drop, or can d = j STILL drop by -O(1/j)?!
   - Test Strengthened Cubic Stage 3 at d = 3 (with c_in = 1 - 2/j): does it also drop by O(1/j)
     when a root of f2 approaches -1/j?
2. FIXED w_j = 1 TURAN RESILIENCE AT k = 50, 60:
   - For single-step D_2^k (where J_i(0) = (8/9) * (i+1)(k-i+1)/(i(k-i)) EXACTLY!!),
     at k = 60, i = 30: J_{30}(0) = (8/9) * (31/30)^2 = 0.9491 < 1 (5% below 1!).
     At k = 100, i = 50: J_{50}(0) = (8/9) * (51/50)^2 = 0.9248 < 1 (7.5% below 1!).
     Can ANY b_{k+1} >= 2 with FIXED w_{k+1} = 1 break U_Turan = (D_b h_i + h_{i-1})^2 - (D_b h_{i-1} + h_{i-2})(D_b h_{i+1} + h_i) < 0
     at t near 0^-?!
   WAIT A SECOND!! WE CAN DO THIS ANALYTICALLY FOR D_2^k (WHERE h_{k, i} = C(k, i) D_2^{k-i}(1)!!)!
   LET'S CHECK IT IN 2 SECONDS!!
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 60

# ============================================================
# 1A. INSPECT THE d=2 EXTREMIZER ON t >= 0 AT j = 80, 160
# ============================================================
print("=== 1A. EXACT d=2 EXTREMIZER ON t >= 0 ===")
def eval_J_d2(a, b1, b2, j, t):
    X  = b1*b2 + (b1+b2)*t + t*t
    Xp = (b1+b2) + 2*t
    M  = a + (2 + a*j)*t + (j + 1)*t*t
    Mp = (2 + a*j) + 2*(j + 1)*t
    Y  = 1 + 3*j*t + j*(j + 1)*t*t
    Yp = 3*j + 2*j*(j + 1)*t
    wXM = X*Mp - Xp*M
    wMY = M*Yp - Mp*Y
    wXY = X*Yp - Xp*Y
    return 4*wXM*wMY / (wXY**2), wXM, wMY, wXY

ts_check = np.concatenate([np.linspace(-5.0, -1.05, 25), np.linspace(-0.99, -0.01, 40), np.linspace(0.0, 5.0, 25), [20.0, 100.0]])
rng = np.random.default_rng(77)
for j in [40, 160]:
    worst = (1e9, None)
    for _ in range(8000):
        r = np.sort(np.exp(rng.uniform(-8, 0, 3)))
        if rng.random() < 0.5:
            r = 1.0 - r[::-1]
        b1, a, b2 = r[0], r[1], r[2]
        ok = True
        for tc in ts_check:
            J, w1, w2, w3 = eval_J_d2(a, b1, b2, j, tc)
            if w1 <= 0 or w2 <= 0 or w3 <= 0 or J < 1.0:
                ok = False
                break
        if not ok:
            continue
        for tc in [0.0, 1e-4, 1e-2, 0.1, 0.5, 2.0, 10.0, 100.0]:
            Ja = eval_J_d2(a, b1, b2, j, tc)[0]
            Jb = eval_J_d2(a, b1, b2, j + 1, tc)[0]
            val = (j**2) * np.log(Jb / Ja)
            if val < worst[0]:
                worst = (val, (a, b1, b2, tc, Ja, Jb))
    a, b1, b2, tc, Ja, Jb = worst[1]
    print(f"  j={j:3d}: j^2*dlnJ={worst[0]:+.4f}, b1={b1:.5f}, a={a:.5f}, b2={b2:.5f}, t={tc}, Ja={Ja:.5f}, Jb={Jb:.5f}, min_J_on_grid={min(eval_J_d2(a,b1,b2,j,t)[0] for t in ts_check):.5f}")

# ============================================================
# 1B. HIGH-DEGREE (d = j) TEST OF STAGE 3 ON t >= 0 AND (-1, 0)!
# ============================================================
print("\n=== 1B. CAN A HIGH-DEGREE (d = j) QUAD STILL DROP BY -O(1/j)? ===")
# Take a degree-d polynomial with (d-2) roots bunched near -1 (so they barely affect t ~ 0!),
# plus 2 roots at (-b2, -b1) near 0!
# Wait: if f2 = (t + 0.99)^{d-2}, f1 = (t + 0.99)^{d-2} (t + a), f0 = (t + 0.99)^{d-2} (t + b1)(t + b2),
# does D_j act on (t+0.99)^{d-2} * g(t) almost like D_{j + d - 2} on g(t) near t = 0?!
# CHECK!! Near t = 0: (t + 0.99)' / (t + 0.99) ~ 1, so:
#   D_j[(t+0.99)^m g(t)] = (t+0.99)^m [ t(1+t) g' + (1 + (j + m(1+t)/(t+0.99)) t) g ]
#                        ~ (t+0.99)^m * D_{j + m} g(t)!!
# AHA!! Let's test if taking b1, a, b2 ~ 1/(j+m) STILL reproduces the EXACT SAME -O(1/j) drop at d = j!!
def eval_J_high_d(m, a, b1, b2, j, t):
    # f2 = (t + r0)^m, f1 = (t + r0)^m (t + a), f0 = (t + r0)^m (t + b1)(t + b2) with r0 = 0.99
    # Wait: compute exact f0, D_j f1, D_j^2 f2 via poly arithmetic (fast in mpmath for m = j <= 60)!
    r0 = mp.mpf('0.95')
    # build (t + r0)^m
    base = [mp.mpf(1)]
    for _ in range(m):
        nb = [mp.mpf(0)] * (len(base) + 1)
        for idx, c in enumerate(base):
            nb[idx] += r0 * c
            nb[idx + 1] += c
        base = nb
    def mul_lin(p, c0):
        nb = [mp.mpf(0)] * (len(p) + 1)
        for idx, c in enumerate(p):
            nb[idx] += c0 * c
            nb[idx + 1] += c
        return nb
    def apply_Dj(p, jj):
        out = [mp.mpf(0)] * (len(p) + 1)
        for idx, c in enumerate(p):
            out[idx] += (1 + idx) * c
            out[idx + 1] += (jj + idx) * c
        return out
    def eval_pd(p, tv):
        v = mp.mpf(0); dv = mp.mpf(0)
        for c in reversed(p):
            dv = dv * tv + v
            v = v * tv + c
        return v, dv
    f2 = base
    f1 = mul_lin(base, mp.mpf(a))
    f0 = mul_lin(mul_lin(base, mp.mpf(b1)), mp.mpf(b2))
    Ma = apply_Dj(f1, j); Ya = apply_Dj(apply_Dj(f2, j), j)
    Mb = apply_Dj(f1, j + 1); Yb = apply_Dj(apply_Dj(f2, j + 1), j + 1)
    def W(u, v): return u[0]*v[1] - u[1]*v[0]
    def getJ(M, Y, tv):
        vX = eval_pd(f0, tv); vM = eval_pd(M, tv); vY = eval_pd(Y, tv)
        return float(4*W(vX, vM)*W(vM, vY)/(W(vX, vY)**2))
    # check min on grid for j and j+1
    grid = [mp.mpf(x) for x in [-2.0, -0.98, -0.9, -0.5, -0.1, -0.02, -0.005, -0.001, 0.0, 0.01, 0.1, 1.0, 10.0]]
    minJa = min(getJ(Ma, Ya, tv) for tv in grid)
    minJb = min(getJ(Mb, Yb, tv) for tv in grid)
    # specific t = 0.01 and t = -0.005
    Ja0 = getJ(Ma, Ya, mp.mpf(t))
    Jb0 = getJ(Mb, Yb, mp.mpf(t))
    return minJa, minJb, Ja0, Jb0

for j in [20, 40, 60]:
    # Scale (b1, a, b2) like 1/(2*j) since j_eff = j + m = 2*j!
    # In d=2 extremizer at j=40, let's use a = 3/(2j), b1 = 0.5/(2j), b2 = 10/(2j)
    je = 2 * j
    minJa, minJb, Ja0, Jb0 = eval_J_high_d(j, 2.5/je, 0.4/je, 15.0/je, j, 0.01)
    print(f"  d={j+2} (m={j}, j={j}): minJa={minJa:.5f}, minJb={minJb:.5f}, at t=0.01: Ja={Ja0:.5f}, Jb={Jb0:.5f}, j*(Jb-Ja)/Ja={j*(Jb0-Ja0)/Ja0:+.4f}, j^2*dlnJ={j*j*np.log(Jb0/Ja0):+.3f}")

# ============================================================
# 2. FIXED w = 1 TURAN BREAK TEST FOR (D_2 + y)^k AT k = 60, 100!
# ============================================================
print("\n=== 2. CAN FIXED w=1 BREAK U_Turan ON h_{k,i} = C(k,i) D_2^{k-i}(1) AT k = 60, 100? ===")
# Remember: for h_k = (D_2 + y)^k(1), h_{k,i}(t) = C(k, i) * D_2^{k-i}(1)!
# At t = 0: D_2^m(1)|_{t=0} = 1, and (D_2^m(1))'|_{t=0} = 2(2^m - 1)!!
# WAIT A SECOND!! CHECK THAT EXACT FORMULA!!
# For p_m(t) = D_2^m(1):
#   p_0 = 1 (p_0(0)=1, p_0'(0)=0)
#   p_m(0) = p_{m-1}(0) = 1
#   p_m'(0) = 2 p_{m-1}'(0) + 2 p_{m-1}(0) = 2 p_{m-1}'(0) + 2 => p_m'(0) = 2(2^m - 1)!!
# Let's check exact Turan of U_i(t) = D_b h_{k,i}(t) + h_{k, i-1}(t) at t = -eps (using exact p_m(t)!)
# Build p_m(t) = D_2^m(1) up to m = 65 in 60-digit mpmath (takes 0.1s!)
p = [[mp.mpf(1)]]
for m in range(1, 65):
    prev = p[-1]
    nxt = [mp.mpf(0)] * (len(prev) + 1)
    for idx, c in enumerate(prev):
        nxt[idx] += (1 + idx) * c
        nxt[idx + 1] += (2 + idx) * c
    p.append(nxt)

def eval_pd_mp(poly, tv):
    v = mp.mpf(0); dv = mp.mpf(0)
    for c in reversed(poly):
        dv = dv * tv + v
        v = v * tv + c
    return v, dv

for k in [40, 60]:
    i = k // 2
    # h_{k, r}(t) = C(k, r) * p[k - r](t)
    c_m1, c_0, c_1, c_2 = mp.binomial(k, i-2), mp.binomial(k, i-1), mp.binomial(k, i), mp.binomial(k, i+1)
    # slope of p[k-i] at 0 is 2(2^{k-i}-1) ~ 2^{k-i+1}. Probe t = -alpha * 2^{-(k-i)}!
    min_U_turan = mp.mpf('1e100')
    min_D_turan = mp.mpf('1e100')
    best_info = None
    for alpha in [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0]:
        t = -mp.mpf(str(alpha)) * (mp.mpf(2) ** (-(k - i)))
        # evaluate h_{k, i-2}, h_{k, i-1}, h_{k, i}, h_{k, i+1}
        vm1, dm1 = eval_pd_mp(p[k - (i - 2)], t); vm1 *= c_m1; dm1 *= c_m1
        v0, d0 = eval_pd_mp(p[k - (i - 1)], t); v0 *= c_0; d0 *= c_0
        v1, d1 = eval_pd_mp(p[k - i], t); v1 *= c_1; d1 *= c_1
        v2, d2 = eval_pd_mp(p[k - (i + 1)], t); v2 *= c_2; d2 *= c_2
        # D_b h_r = t(1+t) d_r + Phi * v_r, where Phi = 1 + b t <= 1 + 2 t (since b >= 2, t < 0)!
        # Minimize U_1^2 - U_0 U_2 over ALL Phi <= 1 + 2t (it is a quadratic in Phi!)
        # U_0 = (u0 + vm1) + Phi * v0
        # U_1 = (u1 + v0)  + Phi * v1
        # U_2 = (u2 + v1)  + Phi * v2
        u0, u1, u2 = t*(1+t)*d0, t*(1+t)*d1, t*(1+t)*d2
        # Pure D_b Turan:
        A_D = u1**2 - u0*u2
        B_D = 2*u1*v1 - (u0*v2 + u2*v0)
        C_D = v1**2 - v0*v2
        if C_D > 0:
            PhiD = min(1 + 2*t, -B_D/(2*C_D))
            valD = (A_D + PhiD*B_D + PhiD**2 * C_D) / (v1**2)
            if valD < min_D_turan: min_D_turan = valD
        # U = D_b h_i + h_{i-1} (FIXED w = 1!) Turan:
        q0, q1, q2 = u0 + vm1, u1 + v0, u2 + v1
        A_U = q1**2 - q0*q2
        B_U = 2*q1*v1 - (q0*v2 + q2*v0)
        C_U = v1**2 - v0*v2
        if C_U > 0:
            PhiU = min(1 + 2*t, -B_U/(2*C_U))
            valU = (A_U + PhiU*B_U + PhiU**2 * C_U) / (v1**2)
            if valU < min_U_turan:
                min_U_turan = valU
                best_info = (alpha, float(PhiU), float((1 - PhiU)/(-t)))
    print(f"  (D_2+y)^{k} (i={i}): min D_Turan/v1^2 = {float(min_D_turan):+.5e}, min U_Turan(w=1)/v1^2 = {float(min_U_turan):+.5e} at {best_info}")
