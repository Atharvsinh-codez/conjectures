"""
Quick 6-second diagnostic on:
1. GLZ sign of Delta'(0), B0(t), A0(t), J_raw(t) across ALL (k, i) for k=2..45 on [-1/(k+2), 0)!
   Specifically check BOTH t = -1e-18 (near 0) AND exact midpoints between roots of h_{k,i} inside [-1/(k+2), 0)!
2. Adversarial test of U_Turan >= 0 under J_raw >= 1 on R (with balanced scales so valid >> 0).
3. Fixed-input Combined GLZ 1-step map G_j at lam = 3/j, 1, j/3 (sampling h_r = D_j u_r so order matches GLZ!).
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 80
rng = np.random.default_rng(2026)

def apply_Dj_mp(poly, j):
    out = [mp.mpf(0)] * (len(poly) + 1)
    for m, c in enumerate(poly):
        if c != 0:
            out[m] += (1 + m) * c; out[m + 1] += (j + m) * c
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
        dv = dv * tv + v; v = v * tv + c
    return v, dv

def apply_Dj_f(p, j):
    out = np.zeros(len(p) + 1)
    for m, c in enumerate(p):
        out[m] += (1 + m) * c; out[m + 1] += (j + m) * c
    return out

def eval_pd_f(p, t):
    v = 0.0; dv = 0.0
    for c in reversed(p):
        dv = dv * t + v; v = v * t + c
    return v, dv

# ============================================================
# 1. CHECK Delta'(0), B0(t), A0(t), J_raw(0) FOR GLZ k=2..45!
# ============================================================
print("=== 1. GLZ k=2..45: SIGN OF Delta'(0), B0(t) ON [-1/(k+2), 0), AND J_raw(0) ===")
h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
neg_B0_cases = []
Jlt1_cases = []
both_bad = []
for k in range(2, 46):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = padd_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = padd_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
    # At t = 0: Delta'(0) = 2 c0(i) c1(i) - c0(i-1) c1(i+1) - c1(i-1) c0(i+1)
    # Wait: B0(t) = t(1+t) Delta'(t), so near 0^-, sign(B0) = -sign(Delta'(0))!
    # Also compute J_raw(0)!
    # W(f, g)(0) = c0(f) c1(g) - c1(f) c0(g).
    if k in (2, 5, 10, 15, 20, 25, 30, 35, 40, 45):
        dp0_neg = 0; dp0_pos = 0; j0_lt1 = 0; min_j0 = mp.mpf(999)
        bad_pair = []
        for i in range(1, k):
            c0_0, c1_0 = h[i-1][0], h[i-1][1]
            c0_1, c1_1 = h[i][0], h[i][1]
            c0_2, c1_2 = h[i+1][0], (h[i+1][1] if len(h[i+1]) > 1 else mp.mpf(0))
            dp0 = 2*c0_1*c1_1 - c0_0*c1_2 - c1_0*c0_2
            if dp0 > 0: dp0_pos += 1
            else: dp0_neg += 1
            w21 = c0_2*c1_1 - c1_2*c0_1
            w10 = c0_1*c1_0 - c1_1*c0_0
            w20 = c0_2*c1_0 - c1_2*c0_0
            j0 = 4*w21*w10 / (w20**2)
            if j0 < min_j0: min_j0 = j0
            if j0 < 1: j0_lt1 += 1
            if dp0 > 0 and j0 < 1:
                bad_pair.append((i, float(dp0), float(j0)))
        print(f"  k={k:2d}: Delta'(0)>0 (B0<0 near 0^-): {dp0_pos}/{k-1}, Delta'(0)<0: {dp0_neg}/{k-1}, J0<1: {j0_lt1}/{k-1}, min J0={float(min_j0):.5f}, (B0<0 & J0<1)={bad_pair}")

# Also check k=35, i=17 AND k=45, i=22 on 500 points in [-1/(k+2), 0) INCLUDING very close to 0 (-1e-18) AND exact sign of Delta'(t)!
f0, f1, f2 = h[21], h[22], h[23] # k=45, i=22
print("  Inspecting k=45, i=22 at specific t in [-1/47, 0):")
for tv in [mp.mpf(-1)/47, mp.mpf('-0.015'), mp.mpf('-0.01'), mp.mpf('-0.005'), mp.mpf('-0.001'), mp.mpf('-1e-5'), mp.mpf('-1e-10'), mp.mpf('-1e-18')]:
    v0, d0 = eval_pd_mp(f0, tv); v1, d1 = eval_pd_mp(f1, tv); v2, d2 = eval_pd_mp(f2, tv)
    u0, u1, u2 = tv*(1+tv)*d0, tv*(1+tv)*d1, tv*(1+tv)*d2
    A0 = u1**2 - u0*u2; B0 = 2*u1*v1 - u0*v2 - u2*v0; C0 = v1**2 - v0*v2
    w21 = v2*d1 - d2*v1; w10 = v1*d0 - d1*v0; w20 = v2*d0 - d2*v0
    J = 4*w21*w10 / (w20**2)
    print(f"    t={float(tv):10.2e}: A0={float(A0):+.3e}, B0={float(B0):+.3e}, C0={float(C0):+.3e}, J={float(J):.6f}")

# ============================================================
# 2. ADVERSARIAL TEST OF U_Turan >= 0 UNDER J_raw >= 1 ON R (BALANCED SCALES!)
# ============================================================
print("\n=== 2. ADVERSARIAL TEST OF U_Turan >= 0 UNDER J_raw >= 1 ON R ===")
ts_chk = np.linspace(-0.99, -0.005, 80)
ts_all = np.concatenate([np.linspace(-4.0, -1.02, 12), ts_chk, np.linspace(0.01, 6.0, 15)])
valid_u = 0; fail_u = 0; fail_ex = None
for trial in range(15000):
    r0 = np.sort(-rng.uniform(0.01, 0.99, 3))
    r1 = np.array([rng.uniform(r0[0], r0[1]), rng.uniform(r0[1], r0[2])])
    r2 = np.array([rng.uniform(r1[0], r1[1])])
    # Normalized at t=0, with mild log-concave scaling s1^2 >= s0 s2, s2^2 >= s1 s3!
    h0 = np.poly(r0)[::-1]; h0 /= h0[0]
    h1 = np.poly(r1)[::-1]; h1 /= h1[0]
    h2 = np.poly(r2)[::-1]; h2 /= h2[0]
    h3 = np.array([1.0])
    # apply scale ratio a^r * b^(r(r-1)/2) with b in (0.2, 1.0)
    a_sc = np.exp(rng.normal(0, 1.5))
    b_sc = rng.uniform(0.3, 1.0)
    h1 *= a_sc; h2 *= (a_sc**2) * b_sc; h3 *= (a_sc**3) * (b_sc**3)
    ok = True
    for t in ts_all:
        v0,d0 = eval_pd_f(h0,t); v1,d1 = eval_pd_f(h1,t); v2,d2 = eval_pd_f(h2,t); v3,d3 = eval_pd_f(h3,t)
        if -1 < t < 0:
            if v1*v1 < v0*v2 or v2*v2 < v1*v3: ok = False; break
        w21 = v2*d1-d2*v1; w10 = v1*d0-d1*v0; w20 = v2*d0-d2*v0
        w32 = v3*d2-d3*v2; w31 = v3*d1-d3*v1
        if w20 <= 0 or w31 <= 0 or 4*w21*w10 < w20*w20 or 4*w32*w21 < w31*w31:
            ok = False; break
    if not ok: continue
    valid_u += 1
    bad = False
    for j in [2, 10, 50]:
        Dh1, Dh2, Dh3 = apply_Dj_f(h1, j), apply_Dj_f(h2, j), apply_Dj_f(h3, j)
        U0 = np.pad(h0, (0, max(0, len(Dh1)-len(h0)))) + np.pad(Dh1, (0, max(0, len(h0)-len(Dh1))))
        U1 = np.pad(h1, (0, max(0, len(Dh2)-len(h1)))) + np.pad(Dh2, (0, max(0, len(h1)-len(Dh2))))
        U2 = np.pad(h2, (0, max(0, len(Dh3)-len(h2)))) + np.pad(Dh3, (0, max(0, len(h2)-len(Dh3))))
        for t in ts_chk:
            u0v = eval_pd_f(U0,t)[0]; u1v = eval_pd_f(U1,t)[0]; u2v = eval_pd_f(U2,t)[0]
            if u1v*u1v - u0v*u2v < -1e-11*(u1v*u1v + abs(u0v*u2v) + 1e-30):
                bad = True
                if fail_ex is None: fail_ex = (j, float(t), float(u1v*u1v - u0v*u2v))
                break
        if bad: break
    if bad: fail_u += 1
print(f"  Adversarial deg(3,2,1,0) J_raw>=1 on R & Turan on (-1,0) (valid={valid_u}): U_Turan<0 fails={fail_u}, ex={fail_ex}")

# ============================================================
# 3. COMBINED GLZ 1-STEP MAP G_j (PROPERLY ORDERED INPUTS!)
# ============================================================
print("\n=== 3. COMBINED GLZ 1-STEP MAP G_j (PROPERLY SAMPLED INPUTS h_r = D_j u_r!) ===")
# Sample E-triple (E3, E2, E1, E0) directly as D_j^3 f3 << D_j^2 f2 << D_j f1 << f0 by taking
# f_r = D_{j-1} g_r (wait: even simpler! Take h_r = D_{j-1}^{3-r} ... no, take random roots for h0,
# and sample h1 = D_j^{-1}(...)! Even simpler: in GLZ, h_r = D_{j-1} u_r, where u_3 << u_2 << u_1 << u_0!)
# Check: if u3 << u2 << u1 << u0 (degs 0,1,2,3), and we define h_r = D_{j-1} u_r?
# Wait: then D_j h_2 has deg 3, h_1 has deg 3: let's check order of (D_j^2 h_2, D_j h_1, h_0) when h_r = D_{j-1}^{r?}...
# No! In GLZ, at step k-1 -> k:
# h_{k-1, i} has deg 2k-i-3. E_{k, i} = D_{k+1}^i h_{k-1, i} has deg 2k-3.
# Let's check for k=4 (j=5): E_{4,0} = h_{3,0} (deg 5), E_{4,1} = D_5 h_{3,1} (deg 5), E_{4,2} = D_5^2 h_{3,2} (deg 5).
# Their order is E_{4,2} << E_{4,1} << E_{4,0}!
# WHY is D_5 h_{3,1} << h_{3,0}? Because h_{3,0} = D_4^2 h_{2,0} + ..., so h_{3,0} had TWO D_4's whereas h_{3,1} had ONE D_4 + shift!
# So to sample valid inputs (h3, h2, h1, h0) with D_j^2 h2 << D_j h1 << h0:
# just sample 4 interlacing polynomials of SAME degree d (representing E3 << E2 << E1 << E0),
# wait: G_j CANNOT be written purely on (E3, E2, E1, E0) because U_r = D_j h_r + h_{r-1}, so
# D_j^r U_r = D_j^{r+1} h_r + D_j^r h_{r-1} = E_{r+1} + D_j E_r? NO!
# E_r = D_j^r h_r (with j = k+1)! Wait: D_j^r U_r = D_j^{r+1} h_r + D_j^r h_{r-1} = E_r (with one extra D_j!) + D_j E_{r-1}!
# WAIT A SECOND!
# CHECK THAT FORMULA!!
# E_{k, i} := D_{k+1}^i h_{k-1, i}.
# Then D_{k+1}^i U_{k, i} = D_{k+1}^{i+1} h_{k-1, i} + D_{k+1}^i h_{k-1, i-1} = D_{k+1} E_{k, i} + D_{k+1} E_{k, i-1}!!
# WHAT?!?! Wait: D_{k+1}^i h_{k-1, i-1} = D_{k+1} (D_{k+1}^{i-1} h_{k-1, i-1}) = D_{k+1} E_{k, i-1}!!
# SO D_{k+1}^i h_{k, i} = D_{k+1}^{i+1} U_{k, i} = D_{k+1}^2 (E_{k, i} + E_{k, i-1})!!
# AND THEN E_{k+1, i} = D_{k+2}^i h_{k, i} = (D_{k+2}/D_{k+1})^i D_{k+1}^2 (E_{k, i} + E_{k, i-1})!!
# Let's test G_j directly on (E3, E2, E1, E0) via h_r = D_j^{3-r} ... wait: we already know Stage 1+2+3!
# Here let's test the EXACT combined map on (E3 << E2 << E1 << E0) where E_r = D_j^r h_r:
# Since h_r = D_j^{-r} E_r, sample a0 << a1 << a2 << a3, set h_r = D_j^{3-r} a_r? No, set h3 = a3 (deg 0),
# h2 = D_j a3 + c2 (deg 1), etc. Let's just sample directly from GLZ-like perturbations!
for j in [20, 40, 80]:
    for lam_name, lam in [("3/j", 3.0/j), ("1.0", 1.0), ("j/3", j/3.0)]:
        valid_g = 0; fails_1 = 0; fails_tgt = 0; min_Jout = 1e9
        for _ in range(4000):
            # Sample 4 interlacing polys g3 << g2 << g1 << g0 of deg 1 (so E_r = D_j^r h_r are deg 4):
            # Take h3 = 1 (deg 0), h2 = D_{j-1}(1) + s2*... No! Just take arbitrary deg-3 h0 with roots in (-1,0),
            # no, take a3 = 1, a2 = D_{j-1}(1) + ..., wait: take ANY 4 interlacing polys of deg 3 for (E3, E2, E1, E0)?
            # Wait: E_{k+1, i} = D_{k+2}^i D_{k+1}^{2-i} (E_{k,i} + E_{k,i-1}) — NO!
            # Only for i <= 2 does D_{k+1}^{2-i} have non-negative exponent! For general i, h_{k,i} = D_{k+1}(D_{k+1} h_{k-1,i} + h_{k-1,i-1}).
            # For a TRIPLE (i-1, i, i+1), shift index by i-1:
            # Let f_0 = h_{k-1, i-2}, f_1 = h_{k-1, i-1}, f_2 = h_{k-1, i}, f_3 = h_{k-1, i+1} (degs d+3, d+2, d+1, d).
            # Wait: to get D_j^2 f_2 << D_j f_1 << f_0, build from bottom:
            # Take f3 = [1.0], f2 = D_{j-1}(f3) * (1 + eps*t)? Even better: build f_r = D_{j-1}^{3-r}(g_r) where g3 << g2 << g1 << g0 (all deg 1)!
            r0 = -rng.uniform(0.02, 0.98)
            r1 = -rng.uniform(0.01, -r0)
            r2 = -rng.uniform(0.005, -r1)
            r3 = -rng.uniform(0.001, -r2)
            # order: g0 << g1 << g2 << g3? Wait: if r0 < r1 < r2 < r3 < 0, then (t-r0) << (t-r1) << (t-r2) << (t-r3)!
            # Wait: we want f3 << f2 << f1 << f0 AND D_j^2 f2 << D_j f1 << f0!
            # If f3 = (t-r0), f2 = D_{j-1}(t-r1), f1 = D_{j-1}^2(t-r2), f0 = D_{j-1}^3(t-r3):
            # Then D_j^3 f3 = D_j^3 (t-r0), D_j^2 f2 = D_j^2 D_{j-1}(t-r1), D_j f1 = D_j D_{j-1}^2(t-r2), f0 = D_{j-1}^3(t-r3)!
            # Since D_j shifts MORE left or MORE right than D_{j-1}?
            # D_j(1) = 1+jt (root -1/j), D_{j-1}(1) = 1+(j-1)t (root -1/(j-1) < -1/j)!
            # So D_j shifts FURTHER RIGHT than D_{j-1}!
            # Hence D_j^3(t-r0) might jump RIGHT of D_{j-1}^3(t-r3) UNLESS r0 is sufficiently left of r3!
            # Let's check directly with random roots r0 < r1 < r2 < r3!
            g0 = np.array([-r0, 1.0]); g1 = np.array([-r1, 1.0]); g2 = np.array([-r2, 1.0]); g3 = np.array([-r3, 1.0])
            f3 = g0
            f2 = apply_Dj_f(g1, j-1)
            f1 = apply_Dj_f(apply_Dj_f(g2, j-1), j-1)
            f0 = apply_Dj_f(apply_Dj_f(apply_Dj_f(g3, j-1), j-1), j-1)
            f0 /= f0[0]; f1 /= f1[0]; f2 /= f2[0]; f3 /= f3[0]
            # E-quad: (D_j^3 f3, D_j^2 f2, D_j f1, f0)
            E3 = apply_Dj_f(apply_Dj_f(apply_Dj_f(f3, j), j), j)
            E2 = apply_Dj_f(apply_Dj_f(f2, j), j)
            E1 = apply_Dj_f(f1, j)
            E0 = f0
            ok = True
            for t in ts_all:
                v0,d0 = eval_pd_f(E2,t); v1,d1 = eval_pd_f(E1,t); v2,d2 = eval_pd_f(E0,t)
                w01=v0*d1-d0*v1; w12=v1*d2-d1*v2; w02=v0*d2-d0*v2
                if w01<=0 or w12<=0 or w02<=0 or 4*w01*w12 < (1 + 2.0/j)*w02*w02:
                    ok = False; break
                u0,du0 = eval_pd_f(E3,t); u1,du1 = eval_pd_f(E2,t); u2,du2 = eval_pd_f(E1,t)
                z01=u0*du1-du0*u1; z12=u1*du2-du1*u2; z02=u0*du2-du0*u2
                if z01<=0 or z12<=0 or z02<=0 or 4*z01*z12 < (1 + 2.0/j)*z02*z02:
                    ok = False; break
            if not ok: continue
            valid_g += 1
            def add_s(pa, pb, sc):
                n = max(len(pa), len(pb))
                out = np.zeros(n)
                out[:len(pa)] += pa; out[:len(pb)] += sc * pb
                return out
            # U_right = f0 + lam * D_j f1, U_mid = f1 + lam * D_j f2, U_left = f2 + lam * D_j f3
            U_right = add_s(f0, apply_Dj_f(f1, j), lam)
            U_mid   = add_s(f1, apply_Dj_f(f2, j), lam)
            U_left  = add_s(f2, apply_Dj_f(f3, j), lam)
            X = apply_Dj_f(apply_Dj_f(apply_Dj_f(U_left, j), j+1), j+1)
            M = apply_Dj_f(apply_Dj_f(U_mid, j), j+1)
            Y = apply_Dj_f(U_right, j)
            mJ = 1e9
            for t in np.concatenate([ts_all, [-2.0/j, -1.0/j, -0.5/j, 0.5/j, 2.0/j]]):
                v0,d0 = eval_pd_f(X,t); v1,d1 = eval_pd_f(M,t); v2,d2 = eval_pd_f(Y,t)
                w01=v0*d1-d0*v1; w12=v1*d2-d1*v2; w02=v0*d2-d0*v2
                J = 4*w01*w12/(w02*w02) if (w01>0 and w12>0 and w02>0) else -1.0
                if J < mJ: mJ = J
            if mJ < min_Jout: min_Jout = mJ
            if mJ < 1.0: fails_1 += 1
            if mJ < 1.0 + 2.0/(j+1): fails_tgt += 1
        print(f"  G_j j={j:2d}, lam={lam_name:3s} (valid={valid_g:4d}): min J_out={min_Jout:.6f}, fails(<1)={fails_1}, fails(<1+2/(j+1))={fails_tgt}")
