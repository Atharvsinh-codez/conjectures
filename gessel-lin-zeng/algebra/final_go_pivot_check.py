"""
Final 5-second Go/Pivot diagnostic before reporting to @coordinator (#33):
1. Adversarial test of U_Turan >= 0 under J_raw >= 1 on R for deg(4,3,2,1) (where all Wronskians exist).
2. Combined GLZ 1-step map G_j at GLZ boundary scale lam = 3/j (and 1, j/3) on perturbed GLZ inputs!
3. Check t -> inf and t = 0 behavior of B_k^(E)(t, y) = sum_{i=0}^{k-1} E_{k,i}(t) y^i (where E_{k,i} = D_{k+1}^i h_{k-1,i}).
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 80
rng = np.random.default_rng(2026)

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

# ============================================================
# 1. ADVERSARIAL U_Turan >= 0 UNDER J_raw >= 1 ON R, deg(4,3,2,1)
# ============================================================
print("=== 1. ADVERSARIAL U_Turan >= 0 UNDER J_raw >= 1 ON R, deg(4,3,2,1) ===")
ts_chk = np.linspace(-0.99, -0.005, 60)
ts_all = np.concatenate([np.linspace(-4.0, -1.05, 10), ts_chk, np.linspace(0.01, 5.0, 12)])
valid_u = 0; fail_u = 0; neg_t3 = 0; fail_ex = None
for trial in range(12000):
    r0 = np.sort(-rng.uniform(0.01, 0.99, 4))
    r1 = np.array([rng.uniform(r0[m], r0[m+1]) for m in range(3)])
    r2 = np.array([rng.uniform(r1[m], r1[m+1]) for m in range(2)])
    r3 = np.array([rng.uniform(r2[0], r2[1])])
    h0 = np.poly(r0)[::-1]; h0 /= h0[0]
    h1 = np.poly(r1)[::-1]; h1 /= h1[0]
    h2 = np.poly(r2)[::-1]; h2 /= h2[0]
    h3 = np.poly(r3)[::-1]; h3 /= h3[0]
    # Extreme skew scales a^r * b^(r(r-1)/2)
    a_sc = np.exp(rng.normal(0, 2.0))
    b_sc = rng.uniform(0.15, 1.0)
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
    bad = False; t3_neg = False
    for j in [2, 10, 50]:
        Dh1, Dh2, Dh3 = apply_Dj_f(h1, j), apply_Dj_f(h2, j), apply_Dj_f(h3, j)
        U0 = np.pad(h0, (0, max(0, len(Dh1)-len(h0)))) + np.pad(Dh1, (0, max(0, len(h0)-len(Dh1))))
        U1 = np.pad(h1, (0, max(0, len(Dh2)-len(h1)))) + np.pad(Dh2, (0, max(0, len(h1)-len(Dh2))))
        U2 = np.pad(h2, (0, max(0, len(Dh3)-len(h2)))) + np.pad(Dh3, (0, max(0, len(h2)-len(Dh3))))
        for t in ts_chk:
            a0=eval_pd_f(h0,t)[0]; a1=eval_pd_f(h1,t)[0]; a2=eval_pd_f(h2,t)[0]
            b1=eval_pd_f(Dh1,t)[0]; b2=eval_pd_f(Dh2,t)[0]; b3=eval_pd_f(Dh3,t)[0]
            if 2*a1*b2 - a0*b3 - a2*b1 < -1e-11: t3_neg = True
            u0v = eval_pd_f(U0,t)[0]; u1v = eval_pd_f(U1,t)[0]; u2v = eval_pd_f(U2,t)[0]
            if u1v*u1v - u0v*u2v < -1e-11*(u1v*u1v + abs(u0v*u2v) + 1e-30):
                bad = True
                if fail_ex is None: fail_ex = (j, float(t), float(u1v*u1v - u0v*u2v), float(a_sc), float(b_sc))
                break
        if bad: break
    if t3_neg: neg_t3 += 1
    if bad: fail_u += 1
print(f"  deg(4,3,2,1) (valid={valid_u}): T3<0={neg_t3}, U_Turan<0 fails={fail_u}, ex={fail_ex}")

# ============================================================
# 2. COMBINED GLZ 1-STEP MAP G_j ON VALID 4-TUPLES AT lam = 3/j, 1, j/3!
# ============================================================
print("\n=== 2. COMBINED GLZ 1-STEP MAP G_j AT lam = 3/j, 1, j/3 ===")
# Build valid (f3, f2, f1, f0) with degs (1,2,3,4) where f_r = D_j^r? NO:
# Start from ANY 4 interlacing polynomials E3 << E2 << E1 << E0 of deg 4? No, we need f_r so D_j^r f_r = E_r!
# Easy exact construction of f_r with D_j^3 f3 << D_j^2 f2 << D_j f1 << f0:
# Take arbitrary deg-1 polynomials g3 << g2 << g1 << g0, and set:
#   f3 = g3,  f2 = D_j(g2) + eps*g3? No, just f3 = g3, f2 = D_{j+1}(g2)? Wait:
#   D_j^3 f3 = D_j^3 g3. If we set f2 = D_j g2 (same j!), then D_j^2 f2 = D_j^3 g2!
#   GENIUS!! If f3 = g3, f2 = D_j g2, f1 = D_j^2 g1, f0 = D_j^3 g0 (ALL with D_j!),
#   then (E3, E2, E1, E0) = (D_j^3 g3, D_j^3 g2, D_j^3 g1, D_j^3 g0)!
#   And by Zero-Loss D_j Preservation, if (g3, g2, g1, g0) has J >= 1 + c/j, SO DOES (E3, E2, E1, E0)!!
#   And we can ALSO perturb f_r slightly (or vary j in f_r to j-delta)!
for j in [20, 40, 80]:
    for lam_name, lam in [("3/j", 3.0/j), ("1.0", 1.0), ("j/3", j/3.0)]:
        valid_g = 0; fails_1 = 0; fails_tgt = 0; min_Jout = 1e9; worst_drop = 0.0
        for _ in range(3000):
            # Take deg-2 interlacing g3 << g2 << g1 << g0 with roots in (-1, 0)
            r0 = np.sort(-rng.uniform(0.01, 0.99, 2))
            # Pick 4 interlaced deg-2 polys: r3[0] < r2[0] < r1[0] < r0[0] < r3[1] < r2[1] < r1[1] < r0[1]
            # Wait: for g3 << g2 << g1 << g0 all of deg 2: need r3_1 < r2_1 < r1_1 < r0_1 < r3_2 < r2_2 < r1_2 < r0_2
            pts = np.sort(rng.uniform(0.005, 0.995, 8))
            if rng.random() < 0.5:
                pts = np.sort(np.exp(rng.uniform(-np.log(j)-1.0, 0, 8)))
            g3 = np.poly([-pts[7], -pts[3]])[::-1]; g3 /= g3[0]
            g2 = np.poly([-pts[6], -pts[2]])[::-1]; g2 /= g2[0]
            g1 = np.poly([-pts[5], -pts[1]])[::-1]; g1 /= g1[0]
            g0 = np.poly([-pts[4], -pts[0]])[::-1]; g0 /= g0[0]
            # Check J(g3, g2, g1) and J(g2, g1, g0) on ts_all (which equals J(E) >= J(g) by zero-loss D_j!)
            # Actually compute J_in directly on (E3, E2, E1) and (E2, E1, E0)!
            f3 = g3
            f2 = apply_Dj_f(g2, j)
            f1 = apply_Dj_f(apply_Dj_f(g1, j), j)
            f0 = apply_Dj_f(apply_Dj_f(apply_Dj_f(g0, j), j), j)
            E3 = apply_Dj_f(apply_Dj_f(apply_Dj_f(f3, j), j), j)
            E2 = apply_Dj_f(apply_Dj_f(f2, j), j)
            E1 = apply_Dj_f(f1, j)
            E0 = f0
            min_Jin = 1e9
            ok = True
            for t in ts_all:
                v0,d0 = eval_pd_f(E2,t); v1,d1 = eval_pd_f(E1,t); v2,d2 = eval_pd_f(E0,t)
                w01=v0*d1-d0*v1; w12=v1*d2-d1*v2; w02=v0*d2-d0*v2
                if w01<=0 or w12<=0 or w02<=0: ok = False; break
                J1 = 4*w01*w12/(w02*w02)
                u0,du0 = eval_pd_f(E3,t); u1,du1 = eval_pd_f(E2,t); u2,du2 = eval_pd_f(E1,t)
                z01=u0*du1-du0*u1; z12=u1*du2-du1*u2; z02=u0*du2-du0*u2
                if z01<=0 or z12<=0 or z02<=0: ok = False; break
                J2 = 4*z01*z12/(z02*z02)
                if min(J1, J2) < min_Jin: min_Jin = min(J1, J2)
            if not ok or min_Jin < 1.0 + 2.0/j: continue
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
            # Next E-triple: (D_{j+1}^2 D_j U_left, D_{j+1} D_j U_mid, D_j U_right)
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
            if min_Jin - mJ > worst_drop: worst_drop = min_Jin - mJ
        print(f"  G_j j={j:2d}, lam={lam_name:3s} (valid={valid_g:4d}): min J_out={min_Jout:.6f}, max(Jin-Jout)={worst_drop:.5f}, fails(<1)={fails_1}, fails(<1+2/(j+1))={fails_tgt}")

# ============================================================
# 3. LEADING-COEFF (t -> inf) CHECK ON B_k^(E)(t, y) = sum E_{k,i}(t) y^i
# ============================================================
print("\n=== 3. B_k^(E)(t, y) = sum_{i=0}^{k-1} (D_{k+1}^i h_{k-1, i}(t)) y^i AT t -> inf ===")
h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
for k in range(2, 41):
    # Compute E_{k, i} = D_{k+1}^i h_{k-1, i}
    if k in (5, 10, 20, 30, 40):
        j = k + 1
        lcs = []
        for i in range(k):
            p = h[i]
            for _ in range(i):
                p = apply_Dj_mp(p, j)
            lcs.append(p[-1])
        # Check Newton ratios lcs[i]^2 / (lcs[i-1]*lcs[i+1]) vs (i+1)(k-i)/(i(k-i-1))
        # Wait: for sum_{i=0}^{k-1} c_i y^i to be real-rooted, need c_i^2 / (c_{i-1} c_{i+1}) >= (i+1)(k-i)/(i(k-i-1))?
        # Specifically for deg m = k-1, Newton inequality is c_i^2/(c_{i-1}c_{i+1}) >= (i+1)(m-i+1)/(i(m-i))!
        m = k - 1
        min_ratio = mp.mpf(999)
        min_newton = mp.mpf(999)
        for i in range(1, m):
            r = (lcs[i]**2) / (lcs[i-1] * lcs[i+1])
            nb = mp.mpf((i+1)*(m-i+1)) / mp.mpf(i*(m-i))
            if r < min_ratio: min_ratio = r
            if r / nb < min_newton: min_newton = r / nb
        print(f"  k={k:2d}: t->inf lc min(c_i^2/(c_{{i-1}}c_{{i+1}})) = {float(min_ratio):.6f}, min(ratio/Newton) = {float(min_newton):.6f}")
    # advance h to k
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = padd_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = padd_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
