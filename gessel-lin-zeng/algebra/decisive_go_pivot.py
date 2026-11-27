"""
Decisive Go/Pivot Check (< 15s total):
1. TWO-REGIME SPLIT AT t = -1/(k+2) (#33 item 1):
   - On GLZ k = 35, 45: check whether B_0(t) = t(1+t) Delta'(t) or A_0(t) EVER goes negative
     inside [-1/(k+2), 0)! (Check fine grid inside [-1/(k+2), 0)!)
   - Check whether plain Turan + J_raw >= 1 on ALL of R EVEN implies T3 >= 0 on (-1, 0) for
     arbitrary 4-tuples (h_{i-2}, h_{i-1}, h_i, h_{i+1})!
     (Remember T3 = 2 h_{i-1}(D_j h_i) - h_{i-2}(D_j h_{i+1}) - h_i(D_j h_{i-1})!)
2. COMBINED GLZ 1-STEP MAP G_j AT ACTUAL GLZ BOUNDARY SCALE lambda = 3/j (#33 item 2):
   - Test G_j(h0, h1, h2, h3) = (D_j(h0 + lam D_j h1), D_{j+1} D_j(h1 + lam D_j h2), D_{j+1}^2 D_j(h2 + lam D_j h3))
     at lam = 3/j (GLZ's exact scale at i = k - 3!) on d=3 tuples with J_in >= 1 + 2/j!
3. FIXED-k DIFFERENTIAL RECURRENCE CHECK (#33 item 2):
   - Check in exact SymPy whether A_{k, i}(x) at k = 6 satisfies ANY linear relation among
     {A_{k, i}, A'_{k, i}, A_{k, i-1}, A'_{k, i-1}, A_{k, i+1}, A'_{k, i+1}} with degree <= 2 poly coeffs in x!
"""
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 80
rng = np.random.default_rng(2026)

# ============================================================
# 1A. CHECK B_0(t) AND A_0(t) ON [-1/(k+2), 0) FOR GLZ k = 35!
# ============================================================
print("=== 1A. TWO-REGIME SPLIT: DO A_0(t) OR B_0(t) GO NEGATIVE INSIDE [-1/(k+2), 0) ON GLZ k=35? ===")
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
        dv = dv * tv + v
        v = v * tv + c
    return v, dv

h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
for k in range(2, 36):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = padd_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = padd_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn

# At k = 35, i = 17, scan 300 log-spaced points inside [-1/37, -1e-18]!
f0, f1, f2 = h[16], h[17], h[18]
neg_A0 = 0; neg_B0 = 0; neg_B0_and_Jlt1 = 0; jlt1_cnt = 0
min_B0_ex = None; min_A0_ex = None
for m in range(300):
    # log-spaced from -1/37 to -1e-18
    t = -(mp.mpf(1)/37) * mp.power(mp.mpf('1e-16'), mp.mpf(m)/299)
    v0, d0 = eval_pd_mp(f0, t); v1, d1 = eval_pd_mp(f1, t); v2, d2 = eval_pd_mp(f2, t)
    u0, u1, u2 = t*(1+t)*d0, t*(1+t)*d1, t*(1+t)*d2
    A0 = u1**2 - u0*u2
    B0 = 2*u1*v1 - (u0*v2 + u2*v0)
    C0 = v1**2 - v0*v2
    w21 = v2*d1 - d2*v1; w10 = v1*d0 - d1*v0; w20 = v2*d0 - d2*v0
    J = 4*w21*w10 / (w20**2)
    if A0 < 0:
        neg_A0 += 1
        if min_A0_ex is None: min_A0_ex = (float(t), float(J))
    if B0 < 0:
        neg_B0 += 1
        if min_B0_ex is None: min_B0_ex = (float(t), float(J))
    if J < 1:
        jlt1_cnt += 1
        if B0 < 0: neg_B0_and_Jlt1 += 1
print(f"  GLZ k=35, i=17 on [-1/37, -1e-18] (300 pts): neg_A0={neg_A0} (ex={min_A0_ex}), neg_B0={neg_B0} (ex={min_B0_ex}), J<1:{jlt1_cnt}, (B0<0 & J<1)={neg_B0_and_Jlt1}")

# ============================================================
# 1B. DOES PLAIN TURAN + J_raw >= 1 ON R IMPLY T3 >= 0 (CROSS TERM OF U-TURAN)?!
# ============================================================
print("\n=== 1B. CAN T3 (CROSS TERM OF U-TURAN) GO NEGATIVE EVEN WHEN J_raw >= 1 ON ALL R? ===")
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

# Generate random 4-tuples (h0, h1, h2, h3) of degs (4, 3, 2, 1) with h3 << h2 << h1 << h0 in (-1,0)
# satisfying J_raw(h0, h1, h2) >= 1 AND J_raw(h1, h2, h3) >= 1 on R AND plain Turan on (-1,0)!
# Check if T3 = 2 h1 (D_j h2) - h0 (D_j h3) - h2 (D_j h1) EVER goes negative on (-1, 0),
# AND if U_Turan = U1^2 - U0 U2 (where U_r = D_j h_{r+1} + h_r) EVER goes negative on (-1, 0)!
ts_chk = np.linspace(-0.99, -0.01, 60)
ts_all = np.concatenate([np.linspace(-3.0, -1.05, 10), ts_chk, np.linspace(0.01, 5.0, 15)])
valid_1b = 0; neg_T3_cnt = 0; neg_U_cnt = 0
for _ in range(5000):
    r0 = np.sort(-rng.uniform(0.02, 0.98, 4))
    r1 = np.array([rng.uniform(r0[m], r0[m+1]) for m in range(3)])
    r2 = np.array([rng.uniform(r1[m], r1[m+1]) for m in range(2)])
    r3 = np.array([rng.uniform(r2[0], r2[1])])
    # random positive scales
    s_fac = np.exp(rng.normal(0, 0.8, 4))
    h0 = s_fac[0]*np.poly(r0)[::-1]
    h1 = s_fac[1]*np.poly(r1)[::-1]
    h2 = s_fac[2]*np.poly(r2)[::-1]
    h3 = s_fac[3]*np.poly(r3)[::-1]
    # Check plain Turan on (-1,0) AND J_raw >= 1 on ts_all for BOTH (h0,h1,h2) and (h1,h2,h3)
    ok = True
    for t in ts_all:
        v0,d0 = eval_pd_f(h0,t); v1,d1 = eval_pd_f(h1,t); v2,d2 = eval_pd_f(h2,t); v3,d3 = eval_pd_f(h3,t)
        if -1 < t < 0:
            if v1*v1 < v0*v2 or v2*v2 < v1*v3:
                ok = False; break
        w21 = v2*d1-d2*v1; w10 = v1*d0-d1*v0; w20 = v2*d0-d2*v0
        w32 = v3*d2-d3*v2; w31 = v3*d1-d3*v1
        if w20 == 0 or w31 == 0 or 4*w21*w10 < w20*w20 or 4*w32*w21 < w31*w31:
            ok = False; break
    if not ok: continue
    valid_1b += 1
    j = 6
    Dh1, Dh2, Dh3 = apply_Dj_f(h1, j), apply_Dj_f(h2, j), apply_Dj_f(h3, j)
    U0 = np.pad(h0, (0, max(0, len(Dh1)-len(h0)))) + np.pad(Dh1, (0, max(0, len(h0)-len(Dh1))))
    U1 = np.pad(h1, (0, max(0, len(Dh2)-len(h1)))) + np.pad(Dh2, (0, max(0, len(h1)-len(Dh2))))
    U2 = np.pad(h2, (0, max(0, len(Dh3)-len(h2)))) + np.pad(Dh3, (0, max(0, len(h2)-len(Dh3))))
    t3_bad = False; u_bad = False
    for t in ts_chk:
        a0 = eval_pd_f(h0,t)[0]; a1 = eval_pd_f(h1,t)[0]; a2 = eval_pd_f(h2,t)[0]
        b1 = eval_pd_f(Dh1,t)[0]; b2 = eval_pd_f(Dh2,t)[0]; b3 = eval_pd_f(Dh3,t)[0]
        T3 = 2*a1*b2 - a0*b3 - a2*b1
        u0v = eval_pd_f(U0,t)[0]; u1v = eval_pd_f(U1,t)[0]; u2v = eval_pd_f(U2,t)[0]
        if T3 < -1e-12: t3_bad = True
        if u1v*u1v - u0v*u2v < -1e-12: u_bad = True
    if t3_bad: neg_T3_cnt += 1
    if u_bad: neg_U_cnt += 1
print(f"  Random 4-tuples with J_raw>=1 on R & Turan on (-1,0) (valid={valid_1b}): T3<0 fails={neg_T3_cnt}, U_Turan<0 fails={neg_U_cnt}")

# ============================================================
# 2. COMBINED GLZ MAP G_j AT ACTUAL GLZ BOUNDARY SCALE lambda = 3/j!
# ============================================================
print("\n=== 2. COMBINED GLZ 1-STEP MAP G_j AT ACTUAL GLZ SCALE lambda = 3/j ===")
for j in [20, 40, 80]:
    lam = 3.0 / j
    min_Jout = 1e9
    valid_g = 0
    fails_1 = 0
    fails_target = 0
    for _ in range(4000):
        # d = (3, 2, 1, 0): h3 = 1, h2 = t+a, h1 = (t+b1)(t+b2), h0 = (t+c1)(t+c2)(t+c3)
        # Put roots near -1/j to probe adversarial corner!
        if rng.random() < 0.5:
            pts = np.sort(np.exp(rng.uniform(-np.log(j)-1.5, 0, 6)))
        else:
            pts = np.sort(rng.uniform(0.005, 0.995, 6))
        c1, b1, a, b2, c2, c3 = pts[0], pts[1], pts[2], pts[3], pts[4], pts[5]
        # Wait: ensure h3 << h2 << h1 << h0: c1 < b1 < c2 < b2 < c3, and b1 < a < b2!
        # Let's sort properly:
        c = np.sort(rng.uniform(0.001, 0.999, 3)) if rng.random() < 0.5 else np.sort(np.exp(rng.uniform(-np.log(j)-1.0, 0, 3)))
        b = np.array([rng.uniform(c[0], c[1]), rng.uniform(c[1], c[2])])
        a_r = rng.uniform(b[0], b[1])
        h3 = np.array([1.0])
        h2 = np.poly([-a_r])[::-1]
        h1 = np.poly(-b)[::-1]
        h0 = np.poly(-c)[::-1]
        # Scale h_r so that h_r + lam * D_j h_{r+1} has ratio lam (normalize at t=0: h_r(0) = 1!)
        h0 = h0 / h0[0]; h1 = h1 / h1[0]; h2 = h2 / h2[0]; h3 = h3 / h3[0]
        # Check input quads (h0, D_j h1, D_j^2 h2) and (h1, D_j h2, D_j^2 h3) have J >= 1 + 2/j on ts_all!
        Dj_h1 = apply_Dj_f(h1, j); Dj2_h2 = apply_Dj_f(apply_Dj_f(h2, j), j)
        Dj_h2 = apply_Dj_f(h2, j); Dj2_h3 = apply_Dj_f(apply_Dj_f(h3, j), j)
        ok = True
        for t in ts_all:
            v0,d0 = eval_pd_f(h0,t); v1,d1 = eval_pd_f(Dj_h1,t); v2,d2 = eval_pd_f(Dj2_h2,t)
            w01=v0*d1-d0*v1; w12=v1*d2-d1*v2; w02=v0*d2-d0*v2
            if w01<=0 or w12<=0 or w02<=0 or 4*w01*w12 < (1 + 2.0/j)*w02*w02:
                ok = False; break
            u0,du0 = eval_pd_f(h1,t); u1,du1 = eval_pd_f(Dj_h2,t); u2,du2 = eval_pd_f(Dj2_h3,t)
            z01=u0*du1-du0*u1; z12=u1*du2-du1*u2; z02=u0*du2-du0*u2
            if z01<=0 or z12<=0 or z02<=0 or 4*z01*z12 < (1 + 2.0/j)*z02*z02:
                ok = False; break
        if not ok: continue
        valid_g += 1
        # Form U0 = h0 + lam*Dj h1, U1 = h1 + lam*Dj h2, U2 = h2 + lam*Dj h3
        def add_s(pa, pb, sc):
            n = max(len(pa), len(pb))
            out = np.zeros(n)
            out[:len(pa)] += pa; out[:len(pb)] += sc * pb
            return out
        U0 = add_s(h0, Dj_h1, lam); U1 = add_s(h1, Dj_h2, lam); U2 = add_s(h2, apply_Dj_f(h3, j), lam)
        # Apply G_j: X = D_j U0, M = D_{j+1} D_j U1, Y = D_{j+1}^2 D_j U2
        X = apply_Dj_f(U0, j)
        M = apply_Dj_f(apply_Dj_f(U1, j), j + 1)
        Y = apply_Dj_f(apply_Dj_f(apply_Dj_f(U2, j), j + 1), j + 1)
        mJ = 1e9
        for t in np.concatenate([ts_all, [-2.0/j, -1.0/j, -0.5/j, 0.5/j, 2.0/j]]):
            v0,d0 = eval_pd_f(X,t); v1,d1 = eval_pd_f(M,t); v2,d2 = eval_pd_f(Y,t)
            w01=v0*d1-d0*v1; w12=v1*d2-d1*v2; w02=v0*d2-d0*v2
            J = 4*w01*w12/(w02*w02) if (w01>0 and w12>0 and w02>0) else -1.0
            if J < mJ: mJ = J
        if mJ < min_Jout: min_Jout = mJ
        if mJ < 1.0: fails_1 += 1
        if mJ < 1.0 + 2.0/(j+1): fails_target += 1
    print(f"  G_j at j={j:2d}, lam=3/j={lam:.4f} (valid={valid_g:4d}): min J_out={min_Jout:.5f}, fails(<1)={fails_1}, fails(<1+2/(j+1))={fails_target}")

# ============================================================
# 3. FIXED-k LINEAR DIFFERENTIAL RECURRENCE CHECK AT k = 6
# ============================================================
print("\n=== 3. FIXED-k DIFFERENTIAL RECURRENCE CHECK AT k = 6 ===")
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
A = {(1,0): [1, 1], (1,1): [1]}
for k in range(2, 7):
    for i in range(k + 1):
        L = 3*k - 2 - i
        u = [0]
        if i <= k - 1: u = padd_int(u, apply_SL_int(A[(k-1, i)], L - 1))
        if i >= 1: u = padd_int(u, A[(k-1, i-1)])
        A[(k, i)] = apply_SL_int(u, L)
x = sp.Symbol('x')
# Check if A_{6,3}, A'_{6,3}, A_{6,2}, A'_{6,2}, A_{6,4}, A'_{6,4} have ANY deg<=2 polynomial syzygy!
polys = []
for idx in [2, 3, 4]:
    P = sum(c * x**m for m, c in enumerate(A[(6, idx)]))
    dP = sp.diff(P, x)
    for r in [0, 1, 2]:
        polys.append(sp.expand(x**r * P))
        polys.append(sp.expand(x**r * dP))
max_d = max(sp.Poly(P, x).degree() for P in polys)
M = sp.Matrix([[sp.Poly(P, x).nth(m) for P in polys] for m in range(max_d + 1)])
print(f"  k=6, i=3 (18 basis funcs x^r A_{{6,j}}, x^r A'_{{6,j}}, max deg={max_d}): rank = {M.rank()} / {len(polys)} (nullity = {len(polys) - M.rank()})")
