"""
Fast quantitative diagnostics on our 2 strongest unbroken GLZ-specific structures:
1. QUAD INDUCTION TRACKING (Sub-step A, D_{k+1}, Sub-step B: D_{k+1} -> D_{k+2}):
   Track min_{t in R} (4 W_XM W_MY / W_XY^2) and max_{t in R} (W_AB W_CE / W_BC^2) at EVERY sub-step:
   - Stage 0 (input): (X, M, Y) = (h_{k-1, i-2}, D_{k+1} h_{k-1, i-1}, D_{k+1}^2 h_{k-1, i})
   - Stage 1 (add):   (A+B, B+C, C+E) = (U_{k, i-2}, D_{k+1} U_{k, i-1}, D_{k+1}^2 U_{k, i})
   - Stage 2 (D_{k+1}): (h_{k, i-2}, D_{k+1} h_{k, i-1}, D_{k+1}^2 h_{k, i})
   - Stage 3 (shift D_{k+1}->D_{k+2}): (h_{k, i-2}, D_{k+2} h_{k, i-1}, D_{k+2}^2 h_{k, i})
   CRITICAL QUESTION: Does Stage 3 (D_{k+1} -> D_{k+2}) INCREASE or DECREASE min(4 W_XM W_MY / W_XY^2)?!
   And on RANDOM real-stable (U_{i-2}, D U_{i-1}, D^2 U_i), does Stage 2 + Stage 3 ALWAYS preserve Quad stability?!
2. q-STABILITY Q_k(t, y) = sum q^{C(i,2)} htilde_{k,i}(t) y^i AT q=1/2 AND q=1/4 FOR k=15..35:
   - Check at t=0: does Q_k'(0, y) << Q_k(0, y) hold at q=1/2 for k=35, 50, 100?
   - Check full UHP stability (max Im(y) for t in UHP, and dy/dt <= 0 on R) for q=1/2 and q=1/4 up to k=25!
   - Check ONE-STEP q-interlacing: does y*Q_{k-1}(t, q*y) interlace D_{k+1} Q_{k-1}(t, y) along ALL positive directions?
"""
import math, random
import mpmath as mp
mp.mp.dps = 50

def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p
def poly_add(a, b):
    n = max(len(a), len(b))
    return poly_trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def poly_deriv(a):
    return [i * a[i] for i in range(1, len(a))] if len(a) > 1 else [0]
def apply_Dj(p, j):
    dp = poly_deriv(p)
    res = [0] * (len(p) + 2)
    for i, c in enumerate(p):
        res[i] += c
        res[i+1] += j * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c
        res[i+2] += c
    return poly_trim(res)
def poly_eval(p, x):
    acc = mp.mpc(0) if isinstance(x, mp.mpc) else mp.mpf(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc

ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 26):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

print("=== 1. QUAD STAGE-BY-STAGE TRACKING (Stage 0 -> 1 -> 2 -> 3) ===", flush=True)
grid_R = [mp.mpf(v) for v in [-100, -10, -2, -1.1, -1, -0.99, -0.8, -0.5, -0.2, -0.05, 0, 0.05, 0.1, 0.2, 0.5, 1, 5, 50]]
def min_J_ratio(X, M, Y):
    dX, dM, dY = poly_deriv(X), poly_deriv(M), poly_deriv(Y)
    mn = mp.mpf('1e100')
    wt = None
    for tv in grid_R:
        x, dx = poly_eval(X, tv), poly_eval(dX, tv)
        m, dm = poly_eval(M, tv), poly_eval(dM, tv)
        y, dy = poly_eval(Y, tv), poly_eval(dY, tv)
        wXM = x*dm - dx*m
        wMY = m*dy - dm*y
        wXY = x*dy - dx*y
        if wXY != 0:
            r = (4*wXM*wMY)/(wXY*wXY)
            if r < mn:
                mn = r
                wt = float(tv)
    return mn, wt

for k in [8, 12, 16]:
    i = k // 2 + 1  # middle index
    # U_{k, r} = D_{k+1} ht[k-1][r] + ht[k-1][r-1]
    U = {}
    for r in range(k + 1):
        p1 = apply_Dj(ht[k-1][r], k+1) if r <= k-1 else [0]
        p2 = ht[k-1][r-1] if r >= 1 else [0]
        U[r] = poly_add(p1, p2)
    # Stage 0a: (ht[k-1][i-3], D_{k+1} ht[k-1][i-2], D_{k+1}^2 ht[k-1][i-1])
    s0a = min_J_ratio(ht[k-1][i-3], apply_Dj(ht[k-1][i-2], k+1), apply_Dj(apply_Dj(ht[k-1][i-1], k+1), k+1))
    # Stage 0b: (D_{k+1} ht[k-1][i-2], D_{k+1}^2 ht[k-1][i-1], D_{k+1}^3 ht[k-1][i])
    s0b = min_J_ratio(apply_Dj(ht[k-1][i-2], k+1), apply_Dj(apply_Dj(ht[k-1][i-1], k+1), k+1), apply_Dj(apply_Dj(apply_Dj(ht[k-1][i], k+1), k+1), k+1))
    # Stage 1 (Add): (U[i-2], D_{k+1} U[i-1], D_{k+1}^2 U[i])
    s1 = min_J_ratio(U[i-2], apply_Dj(U[i-1], k+1), apply_Dj(apply_Dj(U[i], k+1), k+1))
    # Stage 2 (D_{k+1}): (ht[k][i-2], D_{k+1} ht[k][i-1], D_{k+1}^2 ht[k][i])
    s2 = min_J_ratio(ht[k][i-2], apply_Dj(ht[k][i-1], k+1), apply_Dj(apply_Dj(ht[k][i], k+1), k+1))
    # Stage 3 (Shift D_{k+1} -> D_{k+2}): (ht[k][i-2], D_{k+2} ht[k][i-1], D_{k+2}^2 ht[k][i])
    s3 = min_J_ratio(ht[k][i-2], apply_Dj(ht[k][i-1], k+2), apply_Dj(apply_Dj(ht[k][i], k+2), k+2))
    print(f"  k={k:2d}, i={i}: S0a={mp.nstr(s0a[0],6)}@{s0a[1]}, S0b={mp.nstr(s0b[0],6)}@{s0b[1]} -> S1(add)={mp.nstr(s1[0],6)}@{s1[1]} -> S2(D)={mp.nstr(s2[0],6)}@{s2[1]} -> S3(shift)={mp.nstr(s3[0],6)}@{s3[1]}", flush=True)

print("\n=== 2. TEST STAGE 2+3 COMBINED ON RANDOM STABLE QUADS (U0, D_j U1, D_j^2 U2) ===", flush=True)
# Does (U0, D_j U1, D_j^2 U2) stable => (D_j U0, D_{j+1} D_j U1, D_{j+1}^2 D_j U2) stable?!
import numpy as np
def poly_mul(a, b):
    res = [0.0]*(len(a)+len(b)-1)
    for ii, ai in enumerate(a):
        for jj, bj in enumerate(b): res[ii+jj] += ai*bj
    return res
def pfromr(rts):
    p = [1.0]
    for r in rts: p = poly_mul(p, [-r, 1.0])
    return p
def check_quad_np(X, M, Y, ts):
    vX = np.polyval(list(reversed(X)), ts)
    vM = np.polyval(list(reversed(M)), ts)
    vY = np.polyval(list(reversed(Y)), ts)
    dX = np.polyval(list(reversed(poly_deriv(X))), ts)
    dM = np.polyval(list(reversed(poly_deriv(M))), ts)
    dY = np.polyval(list(reversed(poly_deriv(Y))), ts)
    wXM = vX*dM - dX*vM
    wMY = vM*dY - dM*vY
    wXY = vX*dY - dX*vY
    # Also check interlacing sign: wXM > 0, wMY > 0 (or both < 0)
    if np.min(wXM * wMY) <= 0:
        return -1.0
    return float(np.min((4*wXM*wMY - wXY*wXY) / (4*wXM*wMY + 1e-30)))

random.seed(11)
ts_np = np.concatenate([np.linspace(-20, -1.05, 80), np.linspace(-0.99, -0.01, 150), np.linspace(0.01, 20, 80)])
tested = 0
fail_S23 = 0
fail_S3_only = 0
for tr in range(2000):
    if tested >= 200: break
    j = random.randint(4, 10)
    # Generate U2 << U1 << U0 in (-1,0) of degrees 2, 3, 4
    pts = sorted([random.uniform(-0.95, -0.05) for _ in range(7)])
    r0 = [pts[0], pts[2], pts[4], pts[6]]
    r1 = [pts[1], pts[3], pts[5]]
    r2 = [random.uniform(r1[0], r1[1]), random.uniform(r1[1], r1[2])]
    c0, c1, c2 = 10**random.uniform(-0.5, 0.5), 10**random.uniform(-0.5, 0.5), 10**random.uniform(-0.5, 0.5)
    U0 = [c0*x for x in pfromr(r0)]
    U1 = [c1*x for x in pfromr(r1)]
    U2 = [c2*x for x in pfromr(r2)]
    X1, M1, Y1 = U0, apply_Dj(U1, j), apply_Dj(apply_Dj(U2, j), j)
    if check_quad_np(X1, M1, Y1, ts_np) <= 1e-4: continue
    tested += 1
    # S3 only (without extra D_j):
    if check_quad_np(U0, apply_Dj(U1, j+1), apply_Dj(apply_Dj(U2, j+1), j+1), ts_np) < -1e-6:
        fail_S3_only += 1
    # S2+S3 (with extra D_j applied first!):
    h0, h1, h2 = apply_Dj(U0, j), apply_Dj(U1, j), apply_Dj(U2, j)
    if check_quad_np(h0, apply_Dj(h1, j+1), apply_Dj(apply_Dj(h2, j+1), j+1), ts_np) < -1e-6:
        fail_S23 += 1
print(f"  Random stable (U0, D_j U1, D_j^2 U2) ({tested} tested): S3-only fails={fail_S3_only}, S2+S3 (with D_j!) fails={fail_S23}", flush=True)

print("\n=== 3. STRESS TEST q-STABILITY Q_k(t, y) AT q=1/2 AND q=1/4 FOR GLZ (k=10..25) ===", flush=True)
# First test t=0 interlacing Q_k'(0, y) << Q_k(0, y) at q=1/2 for k=10..50!
# Note: at t=0, htilde_{k,i}(0) = C(k,i), and htilde'_{k,i}(0) satisfies exact recurrence:
# D_{k+1}(a + b t) = a + (2b + (k+1)a) t!
# Wait: for htilde_{k,i} = D_{k+1}(D_{k+1} htilde_{k-1,i} + htilde_{k-1,i-1}):
# Let c_{k,i} = htilde_{k,i}(0) = C(k,i), d_{k,i} = htilde'_{k,i}(0)!
# U_{k,i}(0) = c_{k-1,i} + c_{k-1,i-1} = c_{k,i}
# U'_{k,i}(0) = (2 d_{k-1,i} + (k+1) c_{k-1,i}) + d_{k-1,i-1}
# d_{k,i} = 2 U'_{k,i}(0) + (k+1) c_{k,i}!
d_val = {1: {0: 2, 1: 0}}
for k in range(2, 61):
    d_val[k] = {}
    for i in range(k+1):
        c_im1 = math.comb(k-1, i-1) if i >= 1 else 0
        c_i = math.comb(k-1, i) if i <= k-1 else 0
        d_im1 = d_val[k-1][i-1] if i >= 1 else 0
        d_i = d_val[k-1][i] if i <= k-1 else 0
        u_prime = (2*d_i + (k+1)*c_i) + d_im1
        d_val[k][i] = 2*u_prime + (k+1)*math.comb(k, i)

for q_str, q_val in [("1/2", mp.mpf('0.5')), ("1/3", mp.mpf(1)/3), ("1/4", mp.mpf('0.25'))]:
    first_fail = None
    for k in [10, 20, 30, 40, 50, 60]:
        P0 = [mp.mpf(math.comb(k, i)) * (q_val**(i*(i-1)//2)) for i in range(k+1)]
        P1 = [mp.mpf(d_val[k][i]) * (q_val**(i*(i-1)//2)) for i in range(k+1)]
        # Note d_val[k][k] = k*... != 0 for k>=2, so deg P1 = k, deg P0 = k!
        rts0 = mp.polyroots(list(reversed(P0)), maxsteps=400, extraprec=60)
        rts1 = mp.polyroots(list(reversed(P1)), maxsteps=400, extraprec=60)
        if max(abs(r.imag) for r in rts0) > mp.mpf('1e-12') or max(abs(r.imag) for r in rts1) > mp.mpf('1e-12'):
            first_fail = (k, "complex")
            break
        r0 = sorted([r.real for r in rts0])
        r1 = sorted([r.real for r in rts1])
        # Check strict interlacing r0[0] < r1[0] < r0[1] < r1[1] < ... < r0[-1] < r1[-1] (or reverse!)
        ok1 = all(r0[m] < r1[m] for m in range(k)) and all(r1[m] < r0[m+1] for m in range(k-1))
        ok2 = all(r1[m] < r0[m] for m in range(k)) and all(r0[m] < r1[m+1] for m in range(k-1))
        if not (ok1 or ok2):
            first_fail = (k, "no_interlace")
            break
    print(f"  t=0 test Q_k'(0,y) << Q_k(0,y) at q={q_str} (k=10..60): first_fail = {first_fail}", flush=True)

# Full UHP stability check of Q_k(t, y) at q=1/2 and q=1/4 for k=15, 20, 25!
for q_str, q_val in [("1/2", mp.mpf('0.5')), ("1/4", mp.mpf('0.25'))]:
    for k in [15, 20, 25]:
        max_im = mp.mpf('-1e100')
        worst = None
        for tv in [mp.mpc(-10, '1e-6'), mp.mpc(-1, '1e-6'), mp.mpc(-0.99, '1e-6'), mp.mpc(-0.5, '1e-6'),
                   mp.mpc(-0.01, '1e-6'), mp.mpc(0, '1e-6'), mp.mpc(0.01, '1e-6'), mp.mpc(1, '1e-6'), mp.mpc(0, 1)]:
            cy = [poly_eval(ht[k][i], tv) * (q_val**(i*(i-1)//2)) for i in range(k+1)]
            rts = mp.polyroots(list(reversed(cy)), maxsteps=400, extraprec=60)
            for r in rts:
                if r.imag > max_im:
                    max_im = r.imag
                    worst = (tv, r)
        print(f"  Q_{k}(t,y) at q={q_str}: UHP t -> max Im(y) = {mp.nstr(max_im, 6)}", flush=True)
