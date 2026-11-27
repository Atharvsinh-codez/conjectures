"""
Two decisive scripts in one (< 15s total):
PART 1: Settle fixed w_j = 1, b_j >= 2 at k = 120, 160 for (D_2 + y)^k + (D_{b*} + y)!
   - Does U_Turan(w=1) ACTUALLY go negative at k=120, 160?!
   - And if so, does ONE MORE step (b_{k+2} = 1e40, w = 1) break INTERLACING with FIXED w = 1?!
PART 2: Route I Stage 3 closure via the smooth j-pole function r_i(t) = -W(h_{i-1}, D_0 h_i) / W(h_{i-1}, t h_i)!
   - Notice: because h_i << h_{i-1} in (-1,0) and deg h_{i-1} = deg h_i + 1,
     t h_i has SAME degree as h_{i-1} and root at 0 > all roots of h_{i-1}, so h_{i-1} << t h_i
     => W(h_{i-1}, t h_i) > 0 STRICTLY ON ALL OF R!!
   - Check on GLZ (k = 4..16):
     (a) max_{i, t in R} r_i(t) / (k + 1)
     (b) max_{i, t in R} |r_{i+1}(t) - r_i(t)|   (and min/max of r_{i+1}(t) - r_i(t)!)
     (c) Check whether (r_{i+1} - r_i) is BOUNDED on ALL of R AND check how Stage 1 & Stage 2 transform r_i(t)!
   - Also check the EXACT universal Stage 3 bound on d=2,3 triples when BOTH
     r_0(t), r_1(t) <= 0.75 j  AND  |r_1(t) - r_0(t)| <= 1 (or <= 0.5)!
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 120

# ============================================================
# PART 1: TEST (D_2 + y)^k + (D_{b*} + y) WITH FIXED w = 1 AT k = 120, 160!
# ============================================================
print("=== PART 1: FIXED w = 1 ATTACK ON (D_2 + y)^k AT k = 120, 160 ===")
# Build p_m(t) = D_2^m(1) up to m = 85 (takes 0.3s in mpmath!)
p = [[mp.mpf(1)]]
for m in range(1, 86):
    prev = p[-1]
    nxt = [mp.mpf(0)] * (len(prev) + 1)
    for idx, c in enumerate(prev):
        nxt[idx] += (1 + idx) * c
        nxt[idx + 1] += (2 + idx) * c
    p.append(nxt)

def eval_pd(poly, tv):
    v = mp.mpf(0); dv = mp.mpf(0)
    for c in reversed(poly):
        dv = dv * tv + v
        v = v * tv + c
    return v, dv

for k in [80, 120, 160]:
    i = k // 2
    c_m1, c_0, c_1, c_2 = mp.binomial(k, i-2), mp.binomial(k, i-1), mp.binomial(k, i), mp.binomial(k, i+1)
    min_U = mp.mpf('1e100')
    best = None
    for alpha in [0.01, 0.03, 0.05, 0.08, 0.1, 0.12, 0.15, 0.2, 0.25, 0.3]:
        t = -mp.mpf(str(alpha)) * (mp.mpf(2) ** (-(k - i)))
        vm1, dm1 = eval_pd(p[k - (i - 2)], t); vm1 *= c_m1; dm1 *= c_m1
        v0, d0 = eval_pd(p[k - (i - 1)], t); v0 *= c_0; d0 *= c_0
        v1, d1 = eval_pd(p[k - i], t); v1 *= c_1; d1 *= c_1
        v2, d2 = eval_pd(p[k - (i + 1)], t); v2 *= c_2; d2 *= c_2
        u0, u1, u2 = t*(1+t)*d0, t*(1+t)*d1, t*(1+t)*d2
        q0, q1, q2 = u0 + vm1, u1 + v0, u2 + v1
        A_U = q1**2 - q0*q2
        B_U = 2*q1*v1 - (q0*v2 + q2*v0)
        C_U = v1**2 - v0*v2
        if C_U > 0:
            PhiU = min(1 + 2*t, -B_U / (2*C_U))
            valU = (A_U + PhiU*B_U + (PhiU**2)*C_U) / (v1**2)
            if valU < min_U:
                min_U = valU
                best = (alpha, float(PhiU), float((1 - PhiU)/(-t)), t, PhiU)
    print(f"  k={k:3d} (i={i}): min U_Turan(w=1)/v1^2 = {float(min_U):+.5e} at alpha={best[0]}, Phi*={best[1]:.4f}, b*={best[2]:.3e}")
    if min_U < 0:
        # Check if ONE MORE step with FIXED w=1 breaks INTERLACING W(U2_next, U1_next) < 0!
        # Wait: at step k+1 with b* = (1 - Phi*)/(-t), we have U_{k+1, i-1}, U_{k+1, i}, U_{k+1, i+1}
        # with U_i(t)^2 < U_{i-1}(t) U_{i+1}(t)!
        # Can we choose b_{k+2} >= 2 (with w = 1!) so that U_{k+2, i+1}(t) = D_{b_{k+2}} U_{k+1, i+1}(t) + U_{k+1, i}(t) = 0?!
        # That requires D_{b_{k+2}} U_{k+1, i+1}(t) = -U_{k+1, i}(t), i.e.,
        #   t(1+t) U'_{i+1}(t) + (1 + b_{k+2} t) U_{i+1}(t) = -U_i(t)!
        # Let's compute the EXACT b_{k+2} needed and W(U2_next, U1_next)(t)!
        t_s = best[3]
        b_s = (1 - best[4]) / (-t_s)
        def apply_Dj(poly, jv):
            out = [mp.mpf(0)] * (len(poly) + 1)
            for idx, c in enumerate(poly):
                out[idx] += (1 + idx) * c
                out[idx + 1] += (jv + idx) * c
            return out
        def padd(a, b):
            n = max(len(a), len(b))
            out = [mp.mpf(0)] * n
            for idx in range(len(a)): out[idx] += a[idx]
            for idx in range(len(b)): out[idx] += b[idx]
            return out
        def pscale(a, s): return [s*c for c in a]
        h_m1 = pscale(p[k - (i - 2)], c_m1)
        h_0  = pscale(p[k - (i - 1)], c_0)
        h_1  = pscale(p[k - i], c_1)
        h_2  = pscale(p[k - (i + 1)], c_2)
        U0 = padd(apply_Dj(h_0, b_s), h_m1)
        U1 = padd(apply_Dj(h_1, b_s), h_0)
        U2 = padd(apply_Dj(h_2, b_s), h_1)
        u0v, u0d = eval_pd(U0, t_s)
        u1v, u1d = eval_pd(U1, t_s)
        u2v, u2d = eval_pd(U2, t_s)
        # Theorem A at r = t_s (root of D_{b2} U2 + U1 = 0):
        # Need b2 such that t_s*(1+t_s)*u2d + (1 + b2*t_s)*u2v + u1v = 0
        b2 = (-(u1v + t_s*(1+t_s)*u2d)/u2v - 1) / t_s
        W_U2_U1 = u2v*u1d - u2d*u1v
        turan_U = u1v**2 - u0v*u2v
        # Theorem A identity: (D_{b2} U1 + U0)(t_s) * U2(t_s) = t_s*(1+t_s)*W(U2, U1) - 1*(U1^2 - U0 U2)
        # For interlacing (U2_next << U1_next), at a root of U2_next we need U1_next * U2 < 0 (Wait: (-t)(1+t)W > 0 => t(1+t)W < 0!)
        lhs_thmA = t_s*(1+t_s)*W_U2_U1 - turan_U
        print(f"    -> Check step k+2 (FIXED w=1!): b2={float(b2):+.3e}, t(1+t)W(U2,U1)={float(t_s*(1+t_s)*W_U2_U1):+.3e}, -turan_U={float(-turan_U):+.3e}, sum={float(lhs_thmA):+.3e}")

# ============================================================
# PART 2: ROUTE I j-POLE FUNCTION r_i(t) = -W(h_{i-1}, D_0 h_i) / W(h_{i-1}, t h_i) ON GLZ!
# ============================================================
print("\n=== PART 2: SMOOTH j-POLE FUNCTION r_i(t) AND ITS DIFFERENCE r_{i+1}(t) - r_i(t) ON GLZ ===")
mp.mp.dps = 60

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

def mul_t_mp(poly):
    return [mp.mpf(0)] + list(poly)

h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
h_store = {1: h}
for k in range(2, 17):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = padd_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = padd_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
    h_store[k] = h

ts_r = [mp.mpf(x) for x in [-100.0, -5.0, -1.2, -1.01, -0.99, -0.8, -0.5, -0.2, -0.05, -1e-3, -1e-6, 0.0, 1e-4, 0.05, 0.2, 1.0, 5.0, 100.0]]
for k in [4, 6, 8, 10, 12, 14, 16]:
    hk = h_store[k]
    # Compute r_i(t) for i = 1..k (where r_i(t) is the unique j-root of W(h_{k,i-1}, D_j h_{k,i}) = 0)
    r_vals = {}
    max_r = -1e9
    max_r_where = None
    for i in range(1, k + 1):
        f0, f1 = hk[i - 1], hk[i]
        D0_f1 = apply_Dj_mp(f1, 0)
        t_f1 = mul_t_mp(f1)
        r_vals[i] = []
        for t in ts_r:
            v0 = eval_pd(f0, t)
            vD0 = eval_pd(D0_f1, t)
            vt = eval_pd(t_f1, t)
            w_num = v0[0]*vD0[1] - v0[1]*vD0[0]
            w_den = v0[0]*vt[1] - v0[1]*vt[0]
            rv = float(-w_num / w_den)
            r_vals[i].append(rv)
            if rv > max_r:
                max_r = rv
                max_r_where = (i, float(t))
    max_diff = -1e9
    min_diff = 1e9
    max_abs_diff = 0.0
    diff_where = None
    for i in range(1, k):
        for ti, t in enumerate(ts_r):
            df = r_vals[i + 1][ti] - r_vals[i][ti]
            if df > max_diff: max_diff = df
            if df < min_diff: min_diff = df
            if abs(df) > max_abs_diff:
                max_abs_diff = abs(df)
                diff_where = (i, float(t), df)
    print(f"  k={k:2d}: max r_i/(k+1) = {max_r/(k+1):.4f} (max_r={max_r:.2f} at {max_r_where}) | r_{{i+1}}-r_i in [{min_diff:+.4f}, {max_diff:+.4f}], max|diff|={max_abs_diff:.4f} at {diff_where}")
