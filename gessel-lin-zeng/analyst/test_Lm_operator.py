"""Test the shifted-index operator L_m = D_m + t*y*d/dy:
  L_m(X + 2 M y + Y y^2) = (D_m X) + 2 (D_{m+1} M) y + (D_{m+2} Y) y^2 !
1. Does L_m preserve real stability (J >= 1 on R) on random real-stable quadratics X + 2My + Yy^2
   with ALL roots in (-1, 0) (and also check BOTH orientations: X < M < Y and Y < M < X!!)?!
2. WAIT A SECOND: Check BOTH signs of t*y*d/dy:
   (a) (D_m X, D_{m+1} M, D_{m+2} Y)  (m INCREASING with y-power, matching GLZ!)
   (b) (D_m X, D_{m-1} M, D_{m-2} Y)  (m DECREASING with y-power!)
3. Compare GLZ's ACTUAL J_ratio(D_{k+1} U_0, D_{k+2} D_{k+1} U_1, D_{k+2}^2 D_{k+1} U_2)
   with L_{k+1}(U_0, D_{k+1} U_1, D_{k+1}^2 U_2) = (D_{k+1} U_0, D_{k+2} D_{k+1} U_1, D_{k+3} D_{k+1}^2 U_2)!
"""
import numpy as np
import mpmath as mp

def apply_D_flt(m, p):
    res = np.zeros(len(p) + 1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def eval_W(f, g, ts):
    fp = np.arange(1, len(f))*f[1:]; gp = np.arange(1, len(g))*g[1:]
    return np.polyval(fp[::-1], ts)*np.polyval(g[::-1], ts) - np.polyval(f[::-1], ts)*np.polyval(gp[::-1], ts)

print("=== 1. TEST L_m: (X, M, Y) -> (D_m X, D_{m+1} M, D_{m+2} Y) ON RANDOM REAL-STABLE QUADRATICS IN (-1, 0) ===")
rng = np.random.default_rng(2026)
ts = np.concatenate([-10**np.linspace(2, -4, 40), np.linspace(-0.99, -0.01, 100), [0.0], 10**np.linspace(-4, 2, 40)])

for direction in ["X < M < Y (GLZ orientation)", "Y < M < X (reverse)"]:
    inc_fails = 0; dec_fails = 0; valid = 0
    worst_J_inc = 1e9; worst_J_dec = 1e9
    max_J_drop_inc = -1e9
    for _ in range(500):
        # Generate random real-stable quadratic P(t, y) = det(I + t A + y B) with roots in (-1, 0) at y=0!
        # Even simpler: P(t, y) = (1 + u1 t + v1 y)(1 + u2 t + v2 y) * prod_{r=1}^2 (1 + w_r t)
        # Wait: for general non-factorable quadratic in y of degree d in t:
        # Take random 2x2 PSD pencil in y? No, take P(t, y) = f(t) + 2 g(t) y + h(t) y^2 from a 2-variable product:
        # P(t, y) = prod_{j=1}^2 (p_j(t) + q_j(t) y) where p_j < q_j on (-1, 0)!!
        # Let's build p_1, q_1 (deg 2, p_1 < q_1 in (-1, 0)) and p_2, q_2 (deg 2, p_2 < q_2 in (-1, 0))!
        pts1 = -np.sort(rng.uniform(0.02, 0.98, 4))[::-1]  # ascending in (-1, 0)
        pts2 = -np.sort(rng.uniform(0.02, 0.98, 4))[::-1]
        def from_rts(rs):
            p = np.array([1.0])
            for r in rs: p = np.convolve(p, [1.0, -1.0/r])  # (1 - t/r), positive at t=0, roots at r in (-1, 0)
            return p
        p1, q1 = from_rts(pts1[0::2]), from_rts(pts1[1::2])  # p1 < q1 (roots of p1 strictly left of q1!)
        p2, q2 = from_rts(pts2[0::2]), from_rts(pts2[1::2])  # p2 < q2
        # Product (p1 + q1 y)(p2 + q2 y) = (p1 p2) + (p1 q2 + p2 q1) y + (q1 q2) y^2 -> X=p1 p2, 2M=p1 q2+p2 q1, Y=q1 q2!
        X = np.convolve(p1, p2)
        M = 0.5 * (np.convolve(p1, q2) + np.convolve(p2, q1))
        Y = np.convolve(q1, q2)
        if direction.startswith("Y <"):
            X, Y = Y, X
        wxm = eval_W(X, M, ts); wmy = eval_W(M, Y, ts); wxy = eval_W(X, Y, ts)
        j_in = float(np.min(4*wxm*wmy / (wxy*wxy + 1e-300)))
        if j_in < 1.0 - 1e-9: continue
        valid += 1
        m = float(rng.uniform(4.0, 20.0))
        # (a) m, m+1, m+2
        X1, M1, Y1 = apply_D_flt(m, X), apply_D_flt(m+1, M), apply_D_flt(m+2, Y)
        w1 = eval_W(X1, M1, ts); w2 = eval_W(M1, Y1, ts); w12 = eval_W(X1, Y1, ts)
        j_inc = float(np.min(4*w1*w2 / (w12*w12 + 1e-300)))
        if j_inc < 1.0 - 1e-7 or np.min(w1*w2) < -1e-9: inc_fails += 1
        worst_J_inc = min(worst_J_inc, j_inc)
        max_J_drop_inc = max(max_J_drop_inc, j_in - j_inc)
        # (b) m, m-1, m-2
        X2, M2, Y2 = apply_D_flt(m, X), apply_D_flt(m-1, M), apply_D_flt(m-2, Y)
        w1 = eval_W(X2, M2, ts); w2 = eval_W(M2, Y2, ts); w12 = eval_W(X2, Y2, ts)
        j_dec = float(np.min(4*w1*w2 / (w12*w12 + 1e-300)))
        if j_dec < 1.0 - 1e-7 or np.min(w1*w2) < -1e-9: dec_fails += 1
        worst_J_dec = min(worst_J_dec, j_dec)
    print(f"  {direction:28s} (valid={valid}):")
    print(f"    (m, m+1, m+2): fails={inc_fails}/{valid}, worst_J={worst_J_inc:.5f}, max(J_in-J_out)={max_J_drop_inc:+.5f}")
    print(f"    (m, m-1, m-2): fails={dec_fails}/{valid}, worst_J={worst_J_dec:.5f}")

print("\n=== 2. TEST EVEN MORE GENERAL (D_m, D_{m+s}, D_{m+2s}) ON X < M < Y FOR s = 0.5, 1, 5, 50 ===")
# What if X, M, Y are close (like GLZ: M = X + eps*Z, Y = X + 2*eps*Z)?!
close_fails = 0
for s_step in [0.1, 1.0, 5.0, 50.0]:
    fails = 0; min_J = 1e9; valid = 0
    for _ in range(200):
        pts1 = -np.sort(rng.uniform(0.02, 0.98, 4))[::-1]
        pts2 = -np.sort(rng.uniform(0.02, 0.98, 4))[::-1]
        def from_rts(rs):
            p = np.array([1.0])
            for r in rs: p = np.convolve(p, [1.0, -1.0/r])
            return p
        p1, q1 = from_rts(pts1[0::2]), from_rts(pts1[1::2])
        p2, q2 = from_rts(pts2[0::2]), from_rts(pts2[1::2])
        # Make q1 close to p1, q2 close to p2!
        q1 = p1 + 0.05 * q1; q2 = p2 + 0.05 * q2
        X = np.convolve(p1, p2)
        M = 0.5 * (np.convolve(p1, q2) + np.convolve(p2, q1))
        Y = np.convolve(q1, q2)
        valid += 1
        m = 5.0
        X1, M1, Y1 = apply_D_flt(m, X), apply_D_flt(m+s_step, M), apply_D_flt(m+2*s_step, Y)
        w1 = eval_W(X1, M1, ts); w2 = eval_W(M1, Y1, ts); w12 = eval_W(X1, Y1, ts)
        j_out = float(np.min(4*w1*w2 / (w12*w12 + 1e-300)))
        if j_out < 1.0 - 1e-7 or np.min(w1*w2) < -1e-9: fails += 1
        min_J = min(min_J, j_out)
    print(f"  Close X<M<Y, s={s_step:4.1f}: fails={fails}/{valid}, min_J={min_J:.6f}")
