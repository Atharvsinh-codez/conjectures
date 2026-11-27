"""CRITICAL TEST OF L_m = (D_m, D_{m+1}, D_{m+2}) ON GENERAL (NON-FACTORED!) REAL-STABLE QUADRATICS IN (-1, 0):
1. Generate NON-FACTORED real-stable quadratics P(t, y) = X(t) + 2 M(t) y + Y(t) y^2 with ALL t-roots in (-1, 0)
   using P(t, y) = det(I + t*A + y*B) where A > I (eigenvalues > 1 so t = -1/lam in (-1, 0)!) and B >= 0,
   AND by taking convex combos / derivatives of products, AND pushing J_in -> 1.0001!
   Test:
   (a) Does P(t, (1+t)y) stay REAL-STABLE?!
   (b) Does L_m(X, M, Y) = (D_m X, D_{m+1} M, D_{m+2} Y) preserve J >= 1 + eps with ZERO LOSS (J_out >= J_in)?!
2. On GLZ (k = 3..14), compare the 3 triples:
   - T_old = (U_0, D_{k+1} U_1, D_{k+1}^2 U_2)
   - T_Lm  = (D_{k+1} U_0, D_{k+2} D_{k+1} U_1, D_{k+3} D_{k+1}^2 U_2)  = L_{k+1}(T_old)
   - T_GLZ = (D_{k+1} U_0, D_{k+2} D_{k+1} U_1, D_{k+2}^2 D_{k+1} U_2)  (differs ONLY by +t(1+2t) D_{k+1} U_2 in slot 3!)
   - T_alt = (D_{k+1} U_0, D_{k+2} D_{k+1} U_1, D_{k+1} D_{k+3} D_{k+1} U_2)  (differs from T_GLZ ONLY by +t D_{k+1} U_2!)
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

ts = np.concatenate([-10**np.linspace(2, -4, 40), np.linspace(-0.99, -0.01, 100), [0.0], 10**np.linspace(-4, 2, 40)])

print("=== 1. NON-FACTORED REAL-STABLE QUADRATICS IN (-1, 0): TEST P(t, (1+t)y) AND L_m ===")
# Construct general non-factored real-stable quadratic X(t) + 2 M(t) y + Y(t) y^2 of deg d=4 in t:
# Take det(I_4 + t * A + y * B) where A is 4x4 symmetric with eigenvalues > 1 (so roots in (-1, 0)),
# and B is 4x4 rank-2 PSD matrix (so degree in y is EXACTLY 2!!)!
# Wait: check orientation! For A > I, B >= 0, det(I + t A + y B) = 0 with t in H+, y in H+ -> impossible (Im > 0)!
# Wait: check sign of W(X, M): at t=0, det(I + y B) = 1 + tr(B) y + e_2(B) y^2 (positive roots? No, negative roots!).
# If t in H+, y in H+ has no zeros, then as y increases, t-roots DECREASE (move left: Y < M < X)!
# To get GLZ orientation X < M < Y (t in H+, y in H- has no zeros), substitute y -> -y?
# Wait: if y -> -y, roots move RIGHT as y increases, which is X < M < Y — EXCEPT at y>0 small, roots could cross 0!
# Better: take det(I + t A + (1 - y) B) = det((I + B) + t A - y B), NO: take det((I + B) + t(A + B) - y B)?
# Simplest exact parametrization of ALL real-stable quadratics with X < M < Y in (-1, 0):
# Take A > I, B >= 0 (rank 2) with A - B > I (so even at y = -1, roots of I + t A - B are in (-1, 0)? Wait:
# in det(I + t A - y B), at y=0 roots are -1/eig(A) in (-1, 0); at y=inf, coeff of y^2 is e_2(B) * det(I + t A_perp) (deg 2, not deg 4)!
# To have EQUAL DEGREE d=4 for X, M, Y AND X < M < Y in (-1, 0):
# Take P(t, y) = det(I_4 + t * (A_0 + y * A_1)) where A_0 > I_4, A_1 >= 0 (rank 2)!
# CHECK THAT!!
# - At y = 0: det(I + t A_0) has deg 4, roots -1/eig(A_0) in (-1, 0).
# - Wait: as y increases (y > 0), eigenvalues of A_0 + y A_1 INCREASE, so roots -1/eig(A_0 + y A_1) MOVE RIGHT (toward 0^-!!),
#   staying INSIDE (-1, 0) FOR ALL y >= 0!!
# - Coeff of y^2 has deg 4? Wait: rank(A_1) = 2 => t*y*A_1 contributes t^2 y^2 * (deg 2 in t) = deg 4 in t!!
#   Actually Y(t) = [y^2] det(I + t A_0 + t y A_1) has factor t^2 (roots at 0)!
# To avoid factor t^2, take det(I_4 + t A_0 + (1 + t) y A_1)!! No, that ALREADY has (1+t)!
# Just take det( (I + y C_0) + t (A_0 + y A_1) ) where (C_0, A_1) is rank-2 PSD pencil with A_1 > C_0 > 0!!
# CHECK: at y=0, det(I + t A_0) has roots in (-1, 0). For any y >= 0, A_0 + y A_1 > I + y C_0, so ALL roots of
# det(I + y C_0 + t(A_0 + y A_1)) lie in (-1, 0)!!
# And is det( (I + t A_0) + y (C_0 + t A_1) ) real-stable in (t, y)?
# Let u = t in H+, v = y in H+! No, t*y has product of H+!
# Wait: for X < M < Y, P(t, -y) is real-stable in (t, y) in H+ x H+!
# That means P(t, y) = det( (I + t A_0) - y (C_0 + t A_1) )? No, if t, y in H+, Im(I + t A_0 + y B_0) > 0 ALWAYS,
# so det(I + t A_0 + y B_0) has Im(t) and Im(y) OPPOSITE => Y < M < X (roots move LEFT as y increases)!
# Therefore, replacing y -> 1/y (keeping deg 2: y^2 P(t, 1/y) = Y + 2 M y + X y^2) reverses to X < M < Y!!
# Check: if P(t, y) = det(I + t A_0 + y B_0) with rank(B_0) = 2, wait, Y has deg 2 instead of 4!
# To have equal deg 4 in t and deg 2 in y with X < M < Y: just take det(M_0 + t M_1 + y M_2) with M_1 > 0, M_2 >= 0 (rank 2)
# and evaluate at y -> y_0 - y? No, M_0 + t M_1 + y M_2 has deg 2 in y and deg 4 in t ONLY if M_2 is rank 2 (then [y^2] has deg 2 in t).
# Wait: ANY real-stable polynomial of deg 2 in y and equal deg d in t with X < M < Y is:
#   take ANY real-stable Q(t, w) = X + 2 M w + Y w^2 (from a sum of rank-2 or products of (p_j + q_j w) + convex combos in a pencil)!
# Even better: take P(t, y) = det( (1+t) I_2 + t A_2 + y B_2 )? No, just use (I_2 + t A_1 + y B_1) * ... NO, use
# Helton-Vinnikov for deg-2 in y, deg-d in t: take ANY d x d real symmetric H_0, H_1 with H_0 > 0, H_1 >= 0 (rank 2)?
# WAIT: for X, M, Y ALL of degree d=4 and X < M < Y, take:
#   P(t, y) = det( t * I_4 + S_0 - y * S_1 ) where 0 < S_0 < I_4 (roots of X = det(t I + S_0) are -eig(S_0) in (-1, 0)!)
#   and S_1 >= 0 (rank 2!) with S_0 - S_1 > 0!!
# CHECK THAT!!
#   - [y^0]: det(t I + S_0), deg 4, roots in (-1, 0)!
#   - [y^1]: deg 3 in t... wait, to get deg 4 in t for [y^2], S_1 MUST be multiplied by t? No, wait:
#     if P(t, y) = X(t) + 2 M(t) y + Y(t) y^2 has deg(X)=deg(M)=deg(Y)=d, then at t -> inf, P(t, y) ~ t^d (c0 + 2 c1 y + c2 y^2),
#     so the coefficient of t^d MUST have degree 2 in y!!
#     Therefore in det(A_0 + t A_1 + y A_2), to have [t^d y^2] != 0, we need deg_t + deg_y = d + 2, so the matrix size MUST be d + 2 = 6,
#     with rank(A_1) = 4, rank(A_2) = 2!!
# YES!! Take 6x6 PSD matrices A_1 (rank 4) and A_2 (rank 2), and A_0 > 0 (or with roots in (-1,0)):
#   P(t, z) = det(A_0 + t A_1 + z A_2) has deg 4 in t, deg 2 in z, AND [t^4 z^2] != 0 (ALL THREE X, M, Y HAVE DEGREE 4!!)!
#   And with z = 1/y (i.e. reversing y^0, y^1, y^2!), we get X < M < Y AND ALL t-roots of X, M, Y in (-1, 0) whenever A_1 > A_0 > 0!

rng = np.random.default_rng(777)
valid = 0; p1t_fails = 0; Lm_fails = 0; Lm_J_drop_cnt = 0; max_J_drop = -1e9
worst_J_out = 1e9; worst_p1t_im = 0.0

for _ in range(500):
    # Build 6x6 matrix A_0 > 0, A_1 >= 0 (rank 4), A_2 >= 0 (rank 2), with A_1 > A_0 on range(A_1)!
    # Easiest: choose random orthonormal 6x6 U, let A_1 = diag(u1, u2, u3, u4, 0, 0) with u_i in (1.5, 10.0),
    # A_2 = U @ diag(v1, v2, 0, 0, 0, 0) @ U.T with v_i > 0, and A_0 = I_6 (so all t-roots of det(I + t A_1) are -1/u_i in (-1, 0)!!).
    # Check roots of Y(t) = [z^2] det(I + t A_1 + z A_2): they are roots of det(I_4 + t * (compressed A_1)) -> by Cauchy interlacing,
    # eigenvalues of compressed A_1 lie in [min u_i, max u_i] subset (1.5, 10), wait: A_1 has two 0 eigenvalues!
    # Wait: compressing A_1 to kernel(A_2) (dim 4) CAN have eigenvalues < 1 unless A_1|_ker(A_2) > I_4!
    # So just pick u_i in (5.0, 20.0) and verify all roots of X, M, Y lie in (-1, 0)!
    u = rng.uniform(5.0, 20.0, 4)
    A1 = np.diag([u[0], u[1], u[2], u[3], 0.0, 0.0])
    Q, _ = np.linalg.qr(rng.normal(0, 1, (6, 6)))
    v = rng.uniform(0.5, 5.0, 2)
    A2 = Q @ np.diag([v[0], v[1], 0.0, 0.0, 0.0, 0.0]) @ Q.T
    A0 = np.eye(6)
    # Extract coeffs of t^r z^j (r = 0..4, j = 0..2) by interpolation!
    t_pts = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    z_pts = np.array([0.0, 1.0, -1.0])
    grid = np.array([[np.linalg.det(A0 + tv*A1 + zv*A2) for zv in z_pts] for tv in t_pts])
    C = np.linalg.solve(np.vander(t_pts, 5, increasing=True), np.linalg.solve(np.vander(z_pts, 3, increasing=True), grid.T).T)
    # Reverse z -> 1/y so X = [z^2], 2M = [z^1], Y = [z^0] has X < M < Y!
    # Check roots of X, M, Y are ALL in (-1, 0)!
    X = C[:, 2]; M = 0.5 * C[:, 1]; Y = C[:, 0]
    if abs(X[-1]) < 1e-8 or abs(Y[-1]) < 1e-8: continue
    rx = np.roots(X[::-1]); ry = np.roots(Y[::-1])
    if np.max(np.abs(rx.imag)) > 1e-7 or np.max(np.abs(ry.imag)) > 1e-7: continue
    if np.min(rx.real) <= -0.99 or np.max(rx.real) >= -0.01: continue
    if np.min(ry.real) <= -0.99 or np.max(ry.real) >= -0.01: continue
    # Ensure X < M < Y orientation (W(X, M) > 0)!
    wxm = eval_W(X, M, ts); wmy = eval_W(M, Y, ts); wxy = eval_W(X, Y, ts)
    if np.min(wxm) <= 0 or np.min(wmy) <= 0:
        X, Y = Y, X
        wxm = eval_W(X, M, ts); wmy = eval_W(M, Y, ts); wxy = eval_W(X, Y, ts)
    if np.min(wxm) <= 0 or np.min(wmy) <= 0: continue
    j_in = float(np.min(4*wxm*wmy / (wxy*wxy + 1e-300)))
    if j_in < 1.0 - 1e-7: continue
    valid += 1

    # (a) Test if P(t, (1+t)y) = X + 2(1+t)M y + (1+t)^2 Y y^2 is real-stable!
    M1t = np.convolve(M, [1.0, 1.0]); Y1t = np.convolve(Y, [1.0, 2.0, 1.0])
    # Pad X to deg 6
    Xp = np.pad(X, (0, 2)); M1p = np.pad(M1t, (0, 1))
    w1a = eval_W(Xp, M1p, ts); w2a = eval_W(M1p, Y1t, ts); w12a = eval_W(Xp, Y1t, ts)
    if np.min(4*w1a*w2a - w12a*w12a) < -1e-7 * (np.max(w12a*w12a) + 1):
        p1t_fails += 1

    # (b) Test L_m(X, M, Y) = (D_m X, D_{m+1} M, D_{m+2} Y)!
    m = float(rng.uniform(3.0, 20.0))
    Xo, Mo, Yo = apply_D_flt(m, X), apply_D_flt(m+1, M), apply_D_flt(m+2, Y)
    w1 = eval_W(Xo, Mo, ts); w2 = eval_W(Mo, Yo, ts); w12 = eval_W(Xo, Yo, ts)
    j_out = float(np.min(4*w1*w2 / (w12*w12 + 1e-300)))
    if j_out < 1.0 - 1e-6 or np.min(w1*w2) < -1e-7: Lm_fails += 1
    if j_out < j_in - 1e-5: Lm_J_drop_cnt += 1
    worst_J_out = min(worst_J_out, j_out)
    max_J_drop = max(max_J_drop, j_in - j_out)

print(f"  NON-FACTORED 6x6 determinantal real-stable quads in (-1, 0): valid = {valid}")
print(f"    (a) P(t, (1+t)y) real-stable fails: {p1t_fails}/{valid}")
print(f"    (b) L_m = (D_m, D_{{m+1}}, D_{{m+2}}) fails (J<1): {Lm_fails}/{valid}, worst J_out = {worst_J_out:.6f}")
print(f"        J_out < J_in count: {Lm_J_drop_cnt}/{valid}, max(J_in - J_out) = {max_J_drop:+.6f}")

print("\n=== 2. GLZ COMPARISON: T_old -> T_Lm (exact L_{k+1}) vs T_GLZ (D_{k+2}^2) vs T_alt (D_{k+1} D_{k+3}) ===")
mp.mp.dps = 50
def apply_D_int(m, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly): res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv_int(p): return [r*p[r] for r in range(1, len(p))]
def eval_mp(p, tv):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * tv + c
    return acc

KMAX = 12
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
U_all = [None]*(KMAX+1)
for k in range(2, KMAX+1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D_int(m, H[k-1][0])
    for i in range(1, k): U[i] = poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    U_all[k] = U
    for i in range(k+1): H[k][i] = apply_D_int(m, U[i])

ts_mp = [mp.mpf(s) for s in ["-10", "-2", "-1.1", "-1", "-0.95", "-0.5", "-0.1", "-0.01", "0", "0.01", "0.05", "0.2", "1", "10"]]
def min_J_mp(X, M, Y):
    Xp, Mp, Yp = deriv_int(X), deriv_int(M), deriv_int(Y)
    mn = mp.mpf('1e100')
    for tv in ts_mp:
        xv, mv, yv = eval_mp(X, tv), eval_mp(M, tv), eval_mp(Y, tv)
        xp, mpv, yp = eval_mp(Xp, tv), eval_mp(Mp, tv), eval_mp(Yp, tv)
        w1 = xp*mv - xv*mpv; w2 = mpv*yv - mv*yp; w12 = xp*yv - xv*yp
        if w12 != 0:
            r = 4*w1*w2 / (w12*w12)
            if r < mn: mn = r
    return float(mn)

for k in range(3, 12):
    m = k + 1
    i = k // 2 + 1 if k // 2 + 1 <= k else k
    U0, U1, U2 = U_all[k][i-2], U_all[k][i-1], U_all[k][i]
    T_old = (U0, apply_D_int(m, U1), apply_D_int(m, apply_D_int(m, U2)))
    T_Lm  = (apply_D_int(m, T_old[0]), apply_D_int(m+1, T_old[1]), apply_D_int(m+2, T_old[2]))
    T_GLZ = (H[k][i-2], apply_D_int(m+1, H[k][i-1]), apply_D_int(m+1, apply_D_int(m+1, H[k][i])))
    T_alt = (H[k][i-2], apply_D_int(m+1, H[k][i-1]), apply_D_int(m, apply_D_int(m+2, H[k][i])))
    j_old = min_J_mp(*T_old); j_lm = min_J_mp(*T_Lm); j_glz = min_J_mp(*T_GLZ); j_alt = min_J_mp(*T_alt)
    print(f"  k={k:2d}, i={i:2d}: T_old={j_old:.5f} -> T_Lm={j_lm:.5f} | T_GLZ={j_glz:.5f} | T_alt={j_alt:.5f}  (GLZ-Lm={j_glz-j_lm:+.5f})")
