import numpy as np

def apply_D(p, m):
    d = len(p) - 1
    res = np.zeros(d + 2, dtype=np.float64)
    for n in range(d + 1):
        res[n] += (1 + n) * p[n]
        res[n + 1] += (m + n) * p[n]
    return res

def wronskian(f, g):
    # W(f, g) = f g' - f' g  (CORRECT sign: for r(f) < r(g) equal degree monic, W(f, g) > 0!)
    fp = np.arange(1, len(f)) * f[1:] if len(f) > 1 else np.array([0.0])
    gp = np.arange(1, len(g)) * g[1:] if len(g) > 1 else np.array([0.0])
    a = np.convolve(f, gp)
    b = np.convolve(fp, g)
    n = max(len(a), len(b))
    res = np.zeros(n)
    res[:len(a)] += a
    res[:len(b)] -= b
    return res

def eval_poly(p, t):
    v = np.zeros_like(t, dtype=np.float64)
    for c in p[::-1]:
        v = v * t + c
    return v

def min_J_on_R(X, M, Y, ts):
    wxm = eval_poly(wronskian(X, M), ts)
    wmy = eval_poly(wronskian(M, Y), ts)
    wxy = eval_poly(wronskian(X, Y), ts)
    if np.min(wxm) <= 0 or np.min(wmy) <= 0:
        return -1.0
    with np.errstate(divide='ignore', invalid='ignore'):
        J = 4.0 * wxm * wmy / (wxy * wxy)
    J = np.where(np.abs(wxy) < 1e-14 * np.sqrt(wxm * wmy), 1e12, J)
    return float(np.min(J))

ts = np.concatenate([
    -np.geomspace(1e-4, 1e3, 120),
    np.linspace(-0.9999, -0.0001, 400),
    -1.0 - np.geomspace(1e-4, 1e3, 120),
    np.geomspace(1e-4, 1e3, 120),
])

print("=== TEST L_m = (D_m, D_{m+1}, D_{m+2}) WITH CORRECT ORIENTATION r(X) < r(M) < r(Y) ON NON-FACTORED 6x6 DET QUADS ===")
rng = np.random.default_rng(777)
valid = 0; p1t_fails = 0; Lm_fails_below1 = 0; Lm_drop_cnt = 0; max_drop = -1e9; worst_Jout = 1e9
for _ in range(800):
    u = rng.uniform(2.0, 20.0, 4)
    A1 = np.diag([u[0], u[1], u[2], u[3], 0.0, 0.0])
    Q, _ = np.linalg.qr(rng.normal(0, 1, (6, 6)))
    v = rng.uniform(0.2, 5.0, 2)
    A2 = Q @ np.diag([v[0], v[1], 0.0, 0.0, 0.0, 0.0]) @ Q.T
    A0 = np.eye(6)
    t_pts = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    z_pts = np.array([0.0, 1.0, -1.0])
    grid = np.array([[np.linalg.det(A0 + tv*A1 + zv*A2) for zv in z_pts] for tv in t_pts])
    C = np.linalg.solve(np.vander(t_pts, 5, increasing=True), np.linalg.solve(np.vander(z_pts, 3, increasing=True), grid.T).T)
    X, M, Y = C[:, 2], 0.5 * C[:, 1], C[:, 0]
    if abs(X[-1]) < 1e-8 or abs(Y[-1]) < 1e-8: continue
    rx = np.roots(X[::-1]); ry = np.roots(Y[::-1])
    if np.max(np.abs(rx.imag)) > 1e-7 or np.max(np.abs(ry.imag)) > 1e-7: continue
    if np.min(rx.real) <= -0.999 or np.max(rx.real) >= -0.001: continue
    if np.min(ry.real) <= -0.999 or np.max(ry.real) >= -0.001: continue
    # Ensure CORRECT orientation W(X, M) = X M' - X' M > 0 (so r(X) < r(M) < r(Y))!
    if np.min(eval_poly(wronskian(X, M), ts)) <= 0:
        X, Y = Y, X
    Jin = min_J_on_R(X, M, Y, ts)
    if Jin < 1.0 - 1e-6: continue
    # Verify roots really satisfy r(X) < r(Y)
    rx = np.sort(np.roots(X[::-1]).real); ry = np.sort(np.roots(Y[::-1]).real)
    assert all(rx[k] < ry[k] for k in range(len(rx))), f"Orientation mismatch: {rx} vs {ry}"
    valid += 1

    # (a) Test P(t, (1+t)y) -> triple (X, (1+t)M, (1+t)^2 Y)
    Xp = np.pad(X, (0, 2))
    M1 = np.pad(np.convolve(M, [1.0, 1.0]), (0, 1))
    Y2 = np.convolve(Y, [1.0, 2.0, 1.0])
    Jp1t = min_J_on_R(Xp, M1, Y2, ts)
    if Jp1t < Jin - 1e-5:
        p1t_fails += 1

    # (b) Test L_m = (D_m X, D_{m+1} M, D_{m+2} Y) for m in [1.1, 20]
    m = float(rng.uniform(1.1, 20.0))
    Xo, Mo, Yo = apply_D(X, m), apply_D(M, m + 1), apply_D(Y, m + 2)
    Jout = min_J_on_R(Xo, Mo, Yo, ts)
    if Jout < 1.0 - 1e-6: Lm_fails_below1 += 1
    if Jout < Jin - 1e-5: Lm_drop_cnt += 1
    worst_Jout = min(worst_Jout, Jout)
    max_drop = max(max_drop, Jin - Jout)

print(f"  valid = {valid}")
print(f"  (a) (X, (1+t)M, (1+t)^2 Y) J_out < J_in: {p1t_fails}/{valid}")
print(f"  (b) L_m(X, M, Y) J_out < 1: {Lm_fails_below1}/{valid}, worst J_out = {worst_Jout:.6f}")
print(f"      L_m(X, M, Y) J_out < J_in: {Lm_drop_cnt}/{valid}, max(J_in - J_out) = {max_drop:+.6f}")

# (c) ALSO push J_in -> 1.00001 and roots very close to -1 or 0!
print("\n=== EXTREME STRESS TEST OF (X, (1+t)M, (1+t)^2 Y) AND L_m AT J_in ~ 1.0001 AND ROOTS NEAR -1 ===")
valid2 = 0; p1t_fail2 = 0; lm_fail2 = 0; max_drop2 = -1e9; max_drop_p1t = -1e9
for _ in range(3000):
    d = 3
    # Put roots very close to -1 AND very close to 0!
    pts = np.sort(-10**rng.uniform(-4, 0, size=3*d))
    if rng.random() < 0.5:
        pts = np.sort(-1.0 + 10**rng.uniform(-4, 0, size=3*d))
    r0, r1, r2 = pts[0::3], pts[1::3], pts[2::3]
    def mk(rs):
        p = np.array([1.0])
        for r in rs: p = np.convolve(p, [-r, 1.0])
        return p
    X, M, Y = mk(r0), mk(r1), mk(r2)
    # Rescale M so Jin is close to 1.001!
    J0 = min_J_on_R(X, M, Y, ts)
    if J0 <= 0: continue
    c = np.sqrt(1.001 / J0)
    M = M * c
    Jin = min_J_on_R(X, M, Y, ts)
    valid2 += 1
    Xp = np.pad(X, (0, 2)); M1 = np.pad(np.convolve(M, [1.0, 1.0]), (0, 1)); Y2 = np.convolve(Y, [1.0, 2.0, 1.0])
    Jp1t = min_J_on_R(Xp, M1, Y2, ts)
    if Jp1t < Jin - 1e-4: p1t_fail2 += 1
    max_drop_p1t = max(max_drop_p1t, Jin - Jp1t)
    m = 1.01 if rng.random() < 0.5 else 50.0
    Jout = min_J_on_R(apply_D(X, m), apply_D(M, m+1), apply_D(Y, m+2), ts)
    if Jout < Jin - 1e-4: lm_fail2 += 1
    max_drop2 = max(max_drop2, Jin - Jout)
print(f"  extreme d=3: valid={valid2}")
print(f"    (X, (1+t)M, (1+t)^2 Y): J_out < J_in count={p1t_fail2}/{valid2}, max_drop={max_drop_p1t:+.6f}")
print(f"    L_m = (D_m, D_{{m+1}}, D_{{m+2}}): J_out < J_in count={lm_fail2}/{valid2}, max_drop={max_drop2:+.6f}")
