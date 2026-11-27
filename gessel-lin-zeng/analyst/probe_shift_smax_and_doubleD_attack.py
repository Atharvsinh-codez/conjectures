"""Two crucial diagnostics:
1. For Q_m = (htilde_{k-1,i-2}, D_m htilde_{k-1,i-1}, D_m^2 htilde_{k-1,i}), at m = k + 1 + s:
   Find the EXACT critical s_max(k) where min_{t in R} J_ratio(Q_{k+1+s}) first drops below 1!
   (Remember: GLZ only needs s = 1! How big is s_max(k) compared to 1?!)
2. Targeted adversarial attack on GLZ-shaped double-operator product prod_{j=1}^k D_{m_j}(D_{m_j} + y)(1)
   with arbitrary m_j > 1 (including m_j = 1 + 1e-4, 1e6, greedy optimization to minimize Turan on (-1, 0)!)
   Can we break Turan or real-rootedness up to k = 12, or does outer D_{m_j} stop @algebraist's k=9 attack?!
3. Check whether for ANY real-stable cubic P(t, y) = A + 3 B y + 3 C y^2 + E y^3 (equal deg, roots in (-1,0)),
   the Plucker ratio R = W_{AB} W_{CE} / W_{BC}^2 is bounded above (by 1? by 9? unbounded?)!
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 50

def apply_D_mp(m, poly):
    res = [mp.mpf(0)]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def apply_D_int(m, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv(p): return [r*p[r] for r in range(1, len(p))]

def eval_mp(p, tv):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * tv + c
    return acc

KMAX = 14
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D_int(m, H[k-1][0])
    for i in range(1, k): U[i] = poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1): H[k][i] = apply_D_int(m, U[i])

print("=== 1. EXACT s_max(k) WHERE J_ratio(Q_{k+1+s}) FIRST DROPS BELOW 1 ON R (k = 3..14) ===")
ts_probe = [mp.mpf(s) for s in ["-10", "-2", "-1.1", "-1", "-0.99", "-0.9", "-0.5", "-0.1", "-0.01", "-1e-5", "0", "1e-6", "1e-4", "1e-3", "0.01", "0.05", "0.2", "1", "10"]]

def min_J_for_s(k, i, s_val):
    m = mp.mpf(k + 1) + mp.mpf(s_val)
    X = [mp.mpf(c) for c in H[k-1][i-2]]
    M = apply_D_mp(m, [mp.mpf(c) for c in H[k-1][i-1]])
    Y = apply_D_mp(m, apply_D_mp(m, [mp.mpf(c) for c in H[k-1][i]]))
    Xp, Mp, Yp = deriv(X), deriv(M), deriv(Y)
    mn = mp.mpf('1e100'); worst_t = None
    for tv in ts_probe:
        xv, mv, yv = eval_mp(X, tv), eval_mp(M, tv), eval_mp(Y, tv)
        xp, mpv, yp = eval_mp(Xp, tv), eval_mp(Mp, tv), eval_mp(Yp, tv)
        wxm = xp*mv - xv*mpv; wmy = mpv*yv - mv*yp; wxy = xp*yv - xv*yp
        if wxy == 0: continue
        r = 4 * wxm * wmy / (wxy * wxy)
        if r < mn: mn = r; worst_t = tv
    return float(mn), float(worst_t)

for k in [4, 6, 8, 10, 12, 14]:
    i = k // 2 + 1 if k // 2 + 1 < k else k - 1
    res = []
    for s_val in [0, 1, 5, 20, 100, 1000, 10**6]:
        mn, wt = min_J_for_s(k, i, s_val)
        res.append(f"s={s_val}:{mn:.4f}(t={wt:g})")
    print(f"  k={k:2d}, i={i:2d}: " + ", ".join(res), flush=True)

print("\n=== 2. TARGETED ADVERSARIAL ATTACK ON prod_{j=1}^k D_{m_j}(D_{m_j} + y)(1) (m_j > 1) ===")
# Use float64 normalized polynomials + greedy search over m_j in {1.0001, 1.01, 1.5, 2, 5, 20, 100, 1000, 1e5}
# to try to drive Turan f_i^2 - f_{i-1} f_{i+1} < 0 on (-1, 0) OR drive J < 0!
def apply_D_flt(m, p):
    res = np.zeros(len(p) + 1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def step_double_D(polys, m):
    k = len(polys)
    U = [np.zeros(1) for _ in range(k + 1)]
    for i, p in enumerate(polys):
        dp = apply_D_flt(m, p)
        if len(U[i]) < len(dp): U[i] = np.pad(U[i], (0, len(dp) - len(U[i])))
        U[i][:len(dp)] += dp
        if len(U[i+1]) < len(p): U[i+1] = np.pad(U[i+1], (0, len(p) - len(U[i+1])))
        U[i+1][:len(p)] += p
    return [apply_D_flt(m, u) for u in U]

ts_neg = np.concatenate([-10**np.linspace(-9, -1, 40), np.linspace(-0.12, -0.999, 60)])

def min_turan_rel(polys):
    """Min over i=1..k-1 and t in (-1,0) of (f_i^2 - f_{i-1}f_{i+1}) / (f_i^2 + |f_{i-1}f_{i+1}|)."""
    k = len(polys) - 1
    vals = [np.polyval(p[::-1], ts_neg) for p in polys]
    mn = 1e9; info = None
    for i in range(1, k):
        num = vals[i]**2 - vals[i-1]*vals[i+1]
        den = vals[i]**2 + np.abs(vals[i-1]*vals[i+1]) + 1e-300
        rat = num / den
        idx = int(np.argmin(rat))
        if rat[idx] < mn:
            mn = float(rat[idx]); info = (i, float(ts_neg[idx]))
    return mn, info

# Try greedy minimization + extreme bimodal sequences up to k = 12!
cands = [1.0001, 1.01, 1.1, 2.0, 10.0, 100.0, 1e4, 1e6]
polys = [np.array([1.0])]
chosen_ms = []
for step in range(1, 11):
    best_m = None; best_score = 1e9; best_polys = None; best_info = None
    for mv in cands:
        cand_p = step_double_D(polys, mv)
        if step == 1:
            score, info = 1.0, None
        else:
            score, info = min_turan_rel(cand_p)
        if score < best_score:
            best_score = score; best_m = mv; best_polys = cand_p; best_info = info
    # Also try starting from m=1.0001 then 1e6
    polys = best_polys; chosen_ms.append(best_m)
    if step >= 2:
        print(f"  Greedy step k={step:2d}: chosen m={best_m:g}, min_rel_Turan={best_score:+.6f} at {best_info}")

# Also test explicit extreme patterns mimicking @algebraist's attack: tiny (1+eps) vs huge (1e6)!
for label, ms in [
    ("[1.0001]*5 + [1e6]*5", [1.0001]*5 + [1e6]*5),
    ("[1e6]*5 + [1.0001]*5", [1e6]*5 + [1.0001]*5),
    ("[1.0001, 1e6]*5", [1.0001, 1e6]*5),
    ("[1.0001]*8 + [1e6]*2", [1.0001]*8 + [1e6]*2),
    ("[1e6]*2 + [1.0001]*8", [1e6]*2 + [1.0001]*8),
]:
    p = [np.array([1.0])]
    for mv in ms: p = step_double_D(p, mv)
    sc, inf = min_turan_rel(p)
    print(f"  Pattern {label:25s}: min_rel_Turan = {sc:+.6f} at {inf}", flush=True)

print("\n=== 3. PLUCKER RATIO R = W_AB W_CE / W_BC^2 FOR REAL-STABLE CUBIC P(t,y) = A + 3By + 3Cy^2 + Ey^3 ===")
# Generate random real-stable cubics P(t, y) = prod_{j=1}^3 (a_j + b_j t + c_j y) * (or higher deg in t: sum of rank-1?
# Better: P(t, y) = det(I + t*M_1 + y*M_2) where M_1, M_2 are PSD 3x3 matrices! EVERY real-stable cubic in y of deg 3 in t
# has a Helton-Vinnikov determinantal representation det(A_0 + t A_1 + y A_2)!!)
rng = np.random.default_rng(2026)
ts_r = np.linspace(-5, 5, 200)
max_R_unnorm = 0.0; max_R_binom = 0.0
for _ in range(300):
    # Random 3x3 PSD M1, M2
    X1 = rng.normal(0, 1, (3, 3)); M1 = X1 @ X1.T + 0.01*np.eye(3)
    X2 = rng.normal(0, 1, (3, 3)); M2 = X2 @ X2.T + 0.01*np.eye(3)
    # Extract coeffs of y^0, y^1, y^2, y^3 as polys in t (deg 3):
    # Evaluate at 4 y-values and 4 t-values? Even easier: det(I + t M1 + y M2) at t = t_grid
    # Wait: compute A(t), B(t), C(t), E(t) where det(I + t M1 + y M2) = A + 3 B y + 3 C y^2 + E y^3!
    # Use exact interpolation in y at y = 0, 1, -1, 2 for each t, and numerical deriv in t (or exact polyfit deg 3)!
    t_pts = np.array([-2.0, -1.0, 0.0, 1.0])
    y_pts = np.array([0.0, 1.0, -1.0, 2.0])
    V_y = np.vander(y_pts, 4, increasing=True)
    V_t = np.vander(t_pts, 4, increasing=True)
    grid = np.zeros((4, 4))
    for it, tv in enumerate(t_pts):
        for iy, yv in enumerate(y_pts):
            grid[it, iy] = np.linalg.det(np.eye(3) + tv * M1 + yv * M2)
    # Solve coeffs[r, i] for t^r y^i
    C_mat = np.linalg.solve(V_t, np.linalg.solve(V_y, grid.T).T)
    # Unnormalized: p0, p1, p2, p3 = C_mat[:, 0..3]; binom-normalized: A=p0, B=p1/3, C=p2/3, E=p3
    p = [C_mat[:, i] for i in range(4)]
    def eval_W_p(f, g):
        fp = np.arange(1, 4)*f[1:]; gp = np.arange(1, 4)*g[1:]
        return np.polyval(fp[::-1], ts_r)*np.polyval(g[::-1], ts_r) - np.polyval(f[::-1], ts_r)*np.polyval(gp[::-1], ts_r)
    w01 = eval_W_p(p[0], p[1]); w12 = eval_W_p(p[1], p[2]); w23 = eval_W_p(p[2], p[3])
    if np.min(np.abs(w12)) < 1e-9: continue
    r_un = float(np.max((w01 * w23) / (w12 * w12)))
    # Notice for binom-normalized (p0, p1/3, p2/3, p3): W_AB*W_CE / W_BC^2 = (1/3 * 1/3) / (1/9)^2 * r_un = 9 * r_un!
    max_R_unnorm = max(max_R_unnorm, r_un)
    max_R_binom = max(max_R_binom, 9.0 * r_un)
print(f"  Random 3x3 determinantal real-stable cubics (300 trials): max R(p0,p1,p2,p3)={max_R_unnorm:.4f}, max R(binom)={max_R_binom:.4f}")
