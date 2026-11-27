"""Fast diagnostic of the TWO exact closure links for the Plucker Quad route + q<=1/4 stability check:
1. Check Q_{k, q}(t, y) = sum q^{i(i-1)/2} htilde_{k,i}(t) y^i Rayleigh at t=0 for q = 1/4, 1/8, 1/100 up to k = 40!
   (Wait: earlier we tested single-step Herglotz residue for q=1/2; what about Q_{k, q} ITSELF at q = 1/4?!)
2. For the D-lifted quad route (X, M, Y) = (htilde_{k-1,i-2}, D_{k+1} htilde_{k-1,i-1}, D_{k+1}^2 htilde_{k-1,i}):
   Test LINK 1 (index shift D_{k+1} -> D_{k+2}):
   Compare J_ratio = 4 W(X,M) W(M,Y) / W(X,Y)^2 BEFORE and AFTER shifting m = k+1 -> m = k+2!
   Wait: DOES shifting m -> m+1 ALWAYS IMPROVE (or preserve >= 1) J_ratio?!
   And DOES applying D_m ALWAYS IMPROVE (or preserve) the Plucker ratio R = W_{AB} W_{CE} / W_{BC}^2?!
   Test on BOTH GLZ (k=3..14) AND random real-stable quads (A, B, C, E)!!
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 50

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

def deriv_int(p):
    return [r*p[r] for r in range(1, len(p))]

def eval_mp(p, tv):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * tv + c
    return acc

print("=== 1. TEST Q_{k, q}(t, y) = sum q^{i(i-1)/2} htilde_{k,i}(t) y^i AT t=0 FOR q=1/2, 1/4, 1/8, 1/100 (k=2..40) ===")
KMAX = 40
a0 = [[0]*(k+1) for k in range(KMAX+1)]
a1 = [[0]*(k+1) for k in range(KMAX+1)]
a0[1][0], a1[1][0] = 1, 2; a0[1][1], a1[1][1] = 1, 0
for k in range(2, KMAX+1):
    m = k + 1
    for i in range(k+1):
        u0 = (a0[k-1][i] if i<=k-1 else 0) + (a0[k-1][i-1] if i>=1 else 0)
        u1 = ((2*a1[k-1][i] + m*a0[k-1][i]) if i<=k-1 else 0) + (a1[k-1][i-1] if i>=1 else 0)
        a0[k][i] = u0; a1[k][i] = 2*u1 + m*u0

for q_str in ["0.5", "0.25", "0.125"]:
    q = mp.mpf(q_str)
    fails = []
    for k in range(2, 20):
        mp.mp.dps = 120
        p0 = [q**(i*(i-1)//2) * a0[k][i] for i in range(k+1)]
        p1 = [q**(i*(i-1)//2) * a1[k][i] for i in range(k+1)]
        # Initial guesses for roots of sum q^{i(i-1)/2} C(k,i) y^i: y_j ~ -q^{-(j-1)} * j / (k-j+1)
        rts = []
        for j in range(1, k + 1):
            yj = - (q ** (-(j - 1))) * mp.mpf(j) / mp.mpf(k - j + 1)
            for _ in range(30):
                v0 = mp.mpf(0); dv0 = mp.mpf(0)
                for idx in range(k, -1, -1):
                    if idx < k: dv0 = dv0 * yj + v0
                    v0 = v0 * yj + p0[idx]
                if dv0 == 0: break
                step = v0 / dv0
                yj -= step
                if abs(step) < mp.mpf('1e-50') * abs(yj): break
            rts.append(yj)
        # Verify distinctness!
        rts.sort()
        assert all(rts[m+1] - rts[m] > mp.mpf('1e-30') * abs(rts[m]) for m in range(k - 1)), f"dup root at k={k}, q={q_str}"
        signs = []
        for yj in rts:
            v0 = mp.mpf(0); dv0 = mp.mpf(0); v1 = mp.mpf(0)
            for idx in range(k, -1, -1):
                if idx < k: dv0 = dv0 * yj + v0
                v0 = v0 * yj + p0[idx]; v1 = v1 * yj + p1[idx]
            prod = dv0 * v1
            signs.append(1 if prod > 0 else (-1 if prod < 0 else 0))
        if len(set(signs)) > 1 or 0 in signs:
            fails.append((k, signs.count(1), signs.count(-1)))
            if len(fails) >= 3: break
    print(f"  q = {q_str:5s}: Rayleigh at t=0 fails (k=2..19) = {fails}", flush=True)
mp.mp.dps = 50

print("\n=== 2. TEST HOW D_m PRESERVES (OR BREAKS!) PLUCKER RATIO W_{AB} W_{CE} / W_{BC}^2 ON RANDOM REAL-STABLE QUADS ===")
# Take random equal-degree real-stable 4-tuples (A, B, C, E) on (-1, 0) where
# (A, B, C) and (B, C, E) are real-stable (J >= 0 on R) AND W_{AB} W_{CE} / W_{BC}^2 <= rho_in!
# Compute max over t in R of W(DA, DB) W(DC, DE) / W(DB, DC)^2!
def poly_from_roots(rts):
    p = np.array([1.0])
    for r in rts:
        p = np.convolve(p, [-r, 1.0])  # (1 - t/r)? No, (-r + t) => at t=0 positive if r < 0!
    return p  # coeffs in ascending powers of t

def apply_D_flt(m, p):
    res = np.zeros(len(p) + 1)
    for r, c in enumerate(p):
        res[r] += c * (r + 1)
        res[r + 1] += c * (r + m)
    return res

def eval_W(f, g, ts):
    fp = np.arange(1, len(f)) * f[1:]
    gp = np.arange(1, len(g)) * g[1:]
    fv = np.polyval(f[::-1], ts); gv = np.polyval(g[::-1], ts)
    fpv = np.polyval(fp[::-1], ts); gpv = np.polyval(gp[::-1], ts)
    return fpv * gv - fv * gpv

rng = np.random.default_rng(77)
ts_grid = np.concatenate([-10**np.linspace(2, -4, 40), np.linspace(-0.99, -0.01, 100), [0.0], 10**np.linspace(-4, 2, 40)])
max_amplification = 0.0
worse_cases = 0; valid = 0
for _ in range(300):
    d = 4
    # Generate 4 strictly interlacing root sets in (-1, 0): a_j < b_j < c_j < e_j < a_{j+1}
    pts = -np.sort(rng.uniform(0.01, 0.99, 4 * d))[::-1]  # ascending in (-1, 0)
    ra, rb, rc, re = pts[0::4], pts[1::4], pts[2::4], pts[3::4]
    A = poly_from_roots(ra); B = poly_from_roots(rb); C = poly_from_roots(rc); E = poly_from_roots(re)
    # Check J >= 0 on grid for (A,B,C) and (B,C,E)
    wab = eval_W(A, B, ts_grid); wbc = eval_W(B, C, ts_grid); wce = eval_W(C, E, ts_grid)
    wac = eval_W(A, C, ts_grid); wbe = eval_W(B, E, ts_grid)
    if np.min(4*wab*wbc - wac**2) < 0 or np.min(4*wbc*wce - wbe**2) < 0:
        continue
    valid += 1
    rho_in = float(np.max((wab * wce) / (wbc**2 + 1e-300)))
    m = float(10 ** rng.uniform(0.1, 1.5)) + 1.0
    DA, DB, DC, DE = apply_D_flt(m, A), apply_D_flt(m, B), apply_D_flt(m, C), apply_D_flt(m, E)
    dwab = eval_W(DA, DB, ts_grid); dwbc = eval_W(DB, DC, ts_grid); dwce = eval_W(DC, DE, ts_grid)
    rho_out = float(np.max((dwab * dwce) / (dwbc**2 + 1e-300)))
    if rho_out > rho_in + 1e-9:
        worse_cases += 1
    if rho_out / max(rho_in, 1e-12) > max_amplification:
        max_amplification = rho_out / max(rho_in, 1e-12)
print(f"  Random real-stable quads (valid={valid}): D_m increased max(W_AB W_CE / W_BC^2) in {worse_cases}/{valid}, max amp={max_amplification:.4f}")

print("\n=== 3. GLZ STEP-BY-STEP DIAGNOSTIC OF LINK 1 (m=k+1 -> k+2) AND LINK 2 (PLUCKER RATIO) ===")
H = [[None]*(k+1) for k in range(13)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, 13):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D_int(m, H[k-1][0])
    for i in range(1, k): U[i] = poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1): H[k][i] = apply_D_int(m, U[i])

ts_check = [mp.mpf(s) for s in ["-5", "-1.1", "-1", "-0.99", "-0.8", "-0.5", "-0.2", "-0.05", "-0.005", "0", "0.01", "0.5", "5"]]
def min_J_ratio(X, M, Y):
    Xp, Mp, Yp = deriv_int(X), deriv_int(M), deriv_int(Y)
    mn = mp.mpf('1e100')
    for tv in ts_check:
        xv, mv, yv = eval_mp(X, tv), eval_mp(M, tv), eval_mp(Y, tv)
        xp, mpv, yp = eval_mp(Xp, tv), eval_mp(Mp, tv), eval_mp(Yp, tv)
        wxm = xp*mv - xv*mpv; wmy = mpv*yv - mv*yp; wxy = xp*yv - xv*yp
        if wxy == 0: continue
        r = 4 * wxm * wmy / (wxy * wxy)
        if r < mn: mn = r
    return float(mn)

for k in range(3, 12):
    # Compare J_ratio at row k-1 with D_k (previous m!) vs D_{k+1} (shifted m!)
    mn_old = 1e9; mn_new = 1e9; min_diff = 1e9
    for i in range(2, k):
        X = H[k-1][i-2]
        M_old = apply_D_int(k, H[k-1][i-1]); Y_old = apply_D_int(k, apply_D_int(k, H[k-1][i]))
        M_new = apply_D_int(k+1, H[k-1][i-1]); Y_new = apply_D_int(k+1, apply_D_int(k+1, H[k-1][i]))
        r_old = min_J_ratio(X, M_old, Y_old)
        r_new = min_J_ratio(X, M_new, Y_new)
        mn_old = min(mn_old, r_old); mn_new = min(mn_new, r_new)
        min_diff = min(min_diff, r_new - r_old)
    print(f"  k={k:2d}: min J_ratio with D_k = {mn_old:.6f}, with D_{{k+1}} = {mn_new:.6f}, min(new-old) = {min_diff:+.6f}")
