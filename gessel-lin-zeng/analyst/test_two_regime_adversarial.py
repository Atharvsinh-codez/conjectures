"""
CRITICAL ADVERSARIAL TEST OF TWO-REGIME CLOSURE!
Check:
1. On GLZ k=50, 70: densely scan around ALL roots in [-1/(k+2), 0) to see the EXACT
   set where J_raw(t) < 1, and whether B(t) = t(1+t)A'(t) is EVER negative where J_raw(t) < 1!
2. On GENERAL triples (g << f << h with deg g = d+1, deg f = d, deg h = d-1, roots in (-1,0)):
   (a) Does J_raw(g,f,h) >= 1 on (-1, -1/m] imply J_raw(D_m g, D_m f, D_m h) >= 1 on (-1, -1/m]?
   (b) Does it imply J_raw(D_m g, D_m f, D_m h) >= 1 on (-1, -1/(m+1)]?
   (c) What if we ALSO assume A>0, B>=0 whenever J_raw<1 on [-1/m, 0) -- does D_m preserve
       that whole Two-Regime condition AND does the addition step F_i = D_m h_i + h_{i-1} preserve it?!
   Use BOTH random AND adversarial extreme roots (roots clustered near -1, near -1/m, near 0,
   exponentially spaced roots like GLZ, etc.)!
"""
import numpy as np
import mpmath as mp

def eval_Jraw_AB(g, f, h, grid):
    gp, fp, hp = np.polyder(g), np.polyder(f), np.polyder(h)
    gv, fv, hv = np.polyval(g, grid), np.polyval(f, grid), np.polyval(h, grid)
    gd, fd, hd = np.polyval(gp, grid), np.polyval(fp, grid), np.polyval(hp, grid)
    Wgf = gv*fd - gd*fv
    Wfh = fv*hd - fd*hv
    Wgh = gv*hd - gd*hv
    A = fv*fv - gv*hv
    Ap = 2*fv*fd - gd*hv - gv*hd
    B = grid*(1+grid)*Ap
    with np.errstate(divide='ignore', invalid='ignore'):
        J = 4*Wgf*Wfh / (Wgh**2)
    return J, A, B

def apply_Dm_np(m, p):
    c = p[::-1]
    out = np.zeros(len(c) + 1)
    for j, cj in enumerate(c):
        out[j] += (j + 1) * cj
        out[j + 1] += (j + m) * cj
    return out[::-1]

print("=== PART A: Adversarial test of D_m on J_raw >= 1 on (-1, -1/m] ===")
rng = np.random.default_rng(12345)
for mode in ["uniform", "near_zero_exp", "near_split", "near_minus1"]:
    tested = 0
    fail_same = 0    # fails on (-1, -1/m]
    fail_next = 0    # fails on (-1, -1/(m+1)]
    fail_two_reg = 0 # fails full two-regime (J>=1 on RegI, and either J>=1 or B>=0 on RegII)
    min_J_same = 1e9
    min_J_next = 1e9
    worst_case = None
    for trial in range(1500):
        m = int(rng.integers(3, 40))
        d = int(rng.integers(2, 7))
        if mode == "uniform":
            rg = np.sort(-rng.uniform(1e-4, 0.999, d+1))
        elif mode == "near_zero_exp":
            # Exponentially clustered near 0 (like GLZ!)
            base = rng.uniform(1.5, 4.0)
            scale = rng.uniform(0.01, 0.5)
            rg = np.sort(-np.clip(scale * base**(-np.arange(d+1) + rng.uniform(-0.3, 0.3, d+1)), 1e-6, 0.999))
        elif mode == "near_split":
            # Clustered near -1/m
            c0 = -1.0 / m
            spread = rng.uniform(0.01, 0.5) * abs(c0)
            rg = np.sort(np.clip(c0 + rng.uniform(-spread, spread*0.9, d+1), -0.999, -1e-6))
        else: # near_minus1
            rg = np.sort(-1.0 + rng.uniform(1e-4, 0.1, d+1))
        if len(np.unique(rg)) < d+1:
            continue
        rf = np.array([rng.uniform(rg[j] + 0.05*(rg[j+1]-rg[j]), rg[j+1] - 0.05*(rg[j+1]-rg[j])) for j in range(d)])
        rh = np.array([rng.uniform(rf[j] + 0.05*(rf[j+1]-rf[j]), rf[j+1] - 0.05*(rf[j+1]-rf[j])) for j in range(d-1)])
        g = np.poly(rg); g /= np.polyval(g, 0)
        f = np.poly(rf); f /= np.polyval(f, 0)
        h = np.poly(rh); h /= np.polyval(h, 0)
        scale_f = rng.uniform(1.0, 1.3)
        f *= scale_f

        # Build dense grid on (-1, -1/m] and [-1/m, 0)
        gI = np.concatenate([np.linspace(-0.999, -1.0/m, 150),
                             -1.0/m - np.geomspace(1e-6, min(0.1, 1 - 1.0/m), 40)])
        gI = gI[(gI > -1.0) & (gI <= -1.0/m)]
        gII = np.concatenate([np.linspace(-1.0/m, -1e-7, 100),
                              -np.geomspace(1e-7, 1.0/m, 60)])
        Jin_I, Ain_I, Bin_I = eval_Jraw_AB(g, f, h, gI)
        Jin_II, Ain_II, Bin_II = eval_Jraw_AB(g, f, h, gII)
        if np.nanmin(Jin_I) < 1.0 or np.nanmin(Ain_I) <= 0 or np.nanmin(Ain_II) <= 0:
            continue
        # Also check full two-regime input: on gII, wherever J < 1 we have B >= 0
        in_two_reg = np.all((Jin_II >= 1.0) | (Bin_II >= 0.0))
        tested += 1
        gm, fm, hm = apply_Dm_np(m, g), apply_Dm_np(m, f), apply_Dm_np(m, h)
        Jout_same, _, _ = eval_Jraw_AB(gm, fm, hm, gI)
        gI_next = np.concatenate([gI, np.linspace(-1.0/m, -1.0/(m+1), 40)])
        Jout_next, _, _ = eval_Jraw_AB(gm, fm, hm, gI_next)
        ms = float(np.nanmin(Jout_same))
        mn = float(np.nanmin(Jout_next))
        if ms < min_J_same: min_J_same = ms
        if mn < min_J_next:
            min_J_next = mn
            worst_case = (m, d, float(np.nanmin(Jin_I)), ms, mn, in_two_reg)
        if ms < 1.0: fail_same += 1
        if mn < 1.0: fail_next += 1
        if in_two_reg and mn < 1.0:
            fail_two_reg += 1
    print(f"mode={mode:14s}: tested={tested:4d}, fail(-1,-1/m]={fail_same:3d} (min={min_J_same:.5f}), fail(-1,-1/(m+1)]={fail_next:3d} (min={min_J_next:.5f}), fail_from_2reg={fail_two_reg:3d}")
    if worst_case:
        print(f"   worst: m={worst_case[0]}, d={worst_case[1]}, minJin_I={worst_case[2]:.5f}, minJout_same={worst_case[3]:.5f}, minJout_next={worst_case[4]:.5f}, in_2reg={worst_case[5]}")

print("\n=== PART B: Exact GLZ k=50 scan around ALL roots in [-1/(k+2), 0) ===")
mp.mp.dps = 150
def poly_D_mp(m, c):
    out = [mp.mpf(0)] * (len(c) + 1)
    for j, cj in enumerate(c):
        out[j] += (j + 1) * cj
        out[j + 1] += (j + m) * cj
    return out
def poly_add_mp(a, b):
    out = [mp.mpf(0)] * max(len(a), len(b))
    for i, v in enumerate(a): out[i] += v
    for i, v in enumerate(b): out[i] += v
    return out
def eval_mp(c, t):
    v = d = mp.mpf(0)
    for coeff in reversed(c):
        d = d*t + v
        v = v*t + coeff
    return v, d

ht = {(1,0): [mp.mpf(1), mp.mpf(2)], (1,1): [mp.mpf(1)]}
for k in range(2, 51):
    for i in range(k+1):
        t1 = poly_D_mp(k+1, ht[(k-1,i)]) if i <= k-1 else [mp.mpf(0)]
        t2 = ht[(k-1,i-1)] if i >= 1 else [mp.mpf(0)]
        ht[(k,i)] = poly_D_mp(k+1, poly_add_mp(t1, t2))

# Check k=50, i=25 (where J_raw(0)=0.96394 is minimal!) across log-grid from -1/52 to 0!
k, i = 50, 25
m = k + 2
g, f, h = ht[(k, i-1)], ht[(k, i)], ht[(k, i+1)]
s_f = f[1]/f[0]
print(f"k=50, i=25: s_f = {float(s_f):.4e}, -1/s_f = {float(-1/s_f):.4e}, -1/(k+2) = {-1/m:.5f}")
# Scan log-spaced t from -1/m to -1e-25 (and 0)
exps = np.linspace(np.log10(1/m), -25, 400)
t_grid = [-mp.mpf(10)**mp.mpf(e) for e in exps] + [mp.mpf(0)]
min_t_Jlt1 = mp.mpf(0)
max_t_Jlt1 = mp.mpf(-1)
cnt_Bneg_Jlt1 = 0
first_root_approx = None
prev_f = None
for t in reversed(t_grid): # from 0 going left
    gv, gd = eval_mp(g, t)
    fv, fd = eval_mp(f, t)
    hv, hd = eval_mp(h, t)
    if prev_f is not None and fv * prev_f < 0 and first_root_approx is None:
        first_root_approx = float(t)
    prev_f = fv
    Wgf = gv*fd - gd*fv
    Wfh = fv*hd - fd*hv
    Wgh = gv*hd - gd*hv
    J = 4*Wgf*Wfh / (Wgh*Wgh)
    Ap = 2*fv*fd - gd*hv - gv*hd
    B = t*(1+t)*Ap
    if J < 1:
        if t < min_t_Jlt1: min_t_Jlt1 = t
        if t > max_t_Jlt1: max_t_Jlt1 = t
        if B < 0: cnt_Bneg_Jlt1 += 1
print(f"  Approx rightmost root of f: {first_root_approx:.4e}")
print(f"  Region with J_raw < 1: [{float(min_t_Jlt1):.4e}, {float(max_t_Jlt1):.4e}]")
print(f"  Points with J_raw < 1 AND B < 0: {cnt_Bneg_Jlt1}")
