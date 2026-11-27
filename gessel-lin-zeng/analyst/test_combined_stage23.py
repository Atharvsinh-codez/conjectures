"""TEST THE TRUE GLZ STAGE 2+3 COMBINED MAP:
  (U_0, D_j U_1, D_j^2 U_2)  --->  (D_j U_0, D_{j+1} D_j U_1, D_{j+1}^2 D_j U_2) !!

Recall why our previous Stage 3 test lost O(1/j): we tested ARBITRARY (f_0, f_1, f_2), NOT (D_j U_0, D_j U_1, D_j U_2)!
In GLZ, f_r = D_j U_r ALREADY HAS D_j APPLIED TO IT!
Test for d = 2 (U_2 = 1, U_1 = c*(1 + a*t), U_0 = (1+b1*t)(1+b2*t)) AND d = 3 (deg U_2 = 1):
1. Given ONLY J(U_0, D_j U_1, D_j^2 U_2) >= 1 (NO pole-gap assumption at all!!):
   What is max (J_in - J_out) for j = 10, 20, 40, 80, 160?!
2. ALSO test WITH pole-gap r* <= 0.75 * j on (U_0, U_1, U_2)!
3. ALSO test if the combined map (U_0, D_j U_1, D_j^2 U_2) -> (D_j U_0, D_{j+1} D_j U_1, D_{j+1}^2 D_j U_2)
   CAN EVER DROP J BELOW 1 (or if its drop is O(1/j^2)!!)!
"""
import numpy as np

def D(m, p):
    res = np.zeros(len(p) + 1)
    for r, c in enumerate(p):
        res[r] += c * (r + 1)
        res[r + 1] += c * (r + m)
    return res

def W(f, g, ts):
    fp = np.arange(1, len(f)) * f[1:]
    gp = np.arange(1, len(g)) * g[1:]
    return np.polyval(fp[::-1], ts) * np.polyval(g[::-1], ts) - np.polyval(f[::-1], ts) * np.polyval(gp[::-1], ts)

ts = np.concatenate([-10**np.linspace(2, -3, 25), np.linspace(-0.99, -0.01, 50), [0.0], 10**np.linspace(-3, 2, 25)])
rng = np.random.default_rng(2026)

print("=== 1. TRUE GLZ COMBINED STAGE 2+3 MAP AT d = 2: (U0, D_j U1, D_j^2 U2) -> (D_j U0, D_{j+1} D_j U1, D_{j+1}^2 D_j U2) ===")
for j in [10, 20, 40, 80, 160]:
    c_scale = np.sqrt(1.0 + 4.0 / j)
    valid_all = 0; max_drop_all = -1e9; min_Jout_all = 1e9; cfg_all = None
    valid_gap = 0; max_drop_gap = -1e9; min_Jout_gap = 1e9; cfg_gap = None
    for _ in range(20000):
        # U2 = 1, U1 = c*(1 + a*t) with a > 1, U0 = (1-t/x1)(1-t/x2)
        # Try BOTH near-pole a ~ j AND moderate a <= 0.75*j!
        if rng.random() < 0.5:
            a = 1.0 + rng.uniform(0.05, j * 0.99 - 1.0)
        else:
            a = j * rng.uniform(0.85, 0.999)
        U2 = np.array([1.0]); U1 = c_scale * np.array([1.0, a])
        Mj = D(j, U1); Yj = D(j, D(j, U2))
        mr = np.sort(np.roots(Mj[::-1]).real)
        if len(mr) < 2: continue
        u1 = 10**rng.uniform(-4, 0); u2 = 10**rng.uniform(-4, 0)
        x1 = mr[0] - u1*(mr[0] + 1.0)*0.99
        x2 = mr[1] - u2*(mr[1] + 1.0/a)*0.99
        U0 = np.array([1.0, -1.0/x1 - 1.0/x2, 1.0/(x1*x2)])

        w1 = W(U0, Mj, ts); w2 = W(Mj, Yj, ts); w12 = W(U0, Yj, ts)
        if w1.max() >= 0 or w2.max() >= 0: continue
        Jin = float(np.min(4*w1*w2 / (w12*w12 + 1e-300)))
        if Jin < 1.0: continue

        # Output of Stage 2+3:
        Xout = D(j, U0)
        Mout = D(j + 1, Mj)
        Yout = D(j + 1, D(j + 1, D(j, U2)))
        wo1 = W(Xout, Mout, ts); wo2 = W(Mout, Yout, ts); wo12 = W(Xout, Yout, ts)
        if wo1.max() >= 0 or wo2.max() >= 0:
            Jout = -1.0
        else:
            Jout = float(np.min(4*wo1*wo2 / (wo12*wo12 + 1e-300)))
        drop = Jin - Jout
        valid_all += 1
        if drop > max_drop_all: max_drop_all = drop; cfg_all = (round(a,2), round(Jin,5), round(Jout,5))
        if Jout < min_Jout_all: min_Jout_all = Jout

        # Check pole gap at r = 0.75 * j on (U0, U1, U2):
        r = 0.75 * j
        if a < r:
            Mr = D(r, U1); Yr = D(r, D(r, U2))
            if W(U0, Mr, ts).max() < 0 and W(Mr, Yr, ts).max() < 0:
                valid_gap += 1
                if drop > max_drop_gap: max_drop_gap = drop; cfg_gap = (round(a,2), round(Jin,5), round(Jout,5))
                if Jout < min_Jout_gap: min_Jout_gap = Jout

    print(f"  j={j:3d} [ALL   n={valid_all:5d}]: max_drop={max_drop_all:+.6f} (j*d={j*max_drop_all:+.4f}, j^2*d={j*j*max_drop_all:+.2f}), min_Jout={min_Jout_all:.5f}, cfg={cfg_all}")
    print(f"  j={j:3d} [r<.75 n={valid_gap:5d}]: max_drop={max_drop_gap:+.6f} (j*d={j*max_drop_gap:+.4f}, j^2*d={j*j*max_drop_gap:+.2f}), min_Jout={min_Jout_gap:.5f}, cfg={cfg_gap}", flush=True)
