import numpy as np
import mpmath as mp

mp.mp.dps = 60

def apply_D_int(m, p):
    res = [0]*(len(p) + 1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv_int(p):
    return [r*p[r] for r in range(1, len(p))] if len(p) > 1 else [0]

def shift_t(p): return [0] + p

def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v

def eval_W_mp(f, g, t):
    return eval_mp(f, t)*eval_mp(deriv_int(g), t) - eval_mp(deriv_int(f), t)*eval_mp(g, t)

# Build GLZ tilde_h up to k=18
KMAX = 18
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    H[k][0] = apply_D_int(m, apply_D_int(m, H[k-1][0]))
    for i in range(1, k):
        H[k][i] = apply_D_int(m, poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1]))
    H[k][k] = apply_D_int(m, H[k-1][k-1])

print("=== 1. POLE FUNCTION r_i(t) = -W(h_{i-1}, D_0 h_i) / W(h_{i-1}, t h_i) AND DIFFERENCE r_{i+1}(t) - r_i(t) ON GLZ ===")
# Notice W(h_{i-1}, t h_i) > 0 EVERYWHERE on R, so r_i(t) is smooth and finite on ALL R!
ts_grid = [mp.mpf(s) for s in ["-100", "-10", "-2", "-1.1", "-0.99", "-0.9", "-0.7", "-0.5", "-0.3", "-0.1", "-0.01", "-1e-5", "0", "1e-5", "0.01", "0.1", "0.5", "1", "10", "100"]]

for k in [5, 8, 11, 14, 17]:
    j0 = mp.mpf(k + 1)
    max_r = -mp.mpf('1e100'); min_r = mp.mpf('1e100')
    max_diff = -mp.mpf('1e100'); min_diff = mp.mpf('1e100')
    max_abs_ddelta = mp.mpf(0)
    worst_diff_loc = None
    # Precompute r_i(t) for i = 1..k
    Rvals = {}
    for i in range(1, k + 1):
        h0, h1 = H[k][i-1], H[k][i]
        D0h1 = apply_D_int(0, h1)
        th1 = shift_t(h1)
        Rvals[i] = []
        for tv in ts_grid:
            rv = -eval_W_mp(h0, D0h1, tv) / eval_W_mp(h0, th1, tv)
            Rvals[i].append(rv)
            if rv > max_r: max_r = rv
            if rv < min_r: min_r = rv
    for i in range(1, k):
        for idx, tv in enumerate(ts_grid):
            diff = Rvals[i+1][idx] - Rvals[i][idx]
            if diff > max_diff: max_diff = diff
            if diff < min_diff:
                min_diff = diff
                worst_diff_loc = (i, float(tv), float(Rvals[i][idx]), float(Rvals[i+1][idx]))
            dd = abs(1/(j0 - Rvals[i+1][idx]) - 1/(j0 - Rvals[i][idx]))
            if dd > max_abs_ddelta: max_abs_ddelta = dd
    print(f"  k={k:2d} (j={k+1:2d}): r_i in [{float(min_r):+7.3f}, {float(max_r):+7.3f}] | r_{{i+1}}-r_i in [{float(min_diff):+7.4f}, {float(max_diff):+7.4f}] | j^2*|d_delta| <= {float(j0*j0*max_abs_ddelta):.4f}")
    print(f"            worst min(r_{{i+1}}-r_i) at i={worst_diff_loc[0]}, t={worst_diff_loc[1]}: r_i={worst_diff_loc[2]:.3f}, r_{{i+1}}={worst_diff_loc[3]:.3f}")

print("\n=== 2. TEST WHETHER BOUNDING |r_{i+1}(t) - r_i(t)| <= 1 ON RANDOM d=2 TRIPLES FORCES O(1/j^2) STAGE 3 DROP!! ===")
# Recall in task-526, d=2 triple (h0, h1, h2) with r* <= 0.75 j STILL dropped by +0.15/j!
# Now add the constraint |r_{12}(t) - r_{01}(t)| <= 1.0 (or 0.5) on R (and r* <= 0.75 j) and measure j^2 * drop!
def apply_D_np(p, m):
    res = np.zeros(len(p)+1)
    for n, c in enumerate(p):
        res[n] += c*(n+1); res[n+1] += c*(n+m)
    return res

def wronskian_np(f, g):
    fp = np.arange(1, len(f))*f[1:] if len(f) > 1 else np.array([0.0])
    gp = np.arange(1, len(g))*g[1:] if len(g) > 1 else np.array([0.0])
    a = np.convolve(f, gp); b = np.convolve(fp, g)
    n = max(len(a), len(b))
    res = np.zeros(n); res[:len(a)] += a; res[:len(b)] -= b
    return res

def eval_np(p, t):
    v = np.zeros_like(t)
    for c in p[::-1]: v = v*t + c
    return v

def from_roots(rs):
    p = np.array([1.0])
    for r in rs: p = np.convolve(p, [-r, 1.0])
    return p

ts_np = np.concatenate([-np.geomspace(1e-3, 1e2, 40), np.linspace(-0.999, -0.001, 120), np.geomspace(1e-3, 1e2, 40)])
rng = np.random.default_rng(42)

for j in [20, 40, 80, 160]:
    n_gap = 0; max_d_gap = -1e9
    n_matched = 0; max_d_matched = -1e9
    for _ in range(4000):
        # d=2: h2 = [1], h1 = deg 1, h0 = deg 2, or d=3: h2 = deg 1, h1 = deg 2, h0 = deg 3!
        # Let's do h2 deg 1, h1 deg 2, h0 deg 3 (all roots in (-1, 0), interlacing h2 << h1 << h0)!
        r0 = np.sort(-rng.uniform(0.01, 0.99, size=3))
        r1 = np.array([rng.uniform(r0[0], r0[1]), rng.uniform(r0[1], r0[2])])
        r2 = np.array([rng.uniform(r1[0], r1[1])])
        h0, h1, h2 = from_roots(r0), from_roots(r1), from_roots(r2)
        # Compute r01(t) and r12(t)!
        D0h1 = apply_D_np(h1, 0); th1 = np.concatenate([[0.0], h1])
        D0h2 = apply_D_np(h2, 0); th2 = np.concatenate([[0.0], h2])
        r01 = -eval_np(wronskian_np(h0, D0h1), ts_np) / eval_np(wronskian_np(h0, th1), ts_np)
        r12 = -eval_np(wronskian_np(h1, D0h2), ts_np) / eval_np(wronskian_np(h1, th2), ts_np)
        rmax = max(float(np.max(r01)), float(np.max(r12)))
        if rmax > 0.75 * j: continue
        # Scale h1 so Jin ~ 1 + 4/j
        M_j = apply_D_np(h1, j); Y_j = apply_D_np(apply_D_np(h2, j), j)
        w1 = eval_np(wronskian_np(h0, M_j), ts_np)
        w2 = eval_np(wronskian_np(M_j, Y_j), ts_np)
        w12 = eval_np(wronskian_np(h0, Y_j), ts_np)
        if np.min(w1) <= 0 or np.min(w2) <= 0: continue
        Jin_raw = float(np.min(4*w1*w2 / (w12*w12 + 1e-300)))
        if Jin_raw <= 0: continue
        c = np.sqrt((1.0 + 4.0/j) / Jin_raw)
        M_j = M_j * c
        M_jp = apply_D_np(h1, j+1) * c; Y_jp = apply_D_np(apply_D_np(h2, j+1), j+1)
        w1p = eval_np(wronskian_np(h0, M_jp), ts_np)
        w2p = eval_np(wronskian_np(M_jp, Y_jp), ts_np)
        w12p = eval_np(wronskian_np(h0, Y_jp), ts_np)
        if np.min(w1p) <= 0 or np.min(w2p) <= 0: continue
        Jin = 1.0 + 4.0/j
        Jout = float(np.min(4*w1p*w2p / (w12p*w12p + 1e-300)))
        drop = Jin - Jout
        n_gap += 1
        max_d_gap = max(max_d_gap, drop)
        if np.max(np.abs(r12 - r01)) <= 1.0:
            n_matched += 1
            max_d_matched = max(max_d_matched, drop)
    print(f"  j={j:3d}: [r*<.75j n={n_gap:4d}] j*d={j*max_d_gap:+.4f}, j^2*d={j*j*max_d_gap:+7.2f} | [|r12-r01|<=1 n={n_matched:4d}] j*d={j*max_d_matched:+.4f}, j^2*d={j*j*max_d_matched:+7.2f}")

print("\n=== 3. CHECK MAX SINGLE-STEP PHASE ADVANCE Delta_i(t) OF E_{k,i} AT t=-0.95..-0.999 FOR k=10..28! ===")
# Does W(E_{k,i}, E_{k,i+r}) > 0 ALWAYS hold for r = 3, 4 as k grows?!
# Measure max_{i, t} Delta_i(t)/pi where Delta_i = atan2(W(E_i, E_{i+1}), E_i E_{i+1} + E'_i E'_{i+1})?
# Better: scale-invariant phase! Just check directly for r=3, 4, 5, 6, 7 at k=26, 28, 30? (Fast via float64? No, use mpmath at k=18 and compute exact angle!)
# Actually let's compute exact phase step using normalized E_i at k=6, 10, 14, 18!
for k in [6, 10, 14, 18]:
    m = k + 1
    E = []
    for i in range(k):
        p = H[k-1][i]
        for _ in range(i): p = apply_D_int(m, p)
        E.append(p)
    # Find smallest r that fails W(E_i, E_{i+r}) > 0 on ts_grid
    min_fail_r = None
    for r in range(1, k):
        for i in range(k - r):
            for tv in [mp.mpf('-0.999'), mp.mpf('-0.99'), mp.mpf('-0.95'), mp.mpf('-0.9'), mp.mpf('-0.8'), mp.mpf('-0.5')]:
                if eval_W_mp(E[i], E[i+r], tv) <= 0:
                    min_fail_r = (r, i, float(tv))
                    break
            if min_fail_r: break
        if min_fail_r: break
    print(f"  k={k:2d}: smallest failing r for W(E_i, E_{{i+r}}) <= 0: {min_fail_r}")
