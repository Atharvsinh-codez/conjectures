import numpy as np

# Test Stage 1 (Pascal), Stage 2 (D_j), and Stage 3 (j -> j+1) on:
#   - Pole bound:   P(j) := max_{i, t} j * delta_i(j; t), where delta_i(j; t) = W(h_{i-1}, t h_i) / W(h_{i-1}, D_j h_i)
#   - Lipschitz:    L(j) := max_{i, t} j^2 * |delta_{i+1}(j; t) - delta_i(j; t)|
#   - Stage 3 drop: j^2 * (J(j) - J(j+1)) on arbitrary triples satisfying P(j) <= 3.5 and L(j) <= 5!
#
# CRITICAL QUESTIONS:
#   Q1: Does Stage 2 (D_j) PRESERVE or AMPLIFY P(j) and L(j)?
#   Q2: Does Stage 1 (Pascal) PRESERVE or AMPLIFY P(j) and L(j)?
#   Q3: Under P(j) <= 3.5 and L(j) <= 5 (and J(j) >= 1 + 2/j), is Stage 3 drop j^2*(J(j) - J(j+1)) ALWAYS < 1.0?!

def apply_D(b, p):
    res = np.zeros(len(p)+1)
    for n, c in enumerate(p):
        res[n] += c*(n+1); res[n+1] += c*(n+b)
    return res

def add_p(a, b):
    n = max(len(a), len(b))
    r = np.zeros(n); r[:len(a)] += a; r[:len(b)] += b
    return r

def wronskian(f, g):
    fp = np.arange(1, len(f))*f[1:] if len(f) > 1 else np.array([0.0])
    gp = np.arange(1, len(g))*g[1:] if len(g) > 1 else np.array([0.0])
    a = np.convolve(f, gp); b = np.convolve(fp, g)
    n = max(len(a), len(b))
    res = np.zeros(n); res[:len(a)] += a; res[:len(b)] -= b
    return res

def eval_p(p, t):
    v = np.zeros_like(t)
    for c in p[::-1]: v = v*t + c
    return v

def from_roots(rs, scale=1.0):
    p = np.array([scale], dtype=np.float64)
    for r in rs: p = np.convolve(p, [-r, 1.0])
    return p

ts = np.concatenate([
    -np.geomspace(1e-4, 1e2, 40),
    np.linspace(-0.999, -0.001, 140),
    np.geomspace(1e-4, 1e2, 40),
])

def delta_vec(h_hi, h_lo, j, ts):
    th = np.concatenate([[0.0], h_lo])
    Dh = apply_D(j, h_lo)
    w_num = eval_p(wronskian(h_hi, th), ts)
    w_den = eval_p(wronskian(h_hi, Dh), ts)
    if np.min(w_den) <= 0:
        return None
    return w_num / w_den

def min_J(X, M, Y, ts):
    w1 = eval_p(wronskian(X, M), ts)
    w2 = eval_p(wronskian(M, Y), ts)
    w12 = eval_p(wronskian(X, Y), ts)
    if np.min(w1) <= 0 or np.min(w2) <= 0: return -1.0
    with np.errstate(divide='ignore', invalid='ignore'):
        J = 4.0 * w1 * w2 / (w12 * w12)
    J = np.where(np.abs(w12) < 1e-14 * np.sqrt(w1 * w2), 1e12, J)
    return float(np.min(J))

rng = np.random.default_rng(777)

print("=== Q1: HOW DOES STAGE 2 (D_j) ACT ON P = j*max(delta) AND L = j^2*max|delta_2 - delta_1|? ===")
# Apply D_j to (h0, h1, h2) -> (D_j h0, D_j h1, D_j h2) at degrees (3, 2, 1) and (4, 3, 2)!
for j in [15, 30, 60]:
    max_P_ratio = 0.0; max_P_diff = -1e9
    max_L_ratio = 0.0; max_L_diff = -1e9
    n_ok = 0
    for _ in range(3000):
        # Generate h2 << h1 << h0 in (-1, 0) of degs (3, 2, 1), sometimes close to GLZ-like spacing
        r0 = np.sort(-rng.uniform(0.01, 0.99, size=3))
        r1 = np.array([rng.uniform(r0[0], r0[1]), rng.uniform(r0[1], r0[2])])
        r2 = np.array([rng.uniform(r1[0], r1[1])])
        h0, h1, h2 = from_roots(r0), from_roots(r1), from_roots(r2)
        d1 = delta_vec(h0, h1, j, ts)
        d2 = delta_vec(h1, h2, j, ts)
        if d1 is None or d2 is None: continue
        Pin = j * max(float(np.max(d1)), float(np.max(d2)))
        if Pin > 4.0: continue
        Lin = (j**2) * float(np.max(np.abs(d2 - d1)))
        # Apply D_j to all three!
        Dh0, Dh1, Dh2 = apply_D(j, h0), apply_D(j, h1), apply_D(j, h2)
        d1o = delta_vec(Dh0, Dh1, j, ts)
        d2o = delta_vec(Dh1, Dh2, j, ts)
        if d1o is None or d2o is None: continue
        Pout = j * max(float(np.max(d1o)), float(np.max(d2o)))
        Lout = (j**2) * float(np.max(np.abs(d2o - d1o)))
        n_ok += 1
        max_P_ratio = max(max_P_ratio, Pout / Pin)
        max_P_diff = max(max_P_diff, Pout - Pin)
        if Lin > 0.1:
            max_L_ratio = max(max_L_ratio, Lout / Lin)
        max_L_diff = max(max_L_diff, Lout - Lin)
    print(f"  j={j:2d} (n={n_ok:4d}): P_out/P_in <= {max_P_ratio:.4f} (diff <= {max_P_diff:+.4f}) | L_out/L_in <= {max_L_ratio:.4f} (diff <= {max_L_diff:+.4f})")

print("\n=== Q3: STAGE 3 DROP j^2*(J(j) - J(j+1)) UNDER P <= 3.5 AND L <= C (C = 3, 5, 10)! ===")
for j in [20, 40, 80, 160]:
    # Generate triples (h0, h1, h2) of degs (3, 2, 1) with P <= 3.5, L <= C, and J(j) = 1 + 2/j!
    # To get MANY samples with small L = j^2|d2 - d1|, perturb r1, r2 around nearly uniform spacing!
    stats = {3: [0, -1e9], 5: [0, -1e9], 10: [0, -1e9], 1000: [0, -1e9]}
    for _ in range(15000):
        # Mix random + nearly evenly spaced roots to populate small L!
        if rng.random() < 0.6:
            # roots of deg 3, 2, 1 that are close to derivatives / GLZ-like
            a = rng.uniform(0.02, 0.4)
            b = rng.uniform(a + 0.05, 0.95)
            r0 = np.sort([-b, -0.5*(a+b) + rng.uniform(-0.1, 0.1)*(b-a), -a])
            r1 = np.array([0.5*(r0[0]+r0[1]) + rng.uniform(-0.15, 0.15)*(r0[1]-r0[0]),
                           0.5*(r0[1]+r0[2]) + rng.uniform(-0.15, 0.15)*(r0[2]-r0[1])])
            r2 = np.array([0.5*(r1[0]+r1[1]) + rng.uniform(-0.15, 0.15)*(r1[1]-r1[0])])
        else:
            r0 = np.sort(-rng.uniform(0.01, 0.99, size=3))
            r1 = np.array([rng.uniform(r0[0], r0[1]), rng.uniform(r0[1], r0[2])])
            r2 = np.array([rng.uniform(r1[0], r1[1])])
        h0, h1, h2 = from_roots(r0), from_roots(r1), from_roots(r2)
        d1 = delta_vec(h0, h1, j, ts)
        d2 = delta_vec(h1, h2, j, ts)
        if d1 is None or d2 is None: continue
        P = j * max(float(np.max(d1)), float(np.max(d2)))
        if P > 3.5: continue
        L = (j**2) * float(np.max(np.abs(d2 - d1)))
        # Tune h1 scale so J(j) == 1 + 2/j!
        Mj = apply_D(j, h1); Yj = apply_D(j, apply_D(j, h2))
        Ju = min_J(h0, Mj, Yj, ts)
        if Ju <= 0: continue
        c = np.sqrt((1.0 + 2.0/j) / Ju)
        Mjp = apply_D(j+1, h1) * c; Yjp = apply_D(j+1, apply_D(j+1, h2))
        Jout = min_J(h0, Mjp, Yjp, ts)
        drop_j2 = (j**2) * ((1.0 + 2.0/j) - Jout)
        for C in stats:
            if L <= C:
                stats[C][0] += 1
                stats[C][1] = max(stats[C][1], drop_j2)
    s_str = " | ".join(f"L<={C:4d}(n={stats[C][0]:4d}): j^2*d<={stats[C][1]:+7.2f}" for C in [3, 5, 10, 1000])
    print(f"  j={j:3d}: {s_str}")
