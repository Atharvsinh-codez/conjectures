import numpy as np

# TEST THE ACTUAL FULL COMBINED GLZ 1-STEP MAP G_j:
#   Input: 4-tuple (h0, h1, h2, h3) of degrees (d+2, d+1, d, d-1) with roots in (-1, 0), h3 << h2 << h1 << h0,
#          such that BOTH base triples (h0, D_j h1, D_j^2 h2) and (h1, D_j h2, D_j^2 h3) have J >= 1 + 2/(j-1)
#          and W(h0, D_j^3 h3) >= 0 (open half-plane).
#   Output: ( D_j(h0 + D_j h1),  D_{j+1} D_j(h1 + D_j h2),  D_{j+1}^2 D_j(h2 + D_j h3) )
#   Question: Does J_out >= 1 + 2/j HOLD (both with and without pole-gap r* <= 0.75 j)?!
#   Notice: this tests Stage 1 + Stage 2 + Stage 3 COMBINED as a SINGLE step!

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
    np.linspace(-0.999, -0.001, 120),
    np.geomspace(1e-4, 1e2, 40),
])

def min_J(X, M, Y, ts):
    w1 = eval_p(wronskian(X, M), ts)
    w2 = eval_p(wronskian(M, Y), ts)
    w12 = eval_p(wronskian(X, Y), ts)
    if np.min(w1) <= 0 or np.min(w2) <= 0:
        return -1.0
    with np.errstate(divide='ignore', invalid='ignore'):
        J = 4.0 * w1 * w2 / (w12 * w12)
    J = np.where(np.abs(w12) < 1e-14 * np.sqrt(w1 * w2), 1e12, J)
    return float(np.min(J))

def pole_rmax(h_hi, h_lo, ts):
    D0 = apply_D(0.0, h_lo)
    th = np.concatenate([[0.0], h_lo])
    w_num = eval_p(wronskian(h_hi, D0), ts)
    w_den = eval_p(wronskian(h_hi, th), ts)
    return float(np.max(-w_num / w_den))

print("=== FULL COMBINED GLZ STEP G_j AT DEGREES (3, 2, 1, 0) ===")
rng = np.random.default_rng(2026)
for j in [10, 20, 40, 80, 160]:
    eps_in = 2.0 / (j - 1)
    eps_out = 2.0 / j
    n_all = 0; fail_all = 0; min_gap_all = 1e9; worst_all = None
    n_pg = 0; fail_pg = 0; min_gap_pg = 1e9; worst_pg = None
    for _ in range(5000):
        # Generate interlacing h3 << h2 << h1 << h0 in (-1, 0) of degs (3, 2, 1, 0)
        if rng.random() < 0.4:
            # adversarial: put roots near -1/j!
            r0 = np.sort(-np.exp(rng.uniform(-np.log(5*j), -0.02, size=3)))
        else:
            r0 = np.sort(-rng.uniform(0.005, 0.995, size=3))
        r1 = np.array([rng.uniform(r0[0], r0[1]), rng.uniform(r0[1], r0[2])])
        r2 = np.array([rng.uniform(r1[0], r1[1])])
        h0 = from_roots(r0, 1.0)
        h1_u = from_roots(r1, 1.0)
        h2_u = from_roots(r2, 1.0)
        h3_u = np.array([1.0])
        # Compute J of unit-scaled triples (h0, D_j h1, D_j^2 h2) and (h1, D_j h2, D_j^2 h3)
        # Note J(s0 h0, s1 D_j h1, s2 D_j^2 h2) = (s1^2 / (s0 s2)) * J_unit!
        # So we can EXACTLY tune s1, s2, s3 so BOTH J1 and J2 are >= 1 + eps_in (even right AT 1 + eps_in)!
        M1_u = apply_D(j, h1_u); Y1_u = apply_D(j, apply_D(j, h2_u))
        M2_u = apply_D(j, h2_u); Y2_u = apply_D(j, apply_D(j, h3_u))
        J1_u = min_J(h0, M1_u, Y1_u, ts)
        J2_u = min_J(h1_u, M2_u, Y2_u, ts)
        if J1_u <= 0 or J2_u <= 0: continue
        # Pick s1 randomly, then s2 = s1^2 * J1_u / target1, s3 = s2^2 * J2_u / (s1 * target2)
        s1 = np.exp(rng.uniform(-2, 2))
        t1 = (1.0 + eps_in) * (1.0 + rng.uniform(0.0, 0.2/j))
        t2 = (1.0 + eps_in) * (1.0 + rng.uniform(0.0, 0.2/j))
        s2 = s1*s1 * J1_u / t1
        s3 = s2*s2 * J2_u / (s1 * t2)
        h1, h2, h3 = h1_u * s1, h2_u * s2, h3_u * s3
        # Check open half-plane W(h0, D_j^3 h3) >= 0!
        E = apply_D(j, apply_D(j, apply_D(j, h3)))
        if np.min(eval_p(wronskian(h0, E), ts)) < 0: continue

        Jin = min(t1, t2)
        # Apply FULL GLZ step G_j!
        out0 = apply_D(j, add_p(h0, apply_D(j, h1)))
        out1 = apply_D(j+1, apply_D(j, add_p(h1, apply_D(j, h2))))
        out2 = apply_D(j+1, apply_D(j+1, apply_D(j, add_p(h2, apply_D(j, h3)))))
        Jout = min_J(out0, out1, out2, ts)
        margin = Jout - (1.0 + eps_out)

        rstar = max(pole_rmax(h0, h1, ts), pole_rmax(h1, h2, ts), pole_rmax(h2, h3, ts))
        n_all += 1
        if margin < 0: fail_all += 1
        if margin < min_gap_all: min_gap_all = margin; worst_all = (rstar, Jin, Jout)
        if rstar <= 0.75 * j:
            n_pg += 1
            if margin < 0: fail_pg += 1
            if margin < min_gap_pg: min_gap_pg = margin; worst_pg = (rstar, Jin, Jout)
    print(f"  j={j:3d}: [ALL n={n_all:4d}] fails={fail_all:4d}, min(Jout-(1+2/j))={min_gap_all:+.6f} (j^2*m={j*j*min_gap_all:+7.2f}) | [r*<=.75j n={n_pg:4d}] fails={fail_pg:4d}, min={min_gap_pg:+.6f} (j^2*m={j*j*min_gap_pg:+7.2f})")
