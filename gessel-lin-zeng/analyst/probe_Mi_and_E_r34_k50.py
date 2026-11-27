import numpy as np
import mpmath as mp

# PART 1: Brutal adversarial test of whether D_b REALLY preserves M_i(t) > 0 universally!
# Recall M(p0, p1, p2) = (-t)(1+t) W(p2, p1) + p1^2 - p0*p2, where p2 << p1 << p0 in (-1, 0).
# Push p0, p1, p2 to have VERY SMALL J_raw (e.g. J_raw < 0.1!) and min M_i -> 1e-4, then apply D_b!

def apply_D_np(b, p):
    res = np.zeros(len(p)+1)
    for n, c in enumerate(p):
        res[n] += c*(n+1); res[n+1] += c*(n+b)
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

def from_roots(rs, scale=1.0):
    p = np.array([scale], dtype=np.float64)
    for r in rs: p = np.convolve(p, [-r, 1.0])
    return p

ts = np.concatenate([
    -1.0 + np.geomspace(1e-6, 0.49, 120),
    -np.geomspace(1e-6, 0.5, 120),
])

def eval_M(p0, p1, p2, ts):
    w21 = eval_np(wronskian_np(p2, p1), ts)
    v0, v1, v2 = eval_np(p0, ts), eval_np(p1, ts), eval_np(p2, ts)
    return (-ts)*(1.0 + ts)*w21 + v1*v1 - v0*v2

print("=== 1. BRUTAL ADVERSARIAL TEST OF D_b PRESERVING M_i(t) > 0 ON (-1, 0) ===")
# Make p1^2 - p0*p2 VERY NEGATIVE while keeping M_i > 0 ONLY barely (by tuning s2 = scale of p2)!
rng = np.random.default_rng(999)
for b in [1.01, 2.0, 10.0, 1000.0, 1e6]:
    valid = 0; fail = 0; worst_rel = 1e9; worst_cfg = None
    for _ in range(4000):
        # Use extreme roots (including clustered roots near 0 or -1!)
        if rng.random() < 0.5:
            r0 = np.sort(-10**rng.uniform(-4, 0, size=3))
        else:
            r0 = np.sort(-1.0 + 10**rng.uniform(-4, 0, size=3))
        r1 = np.array([rng.uniform(r0[0], r0[1]), rng.uniform(r0[1], r0[2])])
        r2 = np.array([rng.uniform(r1[0], r1[1])])
        p0, p1, p2_unit = from_roots(r0, 1.0), from_roots(r1, 1.0), from_roots(r2, 1.0)
        # Find exact max s2 such that M(p0, p1, s2*p2_unit) >= 0 on ts!
        # M = s2 * [(-t)(1+t)W(p2_u, p1) - p0*p2_u] + p1^2
        # Wait: at points where coeff_s2 < 0, we need s2 <= p1^2 / (-coeff_s2)!
        w21_u = eval_np(wronskian_np(p2_unit, p1), ts)
        v0, v1, v2_u = eval_np(p0, ts), eval_np(p1, ts), eval_np(p2_unit, ts)
        coeff_s2 = (-ts)*(1.0 + ts)*w21_u - v0*v2_u
        neg_mask = coeff_s2 < 0
        if not np.any(neg_mask): continue
        s2_max = float(np.min((v1[neg_mask]**2) / (-coeff_s2[neg_mask])))
        if s2_max <= 0 or not np.isfinite(s2_max): continue
        # Set s2 = 0.999 * s2_max (RAZOR'S EDGE of M_i > 0!)
        p2 = p2_unit * (0.999 * s2_max)
        valid += 1
        q0, q1, q2 = apply_D_np(b, p0), apply_D_np(b, p1), apply_D_np(b, p2)
        Mq = eval_M(q0, q1, q2, ts)
        scale_q = np.abs(eval_np(q1, ts))**2 + np.abs(eval_np(q0, ts)*eval_np(q2, ts)) + 1e-300
        rel = float(np.min(Mq / scale_q))
        if rel < -1e-7:
            fail += 1
        if rel < worst_rel:
            worst_rel = rel; worst_cfg = (r0, r1, r2, s2_max)
    print(f"  b={b:9.2f}: razor-edge valid={valid:4d}, fail={fail:4d}, min_rel_Mq={worst_rel:+.6f}")

print("\n=== 2. CHECK W(E_{k,i}, E_{k,i+r}) > 0 FOR r=3..8 ON GLZ UP TO k = 36 (AT t NEAR -1)! ===")
# Use mpmath dps=80 to avoid any precision issues!
mp.mp.dps = 80
def apply_D_int(m, p):
    res = [0]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv_int(p): return [r*p[r] for r in range(1, len(p))] if len(p) > 1 else [0]
def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v
def eval_W_mp(f, g, t):
    return eval_mp(f, t)*eval_mp(deriv_int(g), t) - eval_mp(deriv_int(f), t)*eval_mp(g, t)

KMAX = 34
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    H[k][0] = apply_D_int(m, apply_D_int(m, H[k-1][0]))
    for i in range(1, k):
        H[k][i] = apply_D_int(m, poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1]))
    H[k][k] = apply_D_int(m, H[k-1][k-1])

ts_near_m1 = [mp.mpf(s) for s in ["-1.01", "-1", "-0.999", "-0.99", "-0.97", "-0.95", "-0.9", "-0.8", "-0.5", "-0.1", "0", "1"]]
for k in [20, 26, 30, 34]:
    m = k + 1
    # Check i = 0..8 against i+r for r = 3, 4, 5, 6, 7, 8! (and i = k//2 against i+r!)
    E = {}
    check_is = sorted(set(list(range(0, min(12, k))) + list(range(k//2 - 2, min(k, k//2 + 9)))))
    for i in check_is:
        p = H[k-1][i]
        for _ in range(i): p = apply_D_int(m, p)
        E[i] = p
    out = []
    for r in [3, 4, 5, 6, 7, 8]:
        fails = 0; first = None
        for i in check_is:
            if (i + r) not in E: continue
            for tv in ts_near_m1:
                if eval_W_mp(E[i], E[i+r], tv) <= 0:
                    fails += 1
                    if first is None: first = (i, float(tv))
                    break
        out.append(f"r={r}:{fails}" + (f"@{first}" if first else ""))
    # Also compute max cross-ratio chi = W01*W23/(W02*W13) for i = 0..4!
    max_chi = mp.mpf(0)
    for i in [0, 1, 2, k//2]:
        if all((i+d) in E for d in [0, 1, 2, 3]):
            for tv in ts_near_m1:
                w01 = eval_W_mp(E[i], E[i+1], tv)
                w23 = eval_W_mp(E[i+2], E[i+3], tv)
                w02 = eval_W_mp(E[i], E[i+2], tv)
                w13 = eval_W_mp(E[i+1], E[i+3], tv)
                if w02 > 0 and w13 > 0:
                    chi = (w01*w23)/(w02*w13)
                    if chi > max_chi: max_chi = chi
    print(f"  k={k:2d}: " + " | ".join(out) + f" || max_chi={float(max_chi):.5f}")
