import numpy as np
import mpmath as mp

# TEST 1: On GLZ (k=3..24) AND extreme b_j >= 2 sequences (single & doubled),
#         what is min_{i, t in R} J_raw(h_{k, i-1}, h_{k, i}, h_{k, i+1})?!
#         Recall J_raw(f0, f1, f2) = 4 W(f0, f1) W(f1, f2) / W(f0, f2)^2,
#         and J_raw >= q <==> q-deformed triple (f0, f1, q f2) is real-stable!
#
# TEST 2: 1-STEP INVARIANCE TEST (ZERO STAGE 3!!):
#         Given arbitrary (h0, h1, h2, h3) of degrees (d+2, d+1, d, d-1) with roots in (-1, 0),
#         h3 << h2 << h1 << h0, and J_raw(h0, h1, h2) >= gamma, J_raw(h1, h2, h3) >= gamma:
#         Does the 1-step map:
#           (Sgl) (D_b h1 + h0, D_b h2 + h1, D_b h3 + h2)
#           (Dbl) D_b of the above
#         PRESERVE J_raw >= gamma for all b >= 2?! (Test gamma = 0.888, 0.5, 0.25, 0.125!)

mp.mp.dps = 60

def apply_D_mp(b, p):
    res = [mp.mpf(0)]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+b)
    return res

def add_mp(a, b):
    res = [mp.mpf(0)]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv_mp(p):
    return [r*p[r] for r in range(1, len(p))] if len(p) > 1 else [mp.mpf(0)]

def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v

def eval_W_mp(f, g, t):
    return eval_mp(f, t)*eval_mp(deriv_mp(g), t) - eval_mp(deriv_mp(f), t)*eval_mp(g, t)

ts_mp = [mp.mpf(s) for s in [
    "-1000", "-10", "-2", "-1.1", "-1.001", "-1", "-0.999", "-0.99", "-0.95", "-0.9",
    "-0.7", "-0.5", "-0.3", "-0.1", "-0.01", "-1e-4", "-1e-8", "0", "1e-6", "0.01", "0.1", "1", "10", "1000"
]]

def min_Jraw_seq(P, ts):
    mn = mp.mpf('1e100')
    loc = None
    for i in range(1, len(P) - 1):
        for tv in ts:
            w1 = eval_W_mp(P[i-1], P[i], tv)
            w2 = eval_W_mp(P[i], P[i+1], tv)
            w12 = eval_W_mp(P[i-1], P[i+1], tv)
            if w1 <= 0 or w2 <= 0:
                return -1.0, (i, float(tv))
            if w12 != 0:
                J = 4*w1*w2/(w12*w12)
                if J < mn:
                    mn = J; loc = (i, float(tv))
    return float(mn), loc

print("=== 1. MIN J_raw ON ALL R FOR GLZ AND EXTREME b_j >= 2 SEQUENCES (SINGLE & DOUBLED) ===")
# Note: in our standard convention, tilde_h_{1,0} = 1+2t, tilde_h_{1,1} = 1;
# for general b_j starting from P = [[1]], step j applies D_{b_j}(D_{b_j} + y) or (D_{b_j} + y)!
# Note starting from [1] with b_1 = 2 gives [1+2t, 1] (single) or [D_2(1+2t), 1+2t] (double).
test_seqs = [
    ("GLZ (k=12, dbl, b=3..13)", [mp.mpf(j+2) for j in range(1, 12)], True, [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]),
    ("GLZ (k=20, dbl, b=3..21)", [mp.mpf(j+2) for j in range(1, 20)], True, [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]),
    ("b=[2]*12 (sgl)",           [mp.mpf(2)]*12, False, [[mp.mpf(1)]]),
    ("b=[2]*12 (dbl)",           [mp.mpf(2)]*12, True,  [[mp.mpf(1)]]),
    ("b=[2]*5+[1e6]*5 (sgl)",    [mp.mpf(2)]*5 + [mp.mpf('1e6')]*5, False, [[mp.mpf(1)]]),
    ("b=[2]*5+[1e6]*5 (dbl)",    [mp.mpf(2)]*5 + [mp.mpf('1e6')]*5, True,  [[mp.mpf(1)]]),
    ("b=[1e6]*5+[2]*5 (sgl)",    [mp.mpf('1e6')]*5 + [mp.mpf(2)]*5, False, [[mp.mpf(1)]]),
    ("b=[1e6]*5+[2]*5 (dbl)",    [mp.mpf('1e6')]*5 + [mp.mpf(2)]*5, True,  [[mp.mpf(1)]]),
    ("b=[2,1e6]*5 (sgl)",        [mp.mpf(2), mp.mpf('1e6')]*5, False, [[mp.mpf(1)]]),
    ("b=[2,1e6]*5 (dbl)",        [mp.mpf(2), mp.mpf('1e6')]*5, True,  [[mp.mpf(1)]]),
    ("b=[1.01]*5+[1e6]*5 (sgl)", [mp.mpf('1.01')]*5 + [mp.mpf('1e6')]*5, False, [[mp.mpf(1)]]),
    ("b=[1.5]*5+[1e6]*5 (sgl)",  [mp.mpf('1.5')]*5 + [mp.mpf('1e6')]*5, False, [[mp.mpf(1)]]),
]

for name, bs, dbl, initP in test_seqs:
    P = initP
    for b in bs:
        k = len(P)
        newP = []
        for i in range(k + 1):
            t2 = apply_D_mp(b, P[i]) if i < k else [mp.mpf(0)]
            t1 = P[i-1] if i >= 1 else [mp.mpf(0)]
            u = add_mp(t2, t1)
            newP.append(apply_D_mp(b, u) if dbl else u)
        P = newP
    mn, loc = min_Jraw_seq(P, ts_mp)
    print(f"  {name:28s}: min J_raw = {mn:+10.6f}  at (i, t)={loc}")

print("\n=== 2. 1-STEP INVARIANCE TEST OF J_raw >= gamma UNDER (D_b h1+h0, D_b h2+h1, D_b h3+h2) ===")
# Fast numpy test with degrees (3, 2, 1, 0) and (4, 3, 2, 1)!
def apply_D_np(b, p):
    res = np.zeros(len(p)+1)
    for n, c in enumerate(p):
        res[n] += c*(n+1); res[n+1] += c*(n+b)
    return res

def add_np(a, b):
    n = max(len(a), len(b))
    r = np.zeros(n); r[:len(a)] += a; r[:len(b)] += b
    return r

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

def min_J_np(f0, f1, f2, ts):
    w1 = eval_np(wronskian_np(f0, f1), ts)
    w2 = eval_np(wronskian_np(f1, f2), ts)
    w12 = eval_np(wronskian_np(f0, f2), ts)
    if np.min(w1) <= 0 or np.min(w2) <= 0:
        return -1.0
    with np.errstate(divide='ignore', invalid='ignore'):
        J = 4.0 * w1 * w2 / (w12 * w12)
    J = np.where(np.abs(w12) < 1e-14 * np.sqrt(w1 * w2), 1e12, J)
    return float(np.min(J))

def from_roots(rs, scale=1.0):
    p = np.array([scale], dtype=np.float64)
    for r in rs: p = np.convolve(p, [-r, 1.0])
    return p

ts_np = np.concatenate([
    -np.geomspace(1e-4, 1e3, 60),
    np.linspace(-0.9999, -0.0001, 200),
    -1.0 - np.geomspace(1e-4, 1e3, 60),
    np.geomspace(1e-4, 1e3, 60),
])

rng = np.random.default_rng(12345)
for degs in [(3, 2, 1, 0), (4, 3, 2, 1)]:
    print(f"\n  --- Degrees {degs} ---")
    for b in [2.0, 5.0, 50.0, 1000.0]:
        for gamma in [0.88, 0.5, 0.25, 0.1]:
            valid = 0; fail_sgl = 0; fail_dbl = 0
            min_sgl = 1e9; min_dbl = 1e9
            for _ in range(3000):
                # Generate interlacing h3 << h2 << h1 << h0 in (-1, 0)
                r0 = np.sort(-rng.uniform(0.01, 0.99, size=degs[0]))
                r1 = np.array([rng.uniform(r0[m], r0[m+1]) for m in range(degs[1])])
                r2 = np.array([rng.uniform(r1[m], r1[m+1]) for m in range(degs[2])])
                r3 = np.array([rng.uniform(r2[m], r2[m+1]) for m in range(degs[3])]) if degs[3] > 0 else np.array([])
                # Random positive scales (try both balanced AND skewed!)
                s0 = 1.0
                s1 = np.exp(rng.uniform(-2, 2))
                s2 = np.exp(rng.uniform(-2, 2))
                s3 = np.exp(rng.uniform(-2, 2))
                h0, h1, h2, h3 = from_roots(r0, s0), from_roots(r1, s1), from_roots(r2, s2), from_roots(r3, s3)
                # Ensure W(h0, h1) > 0? Wait: deg h0 = d+1, deg h1 = d, both positive leading coeffs, r0 surrounds r1:
                # At t -> +inf: h0 ~ t^{d+1}, h1 ~ t^d => W(h0, h1) = h0 h1' - h0' h1 ~ (d - (d+1)) t^{2d} = -t^{2d} < 0!
                # Check sign! In GLZ, h_{1,0} = 1+2t, h_{1,1} = 1 -> W(1+2t, 1) = -2 < 0!
                # So W(h0, h1) < 0 and W(h1, h2) < 0, hence W01 * W12 > 0!
                # Flip sign of W by passing (h2, h1, h0) or checking w1 < 0 and w2 < 0!
                J1 = min_J_np(h2, h1, h0, ts_np)
                if J1 < gamma: continue
                J2 = min_J_np(h3, h2, h1, ts_np)
                if J2 < gamma: continue
                # Also check if W(h3, h0) > 0 (half-plane condition!)
                valid += 1
                g0 = add_np(apply_D_np(b, h1), h0)
                g1 = add_np(apply_D_np(b, h2), h1)
                g2 = add_np(apply_D_np(b, h3), h2)
                Js = min_J_np(g2, g1, g0, ts_np)
                if Js < gamma - 1e-6: fail_sgl += 1
                if Js < min_sgl: min_sgl = Js
                d0, d1, d2 = apply_D_np(b, g0), apply_D_np(b, g1), apply_D_np(b, g2)
                Jd = min_J_np(d2, d1, d0, ts_np)
                if Jd < gamma - 1e-6: fail_dbl += 1
                if Jd < min_dbl: min_dbl = Jd
            print(f"    b={b:6.1f}, gamma={gamma:4.2f}: valid={valid:4d} | sgl_fail={fail_sgl:4d} (min={min_sgl:+7.4f}) | dbl_fail={fail_dbl:4d} (min={min_dbl:+7.4f})")
