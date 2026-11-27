"""Targeted 3-part probe:
1. Does Q_{k, q}(t, y) = sum q^{i(i-1)/2} F_{k,i}(t) y^i at q = 1/8 stay real-stable on WILD m_j double-D products
   prod_{j=1}^k D_{m_j}(D_{m_j} + y)(1)?
2. Does the OUTER operator T_j in prod_{j=1}^k T_j(T_j + w_j y)(1) (general 0 < a_j < b_j, w_j > 0)
   STILL break Turan/real-rootedness like single-step prod (T_j + w_j y)(1), OR does the outer T_j prevent it?!
   (Use exact mpmath to avoid any float64 artifacts!)
3. Check the EXACT quantitative contraction of the Pascal step on J_ratio = 4 W01 W12 / W02^2:
   Given J_1 >= 1 + eps, J_2 >= 1 + eps, and a*b = R <= 1 - delta,
   what is the EXACT infimum of J_new = 4(1+a+u)(1+b+v) / ((1+u)(1+v) - a*b)^2
   over ALL u <= 2*sqrt(a/J_1), v <= 2*sqrt(b/J_2)?!
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 60

def apply_T_mp(a, b, poly):
    """T_{a,b} = t(1+t) d/dt + a + b*t (0 < a < b)."""
    res = [mp.mpf(0)]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r + a)
        res[r+1] += c*(r + b)
    return res

def poly_add_mp(p, q, s=mp.mpf(1)):
    res = [mp.mpf(0)]*max(len(p), len(q))
    for i, v in enumerate(p): res[i] += v
    for i, v in enumerate(q): res[i] += s * v
    return res

def eval_mp(p, tv):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * tv + c
    return acc

def deriv_mp(p): return [r*p[r] for r in range(1, len(p))]

print("=== 1. TEST Q_{k, q}(t, y) AT q = 1/8 ON WILD m_j DOUBLE-D PRODUCTS (EXACT MPMATH) ===")
def test_q_double_D(label, ms, q_val=mp.mpf("0.125")):
    polys = [[mp.mpf(1)]]
    for mv in ms:
        m = mp.mpf(mv)
        U = [[mp.mpf(0)] for _ in range(len(polys) + 1)]
        for i, p in enumerate(polys):
            U[i] = poly_add_mp(U[i], apply_T_mp(mp.mpf(1), m, p))
            U[i+1] = poly_add_mp(U[i+1], p)
        polys = [apply_T_mp(mp.mpf(1), m, u) for u in U]
    k = len(ms)
    # Check 3-term necessary condition: min_{t in R} W(p_{i-1},p_i)*W(p_i,p_{i+1}) / W(p_{i-1},p_{i+1})^2 vs q!
    # (Remember: for Q_{k,q} to be real-stable, q <= 4 * W01 * W12 / W02^2? Wait: in p_{i-1} + p_i y + q p_{i+1} y^2,
    # Rayleigh is W01 + 2 q W02 y + q W12 y^2 => disc = 4 q^2 W02^2 - 4 q W01 W12 <= 0 <=> q <= W01*W12 / W02^2!!)
    ps = [deriv_mp(p) for p in polys]
    ts = [mp.mpf(s) for s in ["-2", "-0.99", "-0.5", "-0.05", "-1e-5", "0", "1e-6", "1e-3", "0.05", "1"]]
    min_wrat = mp.mpf('1e100'); worst = None
    for tv in ts:
        vs = [eval_mp(p, tv) for p in polys]
        vps = [eval_mp(dp, tv) for dp in ps]
        w1 = [vps[i]*vs[i+1] - vs[i]*vps[i+1] for i in range(k)]
        w2 = [vps[i]*vs[i+2] - vs[i]*vps[i+2] for i in range(k-1)]
        for i in range(k-1):
            if w2[i] != 0:
                r = w1[i] * w1[i+1] / (w2[i] * w2[i])
                if r < min_wrat: min_wrat = r; worst = (i+1, float(tv))
    print(f"  {label:32s}: min W01*W12/W02^2 = {mp.nstr(min_wrat, 5)} at {worst}", flush=True)

test_q_double_D("GLZ k=6 (m=2..7)", [2,3,4,5,6,7])
test_q_double_D("step [2]*3+[10^4]*3", [2,2,2,10**4,10**4,10**4])
test_q_double_D("rev [10^4]*3+[2]*3", [10**4,10**4,10**4,2,2,2])
test_q_double_D("alt [2,10^4]*3", [2,10**4,2,10**4,2,10**4])
test_q_double_D("tiny [1.001]*3+[10^4]*3", ["1.001"]*3 + [10**4]*3)
test_q_double_D("rev tiny [10^4]*3+[1.001]*3", [10**4]*3 + ["1.001"]*3)

print("\n=== 2. TEST GENERAL DOUBLE-T: prod_{j=1}^k T_{a_j,b_j}(T_{a_j,b_j} + w_j y)(1) (CAN WE BREAK TURAN?) ===")
# Try to drive J_i(t) < 0 and then Turan < 0 using extreme (a_j, b_j, w_j) WITH the outer T_{a_j, b_j}!
# Note: in single-T, @algebraist broke Turan at k=7 by finding J_i(t) < 0 and then applying a T with Phi = a+bt hitting the negative root of J!
# Let's test 60 random extreme trials at k = 7, 8, 9 (fast float64 + exact mpmath confirm if negative!)
rng = np.random.default_rng(2026)
ts_neg_flt = np.concatenate([-10**np.linspace(-8, -1, 30), np.linspace(-0.12, -0.99, 40)])

def apply_T_flt(a, b, p):
    res = np.zeros(len(p) + 1)
    for r, c in enumerate(p):
        res[r] += c*(r + a); res[r+1] += c*(r + b)
    return res

def run_trial_double_T(params):
    polys = [[mp.mpf(1)]]
    for a, b, w in params:
        a, b, w = mp.mpf(a), mp.mpf(b), mp.mpf(w)
        U = [[mp.mpf(0)] for _ in range(len(polys) + 1)]
        for i, p in enumerate(polys):
            U[i] = poly_add_mp(U[i], apply_T_mp(a, b, p))
            U[i+1] = poly_add_mp(U[i+1], p, w)
        polys = [apply_T_mp(a, b, u) for u in U]
    return polys

# Also test targeted J_i < 0 attack: first build k=5 with extreme params, check if ANY J_i(t) < 0 on (-1, 0)!
neg_J_found = 0; neg_turan_found = 0
for trial in range(40):
    k = 6
    # Extreme bimodal (a, b, w)
    params = []
    for j in range(k):
        a = 10 ** rng.uniform(-3, 2)
        b = a + 10 ** rng.uniform(-3, 3)
        w = 10 ** rng.uniform(-2, 2)
        params.append((a, b, w))
    # Build in mpmath (k=6 is instant!)
    polys = run_trial_double_T(params)
    ps = [deriv_mp(p) for p in polys]
    ts_m = [mp.mpf(s) for s in ["-0.99", "-0.9", "-0.5", "-0.1", "-0.01", "-1e-4", "-1e-7"]]
    has_neg_J = False; has_neg_T = False
    for tv in ts_m:
        vs = [eval_mp(p, tv) for p in polys]
        vps = [eval_mp(dp, tv) for dp in ps]
        for i in range(1, k):
            if vs[i]*vs[i] - vs[i-1]*vs[i+1] < 0:
                has_neg_T = True
            w01 = vps[i-1]*vs[i] - vs[i-1]*vps[i]
            w12 = vps[i]*vs[i+1] - vs[i]*vps[i+1]
            w02 = vps[i-1]*vs[i+1] - vs[i-1]*vps[i+1]
            if 4*w01*w12 - w02*w02 < 0:
                has_neg_J = True
    if has_neg_J: neg_J_found += 1
    if has_neg_T: neg_turan_found += 1
print(f"  General double-T prod T_j(T_j + w_j y)(1) at k=6 (40 extreme trials): neg_J={neg_J_found}/40, neg_Turan={neg_turan_found}/40", flush=True)

print("\n=== 3. EXACT WORST-CASE PASCAL STEP ON J_ratio = 4(1+a+u)(1+b+v) / ((1+u)(1+v) - ab)^2 ===")
# Given J1 = J2 = 1 + eps, and a*b = R < 1, minimize J_new over u in [0, 2*sqrt(a/J1)], v in [0, 2*sqrt(b/J2)]
# and a, b > 0 with a*b = R!
# Check specifically at GLZ scale: a = b = sqrt(R) vs asymmetric a != b!
for k in [10, 20, 50, 100]:
    eps = 4.0 / k
    R = 1.0 - 4.0 / k  # actually in GLZ R ~ 1 - 6/k!
    J_in = 1.0 + eps
    min_J_sym = 1e9; min_J_all = 1e9; best_cfg = None
    for log_ratio in np.linspace(-2, 2, 41):
        a = math.sqrt(R) * math.exp(log_ratio)
        b = R / a
        umax = 2.0 * math.sqrt(a / J_in)
        vmax = 2.0 * math.sqrt(b / J_in)
        # Check boundary u = umax, v = vmax (and interior grid)
        for fu in [0.0, 0.5, 0.9, 1.0]:
            for fv in [0.0, 0.5, 0.9, 1.0]:
                u = fu * umax; v = fv * vmax
                den = (1 + u)*(1 + v) - R
                jn = 4*(1 + a + u)*(1 + b + v) / (den * den)
                if jn < min_J_all:
                    min_J_all = jn; best_cfg = (a, b, fu, fv)
                if abs(log_ratio) < 1e-9 and jn < min_J_sym:
                    min_J_sym = jn
    print(f"  k={k:3d} (J_in={J_in:.5f}, R={R:.4f}): min J_new (sym a=b) = {min_J_sym:.5f}, min J_new (all a*b=R) = {min_J_all:.5f} at fu={best_cfg[2]},fv={best_cfg[3]}")
