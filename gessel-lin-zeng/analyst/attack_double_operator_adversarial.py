import numpy as np
import mpmath as mp

# Adversarial attack on DOUBLE-OPERATOR family:
#   P_k(t, y) = prod_{j=1}^k T_{j,1} (T_{j,2} + w_j y) (1)
# Recall how @algebraist broke SINGLE-operator prod_{j=1}^k (T_j + w_j y)(1) at k=9:
#   Wait: why did single-operator break? Because T_j p_i + w_j p_{i-1} can break the INTERLACING
#   between u_i and u_{i-1}, i.e., W(u_i, u_{i-1}) can change sign, OR Turan p_i^2 - p_{i-1}p_{i+1} < 0!
#   Wait: in Theorem A (verified in #23!), single step (T + wy) PRESERVES p_i << p_{i-1} WHENEVER
#   Turan holds: p_i^2 >= p_{i-1} p_{i+1} on (-1, 0)!!
#   So the ONLY way single-operator failed at k=9 was by FIRST violating Turan p_i^2 < p_{i-1} p_{i+1}
#   at some step k0 <= 8, and THEN applying a huge w_9 (Wait: huge w_9 pushes u_i = T p_i + w p_{i-1} -> w p_{i-1}
#   and u_{i-1} -> w p_{i-2}, wait: when does u_i << u_{i-1} fail? When T p_i and w p_{i-1} are COMPARABLE
#   at the exact t where p_i^2 < p_{i-1} p_{i+1}!).
#
# LET'S TEST THIS EXACT MECHANISM ON DOUBLE-OPERATOR T_{j,1}(T_{j,2} + w_j y)(1)!!
# Step 1: Can we make Turan p_{k,i}^2 - p_{k,i-1} p_{k,i+1} < 0 on (-1, 0) for DOUBLE-OPERATOR at small k?!
#   (We ALREADY know GLZ itself violates raw Turan at k=35! Can extreme a_j, b_j, w_j violate raw Turan at k=4, 5, 6, 7?!)
# Step 2: Once we have a point t0 in (-1, 0) where p_i(t0)^2 < p_{i-1}(t0) p_{i+1}(t0),
#   can we choose the NEXT step T_{k+1}(T_{k+1} + w_{k+1} y) so that q_{i+1} = T(T p_{i+1} + w p_i)
#   and q_i = T(T p_i + w p_{i-1}) VIOLATE INTERLACING, and then the step AFTER breaks REAL-ROOTEDNESS?!

def apply_T_flt(p, a, b):
    d = len(p) - 1
    res = np.zeros(d + 2, dtype=np.float64)
    for n in range(d + 1):
        res[n] += (a + n) * p[n]
        res[n + 1] += (b + n) * p[n]
    return res

def add_p(a, b):
    n = max(len(a), len(b))
    r = np.zeros(n)
    r[:len(a)] += a
    r[:len(b)] += b
    return r

def eval_p(p, ts):
    v = np.zeros_like(ts)
    for c in p[::-1]:
        v = v * ts + c
    return v

def step_double(P, a1, b1, a2, b2, w):
    k = len(P)
    newP = []
    for i in range(k + 1):
        t2 = apply_T_flt(P[i], a2, b2) if i < k else np.array([0.0])
        wy = w * P[i - 1] if i >= 1 else np.array([0.0])
        newP.append(apply_T_flt(add_p(t2, wy), a1, b1))
    return newP

def step_single(P, a, b, w):
    k = len(P)
    newP = []
    for i in range(k + 1):
        t2 = apply_T_flt(P[i], a, b) if i < k else np.array([0.0])
        wy = w * P[i - 1] if i >= 1 else np.array([0.0])
        newP.append(add_p(t2, wy))
    return newP

ts = np.concatenate([
    -10**np.linspace(-6, -0.001, 200),
    -1.0 + 10**np.linspace(-6, -0.001, 200),
    np.linspace(-0.99, -0.01, 200),
])

print("=== 1. CAN DOUBLE-OPERATOR T(T + wy)(1) VIOLATE RAW TURAN p_i^2 >= p_{i-1} p_{i+1} ON (-1, 0) AT SMALL k? ===")
# Recall how GLZ violates raw Turan at k=35: near t = 0^-!
# At t -> 0^-, check the Taylor expansion of p_i(t) = c_{i,0} + c_{i,1} t + ...
# WAIT A SECOND! At t = 0: T_{a,b} f(0) = a * f(0)!
# So at t = 0, P_k(0, y) = prod_{j=1}^k a_{j,1} (a_{j,2} + w_j y), which IS real-rooted in y,
# so at t = 0, p_i(0)^2 > p_{i-1}(0) p_{i+1}(0) ALWAYS (by Newton's inequality!).
# Where did GLZ violate Turan? Near t = -1/k (around the RIGHTMOST root of p_{k,i}!)!
# Let's check min_{t in (-1, 0)} (p_i(t)^2 - p_{i-1}(t) p_{i+1}(t)) / (p_i^2 + |p_{i-1}p_{i+1}|)
# for random extreme (a_j, b_j, w_j) at k = 3, 4, 5, 6, 7, 8!

rng = np.random.default_rng(2026)
for mode in ["single_T", "double_same_T", "double_Dm_w1"]:
    best_turan = (1e9, None)
    neg_count = 0
    for trial in range(400):
        P = [np.array([1.0])]
        params = []
        for j in range(7):
            if mode == "double_Dm_w1":
                a = 1.0
                b = 1.0 + 10**rng.uniform(-3, 4)
                w = 1.0
            else:
                a = 10**rng.uniform(-3, 2)
                b = a + 10**rng.uniform(-3, 4)
                w = 10**rng.uniform(-3, 4)
            params.append((a, b, w))
            if mode == "single_T":
                P = step_single(P, a, b, w)
            else:
                P = step_double(P, a, b, a, b, w)
            # normalize P to prevent overflow
            mx = max(np.max(np.abs(p)) for p in P)
            P = [p / mx for p in P]
        # Check Turan at k=7
        found_neg = False
        for i in range(1, len(P) - 1):
            v0 = eval_p(P[i-1], ts)
            v1 = eval_p(P[i], ts)
            v2 = eval_p(P[i+1], ts)
            denom = v1*v1 + np.abs(v0*v2) + 1e-300
            rel = (v1*v1 - v0*v2) / denom
            mn = float(np.min(rel))
            if mn < best_turan[0]:
                best_turan = (mn, (i, float(ts[np.argmin(rel)]), params))
            if mn < -1e-6:
                found_neg = True
        if found_neg:
            neg_count += 1
    print(f"  {mode:15s} (k=7, 400 trials): neg_turan_trials={neg_count:3d}/400, min_rel_turan={best_turan[0]:+.6f} at i={best_turan[1][0]}, t={best_turan[1][1]:.5e}")

print("\n=== 2. TRY TO BREAK REAL-ROOTEDNESS OF DOUBLE-OPERATOR BY TARGETING NEG-TURAN POINT! ===")
# Take the worst Turan configuration for double_same_T (if negative!), then hit with 1 or 2 adversarial steps!
# Specifically, if p_i(t0)^2 < p_{i-1}(t0) p_{i+1}(t0) at t0 in (-1, 0),
# check if there is (a, b, w) such that u_i = T_{a,b} p_i + w p_{i-1} and u_{i+1} = T_{a,b} p_{i+1} + w p_i
# FAIL TO INTERLACE (which implies after T_{a,b} they STILL fail to interlace, and then next step can break RZ!!)!
# Wait: recall Theorem A! u_{i+1} << u_i fails at a root r of p_i ONLY when:
#   u_{i+1}(r) u_i(r) = (T p_{i+1}(r) + w p_i(r)) * (T p_i(r) + w p_{i-1}(r))
# Wait: at a root r of p_i, p_i(r) = 0, so p_i(r)^2 - p_{i-1}(r)p_{i+1}(r) = -p_{i-1}(r)p_{i+1}(r) > 0 ALWAYS (by interlacing)!
# Where CAN Turan fail? At a point t0 where p_i(t0) != 0 and p_{i-1}(t0), p_{i+1}(t0) have the SAME sign!
# Then choosing w so that u_i(t0) = T p_i(t0) + w p_{i-1}(t0) = 0 (i.e. w = -T p_i(t0) / p_{i-1}(t0) > 0)
# makes u_{i+1}(t0) = T p_{i+1}(t0) + w p_i(t0) = [p_{i-1} T p_{i+1} - p_i T p_i] / p_{i-1}!
# Let's search directly for ANY non-interlacing pair (u_{i+1}, u_i) in double_same_T and double_Dm_w1!!

def check_interlace_T(f_lo, f_hi):
    # deg f_lo = deg f_hi - 1; check all real roots in (-1, 0) and interlacing
    r_lo = np.roots(f_lo[::-1])
    r_hi = np.roots(f_hi[::-1])
    if np.max(np.abs(r_lo.imag)) > 1e-5 or np.max(np.abs(r_hi.imag)) > 1e-5:
        return "COMPLEX", max(float(np.max(np.abs(r_lo.imag))), float(np.max(np.abs(r_hi.imag))))
    r_lo = np.sort(r_lo.real); r_hi = np.sort(r_hi.real)
    ok = all(r_hi[m] < r_lo[m] < r_hi[m+1] for m in range(len(r_lo)))
    return ("OK" if ok else "NO_INTERLACE"), 0.0

for mode in ["double_same_T", "double_Dm_w1"]:
    # First find a k=6 or k=7 sequence with the MOST negative Turan
    best = (1e9, None, None)
    for trial in range(500):
        P = [np.array([1.0])]
        params = []
        for j in range(6):
            if mode == "double_Dm_w1":
                a = 1.0; b = 1.0 + 10**rng.uniform(-3, 5); w = 1.0
            else:
                a = 10**rng.uniform(-4, 2); b = a + 10**rng.uniform(-3, 5); w = 10**rng.uniform(-4, 5)
            params.append((a, b, w))
            P = step_double(P, a, b, a, b, w)
            mx = max(np.max(np.abs(p)) for p in P)
            P = [p / mx for p in P]
        for i in range(1, len(P) - 1):
            v0, v1, v2 = eval_p(P[i-1], ts), eval_p(P[i], ts), eval_p(P[i+1], ts)
            rel = (v1*v1 - v0*v2) / (v1*v1 + np.abs(v0*v2) + 1e-300)
            mn = float(np.min(rel))
            if mn < best[0]:
                best = (mn, P, (i, float(ts[np.argmin(rel)]), params))
    print(f"\n  [{mode}] worst k=6 Turan = {best[0]:+.6f} at i={best[2][0]}, t0={best[2][1]:.5e}")
    # Now sweep 200 extreme (a, b, w) on step 7 + 8 to try to break interlacing or RZ!
    P6 = best[1]
    breaks_interlace = 0
    breaks_rz = 0
    for tr in range(200):
        if mode == "double_Dm_w1":
            a7, b7, w7 = 1.0, 1.0 + 10**rng.uniform(-4, 6), 1.0
            a8, b8, w8 = 1.0, 1.0 + 10**rng.uniform(-4, 6), 1.0
        else:
            a7 = 10**rng.uniform(-5, 2); b7 = a7 + 10**rng.uniform(-3, 6); w7 = 10**rng.uniform(-4, 6)
            a8 = 10**rng.uniform(-5, 2); b8 = a8 + 10**rng.uniform(-3, 6); w8 = 10**rng.uniform(-4, 6)
        P7 = step_double(P6, a7, b7, a7, b7, w7)
        mx = max(np.max(np.abs(p)) for p in P7); P7 = [p/mx for p in P7]
        for i in range(1, len(P7)):
            st, val = check_interlace_T(P7[i], P7[i-1])
            if st == "NO_INTERLACE": breaks_interlace += 1; break
            if st == "COMPLEX": breaks_rz += 1; break
        P8 = step_double(P7, a8, b8, a8, b8, w8)
        mx = max(np.max(np.abs(p)) for p in P8); P8 = [p/mx for p in P8]
        for i in range(1, len(P8)):
            st, val = check_interlace_T(P8[i], P8[i-1])
            if st == "COMPLEX": breaks_rz += 1; break
    print(f"    After 200 extreme step 7+8 attacks: breaks_interlace={breaks_interlace}, breaks_rz(float64)={breaks_rz}")
