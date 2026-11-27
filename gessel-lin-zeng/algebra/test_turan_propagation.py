import random
import math
import numpy as np

def poly_mul(a, b):
    out = np.zeros(len(a) + len(b) - 1, dtype=float)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out

def add_np(a, b):
    n = max(len(a), len(b))
    out = np.zeros(n, dtype=float)
    out[:len(a)] += a
    out[:len(b)] += b
    return out

def apply_T(a, b, p):
    d = len(p) - 1
    out = np.zeros(d + 2, dtype=float)
    for m, cm in enumerate(p):
        out[m] += (m + a) * cm
        out[m + 1] += (m + b) * cm
    return out

def eval_poly(p, x):
    acc = np.zeros_like(x, dtype=float)
    for c in reversed(p):
        acc = acc * x + c
    return acc

def poly_from_roots(rs):
    p = np.array([1.0])
    for r in rs:
        p = poly_mul(p, np.array([-r, 1.0]))
    return p

t_grid = np.linspace(-1.0, 0.0, 601)
t_grid_R = np.concatenate([np.linspace(-20.0, -1.0, 200), t_grid, np.linspace(0.0, 20.0, 200)])

# === TEST 1: Does T preserve p1^2 >= C * p0 * p2 on (-1, 0)?! ===
print("=== TEST 1: Does pure T preserve p1^2 >= C * p0 * p2 on (-1, 0)? ===")
random.seed(11)
for C in [1.0, 1.5, 2.0]:
    fail_T_sameC = 0
    fail_T_C1 = 0
    worst_ratio = 999.0
    for trial in range(500):
        d = random.choice([0, 1, 2, 3])
        r0 = sorted([-random.uniform(0.01, 0.99) for _ in range(d + 2)])
        r1 = [random.uniform(r0[m] + 1e-4, r0[m+1] - 1e-4) for m in range(d + 1)]
        r2 = [random.uniform(r1[m] + 1e-4, r1[m+1] - 1e-4) for m in range(d)]
        p0, p1, p2 = poly_from_roots(r0), poly_from_roots(r1), poly_from_roots(r2)
        v0, v1, v2 = eval_poly(p0, t_grid), eval_poly(p1, t_grid), eval_poly(p2, t_grid)
        pos = (v0 * v2) > 0
        if not np.any(pos): continue
        mr = np.max((v0[pos] * v2[pos]) / (v1[pos] ** 2))
        # Set tight equality at max point!
        p2 = p2 / (C * mr)
        a = random.choice([1.0, 2.0, 5.0])
        b = a + random.choice([0.5, 1.0, 5.0, 50.0, 500.0])
        Tp0, Tp1, Tp2 = apply_T(a, b, p0), apply_T(a, b, p1), apply_T(a, b, p2)
        u0, u1, u2 = eval_poly(Tp0, t_grid), eval_poly(Tp1, t_grid), eval_poly(Tp2, t_grid)
        pos2 = (u0 * u2) > 1e-30
        if np.any(pos2):
            min_out_C = np.min((u1[pos2] ** 2) / (u0[pos2] * u2[pos2]))
            if min_out_C < worst_ratio: worst_ratio = min_out_C
            if min_out_C < C - 1e-6: fail_T_sameC += 1
            if min_out_C < 1.0 - 1e-6: fail_T_C1 += 1
    print(f"  Input C={C:3.1f}: fail(out < C)={fail_T_sameC:3d}/500 | fail(out < 1)={fail_T_C1:3d}/500 | worst out_C = {worst_ratio:.6f}")

# === TEST 2: Check EXACT minimum of htilde_{k,i}(t)^2 / (htilde_{k,i-1}(t) htilde_{k,i+1}(t)) and ULC ratio for GLZ! ===
print("\n=== TEST 2: Check Turan ratio and ULC ratio of GLZ htilde_{k,i}(t) on (-1, 0) and on R ===")
def apply_D_int(j, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return out
def add_int(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return out

ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 13):
    for i in range(k + 1):
        t1 = apply_D_int(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D_int(k + 1, add_int(t1, t2))

for k in range(2, 13):
    min_turan_int = 999.0
    min_ulc_int = 999.0
    min_ulc_at = None
    min_ulc_R = 999.0
    for i in range(1, k):
        p0 = np.array([float(c) for c in ht[(k, i-1)]])
        p1 = np.array([float(c) for c in ht[(k, i)]])
        p2 = np.array([float(c) for c in ht[(k, i+1)]])
        v0, v1, v2 = eval_poly(p0, t_grid), eval_poly(p1, t_grid), eval_poly(p2, t_grid)
        pos = (v0 * v2) > 1e-30
        ulc_fac = (math.comb(k, i)**2) / (math.comb(k, i-1) * math.comb(k, i+1))
        if np.any(pos):
            ratios = (v1[pos]**2) / (v0[pos] * v2[pos])
            mr = float(np.min(ratios))
            if mr < min_turan_int: min_turan_int = mr
            if mr / ulc_fac < min_ulc_int:
                min_ulc_int = mr / ulc_fac
                idx = np.argmin(ratios)
                t_at = float(t_grid[pos][idx])
                min_ulc_at = (i, round(t_at, 4))
        # Also on R
        w0, w1, w2 = eval_poly(p0, t_grid_R), eval_poly(p1, t_grid_R), eval_poly(p2, t_grid_R)
        posR = (w0 * w2) > 1e-30
        if np.any(posR):
            mrR = float(np.min((w1[posR]**2) / (w0[posR] * w2[posR]))) / ulc_fac
            if mrR < min_ulc_R: min_ulc_R = mrR
    print(f"  k={k:2d}: min Turan on (-1,0) = {min_turan_int:.6f} | min ULC on (-1,0) = {min_ulc_int:.6f} at {min_ulc_at} | min ULC on R = {min_ulc_R:.6f}")

# === TEST 3: Check 4-polynomial interior step q_i^2 >= C_out * q_{i-1} q_{i+1} with ULC weights! ===
print("\n=== TEST 3: Test if q_i = T p_i + w p_{i-1} preserves ULC-like / Turan inequality on 4-tuple p3 << p2 << p1 << p0! ===")
# Check what happens to (q_1, q_2, q_3) when p_i = binom(k, i) * bar_p_i with bar_p_1^2 >= bar_p_0 bar_p_2 AND bar_p_2^2 >= bar_p_1 bar_p_3!
for k in [3, 5, 10]:
    # Suppose bar_p_i^2 >= bar_p_{i-1} bar_p_{i+1} on (-1, 0), AND bar_p_1 bar_p_2 >= bar_p_0 bar_p_3!
    # Wait: test on 300 random 4-tuples p3 << p2 << p1 << p0 scaled so p_i^2 >= C_i * p_{i-1} p_{i+1} (with ULC factor)!
    fails_ulc = 0
    fails_turan1 = 0
    worst_ulc_out = 999.0
    for trial in range(300):
        d = random.choice([0, 1, 2])
        r0 = sorted([-random.uniform(0.02, 0.98) for _ in range(d + 3)])
        r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 2)]
        r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d + 1)]
        r3 = [random.uniform(r2[m]+1e-4, r2[m+1]-1e-4) for m in range(d)]
        p0, p1, p2, p3 = poly_from_roots(r0), poly_from_roots(r1), poly_from_roots(r2), poly_from_roots(r3)
        # Scale p0, p1, p2, p3 to satisfy ULC for index (i-1, i, i+1, i+2) = (0, 1, 2, 3) at level k
        # Wait: let's test the SIMPLEST uniform Turan: p1^2 >= c1*p0*p2, p2^2 >= c2*p1*p3, p1*p2 >= c12*p0*p3!
        v = [eval_poly(p, t_grid) for p in [p0, p1, p2, p3]]
        def max_rat(a_idx, b_idx, c_idx, d_idx):
            # max (v[a]*v[b]) / (v[c]*v[d]) where numerator > 0 and denom > 0
            num = v[a_idx] * v[b_idx]
            den = v[c_idx] * v[d_idx]
            mask = (num > 1e-30) & (den > 1e-30)
            return float(np.max(num[mask] / den[mask])) if np.any(mask) else 1e-12
        # For ULC at level k, i=1 and i=2:
        c1_in = (math.comb(k, 1)**2) / (math.comb(k, 0) * math.comb(k, 2))
        c2_in = (math.comb(k, 2)**2) / (math.comb(k, 1) * math.comb(k, 3))
        # Scale p2 so p1^2 >= c1_in * p0 * p2
        s2 = 1.0 / (c1_in * max_rat(0, 2, 1, 1))
        p2 = p2 * s2
        v[2] = v[2] * s2
        # Scale p3 so p2^2 >= c2_in * p1 * p3 AND p1*p2 >= c12 * p0*p3
        c12_in = (math.comb(k, 1) * math.comb(k, 2)) / (math.comb(k, 0) * math.comb(k, 3))
        s3 = min(1.0 / (c2_in * max_rat(1, 3, 2, 2)), 1.0 / (c12_in * max_rat(0, 3, 1, 2)))
        p3 = p3 * s3
        # Apply step q_i = T p_i + w p_{i-1}
        a = random.choice([1.0, 2.0])
        b = a + random.choice([1.0, 5.0, 50.0])
        w = 10 ** random.uniform(-1, 1)
        q1 = add_np(apply_T(a, b, p1), w * p0)
        q2 = add_np(apply_T(a, b, p2), w * p1)
        q3 = add_np(apply_T(a, b, p3), w * p2)
        # Check q2^2 / (q1 * q3) vs ULC factor at level k+1, i=2!
        c2_out = (math.comb(k+1, 2)**2) / (math.comb(k+1, 1) * math.comb(k+1, 3))
        u1, u2, u3 = eval_poly(q1, t_grid), eval_poly(q2, t_grid), eval_poly(q3, t_grid)
        pos = (u1 * u3) > 1e-30
        if np.any(pos):
            mr = float(np.min((u2[pos]**2) / (u1[pos] * u3[pos])))
            if mr / c2_out < worst_ulc_out: worst_ulc_out = mr / c2_out
            if mr < c2_out - 1e-6: fails_ulc += 1
            if mr < 1.0 - 1e-6: fails_turan1 += 1
    print(f"  k={k:2d} (ULC c2_out={c2_out:.4f}): fail(ULC_out)={fails_ulc:3d}/300 | fail(Turan>=1)={fails_turan1:3d}/300 | worst ratio/c2_out = {worst_ulc_out:.6f}")
