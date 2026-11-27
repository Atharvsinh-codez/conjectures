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

def deriv_np(p):
    return np.array([m * p[m] for m in range(1, len(p))], dtype=float)

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

# Dense grid on R (including far left/right!)
t_int = np.linspace(-1.0, 0.0, 501)
t_R = np.concatenate([
    -10.0**np.linspace(3, 0, 100) - 1.0,
    np.linspace(-2.0, -1.0, 100),
    t_int,
    np.linspace(0.0, 1.0, 100),
    10.0**np.linspace(0, 3, 100),
])

# === TEST 1: If p2 << p1 << p0 in (-1, 0) AND p1^2 >= p0*p2 on ALL R, does (Tp1)^2 >= (Tp0)(Tp2) hold on (-1,0) and on R?! ===
print("=== TEST 1: Does p2<<p1<<p0 + (p1^2 >= p0*p2 on ALL R) imply (Tp1)^2 >= (Tp0)(Tp2)? ===")
random.seed(77)
for C in [1.0, 1.2, 1.5, 2.0]:
    fail_int = 0
    fail_R = 0
    fail_deriv = 0
    worst_int = 999.0
    worst_R = 999.0
    worst_deriv = 999.0
    for trial in range(500):
        d = random.choice([0, 1, 2, 3])
        r0 = sorted([-random.uniform(0.02, 0.98) for _ in range(d + 2)])
        r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 1)]
        r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d)]
        p0, p1, p2 = poly_from_roots(r0), poly_from_roots(r1), poly_from_roots(r2)
        v0, v1, v2 = eval_poly(p0, t_R), eval_poly(p1, t_R), eval_poly(p2, t_R)
        pos = (v0 * v2) > 1e-30
        if not np.any(pos): continue
        # Also include t -> infinity ratio (leading coeff ratio = 1 since all monic!)
        mr = max(1.0, float(np.max((v0[pos] * v2[pos]) / (v1[pos]**2))))
        p2 = p2 / (C * mr)

        # Test plain derivative (f1')^2 >= f0' f2' on R!
        if d >= 1:
            dp0, dp1, dp2 = deriv_np(p0), deriv_np(p1), deriv_np(p2)
            w0, w1, w2 = eval_poly(dp0, t_R), eval_poly(dp1, t_R), eval_poly(dp2, t_R)
            posd = (w0 * w2) > 1e-30
            if np.any(posd):
                md = float(np.min((w1[posd]**2) / (w0[posd]*w2[posd])))
                if md < worst_deriv: worst_deriv = md
                if md < C - 1e-5: fail_deriv += 1

        # Test T!
        a = random.choice([1.0, 2.0, 5.0])
        b = a + random.choice([0.5, 1.0, 5.0, 50.0, 500.0])
        Tp0, Tp1, Tp2 = apply_T(a, b, p0), apply_T(a, b, p1), apply_T(a, b, p2)
        u0, u1, u2 = eval_poly(Tp0, t_R), eval_poly(Tp1, t_R), eval_poly(Tp2, t_R)
        posR = (u0 * u2) > 1e-30
        if np.any(posR):
            mR = float(np.min((u1[posR]**2) / (u0[posR]*u2[posR])))
            if mR < worst_R: worst_R = mR
            if mR < C - 1e-5: fail_R += 1
        u0i, u1i, u2i = eval_poly(Tp0, t_int), eval_poly(Tp1, t_int), eval_poly(Tp2, t_int)
        posi = (u0i * u2i) > 1e-30
        if np.any(posi):
            mi = float(np.min((u1i[posi]**2) / (u0i[posi]*u2i[posi])))
            if mi < worst_int: worst_int = mi
            if mi < C - 1e-5: fail_int += 1
    print(f"  C={C:3.1f}: fail_deriv(<C)={fail_deriv:3d} (min={worst_deriv:.4f}) | fail_T_int(<C)={fail_int:3d} (min={worst_int:.4f}) | fail_T_R(<C)={fail_R:3d} (min={worst_R:.4f})")

# === TEST 2: What about the SHIFT-SUM q_i = T p_i + w p_{i-1} when input satisfies p_i^2 >= C p_{i-1} p_{i+1} on ALL R
#             AND p_1 p_2 >= C p_0 p_3 on ALL R (4-tuple)?! ===
print("\n=== TEST 2: 4-tuple interior shift-sum q_i = T p_i + w p_{i-1} with input Turan on ALL R ===")
random.seed(88)
for C in [1.0, 1.1, 1.5, 2.0]:
    fail_q2_R = 0
    fail_q2_int = 0
    worst_q2_R = 999.0
    worst_q2_int = 999.0
    for trial in range(300):
        d = random.choice([0, 1, 2])
        r0 = sorted([-random.uniform(0.02, 0.98) for _ in range(d + 3)])
        r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 2)]
        r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d + 1)]
        r3 = [random.uniform(r2[m]+1e-4, r2[m+1]-1e-4) for m in range(d)]
        p0, p1, p2, p3 = poly_from_roots(r0), poly_from_roots(r1), poly_from_roots(r2), poly_from_roots(r3)
        v = [eval_poly(p, t_R) for p in [p0, p1, p2, p3]]
        def mr_R(a_i, b_i, c_i, d_i):
            num, den = v[a_i]*v[b_i], v[c_i]*v[d_i]
            mask = (num > 1e-30) & (den > 1e-30)
            return max(1.0, float(np.max(num[mask]/den[mask]))) if np.any(mask) else 1.0
        s2 = 1.0 / (C * mr_R(0, 2, 1, 1))
        p2 = p2 * s2; v[2] *= s2
        s3 = min(1.0 / (C * mr_R(1, 3, 2, 2)), 1.0 / (C * mr_R(0, 3, 1, 2)))
        p3 = p3 * s3; v[3] *= s3
        a = random.choice([1.0, 2.0])
        b = a + random.choice([0.5, 2.0, 20.0, 200.0])
        w = 10 ** random.uniform(-2, 2)
        q1 = add_np(apply_T(a, b, p1), w * p0)
        q2 = add_np(apply_T(a, b, p2), w * p1)
        q3 = add_np(apply_T(a, b, p3), w * p2)
        u1, u2, u3 = eval_poly(q1, t_R), eval_poly(q2, t_R), eval_poly(q3, t_R)
        posR = (u1 * u3) > 1e-30
        if np.any(posR):
            mR = float(np.min((u2[posR]**2) / (u1[posR]*u3[posR])))
            if mR < worst_q2_R: worst_q2_R = mR
            if mR < C - 1e-5: fail_q2_R += 1
        u1i, u2i, u3i = eval_poly(q1, t_int), eval_poly(q2, t_int), eval_poly(q3, t_int)
        posi = (u1i * u3i) > 1e-30
        if np.any(posi):
            mi = float(np.min((u2i[posi]**2) / (u1i[posi]*u3i[posi])))
            if mi < worst_q2_int: worst_q2_int = mi
            if mi < C - 1e-5: fail_q2_int += 1
    print(f"  C={C:3.1f}: fail_q2_int(<C)={fail_q2_int:3d} (min={worst_q2_int:.4f}) | fail_q2_R(<C)={fail_q2_R:3d} (min={worst_q2_R:.4f})")

# === TEST 3: Check GLZ at k=35..40: does htilde_{k,i}^2 > htilde_{k,i-1} htilde_{k,i+1} STILL hold on (-1,0) and R?
# (Fast check at t=0 and t=-1 and leading coeffs!) ===
print("\n=== TEST 3: Check GLZ at k=35: leading coeff Turan ratio and t=-1 Turan ratio ===")
# Wait: at t -> infinity, lc(htilde_{k, i})!
# Recurrence for lc(htilde_{k, i}):
# deg(ht_{k, i}) = 2k-i-1. D_{k+1} multiplies lc by (deg + k + 1) = (3k - i)!!
# Let's check: D_{k+1} ht_{k-1, i} has deg 2k-i-2, lc = (2k-i-3 + k+1)*lc(ht_{k-1, i}) = (3k-i-2)*lc(ht_{k-1, i}).
# ht_{k-1, i-1} has deg 2k-i-2!
# Then outer D_{k+1} multiplies by (2k-i-2 + k+1) = (3k-i-1)!
lc = {(1, 0): 2, (1, 1): 1}
for k in range(2, 41):
    for i in range(k + 1):
        t1 = (3*k - i - 2) * lc[(k-1, i)] if i <= k - 1 else 0
        t2 = lc[(k-1, i-1)] if i >= 1 else 0
        lc[(k, i)] = (3*k - i - 1) * (t1 + t2)
for k in [5, 10, 20, 35, 40]:
    min_rat = min((lc[(k, i)]**2) / (lc[(k, i-1)] * lc[(k, i+1)]) for i in range(1, k))
    min_ulc = min(((lc[(k, i)]**2) / (lc[(k, i-1)] * lc[(k, i+1)])) / ((math.comb(k,i)**2)/(math.comb(k,i-1)*math.comb(k,i+1))) for i in range(1, k))
    print(f"  k={k:2d}: min lc Turan = {min_rat:.6f} | min lc ULC = {min_ulc:.6f}")
