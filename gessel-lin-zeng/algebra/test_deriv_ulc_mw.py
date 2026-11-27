import math
import random
import numpy as np
import sympy as sp

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

t_int = np.linspace(-1.0, 0.0, 501)
t_R = np.concatenate([-10.0**np.linspace(2, 0, 60) - 1.0, np.linspace(-2.0, -1.0, 60), t_int, np.linspace(0.0, 1.0, 60), 10.0**np.linspace(0, 2, 60)])

# === TEST 1: Does d/dt preserve p1^2 >= p0*p2 on R even with common root at t=0 (i.e. t*p_i)? ===
print("=== TEST 1: Does d/dt preserve f1^2 >= f0*f2 on R for f_i = p_i vs f_i = t*p_i vs f_i = t*(1+t)*p_i? ===")
random.seed(101)
for mode in ["plain p_i", "t * p_i", "t*(1+t) * p_i", "t*(1+t)^10 * p_i"]:
    fails = 0
    worst = 999.0
    for trial in range(400):
        d = random.choice([1, 2, 3])
        r0 = sorted([-random.uniform(0.02, 0.98) for _ in range(d + 2)])
        r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 1)]
        r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d)]
        p0, p1, p2 = poly_from_roots(r0), poly_from_roots(r1), poly_from_roots(r2)
        v0, v1, v2 = eval_poly(p0, t_R), eval_poly(p1, t_R), eval_poly(p2, t_R)
        pos = (v0 * v2) > 1e-30
        mr = max(1.0, float(np.max((v0[pos]*v2[pos])/(v1[pos]**2)))) if np.any(pos) else 1.0
        p2 = p2 / mr
        if mode == "plain p_i":
            f0, f1, f2 = p0, p1, p2
        elif mode == "t * p_i":
            w = np.array([0.0, 1.0])
            f0, f1, f2 = poly_mul(w, p0), poly_mul(w, p1), poly_mul(w, p2)
        elif mode == "t*(1+t) * p_i":
            w = np.array([0.0, 1.0, 1.0])
            f0, f1, f2 = poly_mul(w, p0), poly_mul(w, p1), poly_mul(w, p2)
        else:
            w = np.array([0.0, 1.0])
            for _ in range(10): w = poly_mul(w, np.array([1.0, 1.0]))
            f0, f1, f2 = poly_mul(w, p0), poly_mul(w, p1), poly_mul(w, p2)
        df0, df1, df2 = deriv_np(f0), deriv_np(f1), deriv_np(f2)
        u0, u1, u2 = eval_poly(df0, t_R), eval_poly(df1, t_R), eval_poly(df2, t_R)
        posd = (u0 * u2) > 1e-30
        if np.any(posd):
            md = float(np.min((u1[posd]**2)/(u0[posd]*u2[posd])))
            if md < worst: worst = md
            if md < 1.0 - 1e-5: fails += 1
    print(f"  {mode:18s}: fails(<1) = {fails:3d}/400 | worst min ratio = {worst:.6f}")

# === TEST 2: Check ULC propagation at MIDDLE index i = k//2 (where c_i ~ 1 + 4/k is SMALL!) ===
print("\n=== TEST 2: Test 4-tuple shift-sum with ULC at MIDDLE index i = k//2 for k=6, 10, 20 ===")
random.seed(202)
for k in [6, 10, 20]:
    i = k // 2
    c1_in = (math.comb(k, i-1)**2) / (math.comb(k, i-2) * math.comb(k, i))
    c2_in = (math.comb(k, i)**2) / (math.comb(k, i-1) * math.comb(k, i+1))
    c12_in = (math.comb(k, i-1)*math.comb(k, i)) / (math.comb(k, i-2)*math.comb(k, i+1))
    c2_out = (math.comb(k+1, i)**2) / (math.comb(k+1, i-1) * math.comb(k+1, i+1))
    fails_ulc = 0
    fails_1 = 0
    worst = 999.0
    for trial in range(300):
        d = random.choice([0, 1, 2])
        r0 = sorted([-random.uniform(0.02, 0.98) for _ in range(d + 3)])
        r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 2)]
        r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d + 1)]
        r3 = [random.uniform(r2[m]+1e-4, r2[m+1]-1e-4) for m in range(d)]
        p0, p1, p2, p3 = poly_from_roots(r0), poly_from_roots(r1), poly_from_roots(r2), poly_from_roots(r3)
        v = [eval_poly(p, t_int) for p in [p0, p1, p2, p3]]
        def mr(ai, bi, ci, di):
            num, den = v[ai]*v[bi], v[ci]*v[di]
            mask = (num > 1e-30) & (den > 1e-30)
            return float(np.max(num[mask]/den[mask])) if np.any(mask) else 1e-12
        s2 = 1.0 / (c1_in * mr(0, 2, 1, 1))
        p2 *= s2; v[2] *= s2
        s3 = min(1.0 / (c2_in * mr(1, 3, 2, 2)), 1.0 / (c12_in * mr(0, 3, 1, 2)))
        p3 *= s3
        a, b = 1.0, 1.0 + random.choice([1.0, 5.0, 50.0])
        w = 10 ** random.uniform(-1, 1)
        q1 = add_np(apply_T(a, b, p1), w * p0)
        q2 = add_np(apply_T(a, b, p2), w * p1)
        q3 = add_np(apply_T(a, b, p3), w * p2)
        u1, u2, u3 = eval_poly(q1, t_int), eval_poly(q2, t_int), eval_poly(q3, t_int)
        pos = (u1 * u3) > 1e-30
        if np.any(pos):
            m = float(np.min((u2[pos]**2)/(u1[pos]*u3[pos])))
            if m / c2_out < worst: worst = m / c2_out
            if m < c2_out - 1e-5: fails_ulc += 1
            if m < 1.0 - 1e-5: fails_1 += 1
    print(f"  k={k:2d}, i={i:2d} (c2_in={c2_in:.4f} -> c2_out={c2_out:.4f}): fail_ULC={fails_ulc:3d}/300 | fail(<1)={fails_1:3d}/300 | worst/c2_out={worst:.4f}")

# === TEST 3: Search for exact differential relations (generalizing Ma-Wang Prop. 1) on htilde_{k, i}! ===
print("\n=== TEST 3: Search for exact differential/intertwining relations among htilde_{k, i} ===")
t = sp.Symbol('t')
def D_sym(j, p):
    return sp.expand(t*(1+t)*sp.diff(p, t) + (1 + j*t)*p)

ht = {(1, 0): 1 + 2*t, (1, 1): sp.Integer(1)}
for k in range(2, 7):
    for i in range(k + 1):
        t1 = D_sym(k + 1, ht[(k-1, i)]) if i <= k - 1 else 0
        t2 = ht[(k-1, i-1)] if i >= 1 else 0
        ht[(k, i)] = D_sym(k + 1, sp.expand(t1 + t2))

# Check Ma-Wang Prop 1 in t!
# A_{k, k-1}(x) = k*A_{k+1, k+1}(x) - (x/2)*A'_{k+1, k+1}(x)
# In t: ht_{k, k-1}(t) has deg k, ht_{k+1, k+1}(t) has deg k!
print("Check Ma-Wang Prop 1 in t (relating ht_{k, k-1} to ht_{k+1, k+1} and its derivative):")
for k in range(1, 6):
    # Solve ht_{k, k-1} = (a + b*t)*ht_{k+1,k+1} + (c + d*t + e*t^2)*ht'_{k+1,k+1}
    a, b, c, d, e = sp.symbols('a b c d e')
    f = ht[(k+1, k+1)]
    expr = sp.expand((a + b*t)*f + (c + d*t + e*t**2)*sp.diff(f, t) - ht[(k, k-1)])
    sol = sp.solve(sp.Poly(expr, t).coeffs(), [a, b, c, d, e])
    print(f"  k={k}: {sol}")

# Check if ht_{k, i-1} is a linear differential combination of ht_{k+1, i+1} or ht_{k, i}!
print("Check if ht_{k, i} = L(ht_{k+1, i+2}) (same degree 2k-i-1 = 2(k+1)-(i+2)-1!):")
for k in range(2, 6):
    for i in range(k - 1):
        a, b, c, d, e = sp.symbols('a b c d e')
        f = ht[(k+1, i+2)]
        expr = sp.expand((a + b*t)*f + (c + d*t + e*t**2)*sp.diff(f, t) - ht[(k, i)])
        sol = sp.solve(sp.Poly(expr, t).coeffs(), [a, b, c, d, e])
        if sol:
            print(f"  MATCH k={k}, i={i}: ht_{{{k},{i}}} = L(ht_{{{k+1},{i+2}}}) -> {sol}")

# Check if ht_{k, i-1} = L(ht_{k, i}) (deg increases by 1!)
print("Check if ht_{k, i-1} = (a + b*t)*ht_{k,i} + (c + d*t + e*t^2)*ht'_{k,i}:")
for k in range(2, 6):
    for i in range(1, k + 1):
        a, b, c, d, e = sp.symbols('a b c d e')
        f = ht[(k, i)]
        expr = sp.expand((a + b*t)*f + (c + d*t + e*t**2)*sp.diff(f, t) - ht[(k, i-1)])
        sol = sp.solve(sp.Poly(expr, t).coeffs(), [a, b, c, d, e])
        if sol:
            print(f"  MATCH k={k}, i={i}: ht_{{{k},{i-1}}} = L(ht_{{{k},{i}}}) -> {sol}")
