import math
import random
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

# === TEST 1: For which k does Q_{k,i}^{(c)}(t, y) = ht_{k,i-1} + c*ht_{k,i}*y + ht_{k,i+1}*y^2 fail at t=0?! ===
print("=== TEST 1: Exact t=0 interlacing test for quadratic ht_{k,i-1} + c*ht_{k,i}*y + ht_{k,i+1}*y^2 up to k=60 ===")
c0 = {(1, 0): 1, (1, 1): 1}
c1 = {(1, 0): 2, (1, 1): 0}
for k in range(2, 61):
    j = k + 1
    for i in range(k + 1):
        a0 = c0[(k-1, i)] if i <= k - 1 else 0
        a1 = (2 * c1[(k-1, i)] + j * c0[(k-1, i)]) if i <= k - 1 else 0
        s0 = a0 + (c0[(k-1, i-1)] if i >= 1 else 0)
        s1 = a1 + (c1[(k-1, i-1)] if i >= 1 else 0)
        c0[(k, i)] = s0
        c1[(k, i)] = 2 * s1 + j * s0

def quad_interlaces(a0, b0, d0, a1, b1, d1):
    # Check if P0(y) = a0 + b0*y + d0*y^2 and P1(y) = a1 + b1*y + d1*y^2 have real roots and P1 << P0
    # Using exact ints or high-precision floats (quadratics!)
    disc0 = b0*b0 - 4*a0*d0
    disc1 = b1*b1 - 4*a1*d1
    if disc0 <= 0 or disc1 <= 0:
        return False, "nonreal"
    # P1 << P0 (same deg 2, positive lc) iff resultant/Bezout condition:
    # r1_small < r0_small < r1_large < r0_large
    # Use float carefully via log or exact sqrt? With Python ints, use math.isqrt or Decimal!
    from decimal import Decimal, getcontext
    getcontext().prec = 80
    A0, B0, D0 = Decimal(a0), Decimal(b0), Decimal(d0)
    A1, B1, D1 = Decimal(a1), Decimal(b1), Decimal(d1)
    s0 = (B0*B0 - 4*A0*D0).sqrt()
    s1 = (B1*B1 - 4*A1*D1).sqrt()
    r0_lo, r0_hi = (-B0 - s0) / (2*D0), (-B0 + s0) / (2*D0)
    r1_lo, r1_hi = (-B1 - s1) / (2*D1), (-B1 + s1) / (2*D1)
    if r1_lo < r0_lo < r1_hi < r0_hi:
        return True, "ok"
    return False, f"order({float(r1_lo):.4f},{float(r0_lo):.4f},{float(r1_hi):.4f},{float(r0_hi):.4f})"

for c_num, c_den in [(2, 1), (21, 10), (22, 10), (5, 2), (3, 1)]:
    first_fail = None
    for k in range(2, 61):
        for i in range(1, k):
            ok, msg = quad_interlaces(
                c_den * c0[(k, i-1)], c_num * c0[(k, i)], c_den * c0[(k, i+1)],
                c_den * c1[(k, i-1)], c_num * c1[(k, i)], c_den * c1[(k, i+1)],
            )
            if not ok:
                first_fail = (k, i, msg)
                break
        if first_fail: break
    print(f"  c={c_num}/{c_den}: quad P1<<P0 up to k=60 -> {'PASS ALL!' if first_fail is None else f'FAIL at k={first_fail[0]}, i={first_fail[1]} ({first_fail[2]})'}")

# === TEST 2: Does Theorem A fail if C < 1 (e.g. C = 0.5, 0.8, 0.95)? ===
print("\n=== TEST 2: Is C=1.0 sharp for Theorem A (p1^2 >= C*p0*p2 on (-1,0) => T p2 + w p1 << T p1 + w p0)? ===")
t_grid = np.linspace(-1.0, 0.0, 501)
def roots_real(p):
    r = np.roots(p[::-1])
    if np.max(np.abs(r.imag)) > 1e-6: return None
    return np.sort(r.real)
def interlaces(p, q):
    rp, rq = roots_real(p), roots_real(q)
    if rp is None or rq is None: return False
    return len(rq) == len(rp) + 1 and all(rq[m] < rp[m] < rq[m+1] for m in range(len(rp)))

random.seed(55)
for C in [0.1, 0.5, 0.8, 0.95, 1.0]:
    fails = 0
    for trial in range(400):
        d = random.choice([0, 1, 2])
        r0 = sorted([-random.uniform(0.01, 0.99) for _ in range(d + 2)])
        r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 1)]
        r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d)]
        p0, p1, p2 = poly_from_roots(r0), poly_from_roots(r1), poly_from_roots(r2)
        v0, v1, v2 = eval_poly(p0, t_grid), eval_poly(p1, t_grid), eval_poly(p2, t_grid)
        pos = (v0 * v2) > 1e-30
        if not np.any(pos): continue
        mr = np.max((v0[pos] * v2[pos]) / (v1[pos]**2))
        p2 = p2 / (C * mr)
        a = 1.0
        b = a + random.choice([0.5, 2.0, 20.0, 500.0])
        w = 10 ** random.uniform(-3, 3)
        q2 = add_np(apply_T(a, b, p2), w * p1)
        q1 = add_np(apply_T(a, b, p1), w * p0)
        if not interlaces(q2, q1):
            fails += 1
    print(f"  C={C:4.2f}: fails(q2<<q1) = {fails:3d}/400")

# === TEST 3: Inspect WHY pure T failed p1^2 >= p0*p2 on (-1,0) in 25/500 trials! ===
print("\n=== TEST 3: Inspect a failure of pure T preserving p1^2 >= p0*p2 on (-1, 0) ===")
random.seed(11)
for trial in range(200):
    d = 1
    r0 = sorted([-random.uniform(0.02, 0.98) for _ in range(d + 2)])
    r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 1)]
    r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d)]
    p0, p1, p2 = poly_from_roots(r0), poly_from_roots(r1), poly_from_roots(r2)
    v0, v1, v2 = eval_poly(p0, t_grid), eval_poly(p1, t_grid), eval_poly(p2, t_grid)
    pos = (v0 * v2) > 1e-30
    if not np.any(pos): continue
    mr = np.max((v0[pos] * v2[pos]) / (v1[pos]**2))
    p2 = p2 / mr # tight C=1 on (-1, 0)
    a, b = 1.0, 50.0
    Tp0, Tp1, Tp2 = apply_T(a, b, p0), apply_T(a, b, p1), apply_T(a, b, p2)
    u0, u1, u2 = eval_poly(Tp0, t_grid), eval_poly(Tp1, t_grid), eval_poly(Tp2, t_grid)
    pos2 = (u0 * u2) > 1e-30
    if np.any(pos2) and np.min((u1[pos2]**2) / (u0[pos2]*u2[pos2])) < 0.8:
        idx = np.argmin((u1[pos2]**2) / (u0[pos2]*u2[pos2]))
        t_fail = t_grid[pos2][idx]
        # Also check if p1^2 >= p0*p2 held on ALL R!
        t_R = np.linspace(-20, 20, 2001)
        min_R = np.min(eval_poly(p1, t_R)**2 - eval_poly(p0, t_R)*eval_poly(p2, t_R))
        # Check if p0 + 2 p1 y + p2 y^2 was real-stable!
        print(f"  Found fail: r0={[round(z,3) for z in r0]}, r1={[round(z,3) for z in r1]}, r2={[round(z,3) for z in r2]}")
        print(f"    out min ratio = {np.min((u1[pos2]**2)/(u0[pos2]*u2[pos2])):.4f} at t={t_fail:.4f} | min(p1^2-p0p2) on R = {min_R:.4e}")
        # Check roots of p0, p1, p2 vs where p1^2 - p0*p2 was tight!
        idx0 = np.argmax((v0[pos]*v2[pos])/(v1[pos]**2))
        print(f"    input tightest at t={t_grid[pos][idx0]:.4f}")
        break

# === TEST 4: Check if Q(t,y) = p_0 + 2*p_1*y + p_2*y^2 real-stable on R? Wait:
# Check for WHICH c_num/c_den in TEST 1 does M_{b} preserve real stability of p_0 + c*p_1*y + p_2*y^2!
# Wait: check at c = 3 (i.e. p_1^2 >= (4/9) p_0 p_2) vs c = 2 (p_1^2 >= p_0 p_2)! ===
# Wait: EVEN BETTER! Check if p_0 + c*p_1*y + p_2*y^2 for GLZ has ALL REAL ROOTS in (-1,0) vs UHP stability!
# Let's check for c = 2.2 (where TEST 1 passes up to k=60!) whether Q_{k,i}^{(c)}(t, y) is REAL STABLE for GLZ k=2..12!
print("\n=== TEST 4: Check UHP real stability of quadratic Q_{k,i}^{(c)}(t, y) = ht_{k,i-1} + c*ht_{k,i}*y + ht_{k,i+1}*y^2 for GLZ ===")
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

def eval_c(p, z):
    acc = 0j
    for coef in reversed(p): acc = acc * z + coef
    return acc

test_ts = [complex(re, im) for re in [-10, -2, -1.05, -0.95, -0.5, -0.1, 0, 0.1, 1, 10] for im in [1e-4, 0.02, 0.5]]
for c_val in [2.0, 2.2, 2.5, 3.0]:
    bad = None
    for k in range(2, 13):
        for i in range(1, k):
            p0, p1, p2 = ht[(k, i-1)], ht[(k, i)], ht[(k, i+1)]
            for tv in test_ts:
                a0, a1, a2 = eval_c(p0, tv), c_val * eval_c(p1, tv), eval_c(p2, tv)
                # roots of a2 y^2 + a1 y + a0 = 0
                r = np.roots([a2, a1, a0])
                if max(z.imag for z in r) > 1e-6:
                    bad = (k, i, tv, max(z.imag for z in r))
                    break
            if bad: break
        if bad: break
    print(f"  c={c_val:3.1f}: Q_{{k,i}}^{{(c)}}(t,y) UHP stable k=2..12 -> {'PASS ALL!' if bad is None else f'FAIL at {bad}'}")
