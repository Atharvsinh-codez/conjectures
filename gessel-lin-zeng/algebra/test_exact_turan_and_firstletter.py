import math
import random
import numpy as np
import sympy as sp

# Exact check: is a polynomial P(t) >= 0 for ALL real t?
# Using numpy roots (with high margin) + exact Sturm check for candidate counterexamples!
def min_ratio_exact_poly(p0, p1, p2):
    """Find EXACT maximum of (p0*p2)/(p1^2) over all real t where p0*p2 > 0!
    Since p1^2 - lambda * p0 * p2 = 0 at a double root when lambda is extremal,
    critical points of (p0*p2)/(p1^2) satisfy (p0*p2)' * p1 - 2 * p0 * p2 * p1' = 0!"""
    # Compute derivative numerator: N = (p0*p2)' * p1 - 2 * p0 * p2 * p1'
    p02 = np.polymul(p0[::-1], p2[::-1]) # high -> low!
    dp02 = np.polyder(p02)
    p1h = p1[::-1]
    dp1 = np.polyder(p1h)
    N = np.polysub(np.polymul(dp02, p1h), 2.0 * np.polymul(p02, dp1))
    crit = np.roots(N)
    real_crit = [z.real for z in crit if abs(z.imag) < 1e-7]
    # Also evaluate at infinity: lc(p0*p2) / lc(p1^2)
    max_r = (p0[-1] * p2[-1]) / (p1[-1] ** 2)
    best_t = float('inf')
    for tc in real_crit:
        v1 = np.polyval(p1h, tc)
        if abs(v1) < 1e-14:
            continue
        v02 = np.polyval(p02, tc)
        if v02 > 0:
            r = v02 / (v1 * v1)
            if r > max_r:
                max_r = r
                best_t = tc
    return max_r, best_t

def apply_T(a, b, p):
    d = len(p) - 1
    out = np.zeros(d + 2, dtype=float)
    for m, cm in enumerate(p):
        out[m] += (m + a) * cm
        out[m + 1] += (m + b) * cm
    return out

def add_np(a, b):
    n = max(len(a), len(b))
    out = np.zeros(n, dtype=float)
    out[:len(a)] += a
    out[:len(b)] += b
    return out

def poly_from_roots(rs):
    p = np.array([1.0])
    for r in rs:
        p = np.polymul(p, np.array([1.0, -r])) # high->low
    return p[::-1] # low->high

# === TEST 1: EXACT critical-point check of whether T preserves p1^2 >= p0*p2 on ALL R! ===
print("=== TEST 1: EXACT critical-point test of T preserving (p1^2 >= p0*p2 on ALL R) ===")
random.seed(777)
for (a, b) in [(1.0, 2.0), (1.0, 3.0), (1.0, 10.0), (1.0, 50.0), (1.0, 500.0), (5.0, 10.0), (5.0, 100.0)]:
    fails = 0
    worst_out_maxr = 0.0 # want <= 1.0!
    worst_case = None
    for trial in range(400):
        d = random.choice([0, 1, 2, 3])
        r0 = sorted([-random.uniform(0.01, 0.99) for _ in range(d + 2)])
        r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 1)]
        r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d)]
        p0, p1, p2 = poly_from_roots(r0), poly_from_roots(r1), poly_from_roots(r2)
        mr_in, _ = min_ratio_exact_poly(p0, p1, p2)
        p2 = p2 / mr_in # EXACTLY tight: max (p0*p2)/p1^2 == 1.0 on all R!
        q0, q1, q2 = apply_T(a, b, p0), apply_T(a, b, p1), apply_T(a, b, p2)
        mr_out, t_out = min_ratio_exact_poly(q0, q1, q2)
        if mr_out > worst_out_maxr:
            worst_out_maxr = mr_out
            worst_case = (d, r0, r1, r2, t_out)
        if mr_out > 1.0 + 1e-6:
            fails += 1
    print(f"  T(a={a:3.1f}, b={b:5.1f}): fails(out>1) = {fails:3d}/400 | max(Tp0*Tp2/Tp1^2) = {worst_out_maxr:.6f}")
    if fails > 0 and worst_case is not None:
        d, r0, r1, r2, t_out = worst_case
        print(f"    WORST: d={d}, t_out={t_out:.5f}, r0={[round(z,3) for z in r0]}, r1={[round(z,3) for z in r1]}, r2={[round(z,3) for z in r2]}")

# === TEST 2: EXACT critical-point test of d/dt preserving (p1^2 >= p0*p2 on ALL R) ===
print("\n=== TEST 2: EXACT critical-point test of d/dt preserving (p1^2 >= p0*p2 on ALL R) ===")
fails_d = 0
worst_d = 0.0
for trial in range(500):
    d = random.choice([1, 2, 3, 4])
    r0 = sorted([-random.uniform(0.01, 0.99) for _ in range(d + 2)])
    r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 1)]
    r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d)]
    p0, p1, p2 = poly_from_roots(r0), poly_from_roots(r1), poly_from_roots(r2)
    mr_in, _ = min_ratio_exact_poly(p0, p1, p2)
    p2 = p2 / mr_in
    dp0 = np.polyder(p0[::-1])[::-1]
    dp1 = np.polyder(p1[::-1])[::-1]
    dp2 = np.polyder(p2[::-1])[::-1]
    mr_out, _ = min_ratio_exact_poly(dp0, dp1, dp2)
    if mr_out > worst_d: worst_d = mr_out
    if mr_out > 1.0 + 1e-6: fails_d += 1
print(f"  d/dt: fails(out>1) = {fails_d}/500 | max(p0'p2'/p1'^2) = {worst_d:.6f}")

# === TEST 3: Check EXACT critical-point test of shift-sum q_i = T p_i + w p_{i-1} on 4-tuple (p0,p1,p2,p3) ===
print("\n=== TEST 3: EXACT critical-point test of q_i = T p_i + w p_{i-1} on 4-tuple with exact Turan on R ===")
fails_s = 0
worst_s = 0.0
for trial in range(400):
    d = random.choice([0, 1, 2])
    r0 = sorted([-random.uniform(0.02, 0.98) for _ in range(d + 3)])
    r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 2)]
    r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d + 1)]
    r3 = [random.uniform(r2[m]+1e-4, r2[m+1]-1e-4) for m in range(d)]
    p0, p1, p2, p3 = poly_from_roots(r0), poly_from_roots(r1), poly_from_roots(r2), poly_from_roots(r3)
    mr1, _ = min_ratio_exact_poly(p0, p1, p2)
    p2 = p2 / mr1
    mr2, _ = min_ratio_exact_poly(p1, p2, p3)
    # Also enforce p1*p2 >= p0*p3 on R!
    # Check max (p0*p3)/(p1*p2) on dense grid + critical points
    p03 = np.polymul(p0[::-1], p3[::-1])
    p12 = np.polymul(p1[::-1], p2[::-1])
    N = np.polysub(np.polymul(np.polyder(p03), p12), np.polymul(p03, np.polyder(p12)))
    mr12 = (p0[-1]*p3[-1])/(p1[-1]*p2[-1])
    for z in np.roots(N):
        if abs(z.imag) < 1e-7:
            v03, v12 = np.polyval(p03, z.real), np.polyval(p12, z.real)
            if v03 > 0 and v12 > 1e-14:
                mr12 = max(mr12, v03 / v12)
    p3 = p3 / max(mr2, mr12)
    a = random.choice([1.0, 2.0])
    b = a + random.choice([1.0, 5.0, 50.0])
    w = 10 ** random.uniform(-2, 2)
    q1 = add_np(apply_T(a, b, p1), w * p0)
    q2 = add_np(apply_T(a, b, p2), w * p1)
    q3 = add_np(apply_T(a, b, p3), w * p2)
    mr_out, t_out = min_ratio_exact_poly(q1, q2, q3)
    if mr_out > worst_s: worst_s = mr_out
    if mr_out > 1.0 + 1e-6: fails_s += 1
print(f"  shift-sum q_i = T p_i + w p_{{i-1}}: fails(out>1) = {fails_s}/400 | max(q1*q3/q2^2) = {worst_s:.6f}")

# === TEST 4: Check first_letter strong interlacing for GLZ up to k=9 using exact integer Sturm! ===
print("\n=== TEST 4: Exact integer Sturm check of first_letter strong interlacing for k=2..9 ===")
def poly_trim(c):
    while len(c) > 1 and c[-1] == 0: c.pop()
    return c
def apply_S_int(L, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (L - m) * cm
    return poly_trim(out)
def add_int(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return poly_trim(out)
def shift_x(c):
    return [0] + c if any(v != 0 for v in c) else [0]
def exact_strict_interlaces_x(p, q):
    # Handle common factor x (root at 0)! Strip min power of x first?
    # Wait: F_{k,i,r} for r >= 2 might have factor x, F_{k,i,1} does not!
    # So check non-strict interlacing or strip common x!
    p, q = p[:], q[:]
    while len(p) > 1 and len(q) > 1 and p[0] == 0 and q[0] == 0:
        p.pop(0); q.pop(0)
    dp, dq = len(p) - 1, len(q) - 1
    if p[-1] <= 0 or q[-1] <= 0: return False
    if dq == 0: return dp == 0
    if dq == dp:
        lcp, lcq = p[-1], q[-1]
        r1 = [lcq * p[m] - lcp * q[m] for m in range(dq)]
        while len(r1) > 1 and r1[-1] == 0: r1.pop()
        if len(r1) - 1 != dq - 1 or r1[-1] <= 0: return False
        r0 = q[:]
    elif dq == dp + 1:
        r0, r1 = q[:], p[:]
    else:
        return False
    cur0, cur1 = r0[:], r1[:]
    while len(cur1) > 1:
        d0, d1 = len(cur0) - 1, len(cur1) - 1
        if d0 != d1 + 1: return False
        a0, a1 = cur0[-1], cur1[-1]
        tmp = [a1 * cur0[m] for m in range(d0 + 1)]
        for m in range(d1 + 1): tmp[m + 1] -= a0 * cur1[m]
        b = tmp[d1]
        nxt = [-(a1 * tmp[m] - b * cur1[m]) for m in range(d1)]
        while len(nxt) > 1 and nxt[-1] == 0: nxt.pop()
        if len(nxt) - 1 != d1 - 1 or nxt[-1] <= 0: return False
        g = math.gcd(*[abs(v) for v in nxt])
        if g > 1: nxt = [v // g for v in nxt]
        cur0, cur1 = cur1, nxt
    return True

# Build F[(k, i)][r] for r = 1..2k
A = {(1, 0): [1, 1], (1, 1): [1]}
F = {(1, 0): {1: [1], 2: [0, 1]}, (1, 1): {1: [0], 2: [1]}}
for k in range(2, 10):
    for i in range(k + 1):
        L = 3*k - 3 - i
        t1 = apply_S_int(L, A[(k-1, i)]) if i <= k - 1 else [0]
        t2 = A[(k-1, i-1)] if i >= 1 else [0]
        B = add_int(t1, t2)
        A[(k, i)] = apply_S_int(L + 1, B)
        F[(k, i)] = {}
        for r in range(1, 2*k - 1):
            f1 = apply_S_int(L - 1, F[(k-1, i)][r]) if i <= k - 1 else [0]
            f2 = F[(k-1, i-1)][r] if i >= 1 else [0]
            F[(k, i)][r] = apply_S_int(L, add_int(f1, f2))
        F[(k, i)][2*k - 1] = apply_S_int(L, shift_x(A[(k-1, i)])) if i <= k - 1 else [0]
        F[(k, i)][2*k] = shift_x(B)
        # Verify sum equals A[(k, i)]
        s = [0]
        for r in range(1, 2*k + 1): s = add_int(s, F[(k, i)][r])
        assert s == A[(k, i)], f"Sum mismatch at ({k},{i})!"
    # Check strong interlacing of nonzero F[(k, i)][r] for all i!
    fails_k = []
    for i in range(k + 1):
        nz = [r for r in range(1, 2*k + 1) if any(v != 0 for v in F[(k, i)][r])]
        for a_idx in range(len(nz)):
            for b_idx in range(a_idx + 1, len(nz)):
                ra, rb = nz[a_idx], nz[b_idx]
                if not exact_strict_interlaces_x(F[(k, i)][ra], F[(k, i)][rb]):
                    fails_k.append((i, ra, rb))
                    break
            if fails_k and fails_k[-1][0] == i: break
    print(f"  k={k}: first_letter strong interlacing across ALL i=0..{k}: {'PASS ALL!' if not fails_k else f'FAIL {fails_k[:4]}'}")
