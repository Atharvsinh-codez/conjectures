import math
import random
import numpy as np

def poly_trim(c):
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c

def apply_D(j, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return poly_trim(out)

def add_poly(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return poly_trim(out)

def exact_strict_interlaces(p, q):
    dp, dq = len(p) - 1, len(q) - 1
    if p[-1] <= 0 or q[-1] <= 0:
        return False, "nonpos_lc"
    if dq == 0:
        return dp == 0, "deg0"
    if dq == dp:
        lcp, lcq = p[-1], q[-1]
        r1 = [lcq * p[m] - lcp * q[m] for m in range(dq)]
        while len(r1) > 1 and r1[-1] == 0:
            r1.pop()
        if len(r1) - 1 != dq - 1 or r1[-1] <= 0:
            return False, "eq_lc_fail"
        r0 = q[:]
    elif dq == dp + 1:
        r0 = q[:]
        r1 = p[:]
    else:
        return False, f"deg_jump({dp}->{dq})"
    cur0, cur1 = r0[:], r1[:]
    while len(cur1) > 1:
        d0, d1 = len(cur0) - 1, len(cur1) - 1
        if d0 != d1 + 1:
            return False, f"deg_drop({d0},{d1})"
        a0, a1 = cur0[-1], cur1[-1]
        tmp = [a1 * cur0[m] for m in range(d0 + 1)]
        for m in range(d1 + 1):
            tmp[m + 1] -= a0 * cur1[m]
        b = tmp[d1]
        nxt = [-(a1 * tmp[m] - b * cur1[m]) for m in range(d1)]
        while len(nxt) > 1 and nxt[-1] == 0:
            nxt.pop()
        if len(nxt) - 1 != d1 - 1 or nxt[-1] <= 0:
            return False, f"sturm_fail_deg{d1-1}"
        g = math.gcd(*[abs(v) for v in nxt])
        if g > 1:
            nxt = [v // g for v in nxt]
        cur0, cur1 = cur1, nxt
    return True, "ok"

# === TEST 1: Exact Sturm check of Q_k(0, y) vs d_t Q_k(0, y) for k up to 30 at q=1/2, 1/3, 1/4! ===
print("=== TEST 1: Necessary stability test at t=0: Q_k(0, y) vs d_t Q_k(0, y) for k up to 30 ===")
# We only need the first TWO coefficients [t^0] and [t^1] of htilde_{k, i}(t)!
# Recurrence for (c0, c1) of htilde_{k, i}(t):
# D_j(c0 + c1*t + ...) = c0 + (2*c1 + j*c0)*t + ...
c0 = {(1, 0): 1, (1, 1): 1}
c1 = {(1, 0): 2, (1, 1): 0}
for k in range(2, 31):
    j = k + 1
    for i in range(k + 1):
        # D_j(ht_{k-1, i})
        a0 = c0[(k-1, i)] if i <= k - 1 else 0
        a1 = (2 * c1[(k-1, i)] + j * c0[(k-1, i)]) if i <= k - 1 else 0
        # + ht_{k-1, i-1}
        s0 = a0 + (c0[(k-1, i-1)] if i >= 1 else 0)
        s1 = a1 + (c1[(k-1, i-1)] if i >= 1 else 0)
        # D_j(S)
        c0[(k, i)] = s0
        c1[(k, i)] = 2 * s1 + j * s0

for q_num, q_den in [(1, 1), (3, 4), (2, 3), (1, 2), (1, 3), (1, 4)]:
    first_fail = None
    for k in range(2, 31):
        E = k * (k - 1) // 2
        P0 = [c0[(k, i)] * (q_num ** (i*(i-1)//2)) * (q_den ** (E - i*(i-1)//2)) for i in range(k + 1)]
        P1 = [c1[(k, i)] * (q_num ** (i*(i-1)//2)) * (q_den ** (E - i*(i-1)//2)) for i in range(k + 1)]
        ok, msg = exact_strict_interlaces(P1, P0)
        if not ok:
            first_fail = (k, msg)
            break
    print(f"  q={q_num}/{q_den}: P1<<P0 up to k=30 -> {'PASS ALL!' if first_fail is None else f'FAIL at k={first_fail[0]} ({first_fail[1]})'}")

# Also check AT t = -1! (i.e. x -> -infinity, leading coefficients of A_{k, i}(x) and second-leading!)
# Remember A_{k, i}(x) has deg 2k-i-1. Look at lc(A_{k,i}) and next coeff!
# In t: ht_{k, i}(-1) = (-1)^{2k-i-1} * lc(A_{k, i})!
# Wait: easier to compute A_{k, i} top two coeffs directly!
# S_L(a_d x^d + a_{d-1} x^{d-1} + ...) = (L-d) a_d x^{d+1} + ((d+1) a_d + (L-d+1) a_{d-1}) x^d + ...
print("\n=== TEST 1b: Necessary stability test at x -> -infinity (top two coeffs of A_{k,i}(x)) ===")
top0 = {(1, 0): 1, (1, 1): 1} # leading coeff [x^d] A_{k,i}
top1 = {(1, 0): 1, (1, 1): 0} # sub-leading coeff [x^{d-1}] A_{k,i}
for k in range(2, 31):
    for i in range(k + 1):
        # S_{3k-3-i} A_{k-1, i}: deg of A_{k-1, i} is d1 = 2k-i-3, L1 = 3k-3-i -> L1 - d1 = k!
        if i <= k - 1:
            d1 = 2*k - i - 3
            L1 = 3*k - 3 - i
            a0 = (L1 - d1) * top0[(k-1, i)]
            a1 = (d1 + 1) * top0[(k-1, i)] + (L1 - d1 + 1) * top1[(k-1, i)]
        else:
            a0, a1 = 0, 0
        # A_{k-1, i-1} has deg 2k-i-2 = d1 + 1!
        s0 = a0 + (top0[(k-1, i-1)] if i >= 1 else 0)
        s1 = a1 + (top1[(k-1, i-1)] if i >= 1 else 0)
        # Apply S_{3k-2-i} to S (deg ds = 2k-i-2, L2 = 3k-2-i -> L2 - ds = k!)
        ds = 2*k - i - 2
        L2 = 3*k - 2 - i
        top0[(k, i)] = (L2 - ds) * s0
        top1[(k, i)] = (ds + 1) * s0 + (L2 - ds + 1) * s1

print("  Check top0[(k, i)] formula for k=1..5:")
for k in range(1, 6):
    print(f"    k={k}: {[top0[(k, i)] for i in range(k+1)]}")
for q_num, q_den in [(1, 1), (3, 4), (2, 3), (1, 2), (1, 3), (1, 4)]:
    first_fail = None
    for k in range(2, 31):
        E = k * (k - 1) // 2
        # For A(x) = sum top0_i x^{2k-1-i} + top1_i x^{2k-2-i}, set z = 1/x, w = y/x -> sum top0_i w^i + z * sum top1_i w^i
        # Wait: in (x, y) near x = -inf, keep y fixed or scale y = x * w!
        # Let's check interlacing of top0 and top1 in w!
        P0 = [top0[(k, i)] * (q_num ** (i*(i-1)//2)) * (q_den ** (E - i*(i-1)//2)) for i in range(k + 1)]
        P1 = [top1[(k, i)] * (q_num ** (i*(i-1)//2)) * (q_den ** (E - i*(i-1)//2)) for i in range(k + 1)]
        ok1, m1 = exact_strict_interlaces(P1, P0)
        ok2, m2 = exact_strict_interlaces(P0, P1)
        if not (ok1 or ok2):
            first_fail = (k, m1, m2)
            break
    print(f"  q={q_num}/{q_den}: top1 vs top0 interlacing up to k=30 -> {'PASS ALL!' if first_fail is None else f'FAIL at k={first_fail[0]} ({first_fail[1]}, {first_fail[2]})'}")

# === TEST 2: Test M_{b, q} on PROPERLY q-transformed F(t, y) = T_q(prod (a_r + b_r t + c_r y))! ===
print("\n=== TEST 2: One step of M_{b, q} on Q = T_q(prod (a_r + b_r t + c_r y)), 0 < a_r < b_r, c_r > 0 ===")
def poly_mul_np(a, b):
    out = np.zeros(len(a) + len(b) - 1, dtype=float)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out
def apply_D_np(b, p):
    d = len(p) - 1
    out = np.zeros(d + 2, dtype=float)
    for m, cm in enumerate(p):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + b) * cm
    return out
def add_np(a, b):
    n = max(len(a), len(b))
    out = np.zeros(n, dtype=float)
    out[:len(a)] += a
    out[:len(b)] += b
    return out
def eval_c(p, z):
    acc = 0j
    for c in reversed(p): acc = acc * z + c
    return acc

test_ts = [complex(re, im) for re in [-10, -2, -1.05, -0.95, -0.5, -0.1, 0, 0.1, 1, 10] for im in [1e-4, 0.02, 0.5]]
random.seed(99)
for q in [0.5, 0.25, 0.1]:
    fails_input = 0
    fails_output = 0
    for trial in range(60):
        m = random.choice([2, 3, 4])
        F = [np.array([1.0])]
        for _ in range(m):
            ar = random.uniform(0.1, 5.0)
            br = ar + random.uniform(0.1, 20.0)
            cr = 10 ** random.uniform(-1, 1)
            newF = [np.zeros(1) for _ in range(len(F) + 1)]
            lin = np.array([ar, br])
            for i in range(len(F)):
                newF[i] = add_np(newF[i], poly_mul_np(F[i], lin))
                newF[i + 1] = add_np(newF[i + 1], F[i] * cr)
            F = newF
        # Apply T_q to F!
        Q = [F[i] * (q ** (i*(i-1)//2)) for i in range(len(F))]
        # Check input Q is stable (sanity check!)
        # Apply M_{b, q}: G[i] = D_b Q[i] + q^{i-1} Q[i-1]
        b = random.choice([2.0, 5.0, 50.0])
        k = len(Q) - 1
        G = [add_np(apply_D_np(b, Q[i]) if i <= k else np.array([0.0]),
                    (q ** (i - 1)) * Q[i - 1] if i >= 1 else np.array([0.0]))
             for i in range(k + 2)]
        max_im = -1e9
        for tv in test_ts:
            yc = [eval_c(G[i], tv) for i in range(k + 1, -1, -1)]
            mx = max(abs(c) for c in yc)
            if mx == 0: continue
            yc = [c / mx for c in yc]
            mi = max(z.imag for z in np.roots(yc))
            if mi > max_im: max_im = mi
        if max_im > 1e-6:
            fails_output += 1
    print(f"  q={q:4.2f} on T_q(prod linear): fails={fails_output:2d}/60")

# === TEST 3: Check 3x3 minors of coefficient matrix [t^m] ht_{k,i} and [x^m] A_{k,i}! ===
print("\n=== TEST 3: Check 3x3 minors of coefficient matrices of htilde and A (k=5, 6) ===")
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 7):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

import sympy as sp
for k in [5, 6]:
    row = [ht[(k, i)] for i in range(k, -1, -1)] # increasing degree
    R = len(row)
    C = max(len(p) for p in row)
    M = [[row[r][c] if c < len(row[r]) else 0 for c in range(C)] for r in range(R)]
    neg3 = 0
    for r1 in range(R):
        for r2 in range(r1+1, R):
            for r3 in range(r2+1, R):
                for c1 in range(C):
                    for c2 in range(c1+1, C):
                        for c3 in range(c2+1, C):
                            d = (M[r1][c1]*(M[r2][c2]*M[r3][c3]-M[r2][c3]*M[r3][c2])
                               - M[r1][c2]*(M[r2][c1]*M[r3][c3]-M[r2][c3]*M[r3][c1])
                               + M[r1][c3]*(M[r2][c1]*M[r3][c2]-M[r2][c2]*M[r3][c1]))
                            if d < 0:
                                neg3 += 1
    print(f"  htilde k={k}: neg 3x3 minors in coeff matrix = {neg3}")
