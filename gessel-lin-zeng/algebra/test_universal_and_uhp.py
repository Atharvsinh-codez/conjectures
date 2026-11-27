import math
import random
import numpy as np

def poly_trim(c):
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c

def apply_D(j, c):
    # D_j(p) = t(1+t)p' + (1 + j*t)p
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return poly_trim(out)

def apply_S(L, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (L - m) * cm
    return poly_trim(out)

def add_poly(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)):
        out[i] += a[i]
    for i in range(len(b)):
        out[i] += b[i]
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

def is_real_rooted(p):
    d = len(p) - 1
    if d <= 1:
        return True
    if p[-1] < 0:
        p = [-v for v in p]
    dp = [m * p[m] for m in range(1, d + 1)]
    return exact_strict_interlaces(dp, p)[0]

print("=== PART A: Test arbitrary sequences c_1,...,c_k >= 2 in prod D_{c_j}(D_{c_j} + y)(1) ===")
def test_seq(name, c_list):
    row = [[1]]
    ok_rr = True
    ok_h = True
    ok_S = True
    ok_BA = True
    ok_AB = True
    ok_AD = True
    first_fail_h = None
    first_fail_S = None
    for step, cj in enumerate(c_list, 1):
        k = step
        S_row = []
        new_row = []
        for i in range(k + 1):
            t1 = apply_D(cj, row[i]) if i <= k - 1 else [0]
            t2 = row[i - 1] if i >= 1 else [0]
            S_i = add_poly(t1, t2)
            S_row.append(S_i)
            new_row.append(apply_D(cj, S_i))
        for i in range(k + 1):
            if not is_real_rooted(new_row[i]):
                ok_rr = False
        for i in range(1, k + 1):
            ok, msg = exact_strict_interlaces(new_row[i], new_row[i-1])
            if not ok:
                ok_h = False
                if first_fail_h is None: first_fail_h = (k, i, msg)
            ok2, msg2 = exact_strict_interlaces(S_row[i], S_row[i-1])
            if not ok2:
                ok_S = False
                if first_fail_S is None: first_fail_S = (k, i, msg2)
        if k >= 2:
            for i in range(1, k):
                if not exact_strict_interlaces(row[i-1], apply_D(cj, row[i]))[0]:
                    ok_BA = False
                if not exact_strict_interlaces(apply_D(cj, row[i]), row[i-1])[0]:
                    ok_AB = False
            for i in range(2, k):
                if not exact_strict_interlaces(apply_D(cj, row[i]), row[i-2])[0]:
                    ok_AD = False
        row = new_row
    print(f"  {name:26s}: RR={str(ok_rr):5s} h_i<<h_{{i-1}}={str(ok_h):5s} S_i<<S_{{i-1}}={str(ok_S):5s} | B<<A={str(ok_BA):5s} A<<B={str(ok_AB):5s} A<<D={str(ok_AD):5s} | fail_h={first_fail_h}")

test_seq("GLZ shifted: [2..9]", [2,3,4,5,6,7,8,9])
test_seq("Const 2: [2]*8", [2]*8)
test_seq("Const 10: [10]*8", [10]*8)
test_seq("Wild: [2,100,2,100,2,100]", [2,100,2,100,2,100])
test_seq("Wild2: [100,2,100,2,100,2]", [100,2,100,2,100,2])
test_seq("Steep: [2,4,8,16,32,64,128]", [2,4,8,16,32,64,128])
test_seq("Drop: [100,50,25,12,6,3,2]", [100,50,25,12,6,3,2])
test_seq("Spike: [2,2,2,1000,2,2]", [2,2,2,1000,2,2])
test_seq("Spike2: [1000,2,2,2,2,2]", [1000,2,2,2,2,2])
test_seq("Two-op? [2,3,4,5,6,7,8]", [2,3,4,5,6,7,8])

# Also test TWO DIFFERENT operators D_a(D_b + y) at each step!
print("\n=== PART B: Test D_{a_j}(D_{b_j} + y) with a_j != b_j ===")
def test_two_op(name, ab_list):
    row = [[1]]
    ok_rr = True
    ok_h = True
    first_fail = None
    for step, (aj, bj) in enumerate(ab_list, 1):
        k = step
        new_row = []
        for i in range(k + 1):
            t1 = apply_D(bj, row[i]) if i <= k - 1 else [0]
            t2 = row[i - 1] if i >= 1 else [0]
            new_row.append(apply_D(aj, add_poly(t1, t2)))
        for i in range(k + 1):
            if not is_real_rooted(new_row[i]):
                ok_rr = False
        for i in range(1, k + 1):
            ok, msg = exact_strict_interlaces(new_row[i], new_row[i-1])
            if not ok:
                ok_h = False
                if first_fail is None: first_fail = (k, i, msg)
        row = new_row
    print(f"  {name:30s}: RR={str(ok_rr):5s} h_i<<h_{{i-1}}={str(ok_h):5s} fail={first_fail}")

test_two_op("a=2, b=100 x6", [(2,100)]*6)
test_two_op("a=100, b=2 x6", [(100,2)]*6)
test_two_op("a=2, b=1000 x6", [(2,1000)]*6)
test_two_op("a=1000, b=2 x6", [(1000,2)]*6)
test_two_op("a=j+1, b=j+1 x6 (GLZ)", [(j+1,j+1) for j in range(1,7)])
test_two_op("a=10j, b=j+1 x6", [(10*j,j+1) for j in range(1,7)])
test_two_op("a=j+1, b=10j x6", [(j+1,10*j) for j in range(1,7)])

# === PART C: Adversarial complex upper-half-plane stability test for B_k(x, y) ===
print("\n=== PART C: Adversarial UHP complex stability test for B_k(x, y) = sum binom(k,i) A_{k,i}(x) y^i ===")
A = {(0, 0): [1]}
for k in range(1, 13):
    for i in range(k + 1):
        t1 = apply_S(3*k - 3 - i, A[(k-1, i)]) if i <= k - 1 else [0]
        t2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S(3*k - 2 - i, add_poly(t1, t2))

def eval_poly_c(coeffs, z):
    # Horner
    acc = 0j
    for c in reversed(coeffs):
        acc = acc * z + c
    return acc

random.seed(42)
# Generate adversarial x in UHP: near negative real axis, near 0, positive real axis, large modulus
xs = []
for re in [-100.0, -10.0, -3.0, -1.5, -1.0, -0.5, -0.1, -0.01, 0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0]:
    for im in [1e-5, 1e-3, 0.01, 0.1, 0.5, 2.0, 10.0]:
        xs.append(complex(re, im))
for _ in range(300):
    re = random.uniform(-20, 20)
    im = 10 ** random.uniform(-4, 1.5)
    xs.append(complex(re, im))

for k in range(2, 13):
    w = [math.comb(k, i) for i in range(k + 1)]
    max_im_y = -1e9
    worst_x = None
    for xv in xs:
        # y-polynomial coefficients [c_0, c_1, ..., c_k]
        # np.roots expects highest degree first
        ycoeffs = [w[i] * eval_poly_c(A[(k, i)], xv) for i in range(k, -1, -1)]
        if abs(ycoeffs[0]) == 0:
            continue
        roots_y = np.roots(ycoeffs)
        mi = max(r.imag for r in roots_y)
        if mi > max_im_y:
            max_im_y = mi
            worst_x = xv
    status = "STABLE (no UHP root)" if max_im_y <= 1e-9 else f"UNSTABLE! max Im(y)={max_im_y:.6e} at x={worst_x}"
    print(f"  k={k:2d}: B_k(x,y) -> {status}")

# Also test binom-weighted htilde_{k,i}(t) in UHP!
print("\n=== PART D: Adversarial UHP complex stability test for sum binom(k,i) htilde_{k,i}(t) y^i ===")
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 13):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

ts = []
for re in [-5.0, -1.5, -1.01, -1.0, -0.99, -0.8, -0.5, -0.2, -0.01, 0.0, 0.01, 0.5, 2.0, 10.0]:
    for im in [1e-5, 1e-3, 0.01, 0.1, 0.5, 2.0, 10.0]:
        ts.append(complex(re, im))
for _ in range(300):
    re = random.uniform(-3, 3)
    im = 10 ** random.uniform(-4, 1.5)
    ts.append(complex(re, im))

for k in range(2, 13):
    w = [math.comb(k, i) for i in range(k + 1)]
    max_im_y = -1e9
    worst_t = None
    for tv in ts:
        ycoeffs = [w[i] * eval_poly_c(ht[(k, i)], tv) for i in range(k, -1, -1)]
        if abs(ycoeffs[0]) == 0:
            continue
        roots_y = np.roots(ycoeffs)
        mi = max(r.imag for r in roots_y)
        if mi > max_im_y:
            max_im_y = mi
            worst_t = tv
    status = "STABLE (no UHP root)" if max_im_y <= 1e-9 else f"UNSTABLE! max Im(y)={max_im_y:.6e} at t={worst_t}"
    print(f"  k={k:2d}: sum binom htilde -> {status}")
