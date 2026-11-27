import math
import random
import numpy as np
import sympy as sp

# --- Exact integer polynomial helpers (coeffs low->high) ---
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

def add_poly(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)):
        out[i] += a[i]
    for i in range(len(b)):
        out[i] += b[i]
    return poly_trim(out)

def scale_poly(c, s):
    return poly_trim([v * s for v in c])

def mul_poly(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
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

def is_real_rooted_sturm(p):
    """Exact Sturm test for all simple real roots."""
    d = len(p) - 1
    if d <= 1:
        return True
    if p[-1] < 0:
        p = [-v for v in p]
    dp = [m * p[m] for m in range(1, d + 1)]
    ok, _ = exact_strict_interlaces(dp, p)
    return ok

# Compute A_{k,i}(x) and htilde_{k,i}(t) for k=1..9
def apply_S(L, c):
    # S_L(f) = (1 + Lx)f + x(1-x)f'
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (L - m) * cm
    return poly_trim(out)

A = {(0, 0): [1]}
for k in range(1, 10):
    for i in range(k + 1):
        t1 = apply_S(3*k - 3 - i, A[(k-1, i)]) if i <= k - 1 else [0]
        t2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S(3*k - 2 - i, add_poly(t1, t2))

# === TEST 1: Real stability of B_k(x, y) = sum_{i=0}^k binom(k, i) A_{k,i}(x) y^i ===
# Recall: f(x, y) is real stable iff for all a > 0, b in R, f(x, a*x + b) is real-rooted in x!
print("=== TEST 1: Real stability of B_k(x, y) = sum binom(k,i) A_{k,i}(x) y^i ===")
def eval_line(k, weights, a_num, a_den, b_num, b_den, poly_dict):
    """Compute (a_den*b_den)^k * sum_i w_i * P_{k,i}(x) * (a*x + b)^i as exact int poly in x."""
    # Scale y = (a_num*b_den * x + b_num*a_den) / (a_den*b_den)
    M = a_den * b_den
    lin = [b_num * a_den, a_num * b_den]
    tot = [0]
    lin_pow = [1]
    for i in range(k + 1):
        w = weights[i] * (M ** (k - i))
        term = scale_poly(mul_poly(poly_dict[(k, i)], lin_pow), w)
        tot = add_poly(tot, term)
        lin_pow = mul_poly(lin_pow, lin)
    return tot

test_lines = [
    (1, 1, 0, 1), (1, 1, 1, 1), (1, 1, -1, 1), (1, 1, 5, 1), (1, 1, -5, 1),
    (1, 1, 100, 1), (1, 1, -100, 1), (1, 10, 1, 1), (1, 10, -1, 1), (10, 1, 1, 1),
    (10, 1, -1, 1), (1, 100, -1, 2), (100, 1, -1, 2), (3, 7, -11, 5), (1, 1000, -1, 1),
    (1000, 1, 1, 1), (1000, 1, -1, 1), (1, 1, -1, 2), (1, 1, 1, 2), (2, 1, -3, 1),
    (1, 5, -2, 1), (5, 1, -2, 1), (1, 20, -1, 3), (20, 1, -1, 3), (1, 1, -20, 1),
]

for k in range(2, 9):
    w_binom = [math.comb(k, i) for i in range(k + 1)]
    failed = None
    for (an, ad, bn, bd) in test_lines:
        p = eval_line(k, w_binom, an, ad, bn, bd, A)
        if not is_real_rooted_sturm(p):
            failed = (an, ad, bn, bd)
            break
    print(f"  k={k}: B_k(x, y) binom-weighted A_{{k,i}}: {'PASS all 25 lines' if failed is None else f'FAIL at a={failed[0]}/{failed[1]}, b={failed[2]}/{failed[3]}'}")

# Also test binom-weighted htilde_{k,i}(t)!
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 9):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

print("\n=== TEST 1b: Real stability of sum binom(k,i) htilde_{k,i}(t) y^i ===")
for k in range(2, 9):
    w_binom = [math.comb(k, i) for i in range(k + 1)]
    failed = None
    for (an, ad, bn, bd) in test_lines:
        p = eval_line(k, w_binom, an, ad, bn, bd, ht)
        if not is_real_rooted_sturm(p):
            failed = (an, ad, bn, bd)
            break
    print(f"  k={k}: binom-weighted htilde_{{k,i}}: {'PASS all 25 lines' if failed is None else f'FAIL at a={failed[0]}/{failed[1]}, b={failed[2]}/{failed[3]}'}")

# Also test pure y-direction: at fixed x or t, is B_k(x, y) real-rooted in y?
print("\n=== TEST 1c: B_k(x, y) in y for fixed real x ===")
for k in range(2, 9):
    w_binom = [math.comb(k, i) for i in range(k + 1)]
    fail_x = None
    for xval in [-100, -10, -2, -1, 0, 1, 2, 10, 100]:
        # Evaluate A_{k,i}(xval)
        y_poly = []
        for i in range(k + 1):
            val = sum(c * (xval ** m) for m, c in enumerate(A[(k, i)]))
            y_poly.append(w_binom[i] * val)
        y_poly = poly_trim(y_poly)
        if not is_real_rooted_sturm(y_poly):
            fail_x = xval
            break
    print(f"  k={k}: B_k(x, y) real-rooted in y for x in test set: {'PASS' if fail_x is None else f'FAIL at x={fail_x}'}")

# === TEST 2: Generalized sequences c_1, ..., c_k in H_k(t, y) = prod D_{c_j}(D_{c_j} + y)(1) ===
print("\n=== TEST 2: Generalized operator sequences D_{c_j}(D_{c_j} + y) ===")
def test_seq(name, c_list):
    # Start from h_{0,0} = [1]
    row = [[1]]
    ok_g = True
    ok_AD = True
    ok_BA = True
    for step, cj in enumerate(c_list, 1):
        # S_i = D_{cj} row[i] + row[i-1]
        k = step
        S_row = []
        new_row = []
        for i in range(k + 1):
            t1 = apply_D(cj, row[i]) if i <= k - 1 else [0]
            t2 = row[i - 1] if i >= 1 else [0]
            S_i = add_poly(t1, t2)
            S_row.append(S_i)
            new_row.append(apply_D(cj, S_i))
        # Factor out (1+t) from every entry? Note D_{cj}(1) = 1 + cj*t, not (1+t) unless cj=1!
        # Check interlacing new_row[i] << new_row[i-1] directly (all roots in (-1,0)?)
        for i in range(1, k + 1):
            ok, _ = exact_strict_interlaces(new_row[i], new_row[i-1])
            if not ok:
                ok_g = False
        if k >= 2:
            for i in range(1, k):
                ok, _ = exact_strict_interlaces(row[i-1], apply_D(cj, row[i]))
                if not ok:
                    ok_BA = False
            for i in range(2, k):
                ok, _ = exact_strict_interlaces(apply_D(cj, row[i]), row[i-2])
                if not ok:
                    ok_AD = False
        row = new_row
    print(f"  {name:30s}: h_i<<h_{{i-1}}={ok_g}, B<<A={ok_BA}, A<<D={ok_AD}")

test_seq("GLZ exact: [1,2,3,4,5,6,7,8]", [1,2,3,4,5,6,7,8])
test_seq("Constant c=1: [1]*8", [1]*8)
test_seq("Constant c=2: [2]*8", [2]*8)
test_seq("Constant c=5: [5]*8", [5]*8)
test_seq("Double: [2,4,6,8,10,12,14]", [2,4,6,8,10,12,14])
test_seq("Decreasing: [10,9,8,7,6,5,4]", [10,9,8,7,6,5,4])
test_seq("Wild: [1,10,2,20,3,30,4]", [1,10,2,20,3,30,4])
test_seq("Steep: [1,2,4,8,16,32,64]", [1,2,4,8,16,32,64])
test_seq("Slow: [1,1,1,2,2,2,3,3]", [1,1,1,2,2,2,3,3])
test_seq("Two-op: D_{a_j}(D_{b_j}+y)?", [3,3,3,3,3,3,3])
