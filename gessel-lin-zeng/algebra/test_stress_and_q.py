import math
import random
import numpy as np

def poly_trim(c):
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c

def apply_gen(a, b, p):
    # T(p) = t(1+t)p' + (a + b*t)*p
    d = len(p) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(p):
        out[m] += (m + a) * cm
        out[m + 1] += (m + b) * cm
    return poly_trim(out)

def apply_D(j, c):
    return apply_gen(1, j, c)

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

def scale_poly(c, s):
    return poly_trim([v * s for v in c])

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

# === PART A: Adversarially try to BREAK prod (T_j + y)(1) with wild integer a_j < b_j! ===
print("=== PART A: Adversarial stress test of prod_{j=1}^k (T_j + y)(1) with T_j = t(1+t)d/dt + (a_j + b_j t) ===")
random.seed(123)
# Test 200 extreme random sequences of length k=6!
fails_rr = 0
fails_int = 0
examples_fail = []
for trial in range(200):
    # Pick wildly varying a_j, b_j > a_j spanning orders of magnitude!
    ab = []
    for j in range(6):
        a = random.choice([1, 2, 10, 100, 10000])
        gap = random.choice([1, 2, 10, 100, 10000])
        ab.append((a, a + gap))
    row = [[1]]
    ok_rr = True
    ok_int = True
    for step, (aj, bj) in enumerate(ab, 1):
        new_row = []
        for i in range(step + 1):
            t1 = apply_gen(aj, bj, row[i]) if i <= step - 1 else [0]
            t2 = row[i - 1] if i >= 1 else [0]
            new_row.append(add_poly(t1, t2))
        for i in range(step + 1):
            if not is_real_rooted(new_row[i]):
                ok_rr = False
        for i in range(1, step + 1):
            if not exact_strict_interlaces(new_row[i], new_row[i-1])[0]:
                ok_int = False
        row = new_row
    if not ok_rr:
        fails_rr += 1
        if len(examples_fail) < 5:
            examples_fail.append(("RR", ab))
    elif not ok_int:
        fails_int += 1
        if len(examples_fail) < 5:
            examples_fail.append(("INT", ab))

print(f"  200 wild (a_j, b_j) trials (k=6): fails_RR={fails_rr}, fails_INT={fails_int}")
for kind, ab in examples_fail:
    print(f"    FAIL {kind}: ab={ab}")

# Also test extreme a_j = 1, b_j arbitrary (pure D_{b_j}) for 200 wild trials!
fails_rr_D = 0
fails_int_D = 0
ex_D = []
for trial in range(200):
    bs = [random.choice([2, 3, 5, 20, 100, 1000, 100000]) for _ in range(6)]
    row = [[1]]
    ok_rr = True
    ok_int = True
    for step, bj in enumerate(bs, 1):
        new_row = []
        for i in range(step + 1):
            t1 = apply_D(bj, row[i]) if i <= step - 1 else [0]
            t2 = row[i - 1] if i >= 1 else [0]
            new_row.append(add_poly(t1, t2))
        for i in range(step + 1):
            if not is_real_rooted(new_row[i]):
                ok_rr = False
        for i in range(1, step + 1):
            if not exact_strict_interlaces(new_row[i], new_row[i-1])[0]:
                ok_int = False
        row = new_row
    if not ok_rr:
        fails_rr_D += 1
        if len(ex_D) < 5: ex_D.append(("RR", bs))
    elif not ok_int:
        fails_int_D += 1
        if len(ex_D) < 5: ex_D.append(("INT", bs))
print(f"  200 wild D_{{b_j}} trials (k=6): fails_RR={fails_rr_D}, fails_INT={fails_int_D}")
for kind, bs in ex_D:
    print(f"    FAIL {kind}: bs={bs}")

# Also test: what if we put weights y_j * c_j, i.e. prod (D_{b_j} + w_j * y)(1) with wild positive w_j?!
fails_rr_W = 0
fails_int_W = 0
ex_W = []
for trial in range(200):
    bw = [(random.choice([2, 10, 1000]), random.choice([1, 10, 1000, 100000])) for _ in range(6)]
    row = [[1]]
    ok_rr = True
    ok_int = True
    for step, (bj, wj) in enumerate(bw, 1):
        new_row = []
        for i in range(step + 1):
            t1 = apply_D(bj, row[i]) if i <= step - 1 else [0]
            t2 = scale_poly(row[i - 1], wj) if i >= 1 else [0]
            new_row.append(add_poly(t1, t2))
        for i in range(step + 1):
            if not is_real_rooted(new_row[i]):
                ok_rr = False
        for i in range(1, step + 1):
            if not exact_strict_interlaces(new_row[i], new_row[i-1])[0]:
                ok_int = False
        row = new_row
    if not ok_rr:
        fails_rr_W += 1
        if len(ex_W) < 5: ex_W.append(("RR", bw))
    elif not ok_int:
        fails_int_W += 1
        if len(ex_W) < 5: ex_W.append(("INT", bw))
print(f"  200 wild (D_{{b_j}} + w_j y) trials (k=6): fails_RR={fails_rr_W}, fails_INT={fails_int_W}")
for kind, bw in ex_W:
    print(f"    FAIL {kind}: bw={bw}")

# === PART B: Check F_{la, mu} << G_{la, mu} at k = 13..20 for GLZ htilde ===
print("\n=== PART B: Check F_{la,mu} << G_{la,mu} for k=13..18 (specifically around la/mu = -1) ===")
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 19):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

for k in range(13, 19):
    # Check Q_{k,i}(0, y) vs Q'_{k,i}(0, y) interlacing first! (Instant necessary condition!)
    # Q_{k,i}(t, y) = ht_{k,i-2}(t) + 2 ht_{k,i-1}(t) y + ht_{k,i}(t) y^2
    fail_d0 = []
    for i in range(2, k + 1):
        p0 = [ht[(k, i-2)][0], 2 * ht[(k, i-1)][0], ht[(k, i)][0]]
        p1 = [ht[(k, i-2)][1], 2 * ht[(k, i-1)][1], ht[(k, i)][1] if len(ht[(k, i)]) > 1 else 0]
        ok1, _ = exact_strict_interlaces(p1, p0)
        ok2, _ = exact_strict_interlaces(p0, p1)
        if not (ok1 or ok2):
            fail_d0.append(i)
    # Also test direct F_{la,mu} << G_{la,mu} at la=-1, mu=1
    fail_fm = []
    for i in range(2, k + 1):
        F = add_poly(scale_poly(ht[(k, i)], -1), ht[(k, i-1)])
        G = add_poly(scale_poly(ht[(k, i-1)], -1), ht[(k, i-2)])
        if not exact_strict_interlaces(F, G)[0]:
            fail_fm.append(i)
    print(f"  k={k:2d}: Q(0,y) vs Q_t(0,y) fails at i={fail_d0[:5]} | F_{{-1,1}}<<G_{{-1,1}} fails at i={fail_fm[:5]}")

# === PART C: Test q-exponential weighted stability sum q^{i(i-1)/2} A_{k,i}(x) y^i ===
print("\n=== PART C: Test q-weighted sum q^{i(i-1)/2} A_{k,i}(x) y^i for real stability ===")
A = {(0, 0): [1]}
for k in range(1, 13):
    for i in range(k + 1):
        t1 = apply_S(3*k - 3 - i, A[(k-1, i)]) if i <= k - 1 else [0]
        t2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S(3*k - 2 - i, add_poly(t1, t2))

def eval_poly_c(coeffs, z):
    acc = 0j
    for c in reversed(coeffs):
        acc = acc * z + c
    return acc

xs = [complex(re, im) for re in [-10.0, -2.0, -0.5, -0.05, 0.0, 0.05, 0.5, 2.0, 10.0] for im in [1e-4, 0.01, 0.2, 2.0]]
for q_name, q_val in [("q=1/2", 0.5), ("q=1/4", 0.25), ("q=1/8", 0.125), ("q=1/16", 0.0625)]:
    bad_k = None
    for k in range(2, 13):
        w = [q_val ** (i * (i - 1) // 2) for i in range(k + 1)]
        max_im = -1e9
        for xv in xs:
            yc = [w[i] * eval_poly_c(A[(k, i)], xv) for i in range(k, -1, -1)]
            if abs(yc[0]) == 0: continue
            mi = max(r.imag for r in np.roots(yc))
            if mi > max_im: max_im = mi
        if max_im > 1e-8:
            bad_k = (k, max_im)
            break
    print(f"  {q_name:8s}: {'ALL PASS k=2..12!' if bad_k is None else f'FAIL at k={bad_k[0]} (max Im y={bad_k[1]:.3e})'}")
