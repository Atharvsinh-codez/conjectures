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

def apply_gen(a, b, c_op, p):
    # T(p) = t(1+t)p' + (a + b*t)*p   (here c_op unused)
    d = len(p) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(p):
        out[m] += (m + a) * cm
        out[m + 1] += (m + b) * cm
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

print("=== 1. Test prod_{j=1}^k (D_{b_j} + y)(1) WITHOUT outer D_{a_j}! ===")
def test_no_outer(name, b_list):
    row = [[1]]
    ok_rr = True
    ok_int = True
    fail_rr = None
    fail_int = None
    for step, bj in enumerate(b_list, 1):
        new_row = []
        for i in range(step + 1):
            t1 = apply_D(bj, row[i]) if i <= step - 1 else [0]
            t2 = row[i - 1] if i >= 1 else [0]
            new_row.append(add_poly(t1, t2))
        for i in range(step + 1):
            if not is_real_rooted(new_row[i]):
                ok_rr = False
                if fail_rr is None: fail_rr = (step, i)
        for i in range(1, step + 1):
            ok, msg = exact_strict_interlaces(new_row[i], new_row[i-1])
            if not ok:
                ok_int = False
                if fail_int is None: fail_int = (step, i, msg)
        row = new_row
    print(f"  {name:30s}: RR={str(ok_rr):5s} (fail={fail_rr}) | q_i<<q_{{i-1}}={str(ok_int):5s} (fail={fail_int})")

test_no_outer("GLZ b=[2..10]", [2,3,4,5,6,7,8,9,10])
test_no_outer("Const b=2 x8", [2]*8)
test_no_outer("Wild [2,100,2,100,2,100]", [2,100,2,100,2,100])
test_no_outer("Wild [100,2,100,2,100,2]", [100,2,100,2,100,2])
test_no_outer("Extreme [2,10000,2,10000,2]", [2,10000,2,10000,2])
test_no_outer("Extreme [10000,2,10000,2,10000]", [10000,2,10000,2,10000])
test_no_outer("Drop [1000,100,10,2,2,2]", [1000,100,10,2,2,2])
test_no_outer("Rise [2,2,2,10,100,1000]", [2,2,2,10,100,1000])

print("\n=== 2. Test prod_{j=1}^k (T_j + y)(1) for OTHER differential operators T_j! ===")
# T_j(p) = t(1+t)p' + (a_j + b_j t) p, with a_j > 0, b_j > a_j
def test_gen_op(name, ab_list):
    row = [[1]]
    ok_rr = True
    ok_int = True
    fail_int = None
    for step, (aj, bj) in enumerate(ab_list, 1):
        new_row = []
        for i in range(step + 1):
            t1 = apply_gen(aj, bj, 0, row[i]) if i <= step - 1 else [0]
            t2 = row[i - 1] if i >= 1 else [0]
            new_row.append(add_poly(t1, t2))
        for i in range(step + 1):
            if not is_real_rooted(new_row[i]):
                ok_rr = False
        for i in range(1, step + 1):
            ok, msg = exact_strict_interlaces(new_row[i], new_row[i-1])
            if not ok:
                ok_int = False
                if fail_int is None: fail_int = (step, i, msg)
        row = new_row
    print(f"  {name:35s}: RR={str(ok_rr):5s} | q_i<<q_{{i-1}}={str(ok_int):5s} (fail={fail_int})")

test_gen_op("a=1, b=j+1 (GLZ inner)", [(1, j+1) for j in range(1, 8)])
test_gen_op("a=j, b=2j", [(j, 2*j) for j in range(1, 8)])
test_gen_op("a=1,100,1,100; b=2,200,2,200", [(1,2),(100,200),(1,2),(100,200),(1,2),(100,200)])
test_gen_op("a=100,1,100,1; b=101,100,101,100", [(100,101),(1,100),(100,101),(1,100),(100,101)])
test_gen_op("a=1,100,1,100; b=1000,101,1000,101", [(1,1000),(100,101),(1,1000),(100,101),(1,1000)])

print("\n=== 3. Test prod (L_j + y)(1) for SIMPLE operators: L_j(p) = (t + c_j)*p + d_j*p' ===")
def apply_L(cj, dj, p):
    # (t + cj)*p + dj*p'
    d = len(p) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(p):
        out[m] += cj * cm
        out[m + 1] += cm
        if m >= 1:
            out[m - 1] += dj * m * cm
    return poly_trim(out)

def test_L(name, cd_list):
    row = [[1]]
    ok_rr = True
    ok_int = True
    fail_int = None
    for step, (cj, dj) in enumerate(cd_list, 1):
        new_row = []
        for i in range(step + 1):
            t1 = apply_L(cj, dj, row[i]) if i <= step - 1 else [0]
            t2 = row[i - 1] if i >= 1 else [0]
            new_row.append(add_poly(t1, t2))
        for i in range(step + 1):
            if not is_real_rooted(new_row[i]):
                ok_rr = False
        for i in range(1, step + 1):
            ok, msg = exact_strict_interlaces(new_row[i], new_row[i-1])
            if not ok:
                ok_int = False
                if fail_int is None: fail_int = (step, i, msg)
        row = new_row
    print(f"  {name:35s}: RR={str(ok_rr):5s} | q_i<<q_{{i-1}}={str(ok_int):5s} (fail={fail_int})")

test_L("c=0, d=-1 (Hermite!)", [(0, -1)]*7)
test_L("c=j, d=-1", [(j, -1) for j in range(1, 8)])
test_L("c=[1,100,1,100], d=-1", [(1,-1),(100,-1),(1,-1),(100,-1),(1,-1),(100,-1)])
test_L("c=[0,0,0,0,0], d=[-1,-100,-1,-100,-1]", [(0,-1),(0,-100),(0,-1),(0,-100),(0,-1)])
test_L("c=[1,100,1,100], d=[-1,-100,-1,-100]", [(1,-1),(100,-100),(1,-1),(100,-100),(1,-1)])
test_L("c=[1,100,1,100], d=[-100,-1,-100,-1]", [(1,-100),(100,-1),(1,-100),(100,-1),(1,-100)])

print("\n=== 4. Test F_{lam,mu} << G_{lam,mu} (i.e. real-rootedness of a*l*p_i + (a*m+b*l)*p_{i-1} + b*m*p_{i-2}) on GLZ htilde ===")
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 9):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

# Test if (lam*p_i + mu*p_{i-1}) << (lam*p_{i-1} + mu*p_{i-2}) for all real lam, mu!
# Note: test exact values of lam/mu in {-100,-10,-3,-1,-0.5,-0.1, 0.1, 1, 10}
ratios = [(-100,1),(-10,1),(-3,1),(-2,1),(-1,1),(-1,2),(-1,3),(-1,10),(-1,100),(1,10),(1,1),(10,1)]
for k in range(2, 8):
    fails = []
    for i in range(2, k + 1):
        pi, pi1, pi2 = ht[(k, i)], ht[(k, i-1)], ht[(k, i-2)]
        for (la, mu) in ratios:
            F = add_poly(scale_poly(pi, la), scale_poly(pi1, mu))
            G = add_poly(scale_poly(pi1, la), scale_poly(pi2, mu))
            # Check if F << G (handle negative lc if mu < 0)
            if mu < 0:
                F = scale_poly(F, -1)
                G = scale_poly(G, -1)
            ok, msg = exact_strict_interlaces(F, G)
            if not ok:
                fails.append((i, la, mu, msg))
                break
    print(f"  k={k}: F_{{la,mu}} << G_{{la,mu}}: {'PASS ALL' if not fails else f'FAIL {fails[:3]}'}")

print("\n=== 5. Test if (htilde_{k,0},...,htilde_{k,k}) is a REAL-ROOTED BASIS (i.e. sum c_i htilde_{k,i} real-rooted whenever sum c_i y^i has all real/neg roots) ===")
# Specifically test if Turan / quadratic c_0 p_{i-2} + c_1 p_{i-1} + c_2 p_i with c_1^2 >= 4 c_0 c_2 is real-rooted!
# Wait: that's IDENTICAL to Test 4! Check if it passes when la, mu have SAME sign vs OPPOSITE sign!
for k in range(2, 8):
    fails_pos = []
    for i in range(2, k + 1):
        pi, pi1, pi2 = ht[(k, i)], ht[(k, i-1)], ht[(k, i-2)]
        for (la, mu) in [(1,100),(1,10),(1,3),(1,1),(3,1),(10,1),(100,1)]:
            F = add_poly(scale_poly(pi, la), scale_poly(pi1, mu))
            G = add_poly(scale_poly(pi1, la), scale_poly(pi2, mu))
            ok, msg = exact_strict_interlaces(F, G)
            if not ok:
                fails_pos.append((i, la, mu, msg))
                break
    print(f"  k={k}: F_{{la,mu}} << G_{{la,mu}} (la,mu > 0): {'PASS ALL' if not fails_pos else f'FAIL {fails_pos[:3]}'}")
