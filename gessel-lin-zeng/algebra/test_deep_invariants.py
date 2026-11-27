import math
import numpy as np
import sympy as sp

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

def is_real_rooted(p):
    d = len(p) - 1
    if d <= 1:
        return True
    if p[-1] < 0:
        p = [-v for v in p]
    dp = [m * p[m] for m in range(1, d + 1)]
    return exact_strict_interlaces(dp, p)[0]

# Compute A_{k,i}(x) and htilde_{k,i}(t) up to k=14
A = {(0, 0): [1]}
for k in range(1, 15):
    for i in range(k + 1):
        t1 = apply_S(3*k - 3 - i, A[(k-1, i)]) if i <= k - 1 else [0]
        t2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S(3*k - 2 - i, add_poly(t1, t2))

# === 1. Exact check of B_k(0, y) vs d_x B_k(0, y) interlacing! ===
print("=== 1. Exact check of B_k(0, y) vs d_x B_k(0, y) interlacing (necessary for real stability!) ===")
for k in range(2, 14):
    P0 = [math.comb(k, i) * A[(k, i)][0] for i in range(k + 1)]
    P1 = [math.comb(k, i) * (A[(k, i)][1] if len(A[(k, i)]) > 1 else 0) for i in range(k + 1)]
    # Check P1 << P0 or P0 << P1!
    ok_10, m10 = exact_strict_interlaces(P1, P0)
    ok_01, m01 = exact_strict_interlaces(P0, P1)
    # Also check exact real-rootedness of B_k(x, a*x + b) at a=1, b=-20? Let's check exact line for k=10!
    print(f"  k={k:2d}: P1<<P0={ok_10} ({m10}) | P0<<P1={ok_01} ({m01}) | deg P0={len(P0)-1}, deg P1={len(P1)-1}")

# Find exact rational witness line for k=10 instability of B_{10}(x, y)!
# Since dy/dx > 0 at some root y_0 of B_{10}(0, y), near (0, y_0) the curve B_{10}(x, y)=0 has POSITIVE slope!
# So a line y = a*x + y_0 with small positive slope a (between 0 and dy/dx) or y_0 slightly shifted misses real roots!
print("\n=== 2. Exact rational witness proving B_{10}(x, y) is NOT real stable ===")
P0_10 = [math.comb(10, i) * A[(10, i)][0] for i in range(11)]
P1_10 = [math.comb(10, i) * A[(10, i)][1] for i in range(11)]
r0 = np.sort(np.roots(list(reversed(P0_10))).real)
for r in r0:
    # compute dy/dx = -P1(r) / P0'(r)
    v1 = sum(c * (r**i) for i, c in enumerate(P1_10))
    vp = sum(i * c * (r**(i-1)) for i, c in enumerate(P0_10) if i >= 1)
    slope = -v1 / vp
    if slope > 0:
        print(f"  POSITIVE SLOPE at root y0={r:.6f}: dy/dx = {slope:+.4f} > 0 !!")

# Check exact line B_{10}(x, a*x + b) near that root!
def eval_line_A(k, a_num, a_den, b_num, b_den):
    M = a_den * b_den
    lin = [b_num * a_den, a_num * b_den]
    tot = [0]
    lin_pow = [1]
    for i in range(k + 1):
        w = math.comb(k, i) * (M ** (k - i))
        term = scale_poly(mul_poly(A[(k, i)], lin_pow), w)
        tot = add_poly(tot, term)
        lin_pow = mul_poly(lin_pow, lin)
    return tot

for (an, ad, bn, bd) in [(1, 1, -20, 1), (10, 1, -15, 1), (100, 1, -15, 1), (1000, 1, -15, 1), (100, 1, -30, 1), (500, 1, -25, 1)]:
    p = eval_line_A(10, an, ad, bn, bd)
    rr = is_real_rooted(p)
    if not rr:
        print(f"  EXACT STURM COUNTEREXAMPLE to B_10(x,y) stability: y = ({an}/{ad})*x + ({bn}/{bd}) -> NOT real-rooted!")

# Check A_{k, i}(0) and A'_{k, i}(0) formula!
print("\n=== 3. Formula for A'_{k,i}(0) / binom(k,i) ===")
for k in range(1, 7):
    vals = [A[(k, i)][1] // math.comb(k, i) if len(A[(k, i)]) > 1 else 0 for i in range(k + 1)]
    print(f"  k={k}: A'_{{k,i}}(0)/binom(k,i) = {vals}")

# === 4. Check if H_k(t, y) = sum htilde_{k,i}(t) y^i has ONLY REAL y-roots for ALL t in (-1, 0)! ===
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 13):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

print("\n=== 4. Test if H_k(t, y) = sum htilde_{k,i}(t) y^i is REAL-ROOTED in y for every t in (-1, 0) ===")
for k in range(2, 12):
    ok_all = True
    fail_t = None
    for num in range(-99, 0):
        # evaluate at t = num/100 exactly!
        # Multiply ht[(k,i)] by 100^{2k-1}
        dmax = 2*k - 1
        y_poly = []
        for i in range(k + 1):
            val = sum(c * (num**m) * (100**(dmax - m)) for m, c in enumerate(ht[(k, i)]))
            y_poly.append(val)
        y_poly = poly_trim(y_poly)
        if not is_real_rooted(y_poly):
            ok_all = False
            fail_t = num / 100
            break
    print(f"  k={k:2d}: H_k(t, y) real-rooted in y for t in (-1, 0): {'PASS (all 99 pts)' if ok_all else f'FAIL at t={fail_t}'}")

# === 5. Test if F_{la, mu} << G_{la, mu} holds at k = 8, 9, 10, 11, 12! ===
print("\n=== 5. Test F_{la,mu} << G_{la,mu} for k = 8..12 (both + and - ratios) ===")
ratios = [(-100,1),(-10,1),(-2,1),(-1,1),(-1,2),(-1,10),(-1,100),(1,10),(1,1),(10,1)]
for k in range(8, 13):
    fails = []
    for i in range(2, k + 1):
        pi, pi1, pi2 = ht[(k, i)], ht[(k, i-1)], ht[(k, i-2)]
        for (la, mu) in ratios:
            F = add_poly(scale_poly(pi, la), scale_poly(pi1, mu))
            G = add_poly(scale_poly(pi1, la), scale_poly(pi2, mu))
            if mu < 0:
                F = scale_poly(F, -1)
                G = scale_poly(G, -1)
            ok, msg = exact_strict_interlaces(F, G)
            if not ok:
                fails.append((i, la, mu, msg))
                break
    print(f"  k={k:2d}: F_{{la,mu}} << G_{{la,mu}}: {'PASS ALL' if not fails else f'FAIL {fails[:3]}'}")
