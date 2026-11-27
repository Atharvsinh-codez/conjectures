import math

def poly_trim(c):
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c

def apply_Dtilde(j, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j + 1) * cm
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

KMAX = 16
g = {(1, 0): [1, 2], (1, 1): [1]}
S = {}
for k in range(2, KMAX + 1):
    for i in range(k + 1):
        t1 = apply_Dtilde(k, g[(k-1, i)]) if i <= k - 1 else [0]
        t2 = g[(k-1, i-1)] if i >= 1 else [0]
        S[(k, i)] = add_poly(t1, t2)
        g[(k, i)] = apply_Dtilde(k, S[(k, i)])

def run_test(name, pairs):
    fails = []
    for tag, p, q in pairs:
        ok, msg = exact_strict_interlaces(p, q)
        if not ok:
            fails.append((tag, msg))
    print(f"{name:44s}: {len(pairs)-len(fails):3d}/{len(pairs):3d} PASS | fails={fails[:5]}")

# Test intermediate S-chain and D_k vs D_{k+1} shifts!
run_test("S_{k,i-1} << g_{k,i} (= D_k S_{k,i})", [((k,i), S[(k,i-1)], g[(k,i)]) for k in range(2,KMAX+1) for i in range(1,k+1)])
run_test("g_{k,i} << S_{k,i-1}", [((k,i), g[(k,i)], S[(k,i-1)]) for k in range(2,KMAX+1) for i in range(1,k+1)])
run_test("g_{k,i} << S_{k,i-2}", [((k,i), g[(k,i)], S[(k,i-2)]) for k in range(2,KMAX+1) for i in range(2,k+1)])
run_test("S_{k,i-1} << D_{k+1} S_{k,i}", [((k,i), S[(k,i-1)], apply_Dtilde(k+1, S[(k,i)])) for k in range(2,KMAX+1) for i in range(1,k+1)])
run_test("D_{k+1} S_{k,i} << S_{k,i-2}", [((k,i), apply_Dtilde(k+1, S[(k,i)]), S[(k,i-2)]) for k in range(2,KMAX+1) for i in range(2,k+1)])
run_test("D_k D_{k+1} S_{k,i} << D_{k+1} D_k S_{k,i}", [((k,i), apply_Dtilde(k, apply_Dtilde(k+1, S[(k,i)])), apply_Dtilde(k+1, apply_Dtilde(k, S[(k,i)]))) for k in range(2,KMAX+1) for i in range(k+1)])
run_test("D_{k+1} D_k S_{k,i} << D_k D_{k+1} S_{k,i}", [((k,i), apply_Dtilde(k+1, apply_Dtilde(k, S[(k,i)])), apply_Dtilde(k, apply_Dtilde(k+1, S[(k,i)]))) for k in range(2,KMAX+1) for i in range(k+1)])

# Also test: g_{k, i-1} << D_k g_{k, i} (with SAME k instead of k+1!)
run_test("g_{k,i-1} << D_k g_{k,i}", [((k,i), g[(k,i-1)], apply_Dtilde(k, g[(k,i)])) for k in range(1,KMAX+1) for i in range(1,k+1)])
run_test("D_k g_{k,i} << g_{k,i-1}", [((k,i), apply_Dtilde(k, g[(k,i)]), g[(k,i-1)]) for k in range(1,KMAX+1) for i in range(1,k+1)])
run_test("D_k g_{k,i} << g_{k,i-2}", [((k,i), apply_Dtilde(k, g[(k,i)]), g[(k,i-2)]) for k in range(2,KMAX+1) for i in range(2,k+1)])

# Test with arbitrary parameter j: for which j does g_{k,i-1} << D_j g_{k,i} << g_{k,i-2}?
print("\n=== Testing parameter j in g_{k,i-1} << D_j g_{k,i} << g_{k,i-2} ===")
for dj in [-2, -1, 0, 1, 2, 3]:
    run_test(f"g_{{k,i-1}} << D_{{k+{dj}}} g_{{k,i}}", [((k,i), g[(k,i-1)], apply_Dtilde(k+dj, g[(k,i)])) for k in range(3,KMAX+1) for i in range(1,k+1)])
    run_test(f"D_{{k+{dj}}} g_{{k,i}} << g_{{k,i-2}}", [((k,i), apply_Dtilde(k+dj, g[(k,i)]), g[(k,i-2)]) for k in range(3,KMAX+1) for i in range(2,k+1)])
