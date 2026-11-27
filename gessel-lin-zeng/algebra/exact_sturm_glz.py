import sympy as sp

t = sp.Symbol('t')
x = sp.Symbol('x')
y = sp.Symbol('y')

def poly_from_coeffs(c):
    # Trim leading zeros; c[0] is constant term, c[-1] is leading term
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c

def apply_Dtilde_coeffs(j, c):
    # Dtilde_j(p) = t(1+t)p' + (1 + (j+1)t)p
    # coeff of t^m in p is c[m]
    # t*p' -> m*c[m] t^m
    # t^2*p' -> m*c[m] t^{m+1}
    # 1*p -> c[m] t^m
    # (j+1)t*p -> (j+1)*c[m] t^{m+1}
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j + 1) * cm
    return poly_from_coeffs(out)

def add_coeffs(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)):
        out[i] += a[i]
    for i in range(len(b)):
        out[i] += b[i]
    return poly_from_coeffs(out)

def exact_strict_interlaces(p, q):
    """Exact Hermite-Biehler / Sturm test for p << q (all simple real roots, lc > 0,
    deg q in {deg p, deg p + 1}, rightmost root of q > rightmost root of p).
    Polynomials represented as list of ints [c_0, c_1, ..., c_d]."""
    dp, dq = len(p) - 1, len(q) - 1
    if p[-1] <= 0 or q[-1] <= 0:
        return False, "nonpos_lc"
    if dq == 0:
        return dp == 0, "deg0"
    # Convert to Fraction for exact remainder sequence (or primitive integer PRS!)
    if dq == dp:
        # r1 = lc(q)*p - lc(p)*q must have deg d-1, lc > 0, and r1 << q
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

    # Signed pseudo-remainder sequence!
    # Since lc(r1) > 0, multiplying r0 by lc(r1)^2 > 0 preserves signs!
    # Wait: deg(r0) = deg(r1) + 1 always if strict interlacing holds!
    cur0 = r0[:]
    cur1 = r1[:]
    import math
    while len(cur1) > 1:
        d0, d1 = len(cur0) - 1, len(cur1) - 1
        if d0 != d1 + 1:
            return False, f"deg_drop({d0},{d1})"
        # Compute exact remainder of cur0 divided by cur1, scaled by positive factor lc(cur1)^2:
        # Step 1: c = lc(cur1)*cur0 - lc(cur0)*t*cur1 (deg <= d0 - 1 = d1)
        a0, a1 = cur0[-1], cur1[-1]
        tmp = [a1 * cur0[m] for m in range(d0 + 1)]
        for m in range(d1 + 1):
            tmp[m + 1] -= a0 * cur1[m]
        # Now tmp[d0] == 0; deg(tmp) <= d1
        # Step 2: rem_scaled = a1*tmp - tmp[d1]*cur1 (deg <= d1 - 1)
        b = tmp[d1]
        rem = [a1 * tmp[m] - b * cur1[m] for m in range(d1)]
        # Negate for Sturm!
        nxt = [-v for v in rem]
        while len(nxt) > 1 and nxt[-1] == 0:
            nxt.pop()
        if len(nxt) - 1 != d1 - 1 or nxt[-1] <= 0:
            return False, f"sturm_fail_at_deg_{d1-1}"
        # Divide by gcd to keep integers tiny!
        g = math.gcd(*[abs(v) for v in nxt])
        if g > 1:
            nxt = [v // g for v in nxt]
        cur0, cur1 = cur1, nxt
    return True, "ok"

KMAX = 20
g = {(1, 0): [1, 2], (1, 1): [1]}
S = {}
Dprev = {}
for k in range(2, KMAX + 1):
    for i in range(k + 1):
        t1 = apply_Dtilde_coeffs(k, g[(k-1, i)]) if i <= k - 1 else [0]
        if i <= k - 1:
            Dprev[(k, i)] = t1
        t2 = g[(k-1, i-1)] if i >= 1 else [0]
        S[(k, i)] = add_coeffs(t1, t2)
        g[(k, i)] = apply_Dtilde_coeffs(k, S[(k, i)])

def run_test(name, gen):
    tot = 0
    fails = []
    for tag, p, q in gen:
        tot += 1
        ok, msg = exact_strict_interlaces(p, q)
        if not ok:
            fails.append((tag, msg))
    print(f"{name:40s}: {tot - len(fails):3d}/{tot:3d} PASS | fails={fails[:5]}")

print(f"=== EXACT INTEGER STURM INTERLACING TESTS (k=1..{KMAX}) ===")
run_test("g_{k,i} << g_{k,i-1}", [((k,i), g[(k,i)], g[(k,i-1)]) for k in range(1,KMAX+1) for i in range(1,k+1)])
run_test("S_{k,i} << S_{k,i-1}", [((k,i), S[(k,i)], S[(k,i-1)]) for k in range(2,KMAX+1) for i in range(1,k+1)])
run_test("D_k g_{k-1,i} << D_k g_{k-1,i-1} (A<<C)", [((k,i), Dprev[(k,i)], Dprev[(k,i-1)]) for k in range(2,KMAX+1) for i in range(1,k)])
run_test("D_k g_{k-1,i} << g_{k-1,i-2} (A<<D)", [((k,i), Dprev[(k,i)], g[(k-1,i-2)]) for k in range(3,KMAX+1) for i in range(2,k)])
run_test("g_{k-1,i-1} << D_k g_{k-1,i} (B<<A)", [((k,i), g[(k-1,i-1)], Dprev[(k,i)]) for k in range(2,KMAX+1) for i in range(1,k)])
run_test("D_k g_{k-1,i} << g_{k-1,i-1} (A<<B)", [((k,i), Dprev[(k,i)], g[(k-1,i-1)]) for k in range(2,KMAX+1) for i in range(1,k)])
run_test("D_k g_{k-1,i} << S_{k,i-1}", [((k,i), Dprev[(k,i)], S[(k,i-1)]) for k in range(2,KMAX+1) for i in range(1,k)])
run_test("g_{k-1,i-1} << S_{k,i-1}", [((k,i), g[(k-1,i-1)], S[(k,i-1)]) for k in range(2,KMAX+1) for i in range(1,k+1)])
run_test("g_{k-1,i} << g_{k,i+1}", [((k,i), g[(k-1,i)], g[(k,i+1)]) for k in range(2,KMAX+1) for i in range(k)])
run_test("g_{k,i+1} << g_{k-1,i}", [((k,i), g[(k,i+1)], g[(k-1,i)]) for k in range(2,KMAX+1) for i in range(k)])
run_test("g_{k,i} << g_{k+1,i+1}", [((k,i), g[(k,i)], g[(k+1,i+1)]) for k in range(1,KMAX) for i in range(k+1)])
run_test("g_{k-1,i} << S_{k,i}", [((k,i), g[(k-1,i)], S[(k,i)]) for k in range(2,KMAX+1) for i in range(k)])
run_test("S_{k,i} << g_{k-1,i-2}", [((k,i), S[(k,i)], g[(k-1,i-2)]) for k in range(3,KMAX+1) for i in range(2,k+1)])
