import math
import numpy as np

def poly_trim(c):
    while len(c) > 1 and c[-1] == 0: c.pop()
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
    if p[-1] <= 0 or q[-1] <= 0: return False, "nonpos_lc"
    if dq == 0: return dp == 0, "deg0"
    if dq == dp:
        lcp, lcq = p[-1], q[-1]
        r1 = [lcq * p[m] - lcp * q[m] for m in range(dq)]
        while len(r1) > 1 and r1[-1] == 0: r1.pop()
        if len(r1) - 1 != dq - 1 or r1[-1] <= 0: return False, "eq_lc_fail"
        r0 = q[:]
    elif dq == dp + 1:
        r0, r1 = q[:], p[:]
    else:
        return False, f"deg_jump({dp}->{dq})"
    cur0, cur1 = r0[:], r1[:]
    while len(cur1) > 1:
        d0, d1 = len(cur0) - 1, len(cur1) - 1
        if d0 != d1 + 1: return False, f"deg_drop({d0},{d1})"
        a0, a1 = cur0[-1], cur1[-1]
        tmp = [a1 * cur0[m] for m in range(d0 + 1)]
        for m in range(d1 + 1): tmp[m + 1] -= a0 * cur1[m]
        b = tmp[d1]
        nxt = [-(a1 * tmp[m] - b * cur1[m]) for m in range(d1)]
        while len(nxt) > 1 and nxt[-1] == 0: nxt.pop()
        if len(nxt) - 1 != d1 - 1 or nxt[-1] <= 0: return False, f"sturm_fail_deg{d1-1}"
        g = math.gcd(*[abs(v) for v in nxt])
        if g > 1: nxt = [v // g for v in nxt]
        cur0, cur1 = cur1, nxt
    return True, "ok"

# Compute htilde_{k, i} for k=1..14
KMAX = 14
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, KMAX + 1):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

# === TEST 1: Is E_{k, i} = D_{k+1}^i ht_{k-1, i} (0 <= i <= k-1) a STRONG (all-pairs!) interlacing sequence?! ===
print("=== TEST 1: Is E_{k, i} = D_{k+1}^i ht_{k-1, i} a STRONG interlacing sequence (all a < b: E_a << E_b)?! ===")
for k in range(2, KMAX + 1):
    E = []
    for i in range(k):
        p = ht[(k-1, i)][:]
        for _ in range(i):
            p = apply_D(k + 1, p)
        E.append(p)
    adj_fails = []
    strong_fails = []
    for a in range(k):
        for b in range(a + 1, k):
            ok, msg = exact_strict_interlaces(E[a], E[b])
            if not ok:
                if b == a + 1:
                    adj_fails.append((a, b, msg))
                strong_fails.append((a, b, msg))
    print(f"  k={k:2d} (deg {2*k-3:2d}): adj fails={len(adj_fails):2d} {adj_fails[:2]} | ALL-PAIRS strong fails={len(strong_fails):2d} {strong_fails[:3]}")

# === TEST 2: Also test F_{k, i} = D_{k+1}^i ht_{k, i} (or D_{k+2}^i ht_{k, i})! ===
print("\n=== TEST 2: Test E_step = [D_{k+1}^i ht_{k, i}] (i=0..k) strong interlacing ===")
for k in range(2, 12):
    F = []
    for i in range(k + 1):
        p = ht[(k, i)][:]
        for _ in range(i):
            p = apply_D(k + 1, p)
        F.append(p)
    strong_fails = [(a, b) for a in range(k+1) for b in range(a+1, k+1) if not exact_strict_interlaces(F[a], F[b])[0]]
    print(f"  k={k:2d}: D_{{k+1}}^i ht_{{k,i}} strong fails = {len(strong_fails)} {strong_fails[:4]}")

# === TEST 3: Check REAL STABILITY of sum E_{k, i}(t) y^i AND sum binom(k-1,i) E_{k,i}(t) y^i AND Q = X + 2My + Yy^2! ===
print("\n=== TEST 3: Check t=0 interlacing P1 << P0 (up to k=45!) for E_{k,i} (unweighted, binom-weighted, and quadratic X+2My+Yy^2) ===")
# Track (c0, c1) of ht_{k, i} up to k=45
c0 = {(1, 0): 1, (1, 1): 1}
c1 = {(1, 0): 2, (1, 1): 0}
for k in range(2, 46):
    j = k + 1
    for i in range(k + 1):
        a0 = c0[(k-1, i)] if i <= k - 1 else 0
        a1 = (2 * c1[(k-1, i)] + j * c0[(k-1, i)]) if i <= k - 1 else 0
        s0 = a0 + (c0[(k-1, i-1)] if i >= 1 else 0)
        s1 = a1 + (c1[(k-1, i-1)] if i >= 1 else 0)
        c0[(k, i)] = s0
        c1[(k, i)] = 2 * s1 + j * s0

def apply_D_c01(j, v0, v1, times):
    for _ in range(times):
        v0, v1 = v0, 2 * v1 + j * v0
    return v0, v1

# Test binom-weighted E_{k,i}(t) at t=0: P0(y) = sum binom(k-1,i)*E0_i y^i, P1(y) = sum binom(k-1,i)*E1_i y^i
for mode in ["unweighted E", "binom E", "1/i! E"]:
    first_fail = None
    for k in range(3, 30):
        E0, E1 = [], []
        for i in range(k):
            v0, v1 = apply_D_c01(k + 1, c0[(k-1, i)], c1[(k-1, i)], i)
            if mode == "binom E":
                w = math.comb(k - 1, i)
            elif mode == "1/i! E":
                w = math.factorial(k - 1) // math.factorial(i)
            else:
                w = 1
            E0.append(w * v0)
            E1.append(w * v1)
        ok1, m1 = exact_strict_interlaces(E1, E0)
        ok2, m2 = exact_strict_interlaces(E0, E1)
        if not (ok1 or ok2):
            first_fail = (k, m1, m2)
            break
    print(f"  {mode:15s}: P1 vs P0 interlacing up to k=29 -> {'PASS ALL!' if first_fail is None else f'FAIL at k={first_fail[0]} ({first_fail[1]}, {first_fail[2]})'}")

# Test quadratic Q_{k,i}(t, y) = X + 2 M y + Y y^2 at t=0 up to k=45!
from decimal import Decimal, getcontext
getcontext().prec = 100
first_fail_Q = None
for k in range(3, 46):
    for i in range(2, k):
        x0, x1 = c0[(k-1, i-2)], c1[(k-1, i-2)]
        m0, m1 = apply_D_c01(k + 1, c0[(k-1, i-1)], c1[(k-1, i-1)], 1)
        y0, y1 = apply_D_c01(k + 1, c0[(k-1, i)], c1[(k-1, i)], 2)
        # Check P0(y) = x0 + 2*m0*y + y0*y^2 and P1(y) = x1 + 2*m1*y + y1*y^2
        # Check if P0 << P1 or P1 << P0!
        ok1, _ = exact_strict_interlaces([x1, 2*m1, y1], [x0, 2*m0, y0])
        ok2, _ = exact_strict_interlaces([x0, 2*m0, y0], [x1, 2*m1, y1])
        if not (ok1 or ok2):
            first_fail_Q = (k, i)
            break
    if first_fail_Q: break
print(f"  Quad X+2My+Yy^2: P1 vs P0 interlacing up to k=45 -> {'PASS ALL!' if first_fail_Q is None else f'FAIL at {first_fail_Q}'}")

# === TEST 4: Adversarial UHP real-stability test of sum binom(k-1,i) E_{k,i}(t) y^i AND X + 2My + Yy^2 for k=3..11! ===
print("\n=== TEST 4: UHP complex real-stability test of binom-weighted E_{k,i} and Quad X+2My+Yy^2 (k=3..11) ===")
def eval_c(p, z):
    acc = 0j
    for c in reversed(p): acc = acc * z + c
    return acc
test_ts = [complex(re, im) for re in [-10, -2, -1.05, -0.95, -0.5, -0.1, 0, 0.1, 1, 10] for im in [1e-4, 0.02, 0.5]]

for k in range(3, 12):
    E = []
    for i in range(k):
        p = ht[(k-1, i)][:]
        for _ in range(i): p = apply_D(k + 1, p)
        E.append(p)
    # Check binom E
    max_im_E = -1e9
    for tv in test_ts:
        yc = [math.comb(k-1, i) * eval_c(E[i], tv) for i in range(k-1, -1, -1)]
        mx = max(abs(c) for c in yc)
        yc = [c / mx for c in yc]
        mi = max(z.imag for z in np.roots(yc))
        if mi > max_im_E: max_im_E = mi
    # Check all quads X + 2My + Yy^2 (which are E[i-2] + 2 E[i-1] y + E[i] y^2, Wait: without extra D^{i-2}!)
    max_im_Q = -1e9
    for i in range(2, k):
        X = ht[(k-1, i-2)]
        M = apply_D(k+1, ht[(k-1, i-1)])
        Y = apply_D(k+1, apply_D(k+1, ht[(k-1, i)]))
        for tv in test_ts:
            yc = [eval_c(Y, tv), 2*eval_c(M, tv), eval_c(X, tv)]
            mx = max(abs(c) for c in yc)
            yc = [c / mx for c in yc]
            mi = max(z.imag for z in np.roots(yc))
            if mi > max_im_Q: max_im_Q = mi
    print(f"  k={k:2d}: binom-E max Im(y)={max_im_E:+.3e} | Quad X+2My+Yy^2 max Im(y)={max_im_Q:+.3e}")
