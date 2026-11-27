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

# === TEST 1: Exact t=0 AND t=-1 AND t=infinity necessary stability checks for binom-E, 1/i!-E, and Quad up to k=50! ===
print("=== TEST 1: Necessary stability checks at t=0, t=-1, and t=inf for binom-E, 1/i!-E, Quad X-2My+Yy^2 ===")
# Track (c0, c1) at t=0, AND (lc, sub_lc) at t=infinity!
c0 = {(1, 0): 1, (1, 1): 1}
c1 = {(1, 0): 2, (1, 1): 0}
lc = {(1, 0): 2, (1, 1): 1}
slc = {(1, 0): 1, (1, 1): 0}
# Also track at t=-1: ht(t) = v0 + v1*(1+t) + ...
# Under D_j(f) = t(1+t)f' + (1+jt)f = ((1+t)-1)(1+t)f' + ((1-j) + j(1+t))f:
# At u = 1+t: D_j(v0 + v1*u + ...) = (u^2 - u)(v1 + ...) + ((1-j) + j*u)(v0 + v1*u + ...)
# = (1-j)*v0 + ((j-2)*v1 + j*v0)*u + ... (WAIT: check sign! -u*v1 + (1-j)*v1*u = -j*v1*u! Let's check carefully!)
# u^0: (1-j)*v0
# u^1: (-v1) + (1-j)*v1 + j*v0 = j*v0 - j*v1!
# Wait: check on f = 1 (v0=1, v1=0): D_j(1) = 1+jt = (1-j) + j*u -> u^0 = 1-j, u^1 = j! Match!
# Check on f = u (v0=0, v1=1): D_j(u) = (u^2-u)(1) + ((1-j)+ju)u = -j*u + (j+1)*u^2 -> u^0=0, u^1=-j! Match!
m0 = {(1, 0): -1, (1, 1): 1} # ht(-1)
m1 = {(1, 0): 2, (1, 1): 0}  # ht'(-1)

for k in range(2, 51):
    j = k + 1
    for i in range(k + 1):
        # t = 0
        a0 = c0[(k-1, i)] if i <= k - 1 else 0
        a1 = (2 * c1[(k-1, i)] + j * c0[(k-1, i)]) if i <= k - 1 else 0
        s0 = a0 + (c0[(k-1, i-1)] if i >= 1 else 0)
        s1 = a1 + (c1[(k-1, i-1)] if i >= 1 else 0)
        c0[(k, i)] = s0
        c1[(k, i)] = 2 * s1 + j * s0
        # t = inf
        if i <= k - 1:
            d1 = 2*k - i - 3
            al = (d1 + j) * lc[(k-1, i)]
            as_ = (d1 + 1) * lc[(k-1, i)] + (d1 + j - 1) * slc[(k-1, i)]
        else:
            al, as_ = 0, 0
        sl = al + (lc[(k-1, i-1)] if i >= 1 else 0)
        ss = as_ + (slc[(k-1, i-1)] if i >= 1 else 0)
        ds = 2*k - i - 2
        lc[(k, i)] = (ds + j) * sl
        slc[(k, i)] = (ds + 1) * sl + (ds + j - 1) * ss
        # t = -1
        if i <= k - 1:
            am0 = (1 - j) * m0[(k-1, i)]
            am1 = j * m0[(k-1, i)] - j * m1[(k-1, i)]
        else:
            am0, am1 = 0, 0
        sm0 = am0 + (m0[(k-1, i-1)] if i >= 1 else 0)
        sm1 = am1 + (m1[(k-1, i-1)] if i >= 1 else 0)
        m0[(k, i)] = (1 - j) * sm0
        m1[(k, i)] = j * sm0 - j * sm1

def apply_D_t0(j, v0, v1, times):
    for _ in range(times): v0, v1 = v0, 2*v1 + j*v0
    return v0, v1
def apply_D_tinf(j, d, l, s, times):
    for _ in range(times):
        l, s = (d + j)*l, (d + 1)*l + (d + j - 1)*s
        d += 1
    return l, s
def apply_D_tm1(j, v0, v1, times):
    for _ in range(times): v0, v1 = (1-j)*v0, j*v0 - j*v1
    return v0, v1

# Remember for P(t, -y) = sum w_i E_i(t) (-y)^i = sum w_i E_i(t) z^i (with z = -y):
# Wait: P(t, y_realstable) has y_realstable = -y, so z = y_realstable!
# In z: P(t, z) = sum w_i E_i(t) z^i!
# For P(t, z) to be real stable in (t, z), at t=0 we need P_t(0, z) << P(0, z) OR P(0, z) << P_t(0, z)!
# At t=inf (set u=1/t): P(1/u, z) ~ u^{-D} (L0(z) + u L1(z) + ...), need L0 and L1 interlace in z!
# At t=-1 (set u=1+t): P(-1+u, z) = M0(z) + u M1(z) + ..., need M0 and M1 interlace in z!
for name, weight_fn in [
    ("binom-E", lambda k, i: math.comb(k-1, i)),
    ("1/i!-E",  lambda k, i: math.factorial(k-1) // math.factorial(i)),
]:
    fail_0, fail_inf, fail_m1 = None, None, None
    for k in range(3, 35):
        j = k + 1
        P0, P1, L0, L1, M0, M1 = [], [], [], [], [], []
        for i in range(k):
            w = weight_fn(k, i)
            v0, v1 = apply_D_t0(j, c0[(k-1, i)], c1[(k-1, i)], i)
            P0.append(w * v0); P1.append(w * v1)
            l, s = apply_D_tinf(j, 2*k - i - 3, lc[(k-1, i)], slc[(k-1, i)], i)
            L0.append(w * l); L1.append(w * s)
            u0, u1 = apply_D_tm1(j, m0[(k-1, i)], m1[(k-1, i)], i)
            # Notice at t=-1, ht_{k-1, i}(-1) has sign (-1)^{2k-i-3}, and D_j multiplies by (1-j) < 0,
            # so E_{k, i}(-1) has CONSTANT sign (-1)^{2k-3} = -1 for ALL i! Make positive!
            M0.append(-w * u0); M1.append(w * u1)
        if fail_0 is None:
            ok = exact_strict_interlaces(P0, P1)[0] or exact_strict_interlaces(P1, P0)[0]
            if not ok: fail_0 = k
        if fail_inf is None:
            ok = exact_strict_interlaces(L0, L1)[0] or exact_strict_interlaces(L1, L0)[0]
            if not ok: fail_inf = k
        if fail_m1 is None:
            # Make M1 positive lc if needed
            if M1[-1] < 0: M1 = [-v for v in M1]
            ok = exact_strict_interlaces(M0, M1)[0] or exact_strict_interlaces(M1, M0)[0]
            if not ok: fail_m1 = k
    print(f"  {name:10s} (k=3..34): t=0 fail={fail_0} | t=inf fail={fail_inf} | t=-1 fail={fail_m1}")

# Also test Quad X - 2My + Yy^2 at t=0, t=inf, t=-1 up to k=50!
f0, finf, fm1 = None, None, None
for k in range(3, 51):
    j = k + 1
    for i in range(2, k):
        # X = ht_{k-1,i-2}, M = D_j ht_{k-1,i-1}, Y = D_j^2 ht_{k-1,i}
        x0, x1 = c0[(k-1, i-2)], c1[(k-1, i-2)]
        m0_, m1_ = apply_D_t0(j, c0[(k-1, i-1)], c1[(k-1, i-1)], 1)
        y0, y1 = apply_D_t0(j, c0[(k-1, i)], c1[(k-1, i)], 2)
        if f0 is None and not (exact_strict_interlaces([x0,2*m0_,y0],[x1,2*m1_,y1])[0] or exact_strict_interlaces([x1,2*m1_,y1],[x0,2*m0_,y0])[0]):
            f0 = (k, i)
        xl, xs = lc[(k-1, i-2)], slc[(k-1, i-2)]
        ml, ms = apply_D_tinf(j, 2*k-i-2, lc[(k-1, i-1)], slc[(k-1, i-1)], 1)
        yl, ys = apply_D_tinf(j, 2*k-i-3, lc[(k-1, i)], slc[(k-1, i)], 2)
        if finf is None and not (exact_strict_interlaces([xl,2*ml,yl],[xs,2*ms,ys])[0] or exact_strict_interlaces([xs,2*ms,ys],[xl,2*ml,yl])[0]):
            finf = (k, i)
        # t=-1: sign of X(-1), M(-1), Y(-1) is (-1)^{2k-i-1}
        sgn = 1 if (2*k - i - 1) % 2 == 0 else -1
        xu0, xu1 = sgn * m0[(k-1, i-2)], -sgn * m1[(k-1, i-2)]
        mu0, mu1 = apply_D_tm1(j, m0[(k-1, i-1)], m1[(k-1, i-1)], 1)
        mu0, mu1 = sgn * mu0, -sgn * mu1
        yu0, yu1 = apply_D_tm1(j, m0[(k-1, i)], m1[(k-1, i)], 2)
        yu0, yu1 = sgn * yu0, -sgn * yu1
        if fm1 is None and not (exact_strict_interlaces([xu0,2*mu0,yu0],[xu1,2*mu1,yu1])[0] or exact_strict_interlaces([xu1,2*mu1,yu1],[xu0,2*mu0,yu0])[0]):
            fm1 = (k, i)
print(f"  Quad X-2My+Yy^2 (k=3..50): t=0 fail={f0} | t=inf fail={finf} | t=-1 fail={fm1}")

# === TEST 2: Test binom-E, 1/i!-E, and Quad X-2My+Yy^2 on WILD [2, 100, 2, 100, 2, 100]! ===
print("\n=== TEST 2: Test binom-E, 1/i!-E, and Quad X-2My+Yy^2 on WILD sequences ===")
def eval_c(p, z):
    acc = 0j
    for c in reversed(p): acc = acc * z + c
    return acc
test_ts = [complex(re, im) for re in [-10, -2, -1.05, -0.95, -0.5, -0.1, 0, 0.1, 1, 10] for im in [1e-4, 0.02, 0.5]]

for name, bs in [
    ("GLZ inner [2,3,4,5,6,7]", [2,3,4,5,6,7]),
    ("Wild [2,100,2,100,2,100]", [2,100,2,100,2,100]),
    ("Drop [100,50,25,12,6,3]", [100,50,25,12,6,3]),
]:
    row = [[1]]
    for bj in bs[:-1]:
        k = len(row)
        row = [add_poly(apply_D(bj, row[i]) if i <= k-1 else [0], row[i-1] if i >= 1 else [0]) for i in range(k+1)]
    b_next = bs[-1]
    k = len(row)
    E = []
    for i in range(k):
        p = row[i][:]
        for _ in range(i): p = apply_D(b_next, p)
        E.append(p)
    # Check min Im(y) for binom-E and 1/i!-E!
    min_b, min_f = 1e9, 1e9
    for tv in test_ts:
        yb = [math.comb(k-1, i) * eval_c(E[i], tv) for i in range(k-1, -1, -1)]
        yb = [c / max(abs(z) for z in yb) for c in yb]
        min_b = min(min_b, min(z.imag for z in np.roots(yb)))
        yf = [(math.factorial(k-1)//math.factorial(i)) * eval_c(E[i], tv) for i in range(k-1, -1, -1)]
        yf = [c / max(abs(z) for z in yf) for c in yf]
        min_f = min(min_f, min(z.imag for z in np.roots(yf)))
    print(f"  {name:25s}: binom-E minIm={min_b:+.3e} | 1/i!-E minIm={min_f:+.3e}")
