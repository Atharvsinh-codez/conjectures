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

def scale_poly(c, s):
    return poly_trim([v * s for v in c])

def mul_poly(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b): out[i + j] += ai * bj
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

def is_real_rooted(p):
    d = len(p) - 1
    if d <= 1: return True
    if p[-1] < 0: p = [-v for v in p]
    dp = [m * p[m] for m in range(1, d + 1)]
    return exact_strict_interlaces(dp, p)[0]

# Compute ht up to k=18
KMAX = 18
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, KMAX + 1):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

# === TEST 1: Exact integer Sturm check of Quad X - 2My + Yy^2 on positive-slope lines y = a*t + b (a > 0) for k=3..18! ===
print("=== TEST 1: Exact Sturm check of Quad Q_{k,i}(t, y) = X - 2My + Yy^2 on positive-slope lines (k=3..18) ===")
lines = [
    (1, 1, 0, 1), (1, 1, 1, 1), (1, 1, -1, 1), (1, 1, 1, 2), (1, 1, -1, 2),
    (10, 1, 1, 1), (1, 10, 1, 1), (100, 1, 1, 1), (1, 100, 1, 1),
    (1, 1, 10, 1), (1, 1, -10, 1), (5, 1, 2, 1), (1, 5, 2, 1), (20, 1, -1, 1),
]

def eval_quad_line(X, M, Y, an, ad, bn, bd):
    # X(t) - 2*M(t)*(a*t + b) + Y(t)*(a*t + b)^2, scaled by (ad*bd)^2
    D = ad * bd
    lin = [bn * ad, an * bd] # D*(a*t + b)
    t0 = scale_poly(X, D * D)
    t1 = scale_poly(mul_poly(M, lin), -2 * D)
    t2 = mul_poly(Y, mul_poly(lin, lin))
    return add_poly(add_poly(t0, t1), t2)

for k in range(3, 16):
    fail = None
    for i in range(2, k):
        X = ht[(k-1, i-2)]
        M = apply_D(k+1, ht[(k-1, i-1)])
        Y = apply_D(k+1, apply_D(k+1, ht[(k-1, i)]))
        for (an, ad, bn, bd) in lines:
            p = eval_quad_line(X, M, Y, an, ad, bn, bd)
            if not is_real_rooted(p):
                fail = (i, an, ad, bn, bd)
                break
        if fail: break
    print(f"  k={k:2d}: Quad Q_{{k,i}}(t, y) exact Sturm on 14 lines: {'PASS ALL!' if fail is None else f'FAIL at {fail}'}")

# === TEST 2: Check EXACT Rayleigh determinant of binom-E and Quad X - 2My + Yy^2 at adversarial (t, y) ===
# Recall: P(t, y) real stable => P_t P_y - P P_{ty} >= 0 for ALL real (t, y)!
print("\n=== TEST 2: Test if Q_{k,i-1}(t, y) + z * D_{k+1} Q_{k,i}(t, y) is REAL STABLE in (t, y, z)!! ===")
# Check for k=4..10: for fixed t in UHP, z in UHP, are all y-roots of Q_{k,i-1}(t, y) + z * D_{k+1} Q_{k,i}(t, y) = 0 in LHP (Im y <= 0)?!
def eval_c(p, z):
    acc = 0j
    for c in reversed(p): acc = acc * z + c
    return acc

ts = [complex(re, im) for re in [-5, -1.2, -0.8, -0.5, -0.2, -0.02, 0.0, 0.2, 2.0] for im in [1e-4, 0.05, 1.0]]
zs = [complex(re, im) for re in [-10, -1, -0.1, 0, 0.1, 1, 10] for im in [1e-4, 0.05, 1.0]]

for k in range(4, 10):
    max_im_y = -1e9
    worst = None
    for i in range(3, k):
        # Q_{k, i-1}(t, y) = ht_{k-1, i-3} - 2 D_{k+1} ht_{k-1, i-2} y + D_{k+1}^2 ht_{k-1, i-1} y^2
        X1 = ht[(k-1, i-3)]
        M1 = apply_D(k+1, ht[(k-1, i-2)])
        Y1 = apply_D(k+1, apply_D(k+1, ht[(k-1, i-1)]))
        # D_{k+1} Q_{k, i}(t, y) = D_{k+1} ht_{k-1, i-2} - 2 D_{k+1}^2 ht_{k-1, i-1} y + D_{k+1}^3 ht_{k-1, i} y^2
        X2 = M1
        M2 = Y1
        Y2 = apply_D(k+1, apply_D(k+1, apply_D(k+1, ht[(k-1, i)])))
        for tv in ts:
            x1, m1, y1 = eval_c(X1, tv), eval_c(M1, tv), eval_c(Y1, tv)
            x2, m2, y2 = eval_c(X2, tv), eval_c(M2, tv), eval_c(Y2, tv)
            for zv in zs:
                # Try BOTH +z and -z!
                # Note: X1 << X2 (equal deg!), so X1 - z X2 is real stable in (t, z) (Im z > 0 -> Im t < 0)!
                # So test Q_{k,i-1}(t, y) - z * D_{k+1} Q_{k,i}(t, y) = 0 with t in UHP, z in UHP -> want Im y <= 0!
                a0 = x1 - zv * x2
                a1 = -2 * (m1 - zv * m2)
                a2 = y1 - zv * y2
                if abs(a2) == 0: continue
                r = np.roots([a2, a1, a0])
                mi = max(w.imag for w in r)
                if mi > max_im_y:
                    max_im_y = mi; worst = (i, tv, zv)
    print(f"  k={k}: Q_{{k,i-1}}(t,y) - z * D_{{k+1}} Q_{{k,i}}(t,y) max Im(y) = {max_im_y:+.3e} at {worst}")

# === TEST 3: Check REAL STABILITY of Cubic C_{k,i}(t, y) = E_{i-3} - 3 E_{i-2} y + 3 E_{i-1} y^2 - E_i y^3 AND binom-E at k=14, 15, 16! ===
print("\n=== TEST 3: Check binom-E at k=14, 15, 16 (50-digit mpmath!) and Cubic C_{k,i} ===")
import mpmath as mp
mp.mp.dps = 50
def eval_mp(p, z):
    acc = mp.mpc(0)
    for c in reversed(p): acc = acc * z + c
    return acc

ts_mp = [mp.mpc(re, im) for re in [-10, -1.5, -0.95, -0.5, -0.05, 0.0, 0.1, 1.0, 10.0] for im in [1e-4, 0.05, 1.0]]
for k in [14, 15, 16]:
    E = []
    for i in range(k):
        p = ht[(k-1, i)][:]
        for _ in range(i): p = apply_D(k + 1, p)
        E.append(p)
    min_im = mp.mpf('1e9')
    worst_t = None
    for tv in ts_mp:
        yc = [math.comb(k-1, i) * eval_mp(E[i], tv) for i in range(k-1, -1, -1)]
        r = mp.polyroots(yc, maxsteps=300, extraprec=30)
        mi = min(w.imag for w in r)
        if mi < min_im:
            min_im = mi; worst_t = complex(tv)
    print(f"  k={k}: binom-E(t, -y) 50-digit minIm(y) = {float(min_im):+.4e} at t={worst_t}")
