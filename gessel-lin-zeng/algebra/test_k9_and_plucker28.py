"""
1. Push our k=8 non-interlacing universal prod(T_j + w_j y)(1) counterexample to k=9 REAL-ROOTEDNESS failure!
2. Also test if fixed a_j = 1, w_j = 1 (varying ONLY b_j > 1) can break universal Turan/interlacing/RR!
3. Check GLZ Quad J^Q = 4 W(X,M) W(M,Y) - W(X,Y)^2 > 0 on ALL of R AND Plucker ratio W_AB W_CE / W_BC^2 up to k=30!
"""
import math, random
import mpmath as mp
mp.mp.dps = 60

def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p
def poly_add(a, b):
    n = max(len(a), len(b))
    return poly_trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def poly_deriv(a):
    return [i * a[i] for i in range(1, len(a))] if len(a) > 1 else [0]
def apply_T(p, a, b):
    dp = poly_deriv(p)
    res = [0] * (len(p) + 2)
    for i, c in enumerate(p):
        res[i] += a * c
        res[i+1] += b * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c
        res[i+2] += c
    return poly_trim(res)
def poly_eval(p, x):
    acc = mp.mpf(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc

print("=== 1. PUSH UNIVERSAL k=8 COUNTEREXAMPLE TO k=9 COMPLEX ROOTS ===", flush=True)
# Rebuild exact trial 1 from seed(7)
random.seed(7)
for tr in range(2):
    k = 6
    polys = {0: [mp.mpf(1)]}
    for j in range(1, k+1):
        a_j = mp.mpf(10 ** random.uniform(-2, 2))
        b_j = a_j + mp.mpf(10 ** random.uniform(-2, 2))
        w_j = mp.mpf(10 ** random.uniform(-2, 2))
        new_p = {}
        for i in range(j+1):
            p1 = apply_T(polys[i], a_j, b_j) if i in polys else [mp.mpf(0)]
            p2 = [w_j * c for c in polys[i-1]] if (i-1) in polys else [mp.mpf(0)]
            new_p[i] = poly_add(p1, p2)
        polys = new_p

i = 2
tv = -mp.mpf('0.001')
dp = {idx: poly_deriv(polys[idx]) for idx in range(7)}
f0, df0 = poly_eval(polys[i-1], tv), poly_eval(dp[i-1], tv)
f1, df1 = poly_eval(polys[i], tv),   poly_eval(dp[i], tv)
f2, df2 = poly_eval(polys[i+1], tv), poly_eval(dp[i+1], tv)
delta = f1*f1 - f0*f2
ddelta_u = tv*(1+tv)*(2*f1*df1 - df0*f2 - f0*df2)
phi_star = -ddelta_u / (2*delta)
a7 = mp.mpf('0.01')
b7 = a7 + (phi_star - a7*(1+tv))/tv
w7 = mp.mpf('1e-10')
p7 = {idx: poly_add(apply_T(polys[idx], a7, b7) if idx in polys else [mp.mpf(0)],
                    [w7*c for c in polys[idx-1]] if (idx-1) in polys else [mp.mpf(0)]) for idx in range(8)}
phi8 = mp.mpf(-100)
a8 = mp.mpf('0.01')
b8 = a8 + (phi8 - a8*(1+tv))/tv
dv_ip1 = poly_eval(poly_deriv(p7[i+1]), tv)
w8 = -(tv*(1+tv)*dv_ip1 + phi8*poly_eval(p7[i+1], tv)) / poly_eval(p7[i], tv)
q8_3 = poly_add(apply_T(p7[3], a8, b8), [w8*c for c in p7[2]])
q8_2 = poly_add(apply_T(p7[2], a8, b8), [w8*c for c in p7[1]])

# Now scan (a9, b9, w9) to get complex roots in q9_3 = T9(q8_3) + w9 * q8_2!
# Note q8_3 has roots at -0.005719, -0.001, -8e-8; q8_2 has TWO roots at -0.00190, -0.001007!
# So we want T9(q8_3) to keep its root between -0.005719 and -0.001 CLOSE to -0.005719, or shift the -0.001 root right,
# so that on (-0.0019, -0.0010) T9(q8_3) and q8_2 push a pair of roots into the complex plane!
found_cx = False
for a9 in [mp.mpf('0.001'), mp.mpf('0.1'), mp.mpf('1'), mp.mpf('10'), mp.mpf('100')]:
    for beta9 in [mp.mpf('0.001'), mp.mpf('1'), mp.mpf('100'), mp.mpf('1e4'), mp.mpf('1e6')]:
        b9 = a9 + beta9
        Tq = apply_T(q8_3, a9, b9)
        for p_w in range(-6, 8):
            w9 = mp.mpf(10)**p_w
            q9 = poly_add(Tq, [w9*c for c in q8_2])
            rts9 = mp.polyroots(list(reversed(q9)), maxsteps=200, extraprec=30)
            mx = max(abs(r.imag) for r in rts9)
            if mx > mp.mpf('1e-10'):
                print(f"  BOOM!! COMPLEX ROOTS AT k=9! a9={a9}, b9={b9}, w9={w9}, max|Im|={mp.nstr(mx, 8)}", flush=True)
                print(f"  Nonreal roots: {[complex(r) for r in rts9 if abs(r.imag) > mp.mpf('1e-10')]}", flush=True)
                found_cx = True
                break
        if found_cx: break
    if found_cx: break
if not found_cx:
    print("  No complex roots found in coarse scan at k=9 (check finer w9!)", flush=True)

print("\n=== 2. CHECK FIXED a_j=1, w_j=1 (VARYING ONLY b_j > 1): CAN J_i < 0 ON (-1,0)? ===", flush=True)
random.seed(42)
tgrid_fast = [-mp.mpf(m)/25 for m in range(1, 25)] + [-mp.mpf('1e-4'), -mp.mpf('0.9999')]
neg_J_fixed = 0
min_J_fixed = mp.mpf('1e100')
for tr in range(40):
    k = 8
    polys = {0: [mp.mpf(1)]}
    for j in range(1, k+1):
        b_j = mp.mpf(1) + mp.mpf(10 ** random.uniform(-1, 3))
        new_p = {}
        for idx in range(j+1):
            p1 = apply_T(polys[idx], mp.mpf(1), b_j) if idx in polys else [mp.mpf(0)]
            p2 = polys[idx-1] if (idx-1) in polys else [mp.mpf(0)]
            new_p[idx] = poly_add(p1, p2)
        polys = new_p
    dp = {idx: poly_deriv(polys[idx]) for idx in range(k+1)}
    bad = False
    for idx in range(1, k):
        for tv in tgrid_fast:
            f0, df0 = poly_eval(polys[idx-1], tv), poly_eval(dp[idx-1], tv)
            f1, df1 = poly_eval(polys[idx], tv),   poly_eval(dp[idx], tv)
            f2, df2 = poly_eval(polys[idx+1], tv), poly_eval(dp[idx+1], tv)
            w01 = f0*df1 - df0*f1
            w12 = f1*df2 - df1*f2
            w02 = f0*df2 - df0*f2
            if 4*w01*w12 - w02*w02 < 0: bad = True
            if w02 != 0:
                r = (4*w01*w12)/(w02*w02)
                if r < min_J_fixed: min_J_fixed = r
    if bad: neg_J_fixed += 1
print(f"  Fixed a_j=1, w_j=1, random b_j in [1.1, 1001] (40 trials, k=8): J_i<0 fails={neg_J_fixed}/40, min 4W01W12/W02^2={mp.nstr(min_J_fixed, 6)}", flush=True)

print("\n=== 3. GLZ QUAD J^Q > 0 ON R AND PLUCKER RATIO W_AB W_CE / W_BC^2 FOR k=12..28 ===", flush=True)
ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 29):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_T(ht[k-1][idx], 1, k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_T(poly_add(p1, p2), 1, k + 1)

grid_R = [mp.mpf(v) for v in [-100, -10, -2, -1.1, -1, -0.99, -0.9, -0.7, -0.5, -0.3, -0.1, -0.01, 0, 0.01, 0.1, 1, 10, 100]]
for k in [12, 16, 20, 24, 28]:
    min_JQ_ratio = mp.mpf('1e100')
    max_plucker = mp.mpf(0)
    worst_jq = None
    worst_pl = None
    # Precompute D_{k+1}^r htilde_{k-1, i} for r = 0, 1, 2, 3
    D0 = ht[k-1]
    D1 = {idx: apply_T(D0[idx], 1, k+1) for idx in range(k)}
    D2 = {idx: apply_T(D1[idx], 1, k+1) for idx in range(k)}
    D3 = {idx: apply_T(D2[idx], 1, k+1) for idx in range(k)}
    dD0 = {idx: poly_deriv(D0[idx]) for idx in range(k)}
    dD1 = {idx: poly_deriv(D1[idx]) for idx in range(k)}
    dD2 = {idx: poly_deriv(D2[idx]) for idx in range(k)}
    dD3 = {idx: poly_deriv(D3[idx]) for idx in range(k)}
    # Check Quad (X=D0[i-2], M=D1[i-1], Y=D2[i]) for i=2..k-1
    for idx in range(2, k):
        for tv in grid_R:
            x, dx = poly_eval(D0[idx-2], tv), poly_eval(dD0[idx-2], tv)
            m, dm = poly_eval(D1[idx-1], tv), poly_eval(dD1[idx-1], tv)
            y, dy = poly_eval(D2[idx], tv),   poly_eval(dD2[idx], tv)
            wXM = x*dm - dx*m
            wMY = m*dy - dm*y
            wXY = x*dy - dx*y
            if wXY != 0:
                r = (4*wXM*wMY)/(wXY*wXY)
                if r < min_JQ_ratio:
                    min_JQ_ratio = r
                    worst_jq = (idx, float(tv))
    # Check Plucker (A=D0[i-3], B=D1[i-2], C=D2[i-1], E=D3[i]) for i=3..k-1
    for idx in range(3, k):
        for tv in grid_R:
            a, da = poly_eval(D0[idx-3], tv), poly_eval(dD0[idx-3], tv)
            b, db = poly_eval(D1[idx-2], tv), poly_eval(dD1[idx-2], tv)
            c, dc = poly_eval(D2[idx-1], tv), poly_eval(dD2[idx-1], tv)
            e, de = poly_eval(D3[idx], tv),   poly_eval(dD3[idx], tv)
            wAB = a*db - da*b
            wBC = b*dc - db*c
            wCE = c*de - dc*e
            r = (wAB*wCE)/(wBC*wBC)
            if r > max_plucker:
                max_plucker = r
                worst_pl = (idx, float(tv))
    print(f"  GLZ k={k:2d}: min 4*W_XM*W_MY/W_XY^2 = {mp.nstr(min_JQ_ratio, 7)} at {worst_jq}; max W_AB*W_CE/W_BC^2 = {mp.nstr(max_plucker, 7)} at {worst_pl}", flush=True)
