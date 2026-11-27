"""
Diagnose EXACTLY why Stage 3 (D_m -> D_{m+1}) drops gamma = 4 W_XM W_MY / W_XY^2 by only ~0.35/k^2 on GLZ,
while random triples in test_S23_and_cubic_stages.py showed larger drops (was it np.roots noise, d=0, or root clustering?).
Also compute exact d/ds ln(gamma(s))|_{s=0} for GLZ k=6..24!
"""
import math, random
import mpmath as mp
mp.mp.dps = 50

def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p
def poly_add(a, b):
    n = max(len(a), len(b))
    return poly_trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def poly_mul(a, b):
    res = [0]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b): res[i+j] += ai*bj
    return poly_trim(res)
def poly_deriv(a):
    return [i * a[i] for i in range(1, len(a))] if len(a) > 1 else [0]
def apply_Dj(p, j):
    dp = poly_deriv(p)
    res = [0] * (len(p) + 2)
    for i, c in enumerate(p):
        res[i] += c; res[i+1] += j * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c; res[i+2] += c
    return poly_trim(res)
def poly_eval(p, x):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc
def pfromr(rts):
    p = [mp.mpf(1)]
    for r in rts: p = poly_mul(p, [-r, mp.mpf(1)])
    return p

tgrid = [mp.mpf(v) for v in [-20, -5, -2, -1.1, -1, -0.95, -0.8, -0.5, -0.2, -0.05, 0, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 5, 20]]
def min_gamma_mp(X, M, Y):
    dX, dM, dY = poly_deriv(X), poly_deriv(M), poly_deriv(Y)
    mn = mp.mpf('1e100')
    wt = None
    for tv in tgrid:
        x, dx = poly_eval(X, tv), poly_eval(dX, tv)
        m, dm = poly_eval(M, tv), poly_eval(dM, tv)
        y, dy = poly_eval(Y, tv), poly_eval(dY, tv)
        w1 = x*dm - dx*m; w2 = m*dy - dm*y; w3 = x*dy - dx*y
        if w1 <= 0 or w2 <= 0:
            return mp.mpf(-1), float(tv)
        if w3 != 0:
            r = (4*w1*w2)/(w3*w3)
            if r < mn:
                mn = r; wt = float(tv)
    return mn, wt

print("=== 1. EXACT 50-DIGIT CHECK OF S3 AND S2+S3 ON RANDOM TRIPLES (j=6, 12, 24) ===", flush=True)
random.seed(111)
for j in [6, 12, 24]:
    tested = 0
    fails_S3 = 0
    fails_S23 = 0
    worst_S3 = (mp.mpf(0), None)
    worst_S23 = (mp.mpf(0), None)
    for tr in range(1000):
        if tested >= 40: break
        d = random.choice([0, 1, 2])
        r2 = sorted([-mp.mpf(random.uniform(0.05, 0.95)) for _ in range(d)])
        p2 = pfromr(r2)
        Dp2 = apply_Dj(p2, j)
        rDp2 = sorted([r.real for r in mp.polyroots(list(reversed(Dp2)), maxsteps=200, extraprec=20)])
        edges = [mp.mpf(-1)] + r2 + [mp.mpf(0)]
        r1 = [edges[m] + mp.mpf(random.uniform(0.02, 0.98))*(rDp2[m]-edges[m]) for m in range(d+1)]
        p1 = pfromr(r1)
        Dp1 = apply_Dj(p1, j)
        rDp1 = sorted([r.real for r in mp.polyroots(list(reversed(Dp1)), maxsteps=200, extraprec=20)])
        edges1 = [mp.mpf(-1)] + r1 + [mp.mpf(0)]
        r0 = [edges1[m] + mp.mpf(random.uniform(0.02, 0.98))*(rDp1[m]-edges1[m]) for m in range(d+2)]
        p0 = pfromr(r0)
        g_in, wt_in = min_gamma_mp(p0, apply_Dj(p1, j), apply_Dj(apply_Dj(p2, j), j))
        if g_in <= mp.mpf('1.001'): continue
        tested += 1
        g_s3, wt_s3 = min_gamma_mp(p0, apply_Dj(p1, j+1), apply_Dj(apply_Dj(p2, j+1), j+1))
        h0, h1, h2 = apply_Dj(p0, j), apply_Dj(p1, j), apply_Dj(p2, j)
        g_s2, wt_s2 = min_gamma_mp(h0, apply_Dj(h1, j), apply_Dj(apply_Dj(h2, j), j))
        g_s23, wt_s23 = min_gamma_mp(h0, apply_Dj(h1, j+1), apply_Dj(apply_Dj(h2, j+1), j+1))
        if g_s3 < 1: fails_S3 += 1
        if g_s23 < 1: fails_S23 += 1
        if g_in - g_s3 > worst_S3[0]: worst_S3 = (g_in - g_s3, (d, float(g_in), float(g_s3), wt_s3))
        if g_s2 - g_s23 > worst_S23[0]: worst_S23 = (g_s2 - g_s23, (d, float(g_s2), float(g_s23), wt_s23))
    print(f"  j={j:2d} ({tested} tested): S3_fails={fails_S3} (max_drop={mp.nstr(worst_S3[0],5)} {worst_S3[1]}); S2+S3_fails={fails_S23} (max S2->S23 drop={mp.nstr(worst_S23[0],5)} {worst_S23[1]})", flush=True)

print("\n=== 2. EXACT d/ds ln(gamma(s))|_{s=0} FOR GLZ AT STAGE 3 (k=6..20, ALL t & ALL i) ===", flush=True)
ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 21):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

# Compute exact derivative d/ds ln(4 W(X, M(s)) W(M(s), Y(s)) / W(X, Y(s))^2) at s=0
# where X = h0, M(s) = D_{m+s} h1 => M'(0) = t h1, Y(s) = D_{m+s}^2 h2 => Y'(0) = 2 t D_m h2 + t(1+t) h2
for k in [6, 10, 14, 18, 20]:
    m = k + 1
    min_dlog = mp.mpf('1e100')
    max_dlog = mp.mpf('-1e100')
    worst_min = None
    for i in range(2, k + 1):
        h0, h1, h2 = ht[k][i-2], ht[k][i-1], ht[k][i]
        M = apply_Dj(h1, m)
        Dh2 = apply_Dj(h2, m)
        Y = apply_Dj(Dh2, m)
        # Mdot = t * h1
        Mdot = [0] + h1
        # Ydot = 2 t Dh2 + (t + t^2) h2
        Ydot = poly_add([0] + [2*c for c in Dh2], poly_add([0] + h2, [0, 0] + h2))
        dX, dM, dY = poly_deriv(h0), poly_deriv(M), poly_deriv(Y)
        dMdot, dYdot = poly_deriv(Mdot), poly_deriv(Ydot)
        for tv in tgrid:
            x, dx = poly_eval(h0, tv), poly_eval(dX, tv)
            mv, dm = poly_eval(M, tv), poly_eval(dM, tv)
            yv, dy = poly_eval(Y, tv), poly_eval(dY, tv)
            md, dmd = poly_eval(Mdot, tv), poly_eval(dMdot, tv)
            yd, dyd = poly_eval(Ydot, tv), poly_eval(dYdot, tv)
            wXM = x*dm - dx*mv
            wMY = mv*dy - dm*yv
            wXY = x*dy - dx*yv
            wXM_dot = x*dmd - dx*md
            wMY_dot = (md*dy - dmd*yv) + (mv*dyd - dm*yd)
            wXY_dot = x*dyd - dx*yd
            dlog = wXM_dot/wXM + wMY_dot/wMY - 2*wXY_dot/wXY
            if dlog < min_dlog:
                min_dlog = dlog
                worst_min = (i, float(tv), float(wXM_dot/wXM), float(wMY_dot/wMY), float(wXY_dot/wXY))
            if dlog > max_dlog:
                max_dlog = dlog
    print(f"  GLZ k={k:2d}: min d(ln gamma)/ds = {mp.nstr(min_dlog, 6)} (k^2*min = {mp.nstr(k*k*min_dlog, 5)}) at i={worst_min[0]}, t={worst_min[1]}; max = {mp.nstr(max_dlog, 6)}", flush=True)
    print(f"     components at worst: dW_XM/W_XM={worst_min[2]:.5f}, dW_MY/W_MY={worst_min[3]:.5f}, dW_XY/W_XY={worst_min[4]:.5f}", flush=True)
