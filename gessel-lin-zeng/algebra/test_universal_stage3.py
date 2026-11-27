"""
CRITICAL TEST OF THE CONJECTURED UNIVERSAL STAGE 3 THEOREM:
For ANY h2 << h1 << h0 in (-1, 0) of degrees d-2, d-1, d (d >= 2) such that
Q_j = (h0, D_j h1, D_j^2 h2) is real-stable (gamma_j(t) >= 1 on R):
What is the WORST POSSIBLE value of:
  (1) (d - 1 + j)^2 * (1 - min_R(gamma_{j+1}) / min_R(gamma_j))
  (2) (d - 1 + j)^2 * max_R( -d/dj ln gamma_j(t) )
  (3) j^2 * (1 - min_R(gamma_{j+1}) / min_R(gamma_j))
across adversarial random triples (including roots VERY close to 0, -1, and each other!) for d=2..6, j=4..30?!
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

tgrid = [mp.mpf(v) for v in [-100, -10, -3, -1.5, -1.1, -1.01, -1, -0.99, -0.9, -0.7, -0.5, -0.3, -0.1, -0.02, -0.002,
                              0, 0.002, 0.02, 0.1, 0.3, 0.5, 1, 3, 10, 100]]

def analyze_triple(h0, h1, h2, j, d):
    M = apply_Dj(h1, j)
    Dh2 = apply_Dj(h2, j)
    Y = apply_Dj(Dh2, j)
    Mdot = [0] + h1
    Ydot = poly_add([0] + [2*c for c in Dh2], poly_add([0] + h2, [0, 0] + h2))
    dX, dM, dY = poly_deriv(h0), poly_deriv(M), poly_deriv(Y)
    dMdot, dYdot = poly_deriv(Mdot), poly_deriv(Ydot)
    min_g = mp.mpf('1e100')
    min_dlog = mp.mpf('1e100')
    wt_g = None; wt_d = None
    for tv in tgrid:
        x, dx = poly_eval(h0, tv), poly_eval(dX, tv)
        mv, dm = poly_eval(M, tv), poly_eval(dM, tv)
        yv, dy = poly_eval(Y, tv), poly_eval(dY, tv)
        wXM = x*dm - dx*mv; wMY = mv*dy - dm*yv; wXY = x*dy - dx*yv
        if wXM <= 0 or wMY <= 0:
            return mp.mpf(-1), mp.mpf(-1e100), None, None
        g = (4*wXM*wMY)/(wXY*wXY) if wXY != 0 else mp.mpf('1e100')
        if g < min_g: min_g = g; wt_g = float(tv)
        md, dmd = poly_eval(Mdot, tv), poly_eval(dMdot, tv)
        yd, dyd = poly_eval(Ydot, tv), poly_eval(dYdot, tv)
        wXM_d = x*dmd - dx*md
        wMY_d = (md*dy - dmd*yv) + (mv*dyd - dm*yd)
        wXY_d = x*dyd - dx*yd
        dlog = wXM_d/wXM + wMY_d/wMY - 2*wXY_d/wXY
        if dlog < min_dlog: min_dlog = dlog; wt_d = float(tv)
    return min_g, min_dlog, wt_g, wt_d

print("=== ADVERSARIAL RANDOM TRIPLES: CHECK (d-1+j)^2 * (-d/dj ln gamma) ===", flush=True)
random.seed(2026)
for (d, j) in [(2, 6), (2, 12), (2, 24), (3, 12), (4, 12), (5, 12), (4, 24)]:
    tested = 0
    max_scaled_dlog = mp.mpf(0)
    max_scaled_drop = mp.mpf(0)
    worst_d = None; worst_dr = None
    for tr in range(1500):
        if tested >= 35: break
        # Include adversarial root distributions (some clustered near 0, some near -1, some uniform)
        mode = tr % 3
        if mode == 0:
            r2 = sorted([-mp.mpf(random.uniform(0.01, 0.99)) for _ in range(d - 2)])
        elif mode == 1:
            r2 = sorted([-mp.mpf(10**random.uniform(-3, -0.02)) for _ in range(d - 2)])
        else:
            r2 = sorted([-mp.mpf(1) + mp.mpf(10**random.uniform(-3, -0.02)) for _ in range(d - 2)])
        h2 = pfromr(r2)
        Dh2 = apply_Dj(h2, j)
        rDh2 = sorted([r.real for r in mp.polyroots(list(reversed(Dh2)), maxsteps=200, extraprec=20)])
        edges = [mp.mpf(-1)] + r2 + [mp.mpf(0)]
        r1 = [edges[m] + mp.mpf(random.uniform(0.01, 0.99))*(rDh2[m]-edges[m]) for m in range(d - 1)]
        h1 = pfromr(r1)
        Dh1 = apply_Dj(h1, j)
        rDh1 = sorted([r.real for r in mp.polyroots(list(reversed(Dh1)), maxsteps=200, extraprec=20)])
        edges1 = [mp.mpf(-1)] + r1 + [mp.mpf(0)]
        r0 = [edges1[m] + mp.mpf(random.uniform(0.01, 0.99))*(rDh1[m]-edges1[m]) for m in range(d)]
        h0 = pfromr(r0)
        g_j, dlog, wt_g, wt_d = analyze_triple(h0, h1, h2, j, d)
        if g_j <= mp.mpf('1.0001'): continue
        tested += 1
        g_jp1, _, wt_g2, _ = analyze_triple(h0, h1, h2, j + 1, d)
        scale = (d - 1 + j)**2
        sc_dlog = -scale * dlog
        sc_drop = scale * (1 - g_jp1 / g_j)
        if sc_dlog > max_scaled_dlog:
            max_scaled_dlog = sc_dlog
            worst_d = (mode, float(g_j), wt_d, float(-dlog))
        if sc_drop > max_scaled_drop:
            max_scaled_drop = sc_drop
            worst_dr = (mode, float(g_j), float(g_jp1), wt_g2)
    print(f"  d={d}, j={j:2d} (N={d-1+j:2d}, {tested:2d} tested): max N^2*(-d ln g/dj) = {mp.nstr(max_scaled_dlog, 5)} {worst_d}; max N^2*(1-g_{{j+1}}/g_j) = {mp.nstr(max_scaled_drop, 5)}", flush=True)
