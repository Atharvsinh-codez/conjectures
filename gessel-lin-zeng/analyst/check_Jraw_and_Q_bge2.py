import mpmath as mp

mp.mp.dps = 80

def apply_D(b, p):
    res = [mp.mpf(0)]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+b)
    return res

def add_mp(a, b):
    res = [mp.mpf(0)]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv_mp(p):
    return [r*p[r] for r in range(1, len(p))] if len(p) > 1 else [mp.mpf(0)]

def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v

def eval_W(f, g, t):
    # W(f, g) = f g' - f' g; for deg f = deg g + 1 with g << f, W(g, f) > 0!
    return eval_mp(f, t)*eval_mp(deriv_mp(g), t) - eval_mp(deriv_mp(f), t)*eval_mp(g, t)

ts = [mp.mpf(s) for s in [
    "-100", "-10", "-2", "-1.1", "-1.01", "-1.0001", "-1", "-0.9999", "-0.99", "-0.95",
    "-0.9", "-0.7", "-0.5", "-0.3", "-0.1", "-0.01", "-1e-4", "-1e-8", "-1e-14", "0",
    "1e-8", "0.01", "0.1", "1", "10", "100"
]]

def check_Jraw_and_Q_rayleigh(name, bs, dbl, initP, qs=(0.5, 0.25)):
    P = initP
    for b in bs:
        k = len(P)
        newP = []
        for i in range(k + 1):
            t2 = apply_D(b, P[i]) if i < k else [mp.mpf(0)]
            t1 = P[i-1] if i >= 1 else [mp.mpf(0)]
            u = add_mp(t2, t1)
            newP.append(apply_D(b, u) if dbl else u)
        P = newP
    # Check min J_raw(P[i+1], P[i], P[i-1]) on ts (note W(P[i], P[i-1]) > 0!)
    mn_J = mp.mpf('1e100'); loc_J = None
    for i in range(1, len(P) - 1):
        for tv in ts:
            w1 = eval_W(P[i], P[i-1], tv)
            w2 = eval_W(P[i+1], P[i], tv)
            w12 = eval_W(P[i+1], P[i-1], tv)
            if w1 <= 0 or w2 <= 0:
                mn_J = mp.mpf(-1); loc_J = (i, float(tv), "W<=0"); break
            if w12 != 0:
                J = 4*w1*w2/(w12*w12)
                if J < mn_J: mn_J = J; loc_J = (i, float(tv))
        if mn_J < 0: break
    # Also check Rayleigh Q_t Q_y - Q Q_{ty} for Q(t, y) = sum q^{i(i-1)/2} P[i](t) y^i at a grid of (t, y)!
    # Wait: check sign of Q_t Q_y - Q Q_{ty} for real-stable Q(t, y) with positive coeffs:
    # For (1+2t) + y: Q_t = 2, Q_y = 1, Q_{ty} = 0 => Q_t Q_y - Q Q_{ty} = +2 > 0!
    # Wait: for deg P[i] decreasing in i, at t -> +inf, P[0] ~ c0 t^d, P[1] ~ c1 t^{d-1},
    # Q_t Q_y - Q Q_{ty} = sum_{i,j} (j - i) ... Wait: for real-stable Q(t, y) in H+ x H+,
    # Delta(t, y) := Q_t Q_y - Q Q_{ty} >= 0 EVERYWHERE on R^2 iff Q is real-stable (Wait: that's Rayleigh, NECESSARY for real stability;
    # check sign for (1+2t) + y: Delta = 2*1 - 0 = +2 > 0! Check for (1+2t)(1+3t) + (1+2.5t)y:
    # Q = (1+5t+6t^2) + (1+2.5t)y => Q_t Q_y - Q Q_{ty} = (5+12t)(1+2.5t) - (1+5t+6t^2)(2.5) = 2.5 + 12t + 15t^2 > 0!).
    # YES! Delta = Q_t Q_y - Q Q_{ty} >= 0 on ALL R^2 is NECESSARY for real stability!
    ys_probe = [mp.mpf(s) for s in ["-1e6", "-1e3", "-10", "-1", "-0.1", "-1e-4", "0", "1e-4", "0.1", "1", "10", "1e3", "1e6"]]
    ray_res = []
    for q in qs:
        qm = mp.mpf(str(q))
        min_rel = mp.mpf('1e100'); bad_loc = None
        # precompute P[i](t) and P'[i](t) on ts
        dP = [deriv_mp(p) for p in P]
        for tv in [mp.mpf(s) for s in ["-10", "-1.1", "-1", "-0.99", "-0.5", "-0.01", "-1e-6", "0", "0.1", "10"]]:
            vals = [qm**(i*(i-1)//2) * eval_mp(P[i], tv) for i in range(len(P))]
            dvals = [qm**(i*(i-1)//2) * eval_mp(dP[i], tv) for i in range(len(P))]
            for yv in ys_probe:
                Q = mp.mpf(0); Qt = mp.mpf(0); Qy = mp.mpf(0); Qty = mp.mpf(0)
                yp = mp.mpf(1)
                for i in range(len(P)):
                    Q += vals[i] * yp
                    Qt += dvals[i] * yp
                    if i >= 1:
                        yp_prev = yp / yv if yv != 0 else (mp.mpf(1) if i == 1 else mp.mpf(0))
                        Qy += i * vals[i] * yp_prev
                        Qty += i * dvals[i] * yp_prev
                    yp *= yv
                delta = Qt * Qy - Q * Qty
                scale = abs(Qt * Qy) + abs(Q * Qty) + mp.mpf('1e-300')
                rel = delta / scale
                if rel < min_rel:
                    min_rel = rel; bad_loc = (float(tv), float(yv))
        ray_res.append(f"q={q}:minRay={float(min_rel):+.4f}@{bad_loc if min_rel < 0 else 'OK'}")
    print(f"  {name:30s}: min_Jraw={float(mn_J):+8.5f}@{loc_J} | " + " | ".join(ray_res))

print("=== FIXED-SIGN J_raw AND RAYLEIGH CHECK FOR b_j >= 2 (UP TO k=36!) ===")
tests = [
    ("GLZ k=12 (dbl)",             [mp.mpf(j+2) for j in range(1, 12)], True,  [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]),
    ("GLZ k=24 (dbl)",             [mp.mpf(j+2) for j in range(1, 24)], True,  [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]),
    ("GLZ k=36 (dbl)",             [mp.mpf(j+2) for j in range(1, 36)], True,  [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]),
    ("b=[2]*20 (sgl)",             [mp.mpf(2)]*20, False, [[mp.mpf(1)]]),
    ("b=[2]*20 (dbl)",             [mp.mpf(2)]*20, True,  [[mp.mpf(1)]]),
    ("b=[2]*18+[1e8]*18 (sgl)",    [mp.mpf(2)]*18 + [mp.mpf('1e8')]*18, False, [[mp.mpf(1)]]),
    ("b=[2]*18+[1e8]*18 (dbl)",    [mp.mpf(2)]*18 + [mp.mpf('1e8')]*18, True,  [[mp.mpf(1)]]),
    ("b=[1e8]*18+[2]*18 (sgl)",    [mp.mpf('1e8')]*18 + [mp.mpf(2)]*18, False, [[mp.mpf(1)]]),
    ("b=[1e8]*18+[2]*18 (dbl)",    [mp.mpf('1e8')]*18 + [mp.mpf(2)]*18, True,  [[mp.mpf(1)]]),
    ("b=[2,1e8]*12 (sgl)",         [mp.mpf(2), mp.mpf('1e8')]*12, False, [[mp.mpf(1)]]),
    ("b=[1.5]*10+[1e6]*10 (sgl)",  [mp.mpf('1.5')]*10 + [mp.mpf('1e6')]*10, False, [[mp.mpf(1)]]),
    ("b=[1.01]*5+[1e6]*5 (sgl)",   [mp.mpf('1.01')]*5 + [mp.mpf('1e6')]*5, False, [[mp.mpf(1)]]),
]
for args in tests:
    check_Jraw_and_Q_rayleigh(*args)
