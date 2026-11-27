import numpy as np
import mpmath as mp

mp.mp.dps = 60

def apply_D_int(m, p):
    res = [0]*(len(p) + 1)
    for r, c in enumerate(p):
        res[r] += c*(r+1)
        res[r+1] += c*(r+m)
    return res

def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv_int(p):
    return [r*p[r] for r in range(1, len(p))] if len(p) > 1 else [0]

def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v

def eval_W_mp(f, g, t):
    return eval_mp(f, t)*eval_mp(deriv_int(g), t) - eval_mp(deriv_int(f), t)*eval_mp(g, t)

# Build GLZ tilde_h up to k = 24
KMAX = 24
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    H[k][0] = apply_D_int(m, apply_D_int(m, H[k-1][0]))
    for i in range(1, k):
        H[k][i] = apply_D_int(m, poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1]))
    H[k][k] = apply_D_int(m, H[k-1][k-1])

print("=== 1. CHECK W(E_{k,i}, E_{k,i+r}) > 0 FOR r = 1..10 UP TO k = 24 ===")
# Check at critical t points including -0.99, -0.95, -0.9, -0.8, -0.5, -0.1, -1e-4, 0, 1, -2!
ts_chk = [mp.mpf(s) for s in ["-10", "-1.5", "-0.999", "-0.99", "-0.95", "-0.9", "-0.8", "-0.5", "-0.1", "-0.01", "-1e-4", "0", "0.1", "10"]]
for k in [10, 14, 18, 22, 24]:
    m = k + 1
    E = []
    for i in range(k):
        p = H[k-1][i]
        for _ in range(i):
            p = apply_D_int(m, p)
        E.append(p)
    r_status = []
    for r in range(1, min(11, k)):
        fails = 0
        first = None
        for i in range(k - r):
            for tv in ts_chk:
                if eval_W_mp(E[i], E[i+r], tv) <= 0:
                    fails += 1
                    if first is None: first = (i, float(tv))
                    break
        r_status.append(f"r={r}:{fails}" + (f"@{first}" if first else ""))
    print(f"  k={k:2d}: " + " | ".join(r_status))

print("\n=== 2. SIGN OF d/dj log J(j; t) ON (-1, 0) vs OUTSIDE (-1, 0) ON GLZ (k=5..16) ===")
def poly_mul_t(p, q):
    res = [0]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            res[i+j] += a*b
    return res

def shift_t(p): return [0] + p

ts_in = [mp.mpf(s) for s in ["-0.999", "-0.99", "-0.95", "-0.9", "-0.8", "-0.7", "-0.5", "-0.3", "-0.1", "-0.05", "-0.01", "-1e-4", "-1e-7"]]
ts_out = [mp.mpf(s) for s in ["-100", "-10", "-2", "-1.1", "-1.01", "0", "1e-4", "0.01", "0.1", "0.5", "1", "2", "10", "100"]]

for k in [6, 10, 14, 16]:
    j0 = mp.mpf(k + 1)
    min_in = (mp.mpf('1e100'), None)
    min_out = (mp.mpf('1e100'), None)
    for i in range(1, k):
        h0, h1, h2 = H[k][i-1], H[k][i], H[k][i+1]
        D0h1 = apply_D_int(0, h1); th1 = shift_t(h1)
        D0h2 = apply_D_int(0, h2)
        Y0 = apply_D_int(0, D0h2)
        Y1 = poly_add_int([2*c for c in shift_t(D0h2)], poly_mul_t([0, 1, 1], h2))
        Y2 = shift_t(shift_t(h2))
        for tag, tlist in [("in", ts_in), ("out", ts_out)]:
            for tv in tlist:
                a0 = eval_W_mp(h0, D0h1, tv); a1 = eval_W_mp(h0, th1, tv)
                b0 = eval_W_mp(h0, Y0, tv); b1 = eval_W_mp(h0, Y1, tv); b2 = eval_W_mp(h0, Y2, tv)
                c0 = eval_W_mp(D0h1, Y0, tv)
                c1 = eval_W_mp(D0h1, Y1, tv) + eval_W_mp(th1, Y0, tv)
                c2 = eval_W_mp(D0h1, Y2, tv) + eval_W_mp(th1, Y1, tv)
                c3 = eval_W_mp(th1, Y2, tv)
                w1 = a0 + a1*j0; w1p = a1
                w2 = c0 + c1*j0 + c2*j0**2 + c3*j0**3; w2p = c1 + 2*c2*j0 + 3*c3*j0**2
                w12 = b0 + b1*j0 + b2*j0**2; w12p = b1 + 2*b2*j0
                dlogJ = w1p/w1 + w2p/w2 - 2*w12p/w12
                if tag == "in" and dlogJ < min_in[0]:
                    min_in = (dlogJ, (i, float(tv)))
                if tag == "out" and dlogJ < min_out[0]:
                    min_out = (dlogJ, (i, float(tv)))
    print(f"  k={k:2d} (j={k+1:2d}): (-1,0) min j^2*dlogJ = {float((k+1)**2*min_in[0]):+8.4f} @{min_in[1]} | outside min j^2*dlogJ = {float((k+1)**2*min_out[0]):+8.4f} @{min_out[1]}")

print("\n=== 3. FULL BIVARIATE REAL-STABILITY CHECK OF Q_{k,q} FOR b_j >= 2 (CHECK FIXED y < 0 TOO!) ===")
# Recall a polynomial Q(t, y) is real-stable iff:
#   - For every real t, Q(t, y) has all real y-roots, AND
#   - For every real y, Q(t, y) has all real t-roots! (Actually just one of these + sign of Wronskian is necessary;
#     checking fixed y in {-100, -10, -1, -0.5, -0.1, -0.01, -1e-4, 0.01, 1, 100} catches complex t-roots!)
def apply_D_mp(m, p):
    res = [mp.mpf(0)]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def add_mp(a, b):
    res = [mp.mpf(0)]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

for label, b_list, doubled in [
    ("[1.01]*4+[1e4]*4 sgl", [mp.mpf('1.01')]*4 + [mp.mpf('1e4')]*4, False),
    ("[2]*4+[1e4]*4 sgl",    [mp.mpf('2')]*4 + [mp.mpf('1e4')]*4, False),
    ("[2]*4+[1e4]*4 dbl",    [mp.mpf('2')]*4 + [mp.mpf('1e4')]*4, True),
    ("[2]*6+[1e6]*6 dbl",    [mp.mpf('2')]*6 + [mp.mpf('1e6')]*6, True),
    ("[2]*3+[1e8]*3 sgl",    [mp.mpf('2')]*3 + [mp.mpf('1e8')]*3, False),
    ("[2]*3+[1e8]*3 dbl",    [mp.mpf('2')]*3 + [mp.mpf('1e8')]*3, True),
    ("[1e6]*4+[2]*4 dbl",    [mp.mpf('1e6')]*4 + [mp.mpf('2')]*4, True),
]:
    P = [[mp.mpf(1)]]
    for b in b_list:
        k = len(P)
        newP = []
        for i in range(k + 1):
            t2 = apply_D_mp(b, P[i]) if i < k else [mp.mpf(0)]
            t1 = P[i-1] if i >= 1 else [mp.mpf(0)]
            u = add_mp(t2, t1)
            newP.append(apply_D_mp(b, u) if doubled else u)
        P = newP
    k = len(b_list)
    # Check both fixed y (real t-roots) and fixed t (real y-roots) for q = 0.5, 0.25, 0.125!
    res = []
    for q in [0.5, 0.25, 0.125]:
        qm = mp.mpf(str(q))
        fail_y = None
        for yv in [mp.mpf(s) for s in ["-1e6", "-1e4", "-100", "-10", "-2", "-1", "-0.5", "-0.1", "-0.01", "-1e-4", "1", "100"]]:
            # poly in t: sum_i q^{i(i-1)/2} y^i P[i](t)
            maxd = max(len(p) for p in P)
            pt = [mp.mpf(0)]*maxd
            for i in range(k+1):
                wt = qm**(i*(i-1)//2) * (yv**i)
                for r, c in enumerate(P[i]):
                    pt[r] += wt * c
            while len(pt) > 1 and pt[-1] == 0: pt.pop()
            if len(pt) <= 2: continue
            try:
                roots = mp.polyroots(pt[::-1], maxsteps=200, extraprec=50)
                mx = float(max(abs(mp.im(z)) for z in roots))
                if mx > 1e-12:
                    fail_y = (float(yv), mx)
                    break
            except Exception:
                pass
        res.append(f"q={q}:{'FAIL@y='+str(fail_y[0]) if fail_y else 'PASS'}")
    print(f"  {label:24s}: " + " | ".join(res))
