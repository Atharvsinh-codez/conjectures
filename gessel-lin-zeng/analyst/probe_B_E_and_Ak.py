"""FAST DECISIVE PROBE:
1. Fixed-t real-rootedness in y (for t in (-1, 0) / t < 0) of:
   - E_N(t, y) = sum (1/i!) F_{N,i}(t) y^i
   - Etilde_N(t, y) = sum (1/(N-i)!) F_{N,i}(t) y^i
   - B_N(t, y) = sum C(N,i) F_{N,i}(t) y^i
   Tested on:
   (a) GLZ k = 2..20 (mpmath at dense t in (-1,0) including -1e-12 and -0.99999!)
   (b) EXTREME wild Atomic D_{c,m} (m_j jumping by 1e6, bimodal lam_j!)
   (c) EXTREME wild Laguerre L_{c,m} (m_j jumping by 1e6, bimodal lam_j!)
2. Check B_N(t, y) and E_N(t, y) at t > 0 (contrast!).
3. ADVERSARIAL ATTACK ON ALIGNED BINOMIAL A_k(t, y) = sum C(k,i) D_{k+1}^i htilde_{k,i}(t) y^i:
   - Check Rayleigh R(t, y) = A_t A_y - A A_{ty} on ALL (t, y) in R^2 for k = 2..15 (especially t > 0 AND t in (-1, 0) AND near t = -1)!
   - Check A_N(t, y) = sum C(N,i) D_m^i F_{N,i}(t) y^i on WILD Atomic parameters!
"""
import math
import numpy as np
import mpmath as mp

def apply_D_int(m, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + 1)
        res[r + 1] += coef * (r + m)
    return res

def poly_add_int(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def eval_mp(p, tv):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * tv + c
    return acc

# Build GLZ H[k][i] up to K=18
KMAX = 18
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D_int(m, H[k-1][0])
    for i in range(1, k):
        U[i] = poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1):
        H[k][i] = apply_D_int(m, U[i])

print("=== 1. GLZ FIXED-t REAL-ROOTEDNESS IN y FOR E_k, Etilde_k, B_k (k=2..18, mpmath) ===")
mp.mp.dps = 50
ts_neg = [mp.mpf(s) for s in ["-1e-12", "-1e-7", "-1e-4", "-0.01", "-0.1", "-0.3", "-0.5", "-0.8", "-0.99", "-0.9999"]]
for name, wfn in [
    ("E_k (1/i!)", lambda k, i: mp.mpf(1) / mp.factorial(i)),
    ("Etilde_k (1/(k-i)!)", lambda k, i: mp.mpf(1) / mp.factorial(k - i)),
    ("B_k (C(k,i))", lambda k, i: mp.mpf(math.comb(k, i))),
]:
    fails = []
    for k in range(2, 19):
        for tv in ts_neg:
            coeffs = [wfn(k, i) * eval_mp(H[k][i], tv) for i in range(k+1)]
            try:
                rts = mp.polyroots(coeffs[::-1], maxsteps=200, extraprec=30)
            except Exception:
                rts = mp.polyroots(coeffs[::-1], maxsteps=500, extraprec=60)
            mi = max(abs(mp.im(r)) for r in rts)
            if mi > mp.mpf('1e-20'):
                fails.append((k, mp.nstr(tv, 4), mp.nstr(mi, 4)))
                break
    print(f"  {name:22s}: fails on (-1,0) for k=2..18 = {fails[:5]}", flush=True)

print("\n=== 2. CONTRAST AT t > 0 (t = +0.01, +0.5) FOR GLZ k=2..10 ===")
for name, wfn in [
    ("E_k (1/i!)", lambda k, i: mp.mpf(1) / mp.factorial(i)),
    ("Etilde_k (1/(k-i)!)", lambda k, i: mp.mpf(1) / mp.factorial(k - i)),
    ("B_k (C(k,i))", lambda k, i: mp.mpf(math.comb(k, i))),
]:
    fails = []
    for k in range(2, 11):
        for tv in [mp.mpf("0.01"), mp.mpf("0.5")]:
            coeffs = [wfn(k, i) * eval_mp(H[k][i], tv) for i in range(k+1)]
            rts = mp.polyroots(coeffs[::-1], maxsteps=200, extraprec=30)
            mi = max(abs(mp.im(r)) for r in rts)
            if mi > mp.mpf('1e-20'):
                fails.append((k, mp.nstr(tv, 3), mp.nstr(mi, 4)))
                break
    print(f"  {name:22s}: fails at t>0 (k=2..10) = {fails[:5]}", flush=True)

print("\n=== 3. EXTREME WILD ATOMIC & LAGUERRE TEST FOR E_N, Etilde_N, B_N ===")
def apply_D_mp(c, m, poly):
    d = len(poly) - 1
    res = [mp.mpf(0)] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * (r + m)
    return res

def apply_L_mp(c, m, poly):
    d = len(poly) - 1
    res = [mp.mpf(0)] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * m
    return res

def poly_add_mp(a, b, s):
    n = max(len(a), len(b))
    res = [mp.mpf(0)] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += s * v
    return res

mp.mp.dps = 50
rng = np.random.default_rng(2026)
for kind, op_fn, ts_list in [
    ("EXTREME Atomic D", apply_D_mp, [mp.mpf(s) for s in ["-1e-9", "-1e-5", "-1e-2", "-0.2", "-0.5", "-0.9", "-0.999"]]),
    ("EXTREME Laguerre L", apply_L_mp, [mp.mpf(s) for s in ["-1e-9", "-1e-4", "-1e-1", "-1", "-10", "-1000"]]),
]:
    res_counts = {"E": 0, "Etilde": 0, "B": 0, "F_plain": 0}
    for trial in range(25):
        N = int(rng.integers(3, 8))
        polys = [[mp.mpf(1)]]
        for j in range(N):
            c = mp.mpf(10 ** rng.uniform(-2, 2))
            m = (c + mp.mpf(10 ** rng.uniform(-2, 4))) if "Atomic" in kind else mp.mpf(10 ** rng.uniform(-2, 4))
            lam = mp.mpf(10 ** rng.uniform(-2, 2))
            new_p = [[mp.mpf(0)]] * (len(polys) + 1)
            for i, p in enumerate(polys):
                new_p[i] = poly_add_mp(new_p[i], op_fn(c, m, p), mp.mpf(1))
                new_p[i+1] = poly_add_mp(new_p[i+1], p, lam)
            polys = new_p
        failed = {k: False for k in res_counts}
        for tv in ts_list:
            vals = [eval_mp(p, tv) for p in polys]
            for tag, wfn in [
                ("F_plain", lambda i: mp.mpf(1)),
                ("E", lambda i: mp.mpf(1)/mp.factorial(i)),
                ("Etilde", lambda i: mp.mpf(1)/mp.factorial(N-i)),
                ("B", lambda i: mp.mpf(math.comb(N, i))),
            ]:
                if failed[tag]: continue
                coeffs = [wfn(i) * vals[i] for i in range(N+1)]
                rts = mp.polyroots(coeffs[::-1], maxsteps=200, extraprec=30)
                if max(abs(mp.im(r)) for r in rts) > mp.mpf('1e-18'):
                    failed[tag] = True
        for tag in res_counts:
            if failed[tag]: res_counts[tag] += 1
    print(f"  {kind:20s} (25 extreme trials): {res_counts}", flush=True)

print("\n=== 4. ADVERSARIAL RAYLEIGH CHECK FOR ALIGNED BINOMIAL A_k(t, y) ON ALL (t, y) IN R^2 (k=2..14) ===")
# A_k(t, y) = sum_{i=0}^k binom(k, i) g_{k,i}(t) y^i where g_{k,i} = D_{k+1}^i htilde_{k,i}
# Check sign of A_t A_y - A A_{ty} across (t, y) in R^2!
ts_ray = [mp.mpf(s) for s in ["-10", "-2", "-1.1", "-1", "-0.99", "-0.9", "-0.5", "-0.1", "-0.01", "-1e-5", "0", "1e-5", "0.01", "0.1", "0.5", "1", "2", "10"]]
ys_ray = [mp.mpf(s) for s in ["-100", "-10", "-3", "-1.5", "-1.05", "-1", "-0.95", "-0.75", "-0.5", "-0.2", "-0.05", "-1e-4", "0", "0.1", "1", "10"]]

def deriv_int(p):
    return [r * p[r] for r in range(1, len(p))]

for k in range(2, 15):
    m = k + 1
    gs = []
    gps = []
    for i in range(k+1):
        p = list(H[k][i])
        for _ in range(i): p = apply_D_int(m, p)
        gs.append(p)
        gps.append(deriv_int(p))
    pos_cnt = 0; neg_cnt = 0
    bad_pt = None
    for tv in ts_ray:
        gv = [mp.mpf(math.comb(k, i)) * eval_mp(gs[i], tv) for i in range(k+1)]
        gpv = [mp.mpf(math.comb(k, i)) * eval_mp(gps[i], tv) for i in range(k+1)]
        for yv in ys_ray:
            # Evaluate A, A_t, A_y, A_{ty}
            A = mp.mpf(0); At = mp.mpf(0); Ay = mp.mpf(0); Aty = mp.mpf(0)
            yp = mp.mpf(1)
            for i in range(k+1):
                A += gv[i] * yp
                At += gpv[i] * yp
                if i > 0:
                    yp_prev = yp / yv if yv != 0 else (mp.mpf(1) if i == 1 else mp.mpf(0))
                    Ay += i * gv[i] * yp_prev
                    Aty += i * gpv[i] * yp_prev
                yp *= yv
            ray = At * Ay - A * Aty
            if ray > mp.mpf('1e-30') * (abs(At*Ay) + abs(A*Aty) + 1):
                pos_cnt += 1
            elif ray < -mp.mpf('1e-30') * (abs(At*Ay) + abs(A*Aty) + 1):
                neg_cnt += 1
                if bad_pt is None:
                    bad_pt = (mp.nstr(tv, 4), mp.nstr(yv, 4), mp.nstr(ray, 4))
    print(f"  A_{k:2d}(t,y) Rayleigh on R^2: pos={pos_cnt}, neg={neg_cnt}, first_neg={bad_pt}", flush=True)
