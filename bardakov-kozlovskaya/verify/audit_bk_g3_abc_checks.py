"""
Exact counterexample & verification script for analyst/certify_G3_ABC_full_iv.py:
  1. Test g_{1,fin}(u0) in Case G3-C at TRUE tau ~ 1.04 (<= 1.05) and increasing ct = 40, 1e4, 1e6, 1e8!
  2. Test c_rein, c_re5, c_re50 at TRUE VALID Piece G3 points (kap ~ 0.448 <= 0.45)!
"""
import mpmath as mp
mp.mp.dps = 35
pi = mp.pi

def make_point(n, c, j, r0, p0, D):
    assert r0 * j == c * p0 + D
    thc = pi * c / n; psij = pi * j / n; al1 = pi * (j + c) / n
    r1, r2, r3 = mp.sin(al1), mp.sin(psij), mp.sin(thc)
    sj = r2 / r1
    kap = 2 * pi * r1 * r2 / r3
    ct = mp.mpf(c) / r0
    tau = c * r0 * kap
    B = 4 * pi * c * r0 * thc
    te_eff = B * sj / (2 * tau * (1 - sj) ** 2)
    def zeta(u):
        w = sj * mp.exp(2j * u * thc)
        return r0 * n * mp.log((r3 / r1) / (1 - w)) + 1j * pi * r0 * j - 2j * pi * c * p0 * (u + mp.mpf('0.5'))
    def d1(u):
        w = sj * mp.exp(2j * u * thc)
        return 2j * pi * c * ((r0 + p0) * w - p0) / (1 - w)
    return dict(n=n, c=c, j=j, r0=r0, p0=p0, D=D, thc=thc, psij=psij, al1=al1,
                sj=sj, kap=kap, ct=ct, tau=tau, B=B, te_eff=te_eff, zeta=zeta, d1=d1)

print("=== TEST 1: Case G3-C (D=+1, TRUE tau ~ 1.04 <= 1.05, ct*tau^2 >= 40) ===")
# kap ~ 2*pi^2 * j*(j+c)/(c*n) => tau = c*kap ~ 2*pi^2 * j*(j+c)/n => n ~ 2*pi^2 * j*(j+c) / 1.04!
c_u = mp.mpf('1.22')
for c in [40, 1000, 100000, 10000000]:
    r0 = 1; p0 = 1; D = 1; j = c + 1
    n = 2 * int(round(pi ** 2 * j * (j + c) / 1.04))
    pt = make_point(n, c, j, r0, p0, D)
    ct = pt['ct']; tau = pt['tau']
    assert tau <= mp.mpf('1.05') and ct * (tau ** 2) >= 40 and c * (pt['kap'] ** 2) < 42
    X = ct * tau / (2 * pi ** 2)
    u0 = min(mp.mpf('0.42'), c_u / ((ct * (tau ** 3)) ** mp.mpf('0.25')))
    g1_par_0 = -sum(mp.sqrt(pi / (m * tau)) * mp.exp(-pi ** 2 * m / tau) for m in range(1, 6))
    Lam = mp.exp(-mp.re(pt['zeta'](u0)))
    def f_int(u):
        z = Lam * mp.exp(pt['zeta'](u))
        return 2 * mp.re(z / (1 + z))
    pts = sorted(set([0, u0 / 2, u0, min(2 * u0, mp.mpf('0.5')), 1, 3, 8, 15]))
    g1_fin = -1 + mp.quad(f_int, pts)
    lin_coeff = ct * abs(g1_fin) / mp.sqrt(X)
    g1_mis_actual = lin_coeff ** 2
    print(f"  c={c:8d}, tau={float(tau):.4f}, X={float(X):10.2f}, u0={float(u0):.6f} | "
          f"g1_par(0+)={float(g1_par_0):+.3e}, g1_fin(u0)={float(g1_fin):+.3e} | "
          f"ct*|g1|/sqrt(X)={float(lin_coeff):.4f}, g1_mis={float(g1_mis_actual):.4f} "
          f"(claimed <= 0.2193 -> {'OK' if g1_mis_actual <= 0.2193 else 'VIOLATED!!'})")

print("\n=== TEST 2: c_rein, c_re5, c_re50 at TRUE VALID Piece G3 points (kap ~ 0.448 <= 0.45) ===")
y = mp.mpf('0.45') / (2 * pi)
al_b = y / (1 - y ** 2 / 5); th_m = al_b / 2
te_lo = (th_m / 2) * mp.cos(th_m / 2) / mp.sin(th_m / 2)
S_b = mp.sin(al_b) / al_b
gam_m = mp.mpf('0.45') ** 2 / (2 * pi ** 2 * S_b ** 2 * mp.cos(al_b) ** 3)
c_re5_s  = te_lo * (mp.sin(mp.mpf('0.080')) / mp.mpf('0.080')) ** 2 / mp.mpf('1.05')
c_rein_s = te_lo * (mp.sin(mp.mpf('0.55') * th_m) / (mp.mpf('0.55') * th_m)) ** 2 / (1 + gam_m * mp.mpf('0.31'))
c_re50_s = te_lo * (mp.log(mp.mpf('1.5')) / mp.mpf('0.5')) * (mp.sin(mp.mpf('0.255')) / mp.mpf('0.255')) ** 2

# kap ~ 2*pi^2 * j*(j+c)/(c*n) = 0.448 => n ~ 2*pi^2 * j*(j+c)/(c*0.448)
for label, (c_, p0_, D_) in [
    ("sj~0.51, kap~0.448", (35, 1, 1)),
    ("sj~0.95, kap~0.448", (11, 20, 1)),
    ("sj~0.99, kap~0.448", (11, 100, 1)),
]:
    r0_ = 1; j_ = c_ * p0_ + D_
    n_ = 2 * int(round(pi ** 2 * j_ * (j_ + c_) / (c_ * 0.448)))
    ptk = make_point(n_, c_, j_, r0_, p0_, D_)
    assert ptk['kap'] <= mp.mpf('0.45') and n_ >= 1007 and c_ * (ptk['kap'] ** 2) < 42
    x5  = mp.asin(mp.sqrt(mp.mpf('0.05') * (1 - ptk['sj']) ** 2 / (4 * ptk['sj']))) / ptk['thc']
    x50 = mp.asin(mp.sqrt(mp.mpf('0.50') * (1 - ptk['sj']) ** 2 / (4 * ptk['sj']))) / ptk['thc']
    print(f"  [{label}] n={n_}, kap={float(ptk['kap']):.5f}, sj={float(ptk['sj']):.5f}, te_eff={float(ptk['te_eff']):.6f} (te_lo={float(te_lo):.6f}), x5={float(x5):.3f}, x50={float(x50):.3f}")
    for nm, u0v, xv, c_s in [
        ("c_rein (u0=0.52, x=0.519)", mp.mpf('0.52'), mp.mpf('0.519'), c_rein_s),
        ("c_re5  (u0=0.52, x=x5)",    mp.mpf('0.52'), x5,              c_re5_s),
        ("c_re5  (u0=0.52, x=0.64)",  mp.mpf('0.52'), mp.mpf('0.64'),  c_re5_s),
        ("c_re50 (u0=0.52, x=x50)",   mp.mpf('0.52'), x50,             c_re50_s),
    ]:
        r = (mp.re(ptk['zeta'](u0v)) - mp.re(ptk['zeta'](xv))) / (ptk['tau'] * (xv ** 2 - u0v ** 2))
        print(f"    {nm:28s}: exact={float(r):.6f} vs script={float(c_s):.6f} ({'OK' if r >= c_s else 'VIOLATED!!'})")
    # Also check global worst case of (Re zeta(u0) - Re zeta(x))/(tau*(x^2-u0^2)) for u0 in [0,0.52], x in (u0, x50]!
    worst_50 = mp.mpf(10)
    worst_pair = None
    for u0f in [0.0, 0.25, 0.52]:
        u0v = mp.mpf(str(u0f))
        for k in range(1, 51):
            xv = u0v + (x50 - u0v) * k / 50
            r = (mp.re(ptk['zeta'](u0v)) - mp.re(ptk['zeta'](xv))) / (ptk['tau'] * (xv ** 2 - u0v ** 2))
            if r < worst_50:
                worst_50 = r; worst_pair = (float(u0v), float(xv))
    print(f"    GLOBAL min on u0<=0.52, T(x)<=0.50: {float(worst_50):.6f} at {worst_pair} vs c_re50={float(c_re50_s):.6f} ({'OK' if worst_50 >= c_re50_s else 'VIOLATED!!'})")
