"""
1. TEST D = 1, SMALL tau = 0.3, 0.5, 1.0, 2.0 ACROSS HUGE c_tilde = 15..100000!
   Find the EXACT clean radius choice delta(tau, c_tilde) that works for ALL c_tilde in [15, infinity)!
   Specifically check: at D = 1, tau = 0.5, c_tilde = 10000, what happens at del = 0.15, 0.25, 0.35, del_*?
2. TEST EXACT FINITE-n ROOT DISTORTION (exact a_k vs model a_k^{mod}) ON HARDEST G3 CORNERS
   (especially large c*kap^2 ~ 10..20, small n ~ 132..500)!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

def Re_Psi(D, tau, lam, th):
    L = sqrt(35.0 / max(tau, 1e-4)) + 1.0
    a0 = 0.25 + log(lam) / tau
    u0 = sqrt(a0) if a0 > 0 else 0.0
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    def f(u):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2.0 * pi * D * u + th + (pi if sgnD < 0 else 0.0)
        if r > 1.0:
            ir = 1.0 / r
            return log(r) + 0.5 * log(max(1.0 + 2*ir*cos(ph) + ir*ir, 1e-300))
        return 0.5 * log(max(1.0 + 2*r*cos(ph) + r*r, 1e-300))
    v, _ = quad(f, -L, L, points=[-u0, u0] if u0 > 1e-12 else [0.0], limit=300)
    return -log(lam) + v

def Im_Psi(D, tau, lam, th):
    L = sqrt(35.0 / max(tau, 1e-4)) + 1.0
    a0 = 0.25 + log(lam) / tau
    u0 = sqrt(a0) if a0 > 0 else 0.0
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    def f(u):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2.0 * pi * D * u + th + (pi if sgnD < 0 else 0.0)
        if r > 1.0:
            ir = 1.0 / r
            return np.arctan2(-ir * sin(ph), 1.0 + ir * cos(ph))
        return np.arctan2(r * sin(ph), 1.0 + r * cos(ph))
    v, _ = quad(f, -L, L, points=[-u0, u0] if u0 > 1e-12 else [0.0], limit=300)
    return (2.0 * u0 - 1.0) * th + v

print("=== 1. TEST D = 1 SMALL tau = 0.5, 1.0, 2.0 FOR HUGE c_tilde = 15..100000 ===")
for tau in [0.5, 1.0, 2.0]:
    # Find exact del_*
    dels_g = np.linspace(0.10, 0.495, 40)
    vs_g = [Re_Psi(1, tau, exp(-tau*(0.25 - (0.5 - d)**2)), 0.0) for d in dels_g]
    d_star = float(dels_g[int(np.argmin(vs_g))])
    # Refine del_* by golden section
    lo, hi = max(0.05, d_star - 0.02), min(0.499, d_star + 0.02)
    gr = (sqrt(5.0) - 1.0) / 2.0
    for _ in range(25):
        c1 = hi - gr*(hi - lo); c2 = lo + gr*(hi - lo)
        if Re_Psi(1, tau, exp(-tau*(0.25-(0.5-c1)**2)), 0.0) < Re_Psi(1, tau, exp(-tau*(0.25-(0.5-c2)**2)), 0.0):
            hi = c2
        else:
            lo = c1
    d_star = 0.5 * (lo + hi)
    print(f"\n  D = 1, tau = {tau:.1f} (del_* = {d_star:.4f}, pi*(1-2*del_*) = {pi*(1-2*d_star):.4f}):")
    for del_use, name in [(0.15, "del=0.15"), (0.25, "del=0.25"), (d_star, "del=del*")]:
        u0 = 0.5 - del_use; lam = exp(-tau * (0.25 - u0 * u0))
        p0 = Re_Psi(1, tau, lam, 0.0)
        ths = np.linspace(-pi, pi, 512, endpoint=False)
        re_v = np.array([Re_Psi(1, tau, lam, t) - p0 for t in ths])
        im_v = np.array([Im_Psi(1, tau, lam, t) for t in ths])
        # Also check drop at th = pi and at th = pi*(1-2*del_use)
        drop_pi = p0 - Re_Psi(1, tau, lam, pi)
        rats = []
        for ct in [15, 100, 1000, 10000, 100000]:
            # Careful: for large ct, use fine grid near 0 + coarse tail!
            # Actually 512 grid has step 0.012, fine up to ct*H <= 5000!
            F = np.exp(ct * (re_v + 1j * im_v))
            rats.append(float(np.mean(F).real) / float(np.mean(np.abs(F))))
        print(f"    {name:8s} (drop@pi={drop_pi:.5f}): rat(ct=15,1e2,1e3,1e4,1e5)=({rats[0]:+.4f},{rats[1]:+.4f},{rats[2]:+.4f},{rats[3]:+.4f},{rats[4]:+.4f})")

print("\n=== 2. TEST EXACT FINITE-n VS MODEL DISCRETE ON HIGH c*kap^2 (SMALL n!) G3 POINTS ===")
for c, j, r0, D, target_r3ck2 in [
    (15, 16, 1, 1, 10.0), (20, 21, 1, 1, 15.0), (20, 21, 1, 1, 22.0),
    (30, 31, 1, 1, 15.0), (60, 62, 1, 2, 15.0), (60, 62, 1, 2, 22.0),
    (31, 47, 2, 1, 10.0), (46, 61, 3, 1, 15.0),
]:
    target_kap = sqrt(target_r3ck2 / ((r0 ** 3) * c))
    if r0 * target_kap > 0.45 or target_kap >= 2.5 or c * target_kap * target_kap >= 42: continue
    n_est = int(round(2.0 * pi * pi * j * (j + c) / (c * target_kap)))
    n = max(132, 2 * ((n_est + 1) // 2))
    if j > n // 2 - (c + 1) // 2 - 1: continue
    kap = 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)
    tau = c * r0 * kap; ct = c / r0
    # Find continuum saddle u0:
    lo_u, hi_u = 0.5 - 0.35 / D, 0.5 + 0.08 / D
    ug = np.linspace(lo_u, hi_u, 25)
    u0 = float(ug[int(np.argmin([Re_Psi(D, tau, exp(-tau*(0.25-u*u)), 0.0) for u in ug]))])
    lam_W = exp(-tau * (0.25 - u0 * u0)); lam = lam_W ** (1.0 / r0)

    # Compare exact finite-n log F vs model discrete log F at peak psi = pi/r0 + th/r0 (th in [-1, 1]):
    rj = sin(j * pi / n) / sin((j + c) * pi / n)
    sj = sin(c * pi / n) / sin((j + c) * pi / n)
    M_k = min(n // 2, int(c + sqrt(45.0 * c / max(kap, 1e-8))) + 15)
    ks = np.arange(-M_k, M_k + 1, dtype=np.float64) if c % 2 == 0 else np.arange(-M_k, M_k, dtype=np.float64)
    etas = (2.0 * ks * pi / n) if c % 2 == 0 else ((2.0 * ks + 1.0) * pi / n)
    la_ex = n * np.log(sj / (1.0 - rj * np.exp(1j * etas))) + (1j * pi if j % 2 == 1 else 0.0)
    # Model: u = eta / (2 * theta_c), so -alpha*m*(c+m) = kap*c*(0.25 - u^2), phase = pi*j*(1 + 2*u)
    us = etas / (2.0 * c * pi / n)
    la_mod = (c * kap * (0.25 - us * us)) + 1j * (pi * j * (1.0 + 2.0 * us))

    def eval_F_arr(la, psis):
        rv = la.real + log(lam); big = rv > 0.0; sm = ~big
        lg = -c * (log(lam) + 1j * psis) + (1j * pi if (j*(c+1)+c) % 2 == 1 else 0.0)
        if np.any(big):
            rb = rv[big]; ib = la.imag[big]
            lg = lg + np.sum(rb) + 1j*(np.sum(ib) + len(rb)*pi) + len(rb)*1j*psis
            ch = np.exp(-rb - 1j*ib); ep = np.exp(-1j*psis)
            for k0 in range(0, len(ch), 200):
                lg += np.sum(np.log(1.0 - ch[k0:k0+200][None,:]*ep[:,None]), axis=1)
        if np.any(sm):
            rs = rv[sm]; is_ = la.imag[sm]
            ch = np.exp(rs + 1j*is_); ep = np.exp(1j*psis)
            for k0 in range(0, len(ch), 200):
                lg += np.sum(np.log(1.0 - ch[k0:k0+200][None,:]*ep[:,None]), axis=1)
        return lg

    psis_all = np.linspace(-pi, pi, 512, endpoint=False)
    lg_ex_all = eval_F_arr(la_ex, psis_all)
    sh = float(np.max(lg_ex_all.real))
    F_ex = np.exp(lg_ex_all - sh)
    rat_ex = float(np.mean(F_ex).real) / float(np.mean(np.abs(F_ex)))

    ths_pk = np.linspace(-1.0, 1.0, 21)
    psis_pk = (pi + ths_pk) / r0
    le = eval_F_arr(la_ex, psis_pk); lm = eval_F_arr(la_mod, psis_pk)
    diff = (le - le[10]) - (lm - lm[10])
    print(f"  n={n:5d} c={c:2d} j={j:2d} r0={r0} D={D} kap={kap:.3f} r0k={r0*kap:.3f} ck2={c*kap*kap:4.1f} tau={tau:4.1f} | rat={rat_ex:+.5f} | fin-mod PK: |dRe|<={np.max(np.abs(diff.real)):.4f} |dIm|<={np.max(np.abs(diff.imag)):.4f}")
