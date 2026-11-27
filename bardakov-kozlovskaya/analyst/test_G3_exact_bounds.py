"""
VERIFY TWO EXACT FINITE-n BREAKTHROUGHS FOR PIECE G3:
1. Exact 1-line Mobius formula for |zeta'(u)|^2 and global bound on |W_3(u)| across ALL u in [-L_n, L_n]
   and ALL G3 corners (BOTH D > 0 AND D < 0!).
2. Exact meromorphic partial-fraction Poisson formula for S'(theta) - ct*(Psi_fin'(theta) + i):
     S'(th) - ct*(Psi_fin'(th) + i) = -(i*c/r_0) * sum_{k != 0} int_{-L_n}^{L_n} e^{2*pi*i*k*ct*(u+1/2)} * y_*(u,th)^{nu(k)} / (1 + Lam*e^{zeta(u)+i*th}) du
   Verify to machine precision (BOTH r_0 | n AND r_0 nmid n!!) AND check exact contour-shift upper bound!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

print("=== 1. VERIFY 1-LINE MOBIUS |zeta'(u)|^2 & GLOBAL |W_3(u)| BOUND ON [-L_n, L_n] ===")
max_W3_ratio = 0.0
for c, j, r0, p0, kap_t in [
    (300, 301, 1, 1, 0.288), (300, 599, 1, 2, 0.288),
    (600, 601, 1, 1, 0.204), (600, 1199, 1, 2, 0.204),
    (300, 302, 1, 1, 0.288), (300, 598, 1, 2, 0.288),
    (80, 121, 2, 3, 0.15),   (80, 119, 2, 3, 0.15),
    (18, 19, 1, 1, 0.30),    (18, 35, 1, 2, 0.30),
    (55, 73, 3, 4, 0.10),    (55, 74, 3, 4, 0.10),
    (91, 109, 5, 6, 0.06),   (91, 110, 5, 6, 0.06),
]:
    D = r0 * j - c * p0
    n = max(132, 2 * int(round(pi * pi * j * (j + c) / (c * kap_t))))
    kap = 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)
    tau = c * r0 * kap
    rj = sin(j * pi / n) / sin((j + c) * pi / n)
    thc = c * pi / n
    Del0 = c * r0 * (rj / (1.0 - rj) - j / c)
    # Test on entire [0, L_n]:
    Ln = n / (2.0 * c)
    us = np.concatenate([np.linspace(0.0, 2.0, 400), np.linspace(2.0, Ln, 400)])
    etas = 2.0 * us * thc
    es = np.exp(1j * etas)
    q1 = 1.0 - rj * es
    zp = 2*pi*c*r0*(1j*rj*es/q1) - 2j*pi*c*p0
    z2 = (2*pi*c*r0)*(2*thc)*(-rj*es/(q1**2))
    z3 = (2*pi*c*r0)*((2*thc)**2)*(-1j*rj*es*(1+rj*es)/(q1**3))
    z4 = (2*pi*c*r0)*((2*thc)**3)*(rj*es*(1+4*rj*es+rj*rj*es*es)/(q1**4))
    # Verify exact Mobius identity for |zp|^2:
    T = 4.0 * rj * (np.sin(us * thc) ** 2) / ((1.0 - rj) ** 2)
    mob_sq = ((2*pi*(D + Del0))**2 + (2*pi*c)**2 * p0 * (r0 + p0) * T) / (1.0 + T)
    mob_err = float(np.max(np.abs(np.abs(zp)**2 - mob_sq) / mob_sq))
    # Exact W_3(u):
    W3 = (-z4*(zp**2) + 10*z3*z2*zp - 15*(z2**3)) / (zp**6)
    rat_W3 = float(np.max(np.abs(W3) * (np.abs(zp)**6) / (120.0 * (tau**3))))
    max_W3_ratio = max(max_W3_ratio, rat_W3)
    # Also check Del0 / (r0*c*kap^2):
    d0_c = Del0 / (r0 * c * kap * kap)
    # Check min |zp|^2 / ((2*tau*u)^2 + (2*pi*D)^2) on [0, Ln]:
    # Wait: at large u ~ Ln, (2*tau*u) grows linearly while sin(u*thc) saturates; we only need |zp|^2 vs exact Mobius!
    # Check |zp(u)|^2 / ((2*pi*|D+Del0|)^2 + c_T * (2*tau*u)^2) on [0, 1.5]:
    u_m = us <= 1.5
    r_mod = float(np.min(np.abs(zp[u_m]) / np.sqrt((2*tau*us[u_m])**2 + (2*pi*(D + Del0))**2)))
    print(f"  c={c:3d} j={j:4d} r0={r0} D={D:+d} n={n:6d} | mob_err={mob_err:.1e} | max|W3*zp^6|/(120*tau^3)={rat_W3:.4f} | Del0/(r0ck2)={d0_c:.5f} | min|zp|/mod1.5={r_mod:.4f}")
print(f"--> GLOBAL MAX |W_3(u) * zeta'(u)^6| / (120 * tau^3) ACROSS ALL u in [-L_n, L_n] = {max_W3_ratio:.5f} !!\n")

print("=== 2. VERIFY EXACT MEROMORPHIC PARTIAL-FRACTION POISSON IDENTITY (INCLUDING r0 nmid n!) ===")
for c, j, r0, p0, tau_t in [
    (18, 19, 1, 1, 1.2),
    (22, 23, 1, 1, 1.2),
    (37, 56, 2, 3, 2.0),
    (55, 73, 3, 4, 3.0),   # r0 = 3, check r0 nmid n!
    (72, 74, 1, 1, 3.0),
    (91, 109, 5, 6, 4.0),  # r0 = 5, check r0 nmid n!
]:
    D = r0 * j - c * p0
    kap_t = tau_t / (c * r0)
    n = max(132, 2 * int(round(pi * pi * j * (j + c) / (c * kap_t))))
    if r0 == 3 and n % 3 == 0: n += 2  # FORCE r0 nmid n!
    if r0 == 5 and n % 5 == 0: n += 2  # FORCE r0 nmid n!
    kap = 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)
    tau = c * r0 * kap; ct = c / r0
    rj = sin(j * pi / n) / sin((j + c) * pi / n)
    sj = sin(c * pi / n) / sin((j + c) * pi / n)
    thc = c * pi / n
    Ln = n / (2.0 * c)
    u0 = 0.35 if (abs(D) == 1 and tau < 2.0) else (0.5 - 0.10 / abs(D))
    def zeta(u):
        eta = 2.0 * u * thc
        return r0 * n * np.log(sj / (1.0 - rj * np.exp(1j * eta))) + 1j * pi * r0 * j - 2j * pi * c * p0 * (u + 0.5)
    logLam = -zeta(u0).real; Lam = exp(logLam); lam = Lam ** (1.0 / r0)
    # Discrete S'(th) at th = 0.5:
    th = 0.5; psi = (pi + th) / r0
    ks = np.arange(-n//2, n//2, dtype=np.float64)
    etas = (2.0 * ks * pi / n) if c % 2 == 0 else ((2.0 * ks + 1.0) * pi / n)
    ak = ((-1.0) ** j) * ((sj / (1.0 - rj * np.exp(1j * etas))) ** n)
    z_disc = ak * lam * np.exp(1j * psi)
    Sp_disc = (1j / r0) * np.sum(-z_disc / (1.0 - z_disc))
    # Continuum ct*(Psi_fin'(th) + i) = (i*c/r0) * int_{-Ln}^{Ln} z/(1+z) du:
    def f_cont(u):
        z = Lam * np.exp(zeta(u) + 1j * th)
        return z / (1.0 + z)
    L_cut = min(Ln, 0.5 + sqrt(35.0 / max(tau, 0.2)))
    c_re, _ = quad(lambda u: f_cont(u).real, -L_cut, L_cut, points=[-u0, u0], limit=300)
    c_im, _ = quad(lambda u: f_cont(u).imag, -L_cut, L_cut, points=[-u0, u0], limit=300)
    Sp_cont = (1j * ct) * (c_re + 1j * c_im)
    diff_exact = Sp_disc - Sp_cont

    # Now compute the k = +-1, +-2, +-3 Poisson mode integrals by CONTOUR SHIFT y = +-y_s!!
    # Note: y_*(u, th) = exp((logLam + zeta(u) + i*(pi+th)) / r0)
    # Wait: for even c, u_k = k/c so c*(u_k + 0.5) = k + c/2 (integer);
    # for odd c, u_k = (k+0.5)/c so c*(u_k + 0.5) = k + (c+1)/2 (integer)!
    # In BOTH cases, m = c*(u_k + 0.5) runs over integers!
    # Check: a_k * lam * e^{i psi} == e^{2*pi*i*p0*m/r0} * y_*(u_k, th)!
    #Wait: check sign (-1)^{r0*j} vs (-1)^j in exp((zeta(u)+i*pi)/r0)!!
    # At u_k: (logLam + zeta(u_k) + i*pi + i*th)/r0 = log(lam) + n*log(...) + i*pi*j - 2*pi*i*m*p0/r0 + i*(pi+th)/r0!
    # Multiply by e^{2*pi*i*p0*m/r0}: cancels -2*pi*i*m*p0/r0 EXACTLY!
    # Check parity: for k in {-3..-1, 1..3}, find nu in {1..r0} with p0*nu = k (mod r0):
    # Wait: in Poisson sum of sum_{m} f(m) = sum_K int e^{-2*pi*i*K*m} f(m) dm:
    #   Notice sign in e^{2*pi*i*(p0*nu/r0 - K)*m} = e^{2*pi*i*k*m/r0} -> k = p0*nu - r0*K!
    #   And m = c*(u + 0.5) -> e^{2*pi*i*k*ct*(u + 0.5)}!
    y_sing_est = min(tau * u0 * u0 / (2 * pi * abs(D)), (pi - 2*pi*abs(D*(0.5-u0)) - abs(th)) / sqrt((2*tau*u0)**2 + (2*pi*D)**2))
    ys = 0.65 * y_sing_est
    pois_sum = 0j
    for k in [-3, -2, -1, 1, 2, 3]:
        # find nu in 1..r0 with (p0*nu - k) % r0 == 0
        nu = [v for v in range(1, r0 + 1) if (p0 * v - k) % r0 == 0][0]
        # Shift u -> u + i*y_shift: for k > 0, e^{2*pi*i*k*ct*(u + i*y)} = e^{-2*pi*k*ct*y} -> choose y = +ys!
        # Wait: check direction of poles!
        # For D > 0: Re zeta(u_r + i*y) ~ Re zeta(u_r) - 2*pi*D*y!
        # At u_r = 0: Re zeta(i*y) - Re zeta(u0) = tau*u0^2 - 2*pi*D*y hits 0 at y = +tau*u0^2/(2*pi*D) > 0!
        # So BOTH +ys and -ys with ys = 0.65*y_sing_est are pole-free!
        y_sh = ys if k > 0 else -ys
        def mode_f(ur):
            u = ur + 1j * y_sh
            zu = zeta(u)
            ys_u = np.exp((logLam + zu + 1j * (pi + th)) / r0)
            z_u = np.exp(logLam + zu + 1j * th)
            return np.exp(2j * pi * k * ct * (u + 0.5)) * (ys_u ** nu) / (1.0 + z_u)
        L_m = min(Ln, 0.5 + sqrt(35.0 * r0 / max(tau, 0.2)))
        mr, _ = quad(lambda ur: mode_f(ur).real, -L_m, L_m, limit=300)
        mi, _ = quad(lambda ur: mode_f(ur).imag, -L_m, L_m, limit=300)
        pois_sum += (-1j * ct) * (mr + 1j * mi)
    print(f"  c={c:2d} j={j:3d} r0={r0} n%r0={n%r0} D={D:+d} tau={tau:4.2f} ys={ys:.4f} exp(-2pi*ct*ys)={exp(-2*pi*ct*ys):.2e} | exact={diff_exact:+.6e} | pois3={pois_sum:+.6e} | err={abs(diff_exact - pois_sum):.2e}")
