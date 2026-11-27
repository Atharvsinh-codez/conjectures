"""
TEST EXACT FINITE-n CONTINUUM Psi_fin(-Lambda e^{i theta}) + 1-LINE EXACT IBP ON PIECE G3
(INCLUDING EXTREME HIGH-c*kap^2 CORNER c = 300, j = 301, r0 = 1, D = 1, kap = 0.25, c*kap^2 = 18.75!!)!

For ANY finite (n, c, j) and (p_0, r_0) with r_0*j = c*p_0 + D:
  Define the EXACT finite-n phase-unwrapped function on u in [-n/(2c), n/(2c)]:
    eta(u) := 2 * u * theta_c,
    zeta(u) := (r_0 * n) * log((r_3/r_1) / (1 - s_j e^{i eta(u)})) + i * pi * r_0 * j - 2*pi*i * c * p_0 * (u + 0.5)!
  Notice:
    1) At every grid point u_k = (k/c for even c, (k+0.5)/c for odd c):
         c * (u_k + 0.5) IS AN INTEGER! So exp(zeta(u_k)) == a_k^{r_0} IDENTICALLY!!
    2) Its EXACT u-derivative is:
         zeta'(u) = (2 * c * r_0 * pi) * (i * s_j e^{i eta} / (1 - s_j e^{i eta})) - 2*pi*i * c * p_0!
         At u = +-0.5 (eta = +-theta_c):
           Re zeta'(+-0.5) = -+ (r_0 * c * kap) * ... = -+ tau * (1 + O(alpha_1^2))!
           Im zeta'(+-0.5) = -2*pi * D_bdy, where D_bdy = c * (r_0 * nu_bdy - p_0) = D - r_0 * c * (j/c - nu_bdy)!!
    3) Define EXACT finite-n continuum:
         Psi_fin(-Lambda e^{i th}) := -log(Lambda) - i*th + int_{-n/(2c)}^{n/(2c)} log(1 + Lambda * exp(zeta(u) + i*th)) du!
    4) EXACT 1-LINE SECOND-DERIVATIVE IBP ON Psi_fin (ZERO PARABOLA APPROXIMATION!!):
         Psi_fin''(th) = int_{-n/(2c)}^{n/(2c)} (z/(1+z)) * (-zeta''(u) / zeta'(u)^2) du
                       = [ 1 / zeta'(u) ]_{-u_0}^{u_0}  +  int E(u, th) * (-zeta''(u) / zeta'(u)^2) du!!
       CHECK THAT ONE-LINE EXACT IDENTITY!!
         d/du (1 / zeta'(u)) = -zeta''(u) / zeta'(u)^2!!
         SO int_{-u_0}^{u_0} (-zeta''(u) / zeta'(u)^2) du = -(1 / zeta'(u_0) - 1 / zeta'(-u_0)) = -2 Re(1 / zeta'(u_0))!!
         WAIT A SECOND!! CHECK -2 Re(1 / zeta'(u_0)) WHEN zeta'(u_0) = -2*tau*u_0 - 2*pi*i*D:
           -2 Re( 1 / (-2*tau*u_0 - 2*pi*i*D) ) = -2 * (-2*tau*u_0) / (4*tau^2*u_0^2 + 4*pi^2*D^2) = +... wait, check sign of i^2 = -1!
           With i^2 = -1: Psi_fin''(th) = +2 Re(1 / zeta'(u_0)) < 0 (EXACTLY -H_0 FOR THE PARABOLA, AND EXACTLY -H_fin FOR FINITE n!!)!

TEST THIS RIGHT NOW ON c = 300, j = 301, r0 = 1, D = 1, kap = 0.25 (ck2 = 18.75, tau = 75) AND OTHER G3 CORNERS!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

def test_exact_fin_G3(c, j, r0, D, target_kap):
    n_est = int(round(2.0 * pi * pi * j * (j + c) / (c * target_kap)))
    n = max(132, 2 * ((n_est + 1) // 2))
    kap = 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)
    tau = c * r0 * kap; ct = c / r0
    p0 = round(r0 * j / c)
    assert r0 * j - c * p0 == D or r0 * j - c * p0 == -D
    rj = sin(j * pi / n) / sin((j + c) * pi / n)
    sj = sin(c * pi / n) / sin((j + c) * pi / n)
    th_c = c * pi / n

    # Exact zeta(u) and zeta'(u) on u in [-n/(2c), n/(2c)]:
    def zeta(u):
        eta = 2.0 * u * th_c
        return r0 * n * np.log(sj / (1.0 - rj * np.exp(1j * eta))) + 1j * pi * r0 * j - 2j * pi * c * p0 * (u + 0.5)
    def dzeta(u):
        eta = 2.0 * u * th_c
        eie = np.exp(1j * eta)
        return (2.0 * c * r0 * pi) * (1j * rj * eie / (1.0 - rj * eie)) - 2j * pi * c * p0
    def d2zeta(u):
        eta = 2.0 * u * th_c
        eie = np.exp(1j * eta)
        # d/deta (i rj e^{i eta}/(1 - rj e^{i eta})) = -rj e^{i eta}/(1 - rj e^{i eta})^2
        return (2.0 * c * r0 * pi) * (2.0 * th_c) * (-rj * eie / ((1.0 - rj * eie) ** 2))

    # Choose u0 ~ 0.5 (at lam = 1, u0 = 0.5 EXACTLY! Check zeta(0.5) = 0!):
    # Let's check exact finite-n continuum saddle u0_* in [0.35, 0.52]:
    L_u = min(n / (2.0 * c), 0.5 + sqrt(35.0 / max(tau, 0.1)))
    def Re_Psi_fin(u0, th):
        lam_W = exp(-zeta(u0).real)
        def f(u):
            z = lam_W * np.exp(zeta(u) + 1j * th)
            az = abs(z)
            if az > 1.0:
                return log(az) + log(max(abs(1.0 + 1.0 / z), 1e-300))
            return log(max(abs(1.0 + z), 1e-300))
        v, _ = quad(f, -L_u, L_u, points=[-u0, u0], limit=300, epsabs=1e-10, epsrel=1e-10)
        return -log(lam_W) + v

    # Find exact finite-n continuum saddle u0:
    ug = np.linspace(0.5 - 0.15 / abs(D), 0.5 + 0.02 / abs(D), 21)
    vs = [Re_Psi_fin(u, 0.0) for u in ug]
    u0 = float(ug[int(np.argmin(vs))])
    lam_W = exp(-zeta(u0).real)
    lam = lam_W ** (1.0 / r0)

    # Check exact H_fin_0 = -2 * Re(1 / dzeta(u0)) vs actual curvature -Re_Psi_fin''(0)!
    # Wait: check sign! At u = +u0 > 0, Re(zeta(u)) is DECREASING, so Re(dzeta(u0)) < 0 => 2 * Re(1/dzeta(u0)) < 0!
    # Wait: check exact IBP identity:
    #   d^2/dth^2 log(1 + z) = i * d/dth (z/(1+z)) = i * (i / dzeta) * d/du (z/(1+z)) = -(1 / dzeta) * d/du (z/(1+z))!
    #   IBP: -int (1/dzeta) d/du(z/(1+z)) du = +int (z/(1+z)) * (-d2zeta / dzeta^2) du!
    #   Step-fn 1_{|u|<u0}: int_{-u0}^{u0} (-d2zeta / dzeta^2) du = [ 1 / dzeta(u) ]_{-u0}^{u0}
    #                      = 1/dzeta(u0) - 1/dzeta(-u0) = 2 * Re(1 / dzeta(u0)) < 0!
    #   Wait: check dzeta(-u0) = -conj(dzeta(u0))!
    #     Re dzeta(-u0) = -Re dzeta(u0), Im dzeta(-u0) = +Im dzeta(u0)!
    #     SO 1/dzeta(u0) - 1/dzeta(-u0) = 1/(a+ib) - 1/(-a+ib) = (-2a)/(a^2+b^2)?? NO:
    #     1/(a+ib) - 1/(-a+ib) = ((-a+ib) - (a+ib)) / (-a^2 - b^2) = (-2a) / (-(a^2+b^2)) = +2a / (a^2+b^2) = 2 Re(1/dzeta(u0))!
    #     Since a = Re dzeta(u0) < 0, 2*Re(1/dzeta(u0)) < 0!! EXACT!!
    H_fin0 = -2.0 * (1.0 / dzeta(u0)).real
    p0_v = Re_Psi_fin(u0, 0.0)
    H_num = -2.0 * (Re_Psi_fin(u0, 0.005) - p0_v) / (0.005 ** 2)

    # Compare EXACT finite-n discrete F_pure vs ct * Psi_fin at th = 0.5!!
    M_k = min(n // 2, int(c + sqrt(45.0 * c / max(kap, 1e-8))) + 15)
    ks = np.arange(-M_k, M_k + 1, dtype=np.float64) if c % 2 == 0 else np.arange(-M_k, M_k, dtype=np.float64)
    etas = (2.0 * ks * pi / n) if c % 2 == 0 else ((2.0 * ks + 1.0) * pi / n)
    la = n * np.log(sj / (1.0 - rj * np.exp(1j * etas))) + (1j * pi if j % 2 == 1 else 0.0)
    rv = la.real + log(lam); big = rv > 0.0; sm = ~big
    def eval_disc(psi):
        v = -c * (log(lam) + 1j * psi) + (1j * pi if (j*(c+1)+c) % 2 == 1 else 0.0)
        if np.any(big):
            v += np.sum(rv[big]) + 1j*(np.sum(la.imag[big]) + np.sum(big)*(pi + psi))
            v += np.sum(np.log(1.0 - np.exp(-rv[big] - 1j*(la.imag[big] + psi))))
        if np.any(sm):
            v += np.sum(np.log(1.0 - np.exp(rv[sm] + 1j*(la.imag[sm] + psi))))
        return v
    d0 = eval_disc(pi / r0); d1 = eval_disc(pi / r0 + 0.5 / r0)
    d_disc = (d1 - d0).real
    d_fin  = ct * (Re_Psi_fin(u0, 0.5) - p0_v)
    # Also exact rat_L1 on |w| = lam!
    psis = np.linspace(-pi, pi, 512, endpoint=False)
    lg = -c * (log(lam) + 1j * psis) + (1j * pi if (j*(c+1)+c) % 2 == 1 else 0.0)
    if np.any(big):
        rb = rv[big]; ib = la.imag[big]
        lg = lg + np.sum(rb) + 1j*(np.sum(ib) + len(rb)*pi) + len(rb)*1j*psis
        ch = np.exp(-rb - 1j*ib); ep = np.exp(-1j*psis)
        for k0 in range(0, len(ch), 200): lg += np.sum(np.log(1.0 - ch[k0:k0+200][None,:]*ep[:,None]), axis=1)
    if np.any(sm):
        rs = rv[sm]; is_ = la.imag[sm]
        ch = np.exp(rs + 1j*is_); ep = np.exp(1j*psis)
        for k0 in range(0, len(ch), 200): lg += np.sum(np.log(1.0 - ch[k0:k0+200][None,:]*ep[:,None]), axis=1)
    sh = float(np.max(lg.real)); F = np.exp(lg - sh)
    rat = float(np.mean(F).real) / float(np.mean(np.abs(F)))
    D_eff = abs(dzeta(u0).imag) / (2 * pi)
    print(f"  n={n:6d} c={c:3d} j={j:3d} r0={r0} D={D} (D_eff={D_eff:.3f}) ck2={c*kap*kap:5.2f} tau={tau:5.1f} | H0={H_fin0:.5f} Hnum={H_num:.5f} | dRe=({d_disc:+.4f},{d_fin:+.4f}) err={abs(d_disc-d_fin):.5f} | rat={rat:+.5f}")

print("=== EXACT FINITE-n CONTINUUM Psi_fin + 1-LINE IBP ON EXTREME G3 CORNERS ===")
for c, j, r0, D, k in [
    (300, 301, 1, 1, 0.25),   # ck2 = 18.75, tau = 75!!
    (300, 301, 1, 1, 0.15),   # ck2 = 6.75, tau = 45
    (200, 202, 1, 2, 0.25),   # ck2 = 12.5, D = 2
    (60, 62, 1, 2, 0.30),     # ck2 = 5.4, ct/D2 = 15
    (15, 16, 1, 1, 0.30),     # ck2 = 1.35, ct/D2 = 15
    (31, 47, 2, 1, 0.15),     # r0 = 2, r0k = 0.30
    (46, 61, 3, 1, 0.10),     # r0 = 3, r0k = 0.30
    (76, 91, 5, 1, 0.06),     # r0 = 5, r0k = 0.30
]:
    test_exact_fin_G3(c, j, r0, D, k)
