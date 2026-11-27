"""
VERIFY THE 1-LINE INTEGRATED MEROMORPHIC POISSON + RESIDUE IDENTITY FOR PIECE G3:
  [log F_pure(th) - log F_pure(0)] - ct*[Psi_fin(th) - Psi_fin(0)]
  = sum_{k != 0} (1 / (2*pi*i*k)) * int_{-Ln + i*yk}^{Ln + i*yk} e^{2*pi*i*k*ct*(u+1/2)}
      * zeta'(u) * [ y_*(u,th)^{nu(k)}/(1+z(u,th)) - y_*(u,0)^{nu(k)}/(1+z(u,0)) ] du
    + EXACT RESIDUES AT POLES CROSSED!!

Check:
1. Exact identity match (to 1e-12!) on ALL corners (D > 0, D < 0, r_0 = 1,2,3,5, r_0 nmid n, small/med/large tau)!
2. Check EXACT upper bound from Residue + Shifted-Line Integral across ALL G3 corners
   (including ct = 15*D^2, ct*tau^2/D^4 = 15, r_0*kap = 0.30, r_0^3*c*kap^2 = 25)!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

def test_integrated_poisson(c, j, r0, p0, tau_t, th_test=0.5):
    D = r0 * j - c * p0
    kap_t = tau_t / (c * r0)
    n = max(132, 2 * int(round(pi * pi * j * (j + c) / (c * kap_t))))
    if r0 == 3 and n % 3 == 0: n += 2
    if r0 == 5 and n % 5 == 0: n += 2
    kap = 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)
    tau = c * r0 * kap; ct = c / r0
    rj = sin(j * pi / n) / sin((j + c) * pi / n)
    sj = sin(c * pi / n) / sin((j + c) * pi / n)
    thc = c * pi / n; Ln = n / (2.0 * c)
    u0 = min(0.35, 1.49 / ((ct * (tau ** 3)) ** 0.25)) if (abs(D) == 1 and tau <= 1.8) else (0.5 - 0.10 / abs(D))

    def zeta(u):
        eta = 2.0 * u * thc
        return r0 * n * np.log(sj / (1.0 - rj * np.exp(1j * eta))) + 1j * pi * r0 * j - 2j * pi * c * p0 * (u + 0.5)
    def dzeta(u):
        eta = 2.0 * u * thc; e = np.exp(1j * eta)
        return 2*pi*c*r0*(1j*rj*e/(1.0 - rj*e)) - 2j*pi*c*p0

    logLam = -zeta(u0).real; Lam = exp(logLam); lam = Lam ** (1.0 / r0)

    # 1. Exact discrete [log F(th) - log F(0)]:
    ks = np.arange(-n//2, n//2, dtype=np.float64)
    etas = (2.0 * ks * pi / n) if c % 2 == 0 else ((2.0 * ks + 1.0) * pi / n)
    la = n * np.log(sj / (1.0 - rj * np.exp(1j * etas))) + (1j * pi if j % 2 == 1 else 0.0)
    rv = la.real + log(lam); big = rv > 0.0; sm = ~big
    def eval_disc(th):
        psi = (pi + th) / r0
        v = -c * (log(lam) + 1j * psi)
        if np.any(big):
            v += np.sum(rv[big]) + 1j * (np.sum(la.imag[big]) + np.sum(big) * (pi + psi))
            v += np.sum(np.log(1.0 - np.exp(-rv[big] - 1j * (la.imag[big] + psi))))
        if np.any(sm):
            v += np.sum(np.log(1.0 - np.exp(rv[sm] + 1j * (la.imag[sm] + psi))))
        return v
    d_disc = eval_disc(th_test) - eval_disc(0.0)

    # 2. Exact continuum ct * [Psi_fin(th) - Psi_fin(0)]:
    L_cut = min(Ln, 0.5 + sqrt(40.0 / max(tau, 0.1)))
    def Psi_fin(th):
        def fr(u):
            z = Lam * np.exp(zeta(u) + 1j * th); az = abs(z)
            return (log(az) + log(abs(1+1/z))) if az > 1 else log(abs(1+z))
        def fi(u):
            z = Lam * np.exp(zeta(u) + 1j * th); az = abs(z)
            return np.angle(1 + 1/z) if az > 1 else np.angle(1 + z)
        vr, _ = quad(fr, -L_cut, L_cut, points=[-u0, u0], limit=300, epsabs=1e-11, epsrel=1e-11)
        vi, _ = quad(fi, -L_cut, L_cut, points=[-u0, u0], limit=300, epsabs=1e-11, epsrel=1e-11)
        return (-logLam + vr) + 1j * ((2.0 * u0 - 1.0) * th + vi)
    d_cont = ct * (Psi_fin(th_test) - Psi_fin(0.0))
    err_exact = d_disc - d_cont
    # Normalize Im mod 2*pi:
    err_exact = err_exact.real + 1j * float(np.angle(np.exp(1j * err_exact.imag)))

    # 3. Compute exact pole(s) in the pole-containing half-plane (sgn(D) > 0 => upper half-plane, sgn(D) < 0 => lower)!
    sgnD = 1 if D > 0 else -1
    # Shift in pole-containing direction PAST all inner poles to Y_past where Re(zeta(u_r + i*sgnD*Y)) - Re(zeta(u0)) <= -1.0!
    # Wait: at u_r = 0, Re(zeta(0)) - Re(zeta(u0)) ~ tau*u0^2, and d/dy Re(zeta) = -Im(zeta') ~ -2*pi*D!
    # So Y_past = (tau*u0^2 + 1.0) / (2*pi*|D|) makes |z| <= e^{-1} everywhere on Im u = sgnD * Y_past!
    # Wait: check if Y_past is small enough that 2*Y_past*thc << 1 (so Im(zeta') stays near 2*pi*D)!
    Y_past = min(0.25, (tau * u0 * u0 + 1.0) / (2 * pi * abs(D)))
    # Check max |z| on Im u = sgnD * Y_past:
    urs = np.linspace(-min(Ln, 2.0), min(Ln, 2.0), 400)
    z_line = np.abs(Lam * np.exp(zeta(urs + 1j * sgnD * Y_past)))
    max_z_past = float(np.max(z_line))

    # Find ALL poles of 1 + Lam * exp(zeta(u) + i*th) = 0 in 0 < sgnD * Im(u) < Y_past via Newton!
    def find_poles(th):
        poles = []
        # Initial guesses along |u_r| <= u0 where Im(zeta(u_r)) + th = (2*m + 1)*pi
        for m in range(-abs(D) - 2, abs(D) + 3):
            u_g = 0.5 + ((2 * m + 1) * pi - th) / (2 * pi * D)
            if abs(u_g) > u0 + 0.15: continue
            y_g = sgnD * max(1e-4, tau * max(0.01, u0 * u0 - u_g * u_g) / (2 * pi * abs(D)))
            u = u_g + 1j * y_g
            for _ in range(25):
                f = logLam + zeta(u) + 1j * th - 1j * (2 * m + 1) * pi
                # wrap f.imag into (-pi, pi]
                f = f.real + 1j * (((f.imag + pi) % (2 * pi)) - pi)
                du = -f / dzeta(u)
                u += du
                if abs(du) < 1e-13: break
            if abs(1.0 + Lam * np.exp(zeta(u) + 1j * th)) < 1e-10 and 0 < sgnD * u.imag < Y_past:
                if not any(abs(u - p) < 1e-6 for p in poles):
                    poles.append(u)
        return poles

    p_th = find_poles(th_test)
    p_0  = find_poles(0.0)
    # Exact residue sum over ALL k with sgn(k) == sgnD:
    # Check residue sign carefully!!
    # For sgnD = +1 (k > 0, shift UP from 0 to +Y_past): contour runs +R then - (R + iY) -> encloses poles COUNTER-CLOCKWISE (+2*pi*i * Res)!
    #Wait: int_R = int_{R+iY} + 2*pi*i * sum Res!
    # Our integrand is (1 / (2*pi*i*k)) * e^{2*pi*i*k*ct*(u+1/2)} * zeta'(u) * [y_*(u,th)^nu/(1+z(u,th)) - y_*(u,0)^nu/(1+z(u,0))]!
    # At pole u_p of 1 + z(u, th) = 0: Res = (1/(2*pi*i*k)) * e^{2*pi*i*k*ct*(u_p+1/2)} * zeta'(u_p) * y_*(u_p,th)^nu / (-zeta'(u_p))
    #   Wait: 2*pi*i * Res = -(1/k) * e^{2*pi*i*k*ct*(u_p+1/2)} * y_*(u_p,th)^{nu(k)}!
    # For sgnD = -1 (k < 0, shift DOWN from 0 to -Y_past): int_R = int_{R-iY} - 2*pi*i * sum Res!
    #   Since k < 0: -2*pi*i * (1/(2*pi*i*k)) * (-1) = +1/k = -1/|k|!!
    # Check sum_{k >= 1} (-1/k) * e^{2*pi*i*sgnD*k*ct*(u_p+1/2)} * y_*(u_p, th)^{nu(sgnD*k)}:
    res_sum = 0j
    for poles, sgn_term, th_v in [(p_th, +1.0, th_test), (p_0, -1.0, 0.0)]:
        for up in poles:
            ys_p = np.exp((logLam + zeta(up) + 1j * (pi + th_v)) / r0)
            for kk in range(1, 20):
                k = sgnD * kk
                nu = [v for v in range(1, r0 + 1) if (p0 * v - k) % r0 == 0][0]
                res_sum += sgn_term * (-1.0 / kk) * np.exp(2j * pi * k * ct * (up + 0.5)) * (ys_p ** nu)

    # Also compute exact analytic upper bound from residues alone: sum_poles -log(1 - exp(-2*pi*ct*|Im u_p|))!
    ub_res = sum(-log(1.0 - exp(-2 * pi * ct * abs(up.imag))) for up in (p_th + p_0))
    # Estimate min pole height among p_th + p_0:
    min_yp = min([abs(up.imag) for up in (p_th + p_0)] + [999.0])
    print(f"  c={c:3d} j={j:4d} r0={r0} n%r0={n%r0} D={D:+d} ct={ct:4.1f} tau={tau:5.2f} u0={u0:.3f} | #p=({len(p_th)},{len(p_0)}) min_yp={min_yp:.4f} max|z|_Y={max_z_past:.3f} | exact={abs(err_exact):.5f} res={abs(res_sum):.5f} |exact-res|={abs(err_exact-res_sum):.2e} | ub_res={ub_res:.5f}")

print("=== TEST 1-LINE INTEGRATED MEROMORPHIC POISSON + RESIDUE IDENTITY ACROSS ALL G3 CORNERS ===")
for c, j, r0, p0, tau_t in [
    (15, 16, 1, 1, 1.0),    # ct = 15, ct*tau^2 = 15!
    (18, 19, 1, 1, 0.92),   # ct = 18, ct*tau^2 = 15!
    (18, 35, 1, 2, 0.92),   # D = -1!
    (18, 19, 1, 1, 1.80),   # tau = 1.8
    (30, 31, 1, 1, 0.71),   # ct = 30, ct*tau^2 = 15!
    (37, 56, 2, 3, 1.50),   # r0 = 2, D = +1
    (55, 73, 3, 4, 2.50),   # r0 = 3, D = -1, r0 nmid n!
    (60, 62, 1, 1, 2.00),   # D = +2, ct/D^2 = 15!
    (72, 74, 1, 1, 3.00),   # D = +2, ct/D^2 = 18!
    (91, 109, 5, 6, 4.00),  # r0 = 5, D = -1, r0 nmid n!
    (100, 101, 1, 1, 25.0), # D = +1, large tau = 25!
    (300, 301, 1, 1, 75.0), # D = +1, extreme tau = 75, ck2 = 18.75!
]:
    test_integrated_poisson(c, j, r0, p0, tau_t, th_test=0.5)
