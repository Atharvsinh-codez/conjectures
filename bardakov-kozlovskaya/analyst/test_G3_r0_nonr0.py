"""
CHECK FOR PIECE G3 (D >= 1, r0 >= 2):
1. Compare exact F_pure(lam e^{i psi}) vs the r0-multisection main term F_main^{(r0)}(lam^{r0} e^{i r0 psi})
   AND vs continuum c_tilde * Psi_{D, tau}(-Lambda e^{i theta}) across ALL theta in [-pi, pi] (peak AND valley!!)!
2. Specifically check:
   (a) On the peak |theta| <= 0.65 * (pi - 2*pi*|del|): max |Re(log F_pure - ct*Psi)| and max |Im(log F_pure - ct*Psi)|!
   (b) On the valley |theta| > 0.65 * (pi - 2*pi*|del|): max_{valley} (Re log F_pure(psi) - Re log F_pure(peak))
       vs continuum valley max!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

def Re_Psi_c(D, tau, lam, th):
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

def Im_Psi_c(D, tau, lam, th):
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

print("=== CHECK PEAK ERROR & VALLEY UPPER BOUND ON EXACT FINITE-n G3 POINTS (r0 = 1, 2, 3, 5) ===")
for c, j, r0, D, target_tau in [
    (15, 16, 1, 1, 1.5), (15, 16, 1, 1, 3.0), (20, 21, 1, 1, 6.0),
    (31, 47, 2, 1, 3.0), (31, 47, 2, 1, 6.0),
    (46, 61, 3, 1, 4.0), (46, 61, 3, 1, 8.0),
    (76, 91, 5, 1, 5.0), (76, 91, 5, 1, 10.0),
    (60, 62, 1, 2, 4.0), (60, 62, 1, 2, 10.0),
    (121, 183, 2, 2, 8.0), # 2*183 - 121*3 = 366 - 363 = 3 (skip if != D)
    (122, 184, 2, 2, 8.0), # 2*184 - 122*3 = 368 - 366 = 2!
    (135, 138, 1, 3, 8.0),
]:
    p0 = round(r0 * j / c)
    if abs(r0 * j - c * p0) != D: continue
    target_kap = target_tau / (c * r0)
    n_est = int(round(2.0 * pi * pi * j * (j + c) / (c * target_kap)))
    n = max(132, 2 * ((n_est + 1) // 2))
    kap = 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)
    tau = c * r0 * kap; ct = c / r0
    if D == 1 and tau < 3.0:
        u0 = 0.35
    else:
        lo_u, hi_u = 0.5 - 0.35 / D, 0.5 + 0.08 / D
        ug = np.linspace(lo_u, hi_u, 25)
        vs = [Re_Psi_c(D, tau, exp(-tau*(0.25 - u*u)), 0.0) for u in ug]
        u0 = float(ug[int(np.argmin(vs))])
    lam_W = exp(-tau * (0.25 - u0 * u0))
    lam = lam_W ** (1.0 / r0)
    delta = D * (0.5 - u0)

    # Exact finite-n roots:
    rj = sin(j * pi / n) / sin((j + c) * pi / n)
    sj = sin(c * pi / n) / sin((j + c) * pi / n)
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

    # Check ALL r0 sectors k = 0..r0-1 at theta in [-pi, pi] (64 points)!
    th_lim = 0.65 * (pi - 2 * pi * abs(delta))
    ths = np.linspace(-pi, pi, 64, endpoint=False)
    re_c = np.array([ct * Re_Psi_c(D, tau, lam_W, t) for t in ths])
    im_c = np.array([ct * Im_Psi_c(D, tau, lam_W, t) for t in ths])
    c0 = ct * Re_Psi_c(D, tau, lam_W, 0.0)

    max_pk_re = 0.0; max_pk_im = 0.0; max_val_excess = -1e9
    for sec in range(r0):
        psi_c = (2 * sec + 1) * pi / r0
        d0 = eval_disc(psi_c)
        for idx, t in enumerate(ths):
            dv = eval_disc(psi_c + t / r0)
            dre = (dv.real - d0.real) - (re_c[idx] - c0)
            if abs(t) <= th_lim:
                dim = float(np.angle(np.exp(1j * (dv.imag - im_c[idx]))))
                max_pk_re = max(max_pk_re, abs(dre))
                max_pk_im = max(max_pk_im, abs(dim))
            else:
                # Valley excess of disc ABOVE continuum:
                max_val_excess = max(max_val_excess, dre)
    print(f"  c={c:3d} j={j:3d} r0={r0} D={D} ct/D2={ct/(D*D):4.1f} tau={tau:4.1f} (r0^2*k={r0*r0*kap:.3f}): PK |dRe|<={max_pk_re:.4f} |dIm|<={max_pk_im:.4f} | VAL max(disc-cont)={max_val_excess:+.4f}")
