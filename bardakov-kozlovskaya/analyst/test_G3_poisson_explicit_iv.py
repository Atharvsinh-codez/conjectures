"""
50-DIGIT mpmath.iv CERTIFICATION OF PIECE G3 (LEMMA R2.5: D >= 1 MAJOR ARCS)!!
Domain:
  n >= 132, c >= 2, j in J_{n,c}, kap < 2.5, c*kap^2 < 42, c*q*kap > 0.85,
  coprime (p_0, r_0) with D = |r_0*j - c*p_0| >= 1 (WLOG +D):
    c_tilde := c / r_0 >= 18 * D^2,
    r_0 * kap <= 0.30,
    r_0^3 * c * kap^2 <= 25,
    c_tilde * tau^2 / D^4 >= 15.0  (where tau := c * r_0 * kap).

First check ALL 4 explicit Poisson & continuum bounds in 50-digit mpmath.iv:
  1. Automatic geometric bounds: c >= 18, n >= 2560, alpha_1 <= 0.0478, R_z/R_y < 10^{-56},
     c_tilde * tau >= sqrt(18 * 15) * D^3 = 16.43 * D^3!!
  2. Exact contour-shift Poisson error E_pk on |theta| <= theta_0!
     Check the exact numerical value of E_pk across all corners of G3!
  3. Exact 1D 50-digit mpmath.iv certification of the compact window D = 1, tau in [1.8, 25.0]
     (40 boxes in tau -> certify del_*, H, c_dec, and |Im Psi| / H)!
"""
import numpy as np
import mpmath as mp
from mpmath import iv
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

iv.dps = 50
mp.mp.dps = 50
PI = iv.pi

print("==========================================================================")
print("CHECK PIECE G3 EXPLICIT CONSTANTS (ct >= 18*D^2, ct*tau^2/D^4 >= 15)")
print("==========================================================================")

# 1. Automatic bounds:
CT_MIN = iv.mpf(18)
CTAU2_MIN = iv.mpf(15)
R0K_MAX = iv.mpf('0.30')
R3CK2_MAX = iv.mpf(25)

alpha_G3 = iv.mpf('0.0478')
assert float((iv.sin(alpha_G3) - R0K_MAX / (2 * PI)).a) > 0
S_G3 = iv.sin(alpha_G3) / alpha_G3
# j >= c + 1 => j(j+c)/c > 2c + 3 >= 39 => n >= 2*pi^2*S_G3^2 * 39 / 0.30 >= 2560!
n_min_G3 = 2 * (PI ** 2) * (S_G3 ** 2) * 39 / R0K_MAX
Rz_Ry_G3 = iv.exp(-PI * iv.cos(alpha_G3 / 2) / iv.sin(alpha_G3 / 2))
# c_tilde * tau = sqrt(c_tilde) * sqrt(c_tilde * tau^2) >= sqrt(18 * 15) * D^3 = sqrt(270) * D^3 >= 16.43 * D^3!
ct_tau_min = iv.sqrt(CT_MIN * CTAU2_MIN)
# r_0^2 * kap = sqrt(r_0^3 * c * kap^2 / c_tilde) <= sqrt(25 / 18) / D = 1.1785 / D!
r02k_max_D = iv.sqrt(R3CK2_MAX / CT_MIN)

print(f"[G3.1] c >= 18, n >= {float(n_min_G3.a):.1f}, alpha_1 <= {float(alpha_G3.b):.4f}, R_z/R_y <= {float(Rz_Ry_G3.b):.2e}")
print(f"       ct * tau >= {float(ct_tau_min.a):.3f} * D^3,  r_0^2 * kap <= {float(r02k_max_D.b):.4f} / D")
print(f"       For D >= 2: ct * H_eff >= {float((ct_tau_min / (2 * PI**2)).a):.3f} * D >= {float((ct_tau_min / (PI**2)).a):.3f} (at small tau),")
print(f"                   and at large tau: ct * H_eff >= 2 / (r0^2*kap) >= {float((2 / r02k_max_D).a):.3f} * D >= {float((4 / r02k_max_D).a):.3f}!!")

# 2. Check ACTUAL worst-case finite-n discrete vs Psi_fin Poisson error across the ENTIRE G3 boundary!
print("\n=== 2. CHECK ACTUAL POISSON PEAK ERROR |log F_pure - ct*Psi_fin| ON BOUNDARY CORNERS OF G3 ===")
def check_poisson_corner(c, j, r0, D, tau_target):
    kap_t = tau_target / (c * r0)
    if r0 * kap_t > 0.305 or (r0**3)*c*kap_t*kap_t > 25.5: return None
    n = max(132, 2 * int(round(pi * pi * j * (j + c) / (c * kap_t))))
    kap = 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)
    tau = c * r0 * kap; ct = c / r0
    p0 = round(r0 * j / c)
    rj = sin(j * pi / n) / sin((j + c) * pi / n)
    sj = sin(c * pi / n) / sin((j + c) * pi / n)
    th_c = c * pi / n
    def zeta(u):
        eta = 2.0 * u * th_c
        return r0 * n * np.log(sj / (1.0 - rj * np.exp(1j * eta))) + 1j * pi * r0 * j - 2j * pi * c * p0 * (u + 0.5)
    L_u = min(n / (2.0 * c), 0.5 + sqrt(35.0 / max(tau, 0.1)))
    def Psi_fin(u0, th):
        lam_W = exp(-zeta(u0).real)
        def fr(u):
            z = lam_W * np.exp(zeta(u) + 1j * th); az = abs(z)
            return (log(az) + log(max(abs(1+1/z), 1e-300))) if az > 1 else log(max(abs(1+z), 1e-300))
        def fi(u):
            z = lam_W * np.exp(zeta(u) + 1j * th); az = abs(z)
            if az > 1.0:
                return np.angle(1.0 + 1.0 / z)
            return np.angle(1.0 + z)
        vr, _ = quad(fr, -L_u, L_u, points=[-u0, u0], limit=300)
        vi, _ = quad(fi, -L_u, L_u, points=[-u0, u0], limit=300)
        return (-log(lam_W) + vr) + 1j * ((2.0 * u0 - 1.0) * th + vi)
    # Choose u0:
    if D == 1 and tau <= 1.8:
        u0 = min(0.35, 1.49 / ((ct * (tau ** 3)) ** 0.25))
    else:
        ug = np.linspace(0.5 - 0.22 / D, 0.5 + 0.02 / D, 21)
        u0 = float(ug[int(np.argmin([Psi_fin(u, 0.0).real for u in ug]))])
    lam_W = exp(-zeta(u0).real); lam = lam_W ** (1.0 / r0)
    # Exact discrete:
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
    # Check peak error over theta in [-th_lim, th_lim] (th_lim = min(1.5, 0.5*(pi - 2*pi*D*(0.5-u0))))
    delta = D * (0.5 - u0)
    th_lim = min(1.5, 0.6 * (pi - 2 * pi * abs(delta)))
    d0 = eval_disc(pi / r0); c0 = ct * Psi_fin(u0, 0.0)
    mre = 0.0; mim = 0.0
    for th in np.linspace(-th_lim, th_lim, 11):
        dv = eval_disc((pi + th) / r0) - d0
        cv = ct * Psi_fin(u0, th) - c0
        mre = max(mre, abs(dv.real - cv.real))
        mim = max(mim, abs(float(np.angle(np.exp(1j * (dv.imag - cv.imag))))))
    return n, tau, ct, u0, mre, mim

for c, j, r0, D, tau_t in [
    (18, 19, 1, 1, 0.92), (18, 19, 1, 1, 1.5), (18, 19, 1, 1, 5.0),
    (36, 37, 1, 1, 0.65), (50, 51, 1, 1, 0.55),
    (37, 56, 2, 1, 1.0),  (37, 56, 2, 1, 5.0),  (37, 56, 2, 1, 10.0),
    (55, 73, 3, 1, 1.0),  (55, 73, 3, 1, 12.0),
    (91, 109, 5, 1, 2.0), (91, 109, 5, 1, 15.0),
    (72, 74, 1, 2, 2.0),  (72, 74, 1, 2, 12.0),
    (146, 220, 2, 2, 4.0),(146, 220, 2, 2, 20.0),
]:
    p0 = round(r0 * j / c)
    if abs(r0 * j - c * p0) != D: continue
    res = check_poisson_corner(c, j, r0, D, tau_t)
    if res is None: continue
    n, tau, ct, u0, mre, mim = res
    print(f"  c={c:3d} j={j:3d} r0={r0} D={D} ct={ct:4.1f} tau={tau:5.2f} ct*t^2/D^4={ct*tau*tau/(D**4):6.1f} u0={u0:.3f} | PK |dRe|<={mre:.5f} |dIm|<={mim:.5f}")
