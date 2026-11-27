"""
CHECK EXACT FOURIER DECOMPOSITION & CORRECTION TERMS OF Psi_{D, tau}(-Lambda e^{i theta})
ACROSS THE ENTIRE (D, tau) PLANE FOR PIECE G3 (D >= 1)!

Specifically check:
1. Exact location of the continuum saddle Lambda_*(D, tau) for D = 1, 2, 3, 5, 10 and tau/D in [0.2, 30]!
   (Including tau/D >= 10 where Lambda_* crosses slightly above 1!)
2. Exact Fourier error ratios:
     err_C(r) := |C_r - C_r^{B2}| / H_0,   err_S(r) := |S_r - S_r^{B1}| / H_0
   and total L1 phase error sum_{r>=1} |S_r - S_r^{B1}| / H_0 and Re error sum_{r>=1} |C_r - C_r^{B2}| / H_0!
3. Also check EXACT FINITE-n vs CONTINUUM difference on actual G3 points (c_tilde/D^2 in [15, 100], r0 in [1..5], r0^2*kap up to 1.0)!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

def fourier_CS_exact(D, tau, u0, r):
    # Exact I_r^- and I_r^+ for lam = exp(-tau*(0.25 - u0^2)) (works for BOTH u0 < 0.5 and u0 >= 0.5!!):
    # Note: delta = D * (0.5 - u0) (can be negative if u0 > 0.5, i.e., lam > 1!).
    delta = D * (0.5 - u0)
    # I_r^- = 2 * int_0^{u0} exp(-2*r*tau*u0*t + r*tau*t^2) cos(2*pi*r*(delta + D*t)) dt
    # I_r^+ = 2 * int_0^{inf} exp(-2*r*tau*u0*t - r*tau*t^2) cos(2*pi*r*(delta - D*t)) dt
    L_out = sqrt(40.0 / max(r * tau, 1e-6)) + 1.0
    Im_val, _ = quad(lambda t: 2.0 * exp(-2*r*tau*u0*t + r*tau*t*t) * cos(2*pi*r*(delta + D*t)),
                     0.0, u0, limit=300, epsabs=1e-11, epsrel=1e-11)
    Ip_val, _ = quad(lambda t: 2.0 * exp(-2*r*tau*u0*t - r*tau*t*t) * cos(2*pi*r*(delta - D*t)),
                     0.0, L_out, limit=300, epsabs=1e-11, epsrel=1e-11)
    sgn = 1.0 if (r % 2 == 1) else -1.0
    C_r = (sgn / r) * (Ip_val + Im_val)
    S_r = (sgn / r) * (Ip_val - Im_val)
    H0 = (tau * u0) / (tau*tau*u0*u0 + pi*pi*D*D)
    C_B2 = (sgn / (r * r)) * 2.0 * cos(2*pi*r*delta) * H0
    # S_B1_saw = (sgn / r) * 2 * sin(2*pi*r*delta) / (pi*D) [the pure sawtooth that cancels (2u0 - 1)*th!]
    # S_B1_full = (sgn / r) * 2 * sin(2*pi*r*delta) * (pi*D / (tau^2*u0^2 + pi^2*D^2))
    S_B1_full = (sgn / (r * r)) * 2.0 * sin(2*pi*r*delta) * (pi * D / (tau*tau*u0*u0 + pi*pi*D*D))
    return C_r, C_B2, S_r, S_B1_full, H0

# Exact g1(u0) = (2*u0 - 1) + sum_{r>=1} r * S_r? Wait: Im Psi = (2*u0 - 1)*th - sum S_r sin(r th)
# So d/dth Im Psi|_{th=0} = (2*u0 - 1) - sum_{r>=1} r * S_r!
# Check directly via g1 = d/d(log lam) Re Psi(-lam)!
def g1_direct(D, tau, u0):
    # d/d(log lam) Re Psi(-lam) = -1 + int_{-inf}^{inf} (z/(1+z)) du  (Re part!)
    lam = exp(-tau * (0.25 - u0 * u0))
    L = sqrt(40.0 / max(tau, 1e-4)) + 1.0
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    def f(u):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2.0 * pi * D * u + (pi if sgnD < 0 else 0.0)
        # Re( z / (1 + z) ) = (r cos ph + r^2) / (1 + 2 r cos ph + r^2)
        if r > 1.0:
            ir = 1.0 / r
            return (1.0 + ir * cos(ph)) / max(1.0 + 2*ir*cos(ph) + ir*ir, 1e-300)
        return (r * cos(ph) + r * r) / max(1.0 + 2*r*cos(ph) + r*r, 1e-300)
    # Only has singularities at u_k where r = 1 and ph = pi; since u0 avoids u_k, smooth!
    val, _ = quad(f, -L, L, points=[-u0, u0], limit=400, epsabs=1e-10, epsrel=1e-10)
    return -1.0 + val

print("=== 1. EXACT FOURIER ERROR RATIOS OF BERNOULLI-B2 APPROXIMATION AT SADDLE u0_* ===")
print(f"{'D':>2} {'tau':>6} {'t/D':>5} | {'del*':>7} {'lam*':>7} {'H0':>8} | {'sum|dC|/H0':>10} {'dC1/H0':>8} {'sum|dS|/H0':>10} {'C1/C1_B2':>8}")
for D in [1, 2, 3, 5]:
    for t_over_D in [0.5, 1.0, 2.0, 4.0, 8.0, 15.0, 25.0]:
        tau = D * t_over_D
        if D == 1 and tau < 1.8: continue
        # Find exact u0_* solving g1(u0_*) = 0 in (0.5 - 0.45/D, 0.5 + 0.1/D)
        lo_u = 0.5 - 0.45 / D
        hi_u = 0.5 + 0.15 / D
        g_lo = g1_direct(D, tau, lo_u)
        g_hi = g1_direct(D, tau, hi_u)
        if g_lo * g_hi > 0:
            # Expand hi_u if needed
            hi_u = 0.5 + 0.4 / D
            g_hi = g1_direct(D, tau, hi_u)
        for _ in range(30):
            mu = 0.5 * (lo_u + hi_u)
            if g1_direct(D, tau, mu) < 0: lo_u = mu
            else: hi_u = mu
        u0 = 0.5 * (lo_u + hi_u)
        del_s = D * (0.5 - u0)
        lam_s = exp(-tau * (0.25 - u0 * u0))
        sum_dC = 0.0; sum_dS = 0.0; dC1 = 0.0; c1_rat = 0.0
        for r in range(1, 16):
            Cr, CB2, Sr, SB1, H0 = fourier_CS_exact(D, tau, u0, r)
            sum_dC += abs(Cr - CB2) / H0
            sum_dS += abs(Sr - SB1) / H0
            if r == 1:
                dC1 = abs(Cr - CB2) / H0
                c1_rat = Cr / CB2 if abs(CB2) > 1e-12 else 0.0
        print(f"{D:2d} {tau:6.2f} {t_over_D:5.1f} | {del_s:+7.4f} {lam_s:7.4f} {H0:8.5f} | {sum_dC:10.4f} {dC1:8.4f} {sum_dS:10.4f} {c1_rat:8.4f}")
