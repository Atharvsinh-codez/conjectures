"""
TEST PIECE G3 (D != 0 MAJOR ARCS):
1. Check D = 1 for SMALL tau in [0.2, 2.5]:
   Why did design_G3_Dne0.py have del* = 0.48 at D = 1, tau = 0.50?
   Test evaluating D = 1 at BOTH:
     (a) Exact saddle lam_* (checking actual integral rat_L1 for c_tilde >= 15!),
     (b) Fixed delta = 0.25 (u_0 = 0.25, lam = exp(-3*tau/16)) where pi - 2*pi*delta = pi/2!
2. Test EXACT CLOSED-FORM BERNOULLI-B2 FOURIER BOUND on Psi_{D, tau}(-lam e^{i theta})!
   Check the exact error |C_r - C_r^main| and |S_r - S_r^main| for ALL r >= 1 across D = 1, 2, 3, 5!
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
    pts = [-u0, u0] if u0 > 1e-12 else [0.0]
    val, _ = quad(f, -L, L, points=pts, limit=400, epsabs=1e-10, epsrel=1e-10)
    return -log(lam) + val

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
    pts = [-u0, u0] if u0 > 1e-12 else [0.0]
    val, _ = quad(f, -L, L, points=pts, limit=400, epsabs=1e-10, epsrel=1e-10)
    return (2.0 * u0 - 1.0) * th + val

print("=== 1. TEST D = 1 SMALL tau IN [0.3, 3.0] AT FIXED delta = 0.20, 0.25, 0.30, 0.35 ===")
for tau in [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0]:
    print(f"\n  D = 1, tau = {tau:.2f} (e^(-tau/4) = {exp(-tau/4):.4f}):")
    for delta in [0.15, 0.20, 0.25, 0.30, 0.35]:
        u0 = 0.5 - delta
        lam = exp(-tau * (0.25 - u0 * u0))
        p0 = Re_Psi(1, tau, lam, 0.0)
        p_dt = Re_Psi(1, tau, lam, 0.01)
        H = -2.0 * (p_dt - p0) / (0.01 ** 2)
        im_dt = Im_Psi(1, tau, lam, 0.01) / 0.01 # linear phase slope g1!
        ths = np.linspace(0.05, pi, 24)
        c_dec = min((p0 - Re_Psi(1, tau, lam, th)) / (H * (sin(th/2)**2)) for th in ths) if H > 0 else -1.0
        # Check actual continuum integral rat for c_tilde = 15, 25, 50:
        th_grid = np.linspace(-pi, pi, 128, endpoint=False)
        re_g = np.array([Re_Psi(1, tau, lam, t) - p0 for t in th_grid])
        im_g = np.array([Im_Psi(1, tau, lam, t) for t in th_grid])
        rats = []
        for ct in [15, 25, 50]:
            F = np.exp(ct * (re_g + 1j * im_g))
            rats.append(float(np.mean(F).real) / float(np.mean(np.abs(F))))
        print(f"    del={delta:.2f} (lam={lam:.4f}): H={H:+.5f} g1={im_dt:+.5f} c_dec={c_dec:+.3f} | rat(ct=15,25,50)=({rats[0]:+.4f},{rats[1]:+.4f},{rats[2]:+.4f})")
