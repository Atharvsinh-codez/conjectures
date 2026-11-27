"""
DESIGN PIECE G3 (D != 0 MAJOR ARCS):
1. Check universal continuum saddle Lambda_*(D, tau) for D in {1, 2, 3, 5, 10} and tau in [0.2, 50]:
   - Compute exact delta_* = D * (0.5 - u_0(Lambda_*)) in (0, 0.5)!
   - Compute exact H(Lambda_*), global decay ratio c_dec = inf_{th != 0} (Psi(0) - Psi(th)) / (H * sin^2(th/2)),
     and cubic phase ratio K_3 = sup_{|th| <= th_0} |Im Psi(th)| / (H * |th|^3)!
   - Discover the EXACT UNIVERSAL SCALING LAW of (delta_*, H, c_dec, K_3) in terms of (D, tau)!!
2. Also check exact finite-n discrete vs continuum error at Lambda_* across G3 corners
   (including c_tilde / D^2 = 12, 15, 22, 50 and r0 = 1, 2, 3, 5 and r0*kap up to 0.25)!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

def Re_Psi_exact(D, tau, lam, theta=0.0):
    L = sqrt(35.0 / max(tau, 1e-4)) + 1.0
    arg_u0 = 0.25 + log(lam) / tau
    u0 = sqrt(arg_u0) if arg_u0 > 0 else 0.0
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    def integ(u):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2.0 * pi * D * u + theta + (pi if sgnD < 0 else 0.0)
        if r > 1.0:
            ir = 1.0 / r
            return log(r) + 0.5 * log(max(1.0 + 2*ir*cos(ph) + ir*ir, 1e-300))
        return 0.5 * log(max(1.0 + 2*r*cos(ph) + r*r, 1e-300))
    pts = [-u0, u0] if u0 > 1e-12 else [0.0]
    val, _ = quad(integ, -L, L, points=pts, limit=400, epsabs=1e-11, epsrel=1e-11)
    return -log(lam) + val

def Im_Psi_exact(D, tau, lam, theta):
    # Im Psi_{D, tau}(-lam e^{i th}) = -th + int_{|u|<u0} (-2*pi*D*u + th + arg(1 + 1/z)) du + int_{|u|>u0} arg(1 + z) du
    # Wait: the -2*pi*D*u integrates to 0 by oddness! So Im Psi = (2*u0 - 1)*th + int_{|u|<u0} atan2(im(1+1/z), re(1+1/z)) + int_{|u|>u0} atan2(...)!
    L = sqrt(35.0 / max(tau, 1e-4)) + 1.0
    arg_u0 = 0.25 + log(lam) / tau
    u0 = sqrt(arg_u0) if arg_u0 > 0 else 0.0
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    def integ(u):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2.0 * pi * D * u + theta + (pi if sgnD < 0 else 0.0)
        if r > 1.0:
            ir = 1.0 / r
            return np.arctan2(-ir * sin(ph), 1.0 + ir * cos(ph))
        return np.arctan2(r * sin(ph), 1.0 + r * cos(ph))
    pts = [-u0, u0] if u0 > 1e-12 else [0.0]
    val, _ = quad(integ, -L, L, points=pts, limit=400, epsabs=1e-10, epsrel=1e-10)
    return (2.0 * u0 - 1.0) * theta + val

print("=== UNIVERSAL SCALING OF CONTINUUM SADDLE FOR D = 1, 2, 3, 5 ===")
print(f"{'D':>2} {'tau':>6} {'t/D':>6} | {'lam*':>8} {'del*':>7} | {'H':>9} {'H*D':>8} | {'c_dec':>7} {'max|Im|/H':>9} {'|Im|@0.3/(H*t^3)':>16}")
for D in [1, 2, 3, 5]:
    for tau in [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]:
        if tau / D > 20.0 or tau / D < 0.25: continue
        # Largest singularity u_max = 0.5 - 0.5/D -> lam_sing = exp(-tau*(0.25 - u_max^2))
        u_sing = 0.5 - 0.5 / D
        lam_sing = exp(-tau * (0.25 - u_sing * u_sing))
        # Minimize Re_Psi(lam, 0) on (lam_sing, 1.05) via golden section!
        # First bracket on delta in [0.005, 0.495]:
        dels = np.linspace(0.01, 0.49, 49)
        lams = [exp(-tau * (0.25 - (0.5 - d / D) ** 2)) for d in dels]
        vals = [Re_Psi_exact(D, tau, x, 0.0) for x in lams]
        idx = int(np.argmin(vals))
        lo = lams[min(len(lams)-1, idx+1)] if idx < len(lams)-1 else lam_sing * 1.001
        hi = lams[max(0, idx-1)] if idx > 0 else 1.02
        if lo > hi: lo, hi = hi, lo
        gr = (sqrt(5.0) - 1.0) / 2.0
        c1 = hi - gr*(hi - lo); c2 = lo + gr*(hi - lo)
        f1 = Re_Psi_exact(D, tau, c1, 0.0); f2 = Re_Psi_exact(D, tau, c2, 0.0)
        for _ in range(25):
            if f1 < f2: hi = c2; c2 = c1; f2 = f1; c1 = hi - gr*(hi - lo); f1 = Re_Psi_exact(D, tau, c1, 0.0)
            else: lo = c1; c1 = c2; f1 = f2; c2 = lo + gr*(hi - lo); f2 = Re_Psi_exact(D, tau, c2, 0.0)
        lam_s = 0.5 * (lo + hi)
        u0 = sqrt(max(0.0, 0.25 + log(lam_s) / tau))
        del_s = D * (0.5 - u0)
        p0 = Re_Psi_exact(D, tau, lam_s, 0.0)
        p_dt = Re_Psi_exact(D, tau, lam_s, 0.01)
        H = -2.0 * (p_dt - p0) / (0.01 ** 2)
        # Check global decay ratio c_dec over theta in (0, pi]:
        ths = np.linspace(0.05, pi, 32)
        c_dec = min((p0 - Re_Psi_exact(D, tau, lam_s, th)) / (H * (sin(th / 2.0) ** 2)) for th in ths)
        # Check cubic phase at th = 0.3:
        im03 = abs(Im_Psi_exact(D, tau, lam_s, 0.3))
        max_im = max(abs(Im_Psi_exact(D, tau, lam_s, th)) for th in [0.1, 0.3, 0.6, 1.0])
        print(f"{D:2d} {tau:6.2f} {tau/D:6.2f} | {lam_s:8.5f} {del_s:7.4f} | {H:9.5f} {H*D:8.4f} | {c_dec:7.4f} {max_im/H:9.4f} {im03/(H*0.027):16.4f}")
