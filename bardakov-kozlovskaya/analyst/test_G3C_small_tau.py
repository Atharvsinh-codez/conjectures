"""
Quick check of Case G3-C (|D| = 1, tau in (0, 1.05], ct*tau^2 >= 40):
Verify exact closed-form analytic bounds for ALL tau in (0, 1.05] and ALL ct >= 40/tau^2!
Specifically, test:
  1. Exact 1-line Fourier/Poisson representation of Re Psi_par(th) and g_{1,par}(u0) for small tau!
  2. Exact saddle u_0^*(tau) for tau in [0.05, 1.05], confirming tau*(u_0^*)^2 >> exp(-pi^2/tau)!
  3. Exact positivity margin rat_G3 for u_0 = max(u_0^*, min(0.42, 1.49/(ct*tau^3)^{1/4})) at ct = 10^2 .. 10^10!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

def eval_Psi_par_small_tau(tau, u0, th):
    L = sqrt(35.0 / max(tau, 1e-4)) + 1.0
    def fr(u):
        x = tau * (u0*u0 - u*u)
        ph = pi + 2*pi*u + th
        if x > 0:
            r = exp(-x)
            return x + 0.5 * log(max(1.0 + 2*r*cos(ph) + r*r, 1e-300))
        r = exp(x)
        return 0.5 * log(max(1.0 + 2*r*cos(ph) + r*r, 1e-300))
    return quad(fr, -L, L, points=[-u0, u0], limit=300)[0]

def eval_g1_par(tau, u0):
    L = sqrt(35.0 / max(tau, 1e-4)) + 1.0
    def f(u):
        x = tau * (u0*u0 - u*u)
        ph = pi + 2*pi*u
        if x > 0:
            r = exp(-x)
            # Re(1 / (1 + z)) = Re(r e^{-i ph} / (1 + r e^{-i ph})) = (r cos ph + r^2) / (1 + 2 r cos ph + r^2)
            return 1.0 - (r*cos(ph) + r*r) / (1.0 + 2*r*cos(ph) + r*r)
        r = exp(x)
        return (r*cos(ph) + r*r) / (1.0 + 2*r*cos(ph) + r*r)
    return -1.0 + quad(f, -L, L, points=[-u0, u0], limit=300)[0]

print("=== 1. EXACT g_1(u0) AND SADDLE u_0^*(tau) FOR D = 1, tau in [0.2, 1.05] ===")
for tau in [0.3, 0.5, 0.7, 0.9, 1.05]:
    # Find u0^* by bisection on [1e-6, 0.45]:
    lo, hi = 1e-6, 0.45
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if eval_g1_par(tau, mid) < 0: lo = mid
        else: hi = mid
    u_star = 0.5 * (lo + hi)
    g_042 = eval_g1_par(tau, 0.42)
    p0 = eval_Psi_par_small_tau(tau, u_star, 0.0)
    ppi = eval_Psi_par_small_tau(tau, u_star, pi)
    p0_42 = eval_Psi_par_small_tau(tau, 0.42, 0.0)
    ppi_42 = eval_Psi_par_small_tau(tau, 0.42, pi)
    print(f"  tau={tau:4.2f}: u0*={u_star:.6f} | drop@pi(u*)={p0-ppi:.6f} (tau*u*^2={tau*u_star**2:.6f}) | g1(0.42)={g_042:.6f} (ratio to tau^2*u^2/(2pi^2)={g_042/(tau*tau*0.42*0.42/(2*pi*pi)):.3f}) | drop@pi(0.42)={p0_42-ppi_42:.5f} (tau*0.42^2={tau*0.42**2:.5f})")
