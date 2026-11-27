"""
CHECK FEASIBILITY OF PIECE G3 AT @algebraist's #66 BOUNDARY:
  ct = c / r_0 >= 10 * D^2,  r_0 * kap <= 0.45,  tau = c * r_0 * kap >= 1.0 * |D|,  r_0^3 * c * kap^2 <= 42!
Check:
  1. D = 1 on tau in [1.0, 25.0]: exact saddle u_0*(tau), H, g_1, drop@pi, Poisson exponent ct*tau*u_0^2!
  2. D = 1 on tau >= 25.0 AND D >= 2 on ALL tau >= 1.0*|D|:
     - Finite-n |R_3|/H_0 with Del0 <= 0.27 (from r0^3*c*kap^2 <= 42)
     - Exact Li_3 IVT saddle bracket del_* in (-0.12/(s^2+1), +0.20)
     - Exact Poisson strip / residue error at ct = 10*D^2 and r0^2*kap <= 2.05/D!
"""
import numpy as np
import mpmath as mp
from mpmath import iv
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

iv.dps = 50
mp.mp.dps = 50
PI = iv.pi

def psi_props(D, tau, u0):
    L = sqrt(35.0 / max(tau, 0.1)) + 1.0
    sgnD = -1.0 if D % 2 == 1 else 1.0
    def R(th):
        def f(u):
            lr = tau * (u0 * u0 - u * u); ph = -2 * pi * D * u + th + (pi if sgnD < 0 else 0)
            if lr > 0: ir = exp(-lr); return lr + 0.5 * log(max(1 + 2 * ir * cos(ph) + ir * ir, 1e-300))
            r = exp(lr); return 0.5 * log(max(1 + 2 * r * cos(ph) + r * r, 1e-300))
        return -tau * (u0 * u0 - 0.25) + quad(f, -L, L, points=[-u0, u0], limit=300)[0]
    def I(th):
        def f(u):
            lr = tau * (u0 * u0 - u * u); ph = -2 * pi * D * u + th + (pi if sgnD < 0 else 0)
            if lr > 0: ir = exp(-lr); return np.arctan2(-ir * sin(ph), 1 + ir * cos(ph))
            r = exp(lr); return np.arctan2(r * sin(ph), 1 + r * cos(ph))
        return (2 * u0 - 1) * th + quad(f, -L, L, points=[-u0, u0], limit=300)[0]
    p0 = R(0.0); ppi = R(pi); H = -2 * (R(0.005) - p0) / 0.005 ** 2; g1 = I(0.005) / 0.005
    return p0 - ppi, H, g1

print("=== 1. D = 1 CONTINUUM SADDLE & EXPONENTS ON tau in [1.0, 25.0] ===")
for tau in [1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0, 15.0, 25.0]:
    ug = np.linspace(0.05, 0.52, 48)
    g1s = [psi_props(1, tau, u)[2] for u in ug]
    u_sad = None
    for k in range(len(ug) - 1):
        if g1s[k] * g1s[k + 1] <= 0:
            a, b = ug[k], ug[k + 1]; ga = g1s[k]
            for _ in range(18):
                m = 0.5 * (a + b); gm = psi_props(1, tau, m)[2]
                if ga * gm <= 0: b = m
                else: a = m; ga = gm
            u_sad = 0.5 * (a + b); break
    u_use = u_sad if u_sad is not None else 0.40
    drop, H, g1 = psi_props(1, tau, u_use)
    us_str = f"{u_sad:.4f}" if u_sad is not None else "NONE"
    print(f"  tau={tau:5.1f}: u_sad={us_str:>7} u_use={u_use:.4f} | H={H:.5f} g1={g1:+.5f} drop@pi={drop:.5f} | @ct=10: 10*H={10*H:.3f} 10*drop={10*drop:.3f} 10*tau*u0^2={10*tau*u_use**2:.2f}")

print("\n=== 2. CHECK D = 1 SMALL tau in [1.0, 2.5] WITH u_0 = max(u_sad, min(0.42, 1.49/(ct*tau^3)^{1/4})) AT ct = 10..10^5 ===")
for tau in [1.0, 1.2, 1.5, 2.0]:
    for ct in [10, 15, 25, 100, 1000, 10000]:
        u0 = min(0.42, 1.49 / ((ct * (tau ** 3)) ** 0.25))
        drop, H, g1 = psi_props(1, tau, u0)
        yp = tau * u0 * u0 / (2 * pi)
        pois_ex = 2 * pi * ct * yp
        print(f"  tau={tau:.1f} ct={ct:5d}: u0={u0:.4f} | ct*H={ct*H:6.2f} ct*g1^2/H={ct*g1*g1/H:.3f} ct*drop@pi={ct*drop:6.2f} | 2pi*ct*yp={pois_ex:5.2f} (res={-log(1-exp(-pois_ex)):.4f})")

print("\n=== 3. CHECK D >= 2 AT ct = 10*D^2, tau = 1.0*D (ct*tau = 10*D^3 -> 80 at D=2!) ===")
for D in [2, 3]:
    tau = 1.0 * D; ct = 10.0 * D * D
    ug = np.linspace(0.5 - 0.20 / D, 0.5 + 0.02 / D, 23)
    g1s = [psi_props(D, tau, u)[2] for u in ug]
    u_sad = 0.5 - 0.13 / D
    for k in range(len(ug) - 1):
        if g1s[k] * g1s[k + 1] <= 0:
            u_sad = 0.5 * (ug[k] + ug[k + 1]); break
    drop, H, g1 = psi_props(D, tau, u_sad)
    # Closest pole at th = 0: height ~ (pi - 2*pi*del) * H / 2 (and at u_r ~ u0 - 1/(2D), height ~ tau*(u0^2 - u_r^2)/(2*pi*D))
    del_s = D * (0.5 - u_sad)
    yp1 = (pi - 2 * pi * abs(del_s)) * H / 2
    ur_near = max(0.0, u_sad - (0.5 - abs(del_s)) / D)
    yp2 = tau * (u_sad ** 2 - ur_near ** 2) / (2 * pi * D)
    yp = min(yp1, yp2)
    print(f"  D={D} ct={ct:4.0f} tau={tau:4.1f}: u_sad={u_sad:.4f} del={del_s:.4f} | ct*H={ct*H:.3f} ct*drop@pi={ct*drop:.2f} | 2pi*ct*yp={2*pi*ct*yp:.2f}")
