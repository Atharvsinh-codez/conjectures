"""
1. TEST THE 1-LINE EXPLICIT RADIUS RULE u_0 = min(0.35, 1.49 / (ct * tau^3)^{1/4})
   FOR D = 1, SMALL tau in [0.1, 2.5] ACROSS ALL c_tilde = 15 .. 1,000,000!!
2. VERIFY THE 1-LINE EXACT Li_2 IBP FORMULA FOR g_1(u_0) TO MACHINE PRECISION
   AND CHECK ITS IVT BRACKET FOR ALL D >= 2, ALL tau > 0!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad
import mpmath as mp

def polylog2_c(w):
    if abs(w) < 0.85:
        s = 0j; p = w
        for k in range(1, 60):
            s += p / (k * k)
            p *= w
            if abs(p) < 1e-16: break
        return s
    return complex(mp.polylog(2, w))

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
    v, _ = quad(f, -L, L, points=[-u0, u0] if u0 > 1e-12 else [0.0], limit=300)
    return -log(lam) + v

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
    v, _ = quad(f, -L, L, points=[-u0, u0] if u0 > 1e-12 else [0.0], limit=300)
    return (2.0 * u0 - 1.0) * th + v

print("=== 1. TEST EXPLICIT RULE u_0 = min(0.35, 1.49 / (ct * tau^3)^{1/4}) FOR D = 1, tau <= 2.5 ===")
for tau in [0.2, 0.5, 1.0, 2.0, 2.5]:
    print(f"\n  D = 1, tau = {tau:.1f}:")
    for ct in [15, 100, 1000, 10000, 100000, 1000000]:
        u0 = min(0.35, 1.49 / ((ct * (tau ** 3)) ** 0.25))
        lam = exp(-tau * (0.25 - u0 * u0))
        p0 = Re_Psi(1, tau, lam, 0.0)
        H = -2.0 * (Re_Psi(1, tau, lam, 0.005) - p0) / (0.005 ** 2)
        g1 = Im_Psi(1, tau, lam, 0.005) / 0.005
        # Adaptive integration over theta in [0, pi]: fine peak [0, min(pi, 6/sqrt(ct*H))] + rest!
        w_pk = min(pi, 6.0 / sqrt(max(ct * H, 0.1)))
        t1 = np.linspace(0.0, w_pk, 200)
        t2 = np.linspace(w_pk, pi, 100) if w_pk < pi else np.array([pi])
        def eval_seg(ts):
            if len(ts) < 2: return 0.0, 0.0
            re_v = np.array([Re_Psi(1, tau, lam, t) - p0 for t in ts])
            im_v = np.array([Im_Psi(1, tau, lam, t) for t in ts])
            vals_re = np.exp(ct * re_v) * np.cos(ct * im_v)
            vals_abs = np.exp(ct * re_v)
            return float(np.trapezoid(vals_re, ts)), float(np.trapezoid(vals_abs, ts))
        r1, a1 = eval_seg(t1); r2, a2 = eval_seg(t2)
        rat = (r1 + r2) / (a1 + a2)
        drop_pi = ct * (p0 - Re_Psi(1, tau, lam, pi))
        lin_mis = ct * g1 * g1 / H
        print(f"    ct={ct:7d}: u0={u0:.4f} ct*u0={ct*u0:8.1f} | ct*g1^2/H={lin_mis:.3f} ct*drop@pi={drop_pi:6.2f} | rat={rat:+.5f}")

print("\n=== 2. VERIFY 1-LINE Li_2 IBP FORMULA FOR g_1(u_0) ===")
for D, tau, delta in [(1, 0.5, 0.15), (1, 1.0, 0.15), (1, 4.0, 0.10), (2, 1.0, 0.13), (2, 8.0, 0.05), (3, 6.0, 0.07), (5, 10.0, 0.04)]:
    u0 = 0.5 - delta / D
    lam = exp(-tau * (0.25 - u0 * u0))
    # Direct numerical g1 via Im_Psi'(0):
    g1_num = Im_Psi(D, tau, lam, 1e-4) / 1e-4
    # Exact 1-line Li_2 IBP formula for g1:
    #   g1 = -(2*del*tau*u0/D)*H0 - Re int L_E(u,0) * tau / (2*(tau*u + i*pi*D)^2) du
    #   Substitute L_E = M_E' / (-2*(tau*u + i*pi*D)):
    #   -int L_E * tau/(2*(tau*u+i*pi*D)^2) = +int M_E' * tau/(4*(tau*u+i*pi*D)^3)
    #   Jump of M_E at +u0 is +(pi^2/6 - 2*pi^2*del^2), at -u0 is -(pi^2/6 - 2*pi^2*del^2)!
    #   IBP boundary: -tau/4 * (pi^2/6 - 2*pi^2*del^2) * [ 1/(tau*u0 + i*pi*D)^3 - 1/(-tau*u0 + i*pi*D)^3 ]
    #   IBP integral: +int M_E(u, 0) * 3*tau^2 / (4*(tau*u + i*pi*D)^4) du!
    H0 = (tau * u0) / (tau*tau*u0*u0 + pi*pi*D*D)
    z1 = tau * u0 + 1j * pi * D; z2 = -tau * u0 + 1j * pi * D
    term1 = -(2 * delta * tau * u0 / D) * H0
    term2 = (-(tau / 4.0) * (pi*pi/6.0 - 2*pi*pi*delta*delta) * (1.0/(z1**3) - 1.0/(z2**3))).real
    L = sqrt(35.0 / max(tau, 1e-4)) + 1.0
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    def rg1_f(u):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2 * pi * D * u + (pi if sgnD < 0 else 0.0)
        z = r * (cos(ph) + 1j * sin(ph))
        ME = polylog2_c(-1.0 / z) if abs(u) < u0 else -polylog2_c(-z)
        return (ME * (3.0 * tau * tau) / (4.0 * ((tau * u + 1j * pi * D) ** 4))).real
    rg1, _ = quad(rg1_f, -L, L, points=[-u0, u0], limit=300)
    g1_id = term1 + term2 + rg1
    print(f"  D={D} tau={tau:4.1f} del={delta:.2f}: num={g1_num:+.8f} | id={g1_id:+.8f} (t1={term1:+.6f} t2={term2:+.6f} r={rg1:+.6f}) | err={abs(g1_num-g1_id):.2e}")
