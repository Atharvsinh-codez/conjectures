"""
CHECK EXACT VALUES OF |R_3_par| / H_0, e_M / H_0, AND e_W3 / H_0 ACROSS ALL OF CASE G3-D (|D| = 1, tau in [1.05, 40.0])!
Also check exact finite-n curvature & saddle bracket across [1.05, 40.0]!
"""
import numpy as np
import mpmath as mp
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

def polylog2_c(w):
    if abs(w) < 0.85:
        s = 0j; p = w
        for k in range(1, 60):
            s += p / (k * k); p *= w
            if abs(p) < 1e-16: break
        return s
    return complex(mp.polylog(2, w))

print("=== CASE G3-D (|D| = 1, tau in [1.05, 40.0]): CHECK |R_3_par|/H0, e_M/H0, e_W3/H0 ===")
print(f"{'tau':>5} {'u0':>5} {'c_min':>5} | {'|R3p|/H':>7} {'e_M/H':>7} {'e_W3/H':>7} | {'TOT/H':>7} | {'|B1+B2|/H':>9}")
cases_D = [
    (1.05, 0.42), (1.3, 0.42), (1.6, 0.42), (2.0, 0.42), (2.5, 0.42),
    (2.5, 0.286), (3.0, 0.335), (4.0, 0.398), (6.0, 0.458), (8.0, 0.482),
    (10.0, 0.494), (14.0, 0.502), (20.0, 0.504), (28.0, 0.503), (40.0, 0.501),
]
for tau, u0 in cases_D:
    # c >= max(10, tau / 0.45, 40 / tau^2) (actually ct >= max(10, 40/tau^2) and c >= tau/0.45!)
    c_min = max(10.0, tau / 0.45, 40.0 / (tau * tau))
    H0 = tau * u0 / (tau * tau * u0 * u0 + pi * pi)
    H_ref = tau / (2 * pi * pi) if (tau <= 2.2 and u0 == 0.42) else H0
    L = sqrt(30.0 / tau) + 1.0
    # 1. Exact max_{th in [0, pi]} |R_3_par(th)| / H_ref:
    r3p_max = 0.0
    for th in [0.0, 0.5, 1.0, 1.5, 2.2, pi]:
        def r3_f(u):
            lr = tau * (u0 * u0 - u * u); ph = pi - 2 * pi * u + th
            z = exp(lr) * (cos(ph) + 1j * sin(ph)) if lr <= 0 else exp(-lr) * (cos(-ph) + 1j * sin(-ph))
            ME = polylog2_c(-z) if lr > 0 else -polylog2_c(-z)
            return ME * (15.0 * (tau ** 3)) / (8.0 * ((tau * u + 1j * pi) ** 6))
        rr, _ = quad(lambda u: r3_f(u).real, -L, L, points=[-u0, u0], limit=200)
        ri, _ = quad(lambda u: r3_f(u).imag, -L, L, points=[-u0, u0], limit=200)
        r3p_max = max(r3p_max, abs(rr + 1j * ri) / H_ref)
    # 2. Exact finite-n perturbation bounds at r0*kap <= 0.45 and c >= c_min:
    #    r0*c*kap^2 <= 0.45 * tau (and <= 42)!
    r0ck2 = min(42.0, 0.45 * tau)
    Del0 = 0.00638 * r0ck2
    #    |dz(u)| = |(zeta_fin(u) - zeta_fin(u0)) - (zeta_par(u) - zeta_par(u0))|
    #    Check exact derivative |dz'(t)|:
    #      Re: |Re(z_fin'(t)) - (-2*tau*t)| <= (0.006 + 0.0104*t^2) * 2*tau*|t|
    #      Im: |Im(z_fin'(t)) - 2*pi*D| = 2*pi*|Del0 - Del2(t)| <= 2*pi*(Del0 + 0.0765*r0ck2*t^2)
    #    Wait: for |t| <= 0.5, Del0 and Del2(t) satisfy int_0^{0.5} (Del0 - Del2(t)) dt == 0!!
    #    So |Im_err(u) - Im_err(u0)| <= min(2*pi*0.00638*r0ck2*0.385, 2*pi*(Del0 + 0.0765*r0ck2*max(u^2,u0^2))*|u-u0|) on |u|<=0.5!
    def eM_f(u):
        x = 0.95 * tau * abs(u * u - u0 * u0)
        if x < 1e-12: return 0.0
        mlog = -log(1.0 - exp(-x)) if x < 30 else exp(-x)
        # Re err:
        dre = (0.006 + 0.0104 * max(u * u, u0 * u0)) * tau * abs(u * u - u0 * u0)
        # Im err:
        dim_slope = 2 * pi * (Del0 + 0.077 * r0ck2 * max(u * u, u0 * u0)) * abs(u - u0)
        dim = min(dim_slope, 2 * pi * Del0 * 0.39) if (abs(u) <= 0.55 and u0 <= 0.55) else dim_slope
        dz = sqrt(dre * dre + dim * dim)
        W3p = 15.0 * (2 * tau) ** 3 / (((2 * tau * u) ** 2 + (2 * pi) ** 2) ** 3)
        return mlog * dz * W3p
    eM, _ = quad(eM_f, -L, L, points=[-u0, u0], limit=200)
    # 3. Exact e_W3 = int |M_E_fin| * |W_3_fin - W_3_par| du:
    #    |M_E_fin| <= (pi^2/6) * exp(-0.95*tau*|u^2 - u0^2|)
    #    |W_3_fin - W_3_par| <= |15*(-2*tau)^3/z_par'^6 - 15*(z_fin'')^3/z_fin'^6| + |10*z'*z''*z''' - z'^2*z''''|/|z_fin'|^6
    #    Check |z_fin'(u) - z_par'(u)| / |z_par'(u)|:
    #      At u: |dz_p(u)| <= sqrt(((0.006 + 0.0104*u^2)*2*tau*|u|)^2 + (2*pi*max(Del0, 0.077*r0ck2*u^2))^2)
    #      |z_par'(u)| = sqrt((2*tau*u)^2 + (2*pi)^2)
    #      |z_fin''(u) - (-2*tau)| <= (0.006 + 0.0104*u^2 + 0.072*2*0.072) ... actually |z_fin''(u) - (-2*tau)| / (2*tau) <= 0.006 + 0.011*u^2 + 0.144*abs(u)*th?
    #      Wait: z_fin''(u) = -2*tau_eff * e^{i eta} / (1 + (1-e^{i eta})*s_j/(1-s_j))^2 -> phase of z_fin'' is ONLY 2*eta/(1-s_j) ~ 2*u*tan(al) <= 0.144*|u|!
    #      WAIT A SECOND!! CHECK arg(z_fin''(u))!!
    #      z_fin''(u) = -2*tau_eff * e^{2i u th_c} / (1 - s_j e^{2i u th_c})^2!
    #      Since arg(1 - s_j e^{i eta}) ~ -s_j/(1-s_j) * eta ~ -2*u*tan(al), arg(z_fin''(u)) ~ +4*u*tan(al) <= 0.288*u!!
    #      THANK GOD WE CHECKED THE PHASE OF z_fin''(u)!!
    #      Let's compute the ACTUAL worst-case |W_3_fin(u) - W_3_par(u)| directly on real finite-n parameters (c=ceil(c_min), kap=tau/c) to see its true size!
    c_ex = max(10, int(np.ceil(c_min)))
    j_ex = c_ex + 1
    kap_ex = tau / c_ex
    n_ex = max(132, 2 * int(round(pi * pi * j_ex * (j_ex + c_ex) / (c_ex * kap_ex))))
    kap_a = 2 * pi * sin(j_ex * pi / n_ex) * sin((j_ex + c_ex) * pi / n_ex) / sin(c_ex * pi / n_ex)
    rj = sin(j_ex * pi / n_ex) / sin((j_ex + c_ex) * pi / n_ex); thc = c_ex * pi / n_ex
    tau_a = c_ex * kap_a
    def eW3_exact(u):
        e = np.exp(2j * u * thc); q1 = 1 - rj * e
        zp = 2 * pi * c_ex * (1j * rj * e / q1) - 2j * pi * c_ex
        z2 = (2 * pi * c_ex) * (2 * thc) * (-rj * e / (q1 ** 2))
        z3 = (2 * pi * c_ex) * ((2 * thc) ** 2) * (-1j * rj * e * (1 + rj * e) / (q1 ** 3))
        z4 = (2 * pi * c_ex) * ((2 * thc) ** 3) * (rj * e * (1 + 4 * rj * e + rj * rj * e * e) / (q1 ** 4))
        W3f = (-z4 * zp ** 2 + 10 * z3 * z2 * zp - 15 * z2 ** 3) / (zp ** 6)
        W3p = 15.0 * (2 * tau_a) ** 3 / ((-2 * tau_a * u + 2j * pi) ** 6)
        return (pi * pi / 6.0) * exp(-0.95 * tau_a * abs(u * u - u0 * u0)) * abs(W3f - W3p)
    eW3, _ = quad(eW3_exact, -L, L, points=[-u0, u0], limit=200)
    # Also check |B1+B2| / H_ref at th = 0:
    z1 = tau * u0 + 1j * pi; z2 = -tau * u0 + 1j * pi; delta = 0.5 - u0
    B1_0 = (1j * tau / 4.0) * (2 * pi * delta / (z1 ** 3) + 2 * pi * delta / (z2 ** 3))
    B2_0 = (3.0 * tau * tau / 8.0) * ((2 * pi * pi * delta * delta - pi * pi / 6.0) * (1 / (z1 ** 5) - 1 / (z2 ** 5)))
    print(f"{tau:5.2f} {u0:5.3f} {c_min:5.1f} | {r3p_max:7.4f} {eM/H_ref:7.4f} {eW3/H_ref:7.4f} | {r3p_max+eM/H_ref+eW3/H_ref:7.4f} | {abs((B1_0+B2_0).real)/H_ref:9.4f}")
