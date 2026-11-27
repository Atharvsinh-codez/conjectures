"""
VERIFY THE 2-LINE DOUBLE-IBP EXACT IDENTITY FOR PIECE G3 (D >= 1):
  Psi''(theta) = -H_0 + B(theta) + R_2(theta),
where:
  H_0      = tau * u_0 / (tau^2 u_0^2 + pi^2 D^2),
  B(theta) = (i*tau/4) * [ (2*pi*del + th)/(tau*u_0 + i*pi*D)^3 + (2*pi*del - th)/(-tau*u_0 + i*pi*D)^3 ],
  R_2(th)  = -int_{-inf}^{inf} L_E(u, th) * (3 * tau^2) / (4 * (tau*u + i*pi*D)^4) du,
  L_E(u,th)= log(1 + z^{-1}) for |u| < u_0,  log(1 + z) for |u| > u_0!

TEST:
1. Machine-precision verification of Psi''(th) = -H_0 + B(th) + R_2(th)!
2. Check 1-line ELEMENTARY TRIANGLE BOUND on |R_2(th)| / H_eff across ALL (D, tau) for |th| <= th_lim!
3. Check the EXACT SAME 1-line IBP formula for g_1(Lambda) = -(2*del*tau*u_0/D)*H_0 - Re int L_E(u,0)*tau/(2*(tau*u+i*pi*D)^2) du!
4. Check GLOBAL DECAY on [th_lim, pi] (remember on [0, pi - 2*pi*del], Re Psi''(th) <= -c * H_eff ALREADY PROVES
   Re Psi(th) - Re Psi(0) <= -(c/2)*H_eff*th^2 ON THE ENTIRE INTERVAL |th| <= pi - 2*pi*del!!
   Only the tiny corner [pi - 2*pi*del, pi] near the valley theta = pi needs a 1-line upper bound!!)!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

def L_E_val(D, tau, u0, th, u):
    lam = exp(-tau * (0.25 - u0 * u0))
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    r = lam * exp(tau * (0.25 - u * u))
    ph = -2.0 * pi * D * u + th + (pi if sgnD < 0 else 0.0)
    z = r * (cos(ph) + 1j * sin(ph))
    if abs(u) < u0:
        return np.log(1.0 + 1.0 / z)
    return np.log(1.0 + z)

def R2_and_tri(D, tau, u0, th):
    L = sqrt(40.0 / max(tau, 1e-4)) + 1.0
    pts = [-u0, u0]
    def f_c(u):
        return -L_E_val(D, tau, u0, th, u) * (3.0 * tau * tau) / (4.0 * ((tau * u + 1j * pi * D) ** 4))
    re_v, _ = quad(lambda u: f_c(u).real, -L, L, points=pts, limit=400, epsabs=1e-11, epsrel=1e-11)
    im_v, _ = quad(lambda u: f_c(u).imag, -L, L, points=pts, limit=400, epsabs=1e-11, epsrel=1e-11)
    # 1-line triangle upper bound: int |L_E(u, th)| * 3*tau^2 / (4 * (tau^2 u^2 + pi^2 D^2)^2) du
    tri_v, _ = quad(lambda u: abs(L_E_val(D, tau, u0, th, u)) * (3.0 * tau * tau) / (4.0 * ((tau*tau*u*u + pi*pi*D*D) ** 2)),
                    -L, L, points=pts, limit=400, epsabs=1e-10, epsrel=1e-10)
    return re_v + 1j * im_v, tri_v

def B_exact(D, tau, u0, th):
    delta = D * (0.5 - u0)
    z1 = (tau * u0 + 1j * pi * D) ** 3
    z2 = (-tau * u0 + 1j * pi * D) ** 3
    return (1j * tau / 4.0) * ((2 * pi * delta + th) / z1 + (2 * pi * delta - th) / z2)

print("=== 1. VERIFY 2-LINE DOUBLE-IBP EXACT IDENTITY & 1-LINE TRIANGLE BOUND ON |R_2| ===")
print(f"{'D':>2} {'tau':>5} {'del':>5} | {'H_eff':>8} {'|R2_ex|/H':>9} {'tri_R2/H':>9} {'|ImB|/H@1':>9} | {'id_err':>8}")
for D, tau, delta in [
    (1, 0.5, 0.15), (1, 1.0, 0.15), (1, 2.0, 0.15), (1, 4.0, 0.10), (1, 8.0, 0.02), (1, 16.0, 0.00),
    (2, 1.0, 0.13), (2, 4.0, 0.11), (2, 8.0, 0.05), (2, 16.0, 0.01), (2, 32.0, 0.00),
    (3, 1.5, 0.08), (3, 6.0, 0.07), (3, 12.0, 0.04), (3, 24.0, 0.01),
    (5, 2.5, 0.05), (5, 10.0, 0.04), (5, 20.0, 0.02), (5, 40.0, 0.00),
]:
    u0 = 0.5 - delta / D
    H0 = (tau * u0) / (tau*tau*u0*u0 + pi*pi*D*D)
    B0 = B_exact(D, tau, u0, 0.0).real
    H_eff = H0 - B0 # exact main curvature = -Re(-H0 + B(0))!
    th = 0.5
    R2, tri = R2_and_tri(D, tau, u0, th)
    # Check against E_int_exact at th = 0.5:
    lam = exp(-tau * (0.25 - u0 * u0))
    L = sqrt(40.0 / max(tau, 1e-4)) + 1.0
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    def e_c(u):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2.0 * pi * D * u + th + (pi if sgnD < 0 else 0.0)
        z = r * (cos(ph) + 1j * sin(ph))
        E = (-1.0 / (1.0 + z)) if abs(u) < u0 else (z / (1.0 + z))
        return E * tau / (2.0 * ((tau * u + 1j * pi * D) ** 2))
    er, _ = quad(lambda u: e_c(u).real, -L, L, points=[-u0, u0], limit=300)
    ei, _ = quad(lambda u: e_c(u).imag, -L, L, points=[-u0, u0], limit=300)
    id_err = abs((er + 1j * ei) - (B_exact(D, tau, u0, th) + R2))
    imB1 = abs(B_exact(D, tau, u0, 1.0).imag) / H_eff
    print(f"{D:2d} {tau:5.1f} {delta:5.2f} | {H_eff:8.5f} {abs(R2)/H_eff:9.4f} {tri/H_eff:9.4f} {imB1:9.4f} | {id_err:8.1e}")
