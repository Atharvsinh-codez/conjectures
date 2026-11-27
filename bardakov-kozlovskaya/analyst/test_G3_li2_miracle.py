"""
VERIFY THE 3-STEP POLYLOGARITHM IBP LADDER (6TH-POWER DENOMINATOR!) FOR PIECE G3 (ALL D >= 1!!):
  Psi''(theta) = -H_0 + B_1(theta) + B_2(theta) + R_3(theta),
where:
  H_0 = tau * u_0 / (tau^2 u_0^2 + pi^2 D^2),
  B_1(th) = (i*tau/4) * [ alpha_+ / (tau*u_0 + i*pi*D)^3 - alpha_- / (-tau*u_0 + i*pi*D)^3 ],
  B_2(th) = (3*tau^2/8) * [ (alpha_+^2/2 - pi^2/6) / (tau*u_0 + i*pi*D)^5 - (alpha_-^2/2 - pi^2/6) / (-tau*u_0 + i*pi*D)^5 ],
  R_3(th) = int_{-inf}^{inf} M_E(u, th) * (15 * tau^3) / (8 * (tau*u + i*pi*D)^6) du,
with alpha_+ = 2*pi*del + th, alpha_- = -2*pi*del + th  (|alpha_+|, |alpha_-| < pi),
and  |M_E(u, th)| <= Li_2(r) <= (pi^2 / 6) * exp(-tau * |u^2 - u_0^2|)  FOR ALL theta!

Check:
1. Machine-precision verification of R_2(th) == B_2(th) + R_3(th)!
2. 1-line CLOSED-FORM elementary bound on |R_3(th)| / H_eff across ALL D >= 1 and ALL tau > 0!!
   (ZERO log singularity! PURE exponential integral int exp(-tau*|u^2 - u_0^2|) du!)
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad
import mpmath as mp

def polylog2_c(w):
    # Fast exact series for |w| <= 1 (or mpmath for |w| close to 1)
    if abs(w) < 0.85:
        s = 0j; p = w
        for k in range(1, 60):
            s += p / (k * k)
            p *= w
            if abs(p) < 1e-16: break
        return s
    return complex(mp.polylog(2, w))

def check_li2_ladder(D, tau, delta, th):
    u0 = 0.5 - delta / D
    lam = exp(-tau * (0.25 - u0 * u0))
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    ap = 2 * pi * delta + th
    am = -2 * pi * delta + th
    z1 = tau * u0 + 1j * pi * D
    z2 = -tau * u0 + 1j * pi * D
    H0 = (tau * u0) / (tau*tau*u0*u0 + pi*pi*D*D)
    B1 = (1j * tau / 4.0) * (ap / (z1 ** 3) - am / (z2 ** 3))
    B2 = (3.0 * tau * tau / 8.0) * ((0.5*ap*ap - pi*pi/6.0) / (z1 ** 5) - (0.5*am*am - pi*pi/6.0) / (z2 ** 5))

    L = sqrt(35.0 / max(tau, 1e-4)) + 1.0
    def R2_f(u):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2 * pi * D * u + th + (pi if sgnD < 0 else 0.0)
        z = r * (cos(ph) + 1j * sin(ph))
        LE = np.log(1.0 + 1.0 / z) if abs(u) < u0 else np.log(1.0 + z)
        return -LE * (3.0 * tau * tau) / (4.0 * ((tau * u + 1j * pi * D) ** 4))
    def R3_f(u):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2 * pi * D * u + th + (pi if sgnD < 0 else 0.0)
        z = r * (cos(ph) + 1j * sin(ph))
        ME = polylog2_c(-1.0 / z) if abs(u) < u0 else -polylog2_c(-z)
        return ME * (15.0 * (tau ** 3)) / (8.0 * ((tau * u + 1j * pi * D) ** 6))

    r2r, _ = quad(lambda u: R2_f(u).real, -L, L, points=[-u0, u0], limit=300)
    r2i, _ = quad(lambda u: R2_f(u).imag, -L, L, points=[-u0, u0], limit=300)
    r3r, _ = quad(lambda u: R3_f(u).real, -L, L, points=[-u0, u0], limit=300)
    r3i, _ = quad(lambda u: R3_f(u).imag, -L, L, points=[-u0, u0], limit=300)
    R2 = r2r + 1j * r2i
    R3 = r3r + 1j * r3i
    err = abs(R2 - (B2 + R3))

    # 1-LINE ELEMENTARY CLOSED-FORM BOUND ON |R_3(th)| (VALID FOR ALL theta!!):
    #   |M_E(u, th)| <= (pi^2 / 6) * exp(-tau * |u^2 - u_0^2|)
    #   Outer (|u| > u_0): denom >= (tau^2 u_0^2 + pi^2 D^2)^3,
    #     2 * int_{u_0}^inf exp(-tau*(u^2 - u_0^2)) du <= 2 * min(1/(2*tau*u_0), sqrt(pi)/(2*sqrt(tau)))!
    #   Inner (|u| < u_0): split at u_1 = u_0 / 2 (or 0.6 * u_0):
    #     - [u_1, u_0]: denom >= (tau^2 u_1^2 + pi^2 D^2)^3, int <= min(1/(tau*(u_0+u_1)), u_0 - u_1)
    #     - [0, u_1]:   denom >= (pi^2 D^2)^3, int <= u_1 * exp(-tau*(u_0^2 - u_1^2))
    #   Multiply by (pi^2 / 6) * (15 * tau^3 / 8) = (5 * pi^2 / 16) * tau^3!
    pref = (5.0 * pi * pi / 16.0) * (tau ** 3)
    I_out = 2.0 * min(1.0 / (2 * tau * u0), sqrt(pi) / (2 * sqrt(tau))) / ((tau*tau*u0*u0 + pi*pi*D*D) ** 3)
    u1 = 0.6 * u0
    I_in1 = 2.0 * min(1.0 / (tau * (u0 + u1)), u0 - u1) / ((tau*tau*u1*u1 + pi*pi*D*D) ** 3)
    I_in2 = 2.0 * u1 * exp(-tau * (u0*u0 - u1*u1)) / ((pi * D) ** 6)
    cf_R3 = pref * (I_out + I_in1 + I_in2)

    # Exact triangle bound of (pi^2/6)*exp(-tau|u^2-u0^2|):
    tri_R3, _ = quad(lambda u: pref * exp(-tau * abs(u*u - u0*u0)) / ((tau*tau*u*u + pi*pi*D*D) ** 3),
                     -L, L, points=[-u0, u0], limit=300)
    H_tot0 = H0 - B1.real - B2.real # at th = 0? Let's compute H_main at th=0:
    B1_0 = (1j * tau / 4.0) * (2*pi*delta / (z1**3) + 2*pi*delta / (z2**3))
    B2_0 = (3.0 * tau * tau / 8.0) * ((2*pi*pi*delta*delta - pi*pi/6.0) * (1/(z1**5) - 1/(z2**5)))
    H_main0 = H0 - B1_0.real - B2_0.real
    return H_main0, abs(R3)/H_main0, tri_R3/H_main0, cf_R3/H_main0, abs(B2.real)/H_main0, err

print("=== VERIFY 3-STEP Li_2 IBP LADDER (6TH-POWER DENOMINATOR) FOR ALL D >= 1 ===")
print(f"{'D':>2} {'tau':>5} {'del':>5} | {'H_main0':>8} | {'|R3|/H':>7} {'tri3/H':>7} {'CF3/H':>7} {'|ReB2|/H':>8} | {'id_err':>8}")
for D, tau, delta in [
    (1, 0.5, 0.15), (1, 1.0, 0.15), (1, 2.0, 0.15), (1, 4.0, 0.10), (1, 8.0, 0.02), (1, 16.0, 0.00), (1, 32.0, 0.00),
    (2, 1.0, 0.13), (2, 4.0, 0.11), (2, 8.0, 0.05), (2, 16.0, 0.01), (2, 32.0, 0.00),
    (3, 1.5, 0.08), (3, 6.0, 0.07), (3, 12.0, 0.04), (3, 24.0, 0.01),
    (5, 2.5, 0.05), (5, 10.0, 0.04), (5, 20.0, 0.02), (5, 40.0, 0.00),
]:
    Hm, r3, t3, cf3, rb2, err = check_li2_ladder(D, tau, delta, 0.5)
    print(f"{D:2d} {tau:5.1f} {delta:5.2f} | {Hm:8.5f} | {r3:7.4f} {t3:7.4f} {cf3:7.4f} {rb2:8.4f} | {err:8.1e}")
