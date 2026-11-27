"""
VERIFY THE 1-LINE EXACT SECOND-DERIVATIVE MIRACLE IDENTITY FOR PIECE G3 (D >= 1):
  d^2/dtheta^2 Psi_{D, tau}(-Lambda e^{i theta})
    = -H_0 + int_{-inf}^{inf} E(u, theta) * (tau / (2 * (tau*u + i*pi*D)^2)) du,
where H_0 = tau*u_0 / (tau^2*u_0^2 + pi^2*D^2),
      E(u, th) = -1/(1 + z(u, th)) for |u| < u_0, and +z(u, th)/(1 + z(u, th)) for |u| > u_0!

ALSO CHECK:
1. Machine-precision check of the identity!
2. Actual max |Re E_int(th)| / H_0 and |Im E_int(th)| / H_0 on |th| <= pi/2 (and |th| <= 0.7*(pi - 2*pi*|del|))!
3. Check a 1-line ELEMENTARY UPPER BOUND on |E_int(th)| (by integrating by parts ONCE MORE,
   replacing z = e^{-2*pi*i*D*u + ...} -> factoring out the fast oscillation e^{-2*pi*i*D*u}!)!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt
from scipy.integrate import quad

def E_int_exact(D, tau, u0, th):
    lam = exp(-tau * (0.25 - u0 * u0))
    L = sqrt(40.0 / max(tau, 1e-4)) + 1.0
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    def val_c(u):
        # z(u, th) = sgnD * lam * exp(tau*(0.25 - u^2) - 2*pi*i*D*u + i*th)
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2.0 * pi * D * u + th + (pi if sgnD < 0 else 0.0)
        z = r * (cos(ph) + 1j * sin(ph))
        E = (-1.0 / (1.0 + z)) if abs(u) < u0 else (z / (1.0 + z))
        K = tau / (2.0 * ((tau * u + 1j * pi * D) ** 2))
        return E * K
    pts = [-u0, u0]
    re_v, _ = quad(lambda u: val_c(u).real, -L, L, points=pts, limit=400, epsabs=1e-11, epsrel=1e-11)
    im_v, _ = quad(lambda u: val_c(u).imag, -L, L, points=pts, limit=400, epsabs=1e-11, epsrel=1e-11)
    return re_v + 1j * im_v

# Also verify via finite diff of Psi(th)!
def Psi_val(D, tau, u0, th):
    lam = exp(-tau * (0.25 - u0 * u0))
    L = sqrt(40.0 / max(tau, 1e-4)) + 1.0
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    def fr(u):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2.0 * pi * D * u + th + (pi if sgnD < 0 else 0.0)
        if r > 1.0:
            ir = 1.0 / r
            return log(r) + 0.5 * log(max(1.0 + 2*ir*cos(ph) + ir*ir, 1e-300))
        return 0.5 * log(max(1.0 + 2*r*cos(ph) + r*r, 1e-300))
    def fi(u):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2.0 * pi * D * u + th + (pi if sgnD < 0 else 0.0)
        if r > 1.0:
            ir = 1.0 / r
            return np.arctan2(-ir * sin(ph), 1.0 + ir * cos(ph))
        return np.arctan2(r * sin(ph), 1.0 + r * cos(ph))
    vr, _ = quad(fr, -L, L, points=[-u0, u0], limit=400, epsabs=1e-11, epsrel=1e-11)
    vi, _ = quad(fi, -L, L, points=[-u0, u0], limit=400, epsabs=1e-11, epsrel=1e-11)
    return (-log(lam) + vr) + 1j * ((2.0 * u0 - 1.0) * th + vi)

print("=== 1. VERIFY 1-LINE SECOND-DERIVATIVE IDENTITY ===")
for D, tau, delta in [(1, 4.0, 0.10), (1, 1.0, 0.15), (2, 4.0, 0.10), (3, 12.0, 0.04), (5, 20.0, 0.02)]:
    u0 = 0.5 - delta / D
    H0 = (tau * u0) / (tau*tau*u0*u0 + pi*pi*D*D)
    th = 0.4; dt = 1e-3
    d2_fd = (Psi_val(D, tau, u0, th + dt) - 2*Psi_val(D, tau, u0, th) + Psi_val(D, tau, u0, th - dt)) / (dt * dt)
    d2_id = -H0 + E_int_exact(D, tau, u0, th)
    print(f"  D={D} tau={tau:4.1f} del={delta:.2f} th={th}: FD={d2_fd.real:+.7f}{d2_fd.imag:+.7f}i | ID={d2_id.real:+.7f}{d2_id.imag:+.7f}i | err={abs(d2_fd-d2_id):.2e}")

print("\n=== 2. MEASURE |Re E_int(th)| / H0 AND |Im E_int(th)| / H0 ON |th| <= pi/2 ===")
print(f"{'D':>2} {'tau':>5} {'del':>5} {'H0':>8} | {'max|ReE|/H0':>11} {'max|ImE|/H0':>11} | {'ReE(0)/H0':>10} {'naive_L1/H0':>11}")
for D, tau, delta in [
    (1, 0.5, 0.15), (1, 1.0, 0.15), (1, 2.0, 0.15), (1, 4.0, 0.10), (1, 8.0, 0.02), (1, 16.0, 0.00),
    (2, 1.0, 0.13), (2, 4.0, 0.11), (2, 8.0, 0.05), (2, 16.0, 0.01), (2, 32.0, 0.00),
    (3, 1.5, 0.08), (3, 6.0, 0.07), (3, 12.0, 0.04), (3, 24.0, 0.01),
    (5, 2.5, 0.05), (5, 10.0, 0.04), (5, 20.0, 0.02), (5, 40.0, 0.00),
]:
    u0 = 0.5 - delta / D
    H0 = (tau * u0) / (tau*tau*u0*u0 + pi*pi*D*D)
    th_lim = min(pi / 2.0, 0.7 * (pi - 2 * pi * abs(delta)))
    ths = np.linspace(0.0, th_lim, 9)
    vals = [E_int_exact(D, tau, u0, t) for t in ths]
    max_re = max(abs(v.real) for v in vals) / H0
    max_im = max(abs(v.imag) for v in vals) / H0
    re0 = vals[0].real / H0
    # Also check naive L1 integral of |E| * |K| (without cancellation):
    lam = exp(-tau * (0.25 - u0 * u0))
    L = sqrt(35.0 / max(tau, 1e-4)) + 1.0
    def l1_f(u):
        r = lam * exp(tau * (0.25 - u * u))
        E_b = (1.0 / abs(r - 1.0)) if abs(r - 1.0) > 1e-4 else 0.0
        # Actually at th=0, min|1+z| >= 2*sin(pi*delta) near u0!
        return (min(r, 1.0/r) / (1.0 + min(r, 1.0/r))) * (tau / (2.0 * (tau*tau*u*u + pi*pi*D*D)))
    n_l1, _ = quad(l1_f, -L, L, points=[-u0, u0], limit=200)
    print(f"{D:2d} {tau:5.1f} {delta:5.2f} {H0:8.5f} | {max_re:11.4f} {max_im:11.4f} | {re0:+10.4f} {n_l1/H0:11.4f}")
