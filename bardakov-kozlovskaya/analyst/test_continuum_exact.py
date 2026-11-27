"""
EXACT (pole-safe!) evaluation of the universal continuum function:
  Psi_{D, tau}(W) = -log(W) + int_{-inf}^{inf} log(1 - W * exp(-2*pi*i*D*x - tau*x*(1+x))) dx
and its derivatives g1 = W * Psi'(W) and H = -d^2/dtheta^2 Re Psi(-lam e^{i theta})|_{theta=0}
using EXACT Gaussian-Poisson / exact branch-free integration!

Key trick to evaluate Psi_{D, tau}(-lam * e^{i theta}) to 1e-13 accuracy with ZERO pole errors:
  Set x = -1/2 + u, so z(u) = (-1)^D * lam * exp(tau/4 - tau*u^2 - 2*pi*i*D*u + i*theta).
  We want I(lam, theta) = int_{-inf}^{inf} log(1 + z(u)) du.
  When |z(u)| > 1 (i.e. |u| < u_0 where tau*u_0^2 = tau/4 + log(lam)):
    log(1 + z(u)) = log(z(u)) + log(1 + 1/z(u))
    Here Re(log(z(u))) = log(lam) + tau*(1/4 - u^2).
    And |1/z(u)| < 1 everywhere on (-u_0, u_0) EXCEPT at the endpoints +-u_0!
  Even cleaner for Re Psi(-lam * e^{i theta}) = -log(lam) + int_{-inf}^{inf} (1/2) log(|1 + z(u)|^2) du:
    Use the EXACT identity:
      int_{-inf}^{inf} log|1 + z(u)| du = int_{|u|<u_0} (log(lam) + tau*(1/4 - u^2)) du
                                          + Re int_{|u|<u_0} log(1 + 1/z(u)) du
                                          + Re int_{|u|>u_0} log(1 + z(u)) du!
  Look! Even if |z(u)| = 1 at u = +-u_0, log|1 + z| is ONLY logarithmically singular (integrable!)
  and away from u = +-u_0, |1/z| < 1 (for |u|<u_0) and |z| < 1 (for |u|>u_0) STRICTLY —
  so there are ZERO interior poles anywhere on (-u_0, u_0)!!
  WAIT A SECOND!! WHY DID test_universal_psi.py HIT A POLE AT u = +-1/4 FOR D = 2, tau = 0.80?!
  BECAUSE u_0 WAS EQUAL TO 1/4 WHEN lam = exp(-3*tau/16)!!
  So the endpoint +-u_0 ITSELF hit z(+-u_0) = -1!!
  Let's compute the EXACT continuous Re Psi(-lam) as a function of lam for D = 1, 2, 3 and tau = 0.2..10
  using adaptive scipy.integrate.quad (splitting at +-u_0)!
"""
import numpy as np
from scipy.integrate import quad

def Re_Psi_exact(D, tau, lam, theta=0.0):
    # Compute Re Psi_{D, tau}(-lam * e^{i theta}) = -log(lam) + int_{-inf}^{inf} log|1 + z(u)| du
    # where z(u) = (-1)^D * lam * exp(tau/4 - tau*u^2 - 2*pi*i*D*u + i*theta).
    # Note: log|1 + z(u)| is EVEN in u when theta = 0!
    # For general theta, integrate over [-L, L] splitting at -u0 and +u0!
    L = np.sqrt(35.0 / max(tau, 1e-4)) + 1.0
    arg_u0 = 0.25 + np.log(lam) / tau
    u0 = np.sqrt(arg_u0) if arg_u0 > 0 else 0.0
    sgnD = -1.0 if (D % 2 == 1) else 1.0

    def integrand(u):
        r = lam * np.exp(tau * (0.25 - u * u))
        ph = -2.0 * np.pi * D * u + theta
        if sgnD < 0:
            ph += np.pi
        # log|1 + r e^{i ph}| = 0.5 * log(1 + 2*r*cos(ph) + r^2)
        if r > 1.0:
            inv_r = 1.0 / r
            v = 1.0 + 2.0 * inv_r * np.cos(ph) + inv_r * inv_r
            return np.log(r) + 0.5 * np.log(max(v, 1e-300))
        else:
            v = 1.0 + 2.0 * r * np.cos(ph) + r * r
            return 0.5 * np.log(max(v, 1e-300))

    pts = [-u0, u0] if u0 > 1e-12 else [0.0]
    # Note: since integrand oscillates with freq D, limit subintervals
    val, _ = quad(integrand, -L, L, points=pts, limit=400, epsabs=1e-11, epsrel=1e-11)
    return -np.log(lam) + val

def find_min_lam_exact(D, tau):
    # Minimize Re_Psi_exact(D, tau, lam, theta=0.0) over lam in [0.1, 3.0]!
    # Note: at the saddle point W = -lam_*, d/d(log lam) Re Psi(-lam) = 0!
    # Wait: is it a MINIMUM of Re Psi(-lam) in lam? YES! (Saddle point: min along real axis, max along circle!)
    em = np.exp(-tau / 4.0)
    # Scan grid first to see shape & singularities!
    lams = np.linspace(max(0.05, em * 0.5), 2.0, 80)
    vals = [Re_Psi_exact(D, tau, x, 0.0) for x in lams]
    idx = int(np.argmin(vals))
    lo = lams[max(0, idx - 2)]
    hi = lams[min(len(lams) - 1, idx + 2)]
    gr = (np.sqrt(5.0) - 1.0) / 2.0
    c1 = hi - gr * (hi - lo); c2 = lo + gr * (hi - lo)
    f1 = Re_Psi_exact(D, tau, c1, 0.0); f2 = Re_Psi_exact(D, tau, c2, 0.0)
    for _ in range(30):
        if f1 < f2:
            hi = c2; c2 = c1; f2 = f1; c1 = hi - gr * (hi - lo)
            f1 = Re_Psi_exact(D, tau, c1, 0.0)
        else:
            lo = c1; c1 = c2; f1 = f2; c2 = lo + gr * (hi - lo)
            f2 = Re_Psi_exact(D, tau, c2, 0.0)
    lam_s = 0.5 * (lo + hi)
    # Check global maximum in theta on [-pi, pi] at lam = lam_s!
    p0 = Re_Psi_exact(D, tau, lam_s, 0.0)
    ths = np.linspace(0.0, np.pi, 37)
    pths = [Re_Psi_exact(D, tau, lam_s, th) for th in ths]
    max_diff = max(pths[1:]) - p0 # should be < 0 if theta=0 is strict global max!
    # Curvature via finite diff at small dth = 0.01
    p_dt = Re_Psi_exact(D, tau, lam_s, 0.01)
    H = -2.0 * (p_dt - p0) / (0.01 ** 2)
    # Also compute g1 at lam = 1.0 via finite diff in log(lam)
    g1_at_1 = (Re_Psi_exact(D, tau, np.exp(1e-4), 0.0) - Re_Psi_exact(D, tau, np.exp(-1e-4), 0.0)) / 2e-4
    # Singularities lam_k = exp(-tau*(0.25 - u_k^2)) for |u_k| < 0.5
    sings = []
    for k in range(-D, D + 1):
        uk = (k + (0.0 if D % 2 == 1 else 0.5)) / D
        if abs(uk) < 0.5:
            sings.append(np.exp(-tau * (0.25 - uk * uk)))
    max_sing = max(sings) if sings else 0.0
    return lam_s, max_sing, g1_at_1, H, max_diff

print("=== EXACT POLE-SAFE CONTINUUM SADDLE ANALYSIS FOR D = 1, 2, 3 ===")
for D in [1, 2, 3]:
    print(f"\n--- D = {D} ---")
    print(f"{'tau':>5} {'maxSing':>8} {'lam_*':>8} {'g1@1':>10} {'H(lam*)':>9} {'max_th-p0':>10}")
    for tau in [0.4, 0.8, 1.5, 2.0, 3.0, 5.0, 10.0]:
        ls, ms, g1, H, md = find_min_lam_exact(D, tau)
        print(f"{tau:5.2f} {ms:8.5f} {ls:8.5f} {g1:10.6f} {H:9.4f} {md:10.6f}")
