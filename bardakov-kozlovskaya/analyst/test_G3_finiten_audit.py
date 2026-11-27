"""
EXACT FINITE-n STRESS AUDIT OF PIECE G3 (D >= 1 MAJOR ARCS) ACROSS ALL EXTREME CORNERS:
Test exact finite-n F_pure(lam e^{i psi}) vs continuum c_tilde * Psi_{D, tau}(-Lambda e^{i theta})
across:
  - D in {1, 2, 3}
  - r0 in {1, 2, 3, 5}
  - c_tilde / D^2 in {15, 22, 50}
  - tau = c*r0*kap across small, medium, large (keeping r0*kap <= 0.25, r0^3*c*kap^2 <= 20, c*q*kap > 0.85)

Measure:
1. Exact finite-n rat_L1 on |w| = lam!
2. Exact |log F_exact - c_tilde * Psi_{D, tau}| on the peak |theta| <= 0.6*(pi - 2*pi*|del|)!
3. Exact closed-form -log(1 - r) bound on |R_2(th)| / H_eff across ALL (D, tau)!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd
from scipy.integrate import quad

# =============================================================================
# 1. CHECK 1-LINE CLOSED-FORM -log(1 - r) BOUND ON |R_2(th)| / H_eff!
# =============================================================================
print("=== 1. CHECK UNIVERSAL -log(1 - r) TRIANGLE BOUND ON |R_2(th)| / H_eff (INDEPENDENT OF theta!) ===")
print(f"{'D':>2} {'tau':>5} {'del':>5} | {'H_eff':>8} | {'log_int/H':>9} {'cf_bound/H':>10}")
for D, tau, delta in [
    (1, 0.5, 0.15), (1, 1.0, 0.15), (1, 2.0, 0.15), (1, 4.0, 0.10), (1, 8.0, 0.02), (1, 16.0, 0.00),
    (2, 1.0, 0.13), (2, 4.0, 0.11), (2, 8.0, 0.05), (2, 16.0, 0.01), (2, 32.0, 0.00),
    (3, 1.5, 0.08), (3, 6.0, 0.07), (3, 12.0, 0.04), (3, 24.0, 0.01),
    (5, 2.5, 0.05), (5, 10.0, 0.04), (5, 20.0, 0.02), (5, 40.0, 0.00),
]:
    u0 = 0.5 - delta / D
    H0 = (tau * u0) / (tau*tau*u0*u0 + pi*pi*D*D)
    z1 = (tau * u0 + 1j * pi * D) ** 3
    z2 = (-tau * u0 + 1j * pi * D) ** 3
    B0 = ((1j * tau / 4.0) * (2 * pi * delta / z1 + 2 * pi * delta / z2)).real
    H_eff = H0 - B0
    L = sqrt(40.0 / max(tau, 1e-4)) + 1.0
    # Exact integral of (-log(1 - exp(-tau*|u^2 - u0^2|))) * 3*tau^2 / (4*(tau^2 u^2 + pi^2 D^2)^2)
    def f_log(u):
        d = tau * abs(u * u - u0 * u0)
        if d < 1e-12: return 0.0
        v = -log(1.0 - exp(-d)) if d < 35.0 else exp(-d)
        return v * (3.0 * tau * tau) / (4.0 * ((tau*tau*u*u + pi*pi*D*D) ** 2))
    val, _ = quad(f_log, -L, L, points=[-u0, u0], limit=400, epsabs=1e-10, epsrel=1e-10)
    # 1-line closed-form bound:
    # For |u| <= u0: u0^2 - u^2 >= u0*(u0 - |u|), and denom >= (pi*D)^4 (or split at u0/2!)
    # For |u| >= u0: u^2 - u0^2 >= 2*u0*(|u| - u0) AND denom >= (tau^2 u0^2 + pi^2 D^2)^2!
    # Wait: even simpler! On [0, u0]: int_0^{u0} (-log(1 - e^{-tau*u0*t})) dt <= pi^2 / (6*tau*u0),
    #   AND on [u0, inf): int_0^inf (-log(1 - e^{-2*tau*u0*t})) dt = pi^2 / (12*tau*u0)!
    #   Also for small tau, on [0, u0]: int_0^{u0} (-log(1 - e^{-tau*(u0^2 - u^2)})) du <= u0 * (1 + log(1/(tau*u0^2))...)
    # Let's check the simple 2-piece closed form:
    #   Outer (|u| > u0): denom >= (tau^2 u0^2 + pi^2 D^2)^2, int <= 2 * min(pi^2/(12*tau*u0), sqrt(pi)*2.613/(2*sqrt(tau)))
    #   Inner (|u| < u0): denom >= (pi^2 D^2)^2, int <= 2 * min(pi^2/(6*tau*u0), u0*(2 + max(0, log(1/(tau*u0*u0)))))
    I_out = 2.0 * min((pi*pi)/(12*tau*u0), sqrt(pi)*2.613/(2*sqrt(tau))) * (3*tau*tau) / (4*((tau*tau*u0*u0 + pi*pi*D*D)**2))
    I_in  = 2.0 * min((pi*pi)/(6*tau*u0), u0*(2.0 + max(0.0, log(1.0/(tau*u0*u0))))) * (3*tau*tau) / (4*(pi**4)*(D**4))
    cf = I_out + I_in
    print(f"{D:2d} {tau:5.1f} {delta:5.2f} | {H_eff:8.5f} | {val/H_eff:9.4f} {cf/H_eff:10.4f}")

# =============================================================================
# 2. EXACT FINITE-n STRESS AUDIT ON HARDEST PIECE G3 CORNERS!
# =============================================================================
print("\n=== 2. EXACT FINITE-n STRESS AUDIT ON HARDEST G3 CORNERS (c_tilde/D^2 >= 15, r0*kap <= 0.25) ===")
def Re_Psi_c(D, tau, lam, th):
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

def audit_finiten_G3(c, j, r0, D, target_tau):
    target_kap = target_tau / (c * r0)
    n_est = int(round(2.0 * pi * pi * j * (j + c) / (c * target_kap)))
    n = max(132, 2 * ((n_est + 1) // 2))
    if j > n // 2 - (c + 1) // 2 - 1: return None
    kap = 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)
    tau = c * r0 * kap
    ct = c / r0
    # Choose delta: 0.15 for D=1 & tau < 3.0; else exact continuum saddle!
    if D == 1 and tau < 3.0:
        u0 = 0.35
    else:
        # Quick golden section for continuum saddle
        lo_u, hi_u = 0.5 - 0.35 / D, 0.5 + 0.08 / D
        lams_g = [exp(-tau*(0.25 - u*u)) for u in np.linspace(lo_u, hi_u, 25)]
        vs_g = [Re_Psi_c(D, tau, lm, 0.0) for lm in lams_g]
        u0 = float(np.linspace(lo_u, hi_u, 25)[int(np.argmin(vs_g))])
    lam_W = exp(-tau * (0.25 - u0 * u0))
    lam = lam_W ** (1.0 / r0)

    # Exact finite-n roots a_k:
    rj = sin(j * pi / n) / sin((j + c) * pi / n)
    sj = sin(c * pi / n) / sin((j + c) * pi / n)
    M_k = min(n // 2, int(c + sqrt(45.0 * c / max(kap, 1e-8))) + 15)
    if c % 2 == 0:
        ks = np.arange(-M_k, M_k + 1, dtype=np.float64)
        etas = 2.0 * ks * pi / n
    else:
        ks = np.arange(-M_k, M_k, dtype=np.float64)
        etas = (2.0 * ks + 1.0) * pi / n
    la = n * np.log(sj / (1.0 - rj * np.exp(1j * etas))) + (1j * pi if j % 2 == 1 else 0.0)

    # Evaluate EXACT finite-n F_pure(lam e^{i psi}) on 512 grid!
    psis = np.linspace(-pi, pi, 512, endpoint=False)
    rv = la.real + log(lam)
    big = rv > 0.0; sm = ~big
    lg = -c * (log(lam) + 1j * psis) + (1j * pi if (j*(c+1)+c) % 2 == 1 else 0.0)
    if np.any(big):
        rb = rv[big]; ib = la.imag[big]
        lg = lg + np.sum(rb) + 1j*(np.sum(ib) + len(rb)*pi) + len(rb)*1j*psis
        ch = np.exp(-rb - 1j*ib); ep = np.exp(-1j*psis)
        for k0 in range(0, len(ch), 200):
            lg += np.sum(np.log(1.0 - ch[k0:k0+200][None,:]*ep[:,None]), axis=1)
    if np.any(sm):
        rs = rv[sm]; is_ = la.imag[sm]
        ch = np.exp(rs + 1j*is_); ep = np.exp(1j*psis)
        for k0 in range(0, len(ch), 200):
            lg += np.sum(np.log(1.0 - ch[k0:k0+200][None,:]*ep[:,None]), axis=1)
    sh = float(np.max(lg.real))
    F = np.exp(lg - sh)
    rat_exact = float(np.mean(F).real) / float(np.mean(np.abs(F)))

    # Compare Re(lg(pi/r0 + th/r0) - lg(pi/r0)) vs ct*(Re_Psi_c(th) - Re_Psi_c(0)) at th = 0.5, 1.0!
    def eval_single(psi):
        v = -c * (log(lam) + 1j * psi) + (1j * pi if (j*(c+1)+c) % 2 == 1 else 0.0)
        if np.any(big):
            v += np.sum(rv[big]) + 1j*(np.sum(la.imag[big]) + np.sum(big)*(pi + psi))
            v += np.sum(np.log(1.0 - np.exp(-rv[big] - 1j*(la.imag[big] + psi))))
        if np.any(sm):
            v += np.sum(np.log(1.0 - np.exp(rv[sm] + 1j*(la.imag[sm] + psi))))
        return v
    v0 = eval_single(pi / r0)
    v1 = eval_single(pi / r0 + 0.5 / r0)
    d_ex = (v1 - v0).real
    d_ct = ct * (Re_Psi_c(D, tau, lam_W, 0.5) - Re_Psi_c(D, tau, lam_W, 0.0))
    im0 = float(np.angle(np.exp(1j * v0.imag)))
    return n, kap, tau, r0*kap, (r0**3)*c*kap*kap, lam_W, rat_exact, im0, d_ex, d_ct, abs(d_ex - d_ct)

for c, j, r0, D, t_tau in [
    (15, 16, 1, 1, 1.5), (15, 16, 1, 1, 3.0), (20, 21, 1, 1, 4.0), (20, 21, 1, 1, 8.0),
    (30, 31, 1, 1, 1.0), (30, 31, 1, 1, 5.0), (30, 46, 2, 2, 4.0), (45, 61, 3, 3, 6.0),
    (30, 46, 2, 1, 4.0), # Wait: 2*46 - 30*3 = 92 - 90 = 2 (D=2)! Let's do true D=1, r0=2: c=31, j=47 -> 2*47-31*3=1!
    (31, 47, 2, 1, 3.0), (31, 47, 2, 1, 6.0),
    (46, 61, 3, 1, 4.0), (46, 61, 3, 1, 8.0), # 3*61 - 46*4 = 183 - 184 = -1 (|D|=1!)
    (60, 62, 1, 2, 4.0), (60, 62, 1, 2, 10.0),
    (120, 122, 1, 2, 6.0), (135, 138, 1, 3, 8.0),
    (76, 91, 5, 1, 5.0), # 5*91 - 76*6 = 455 - 456 = -1 (|D|=1, ct=15.2!)
]:
    # Verify D:
    p0 = round(r0 * j / c)
    D_act = abs(r0 * j - c * p0)
    if D_act != D: continue
    res = audit_finiten_G3(c, j, r0, D, t_tau)
    if res is None: continue
    n, kap, tau, r0k, r3ck2, lamW, rat, im0, dex, dct, err = res
    print(f"  n={n:6d} c={c:3d} j={j:3d} r0={r0} D={D} ct/D2={c/(r0*D*D):4.1f} tau={tau:4.1f} r0k={r0k:.3f} r3ck2={r3ck2:4.1f} | rat={rat:+.5f} im0={im0:+.3f} | dRe=({dex:+.3f},{dct:+.3f}) err={err:.4f}")
