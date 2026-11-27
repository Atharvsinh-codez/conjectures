"""
CHECK TWO CRITICAL QUANTITIES FOR LEMMA R2.5 (MAJOR ARCS):
1. For q = 1 (c | j) with c*kap^2 < 42, kap < 2.5:
   What is the ACTUAL maximum perturbation of log F(-lam_* e^{i theta}) caused by the finite-n
   imaginary parts arg(a_m) != 0 (compared to the real-positive-root polynomial prod(1 + |a_m| lam_* e^{i theta}))?!
2. For Major Arcs with r_0 >= 1 and D != 0 (e.g. (200, 201, kap=0.01), (200, 267, kap=0.01),
   (200, 299, kap=0.01), (200, 323, kap=0.001), (50, 51, kap=0.10), (500, 501, kap=0.005)):
   Compute the EXACT difference between:
     (a) Exact finite-n log F_pure(lam e^{i psi})
     (b) Continuum approximation (c / r_0) * Psi_{D, tau}(lam^{r_0} e^{i r_0 psi})
   at the saddle peak psi = pi/r_0 (and theta = 0, 0.05, 0.1, 0.2, 0.5)!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd
from scipy.integrate import quad

# =============================================================================
# 1. EXACT FINITE-n PHASE PERTURBATION FOR q = 1 (c | j)
# =============================================================================
print("=== 1. EXACT FINITE-n VS POSITIVE-MODULUS INTEGRAL FOR q = 1 (c*kap^2 < 42) ===")
max_phase_pert = 0.0; worst_q1 = None
for n in [132, 320, 1000, 5000]:
    for c in [2, 3, 4, 6, 10, 20, 40, 80]:
        if c > 0.065 * n: continue
        for mult in [2, 3, 4, 6]:
            j = mult * c
            if j > (n - c) // 2: break
            kap = 2*pi*sin(j*pi/n)*sin((j+c)*pi/n)/sin(c*pi/n)
            if kap >= 2.5 or c * kap * kap >= 42.0 or c * kap <= 0.16: continue
            rj = sin(j*pi/n)/sin((j+c)*pi/n)
            sj = sin(c*pi/n)/sin((j+c)*pi/n)
            if c % 2 == 1:
                ths = ((2*np.arange(n) + 1) * pi / n + pi) % (2*pi) - pi
            else:
                ths = (2*np.arange(n) * pi / n + pi) % (2*pi) - pi
            # a_k = ((-1)^j * sj^n) / (1 + rj e^{i th})^n
            la = n * np.log(sj / (1.0 + rj * np.exp(1j * ths))) + (1j*pi if j % 2 == 1 else 0.0)
            # Solve saddle for POSITIVE moduli |a_k| = exp(la.real):
            lo, hi = -10.0, 10.0
            for _ in range(40):
                vm = 0.5*(lo + hi)
                z = np.clip(la.real + vm, -60.0, 60.0)
                if float(np.sum(1.0 / (1.0 + np.exp(-z)))) > c: hi = vm
                else: lo = vm
            v = 0.5*(lo + hi)
            lam = exp(v)
            # Compute H_pos and max_k |arg(a_k)| weighted by |a_k|lam/(1+|a_k|lam)^2!
            z_re = np.clip(la.real + v, -60.0, 60.0)
            p = 1.0 / (1.0 + np.exp(-z_re))
            H_pos = float(np.sum(p * (1.0 - p)))
            # Compute EXACT ratio F_exact(-lam e^{i theta}) / F_pos(-lam e^{i theta}) at theta = 0!
            # For each pair (a_k, conj(a_k)), (1 + a_k lam) * (1 + conj(a_k) lam) = (1 + |a_k|lam)^2 - 2|a_k|lam(1 - cos(arg a_k))!
            # So log(F_exact(-lam) / F_pos(-lam)) = sum_k (1/2) log(1 - 4*p_k*(1-p_k)*sin^2(arg_k / 2))!
            angs = (la.imag + pi) % (2*pi) - pi
            # Check worst weighted phase: sum_k p_k*(1-p_k) * |ang_k|  (which bounds |g1_exact - g1_pos| AND |Delta H|!)
            w_ang = float(np.sum(p * (1.0 - p) * np.abs(angs)))
            # Check EXACT integral rat on |w| = lam!
            M_psi = 256
            psis = np.linspace(-pi, pi, M_psi, endpoint=False)
            # Exact log F(-lam e^{i psi}) for big/small roots:
            big = (la.real + v) > 0.0
            log_F = -c * (v + 1j * (pi + psis)) + (1j * pi if (j*(c+1)+c)%2==1 else 0.0)
            # Wait: check exact F(-lam e^{i psi}) = (-1)^{j(c+1)+c} (-lam e^{i psi})^{-c} prod (1 + a_k lam e^{i psi})
            if np.any(big):
                rb = la.real[big] + v; ib = la.imag[big]
                log_F = log_F + np.sum(rb) + 1j * np.sum(ib) + len(rb) * 1j * psis
                inv_u = np.exp(-rb - 1j * ib)
                epm = np.exp(-1j * psis)
                keep = np.abs(inv_u) > 1e-14
                if np.any(keep):
                    log_F += np.sum(np.log(1.0 + inv_u[keep][None, :] * epm[:, None]), axis=1)
            sm = ~big
            if np.any(sm):
                rs = la.real[sm] + v; is_ = la.imag[sm]
                u = np.exp(rs + 1j * is_)
                epp = np.exp(1j * psis)
                keep = np.abs(u) > 1e-14
                if np.any(keep):
                    log_F += np.sum(np.log(1.0 + u[keep][None, :] * epp[:, None]), axis=1)
            sh = float(np.max(log_F.real))
            F = np.exp(log_F - sh)
            rat = float(np.mean(F).real) / float(np.mean(np.abs(F)))
            if w_ang / sqrt(H_pos) > max_phase_pert:
                max_phase_pert = w_ang / sqrt(H_pos)
                worst_q1 = (n, c, j, kap, c*kap*kap, H_pos, w_ang, w_ang/sqrt(H_pos), rat)
print(f"  Worst weighted finite-n phase error sum p(1-p)|ang|/sqrt(H) = {max_phase_pert:.5f}")
print(f"  at (n,c,j,kap,ck2,H,w_ang,w_ang/sH,rat) = {worst_q1}")

# =============================================================================
# 2. CHECK EXACT DISCRETE-VS-CONTINUUM ERROR FOR MAJOR ARCS (r_0 >= 1, D != 0)
# =============================================================================
print("\n=== 2. DISCRETE FINITE-n VS CONTINUUM (c/r0)*Psi_{D, tau}(W) ON MAJOR ARCS ===")
def Re_Psi_cont(D, tau, lam_W, theta_W):
    L = sqrt(35.0 / max(tau, 1e-4)) + 1.0
    arg_u0 = 0.25 + log(lam_W) / tau
    u0 = sqrt(arg_u0) if arg_u0 > 0 else 0.0
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    def integ(u):
        r = lam_W * exp(tau * (0.25 - u * u))
        ph = -2.0 * pi * D * u + theta_W + (pi if sgnD < 0 else 0.0)
        if r > 1.0:
            ir = 1.0 / r
            return log(r) + 0.5 * log(max(1.0 + 2*ir*cos(ph) + ir*ir, 1e-300))
        return 0.5 * log(max(1.0 + 2*r*cos(ph) + r*r, 1e-300))
    pts = [-u0, u0] if u0 > 1e-12 else [0.0]
    val, _ = quad(integ, -L, L, points=pts, limit=300, epsabs=1e-10, epsrel=1e-10)
    return -log(lam_W) + val

for c, j, kap, r0, D, lam in [
    (200, 201, 0.01, 1, 1, 0.6705),
    (200, 267, 0.01, 3, 1, 0.9238),
    (200, 299, 0.01, 2, -2, 0.9071),
    (200, 323, 0.001, 13, -1, 0.968), # lam^13 ~ 0.655
]:
    alpha = kap / c
    tau = c * r0 * kap
    c_tilde = c / r0
    lam_W = lam ** r0
    M_cut = int(np.ceil(sqrt(45.0 / alpha))) + 10
    ms = np.arange(-(c + M_cut), M_cut + 1, dtype=np.float64)
    re_a = -alpha * ms * (c + ms) + log(lam)
    im_a = -2.0 * pi * ms * j / c
    keep = re_a > -35.0
    re_a = re_a[keep]; im_a = im_a[keep]
    # Compare Re log F(lam e^{i (pi/r0 + dth)}) - Re log F(lam e^{i pi/r0}) vs c_tilde*(Psi(r0*dth) - Psi(0))!
    diffs = []
    def re_logF_disc(psi):
        big = re_a > 0.0
        tot = -c * log(lam)
        if np.any(big):
            tot += float(np.sum(re_a[big])) + float(np.sum(np.log(np.abs(1.0 - np.exp(-re_a[big] - 1j*(im_a[big] + psi))))))
        sm = ~big
        if np.any(sm):
            tot += float(np.sum(np.log(np.abs(1.0 - np.exp(re_a[sm] + 1j*(im_a[sm] + psi))))))
        return tot
    p0_d = re_logF_disc(pi / r0)
    p0_c = c_tilde * Re_Psi_cont(abs(D), tau, lam_W, 0.0)
    for th_W in [0.2, 0.5, 1.0, 2.0, pi]:
        d_disc = re_logF_disc(pi / r0 + th_W / r0) - p0_d
        d_cont = c_tilde * Re_Psi_cont(abs(D), tau, lam_W, th_W) - p0_c
        diffs.append((th_W, d_disc, d_cont, abs(d_disc - d_cont)))
    s_d = " | ".join(f"thW={t:.1f}: disc={dd:+.3f} cont={dc:+.3f} err={e:.4f}" for t, dd, dc, e in diffs[:4])
    print(f"  c={c} j={j} r0={r0:2d} D={D:+2d} tau={tau:.2f} lamW={lam_W:.3f}:")
    print(f"     {s_d}")
    print(f"     at thW=pi: disc={diffs[4][1]:+.3f} cont={diffs[4][2]:+.3f} err={diffs[4][3]:.4f}")
