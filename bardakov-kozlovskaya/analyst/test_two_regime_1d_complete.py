"""
COMPLETE UNIFORM DIAGNOSTIC FOR THE TWO-REGIME 1D PROOF (Regime 1: kappa >= kappa0; Regime 2: kappa < kappa0):
Goal: Find the cleanest explicit threshold kappa0 and check ALL worst-case constants across
      EVERY (c, j) with 2 <= c <= floor((n-2)/3) and c+1 <= j <= (n-c)//2 for n = 20, 40, 80, 120, 200!

Check BOTH:
1. In Regime 1 (kappa >= kappa0):
   - What is the EXACT maximum over ALL (c, j) of:
     (a) max_k |w_ext|^n and max_k |w_int|^n (adjacent-root ratio, should be <= exp(-c0 * kappa0)),
     (b) sum_k |w_k|^n (total tail sum! Notice if sum |w_k|^n < log(2), then |K(psi) - 1| is TINY!),
     (c) cubic-root perturbation error |f_{n,c}(u_j) / E_{n,c,j} - 1| (separate u_j <= 0.5 vs u_j > 0.5 near u=+1!),
     (d) I_full lower bound!
2. In Regime 2 (kappa < kappa0):
   - What is the MAXIMUM of |u_j| = r_j^n and (s_j/r_j)^n over ALL (c, j) with kappa < kappa0?!
   - What is the MINIMUM of ratio = I_val / avg|F| with our rescaled contour lam = (1 + sqrt(pi/(q*alpha))/c)^{-1/q}
     across ALL (c, j) with kappa < kappa0 (including the transition kappa ~ kappa0!) up to n = 1000?!
"""
import numpy as np
from math import gcd

def kappa_val(n, c, j):
    return 2.0 * np.pi * np.sin(j * np.pi / n) * np.sin((j + c) * np.pi / n) / np.sin(c * np.pi / n)

print("=== PART 1: Map exact tail sum S_tail = sum |w_ext|^n + sum |w_int|^n vs kappa across ALL (c,j) ===")
for n in [40, 80, 160, 320]:
    # Check for several kappa0 cutoffs: 1.5, 2.0, 2.5, 3.0, 4.0
    for kappa0 in [2.0, 3.0, 4.0]:
        max_rho = 0.0
        max_sum = 0.0
        worst_r1 = None
        # Regime 2 stats
        max_uj_r2 = 0.0
        max_sr_r2 = 0.0
        cnt_r1 = 0; cnt_r2 = 0
        for c in range(2, (n - 2) // 3 + 1):
            if c % 2 == 1:
                ths = ((2*np.arange(n) + 1) * np.pi / n + np.pi) % (2*np.pi) - np.pi
            else:
                ths = (2*np.arange(n) * np.pi / n + np.pi) % (2*np.pi) - np.pi
            th_star = np.pi * (n - c) / n
            ext_mask = np.abs(ths) > th_star + 1e-12
            int_mask = np.abs(ths) < th_star - 1e-12
            for j in range(c + 1, (n - c) // 2 + 1):
                kap = kappa_val(n, c, j)
                rj = np.sin(j * np.pi / n) / np.sin((j + c) * np.pi / n)
                sj = np.sin(c * np.pi / n) / np.sin((j + c) * np.pi / n)
                if kap >= kappa0:
                    cnt_r1 += 1
                    # Compute max_rho and sum |w|^n safely via log!
                    s_tot = 0.0; m_rho = 0.0
                    if np.any(ext_mask):
                        a_ext = np.abs(1.0 + rj * np.exp(1j * ths[ext_mask])) / sj
                        nz = a_ext > 0
                        if np.any(nz):
                            pows = np.exp(n * np.log(a_ext[nz]))
                            m_rho = max(m_rho, float(np.max(pows)))
                            s_tot += float(np.sum(pows))
                    a_int = sj / np.abs(1.0 + rj * np.exp(1j * ths[int_mask]))
                    pows = np.exp(n * np.log(a_int))
                    m_rho = max(m_rho, float(np.max(pows)))
                    s_tot += float(np.sum(pows))
                    if m_rho > max_rho: max_rho = m_rho
                    if s_tot > max_sum:
                        max_sum = s_tot
                        worst_r1 = (c, j, kap, m_rho, s_tot)
                else:
                    cnt_r2 += 1
                    uj_abs = float(np.exp(n * np.log(rj)))
                    sr_n = float(np.exp(n * np.log(sj / rj)))
                    if uj_abs > max_uj_r2: max_uj_r2 = uj_abs
                    if sr_n > max_sr_r2: max_sr_r2 = sr_n
        if n in [80, 320]:
            print(f"n={n:3d}, kap0={kappa0:.1f} | R1({cnt_r1:5d} pts): max|w|^n={max_rho:.4f}, max_sum|w|^n={max_sum:.4f} at {worst_r1[:3] if worst_r1 else None} | R2({cnt_r2:4d} pts): max|u_j|={max_uj_r2:.2e}, max(s/r)^n={max_sr_r2:.2e}")

print("\n=== PART 2: Check Regime 1 cubic-root shift |I_full - I_E| across ALL R1 points at n = 80, 160 ===")
# Specifically inspect EVERY point where |I_full - I_E| > 0.1 (where does cubic-root shift matter?)!
def eval_both_1d(n, c, j, M_psi=64):
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    eip = np.exp(1j * psis)
    sin2 = 4.0 * np.sin(psis / 2.0)**2
    zetas = np.exp(2j * np.pi * np.arange(n) / n)
    if c % 2 == 1:
        ths = ((2*np.arange(n) + 1) * np.pi / n + np.pi) % (2*np.pi) - np.pi
    else:
        ths = (2*np.arange(n) * np.pi / n + np.pi) % (2*np.pi) - np.pi
    th_star = np.pi * (n - c) / n
    ext_mask = np.abs(ths) > th_star + 1e-12
    int_mask = np.abs(ths) < th_star - 1e-12
    rj = np.sin(j * np.pi / n) / np.sin((j + c) * np.pi / n)
    sj = np.sin(c * np.pi / n) / np.sin((j + c) * np.pi / n)
    uj = ((-1)**c) * (rj**n)
    sgn_j = 1.0 if (j % 2 == 0) else -1.0
    # 1. Pure E integral I_E:
    log_K = np.zeros(M_psi, dtype=np.complex128)
    if np.any(ext_mask):
        w_ext = (1.0 + rj * np.exp(1j * ths[ext_mask])) / sj
        nz = np.abs(w_ext) > 0
        if np.any(nz):
            lw = n * np.log(w_ext[nz]) + (1j*np.pi if sgn_j < 0 else 0.0)
            log_K += np.sum(np.log(1.0 - np.exp(lw[None, :] - 1j*psis[:, None])), axis=1)
    w_int = sj / (1.0 + rj * np.exp(1j * ths[int_mask]))
    lw = n * np.log(w_int) + (1j*np.pi if sgn_j < 0 else 0.0)
    keep = lw.real > -40.0
    if np.any(keep):
        log_K += np.sum(np.log(1.0 - np.exp(lw[keep][None, :] + 1j*psis[:, None])), axis=1)
    I_E = float(np.mean(sin2 * np.exp(log_K)).real)
    # 2. Full f integral I_full (only compute if |uj| > 1e-12 or (sj/rj)^n > 1e-12, else I_full == I_E to 1e-10!):
    sr_n = float(np.exp(n * np.log(sj / rj)))
    if abs(uj) < 1e-10 and sr_n < 1e-10:
        return I_E, I_E, 0.0
    z_v = sgn_j * (sj**n) * eip
    b_q = 1.0 + uj - z_v
    c_q = (1.0 - z_v) * (uj - z_v)
    disc = np.sqrt((1.0 - uj)**2 + 2.0*(1.0 + uj)*z_v - 3.0*z_v*z_v)
    x_v = 0.5 * (b_q + disc)
    swap = np.abs(x_v) < np.abs(0.5 * (b_q - disc))
    x_v = np.where(swap, 0.5 * (b_q - disc), x_v)
    y_v = c_q / x_v
    ax = np.exp(np.log(x_v) / n)
    by = np.exp(np.log(y_v) / n)
    w = ax[:, None] + by[:, None] * zetas[None, :]
    log_w = n * np.log(w)
    log_z = np.log(z_v)[:, None]
    use_z = log_z.real > log_w.real
    log_fac = np.where(use_z,
                       log_z + np.log(1.0 - np.exp(np.where(use_z, log_w - log_z, -100.0))),
                       log_w + 1j*np.pi + np.log(1.0 - np.exp(np.where(~use_z, log_z - log_w, -100.0))))
    log_pn = np.sum(log_fac, axis=1)
    e3 = z_v * (1.0 - z_v) * (uj - z_v)
    e3p = (1.0 - z_v)*(uj - z_v) - z_v*(uj - z_v) - z_v*(1.0 - z_v)
    log_meas = np.log(z_v) + np.log(e3p) - (c + 1) * np.log(e3)
    log_ref = np.sum(n * np.log(np.abs(1.0 + rj * np.exp(1j * ths[int_mask])))) + n * np.log(sj) - c * n * np.log(rj)
    exp_sgn = 1.0 if (j % 2 == 0) else -1.0
    integ = np.exp(log_pn + log_meas - log_ref + (1j*np.pi if exp_sgn < 0 else 0.0))
    I_full = float(np.mean(integ).real)
    # Also compute max_psi |integ_full - integ_E|!
    integ_E = sin2 * np.exp(log_K)
    max_diff = float(np.max(np.abs(integ - integ_E)))
    return I_E, I_full, max_diff

for n in [80, 120]:
    worst_diff = 0.0; worst_d_info = None
    worst_I = 1e9; worst_I_info = None
    for c in range(2, (n - 2) // 3 + 1):
        for j in range(c + 1, (n - c) // 2 + 1):
            if kappa_val(n, c, j) < 3.0: continue
            I_E, I_full, mdiff = eval_both_1d(n, c, j, M_psi=48)
            if mdiff > worst_diff:
                worst_diff = mdiff
                worst_d_info = (c, j, n-3*c, I_E, I_full, mdiff)
            if I_full < worst_I:
                worst_I = I_full
                worst_I_info = (c, j, n-3*c, I_E, I_full, mdiff)
    print(f"n={n:3d} (kappa>=3.0): min I_full={worst_I:.4f} at (c,j,d,IE,Ifull,mdiff)={worst_I_info} | max_psi|F-FE|={worst_diff:.4f} at {worst_d_info}")
