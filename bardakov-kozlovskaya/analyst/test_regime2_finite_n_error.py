"""
TEST THE COMPLETE REGIME 2 (kappa < 3.5) ARCHITECTURE:
1. Check how close the exact finite-n integral I_exact(n, c, j) is to the pure Gaussian-phase
   model integral I_model(alpha, c, j) = (-1)^j * T(alpha; j, c) * (normalization) across ALL (c, j)
   in Regime 2 (kappa < 3.5) at n = 80, 160, 320, 1000!
   (Wait: check the exact normalization: prod_{m=-c+1}^{-1} |a_m^{model}| vs I_exact!)
2. Scan EVERY SINGLE (c, j) in Regime 2 (kappa < 3.5) at n = 80, 160, 320, 1000 with
   lam = min(1.0, lam_*) and find the GLOBAL MINIMUM of rat = I / avg(|F|)!!
3. Check the exact formula for lam_* vs a CLOSED-FORM elementary lam_elem!
   Notice: earlier we found lam_* < 1 iff c * q * kappa < 1.15!
   Test closed-form:
     If c * q * kappa >= 1.0:  lam = 1.0
     If c * q * kappa <  1.0:  lam = exp(-c*kap/4) * (c * q * kappa / 1.15)^{1/q}  (or exact saddle of Phi_main!)
   Check if this closed-form lam gives rat >= 0.70 EVERYWHERE in Regime 2!!
"""
import numpy as np
from math import gcd

def kappa_val(n, c, j):
    return 2.0 * np.pi * np.sin(j * np.pi / n) * np.sin((j + c) * np.pi / n) / np.sin(c * np.pi / n)

def solve_lam_star(alpha, c, j):
    g = gcd(c, j); q = c // g
    M = int(np.ceil(c / (2*q) + np.sqrt(50.0 / max(q * alpha, 1e-8)) / q)) + 5
    ss = np.arange(-M, M + 1, dtype=np.float64)
    mid = q * ss + 0.5 * (q - 1)
    ex = -alpha * q * (mid * (c + mid) + (q*q - 1.0) / 12.0)
    lo, hi = -20.0, 0.0 # capped at 1.0 (hi = 0.0)!
    # Check value at hi = 0.0 (lam = 1.0):
    z0 = np.clip(ex, -60.0, 60.0)
    if float(np.sum(1.0 / (1.0 + np.exp(-z0)))) <= g:
        return 1.0
    for _ in range(40):
        vm = 0.5 * (lo + hi)
        z = np.clip(ex + q * vm, -60.0, 60.0)
        if float(np.sum(1.0 / (1.0 + np.exp(-z)))) > g:
            hi = vm
        else:
            lo = vm
    return float(np.exp(0.5 * (lo + hi)))

def eval_exact_and_model(n, c, j, lam, M_psi=128):
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    eip = np.exp(1j * psis)
    mid_fac = 2.0 - lam * eip - (1.0 / lam) * np.conj(eip)
    # 1. Exact finite-n integrand:
    if c % 2 == 1:
        ths = ((2*np.arange(n) + 1) * np.pi / n + np.pi) % (2*np.pi) - np.pi
    else:
        ths = (2*np.arange(n) * np.pi / n + np.pi) % (2*np.pi) - np.pi
    th_star = np.pi * (n - c) / n
    ext_mask = np.abs(ths) > th_star + 1e-12
    int_mask = np.abs(ths) < th_star - 1e-12
    rj = np.sin(j * np.pi / n) / np.sin((j + c) * np.pi / n)
    sj = np.sin(c * np.pi / n) / np.sin((j + c) * np.pi / n)
    sgn_j = 1.0 if (j % 2 == 0) else -1.0
    log_K = np.zeros(M_psi, dtype=np.complex128)
    if np.any(ext_mask):
        w_ext = (1.0 + rj * np.exp(1j * ths[ext_mask])) / sj
        nz = np.abs(w_ext) > 0
        if np.any(nz):
            lw = n * np.log(w_ext[nz]) - np.log(lam) + (1j*np.pi if sgn_j < 0 else 0.0)
            log_K += np.sum(np.log(1.0 - np.exp(lw[None, :] - 1j*psis[:, None])), axis=1)
    w_int = sj / (1.0 + rj * np.exp(1j * ths[int_mask]))
    lw = n * np.log(w_int) + np.log(lam) + (1j*np.pi if sgn_j < 0 else 0.0)
    keep = lw.real > -40.0
    if np.any(keep):
        log_K += np.sum(np.log(1.0 - np.exp(lw[keep][None, :] + 1j*psis[:, None])), axis=1)

    # 2. Pure Gaussian-phase model integrand (same normalization!):
    # Exact alpha = 2*pi^2 * j*(j+c) / (c^2 * n)  (or using exact kap / c!)
    kap = kappa_val(n, c, j)
    alpha = kap / c
    # Exterior m = -(c-1)..-1: w_ext^{model} = exp(+2pi i m j/c + alpha*m*(c+m))
    # Wait: check sign! Earlier we had a_m = (-1)^j s_j^n / R_m = exp(-2pi i m j/c - alpha*m*(c+m)),
    # so R_m / ((-1)^j s_j^n) = 1/a_m = exp(+2pi i m j/c + alpha*m*(c+m))!
    # And for interior m >= 1 or m <= -(c+1): ((-1)^j s_j^n) / R_m = a_m = exp(-2pi i m j/c - alpha*m*(c+m))!
    log_Km = np.zeros(M_psi, dtype=np.complex128)
    if c > 1:
        ms_ext = np.arange(-(c - 1), 0, dtype=np.float64)
        lw_m = alpha * ms_ext * (c + ms_ext) + 2j * np.pi * ms_ext * j / c - np.log(lam)
        log_Km += np.sum(np.log(1.0 - np.exp(lw_m[None, :] - 1j*psis[:, None])), axis=1)
    M_cut = int(np.ceil(np.sqrt(40.0 / max(alpha, 1e-8)))) + 2
    ms_int = np.concatenate([np.arange(1, M_cut + 1, dtype=np.float64),
                             np.arange(-(c + M_cut), -c, dtype=np.float64)])
    lw_mi = -alpha * ms_int * (c + ms_int) - 2j * np.pi * ms_int * j / c + np.log(lam)
    keep_m = lw_mi.real > -40.0
    if np.any(keep_m):
        log_Km += np.sum(np.log(1.0 - np.exp(lw_mi[keep_m][None, :] + 1j*psis[:, None])), axis=1)

    sh = max(float(np.max(log_K.real)), float(np.max(log_Km.real)))
    integ_ex = mid_fac * np.exp(log_K - sh)
    integ_mo = mid_fac * np.exp(log_Km - sh)
    val_ex = float(np.mean(integ_ex).real)
    val_mo = float(np.mean(integ_mo).real)
    rat_ex = val_ex / float(np.mean(np.abs(integ_ex)))
    rat_mo = val_mo / float(np.mean(np.abs(integ_mo)))
    rel_diff_I = abs(val_ex - val_mo) / abs(val_mo)
    rel_diff_L1 = float(np.mean(np.abs(integ_ex - integ_mo))) / abs(val_mo)
    return rat_ex, rat_mo, rel_diff_I, rel_diff_L1

print("=== FULL SCAN OF ALL REGIME-2 POINTS (kappa < 3.5) AT n = 80, 160, 320, 1000 ===")
for n in [80, 160, 320, 1000]:
    cnt = 0
    min_rat = 1e9; worst_rat = None
    max_Ierr = 0.0; worst_Ierr = None
    max_L1err = 0.0; worst_L1err = None
    for c in range(2, int(0.12 * n) + 2):
        for j in range(c + 1, int(0.25 * n) + 2):
            if j > (n - c) // 2: continue
            kap = kappa_val(n, c, j)
            if kap >= 3.5: continue
            cnt += 1
            alpha = kap / c
            lam = solve_lam_star(alpha, c, j)
            rat_ex, rat_mo, di, dl1 = eval_exact_and_model(n, c, j, lam, M_psi=64)
            if rat_ex < min_rat: min_rat = rat_ex; worst_rat = (c, j, kap, lam, rat_ex, rat_mo)
            if di > max_Ierr: max_Ierr = di; worst_Ierr = (c, j, kap, lam, di, dl1)
            if dl1 > max_L1err: max_L1err = dl1; worst_L1err = (c, j, kap, lam, di, dl1)
    print(f"n={n:4d} ({cnt:4d} R2 pts): min_rat={min_rat:.4f} at {worst_rat[:4] if worst_rat else None} | max |I_ex/I_mo - 1|={max_Ierr:.4f} at {worst_Ierr[:4] if worst_Ierr else None} | max L1_err/|I_mo|={max_L1err:.4f} at {worst_L1err[:4] if worst_L1err else None}")
