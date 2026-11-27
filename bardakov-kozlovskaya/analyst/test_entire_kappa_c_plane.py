"""
CRITICAL CHECK OF BOTH REGIMES BEFORE LOCKING IN THE FINAL PROOF ARCHITECTURE:
1. Regime 1 (kappa >= kappa0) on STRICTLY INTERIOR points (dropping j=c AND dropping j=m-c/2 for even c!):
   Check min I_full, max |I_full - I_E|, max_sum |w|^n at n = 40, 80, 120, 200!
2. Regime 2 (kappa < kappa0) across ALL (c, j) INCLUDING LARGE c = 10, 20, 50 with small kappa
   at n = 320, 1000, 2000, 5000!
   Compare unscaled lam = 1 vs rescaled lam = (1 + sqrt(pi/(q*alpha))/c)^{-1/q} vs optimal lam!
"""
import numpy as np
from math import gcd

def kappa_val(n, c, j):
    return 2.0 * np.pi * np.sin(j * np.pi / n) * np.sin((j + c) * np.pi / n) / np.sin(c * np.pi / n)

def test_I_rescaled(n, c, j, lam=1.0, M_psi=128):
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    eip = np.exp(1j * psis)
    mid_fac = 2.0 - lam * eip - (1.0 / lam) * np.conj(eip)
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
            fac = 1.0 - np.exp(lw[None, :] - 1j*psis[:, None])
            log_K += np.sum(np.log(fac), axis=1)
    w_int = sj / (1.0 + rj * np.exp(1j * ths[int_mask]))
    lw = n * np.log(w_int) + np.log(lam) + (1j*np.pi if sgn_j < 0 else 0.0)
    keep = lw.real > -40.0
    if np.any(keep):
        fac = 1.0 - np.exp(lw[keep][None, :] + 1j*psis[:, None])
        log_K += np.sum(np.log(fac), axis=1)
    # Shift log_K to avoid overflow!
    sh = np.max(log_K.real)
    integ = mid_fac * np.exp(log_K - sh)
    val = np.mean(integ)
    abs_val = np.mean(np.abs(integ))
    return float(val.real * np.exp(min(sh, 50.0))), float(val.real / abs_val), float(sh)

print("=== PART 1: Regime 1 (kappa >= 3.0) on STRICTLY INTERIOR points (drop j=c AND j=m-c/2 for even c) ===")
def eval_both_strict(n, c, j, M_psi=48):
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
    integ_E = sin2 * np.exp(log_K)
    I_E = float(np.mean(integ_E).real)
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
    return I_E, I_full, float(np.max(np.abs(integ - integ_E)))

for n in [40, 80, 120, 160]:
    m = n // 2
    min_If = 1e9; max_diff = 0.0; w_I = None; w_d = None
    for c in range(2, (n - 3) // 3 + 1):
        j_hi = m - (c + 1)//2 - (1 if c % 2 == 0 else 0)
        for j in range(c + 1, j_hi + 1):
            if kappa_val(n, c, j) < 3.0: continue
            I_E, I_f, md = eval_both_strict(n, c, j, M_psi=48)
            if I_f < min_If: min_If = I_f; w_I = (c, j, I_E, I_f, md)
            if md > max_diff: max_diff = md; w_d = (c, j, I_E, I_f, md)
    print(f"n={n:3d} strict R1(kap>=3): min I_full={min_If:.4f} at {w_I[:2]} (IE={w_I[2]:.4f}, md={w_I[4]:.4f}) | max|F-FE|={max_diff:.4f} at {w_d[:2]} (IE={w_d[2]:.4f}, If={w_d[3]:.4f})")

print("\n=== PART 2: Regime 2 (kappa < 4.0) stress test INCLUDING LARGE c = 10, 20, 50 at n = 1000, 5000 ===")
for n in [1000, 5000]:
    print(f"\n--- n = {n} ---")
    cases = [
        (10, 11), (10, 15), (10, 20), (10, 30),
        (20, 21), (20, 25), (20, 30), (20, 40),
        (50, 51), (50, 60), (50, 75), (50, 100),
    ]
    for c, j in cases:
        if j > (n - c)//2: continue
        kap = kappa_val(n, c, j)
        if kap > 5.0: continue
        g = gcd(c, j); q = c // g; alpha = kap / c
        lam = (1.0 + np.sqrt(np.pi / (q * alpha)) / c) ** (-1.0 / q)
        _, rat1, sh1 = test_I_rescaled(n, c, j, lam=1.0, M_psi=256)
        _, rat_l, sh_l = test_I_rescaled(n, c, j, lam=lam, M_psi=256)
        print(f"  c={c:2d}, j={j:3d} (g={g:2d},q={q:2d}) kap={kap:.3f}, c*kap={c*kap:6.2f}: unscaled rat={rat1:+.4f}(sh={sh1:5.1f}) | rescaled lam={lam:.3f} rat={rat_l:+.4f}(sh={sh_l:5.1f})")
