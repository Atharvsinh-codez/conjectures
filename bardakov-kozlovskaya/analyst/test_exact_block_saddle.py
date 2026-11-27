"""
DESIGN & STRESS-TEST THE EXACT REGIME 2 (kappa < 4.0) PROOF ENGINE!
Test across ALL (c, j) with kappa < 4.0 at n = 80, 320, 1000, 5000:
1. Check how closely the exact finite-n roots a_m = (-1)^j s_j^n / (1 + r_j e^{i theta_m})^n
   match the Gaussian-phase model a_m^{model} = e^{-2pi i m j / c} e^{-alpha m(c + m)}!
2. Compute the EXACT saddle radius lam_* solving:
     sum_{s in Z} (B_s * lam^q) / (1 + B_s * lam^q) = g
   where B_s = exp(-q * alpha * (q*s + (q-1)/2) * (c + q*s + (q-1)/2))!
   Test the 1D contour integral on |w| = lam_* across EVERY tricky (c, j) in Regime 2 up to n = 5000!
3. ALSO test a SUPER-SIMPLE 2-line explicit formula for lam_* (no root-finding needed!):
   Notice:
   - When q * kappa >= 0.5: earlier we saw lam = 1.0 ALREADY gives rat >= 0.94!
   - When q * kappa < 0.5: how does lam_* behave?! Check lam_* vs 1.0!
4. ALSO check: For q = 1 (c | j), what is max |arg(a_m)| across all significant roots?
   Wait: if |arg(a_m)| < pi/(2c) for all significant roots when q=1? No, c can be 50!
   Wait: when we integrate on |w| = lam_*, the product prod_m (1 + a_m w) concentrates near psi = 0
   (width O(1/sqrt(g))) where EVERY factor (1 + a_m lam_*) is positive real!
"""
import numpy as np
from math import gcd

def kappa_val(n, c, j):
    return 2.0 * np.pi * np.sin(j * np.pi / n) * np.sin((j + c) * np.pi / n) / np.sin(c * np.pi / n)

def solve_lam_star(alpha, c, j):
    g = gcd(c, j); q = c // g
    # B_s for block s: product of q consecutive roots m = q*s + r (r=0..q-1)
    # sum_{r=0}^{q-1} m(c+m) = q * [ (q*s + (q-1)/2)*(c + q*s + (q-1)/2) + (q^2 - 1)/12 ]
    # Shift s so center is at s0 = -(c + q - 1)/(2q)
    M = int(np.ceil(c / (2*q) + np.sqrt(50.0 / max(q * alpha, 1e-8)) / q)) + 5
    ss = np.arange(-M, M + 1, dtype=np.float64)
    mid = q * ss + 0.5 * (q - 1)
    ex = -alpha * q * (mid * (c + mid) + (q*q - 1.0) / 12.0)
    # Solve sum_s 1 / (1 + exp(-ex - q*u)) = g for u = log(lam) via bisection!
    lo, hi = -20.0, 2.0
    for _ in range(50):
        vm = 0.5 * (lo + hi)
        # sigmoid(ex + q*vm)
        z = np.clip(ex + q * vm, -60.0, 60.0)
        val = float(np.sum(1.0 / (1.0 + np.exp(-z))))
        if val > g:
            hi = vm
        else:
            lo = vm
    return float(np.exp(0.5 * (lo + hi)))

def test_contour_lam(n, c, j, lam, M_psi=256):
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    # Evaluate (-1)^j * (1/(2pi i)) oint_{|w|=lam} prod_k (1 + (-1)^j s_j^n lam_k w) / w^{c+1} dw
    # Wait: in terms of p_n(1, u_j, z) at z = (-1)^j * s_j^n * e^{i psi} / lam?
    # Careful! w = (-1)^j / z * s_j^n -> 1 + (-1)^j s_j^n lam_k w = (R_k + (-1)^j s_j^n w) / R_k!
    # Wait: earlier we used z = (-1)^j * lam * s_j^n * e^{i psi}, which corresponds to |w| = 1/lam!
    # Let's check: [z^c] prod (z - R_k) = (-1)^{n-c} (prod R_k) [z^c] prod (1 - z / R_k)!
    # So w = -z / ((-1)^j s_j^n) has |w| = |z| / s_j^n = lam!
    # YES! Putting z = (-1)^j * lam * s_j^n * e^{i psi} IS literally |w| = lam!
    # Check: at lam -> 0, [z^c] prod(1 - z/R_k) = (-1)^c e_c(1/R_k), and z^{-c} has lam^{-c}!
    # So our previous test_I_rescaled(n, c, j, lam) ALREADY uses |w| = lam!
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
            log_K += np.sum(np.log(1.0 - np.exp(lw[None, :] - 1j*psis[:, None])), axis=1)
    w_int = sj / (1.0 + rj * np.exp(1j * ths[int_mask]))
    lw = n * np.log(w_int) + np.log(lam) + (1j*np.pi if sgn_j < 0 else 0.0)
    keep = lw.real > -40.0
    if np.any(keep):
        log_K += np.sum(np.log(1.0 - np.exp(lw[keep][None, :] + 1j*psis[:, None])), axis=1)
    sh = np.max(log_K.real)
    integ = mid_fac * np.exp(log_K - sh)
    val = np.mean(integ)
    abs_val = np.mean(np.abs(integ))
    return float(val.real / abs_val), float(np.min(integ.real) / val.real)

print("=== TEST EXACT BLOCK-SADDLE RADIUS lam_* ACROSS ALL HARD REGIME-2 POINTS UP TO n = 5000 ===")
for n in [1000, 5000]:
    print(f"\n--- n = {n} ---")
    cases = [
        (2, 3), (2, 4), (2, 5), (3, 4), (4, 5), (4, 6), (6, 9),
        (10, 11), (10, 15), (10, 20), (10, 30),
        (20, 21), (20, 25), (20, 30), (20, 40),
        (50, 51), (50, 60), (50, 75), (50, 100),
    ]
    for c, j in cases:
        kap = kappa_val(n, c, j)
        if kap > 4.5: continue
        g = gcd(c, j); q = c // g; alpha = kap / c
        lam_s = solve_lam_star(alpha, c, j)
        rat1, _ = test_contour_lam(n, c, j, 1.0, M_psi=256)
        rats, minrs = test_contour_lam(n, c, j, lam_s, M_psi=256)
        print(f"  c={c:2d}, j={j:3d} (g={g:2d},q={q:2d}) kap={kap:.3f}, q*kap={q*kap:5.2f}: unscaled={rat1:+.4f} | lam_*={lam_s:.3f} -> rat={rats:+.4f} (minRe/avg={minrs:+.3f})")
