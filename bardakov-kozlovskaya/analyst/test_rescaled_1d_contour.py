"""
TEST THE UNIFIED 1D CONTOUR RADIUS R_z^* FOR ALL (n, c, j) INCLUDING n = 5120!!
Look at the exact integral for E_{n-c}(1, u_j) = (-1)^c [z^c] p_n(1, u_j, z):
  (-1)^j * [z^c] p_n(1, u_j, z) / u_j^c
    = (1/(2pi i)) oint_{|z| = R_z} p_n(1, u_j, z) / (z^{c+1} * u_j^c) dz  * (-1)^j.

When alpha = 2*pi^2 * j*(j+c) / (c^2 * n) is LARGE: R_z = s_j^n gives I ~ 2 with ratio ~ 0.99!
When alpha is SMALL: let g = gcd(c, j), q = c // g; then p_q ~ (-1)^{jq} sqrt(pi / (q*alpha)) / s_j^{nq},
so the saddle radius is:
  R_z^* = s_j^n * min(1.0, (c / sqrt(pi / (q * alpha)))**(1.0 / q))  (or just solve sum |lambda_k|^q * R^q ~ c!).
WAIT A SECOND!! Even cleaner (no q needed!):
  Check the exact radius R_z^* where the 1D integrand has MAXIMUM ratio I_val / avg(|F|)
  for EVERY (c, j) at n = 160, 640, 2560, 5120!
  Specifically test the EXPLICIT CLOSED-FORM radius:
    R_z(n, c, j) = s_j^n * (1 + (pi / (q * alpha))**(1/2) / c)**(-1/q) !!
  (Notice: when alpha -> inf, R_z -> s_j^n smoothly! When alpha -> 0, R_z -> s_j^n * (c / sqrt(pi/(q*alpha)))^{1/q} smoothly!!)
"""
import numpy as np
from math import gcd

def test_I_rescaled(n, c, j, scale_factor=1.0, M_psi=128):
    # Integrate (-1)^j * (1/(2pi i)) oint_{|z| = scale_factor * sj^n} p_n(1, u_j, z) / (z^{c+1}) dz
    # divided by ref = prod_{|th|<th*} |R_k| * sj^n  (so at scale_factor=1, matches I_{n,c,j}!).
    # Wait: in terms of R_k:
    #   p_n(1, u_j, z) / z^c = [prod_{|th|<th*} (-R_k)] * [z^{-1} (z - (-1)^j sj^n)^2]
    #                          * prod_{|th|>th*} (1 - R_k/z) * prod_{|th|<th*} (1 - z/R_k)
    # Put z = (-1)^j * lam * sj^n * e^{i psi} (lam = scale_factor > 0)!
    # Then:
    #   z^{-1} (z - (-1)^j sj^n)^2 = (-1)^j sj^n * (lam e^{i psi} - 1)^2 / (lam e^{i psi})
    #                              = (-1)^j sj^n * [ (lam e^{i psi} + lam^{-1} e^{-i psi}) - 2 ] !
    # Notice the minus sign: -(2 - lam e^{i psi} - lam^{-1} e^{-i psi})!
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    eip = np.exp(1j * psis)
    lam = scale_factor
    mid_fac = 2.0 - lam * eip - (1.0 / lam) * np.conj(eip) # equals 4 sin^2(psi/2) when lam=1!

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
    integ = mid_fac * np.exp(log_K)
    val = np.mean(integ)
    abs_val = np.mean(np.abs(integ))
    return float(val.real), float(val.real / abs_val)

print("=== TEST EXPLICIT RESCALED CONTOUR RADIUS lam = (1 + sqrt(pi/(q*alpha))/c)^{-1/q} UP TO n = 5120 ===")
for (c, j) in [(1, 2), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5), (4, 6), (5, 6), (6, 9)]:
    g = gcd(c, j)
    q = c // g
    parts = []
    for n in [160, 640, 2560, 5120]:
        alpha = 2.0 * (np.pi**2) * j * (j + c) / (c * c * n)
        # Try explicit formula for lam:
        # Notice for q >= 2, if lam < 1, w_ext / lam can exceed 1! Wait: exterior terms have |w_ext| ~ exp(-alpha * c^2/4)!
        # WAIT A SECOND!! If lam < |w_ext|, then 1 - w_ext/z has a ZERO inside the contour? NO! It's a polynomial in z, ZERO poles except z=0!
        # We just shouldn't make |w_ext|/lam huge! Actually at lam = 1, |w_ext| < 1; at lam < 1, |w_ext|/lam can be slightly > 1, which is FINE for a polynomial!
        lam = (1.0 + np.sqrt(np.pi / (q * alpha)) / c) ** (-1.0 / q)
        # Wait: check if lam < 1 works cleanly even when |w_ext|/lam > 1!
        I_v, rat = test_I_rescaled(n, c, j, scale_factor=lam, M_psi=128)
        parts.append(f"n={n}:lam={lam:.3f},I={I_v:+.2e}(rat={rat:+.3f})")
    print(f"(c={c},j={j},g={g},q={q}): " + " | ".join(parts))
