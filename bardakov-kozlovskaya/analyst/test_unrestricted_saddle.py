"""
TEST UNRESTRICTED SADDLE RADIUS lam_* in (0, infinity) (NO CAP AT 1.0!!)
on ALL extreme scaling corners of (c, q, kappa) — including:
  c=500, j=1000, q=1, kap=0.01 (where lam=1 had rat=0.0012!)
  c=200, j=201, q=200, kap=0.01 (where lam=1 had rat=0.2157!)
Check BOTH:
1. Pure block saddle lam_block (solving sum_s B_s lam^q / (1 + B_s lam^q) = g, unrestricted!)
2. Exact saddle of the real positive product:
   At psi = pi/q, wait: check rat = I / avg|F| on |w| = lam_block!!
"""
import numpy as np
from math import gcd

def solve_lam_unrestricted(alpha, c, j):
    g = gcd(c, j); q = c // g
    M = int(np.ceil(c / (2*q) + np.sqrt(60.0 / max(q * alpha, 1e-12)) / q)) + 10
    ss = np.arange(-M, M + 1, dtype=np.float64)
    mid = q * ss + 0.5 * (q - 1)
    ex = -alpha * q * (mid * (c + mid) + (q*q - 1.0) / 12.0)
    # Solve sum_s 1 / (1 + exp(-(ex + q*v))) = g for v = log(lam) in (-30, +30)
    lo, hi = -30.0, 30.0
    for _ in range(60):
        vm = 0.5 * (lo + hi)
        z = np.clip(ex + q * vm, -60.0, 60.0)
        val = float(np.sum(1.0 / (1.0 + np.exp(-z))))
        if val > g:
            hi = vm
        else:
            lo = vm
    return float(np.exp(0.5 * (lo + hi)))

def eval_model_at_lam(alpha, c, j, lam, M_psi=512):
    # Evaluate exact Gaussian-phase model integrand on |w| = lam!
    # Formula: F(psi) = (-1)^j * (lam*e^{i psi})^{-c} * prod_{m in Z} (1 - a_m * lam * e^{i psi})
    # Normalized by prod_{m=-(c-1)}^{-1} |a_m|:
    # Equivalently, keep all m in [-c-M_cut, M_cut] directly without splitting at |a_m|=1!
    # Wait: for stability, factor out a_m*lam*e^{i psi} whenever |a_m|*lam > 1!
    M_cut = int(np.ceil(np.sqrt(60.0 / max(alpha, 1e-12)))) + 10
    ms = np.arange(-(c + M_cut), M_cut + 1, dtype=np.float64)
    # log(a_m * lam) = -alpha*m*(c+m) + log(lam) - 2pi i m j/c
    re_a = -alpha * ms * (c + ms) + np.log(lam)
    im_a = -2.0 * np.pi * ms * j / c
    # Keep terms with re_a > -45
    keep = re_a > -45.0
    ms = ms[keep]; re_a = re_a[keep]; im_a = im_a[keep]
    big = re_a > 0.0
    # Wait: for |z| > 1, 1 - z = (-z) * (1 - 1/z)
    # Sum of log(1 - a_m * lam * e^{i psi}) - c * (log(lam) + i*psi):
    # For big terms: log(-a_m * lam * e^{i psi}) + log(1 - a_m^{-1} lam^{-1} e^{-i psi})
    # Notice overall constant phase from (-1)^j * prod_{big} (-a_m) is purely real sign!
    # Let's compute exact log_F(psi) directly in complex arithmetic (only ~1000 terms):
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    log_F = -c * (np.log(lam) + 1j * psis)
    if j % 2 == 1:
        log_F = log_F + 1j * np.pi

    # Big terms: re_a > 0
    if np.any(big):
        rb = re_a[big]; ib = im_a[big]
        # log(-a_m * lam) = rb + i*(ib + pi)
        # sum over big:
        log_F = log_F + np.sum(rb) + 1j * (np.sum(ib) + len(rb) * np.pi) + len(rb) * 1j * psis
        inv_u = np.exp(-rb - 1j * ib) # |inv_u| < 1
        ep_neg = np.exp(-1j * psis)
        for k0 in range(0, len(inv_u), 200):
            ch = inv_u[k0:k0+200]
            log_F += np.sum(np.log(1.0 - ch[None, :] * ep_neg[:, None]), axis=1)

    # Small/equal terms: re_a <= 0
    sm = ~big
    if np.any(sm):
        rs = re_a[sm]; is_ = im_a[sm]
        u = np.exp(rs + 1j * is_) # |u| <= 1
        ep_pos = np.exp(1j * psis)
        for k0 in range(0, len(u), 200):
            ch = u[k0:k0+200]
            # handle exact |u|=1 safely (add tiny jitter? No, M_psi doesn't hit exact 0 unless psi=0,
            # use direct product 1 - u*e^{i psi} via complex log or direct mult for |u| > 0.999999)
            near1 = np.abs(ch) > 1.0 - 1e-9
            if np.any(near1):
                # exact factor
                prod_n1 = np.prod(1.0 - ch[near1][None, :] * ep_pos[:, None], axis=1)
                # avoid log(0)
                with np.errstate(divide='ignore'):
                    log_F += np.where(np.abs(prod_n1) > 0, np.log(prod_n1), -1e9)
            far = ~near1
            if np.any(far):
                log_F += np.sum(np.log(1.0 - ch[far][None, :] * ep_pos[:, None]), axis=1)

    sh = float(np.max(log_F.real))
    F = np.exp(log_F - sh)
    val = float(np.mean(F).real)
    im_err = abs(float(np.mean(F).imag))
    rat = val / float(np.mean(np.abs(F)))
    return rat, im_err, int(np.sum(big))

print("=== UNRESTRICTED SADDLE RADIUS lam_* (NO CAP AT 1.0!!) ON ALL EXTREME CORNERS ===")
print(f"{'c':>5} {'j':>6} {'q':>4} {'kap':>7} {'c*q*k':>8} {'lam_*':>10} {'#big':>5} {'rat@lam*':>10} {'im_err':>9}")
cases = [
    (2, 3, 0.5), (2, 3, 0.2), (2, 3, 0.05),
    (2, 4, 1.5), (2, 4, 0.5), (2, 4, 0.1), (2, 4, 0.01),
    (3, 4, 0.5), (3, 4, 0.1), (3, 6, 1.0), (3, 6, 0.3),
    (4, 8, 1.5), (4, 8, 0.5), (4, 8, 0.1),
    (20, 40, 1.0), (20, 40, 0.3), (20, 40, 0.1), (20, 40, 0.05),
    (100, 200, 1.0), (100, 200, 0.3), (100, 200, 0.1), (100, 200, 0.03), (100, 200, 0.01),
    (500, 1000, 0.1), (500, 1000, 0.03), (500, 1000, 0.01), (500, 1000, 0.003),
    (100, 150, 0.3), (100, 150, 0.1), (100, 150, 0.03), (100, 150, 0.01),
    (50, 51, 1.0), (50, 51, 0.3), (50, 51, 0.1), (50, 51, 0.03), (50, 51, 0.01),
    (200, 201, 0.3), (200, 201, 0.05), (200, 201, 0.01),
]
for c, j, kap in cases:
    g = gcd(c, j); q = c // g
    alpha = kap / c
    lam = solve_lam_unrestricted(alpha, c, j)
    rat, ime, nbig = eval_model_at_lam(alpha, c, j, lam, M_psi=2048)
    print(f"{c:5d} {j:6d} {q:4d} {kap:7.4f} {c*q*kap:8.3f} {lam:10.6f} {nbig:5d} {rat:10.5f} {ime:9.1e}")
