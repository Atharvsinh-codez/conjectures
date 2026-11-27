"""
TEST THE FULL 1D CONTOUR INTEGRAL FOR f_{n, c}(u_j) DIRECTLY (INCLUDING x(z), y(z), e3(z))!!
For ANY even n = 2m, ANY 1 <= c <= floor(n/3) with d = n - 3c >= 2, and ANY j in {c+1, ..., (n-c)//2}:
  u_j = (-1)^c * r_j^n,  r_j = sin(j*pi/n)/sin((j+c)*pi/n) <= 1,  s_j = sin(c*pi/n)/sin((j+c)*pi/n) < r_j!
By Cauchy's integral formula in z around |z| = s_j^n:
  f_{n, c}(u_j) = (1/(2pi i)) oint_{|z| = s_j^n} p_n(x(z), y(z), z) * e3'(z) / e3(z)^{c+1} dz
where:
  e3(z) = z * (1 - z) * (u_j - z),
  x(z) + y(z) = 1 + u_j - z,  x(z) * y(z) = (1 - z) * (u_j - z)!
WAIT A SECOND!!
Check whether the pole z = u_j lies OUTSIDE |z| = s_j^n!
Since j >= c + 1:
  |z| / |u_j| = (s_j / r_j)^n = (sin(c*pi/n) / sin(j*pi/n))^n!
  Check: sin(j*pi/n) > sin(c*pi/n) holds WHENEVER j + c < n (ALWAYS true since j <= (n-c)/2 < n - c!)
  AND j > c (true since j >= c + 1!!)!
  SO |z| < |u_j| <= 1 STRICTLY AT EVERY SINGLE INTERIOR SAMPLE POINT j = c+1 .. (n-c)//2!!
  NO POLES OF e3(z) OTHER THAN z = 0 LIE INSIDE |z| = s_j^n!!
  Wait: WHAT ABOUT THE BRANCH POINTS OF x(z), y(z)?!
  x(z), y(z) satisfy W^2 - (1 + u_j - z) W + (1 - z)(u_j - z) = 0,
  whose discriminant is:
    Delta(z) = (1 + u_j - z)^2 - 4(1 - z)(u_j - z) = (1 - u_j)^2 + 2(1 + u_j) z - 3 z^2!!
  WAIT: p_n(x, y, z) IS SYMMETRIC IN (x, y)!! SO IT HAS **ZERO BRANCH CUTS** IN z — IT IS A PURE POLYNOMIAL IN z!!
  THEREFORE THE CAUCHY INTEGRAL FORMULA AROUND |z| = s_j^n IS **100% EXACT FOR EVERY SINGLE j IN {c+1, ..., (n-c)//2}** (INCLUDING j = (n-c)/2 WHERE u_j = +1 AND x(z), y(z) MEET AT z = 0!!)!!

TEST THIS EXACT 1D INTEGRAL AT EVERY SINGLE (c, j) FOR n = 12, 24, 40, 80, 160!!
Compare the FULL integrand for f_{n, c}(u_j) with the pure E_{n-c}(1, u_j)/u_j^c integrand!
Check BOTH:
  (A) exact sign match,
  (B) the normalized ratio: f_{n,c}(u_j) / ( (-1)^{c+j} * |const| )  [i.e. I_full vs I_E!],
  (C) the minimum over ALL (c, j) of I_full (in units of the main positive term)!!
"""
import numpy as np

def test_1d_full_f(n, M_psi=64):
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    eip = np.exp(1j * psis)
    zetas = np.exp(2j * np.pi * np.arange(n) / n)

    min_Ifull = 1e9
    max_Ifull = -1e9
    worst_info = None
    all_pos = True
    n_pts = 0
    for c in range(1, n // 3 + 1):
        d = n - 3 * c
        if d < 2: continue # d=0 is constant (b=0 proved), d=1 has only 1 root at u=-1 (odd d, palindromic => u=-1 is exact root, sign at 0 proved!)
        j_start = c + 1
        j_end = (n - c) // 2
        if c % 2 == 1:
            ths = ((2*np.arange(n) + 1) * np.pi / n + np.pi) % (2*np.pi) - np.pi
        else:
            ths = (2*np.arange(n) * np.pi / n + np.pi) % (2*np.pi) - np.pi
        th_star = np.pi * (n - c) / n
        int_mask = np.abs(ths) < th_star - 1e-12
        for j in range(j_start, j_end + 1):
            n_pts += 1
            rj = np.sin(j * np.pi / n) / np.sin((j + c) * np.pi / n)
            sj = np.sin(c * np.pi / n) / np.sin((j + c) * np.pi / n)
            uj = ((-1)**c) * (rj**n)
            sgn_j = 1.0 if (j % 2 == 0) else -1.0
            # z on circle |z| = sj^n centered at 0, phased by (-1)^j
            z_v = sgn_j * (sj**n) * eip # (M,)
            # Solve W^2 - (1 + uj - z) W + (1 - z)(uj - z) = 0 cleanly via Vieta (NO float64 cancellation!)
            # x_v = 1 + O(z), y_v = (1 - z)(uj - z) / x_v!
            # Better: disc = sqrt((1 - uj)^2 + 2*(1 + uj)*z - 3*z^2)
            # Notice at uj = +1 (j = (n-c)/2, even c), disc = sqrt(4*z - 3*z^2), and x_v, y_v are both O(1).
            # When |uj| < 0.5, x_v = 0.5*(1 + uj - z + disc) ~ 1, y_v = c_q / x_v has ZERO cancellation!
            b_q = 1.0 + uj - z_v
            c_q = (1.0 - z_v) * (uj - z_v)
            disc_in = (1.0 - uj)**2 + 2.0*(1.0 + uj)*z_v - 3.0*z_v*z_v
            disc = np.sqrt(disc_in)
            x_v = 0.5 * (b_q + disc)
            # Ensure we pick the larger root for x_v so c_q / x_v is accurate
            swap = np.abs(x_v) < np.abs(0.5 * (b_q - disc))
            x_v = np.where(swap, 0.5 * (b_q - disc), x_v)
            y_v = c_q / x_v
            # Compute log(p_n(x_v, y_v, z_v)) accurately using log-sum!
            # Wait: p_n(x, y, z) = prod_{k=0}^{n-1} (z - (x^{1/n} + zeta^k y^{1/n})^n)
            # Use principal n-th roots of x_v, y_v (symmetric, branch doesn't matter!)
            ax = np.exp(np.log(x_v) / n)
            by = np.exp(np.log(y_v) / n)
            w = ax[:, None] + by[:, None] * zetas[None, :] # (M, n)
            log_w = n * np.log(w)
            log_z = np.log(z_v)[:, None]
            use_z = log_z.real > log_w.real
            log_fac = np.where(use_z,
                               log_z + np.log(1.0 - np.exp(np.where(use_z, log_w - log_z, -100.0))),
                               log_w + 1j*np.pi + np.log(1.0 - np.exp(np.where(~use_z, log_z - log_w, -100.0))))
            log_pn = np.sum(log_fac, axis=1)
            # Cauchy measure: (1/(2pi i)) * (e3'(z) / e3(z)^{c+1}) dz
            # With z = z_v(psi), dz/(2pi i) = z_v d(psi)/(2pi)!
            # e3(z) = z*(1-z)*(uj-z), e3'(z) = (1-z)*(uj-z) - z*(uj-z) - z*(1-z)
            e3 = z_v * (1.0 - z_v) * (uj - z_v)
            e3p = (1.0 - z_v)*(uj - z_v) - z_v*(uj - z_v) - z_v*(1.0 - z_v)
            log_meas = np.log(z_v) + np.log(e3p) - (c + 1) * np.log(e3)
            # Reference positive scale from E_{n-c}(1, u_j) / u_j^c:
            # log_ref = sum_{int_mask} n*log|1 + rj e^{i th_k}| + log(sj^n) - c*log(rj^n)
            # Wait: [z^c] p_n / u_j^c has prefactor:
            #   prod_{|th|<th*} |R_k| * s_j^n / |u_j|^c  (Wait: earlier z^{c-1} * (s_j^n * z) / z^c = s_j^n!)
            log_ref = np.sum(n * np.log(np.abs(1.0 + rj * np.exp(1j * ths[int_mask])))) + n * np.log(sj) - c * n * np.log(rj)
            # Expected sign of f_{n, c}(u_j):
            # At u=0, sign is (-1)^c; at j=c, same sign (-1)^c; at j=c+1, opposite (-1)^{c+1}; so at j, (-1)^j !!
            # WAIT: check sign at u_j for even/odd c!
            # Recall (-1)^c E_{n-c}(1, u_j) had sign (-1)^{c+j}, and dividing by u_j^c = (-1)^{c^2} r_j^{nc} = (-1)^c r_j^{nc}
            # gives (-1)^{c+j} / (-1)^c = (-1)^j !!
            exp_sgn = 1.0 if (j % 2 == 0) else -1.0
            log_integrand = log_pn + log_meas - log_ref
            if exp_sgn < 0: log_integrand += 1j * np.pi
            integrand = np.exp(log_integrand)
            I_full = float(np.mean(integrand).real)
            if I_full <= 0: all_pos = False
            if I_full < min_Ifull:
                min_Ifull = I_full
                worst_info = (c, j, I_full, float(np.mean(np.abs(integrand))), float((sj/rj)**n), float(rj**n))
            if I_full > max_Ifull:
                max_Ifull = I_full
    print(f"n={n:3d} ({n_pts:4d} sample pts): ALL_POS={all_pos} | I_full in [{min_Ifull:+.5f}, {max_Ifull:+.5f}] | worst (c,j,I,avg|F|,(s/r)^n,r^n)={worst_info}")

for n in [12, 24, 40, 80, 120]:
    test_1d_full_f(n, M_psi=64)
