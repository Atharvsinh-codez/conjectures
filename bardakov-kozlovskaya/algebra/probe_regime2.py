"""
Explore Regime 2 (kappa < 2.5, even n >= 132, strictly interior c >= 2, j >= c + 1):
1. Check exact size of R_y, R_z, R_z/R_y, c/n, j/n in Regime 2.
2. Inspect the exact frozen 1D contour integral:
     I(rho) = (1/(2*pi)) int_{-pi}^{pi} (-1)^j * [prod_{k=0}^{n-1} (1 - (rho e^{i phi}) / w_k)] / (rho e^{i phi})^c dphi
   at BOTH:
     (a) scaled radius rho = lam * R_z (testing optimal lam as a function of kappa, c, d=gcd(c,j)), AND
     (b) the exact real saddle point of F(z) = prod_{k=0}^{n-1} (1 - z / w_k) / z^c!
3. Check the EXACT Gaussian / Poisson-summation approximation to log prod(1 - z/w_k) on |z| = rho < (1-s)^n!
4. Also test how low Regime 1's kappa_0 can go if we tighten S_small / S_large slightly!
"""
import math
import cmath

def kappa_val(n, c, j):
    return 2.0 * math.pi * math.sin(j * math.pi / n) * math.sin((j + c) * math.pi / n) / math.sin(c * math.pi / n)

def eval_frozen_integrand(n, c, j, rho_over_R0, phi, use_R0=True):
    """Evaluate (-1)^j * (rho e^{i phi})^{-c} * prod_{k=0}^{n-1} (1 - (rho e^{i phi}) / w_k)
    where R0 = (1 - s)^n (or R_z = (r3/r1)^n).
    Return (val, log_abs_integrand, phase).
    """
    th_c = c * math.pi / n
    psi_j = j * math.pi / n
    alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1
    # Work in normalized units u_k = w_k / R0 = ((1 - s e^{i eta_k}) / (1 - s))^n!
    # Wait: Rz / R0 = ((r3/r1) / (1 - s))^n!
    R_scale_1 = (1.0 - s) if use_R0 else (r3 / r1)
    # z_hat = z / R_scale = rho_over_R0 * e^{i phi}
    # log_val = -c * log(rho_over_R0) - i*c*phi + sum_k log(1 - z_hat / u_k)
    # where u_k = ((1 - s e^{i eta_k}) / R_scale_1)^n!
    z_hat = rho_over_R0 * cmath.exp(1j * phi)
    prod = complex(1.0, 0.0)
    log_scale = 0.0
    for k in range(n):
        eta = (2 * k * math.pi / n) if (c % 2 == 0) else ((2 * k + 1) * math.pi / n)
        if eta > math.pi:
            eta -= 2.0 * math.pi
        z_k = (1.0 - s * cmath.exp(1j * eta)) / R_scale_1
        # Compute u_k = z_k^n carefully using log:
        log_r = n * math.log(abs(z_k))
        ang = n * cmath.phase(z_k)
        if log_r > 500:
            # 1 - z_hat / u_k = 1
            continue
        elif log_r < -500:
            # 1 - z_hat / u_k ~ -z_hat / u_k
            prod *= -z_hat * cmath.exp(-1j * ang)
            log_scale += -log_r
        else:
            u_k = math.exp(log_r) * cmath.exp(1j * ang)
            fac = 1.0 - z_hat / u_k
            prod *= fac
            ap = abs(prod)
            if ap > 1e100 or (0 < ap < 1e-100):
                log_scale += math.log(ap)
                prod /= ap
    log_scale -= c * math.log(rho_over_R0)
    prod *= ((-1) ** j) * cmath.exp(-1j * c * phi)
    return prod, log_scale

def contour_quality(n, c, j, rho_over_R0, M=256, use_R0=True):
    """Compute exact ratio: Re(int I(phi) dphi) / int |I(phi)| dphi!
    If this ratio is > 0.5 (say +0.7), the integrand is overwhelmingly positive on the circle!
    """
    vals = []
    max_ls = -1e300
    for q in range(M):
        phi = 2.0 * math.pi * (q + 0.5) / M
        p, ls = eval_frozen_integrand(n, c, j, rho_over_R0, phi, use_R0=use_R0)
        vals.append((p, ls))
        if ls > max_ls:
            max_ls = ls
    s_re = 0.0
    s_abs = 0.0
    for p, ls in vals:
        w = math.exp(ls - max_ls) if (ls - max_ls > -700) else 0.0
        v = p * w
        s_re += v.real
        s_abs += abs(v)
    return s_re / s_abs if s_abs > 0 else 0.0, s_re / M, max_ls

print("=== REGIME 2 PROBE: CONTOUR POSITIVITY RATIO Re(int)/int|I| ACROSS (n, c, j) ===")
cases = [
    # Tiny kappa (n=5120, fixed c, j):
    (5120, 2, 3), (5120, 2, 4), (5120, 3, 4), (5120, 3, 5), (5120, 3, 6),
    (5120, 4, 5), (5120, 4, 6), (5120, 4, 8), (5120, 6, 9), (5120, 6, 12),
    # Moderate kappa (0.5 .. 2.49) at n = 320, 1000:
    (320, 2, 7), (320, 3, 8), (320, 4, 10), (320, 5, 11), (320, 6, 12),
    (320, 8, 12), (320, 10, 15), (320, 12, 16), (320, 15, 18),
    (1000, 10, 20), (1000, 15, 25), (1000, 20, 28), (1000, 30, 33),
]
for n, c, j in cases:
    kap = kappa_val(n, c, j)
    d = math.gcd(c, j)
    # Test several lam = rho / Rz values (or rho / R0):
    best_r = -2.0
    best_lam = None
    ratios = []
    for lam in [0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
        r, _, _ = contour_quality(n, c, j, lam, M=128, use_R0=False)
        ratios.append((lam, r))
        if r > best_r:
            best_r = r
            best_lam = lam
    # Also check exact Rz/R0 ratio:
    th_c, psi_j, alp1 = c*math.pi/n, j*math.pi/n, (j+c)*math.pi/n
    s = math.sin(psi_j) / math.sin(alp1)
    Rz_over_R0 = ((math.sin(th_c)/math.sin(alp1)) / (1.0 - s)) ** n
    r_at_09R0, _, _ = contour_quality(n, c, j, 0.9, M=128, use_R0=True)
    r_at_05R0, _, _ = contour_quality(n, c, j, 0.5, M=128, use_R0=True)
    r_at_15Rz, _, _ = contour_quality(n, c, j, 1.5, M=128, use_R0=False)
    print(f"n={n:4d} c={c:2d} j={j:2d} d={d:2d} kap={kap:.3f} Rz/R0={Rz_over_R0:6.2f} | "
          f"r(.5R0)={r_at_05R0:+.3f} r(.9R0)={r_at_09R0:+.3f} r(1.5Rz)={r_at_15Rz:+.3f} | "
          f"BEST(lam={best_lam:4.1f}*Rz)={best_r:+.4f}")
