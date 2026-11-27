"""
Test the UNIFIED SADDLE-POINT (Hayman) mechanism on ALL Major Arcs (BOTH D = 0 AND D != 0, including r0 nmid c)!
For each major arc (c, j, kap, r0, D):
  1. Solve the 1D real saddle equation S(lam) := Re sum_k (a_k w0 / (a_k w0 - 1)) = c at w0 = lam * exp(i * phi0),
     where phi0 = pi / r0 (or pi if r0 odd — check BOTH!).
  2. Compute:
     - exact curvature H = -d^2/dtheta^2 log|I(phi0 + theta)| at theta = 0,
     - exact |Phi(phi_nu)| at all r0 peaks phi_nu = (2*nu + 1)*pi/r0,
     - exact |Phi'(phi_nu)| and |Phi''(phi_nu)| at all r0 peaks,
     - exact third derivative bound M3 and Hayman ratio M3 / H^{3/2},
     - exact tail mass outside |theta| <= theta_cut (where |Phi| <= 1.2) vs main peak mass!
"""
import math
import cmath

def gaussian_a_list(c, j, kap):
    alpha = math.pi * kap / c
    x_max = int(math.sqrt(0.25 * c * c + 45.0 / alpha)) + 5
    shift = 0.0 if (c % 2 == 0) else 0.5
    roots = []
    sgn = 1.0 if (j % 2 == 0) else -1.0
    for k in range(-x_max, x_max + 1):
        x = k + shift
        log_r = alpha * (0.25 * c * c - x * x)
        ang = 2.0 * math.pi * x * j / c + (0.0 if sgn > 0 else math.pi)
        roots.append((log_r, ang))
    return roots

def log_F_derivs(roots, c, lam, phi):
    """Exact log|I|, Phi, d/dphi log|I|, dPhi/dphi, d^2/dphi^2 log|I|, d^2Phi/dphi^2, d^3Phi/dphi^3 at w = lam*e^{i phi}.
    Using z = a_k * w:
      f(phi) = log(1 - z)
      f'(phi) = -i * z / (1 - z)
      f''(phi) = z / (1 - z)^2
      f'''(phi) = i * z * (1 + z) / (1 - z)^3
    Since log(I) = -c*log(lam) - i*c*phi + sum f(phi):
      Re(f') = d/dphi log|I|,   Im(f') - c = dPhi/dphi
      Re(f'') = d^2/dphi^2 log|I| = -H,  Im(f'') = d^2Phi/dphi^2
      Im(f''') = d^3Phi/dphi^3!
    """
    log_lam = math.log(lam)
    log_abs = -c * log_lam
    tot_ang = -c * phi
    f1 = 0j
    f2 = 0j
    f3 = 0j
    for lr, ang in roots:
        lm = lr + log_lam
        th = ang + phi
        if lm < -40.0:
            continue
        elif lm > 40.0:
            # 1 - z ~ -z -> log(1-z) = log(-z), f' = i, f'' = 0, f''' = 0
            log_abs += lm
            tot_ang += th + math.pi
            f1 += 1j
        else:
            r = math.exp(lm)
            z = complex(r * math.cos(th), r * math.sin(th))
            omz = 1.0 - z
            log_abs += math.log(abs(omz))
            tot_ang += math.atan2(omz.imag, omz.real)
            u = z / omz
            f1 += -1j * u
            f2 += u / omz
            f3 += 1j * u * (1.0 + z) / (omz * omz)
    d1_re = f1.real
    d1_im = f1.imag - c
    H = -f2.real
    d2_im = f2.imag
    d3_im = f3.imag
    return log_abs, tot_ang, d1_re, d1_im, H, d2_im, d3_im

def solve_saddle_lam(roots, c, phi0):
    """Solve dPhi/dphi(lam, phi0) = 0 for lam near 1 (specifically where H > 0!)."""
    # Recall dPhi/dphi = Im(sum -i*z/(1-z)) - c = Re(sum z/(z-1)) - c!
    # Wait: d/d(log lam) [dPhi/dphi] = Re(sum -z/(1-z)^2) = H! So when H > 0, dPhi/dphi is STRICTLY INCREASING in log(lam)!!
    # Check: z = a*lam*e^{i phi0}, dz/d(log lam) = z, d/d(log lam) [z/(z-1)] = -z/(z-1)^2 = -z/(1-z)^2 = -f2!
    # YES! d/d(log lam) (dPhi/dphi) = -Re(f2) = +H!!
    # Start near lam = 1.0 (or scan grid to avoid any pole if D odd and |a0|*lam ~ 1):
    best_x = 0.0
    best_val = 1e9
    for k in range(-40, 41):
        x = k * 0.02
        _, _, _, d1_im, H, _, _ = log_F_derivs(roots, c, math.exp(x), phi0)
        if H > 0 and abs(d1_im) < best_val:
            best_val = abs(d1_im)
            best_x = x
    x = best_x
    for _ in range(25):
        _, _, _, d1_im, H, _, _ = log_F_derivs(roots, c, math.exp(x), phi0)
        if abs(d1_im) < 1e-11:
            break
        if H <= 0:
            break
        step = max(-0.2, min(0.2, d1_im / H))
        x -= step
    return math.exp(x)

def test_case(c, j, kap, r0, D):
    roots = gaussian_a_list(c, j, kap)
    sgn_pref = 1.0 if ((j * (c + 1) + c) % 2 == 0) else -1.0
    # Use phi0 = pi if r0 is odd, else pi/r0:
    phi0 = math.pi if (r0 % 2 == 1) else (math.pi / r0)
    lam = solve_saddle_lam(roots, c, phi0)
    la0, ph0, d1r, d1i, H, d2i, d3i = log_F_derivs(roots, c, lam, phi0)
    eff_ph0 = (ph0 + (0.0 if sgn_pref > 0 else math.pi) + math.pi) % (2.0 * math.pi) - math.pi

    # Check ALL r0 peaks phi_nu = (2*nu + 1)*pi/r0:
    max_pk_ph = 0.0
    max_pk_d1i = 0.0
    min_pk_H = 1e9
    for nu in range(r0):
        pnu = (2 * nu + 1) * math.pi / r0
        _, phn, _, d1in, Hn, _, _ = log_F_derivs(roots, c, lam, pnu)
        ep = (phn + (0.0 if sgn_pref > 0 else math.pi) + math.pi) % (2.0 * math.pi) - math.pi
        max_pk_ph = max(max_pk_ph, abs(ep))
        max_pk_d1i = max(max_pk_d1i, abs(d1in))
        min_pk_H = min(min_pk_H, Hn)

    # Check exact integral ratio at THIS explicit saddle lam!
    M = max(256, 64 * r0)
    vals = []
    mx = -1e300
    for q in range(M):
        phi = 2.0 * math.pi * (q + 0.5) / M
        la, ph, *_ = log_F_derivs(roots, c, lam, phi)
        ep = (ph + (0.0 if sgn_pref > 0 else math.pi) + math.pi) % (2.0 * math.pi) - math.pi
        vals.append((la, ep))
        if la > mx: mx = la
    s_re = sum(math.exp(la - mx) * math.cos(ep) for la, ep in vals if la - mx > -700)
    s_abs = sum(math.exp(la - mx) for la, ep in vals if la - mx > -700)
    rat = s_re / s_abs
    # Hayman ratio |d3i| / H^{3/2} and peak width r0 * sqrt(H):
    hayman = abs(d3i) / (H ** 1.5) if H > 0 else 999.0
    # Scaled curvature H / r0^2:
    H_red = H / (r0 * r0)
    print(f"c={c:4d} j={j:4d} kap={kap:.4f} r0={r0:2d} D={D:+2d} | lam*={lam:.4f} | "
          f"H/r0^2={H_red:7.2f} |d3|/H^1.5={hayman:.4f} | max|Phi_pk|={max_pk_ph:.4f} max|d1i_pk|={max_pk_d1i:.4f} | rat={rat:+.5f}")

print("=== UNIFIED SADDLE-POINT (HAYMAN) TEST ACROSS MAJOR ARCS (D=0 AND D!=0) ===")
tests = [
    # D = 0 (smallest c/r0 and largest c):
    (10, 20, 0.100, 1, 0),
    (12, 24, 0.200, 1, 0),
    (20, 40, 0.050, 1, 0),
    (30, 45, 0.020, 2, 0),
    (30, 40, 0.010, 3, 0),
    (50, 80, 0.005, 5, 0),
    (200, 300, 0.005, 2, 0),
    (200, 400, 0.005, 1, 0),
    (500, 1000, 0.010, 1, 0),
    # D != 0 (smallest c/(r0*D^2), r0 nmid c, and tiny c*kap):
    (30, 31, 0.050, 1, 1),
    (50, 51, 0.020, 1, 1),
    (80, 107, 0.005, 3, 1),
    (80, 107, 0.010, 3, 1),
    (100, 134, 0.005, 3, 2),
    (120, 122, 0.020, 1, 2),
    (200, 203, 0.020, 1, 3),
    (200, 201, 0.010, 1, 1),
    (200, 201, 0.002, 1, 1),
    (200, 301, 0.005, 2, 2),
    (300, 401, 0.002, 3, 1),
    (500, 701, 0.001, 5, 1),
]
for args in tests:
    test_case(*args)
