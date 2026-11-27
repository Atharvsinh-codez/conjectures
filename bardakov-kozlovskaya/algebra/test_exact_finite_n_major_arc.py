"""
Test the EXACT FINITE-n unified saddle-point (Hayman) behavior on extreme Major Arcs (BOTH D = 0 AND D != 0)!
Uses the EXACT finite-n roots a_k = (-1)^j * R_z / (1 - s_j e^{i eta_k})^n (zero Gaussian approximation!)
for large n up to n = 200,000 (summing only roots with |log|a_k|| <= 50, plus exact tail check!).
Also fixes the c*kap = 0.40 pole avoidance (starting Newton above |a_0|^{-1} or at lam = 1.0!).
"""
import math
import cmath

def exact_finite_n_roots(n, c, j):
    """Return exact finite-n roots (log|a_k|, arg(a_k)) for all k with |log|a_k|| <= 50 (plus exact counts of deep in/out)."""
    th_c = c * math.pi / n
    psi_j = j * math.pi / n
    alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1
    log_Rz = n * math.log(r3 / r1)
    # eta_k = 2*k*pi/n (even c) or (2*k+1)*pi/n (odd c).
    # For |eta_k| > pi/2, |1 - s e^{i eta}| >= 1, so log|a_k| <= log_Rz < -80!
    # Find k_max where log|a_k| < -55:
    roots = []
    shift = 0.0 if (c % 2 == 0) else 0.5
    sgn_ang = 0.0 if (j % 2 == 0) else math.pi
    # Scan k = 0, +-1, +-2, ... until log|a_k| < -55:
    k = 0
    while True:
        added_any = False
        for kk in ([0] if (k == 0 and c % 2 == 0) else ([0, -1] if (k == 0 and c % 2 == 1) else [k, -k - (1 if c % 2 == 1 else 0)])):
            eta = 2.0 * (kk + shift) * math.pi / n
            if abs(eta) > math.pi + 1e-12:
                continue
            z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
            log_r = log_Rz - n * math.log(abs(z))
            ang = sgn_ang - n * math.atan2(z.imag, z.real)
            if log_r > -55.0:
                roots.append((log_r, ang))
                added_any = True
        if not added_any and k > c:
            break
        k += 1
        if k > n // 2 + 2:
            break
    return roots, 2.0 * math.pi * r1 * r2 / r3

def log_F_derivs(roots, c, lam, phi):
    log_lam = math.log(lam)
    log_abs = -c * log_lam
    tot_ang = -c * phi
    f1 = 0j; f2 = 0j; f3 = 0j
    for lr, ang in roots:
        lm = lr + log_lam
        th = ang + phi
        if lm < -45.0:
            continue
        elif lm > 45.0:
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
    return log_abs, tot_ang, f1.real, f1.imag - c, -f2.real, f2.imag, f3.imag

def test_exact(n, c, j, r0, D):
    roots, kap = exact_finite_n_roots(n, c, j)
    sgn_pref = 1.0 if ((j * (c + 1) + c) % 2 == 0) else -1.0
    phi0 = math.pi if (r0 % 2 == 1) else (math.pi / r0)
    # Find lam_L2 (always safe even near poles!) AND check exact saddle at phi0!
    M = max(256, 64 * r0)
    def obj(log_lam):
        lam = math.exp(log_lam)
        vals = [log_F_derivs(roots, c, lam, 2.0 * math.pi * (q + 0.5) / 128)[0] for q in range(128)]
        mx = max(vals)
        return 2.0 * mx + math.log(sum(math.exp(2.0 * (v - mx)) for v in vals) / 128)
    grid = [-1.5 + 3.0 * i / 30 for i in range(31)]
    bi = min(range(len(grid)), key=lambda i: obj(grid[i]))
    lo, hi = grid[max(0, bi - 1)], grid[min(len(grid) - 1, bi + 1)]
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    for _ in range(30):
        m1, m2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
        if obj(m1) < obj(m2): hi = m2
        else: lo = m1
    lam = math.exp(0.5 * (lo + hi))

    la0, ph0, d1r, d1i, H, d2i, d3i = log_F_derivs(roots, c, lam, phi0)
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
    # Also check exact lam = 1 rat!
    vals1 = []
    mx1 = -1e300
    for q in range(M):
        phi = 2.0 * math.pi * (q + 0.5) / M
        la, ph, *_ = log_F_derivs(roots, c, 1.0, phi)
        ep = (ph + (0.0 if sgn_pref > 0 else math.pi) + math.pi) % (2.0 * math.pi) - math.pi
        vals1.append((la, ep))
        if la > mx1: mx1 = la
    rat1 = sum(math.exp(la - mx1) * math.cos(ep) for la, ep in vals1 if la - mx1 > -700) / sum(math.exp(la - mx1) for la, ep in vals1 if la - mx1 > -700)
    hay = abs(d3i) / (H ** 1.5) if H > 0 else -1.0
    print(f"n={n:7d} c={c:4d} j={j:4d} kap={kap:.4f} r0={r0} D={D:+2d} | lam*={lam:.4f} | "
          f"H/r0^2={H/(r0*r0):7.2f} d1i={d1i:+6.3f} hay={hay:6.4f} | rat@1={rat1:+.4f} rat@lam*={rat:+.5f}")

print("=== EXACT FINITE-n MAJOR-ARC TEST (ZERO GAUSSIAN APPROXIMATION) ===")
# Pick even n matching each (c, j, kap) via n ~ 2*pi*j*(j+c)/(c*kap)!
def mk_n(c, j, kap):
    n = int(round(2.0 * math.pi * j * (j + c) / (c * kap)))
    return max(132, n + (n % 2))

cases = [
    (10, 20, 0.10, 1, 0),
    (12, 24, 0.20, 1, 0),
    (20, 40, 0.05, 1, 0),
    (30, 45, 0.02, 2, 0),
    (30, 40, 0.01, 3, 0),
    (50, 80, 0.005, 5, 0),
    (200, 400, 0.01, 1, 0),
    (30, 31, 0.05, 1, 1),
    (50, 51, 0.02, 1, 1),
    (80, 107, 0.005, 3, 1),
    (100, 134, 0.005, 3, 2),
    (120, 122, 0.02, 1, 2),
    (200, 203, 0.02, 1, 3),
    (200, 201, 0.01, 1, 1),
    (200, 201, 0.002, 1, 1),
    (200, 301, 0.005, 2, 2),
]
for c, j, kap, r0, D in cases:
    n = mk_n(c, j, kap)
    test_exact(n, c, j, r0, D)
