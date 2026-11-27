"""
Test 100% EXACT FINITE-n Hayman saddle-point bounds for ALL D = 0 Major Arcs:
  c = d * q, j = d * p0 (gcd(p0, q) = 1), d >= 6, q*kap <= 0.70, c*q*kap >= 0.8, c*kap^2 <= 2.0, n >= 132.
Specifically tests:
  1. Exact finite-n peak phase |Phi(phi_nu)| across ALL q peaks phi_nu = (2*nu+1)*pi/q (remember finite-n cubic phase!).
  2. Exact global modulus decay ratio log|I(phi_nu + theta)| - log|I(phi_nu)| <= -c_dec * H * (1 - cos(q*theta))/q^2!
  3. Exact Hayman ratio M3 / H^{3/2} on the peak support!
  4. Exact finite-n positivity ratio Re(int) / int|I| AND explicit Hayman lower bound!
"""
import math

def exact_finite_n_roots(n, c, j):
    th_c = c * math.pi / n
    psi_j = j * math.pi / n
    alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1
    log_Rz = n * math.log(r3 / r1)
    roots = []
    shift = 0.0 if (c % 2 == 0) else 0.5
    sgn_ang = 0.0 if (j % 2 == 0) else math.pi
    k = 0
    while True:
        added = False
        klist = [0] if (k == 0 and c % 2 == 0) else ([0, -1] if (k == 0 and c % 2 == 1) else [k, -k - (1 if c % 2 == 1 else 0)])
        for kk in klist:
            eta = 2.0 * (kk + shift) * math.pi / n
            if abs(eta) > math.pi + 1e-12:
                continue
            z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
            log_r = log_Rz - n * math.log(abs(z))
            ang = sgn_ang - n * math.atan2(z.imag, z.real)
            if log_r > -55.0:
                roots.append((log_r, ang))
                added = True
        if not added and k > c:
            break
        k += 1
        if k > n // 2 + 2:
            break
    return roots, 2.0 * math.pi * r1 * r2 / r3

def derivs(roots, c, lam, phi):
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

def solve_saddle(roots, c, phi0):
    x = 0.0
    for _ in range(30):
        _, _, _, d1i, H, _, _ = derivs(roots, c, math.exp(x), phi0)
        if abs(d1i) < 1e-11 or H <= 0:
            break
        x -= max(-0.2, min(0.2, d1i / H))
    return math.exp(x)

def mk_n_exact(c, j, kap_target):
    # kap ~ 2*pi^2 * j*(j+c)/(c*n) -> n ~ 2*pi^2 * j*(j+c)/(c*kap_target)
    n = int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kap_target)))
    return max(132, n + (n % 2))

print("=== EXACT FINITE-n D=0 MAJOR-ARC HAYMAN STRESS TEST ===")
grid = [
    # (q, p0, d, kap_target)
    (1, 2, 7, 0.50),
    (1, 2, 7, 0.20),
    (1, 2, 10, 0.30),
    (1, 2, 12, 0.60),
    (1, 2, 20, 0.15),
    (1, 2, 50, 0.05),
    (1, 2, 200, 0.02),
    (2, 3, 7, 0.25),
    (2, 3, 10, 0.15),
    (2, 3, 30, 0.05),
    (2, 3, 100, 0.02),
    (3, 4, 7, 0.15),
    (3, 4, 15, 0.08),
    (3, 4, 50, 0.02),
    (5, 7, 7, 0.10),
    (5, 7, 20, 0.03),
    (7, 9, 7, 0.08),
    (10, 11, 7, 0.05),
]
worst_rat = 1e9
worst_pk_ph = 0.0
worst_hay = 0.0
min_Hred = 1e9
for q, p0, d, kt in grid:
    c = d * q
    j = d * p0
    n = mk_n_exact(c, j, kt)
    roots, kap = exact_finite_n_roots(n, c, j)
    sgn_pref = 1.0 if ((j * (c + 1) + c) % 2 == 0) else -1.0
    phi0 = math.pi if (q % 2 == 1) else (math.pi / q)
    lam = solve_saddle(roots, c, phi0)

    # Check ALL q peaks phi_nu = (2*nu+1)*pi/q:
    max_pk_ph = 0.0
    max_pk_d1r = 0.0
    max_pk_d1i = 0.0
    min_H = 1e9
    max_hay = 0.0
    for nu in range(q):
        pnu = (2 * nu + 1) * math.pi / q
        la, ph, d1r, d1i, H, d2i, d3i = derivs(roots, c, lam, pnu)
        ep = (ph + (0.0 if sgn_pref > 0 else math.pi) + math.pi) % (2.0 * math.pi) - math.pi
        max_pk_ph = max(max_pk_ph, abs(ep))
        max_pk_d1r = max(max_pk_d1r, abs(d1r))
        max_pk_d1i = max(max_pk_d1i, abs(d1i))
        min_H = min(min_H, H)
        max_hay = max(max_hay, abs(d3i) / (H ** 1.5))

    # Check exact integral ratio AND global decay constant c_dec on [0, pi/q]:
    M = max(256, 64 * q)
    vals = []
    mx = -1e300
    for m_i in range(M):
        phi = 2.0 * math.pi * (m_i + 0.5) / M
        la, ph, *_ = derivs(roots, c, lam, phi)
        ep = (ph + (0.0 if sgn_pref > 0 else math.pi) + math.pi) % (2.0 * math.pi) - math.pi
        vals.append((la, ep))
        if la > mx: mx = la
    rat = sum(math.exp(la - mx) * math.cos(ep) for la, ep in vals if la - mx > -700) / sum(math.exp(la - mx) for la, ep in vals if la - mx > -700)
    Hred = min_H / (q * q)
    # Check global decay on one sector [-pi/q, pi/q] around phi0:
    la0 = derivs(roots, c, lam, phi0)[0]
    min_c_dec = 1e9
    for t_i in range(1, 50):
        th = (t_i / 50.0) * (math.pi / q)
        lat = derivs(roots, c, lam, phi0 + th)[0]
        denom = (min_H / (q * q)) * 2.0 * (math.sin(0.5 * q * th) ** 2)
        if denom > 1e-9:
            min_c_dec = min(min_c_dec, (la0 - lat) / denom)

    worst_rat = min(worst_rat, rat)
    worst_pk_ph = max(worst_pk_ph, max_pk_ph)
    worst_hay = max(worst_hay, max_hay)
    min_Hred = min(min_Hred, Hred)
    print(f"q={q:2d} d={d:3d} n={n:6d} kap={kap:.4f} qk={q*kap:.3f} ck^2={c*kap*kap:.3f} | "
          f"lam={lam:.4f} H/q^2={Hred:6.2f} c_dec={min_c_dec:.3f} hay={max_hay:.4f} | "
          f"pk_ph={max_pk_ph:.4f} d1i={max_pk_d1i:.4f} | rat={rat:+.5f}")

print(f"\nSUMMARY: worst_rat={worst_rat:+.5f}, worst_pk_ph={worst_pk_ph:.5f}, worst_hay={worst_hay:.4f}, min_H/q^2={min_Hred:.2f}")
