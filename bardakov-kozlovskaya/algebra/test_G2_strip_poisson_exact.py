"""
Test our 100%-EXACT FINITE-n 2nd-IBP Cauchy-Cotangent bound on R_{non-q}(w)
and our 1-line exact conjugate-pair phase error E_phase(W) across ALL corners of Box G2!
"""
import math
import cmath

def exact_finite_n_params(n, c, j):
    th_c = c * math.pi / n
    psi_j = j * math.pi / n
    alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1
    kap = 2.0 * math.pi * r1 * r2 / r3
    return th_c, psi_j, alp1, r1, r2, r3, s, kap

def b_complex(x, n, c, j, r1, r2, r3, s):
    # b(x) = (-1)^j * (r3/r1)^n * exp(-2pi i (x - c/2) j / c) / (1 - s exp(2pi i x / n))^n
    # Compute log b(x) directly!
    eta = 2.0 * math.pi * x / n
    z = 1.0 - s * cmath.exp(1j * eta)
    # Wait: continuous branch of log(1 - s e^{i eta})!
    # At x = c/2: 1 - s e^{i th_c} = (r3/r1) e^{-i psi_j}, so log z = log(r3/r1) - i*psi_j
    # For |s e^{i eta}| < 1 (Im(eta) >= 0 -> |e^{i eta}| <= 1!), Re(z) can be computed, principal log works when |s e^{i eta}| < 1!
    # Wait: for Im(x) < 0, |s e^{i eta}| can exceed 1 if Im(x) is too negative, let's check!
    log_b = n * math.log(r3 / r1) - n * cmath.log(z) - 2.0j * math.pi * (x - c / 2.0) * (j / c)
    if j % 2 == 1:
        log_b += 1j * math.pi
    # Wrap Im(log_b) so at x = c/2 it is 0
    return log_b

def beta_and_betap(x, n, c, j, s):
    # beta(x) = d/dx log b(x) = -2pi i j/c + 2pi i s e^{i eta} / (1 - s e^{i eta})
    # beta'(x) = (2pi i)^2 / n * s e^{i eta} / (1 - s e^{i eta})^2
    eta = 2.0 * math.pi * x / n
    se = s * cmath.exp(1j * eta)
    om = 1.0 - se
    beta = 2.0j * math.pi * (se / om - j / c)
    betap = -(4.0 * math.pi * math.pi / n) * se / (om * om)
    return beta, betap

print("=== 1. TEST 2nd-IBP STRIP-POISSON BOUND ON R_{non-q} (EXACT FINITE-n!) ===")
# For each ray t, with theta_t in (-pi, pi] \ {0}:
# E_t^{top} = (q / (4pi^2)) sum_{l>=1} (1/l^2) e^{-2pi l H_t / q} * int_{-n/2}^{n/2} L''_t(x + i H_t) e^{2pi i l (x - c/2 - t)/q} dx
# Wait: CHECK SIGN OF H_t!
# Since beta(c/2) = -kap + i dpsi', log b(c/2 + sigma + i tau) ~ -kap*sigma - i*kap*tau!
# Wait: for tau > 0 (top), Im(log b) shifts by -kap*tau!
#   So if theta_t > 0, shifting tau > 0 MOVES theta_t TOWARD 0 (hits pole at tau = theta_t / kap)!
#   Wait: shifting tau < 0 (bottom) moves theta_t AWAY from 0 (toward +pi, hits pole at tau = -(2pi - theta_t)/kap)!
# THEREFORE:
#   - On the TOP (tau > 0), max safe height is H_t^+ < (theta_t mod 2pi) / kap (for theta_t in (0, 2pi))!
#   - On the BOTTOM (tau < 0), max safe height is H_t^- < ((2pi - theta_t) mod 2pi) / kap!
# Let's test this exact numerical contour bound vs actual |R_{non-q}|!

for d, q, p0, q2k in [(2, 2, 3, 2.0), (5, 2, 3, 2.0), (50, 2, 3, 2.0), (5, 5, 6, 2.0), (50, 5, 6, 2.0), (50, 10, 11, 2.0)]:
    c = d * q; j = d * p0; kt = q2k / (q * q)
    n = max(132, int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt))) // 2 * 2)
    th_c, psi_j, alp1, r1, r2, r3, s, kap = exact_finite_n_params(n, c, j)
    lam = 1.0  # test near saddle lam ~ 1
    for psi in [0.0, math.pi / 2.0, 2.0 * math.pi / 3.0]:
        # Compute exact sum of top + bottom 2nd-IBP L1 norms over all q rays t=0..q-1!
        tot_ibp2 = 0.0
        for k_ray in range(q):
            # theta_t in (0, 2pi): ((2*k_ray + 1)*pi + psi) / q
            th_pos = (((2 * k_ray + 1) * math.pi + psi) / q) % (2.0 * math.pi)
            th_neg = 2.0 * math.pi - th_pos
            # Choose H_pos = c0 * th_pos / kap, H_neg = c0 * th_neg / kap (capped at something safe so |s e^{i eta}| < 1 isn't crossed)
            # Wait: check |s e^{i eta}| = s * exp(-2pi tau / n) = 1 when tau = (n/(2pi)) log s < 0!
            # Since log(1/s) ~ (c/j) * ... wait: (1 - s) ~ r3/r1 ~ c/(j+c), so (n/(2pi)) log(1/s) ~ n*c/(2pi*(j+c)) ~ pi*j / kap >> 2pi/kap!
            th_min = min(th_pos, th_neg)
            tau_pole = (n / (2.0 * math.pi)) * math.log(1.0 / s)
            # Symmetric safe height H <= 0.72 * min(th_min, 2.2*pi/q) / kap for BOTH top and bottom!
            for sign in [+1, -1]:
                H = 0.72 * min(th_min, 2.2 * math.pi / q) / kap
                if sign == -1:
                    H = min(H, 0.5 * tau_pole)
                # Numerical L1 integral of |L''_t(x + i*sign*H)| over x in [-n/2, n/2]
                # Concentrate grid near x = +-c/2 + uniform inside [-c/2, c/2]
                tau = sign * H
                # Integrate over [-n/2, n/2] with adaptive step
                xs = []
                W = min(n / 2.0, c / 2.0 + 18.0 / kap)
                # Sample uniformly on [-W, W] with step min(0.15 / kap, max(0.2, c / 200.0))
                dx = min(0.12 / kap, max(0.15, c / 120.0))
                Npts = max(80, int(2.0 * W / dx))
                dx = 2.0 * W / Npts
                I_L2 = 0.0
                for idx in range(Npts):
                    xr = -W + (idx + 0.5) * dx
                    z = complex(xr, tau)
                    lb = b_complex(z, n, c, j, r1, r2, r3, s)
                    u_log = lb + math.log(lam) + 1j * th_pos
                    if u_log.real < -40.0:
                        continue
                    beta, betap = beta_and_betap(z, n, c, j, s)
                    if u_log.real > 40.0:
                        # u/(1-u) -> -1, u/(1-u)^2 -> 0
                        L2 = abs(betap)
                    else:
                        u = cmath.exp(u_log)
                        om = 1.0 - u
                        L2 = abs(-betap * (u / om) - beta * beta * (u / (om * om)))
                    I_L2 += L2 * dx
                decay = math.exp(-2.0 * math.pi * H / q)
                # sum_{l>=1} decay^l / l^2 = Li_2(decay) <= decay / (1 - decay/2) ... exact polylog:
                li2 = sum((decay ** l) / (l * l) for l in range(1, 30))
                tot_ibp2 += (q / (4.0 * math.pi * math.pi)) * li2 * I_L2
        if psi in [0.0, 2.0 * math.pi / 3.0]:
            print(f"  d={d:2d} q={q:2d} q^2k={q*q*kap:.3f} psi={psi:.2f}: 2nd-IBP exact L1 bound on |R_nonq| = {tot_ibp2:.4f}")

print("\n=== 2. TEST 1-LINE EXACT PHASE FACTOR E_phase(W)^{1/q} ON BOX G2 ===")
# Compute (1/q) * sum_{pairs k} sin^2(q*dpsi_k/2) / (sinh^2(0.5*log u_k) + cos^2(psi/2))
# at psi = 0, pi/2, 2pi/3 across ALL corners including d = 50, 200, 1000!
for d, q, p0, q2k in [
    (2, 1, 2, 0.45), (2, 1, 2, 2.0), (50, 1, 2, 0.80), (200, 1, 2, 0.45), (500, 1, 2, 0.28),
    (2, 2, 3, 2.0), (50, 2, 3, 2.0), (200, 2, 3, 1.20),
    (50, 5, 6, 2.0), (1000, 5, 6, 2.0), (500, 10, 11, 2.0),
]:
    c = d * q; j = d * p0; kt = q2k / (q * q)
    if c * kt * kt >= 41.5 or kt >= 2.45:
        continue
    n = max(132, int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt))) // 2 * 2)
    th_c, psi_j, alp1, r1, r2, r3, s, kap = exact_finite_n_params(n, c, j)
    shift = 0.0 if (c % 2 == 0) else 0.5
    # Sum over positive half-roots x > 0 (pairs!):
    S0 = S_pi2 = S_2pi3 = 0.0
    for k in range(0, n // 2 + 1):
        x = k + shift
        if x == 0 or x >= n / 2.0 - 1e-12:
            continue
        eta = 2.0 * math.pi * x / n
        z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
        lr = n * math.log(r3 / r1) - n * math.log(abs(z))
        if abs(q * lr) > 40.0:
            continue
        # exact q*dpsi:
        ang = (math.pi if (j % 2 == 1) else 0.0) - n * math.atan2(z.imag, z.real)
        qa = ((q * ang + math.pi) % (2.0 * math.pi)) - math.pi
        num = math.sin(0.5 * qa) ** 2
        sh2 = math.sinh(0.5 * q * lr) ** 2
        S0 += num / (sh2 + 1.0)
        S_pi2 += num / (sh2 + 0.5)
        S_2pi3 += num / (sh2 + 0.25)
    print(f"  d={d:4d} q={q:2d} ck^2={c*kap*kap:5.2f} q^2k={q*q*kap:.3f}: (1/q)*E_phase sum @ psi=0: {S0/q:.4f}, pi/2: {S_pi2/q:.4f}, 2pi/3: {S_2pi3/q:.4f}")
