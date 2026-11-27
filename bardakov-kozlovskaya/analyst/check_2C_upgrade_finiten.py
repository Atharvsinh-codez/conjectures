"""
1. EXACT FINITE-n STRESS TEST OF SUB-REGIME 2C-UPGRADE (c*q*kap <= C2 = 0.85 & 0.90, n >= 132):
   Verify EVERY component of the 1-page Exact q-Multisection + 3-Line Polya-Szego proof:
   (a) Automatic lower bound on n, upper bound on kap, upper bound on alpha_1, upper bound on (2j+c)/n!
   (b) Exact positive-root saddle Lambda_* = lam_*^q:
       - Verify max_k |a_k| * lam_* < 1 STRICTLY (ZERO exterior roots!)
       - Compute exact W_0 = |a_{-c/2}|^q * Lambda_* < 1!
       - Compute exact H_q(Lambda_*) / g >= 0.65 (so H_q >= 0.65 > 0.50 for ALL g >= 1!)
   (c) Exact non-q power-sum tail sum_{q!|r} |tilde_p_r| * lam_*^r / r on finite n!
   (d) Exact finite-n phase perturbation on the q-multiple part (since a_k^q has tiny finite-n arg(a_k^q) != 0)!
2. Also inspect the exact borderline survivors of B_minor_3way (checking d_min for D=0 and min(ct/D^2) for D!=0)!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd, lgamma

# =============================================================================
# 1. FINITE-n STRESS TEST OF SUB-REGIME 2C-UPGRADE (c*q*kap <= 0.90)
# =============================================================================
print("=== 1. EXACT FINITE-n STRESS TEST OF SUB-REGIME 2C-UPGRADE (c*q*kap <= 0.90) ===")
worst_W0 = 0.0; worst_Hq_g = 1e9; worst_L1_nonq = 0.0; worst_q_ang = 0.0
worst_cases = {}

for c in [2, 3, 4, 5, 6, 8, 10, 15, 20, 30]:
    for j_off in range(1, min(c + 4, 25)):
        j = c + j_off
        g = gcd(c, j); q = c // g
        for target_cqk in [0.15, 0.40, 0.70, 0.80, 0.85, 0.90]:
            target_kap = target_cqk / (c * q)
            # Pick even n >= 132 achieving approximately target_kap:
            # kap ~ 2*pi^2 * j*(j+c) / (c * n) => n ~ 2*pi^2 * j*(j+c) / (c * target_kap)
            n_est = int(round(2.0 * pi * pi * j * (j + c) / (c * target_kap)))
            n = max(132, 2 * ((n_est + 1) // 2))
            if j > n // 2 - (c + 1) // 2 - 1: continue
            kap = 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)
            cqk = c * q * kap
            if cqk > 0.905 or kap >= 2.5: continue

            rj = sin(j * pi / n) / sin((j + c) * pi / n)
            sj = sin(c * pi / n) / sin((j + c) * pi / n)
            # Compute log |a_k| and arg(a_k^q) for the active roots around theta = pi (or near -c/2):
            # Roots are at eta_k = (2k+1)pi/n (c odd) or 2k pi/n (c even), near eta = 0!
            # Recall a_k = (-1)^j * sj^n / (1 - rj e^{i eta_k})^n!
            # Wait: check sign convention! At eta_k = +-theta_c: |1 - rj e^{+-i theta_c}| = sj, so |a_k| = 1!
            # Peak modulus is at eta = 0: |1 - rj| = (sin((j+c)pi/n) - sin(j pi/n))/sin((j+c)pi/n) < sj!
            M_k = min(n // 2, int(c + sqrt(45.0 * c / max(kap, 1e-8))) + 10)
            if c % 2 == 0:
                ks = np.arange(-M_k, M_k + 1, dtype=np.float64)
                etas = 2.0 * ks * pi / n
            else:
                ks = np.arange(-M_k, M_k, dtype=np.float64)
                etas = (2.0 * ks + 1.0) * pi / n
            la = n * np.log(sj / (1.0 - rj * np.exp(1j * etas))) + (1j * pi if j % 2 == 1 else 0.0)
            re_a = la.real
            # Note: exact peak modulus |a_max| = exp(max(re_a))
            max_re_a = float(np.max(re_a))
            # Solve exact positive-root saddle: sum_k 1 / (1 + exp(-(q*re_a + V))) = c!
            lo, hi = -80.0, 20.0
            for _ in range(60):
                vm = 0.5 * (lo + hi)
                z = np.clip(q * re_a + vm, -70.0, 70.0)
                if float(np.sum(1.0 / (1.0 + np.exp(-z)))) > c: hi = vm
                else: lo = vm
            V = 0.5 * (lo + hi)
            lam = exp(V / q)
            W0 = exp(q * max_re_a + V) # MUST BE < 1 for ZERO exterior roots!
            z = np.clip(q * re_a + V, -70.0, 70.0)
            p = 1.0 / (1.0 + np.exp(-z))
            H_q = float(np.sum(p * (1.0 - p))) / q

            # Check non-q L1 tail sum_{r=1..50, q!|r} |tilde_p_r| * lam^r / r!
            # Use exact roots u_k = exp(re_a + log(lam) + i*la.imag):
            L1_nonq = 0.0
            if W0 < 1.0 and q >= 2:
                u = np.exp(la + log(lam))
                pu = np.ones_like(u)
                for r in range(1, min(60, 6 * q + 10)):
                    pu *= u
                    if r % q != 0:
                        L1_nonq += abs(float(np.sum(pu).real)) / r

            # Check finite-n phase error on q-multiples: arg(a_k^q) weighted by p*(1-p)/q / sqrt(H_q)!
            # Note: at the saddle W = -lam^q, the positive-root factor is (1 + |a_k|^q lam^q e^{i theta})^{1/q},
            # while actual factor is (1 + a_k^q lam^q e^{i theta})^{1/q} where a_k^q = |a_k|^q e^{i q arg(a_k)}!
            # Wait: check sign of a_k^q!
            # At eta_k = +-theta_c: a_{+-theta_c} = (-1)^j * sj^n / (sj e^{+-i j pi/n})^n = (-1)^j * (-1)^j = +1!
            # Wait: what is arg(a_k^q) at eta = 0 (m = -c/2)?
            # At eta = 0: 1 - rj > 0, so a_k = (-1)^j * |a_k|!
            # Wait: (-1)^{j*q}: if j*q is ODD (e.g. c even, g even, q odd? No: q = c/g, j = g * p0, so j*q = c * p0!
            # Wait: if c is odd and p0 is odd, j*q is ODD, so a_k^q has phase (-1)^{c*p0} = -1?!)
            # CHECK THAT! In our model earlier, im_a = -2*pi*m*j/c had m = 0 at boundary, so m = -c/2 had -2*pi*(-c/2)*j/c = +pi*j!
            # Multiply by q: +pi*j*q = pi * c * p0! Wait: at m = 0, a_0 = +1, so a_0^q = +1!
            # Why did eta = 0 give (-1)^j while m = 0 gave +1?
            # Because for c odd, eta = 0 is NOT a grid point (ks are half-integers, eta_k = +-pi/n, ...)!
            # Let's check actual arg(a_k^q) modulo 2*pi directly!
            ang_q = (q * la.imag + pi) % (2 * pi) - pi
            # Wait: check if a_k^q is near +1 or -1!
            # At boundary eta = theta_c: a = +1 => a^q = +1 ALWAYS!
            w_ang = float(np.sum(p * (1.0 - p) * np.abs(ang_q))) / (q * sqrt(H_q))

            key = round(target_cqk, 2)
            if key not in worst_cases:
                worst_cases[key] = [0.0, 1e9, 0.0, 0.0, None]
            wc = worst_cases[key]
            if W0 > wc[0]: wc[0] = W0
            if H_q / g < wc[1]: wc[1] = H_q / g
            if L1_nonq > wc[2]: wc[2] = L1_nonq; wc[4] = (n, c, j, q, g, cqk)
            if w_ang > wc[3]: wc[3] = w_ang

print(f"{'cqk':>5} | {'max W0':>8} | {'min Hq/g':>8} | {'max L1_nonq':>11} | {'max w_ang/sH':>12} | worst L1 (n,c,j,q,g,cqk)")
for k in sorted(worst_cases):
    w0, hqg, l1, wa, info = worst_cases[k]
    print(f"{k:5.2f} | {w0:8.5f} | {hqg:8.5f} | {l1:11.6f} | {wa:12.6f} | {info}")
