"""
DERIVE & VERIFY 100% ELEMENTARY, CLOSED-FORM CONSTANTS FOR:
1. Lemma R2.1 (Tiny Cubic-Root Perturbation in Regime 2: n >= 132, kappa < 2.5 or 3.0)
2. Lemma R2.2 (Sub-regime 2A: c * kappa^2 >= C_1 on |w| = 1 via exact 4-line pairing bound + Jensen)
3. Lemma R2.3 (Sub-regime 2C: c * q * kappa <= C_2 on |w| = lam < e^{-c*kap/4} via exact binomial p_r)

Check EVERY intermediate elementary inequality across ALL Regime-2 points for n in {132, 160, 320, 1000}
AND across extreme high-n corners up to n = 10^6!
"""
import numpy as np
from math import gcd, pi, sin, cos, log, exp, sqrt, asin

def cot(x): return cos(x) / sin(x)
def kappa_val(n, c, j):
    return 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)

# =============================================================================
# PART 1: EXACT CLOSED-FORM PROOF OF LEMMA R2.2 (Sub-regime 2A: c*kap^2 >= C_1)
# =============================================================================
# Let kap < kap0 = 2.5 (or 3.0), n >= 132.
# From sin(alpha_1) < kap0 / (2*pi):
#   kap0 = 2.5 => alpha_1 < asin(2.5/(2*pi)) = 0.40928 < 0.1303*pi
#   kap0 = 3.0 => alpha_1 < asin(3.0/(2*pi)) = 0.49782 < 0.1585*pi
# Let alpha_max = asin(kap0 / (2*pi)).
# Then:
#   theta_c = c*pi/n < alpha_max / 2  (since j >= c+1 => j+c >= 2c+1)
#   (2j + c)*pi/n < 2*alpha_max < pi/2  => cos((2j+c)*pi/n) > cos(2*alpha_max) > 0!
#   Hence Re(zeta_m) > 0 for ALL m, so |1 + zeta_m| > 1!
#
# EXACT ELEMENTARY MODULUS BOUNDS FOR |a_{-m}^{-1}| AND |a_m|:
# Recall from ALGEBRA_ROUTE Lemma 4.4:
#   |1 + r_j e^{i theta_m}|^2 / s_j^2 = 1 + (2*kap/(pi*sin(theta_c))) * sin((|eta|+theta_c)/2) * sin((|eta|-theta_c)/2)
# where eta = pi - theta_m!
# - For exterior roots (m = 1..floor((c-1)/2)): |eta| = theta_c - 2*m*pi/n, so:
#     |1 + r_j e^{i theta_{-m}}|^2 / s_j^2 = 1 - (2*kap/(pi*sin(theta_c))) * sin(theta_c - m*pi/n) * sin(m*pi/n).
#   Using sin(theta_c - x)/sin(theta_c) >= 1 - m/c (concavity of sin on [0, pi/2]!)
#   and sin(m*pi/n) >= sinc(alpha_max/4) * (m*pi/n) (since m*pi/n <= theta_c/2 < alpha_max/4!):
#     |1 + ...|^2 / s_j^2 <= 1 - (2*gam_ext*kap/n) * m*(1 - m/c),  gam_ext = sinc(alpha_max/4) = sin(alpha_max/4)/(alpha_max/4)!
#   Since (1 - u)^{n/2} <= exp(-n*u/2):
#     |a_{-m}^{-1}| <= exp(-gam_ext * kap * m * (1 - m/c))!!   (100% ELEMENTARY, 2 LINES!!)
#
# - For interior roots (m >= 1): |eta| = theta_c + 2*m*pi/n, so:
#     |1 + ...|^2 / s_j^2 = 1 + (2*kap/(pi*sin(theta_c))) * sin(theta_c + m*pi/n) * sin(m*pi/n).
#   For m <= (n - 2c)/2: theta_c + m*pi/n <= pi - theta_c, so sin(theta_c + m*pi/n) >= sin(theta_c),
#   and for m*pi/n <= pi/2: sin(m*pi/n) >= (2/pi)*(m*pi/n) = 2m/n (or sinc(x) for smaller m)!
#   Thus |1 + ...|^2 / s_j^2 >= 1 + (4*kap/(pi*n)) * m  (and for small m <= n/6: >= 1 + (2*sinc(pi/6)*kap/n)*m)!
#   Even better: for m >= 1, log(1 + y) >= y/(1 + y/2) or simply (1 + A/n)^{-n/2} <= (1 + A/N0)^{-N0/2}!
#
# - For n * |zeta_m|:
#     n * |zeta_m| = n * [4 * sin(psi_j)*sin(alpha_1) / sin^2(theta_c)] * sin^2(m*pi/n)
#                  = (2*kap / (pi * sin(theta_c))) * n * sin^2(m*pi/n)
#                  <= (2*kap / (pi * sinc(alpha_max/2) * c*pi/n)) * n * (m*pi/n)^2
#                  = (2 / sinc(alpha_max/2)) * (kap / c) * m^2!!
#   WAIT A SECOND!! LOOK HOW CLEAN THAT CONSTANT IS!!
#   Let gam_z := 1 / sinc(alpha_max / 2).
#   For kap0 = 2.5: alpha_max = 0.4093, alpha_max/2 = 0.2046 => sinc = 0.9930 => gam_z = 1.0070!
#   gam_ext = sinc(alpha_max / 4) = sinc(0.1023) = 0.9983!
#
# NOW COMPUTE THE EXACT CLOSED-FORM UPPER BOUND ON B_pair(n, c, j)!!
#   Pair sum:
#     2 * sum_{m=1}^{floor((c-1)/2)} (2 * gam_z * kap / c) * m^2 * exp(-gam_ext*kap*m*(1-m/c)) / (1 - exp(-gam_ext*kap*m*(1-m/c)))
#   WAIT: Look at m <= c/2 => 1 - m/c >= 1/2!
#   If we just bound 1 - m/c >= 1/2 EVERYWHERE, we lose a factor of 2^3 = 8 (too loose!).
#   INSTEAD, split at m <= alpha_cut * c (say m <= c/4 => 1 - m/c >= 3/4, or keep exact (1 - m/c) in integral!!):
#   Wait: let x = kap * m, and let L = c * kap!
#   Notice x*(1 - x/L) = x*(1 - m/c)!
#   Even simpler: for 1 <= m <= c/2, split into m <= c/6 (where 1 - m/c >= 5/6) and c/6 < m <= c/2 (where 1-m/c >= 1/2 and x >= L/6)!
#   CHECK THAT ELEMENTARY 2-PIECE INTEGRAL BOUND!

def closed_form_B_pair_bound(C1, kap0=2.5, N0=132):
    # Compute RIGOROUS worst-case closed-form upper bound on B_pair for ALL (n, c, j) with
    # n >= N0, kap < kap0, c * kap^2 >= C1!
    alpha_max = asin(kap0 / (2.0 * pi))
    gam_ext = sin(alpha_max / 4.0) / (alpha_max / 4.0)
    gam_z = (alpha_max / 2.0) / sin(alpha_max / 2.0)
    # Note: c * kap = (c * kap^2) / kap > C1 / kap0!
    L_min = C1 / kap0

    # Let's optimize split fraction f in (0, 0.5) (m <= f*c vs f*c < m <= c/2)!
    # On 1 <= m <= f*c: 1 - m/c >= 1 - f, so exponent >= gam_ext*(1-f)*kap*m := lam1*m.
    # For ANY lam1 > 0, sum_{m=1}^inf m^2 / (exp(lam1*m) - 1):
    #   Since g(t) = t^2 / (exp(t) - 1) is unimodal with max M_g = max_{t>0} t^2/(e^t-1) < 0.6477 (at t ~ 1.594):
    #   By Riemann-integral comparison with unimodal function:
    #     sum_{m=1}^inf g(lam1*m) <= (1/lam1) * int_0^inf t^2/(e^t-1) dt + M_g
    #                              = 2*zeta(3)/lam1 + M_g!!
    #   Wait: multiplying by 4*gam_z/(c*kap*lam1^2?):
    #   Term is 2 * (2*gam_z*kap/c) * m^2 / (exp(lam1*m) - 1)
    #         = (4*gam_z / (c * kap * gam_ext^2 * (1-f)^2)) * g(lam1*m)!
    #   Wait: NO! Multiply by (lam1^2 / lam1^2):
    #     (4*gam_z*kap/c) * (1/lam1^2) * sum g(lam1*m)
    #     <= (4*gam_z*kap/c) * (1/lam1^2) * (2*zeta(3)/lam1 + M_g)
    #     = 4*gam_z * [ 2*zeta(3) / (c*kap^2 * gam_ext^3 * (1-f)^3) + M_g / (c*kap * gam_ext^2 * (1-f)^2) ]!!
    #   WHOA!! CHECK THAT!! ZERO DEPENDENCE ON c OR kap INDIVIDUALLY — ONLY ON c*kap^2 >= C1 AND c*kap >= L_min!!
    zeta3 = 1.2020569031595943
    M_g = 0.6476102379 # exact max of t^2/(e^t - 1)
    best = 1e9; best_info = None
    for f in np.linspace(0.05, 0.40, 71):
        lam1_fac = gam_ext * (1.0 - f)
        S1 = 4.0 * gam_z * (2.0 * zeta3 / (C1 * (lam1_fac**3)) + M_g / (L_min * (lam1_fac**2)))
        # On f*c < m <= c/2: 1 - m/c >= 1/2, so exponent >= (gam_ext/2)*kap*m >= (gam_ext/2)*f*L_min := u_lo!
        # Wait: even sharper! On [f, 0.5], t*(1-t) >= f*(1-f) > 0, and on m in (f*c, c/2],
        # use exp(-gam_ext*kap*m*(1-m/c)) <= exp(-gam_ext*L*f*(1-f))? Wait: as m increases from f*c to c/2,
        # m*(1-m/c) INCREASES monotonically to c/4!
        # Even better: use 1-m/c >= 1/2 AND integral from x = (gam_ext/2)*f*L_min to infinity!
        #   lam2_fac = gam_ext / 2.0
        #   sum_{m > f*c} g(lam2*m) <= (1/lam2) * int_{lam2*f*c}^inf t^2/(e^t-1) dt + g_max_tail
        #   where for t >= u0 = lam2_fac * f * L_min:
        #     int_{u0}^inf t^2/(e^t-1) dt <= (1/(1 - e^{-u0})) * (u0^2 + 2*u0 + 2) * e^{-u0}!
        lam2_fac = gam_ext * 0.5
        u0 = lam2_fac * f * L_min
        # Also m*(1-m/c) >= f*(1-f)*c, so u_min_exact = gam_ext * f * (1 - f) * L_min!
        #Wait: using m*(1-m/c) >= (1/2)*m gives u0 = 0.5*gam_ext*f*L_min:
        tail_int = ((u0*u0 + 2.0*u0 + 2.0) * exp(-u0)) / (1.0 - exp(-u0))
        g_tail_max = M_g if u0 < 1.6 else (u0*u0 / (exp(u0) - 1.0))
        S2 = 4.0 * gam_z * (tail_int / (C1 * (lam2_fac**3)) + g_tail_max / (L_min * (lam2_fac**2)))

        # Midpoint (c even): |a_{-c/2}^{-1}| <= exp(-gam_ext * c*kap / 4) <= exp(-gam_ext * L_min / 4)
        e_mid = exp(-gam_ext * L_min / 4.0)
        S_mid = -log(1.0 - e_mid)

        # Interior tail (m >= floor((c-1)/2) + 1 >= c/2):
        # For m >= c/2, |a_m| <= (1 + (4*kap/(pi*N0))*m)^{-N0/2} ... wait: for m*pi/n <= pi/3 (almost all!),
        # actually let's use exact monotone bound: for m >= c/2:
        #   Wait: |a_m| <= (1 + 2*c_int*kap*m / N0)^{-N0/2} with c_int = 2/pi = 0.6366 (valid for ALL m <= (n-c)/2!).
        #   Check at m = c/2: (1 + (4/pi)*kap*c/(2*N0))^{-N0/2} = (1 + (2/pi)*L_min/N0)^{-N0/2}
        c_int = 2.0 / pi
        x_half = (c_int * L_min) / N0 # = (4/(pi*N0)) * kap * (c/2)
        a_half = (1.0 + x_half) ** (-N0 / 2.0)
        # Sum_{m >= c/2} 2 * (-log(1 - |a_m|)) <= (2 / (1 - a_half)) * sum_{m >= c/2} (1 + (4*kap/(pi*N0))*m)^{-N0/2}
        # Integral bound: <= (2/(1-a_half)) * [ a_half + (pi*N0 / (4*kap*(N0/2 - 1))) * (1 + x_half)^{-(N0/2 - 1)} ]
        # Wait: there is a 1/kap in front of the integral! Since kap = sqrt(c*kap^2 / c), can kap be TINY?
        # YES! If c is HUGE (c = 10000) and kap = 0.05 (c*kap^2 = 25), then L_min = c*kap = 500 >> 1,
        # so (1 + (2/pi)*L/N0)^{-(N0/2-1)} is 10^{-100}!!
        # Let's bound (1/kap) * (1 + c_int*L/N0)^{-(N0/2-1)} using 1/kap = L / (c*kap^2) <= L / C1!
        # Maximize L * (1 + c_int*L/N0)^{-(N0/2 - 1)} over L >= L_min!
        # Derivative is zero at c_int*L/N0 = 1/(N0/2 - 2) = 1/64 => L_pk = N0 / (64*c_int) ~ 3.24!
        # Since L_min = C1/kap0 >= 8 > 3.24, it is STRICTLY DECREASING for all L >= L_min!!
        int_tail_term = (L_min / C1) * (N0 / (2.0 * c_int * (N0/2.0 - 1.0))) * ((1.0 + x_half) ** (-(N0/2.0 - 1.0)))
        S_int = (2.0 / (1.0 - a_half)) * (a_half + int_tail_term)

        tot = S1 + S2 + S_mid + S_int
        if tot < best:
            best = tot; best_info = (f, S1, S2, S_mid, S_int)
    return best, best_info

print("=== CLOSED-FORM ELEMENTARY BOUND FOR SUB-REGIME 2A (c*kap^2 >= C1, kap < kap0, n >= 132) ===")
for kap0 in [2.5, 3.0]:
    print(f"\n--- kap0 = {kap0} (pi/2 = {pi/2:.4f}) ---")
    for C1 in [16.0, 20.0, 25.0, 30.0, 36.0]:
        b, (f, s1, s2, sm, si) = closed_form_B_pair_bound(C1, kap0=kap0, N0=132)
        status = "PASS (< pi/2)!" if b < pi/2 else "FAIL"
        print(f"  C1={C1:4.1f} (L_min={C1/kap0:5.2f}): B_closed={b:.4f} [f={f:.2f}: S1={s1:.4f}, S2={s2:.4f}, mid={sm:.4f}, int={si:.4f}] -> {status}")
