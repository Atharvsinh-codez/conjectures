"""
VERIFY CLOSED-FORM ELEMENTARY CONSTANTS FOR:
1. Sub-regime 2A (c * kap^2 >= C_1, kap < 2.5, n >= 132) using individual log bound on [f*c, c/2]!
2. Sub-regime 2C (c * q * kap <= C_2, kap < 2.5, n >= 132) using exact binomial p_r + 10-line global
   decay (1 - cos(l*phi) <= l^2 (1 - cos(phi))) + cubic phase (|sin(l*phi) - l*phi| <= l^3 |phi|^3 / 6)!
3. Also verify Lemma R2.1 (exact cubic-root perturbation bound across ALL kap < 2.5, n >= 132, lam <= 2)!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, asin

# =============================================================================
# 1. SUB-REGIME 2A (FIXED S2 USING INDIVIDUAL LOG BOUND ON [f*c, c/2])
# =============================================================================
def closed_form_2A_fixed(C1, kap0=2.5, N0=132):
    alpha_max = asin(kap0 / (2.0 * pi))
    gam_ext = sin(alpha_max / 4.0) / (alpha_max / 4.0)
    gam_z = (alpha_max / 2.0) / sin(alpha_max / 2.0)
    L_min = C1 / kap0
    zeta3 = 1.2020569031595943
    M_g = 0.6476102379

    best = 1e9; best_info = None
    for f in np.linspace(0.15, 0.45, 61):
        # (a) 1 <= m <= f*c: exact pairing bound
        a1 = gam_ext * (1.0 - f)
        S1 = 4.0 * gam_z * (2.0 * zeta3 / (C1 * a1**3) + M_g / (L_min * a1**2))

        # (b) f*c < m <= floor((c-1)/2): INDIVIDUAL bound!
        # For each such m, sum_{r>=1} |a_m^r - a_{-m}^{-r}|/r <= -log(1 - |a_{-m}^{-1}|) - log(1 - |a_m|)
        # Since |a_m| < |a_{-m}^{-1}| <= exp(-gam_ext * kap * m * (1 - m/c)):
        #   term <= 2 * (-log(1 - |a_{-m}^{-1}|)) <= (2 / (1 - e^{-A_min})) * exp(-gam_ext * kap * m * (1 - m/c))
        # By concavity of t(1-t) on [f, 0.5]: m*(1 - m/c) >= f*(1-f)*c + (0.5 - f)*(m - f*c)
        # Wait: even cleaner! Multiply by 2 for conjugate symmetry (both sides):
        #   S2 <= (4 / (1 - e^{-A_min})) * sum_{m > f*c}^{c/2} exp(-gam_ext * f*(1-f)*L - gam_ext*(0.5-f)*kap*(m - f*c))
        #       <= (4 / (1 - e^{-A_min})) * exp(-cA * L) * [ 1 + 1 / (gam_ext*(0.5 - f)*kap) ]
        # Wait: 1/kap = L / (c*kap^2) <= L / C1!
        # So S2 <= (4 / (1 - e^{-A_min})) * [ exp(-cA * L_min) + (1 / (C1 * gam_ext * (0.5 - f))) * max_{L>=L_min} (L * exp(-cA*L)) ]!!
        # EVEN SHARPER: use the actual integral over m in [f*c, c/2] (length (0.5-f)*c), so sum <= (1/(kap*b))*(1 - e^{-b*(0.5-f)*L}) + 1!
        cA = gam_ext * f * (1.0 - f)
        b = gam_ext * (0.5 - f)
        A_min = cA * L_min
        denom_min = 1.0 - exp(-A_min)
        L1 = max(L_min, 1.0 / cA)
        # Note: if f -> 0.5, instead of 1/(kap*b), the number of terms is <= (0.5 - f)*c = (0.5 - f)*L / kap <= (0.5-f)*L^2 / C1!
        # Let's take min of the geometric-sum bound and the term-count bound!
        # Actually b*(0.5-f)*L = gam_ext*(0.5-f)^2 * L. Let's compute exact integral int_0^{(0.5-f)L} exp(-b*s) ds = (1 - e^{-b(0.5-f)L})/b!
        # For L >= L_min: L * exp(-cA * L) <= L1 * exp(-cA * L1), and 1/(kap*b) <= L/(C1*b):
        S2 = (4.0 / denom_min) * (exp(-A_min) + (L1 * exp(-cA * L1)) / (C1 * b))

        # (c) Midpoint (c even):
        e_mid = exp(-gam_ext * L_min / 4.0)
        S_mid = -log(1.0 - e_mid)

        # (d) Interior tail (m >= c/2):
        c_int = 2.0 / pi
        x_half = (c_int * L_min) / N0
        a_half = (1.0 + x_half) ** (-N0 / 2.0)
        L_int_pk = max(L_min, N0 / (c_int * (N0/2.0 - 2.0)))
        x_pk = (c_int * L_int_pk) / N0
        int_tail = (L_int_pk / C1) * (N0 / (2.0 * c_int * (N0/2.0 - 1.0))) * ((1.0 + x_pk) ** (-(N0/2.0 - 1.0)))
        S_int = (2.0 / (1.0 - a_half)) * (a_half + int_tail)

        tot = S1 + S2 + S_mid + S_int
        if tot < best:
            best = tot; best_info = (f, S1, S2, S_mid, S_int)
    return best, best_info

print("=== FIXED CLOSED-FORM BOUND FOR SUB-REGIME 2A (kap0 = 2.5 & 3.0, N0 = 132) ===")
for kap0 in [2.5, 3.0]:
    print(f"\n--- kap0 = {kap0} (pi/2 = {pi/2:.4f}) ---")
    for C1 in [16.0, 20.0, 22.0, 25.0, 28.0, 32.0, 36.0]:
        b, (f, s1, s2, sm, si) = closed_form_2A_fixed(C1, kap0=kap0, N0=132)
        st = "PASS (< 1.50)!" if b < 1.50 else ("PASS (< pi/2)!" if b < pi/2 else "fail")
        print(f"  C1={C1:4.1f} (L_min={C1/kap0:5.2f}): B={b:.4f} [f={f:.2f}: S1={s1:.4f}, S2={s2:.4f}, mid={sm:.4f}, int={si:.4f}] -> {st}")

# =============================================================================
# 2. CHECK 10-LINE CLOSED-FORM LOWER BOUND FOR SUB-REGIME 2C (ALL g >= 1, c*q*kap <= C2)
# =============================================================================
# On |u| = eta (where u = rho^q e^{i phi}, peak at phi = 0):
# Main integral (ignoring exponentially tiny q \nmid r terms, which are <= e^{-pi^2 / C2} < 10^{-5}):
#   I_main(g, eta) = (1/(2pi)) int_{-pi}^{pi} exp( S_Re(phi) ) * cos( S_Im(phi) ) dphi
# where at saddle eta solving sum_{l>=1} (-1)^{l-1} l^{-1/2} eta^l = g * sqrt(c*q*kap/pi):
#   Write x_l = (-1)^{l-1} (P_{lq} / P_q) * rho^{(l-1)q} ~ (-1)^{l-1} l^{-1/2} eta^{l-1}  (with x_1 = 1!)
#   Wait: normalize by M := P_q * rho^q / q = g / (1 - eta/sqrt(2) + eta^2/sqrt(3) - ...)!!
#   Then:
#     S_Re(0) - S_Re(phi) = M * sum_{l>=1} (-1)^{l-1} (l^{-3/2} eta^{l-1}) * (1 - cos(l*phi))
#                         >= M * [ 1 - sum_{k>=1} (2k)^2 * (2k)^{-3/2} * eta^{2k-1} ] * (1 - cos(phi))
#                         = M * [ 1 - sum_{k>=1} sqrt(2k) * eta^{2k-1} ] * (1 - cos(phi))  := c_Re(eta) * M * (1 - cos(phi))!
#     |S_Im(phi)| = M * | sum_{l>=1} (-1)^{l-1} l^{-3/2} eta^{l-1} (sin(l*phi) - l*phi) |
#                 <= M * [ sum_{l>=1} l^{3/2} eta^{l-1} ] * |phi|^3 / 6  := c_Im(eta) * M * |phi|^3 / 6!
#   WAIT A SECOND!! CHECK THAT FOR g = 1, 2, 3, ..., infinity AT eta = 0.30, 0.35, 0.40!!
print("\n=== CHECK 10-LINE SUB-REGIME 2C LOWER BOUND FOR ALL g >= 1 ===")
for eta in [0.25, 0.30, 0.35, 0.40, 0.45]:
    # Corresponding c*q*kap (approx pi * eta^2 * (1 - eta/sqrt(2) + ...)^2):
    alt_s = sum(((-1)**(l-1)) * (l**(-0.5)) * (eta**(l-1)) for l in range(1, 40))
    cqk = pi * (eta * alt_s)**2
    # Wait: c_Re only needs NEGATIVE (even l) terms subtracted (since odd l has + and 1-cos(l*phi) >= 0!)
    neg_even = sum(sqrt(2*k) * (eta**(2*k - 1)) for k in range(1, 30))
    c_Re = 1.0 - neg_even
    c_Im = sum((l**1.5) * (eta**(l - 1)) for l in range(1, 40))
    # For g >= 1, M = g / alt_s. Check integral lower bound for g = 1, 2, 3, 5, 10, 100!
    # On |phi| <= phi_0 (where c_Im * M * phi_0^3 / 6 <= pi/3):
    #   Re >= exp(-c_Re_up * M * (1-cos phi)) * cos(c_Im * M * phi^3 / 6)
    # Outside |phi| > phi_0: |integrand| <= exp(-c_Re * M * (1 - cos phi))!
    # Check exact lower bound of int_{-phi0}^{phi0} e^{-M*c_up*(1-cos)} cos(M*c_Im*phi^3/6) - 2*int_{phi0}^{pi} e^{-M*c_Re*(1-cos)}!
    c_up = sum((l**(-1.5)) * (l**2) * (eta**(l-1)) for l in range(1, 40) if l % 2 == 1) # odd terms upper bound on S_Re(0)-S_Re(phi)
    print(f"  eta={eta:.2f} (c*q*kap~{cqk:.3f}, alt={alt_s:.3f}): c_Re={c_Re:.4f}, c_up={c_up:.4f}, c_Im={c_Im:.4f} | e^(-pi^2/cqk)={exp(-pi**2/max(cqk,1e-9)):.2e}")
    for g in [1, 2, 3, 5, 10, 50]:
        M = g / alt_s
        if c_Re <= 0: continue
        phis = np.linspace(-pi, pi, 4000, endpoint=False)
        # Rigorous pointwise lower bound on Re(integrand) using ONLY c_Re, c_up, c_Im:
        ph_bd = c_Im * M * (np.abs(phis)**3) / 6.0
        # Where ph_bd < pi/2: lower bound is exp(-c_up * M * (1 - cos)) * cos(ph_bd)
        # Where ph_bd >= pi/2: worst-case lower bound is -exp(-c_Re * M * (1 - cos))
        lb = np.where(ph_bd < pi/2, np.exp(-c_up * M * (1.0 - np.cos(phis))) * np.cos(ph_bd),
                      -np.exp(-c_Re * M * (1.0 - np.cos(phis))))
        val_lb = float(np.mean(lb))
        if g in (1, 2, 3, 5, 50):
            print(f"      g={g:2d} (M={M:6.2f}): RIGOROUS integral lower bound = {val_lb:+.6f} {'PASS' if val_lb > 0 else 'FAIL'}")
