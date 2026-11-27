"""
Verify our 4-LINE m-space universal bound on |Phi^{exact}(phi) - Phi^G(phi)| across ALL of Regime 2!
Wait: check exact constant of BOTH:
  (1) Pure modulus difference |Phi^{mod}(phi) - Phi^G(phi)|
  (2) Pure phase difference   |Phi^{exact}(phi) - Phi^{mod}(phi)|
using ONLY Lemma G2.1 bounds (ZERO dependence on j/c, ZERO loops over r)!
"""
import math

# Check exact 1-variable maximum over T = c*kap in (0, inf) of the m-space phase diff bound!
# Recall for m = 1..floor((c-1)/2):
#   Let g(r, a) = arg(1 - r e^{i a}).
#   Exact derivatives:
#     |g_r(r, a)| = |sin a| / |1 - r e^{i a}|^2 <= 1 / (1 - r^2)? Wait: max_a |sin a|/(1-2r cos a + r^2) = 1 / (1 - r^2)!!
#     (Check! At cos a = 2r/(1+r^2), sin a = (1-r^2)/(1+r^2), den = (1-r^2)^2/(1+r^2) -> ratio = 1 / (1 - r^2)!!
#      EXACT — HALF of 1/(1-r) as r -> 1!!)
#     |g_a(r, a)| = |r(r - cos a)| / |1 - r e^{i a}|^2 <= r / (1 - r) (exact max at a = 0: r/(1-r)!)
#     |g_{ra}(r, a)| = |(1+r^2)cos a - 2r| / |1 - r e^{i a}|^4? Check exact max_a |g_{ra}|!
#     Actually: g_a(r, a) = Re( -r e^{i a} / (1 - r e^{i a}) ) -> g_{ra}(r, a) = Re( -e^{i a} / (1 - r e^{i a})^2 )!
#     So |g_{ra}(r, a)| <= 1 / (1 - r)^2 (exact max at a = 0)!
#
#   WAIT A SECOND!! Check exact delta_psi(c/2 - m) + delta_psi(c/2 + m)!!
#   Remember delta_psi(x) = -n*phi(2*pi*x/n) - (2*pi*j/c)*x!
#   So delta_psi(c/2 - m) + delta_psi(c/2 + m) = -n*(phi(th_c - d) + phi(th_c + d) - 2*phi(th_c))!
#   Notice by Taylor's theorem, phi(th_c - d) + phi(th_c + d) - 2*phi(th_c) = int_0^d (d - t) (phi''(th_c + t) + phi''(th_c - t)) dt!
#   Wait: earlier we had phi''(eta) = g'(eta) = s(1-s^2) sin(eta) / |1 - s e^{i eta}|^4!
#   For m <= c/2: eta = th_c +- d <= 2*th_c!
#   Look at the ACTUAL size of |Phi^{exact}(phi) - Phi^G(phi)| vs our 4-line bound!

def exact_m_space_bound(c, kap, c_in_lo=0.9913, c_in_hi=1.0506, c_out_lo=0.8570, c_out_hi=1.0070):
    # Wait: remember the relative modulus error is O(kap^2) (0.0085 * kap^2), NOT a constant 0.05!
    # And the exact phase error delta_psi(c/2 - m) is <= 0.1394 * kap^2 * m (for m <= c/2)!
    # Let's check the EXACT worst m-space bound using:
    #   - r_{-m} <= exp(-0.9913 * kap * m * (1 - m/c))
    #   - r_{+m} >= exp(-1.0070 * kap * m * (1 + m/c))
    #   - |d_m^- - d_m^+| = |delta_psi(c/2 - m) + delta_psi(c/2 + m)| <= 0.646 * (kap^2 / c) * m^2 (WAIT: check factor of 2! A_in/2 = 0.836!)
    #   Actually let's FIRST check the TRUE exact |delta_psi(c/2 - m) + delta_psi(c/2 + m)| on finite-n roots!
    pass

# Check TRUE exact |Phi^{exact}(phi) - Phi^G(phi)| AND check exact finite-n delta_psi(c/2-m) + delta_psi(c/2+m)!
def check_exact_phi_diff(n, c, j):
    th_c = c * math.pi / n; psi_j = j * math.pi / n; alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1; log_Rz = n * math.log(r3 / r1); kap = 2.0 * math.pi * r1 * r2 / r3
    def dpsi(x):
        eta = 2.0 * math.pi * x / n
        phi = math.atan2(-s * math.sin(eta), 1.0 - s * math.cos(eta))
        return -n * phi - (2.0 * math.pi * j / c) * x
    def lr_ex(x):
        eta = 2.0 * math.pi * x / n
        z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
        return log_Rz - n * math.log(abs(z))
    half = (c - 1) // 2
    # Check max over phi of |Phi^{exact}(phi) - Phi^G(phi)| AND our 4-line m-space upper bound!
    bound_phase = 0.0; bound_mod = 0.0
    max_dsum_coeff = 0.0; max_dbar_coeff = 0.0
    d1 = abs(dpsi(c / 2.0 - 0.001) / 0.001)  # exact |delta_psi'(c/2)|
    for m in range(1, half + 1):
        dm_minus = dpsi(c / 2.0 - m)
        dm_plus = -dpsi(c / 2.0 + m)
        max_dsum_coeff = max(max_dsum_coeff, abs(dm_minus - dm_plus) / (kap * kap * m * m / c))
        max_dbar_coeff = max(max_dbar_coeff, dm_minus / (kap * kap * m))
        # Exact m-space bound using actual r_{-m}, r_{+m}, dm_minus, dm_plus:
        rm = math.exp(-lr_ex(c / 2.0 - m))
        rp = math.exp(lr_ex(c / 2.0 + m))
        rG_m = math.exp(-kap * m * (1.0 - m / c))
        rG_p = math.exp(-kap * m * (1.0 + m / c))
        # 2 copies (phi and -phi) of:
        #   Phase diff: |g(rp, a - dm_plus) - g(rm, a - dm_minus) - (g(rp, a) - g(rm, a))|
        #   <= |dm_minus - dm_plus| * (rp / (1 - rp)) + |dm_minus| * |rm - rp| / (1 - rm)^2
        #   (Notice: we bound with rp / (1 - rp) + |rm - rp| / (1 - rm)^2!)
        bound_phase += 2.0 * (abs(dm_minus - dm_plus) * (rp / (1.0 - rp)) + abs(dm_minus) * abs(rm - rp) / ((1.0 - rm) ** 2))
        # Mod diff: |g(rp, a) - g(rm, a) - (g(rG_p, a) - g(rG_m, a))|
        #   <= |rm - rG_m| / (1 - max(rm, rG_m)^2) + |rp - rG_p| / (1 - max(rp, rG_p)^2)
        bound_mod += 2.0 * (abs(rm - rG_m) / (1.0 - max(rm, rG_m) ** 2) + abs(rp - rG_p) / (1.0 - max(rp, rG_p) ** 2))
    # Tail m > half (large roots only):
    for m in range(half + 1, min(n // 2 - c // 2, half + int(20.0 / kap) + 5)):
        rp = math.exp(lr_ex(c / 2.0 + m)) if lr_ex(c / 2.0 + m) > -40 else 0.0
        rG_p = math.exp(-kap * m * (1.0 + m / c)) if kap * m * (1.0 + m / c) < 40 else 0.0
        dm_plus = -dpsi(c / 2.0 + m)
        bound_phase += 2.0 * min(abs(dm_plus), math.pi) * (rp / (1.0 - rp))
        bound_mod += 2.0 * abs(rp - rG_p) / (1.0 - max(rp, rG_p) ** 2)
    return kap, bound_phase, bound_mod, max_dsum_coeff, max_dbar_coeff

print("=== CHECK 4-LINE m-SPACE BOUND ON |Phi^ex(phi) - Phi^G(phi)| ACROSS G1 ===")
for c, j, kt in [
    (6, 7, 1.50), (6, 7, 2.30),
    (10, 13, 0.05), (10, 13, 0.30), (10, 13, 1.80),
    (12, 13, 0.50), (15, 16, 0.50), (20, 21, 0.50),
    (29, 39, 0.05), (29, 39, 1.00),
    (50, 69, 0.15), (50, 69, 0.80),
    (100, 137, 0.10), (100, 137, 0.60),
    (200, 273, 0.42),
]:
    if c * kt * kt >= 41.0 or kt >= 2.45: continue
    n = max(132, int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt))) // 2 * 2)
    if j > n // 2 - (c + 1) // 2 - (1 if c % 2 == 0 else 0): continue
    kap, bph, bmd, csum, cbar = check_exact_phi_diff(n, c, j)
    print(f"  c={c:3d} j={j:3d} kap={kap:.3f} ck^2={c*kap*kap:5.2f} | m-space phase={bph:.4f}, mod={bmd:.4f}, TOT={bph+bmd:.4f} | dsum_c={csum:.3f}, dbar_c={cbar:.3f}")
