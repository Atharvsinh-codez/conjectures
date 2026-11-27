"""
Independent 100% RIGOROUS OUTWARD-ROUNDED INTERVAL ARITHMETIC (mpmath.iv, 50-digit precision)
audit of Regime 1 (ALGEBRA_ROUTE.md Section 4, ALL even n >= 84, kappa >= 3.0 and kappa >= 2.5).

Verifies from scratch:
1. All 6 closed-form margins using rigorous interval arithmetic (mp.iv) so ZERO floating-point
   rounding error can hide anywhere!
2. Explicit verification of every side hypothesis (strict Rouché inequalities, eta_pert < 1 for
   small-root tail monotonicity, R_y > 1/2 subcase with corrected rigorous bound r3/r2 < 0.4814).
3. Independent 80-digit mpmath verification of Theorem 3.2 (4-factor identity) and Lemma 4.5
   (exact Cauchy integral B_avg) against exact integer engine.py coefficients!
"""

from __future__ import annotations
import sys
import mpmath as mp
from mpmath import iv

sys.path.insert(0, r"d:\solving-conjectures\work\bardakov-kozlovskaya\verify")
import engine as bk


def iv_cot(x):
    return iv.cos(x) / iv.sin(x)


def S_off_iv(k_min, c_min, theta_max, N0=84):
    k_min = iv.mpf(k_min)
    theta_max = iv.mpf(theta_max)
    N0_iv = iv.mpf(N0)

    # Large factors upper bound
    s_lg = iv.mpf(0)
    for r in (1, 2):
        r_iv = iv.mpf(r)
        fac = iv.mpf(1) - ((r_iv * iv.pi / N0_iv) ** 2) / 6
        s_lg += 2 * ((1 + (2 * k_min * r_iv * fac) / N0_iv) ** (-N0_iv / 2))
    c_J = 4 * k_min / iv.pi
    for r in range(3, 21):
        r_iv = iv.mpf(r)
        s_lg += 2 * ((1 + c_J * r_iv / N0_iv) ** (-N0_iv / 2))
    s_lg += 2 * (N0_iv / (c_J * (N0_iv / 2 - 1))) * ((1 + 20 * c_J / N0_iv) ** (-(N0_iv / 2 - 1)))

    # Small factors upper bound
    half_th = theta_max / 2
    lam = (iv.sin(half_th) / half_th) * k_min
    cands_hi = []
    if c_min <= 2:
        cands_hi.append(iv.exp(-lam * iv.mpf('0.5')).b)
    if c_min <= 3:
        cands_hi.append((2 * iv.exp(-lam * (iv.mpf(2) / 3))).b)
    if c_min <= 4:
        cands_hi.append((2 * iv.exp(-lam * iv.mpf('0.75')) + iv.exp(-lam)).b)
    if c_min <= 5:
        cands_hi.append((2 * iv.exp(-lam * iv.mpf('0.8')) + 2 * iv.exp(-lam * iv.mpf('1.2'))).b)
    cm = iv.mpf(max(c_min, 6))
    r1_b = 2 * iv.exp(-lam * (1 - 1 / cm))
    r2_b = 2 * iv.exp(-lam * 2 * (1 - 2 / cm))
    r3_b = 2 * iv.exp(-lam * iv.mpf('1.5')) / (1 - iv.exp(-lam * iv.mpf('0.5')))
    cands_hi.append((r1_b + r2_b + r3_b).b)
    s_sm_ub = iv.mpf(max(cands_hi))
    s_lg_ub = iv.mpf(s_lg.b)
    return s_sm_ub, s_lg_ub


def solve_rouche_iv(p_mod_ub, Rz_Ry_ub, Ry_omt_ub, omt_lb):
    """Find explicit rational/decimal rho_x, u_y and CERTIFY STRICT Rouché inequalities via iv!"""
    # Use mp.mpf to find a tight candidate, nudge up by 1e-8, then CERTIFY with iv!
    px = float(p_mod_ub.b) / float(omt_lb.a)
    for _ in range(100):
        px = (1 + 1e-9) * float(p_mod_ub.b) / ((1 - px) * (float(omt_lb.a) - px))
    rho_x = iv.mpf(str(px * (1 + 1e-7)))
    lhs_x = rho_x * (1 - rho_x) * (omt_lb - rho_x)
    assert rho_x.b < min(0.5, 0.5 * float(omt_lb.a))
    assert lhs_x.a > p_mod_ub.b, f"Strict Rouché x FAILED: {lhs_x.a} <= {p_mod_ub.b}"

    uy = float(Rz_Ry_ub.b) / float(omt_lb.a)
    for _ in range(100):
        uy = (1 + 1e-9) * float(Rz_Ry_ub.b) / ((1 - uy) * (float(omt_lb.a) - uy * float(Ry_omt_ub.b)))
    u_y = iv.mpf(str(uy * (1 + 1e-7)))
    lhs_y = u_y * (1 - u_y) * (omt_lb - u_y * Ry_omt_ub)
    assert u_y.b < 0.5 and (u_y * Ry_omt_ub).b < 0.5 * float(omt_lb.a)
    assert lhs_y.a > Rz_Ry_ub.b, f"Strict Rouché y FAILED: {lhs_y.a} <= {Rz_Ry_ub.b}"

    dy = iv.mpf(((u_y + rho_x) / (1 - rho_x)).b)
    dz = iv.mpf((1 / (((1 - rho_x) ** 3) * (1 - dy)) - 1).b)
    return rho_x, dy, dz


def eval_closed_iv(label, k_min, c_min, theta_max, dx, dy, dz, t1_ub, LmC_Ldy_ub, eta_pert_ub, t2_quad_ub, Bmin_lb, N0=84):
    N0_iv = iv.mpf(N0)
    # Verify side-condition for small-root monotonicity: 1/N0 + log(0.5 + eta_pert/N0) < 0
    assert (1 / N0_iv + iv.log(iv.mpf('0.5') + eta_pert_ub / N0_iv)).b < 0, f"Small-root monotonicity failed for {label}!"
    assert eta_pert_ub.b < N0

    s_sm, s_lg = S_off_iv(k_min, c_min, theta_max, N0=N0)
    inv_1pdz = ((1 + dx) ** 3) * (1 + dy)
    s_lg_t = s_lg * (1 + dz) * ((1 - eta_pert_ub / N0_iv) ** (-N0_iv))
    s_sm_t = (s_sm * iv.exp(2 * eta_pert_ub) + (N0_iv / 3) * ((iv.mpf('0.5') + eta_pert_ub / N0_iv) ** N0_iv)) * inv_1pdz
    s_off_t = s_sm_t + s_lg_t

    t2 = LmC_Ldy_ub + t2_quad_ub
    M1 = iv.exp(t1_ub + t2 + s_off_t) - 1
    ep_bdy = iv.exp(eta_pert_ub) - 1
    R_bdy = dz + 2 * ep_bdy + (((1 + ep_bdy) ** 2) * inv_1pdz - 1)
    margin = Bmin_lb - M1 * (2 + R_bdy)
    print(f"  {label:22s}: eta_pert<={float(eta_pert_ub.b):.4f}, M1<={float(M1.b):.4f}, Rbdy<={float(R_bdy.b):.4f}, Bmin>={float(Bmin_lb.a):.4f} ==> RIGOROUS MARGIN_LB = {float(margin.a):+.6f}")
    assert margin.a > 0, f"RIGOROUS MARGIN FAILED FOR {label}: {margin}"
    return float(margin.a)


def audit_regime1_iv(k0_str="3.0", N0=84):
    iv.dps = 50
    k0 = iv.mpf(k0_str)
    N0_iv = iv.mpf(N0)
    print(f"\n=== RIGOROUS OUTWARD-ROUNDED INTERVAL CERTIFICATION OF REGIME 1 (n >= {N0}, kappa >= {k0_str}) ===")

    # Case 1A (Subcase R_y > 1/2 check first!)
    # Rigorous bound: cot(psi_j) <= (2^{1/N0}-1)/sin(2*pi/N0) + tan(0.075*pi)
    cot_max = (iv.mpf(2) ** (1 / N0_iv) - 1) / iv.sin(2 * iv.pi / N0_iv) + iv.tan(iv.mpf('0.075') * iv.pi)
    sin_psi_min = 1 / iv.sqrt(1 + cot_max ** 2)
    r3_r2_max = iv.sin(iv.mpf('0.15') * iv.pi) / sin_psi_min
    Rz_Ry_Rygt = r3_r2_max ** N0_iv
    print(f"  [1A R_y>1/2 fix] rigorous r3/r2 <= {float(r3_r2_max.b):.5f} => Rz/Ry <= {float(Rz_Ry_Rygt.b):.3e} (<< 4/n^2 * 10^-20!)")
    assert (N0_iv * N0_iv * Rz_Ry_Rygt).b < 1e-22

    # Case 1A-i (R_y <= 1/2)
    Rz_Ry_1A = iv.exp(-iv.pi * iv_cot(iv.mpf('0.15') * iv.pi + iv.pi / (2 * N0_iv)))
    sinc_15 = iv.sin(iv.mpf('0.15') * iv.pi) / (iv.mpf('0.15') * iv.pi)
    c_j = iv.sqrt(sinc_15 * k0 / (2 * iv.pi * iv.pi))
    cot_247 = iv.sqrt(1 - iv.mpf('0.70') ** 2) / iv.mpf('0.70')
    Rz_1A_a = (iv.sin(iv.mpf('0.15') * iv.pi) / iv.mpf('0.70')) ** N0_iv
    Rz_1A_b = iv.exp(-iv.pi * cot_247 * c_j * iv.sqrt(N0_iv))
    Rz_1A = iv.mpf(max(Rz_1A_a.b, Rz_1A_b.b))
    nRz_1A = N0_iv * Rz_1A
    dx_1A, dy_1A, dz_1A = solve_rouche_iv(iv.mpf('0.5') * Rz_1A, Rz_Ry_1A, iv.mpf('0.5'), iv.mpf('0.5'))
    L_dy_1A = -iv.log(1 - dy_1A)
    c_dx_1A = iv.mpf('0.25') / ((1 - dx_1A) * (iv.mpf('0.5') - dx_1A))
    c_uy_1A = 1 / ((1 - dy_1A) * (iv.mpf('0.5') - iv.mpf('0.5') * dy_1A))
    C_Ldy_1A = (L_dy_1A / dy_1A) * ((c_uy_1A + c_dx_1A) / (1 - dx_1A))
    r2_r3_Ldy_1A = C_Ldy_1A * (Rz_Ry_1A ** ((N0_iv - 1) / N0_iv))
    k_Ldy_cands = [
        (iv.mpf('1.5') * Rz_Ry_1A).b,
        (iv.mpf('2.5') * iv.exp(-2 * iv.pi * iv_cot(iv.mpf('0.20') * iv.pi))).b,
        (iv.mpf('0.5') * N0_iv * ((iv.sin(iv.mpf('0.15') * iv.pi) / iv.sin(iv.mpf('0.25') * iv.pi)) ** N0_iv)).b,
    ]
    LmC_Ldy_1A = C_Ldy_1A * iv.mpf(max(k_Ldy_cands)) + 2 * r2_r3_Ldy_1A
    t1_1A = (iv.mpf('0.5') / ((1 - dx_1A) * (iv.mpf('0.5') - dx_1A))) * nRz_1A
    t2_quad_1A = (r2_r3_Ldy_1A * L_dy_1A) / (2 * (1 - L_dy_1A / N0_iv)) + (r2_r3_Ldy_1A ** 2) / (2 * (1 - r2_r3_Ldy_1A / N0_iv))
    Bmin_1A = 2 - (5 * (Rz_Ry_1A ** ((N0_iv - 1) / N0_iv)) + 2 * Rz_Ry_1A)
    m1A = eval_closed_iv("Case 1A (c<=0.15n)", k0, 2, iv.mpf('0.15') * iv.pi, dx_1A, dy_1A, dz_1A,
                         t1_1A, LmC_Ldy_1A, r2_r3_Ldy_1A, t2_quad_1A, Bmin_1A, N0=N0)

    # Helper for sa1, sa2
    def sa_iv(th):
        return iv.log(iv_cot(th / 4)) / iv.pi, iv_cot(th / 2) / (2 * iv.pi)

    # Case 1B-I
    k_1B = 2 * iv.pi * iv.sin(iv.mpf('0.30') * iv.pi)
    Rz_Ry_1BI = iv.exp(-2 * iv.pi / iv.sqrt(3))
    Ry_1BI_even = iv.exp(-2 * iv.pi * iv_cot(iv.mpf('0.425') * iv.pi))
    omt_1BI = 1 - Ry_1BI_even
    dx_1BI, dy_1BI, dz_1BI = solve_rouche_iv(Ry_1BI_even * Rz_Ry_1BI, Rz_Ry_1BI, Ry_1BI_even, omt_1BI)
    L_dy_1BI = -iv.log(1 - dy_1BI)
    c_1BI_exp = iv.pi / (2 * iv.sqrt(3))
    t1_1BI_o = (5 * iv.exp(-5 * c_1BI_exp)) / ((1 - dx_1BI) ** 2)
    t1_1BI_e = (6 * iv.exp(-6 * c_1BI_exp)) / ((1 - dx_1BI) * (omt_1BI - dx_1BI))
    t1_1BI = iv.mpf(max(t1_1BI_o.b, t1_1BI_e.b))
    r2r3_1BI = 1 / iv.sin(iv.mpf('0.15') * iv.pi)
    sa1_1B, sa2_1B = sa_iv(iv.mpf('0.15') * iv.pi)
    g_abs_1BI = r2r3_1BI * iv.cos(iv.mpf('0.30') * iv.pi)
    C_Ldy_1BI = (L_dy_1BI / dy_1BI) * ((1 / ((1 - dy_1BI) * (omt_1BI - dy_1BI * Ry_1BI_even)) + Ry_1BI_even / ((1 - dx_1BI) * (omt_1BI - dx_1BI))) / (1 - dx_1BI))
    LmC_Ldy_1BI = iv.mpf(max(
        ((iv.mpf('2.5') + 2 * g_abs_1BI) * L_dy_1BI).b,
        ((iv.mpf('3.5') + 2 * g_abs_1BI) * C_Ldy_1BI * iv.exp(-3 * iv.pi / iv.sqrt(3))).b,
    ))
    r2_r3_Ldy_1BI = r2r3_1BI * L_dy_1BI
    t2_quad_1BI = (sa1_1B * (L_dy_1BI ** 2) / (2 * (1 - L_dy_1BI / N0_iv)) + sa2_1B * (L_dy_1BI ** 2) / (2 * (1 - r2_r3_Ldy_1BI / N0_iv)))
    Bmin_1BI = 2 - (2 * g_abs_1BI * (1 + Ry_1BI_even ** 2) + 1) * (Rz_Ry_1BI / omt_1BI)
    m1BI = eval_closed_iv("Case 1B-I (j>=c+2)", k_1B, 13, iv.pi / 3, dx_1BI, dy_1BI, dz_1BI,
                          t1_1BI, LmC_Ldy_1BI, r2_r3_Ldy_1BI, t2_quad_1BI, Bmin_1BI, N0=N0)

    # Case 1B-IIa
    Rz_Ry_IIa = iv.exp(-iv.pi * iv_cot(iv.mpf('0.25') * iv.pi + iv.pi / (2 * N0_iv)))
    Ry_IIa = iv.exp(-iv.mpf('0.15') * N0_iv * iv.pi * iv_cot(iv.mpf('0.425') * iv.pi))
    omt_IIa = 1 - Ry_IIa
    dx_IIa, dy_IIa, dz_IIa = solve_rouche_iv(Ry_IIa * Ry_IIa * Rz_Ry_IIa, Rz_Ry_IIa, Ry_IIa, omt_IIa)
    L_dy_IIa = -iv.log(1 - dy_IIa)
    t1_IIa = (iv.mpf('0.55') * N0_iv * Ry_IIa * Ry_IIa * Rz_Ry_IIa) / ((1 - dx_IIa) * (omt_IIa - dx_IIa))
    r2r3_IIa = iv.cos(iv.pi / N0_iv) + iv_cot(iv.mpf('0.15') * iv.pi) * iv.sin(iv.pi / N0_iv)
    g_abs_IIa = r2r3_IIa * iv.cos(iv.mpf('0.30') * iv.pi)
    LmC_Ldy_IIa = (iv.mpf('1.5') + 2 * g_abs_IIa) * L_dy_IIa
    r2_r3_Ldy_IIa = r2r3_IIa * L_dy_IIa
    t2_quad_IIa = (sa1_1B * (L_dy_IIa ** 2) / (2 * (1 - L_dy_IIa / N0_iv)) + sa2_1B * (L_dy_IIa ** 2) / (2 * (1 - r2_r3_Ldy_IIa / N0_iv)))
    Bmin_IIa = 2 - (2 * g_abs_IIa * (1 + Ry_IIa ** 2) + 1) * (Rz_Ry_IIa / omt_IIa)
    mIIa = eval_closed_iv("Case 1B-IIa (.15-.25n)", k_1B, 13, iv.mpf('0.25') * iv.pi, dx_IIa, dy_IIa, dz_IIa,
                          t1_IIa, LmC_Ldy_IIa, r2_r3_Ldy_IIa, t2_quad_IIa, Bmin_IIa, N0=N0)

    # Case 1B-IIb
    k_IIb = iv.pi * iv.sqrt(3)
    Rz_Ry_IIb = iv.exp(-iv.pi * iv_cot(iv.mpf('0.31') * iv.pi + iv.pi / (2 * N0_iv)))
    Ry_IIb = iv.exp(-(iv.mpf('0.07') * N0_iv - 2) * iv.pi * (iv.sqrt(2) - 1))
    omt_IIb = 1 - Ry_IIb
    dx_IIb, dy_IIb, dz_IIb = solve_rouche_iv(Ry_IIb * Ry_IIb * Rz_Ry_IIb, Rz_Ry_IIb, Ry_IIb, omt_IIb)
    L_dy_IIb = -iv.log(1 - dy_IIb)
    t1_IIb = (iv.mpf('0.25') * N0_iv * Ry_IIb * Ry_IIb * Rz_Ry_IIb) / ((1 - dx_IIb) * (omt_IIb - dx_IIb))
    r2r3_IIb = iv.cos(iv.pi / N0_iv) + iv_cot(iv.mpf('0.25') * iv.pi) * iv.sin(iv.pi / N0_iv)
    sa1_IIb, sa2_IIb = sa_iv(iv.mpf('0.25') * iv.pi)
    LmC_Ldy_IIb = iv.mpf('1.5') * L_dy_IIb
    r2_r3_Ldy_IIb = r2r3_IIb * L_dy_IIb
    t2_quad_IIb = (sa1_IIb * (L_dy_IIb ** 2) / (2 * (1 - L_dy_IIb / N0_iv)) + sa2_IIb * (L_dy_IIb ** 2) / (2 * (1 - r2_r3_Ldy_IIb / N0_iv)))
    g_max_IIb = r2r3_IIb * (-iv.cos(iv.mpf('0.62') * iv.pi + iv.pi / N0_iv))
    Bmin_IIb = 2 - (2 * g_max_IIb * (1 + Ry_IIb ** 2) + 1) * (Rz_Ry_IIb / omt_IIb)
    mIIb = eval_closed_iv("Case 1B-IIb (.25-.31n)", k_IIb, 22, iv.mpf('0.31') * iv.pi, dx_IIb, dy_IIb, dz_IIb,
                          t1_IIb, LmC_Ldy_IIb, r2_r3_Ldy_IIb, t2_quad_IIb, Bmin_IIb, N0=N0)

    # Case 1B-III even & odd
    k_III = iv.pi * iv.sqrt(3)
    r2r3_III = iv.cos(iv.pi / N0_iv) + iv_cot(iv.mpf('0.31') * iv.pi) * iv.sin(iv.pi / N0_iv)
    sa1_III, sa2_III = sa_iv(iv.mpf('0.31') * iv.pi)
    Rz_Ry_III = iv.exp(-iv.pi / iv.sqrt(3))
    LmC_III = iv.mpf('1.5') + 2 * iv.cos(iv.mpf('0.62') * iv.pi)
    g_max_III = iv.mpf('0.5') * r2r3_III
    cot_345 = iv_cot(iv.mpf('0.345') * iv.pi)

    # Even
    Ry_e = iv.exp(-2 * iv.pi * cot_345)
    Rz_e = iv.exp(-3 * iv.pi * iv_cot((iv.mpf(1) / 3 + iv.mpf(1) / (6 * N0_iv)) * iv.pi))
    omt_e = 1 - Ry_e
    dx_e, dy_e, dz_e = solve_rouche_iv(Ry_e * Rz_e, Rz_Ry_III, Ry_e, omt_e)
    L_dy_e = -iv.log(1 - dy_e)
    t1_e = (4 * iv.exp(-5 * iv.pi * cot_345)) / ((1 - dx_e) * (omt_e - dx_e))
    r2_r3_Ldy_e = r2r3_III * L_dy_e
    t2_quad_e = (sa1_III * (L_dy_e ** 2) / (2 * (1 - L_dy_e / N0_iv)) + sa2_III * (L_dy_e ** 2) / (2 * (1 - r2_r3_Ldy_e / N0_iv)))
    Bmin_e = 2 - (2 * g_max_III * (1 + Ry_e ** 2) + 1) * (Rz_Ry_III / omt_e)
    mIIIe = eval_closed_iv("Case 1B-III even", k_III, 27, iv.pi / 3, dx_e, dy_e, dz_e,
                           t1_e, LmC_III * L_dy_e, r2_r3_Ldy_e, t2_quad_e, Bmin_e, N0=N0)

    # Odd
    Ry_o = iv.exp(-iv.pi * cot_345)
    Rz_o = iv.exp(-2 * iv.pi / iv.sqrt(3))
    dx_o, dy_o, dz_o = solve_rouche_iv(Ry_o * Rz_o, Rz_Ry_III, iv.mpf(0), iv.mpf(1))
    L_dy_o = -iv.log(1 - dy_o)
    t1_o = (3 * iv.exp(-3 * iv.pi * cot_345)) / ((1 - dx_o) ** 2)
    r2_r3_Ldy_o = r2r3_III * L_dy_o
    t2_quad_o = (sa1_III * (L_dy_o ** 2) / (2 * (1 - L_dy_o / N0_iv)) + sa2_III * (L_dy_o ** 2) / (2 * (1 - r2_r3_Ldy_o / N0_iv)))
    Bmin_o = 2 - (2 * g_max_III + 1) * Rz_Ry_III
    mIIIo = eval_closed_iv("Case 1B-III odd", k_III, 27, iv.pi / 3, dx_o, dy_o, dz_o,
                           t1_o, LmC_III * L_dy_o, r2_r3_Ldy_o, t2_quad_o, Bmin_o, N0=N0)

    print(f"  >>> RIGOROUS INTERVAL GLOBAL WORST MARGIN (k0={k0_str}): {min(m1A, m1BI, mIIa, mIIb, mIIIe, mIIIo):+.6f} <<<")


def audit_theorem32_and_lemma45():
    print("\n=== INDEPENDENT 80-DPS CHECK OF THEOREM 3.2 (4-FACTOR IDENTITY) & LEMMA 4.5 (B_avg) ===")
    mp.mp.dps = 80
    for n, cases in [(16, [(2, 3), (3, 4), (4, 5)]), (20, [(2, 3), (5, 7), (4, 6)])]:
        coeffs = bk.compute_elementary_coeffs(n)
        for c, j in cases:
            r1 = mp.sin((j + c) * mp.pi / n)
            r2 = mp.sin(j * mp.pi / n)
            r3 = mp.sin(c * mp.pi / n)
            s = r2 / r1
            Ry = s ** n
            Rz = (r3 / r1) ** n
            t = Ry if (c % 2 == 0) else -Ry
            u = t / ((1 + t) ** 2)

            # Exact LHS from engine.py
            d = (n - 3 * c) // 2
            Q_val = mp.mpf(0)
            for b in range(d, -1, -1):
                Q_val = Q_val * u + mp.mpf(coeffs[(n - 3 * c - 2 * b, b, c)])
            lhs = ((-1) ** j) * ((1 + t) ** (n - 3 * c)) * Q_val

            # Compute Theorem 3.2 4-factor integral & Lemma 4.5 B_avg on 64 points
            th_c = c * mp.pi / n
            etas = []
            for k in range(n):
                et = (2 * k * mp.pi / n) if (c % 2 == 0) else ((2 * k + 1) * mp.pi / n)
                if et > mp.pi:
                    et -= 2 * mp.pi
                etas.append(et)
            small_e = [et for et in etas if abs(et) < th_c - mp.mpf('1e-20')]
            large_e = [et for et in etas if abs(et) > th_c + mp.mpf('1e-20')]
            assert len(small_e) == c - 1 and len(large_e) == n - c - 1

            prod_lg_abs = mp.mpf(1)
            for et in large_e:
                prod_lg_abs *= abs(1 - s * mp.exp(1j * et)) ** n
            C_mcj = Rz * prod_lg_abs / (Ry ** c)

            N_pts = 160
            acc_full = mp.mpc(0)
            acc_B = mp.mpc(0)
            for idx in range(N_pts):
                phi = -mp.pi + 2 * mp.pi * idx / N_pts
                w = mp.exp(1j * phi)
                # Solve T^3 - (1+t)T^2 + t*T - t*Rz*w = 0 by Newton starting from 1 and t
                def newton_cubic(T0):
                    T = mp.mpc(T0)
                    for _ in range(30):
                        f = T * (T - 1) * (T - t) - t * Rz * w
                        df = 3 * T * T - 2 * (1 + t) * T + t
                        T -= f / df
                    return T
                xv = newton_cubic(1)
                yv = newton_cubic(t)
                zv = (1 + t) - xv - yv
                dx = xv - 1
                dy = yv / (t * xv) - 1
                dz = zv / (xv * Rz * w) - 1
                # Verify exact Vieta inversion
                assert abs((1 + dz) * ((1 + dx) ** 3) * (1 + dy) - 1) < mp.mpf('1e-65')
                eps_s = (1 + dy) ** (mp.mpf(1) / n) - 1
                F1 = (1 + dx) ** (n - 3 * c)
                F2 = mp.mpc(1)
                for et in large_e:
                    ak = -s * mp.exp(1j * et) / (1 - s * mp.exp(1j * et))
                    F2 *= (1 + ak * eps_s) ** n
                F2 /= (1 + dy) ** c
                F3 = mp.mpc(1)
                z_over_x = zv / xv
                for et in small_e:
                    wk_t = (1 - s * (1 + eps_s) * mp.exp(1j * et)) ** n
                    F3 *= (1 - wk_t / z_over_x)
                for et in large_e:
                    wk_t = (1 - s * (1 + eps_s) * mp.exp(1j * et)) ** n
                    F3 *= (1 - z_over_x / wk_t)
                ep_p = (1 + (-s * mp.exp(1j * th_c) / (1 - s * mp.exp(1j * th_c))) * eps_s) ** n - 1
                ep_m = (1 + (-s * mp.exp(-1j * th_c) / (1 - s * mp.exp(-1j * th_c))) * eps_s) ** n - 1
                Z = ((-1) ** j) / w
                B = (2 - Z - 1 / Z) - (1 / Z) * dz + (ep_p + ep_m) - Z * ((1 + ep_p) * (1 + ep_m) * ((1 + dx) ** 3) * (1 + dy) - 1)
                acc_B += B
                acc_full += F1 * F2 * F3 * B

            rhs = (C_mcj * (acc_full / N_pts)).real
            rel_err = abs(rhs - lhs) / abs(lhs)
            B_num = (acc_B / N_pts).real
            dy_p0 = -(1 + t * t) * Rz / (t * (1 - t))
            dz_p0 = (1 - 2 * t * t) * Rz / (t * (1 - t))
            g_thc = -r2 * mp.cos((j + c) * mp.pi / n) / r3
            B_formula = 2 - ((-1) ** j) * (2 * g_thc * dy_p0 - dz_p0)
            b_err = abs(B_num - B_formula)
            print(f"  n={n:2d}, c={c}, j={j}: Thm3.2 rel_err = {float(rel_err):.2e}, Lemma4.5 B_avg err = {float(b_err):.2e}")
            assert rel_err < mp.mpf('1e-40') and b_err < mp.mpf('1e-40')


def main():
    audit_regime1_iv("3.0", N0=84)
    audit_regime1_iv("2.5", N0=84)
    audit_theorem32_and_lemma45()
    print("\n=== ALL REGIME 1 OUTWARD-ROUNDED INTERVAL BOUNDS & EXACT IDENTITIES CERTIFIED! ===")


if __name__ == "__main__":
    main()
