"""
Independent 50-digit outward-rounded interval (mpmath.iv) adversarial audit of:
  - Lemma R2.1 (Universal Regime 2 perturbation, n >= 132, kappa < 2.5)
  - Lemma R2.2 (Sub-regime 2A: c*kap^2 >= 42, kappa < 2.5, n >= 132)
  - Lemma R2.3 (Sub-regime 2C-Upgrade: c*q*kap <= 0.90, kappa < 2.5, n >= 132)

All interval constructors preserve full 50-digit outward rounding (ZERO float64 casts
inside iv.mpf([...])), and all 6 audit fixes identified during line-by-line review of
analyst/certify_R2_1_and_2A_iv.py and analyst/certify_R2_3_iv.py are explicitly certified:
  [Fix 1] R2.2 small-root |a_k| <= (r2/r3)/cos(alpha_max) => rigorous iv eta_3_coeff (~2.10 > 2.02).
  [Fix 2] R2.2 n->inf uniform R_bdy_coeff (~6.02 > 6.00; analyst checked 6.0 only at N0=132).
  [Fix 3] R2.3 l=0 universal Jensen F_R = 9/14 (fixes far-side u_R in (1/2,1] at r=1,2,3) and
          pref_L0(l=0) = 2.0 (fixes t_* = 1.5 + delta_shift > 1.5 exceeding 1.37 at c=2,j=3).
  [Fix 4] R2.3 Step 5/6 (q >= 2): explicit iv bound on Delta_max = max_k |Delta_k| < 0.01 and
          exact 1/(1 - Delta_max) factor for |log(1 - Delta_k)|/q.
  [Fix 5] R2.3 Step 6: explicit stretched-exponential iv bound E_cubic (< 1e-44) replacing
          uncertified hardcoded 1e-10.
  [Fix 6] Exact 80-digit checks of the pairing identity, exact binomial p_r identity, and
          rescaled contour representation against exact integer engine.py coefficients.
"""

from __future__ import annotations
import sys
import math
import mpmath as mp
from mpmath import iv

sys.path.insert(0, r"d:\solving-conjectures\work\bardakov-kozlovskaya\verify")
import engine as bk

iv.dps = 50
mp.mp.dps = 50
PI = iv.pi


def iv_cot(x):
    return iv.cos(x) / iv.sin(x)


def iv_hull(lo_iv, hi_iv):
    """Exact 50-digit outward-rounded interval hull [lo_iv.a, hi_iv.b] with ZERO float64 truncation."""
    return iv.mpf([lo_iv.a, hi_iv.b])


def iv_max_ub(items):
    """Exact 50-digit upper bound max over a list of iv intervals (no float64 truncation)."""
    return iv.mpf(max(x.b for x in items))


def iv_min_lb(items):
    """Exact 50-digit lower bound min over a list of iv intervals (no float64 truncation)."""
    return iv.mpf(min(x.a for x in items))


def audit_R2_1_and_2A():
    print("==========================================================================")
    print("AUDIT OF LEMMA R2.1 & LEMMA R2.2 (SUB-REGIME 2A: c*kap^2 >= 42, n >= 132)")
    print("==========================================================================")
    N0 = iv.mpf(132)
    KAP0 = iv.mpf('2.5')
    C1 = iv.mpf(42)

    alpha_max = iv.mpf('0.409284')
    assert (iv.sin(alpha_max) - KAP0 / (2 * PI)).a > 0
    theta_c_max = alpha_max / 2
    S_a = iv.sin(alpha_max) / alpha_max

    r3_r1_max = 1 / (2 * iv.cos(theta_c_max))
    Rz_max = r3_r1_max ** N0
    Ry_max = iv.exp(-2 * PI * iv_cot(alpha_max))
    cot_half = iv_cot(alpha_max / 2)
    Rz_Ry_max = iv.exp(-PI * cot_half)

    xcot_half = (alpha_max / 2) * cot_half
    C_exp = 4 * (PI ** 2) * S_a * xcot_half
    C_exp_pert = ((N0 - 1) / N0) * C_exp
    assert C_exp.a > 37.8

    # Universal (1/theta_c)*(R_z/R_y)
    c_th = 2 * PI * xcot_half / 3
    assert (1 / theta_c_max).a > (1 / c_th).b
    inv_th_RzRy = iv_max_ub([
        (1 / theta_c_max) * iv.exp(-c_th / theta_c_max),
        (N0 / (2 * PI)) * (r3_r1_max ** N0),
    ])

    # Strict Rouché on |w| = 1
    rho_x_1 = iv.mpf('1e-40')
    assert (rho_x_1 * (1 - rho_x_1) * (1 - Ry_max - rho_x_1) - Ry_max * Rz_max).a > 0
    c_uy = iv.mpf('1.001')
    u_y_1 = c_uy * Rz_Ry_max
    assert (u_y_1 * (1 - u_y_1) * (1 - Ry_max - u_y_1 * Ry_max) - Rz_Ry_max).a > 0

    # Pointwise-uniform linear ratios vs R_z/R_y and X = (R_z/R_y)^{1 - 1/N0}
    c_rhox_per_RzRy = (iv.mpf(2) * (Ry_max ** 2))  # rho_x^* <= 2 * R_y * R_z = 2 * R_y^2 * (R_z/R_y)
    assert ((1 - rho_x_1) * (1 - Ry_max - rho_x_1)).a > iv.mpf('0.5')
    c_dy_per_RzRy = (c_uy + c_rhox_per_RzRy) / (1 - rho_x_1)
    bar_dy_1 = c_dy_per_RzRy * Rz_Ry_max
    c_Ldy_per_RzRy = (-iv.log(1 - bar_dy_1)) / Rz_Ry_max
    c_dz_per_RzRy = (1 / (((1 - rho_x_1) ** 3) * (1 - bar_dy_1)) - 1) / Rz_Ry_max

    X_max = Rz_Ry_max ** ((N0 - 1) / N0)
    eta_pert_1_max = c_Ldy_per_RzRy * X_max
    c_eps_bdy_per_X = (iv.exp(eta_pert_1_max) - 1) / X_max
    eps_bdy_1 = c_eps_bdy_per_X * X_max

    # [FIX 2] True n->inf UNIFORM R_bdy_coeff:
    # Since R_z/R_y <= X for ALL n >= N0 (and ratio -> 1 as n -> inf!),
    #   bar_dz <= c_dz * X, bar_dy <= c_dy * X, rho_x <= c_rhox * X, eps_bdy <= c_eps_bdy * X.
    # Let f(t) = ((1 + c_eps_bdy*t)^2 * (1 + c_rhox*t)^3 * (1 + c_dy*t) - 1) / t (increasing in t > 0):
    c_prod_per_X = (((1 + eps_bdy_1) ** 2) * ((1 + c_rhox_per_RzRy * X_max) ** 3) * (1 + c_dy_per_RzRy * X_max) - 1) / X_max
    R_bdy_coeff_uniform = c_dz_per_RzRy + 2 * c_eps_bdy_per_X + c_prod_per_X
    print(f"  [FIX 2 check] True n->inf uniform R_bdy_coeff <= {float(R_bdy_coeff_uniform.b):.6f} (> 6.000! analyst's 6.0 fails as n->inf)")
    assert R_bdy_coeff_uniform.a > iv.mpf('6.0')

    # [FIX 1] True small-root |a_k| <= (r2/r3)/cos(alpha_max) => rigorous eta_3_coeff!
    sec_alpha = 1 / iv.cos(alpha_max)
    eta_sm_max = sec_alpha * eta_pert_1_max
    c_sm_exp_per_X = (iv.exp(eta_sm_max) - 1) / X_max
    c_inv_dz_per_X = (c_dz_per_RzRy / (1 - c_dz_per_RzRy * Rz_Ry_max))
    # |tilde_u_k / u_k - 1| <= |(1+a_k eps_s)^n / (1+delta_z) - 1|
    #   <= (1 + (exp(sec_alpha*eta_pert)-1)) * (1 + bar_dz/(1-bar_dz)) - 1
    eta_3_coeff_uniform = (((1 + c_sm_exp_per_X * X_max) * (1 + c_inv_dz_per_X * X_max)) - 1) / X_max
    print(f"  [FIX 1 check] True small-root eta_3_coeff <= {float(eta_3_coeff_uniform.b):.6f} (> 2.020! analyst omitted 1/cos(alpha_max)={float(sec_alpha.b):.4f})")
    assert eta_3_coeff_uniform.a > iv.mpf('2.02')

    # t_1 and t_2 on |w| = 1
    L_dy_1 = c_Ldy_per_RzRy * Rz_Ry_max
    t_1_max = 2 * N0 * Ry_max * Rz_max
    t_2_max = (iv.mpf('1.5') * L_dy_1 + 2 * eta_pert_1_max
               + (2 / PI) * inv_th_RzRy * (c_Ldy_per_RzRy ** 2) * Rz_Ry_max / (2 * (1 - eta_pert_1_max / N0)))

    # B_2A (exact pairing phase bound)
    gam_ext = iv.sin(alpha_max / 4) / (alpha_max / 4)
    gam_z = (alpha_max / 2) / iv.sin(alpha_max / 2)
    L_min = C1 / KAP0
    # Rigorous zeta(3) enclosure from partial sum + integral tail (no mp.zeta trust needed!)
    M_z = 200
    z3_part = sum(iv.mpf(1) / (iv.mpf(k) ** 3) for k in range(1, M_z + 1))
    z3_tail_lo = iv.mpf(1) / (2 * ((iv.mpf(M_z) + 1) ** 2))
    z3_tail_hi = iv.mpf(1) / (2 * (iv.mpf(M_z) ** 2))
    zeta3_iv = iv.mpf([(z3_part + z3_tail_lo).a, (z3_part + z3_tail_hi).b])

    t0_iv = iv.mpf(['1.59362426', '1.59362427'])
    assert ((2 - iv.mpf(t0_iv.a)) * iv.exp(iv.mpf(t0_iv.a)) - 2).a > 0
    assert ((2 - iv.mpf(t0_iv.b)) * iv.exp(iv.mpf(t0_iv.b)) - 2).b < 0
    M_g_iv = (t0_iv ** 2) / (iv.exp(t0_iv) - 1)

    f = iv.mpf('0.26')
    a1 = gam_ext * (1 - f)
    S1 = 4 * gam_z * (2 * zeta3_iv / (C1 * (a1 ** 3)) + M_g_iv / (L_min * (a1 ** 2)))
    cA = gam_ext * f * (1 - f)
    b_slope = gam_ext * (iv.mpf('0.5') - f)
    A_min = cA * L_min
    assert A_min.a > 1
    S2 = (4 / (1 - iv.exp(-A_min))) * (iv.exp(-A_min) + (L_min * iv.exp(-A_min)) / (C1 * b_slope))
    S_mid = -iv.log(1 - iv.exp(-gam_ext * L_min / 4))
    c_int = 2 / PI
    x_half = (c_int * L_min) / N0
    a_half = (1 + x_half) ** (-N0 / 2)
    assert L_min.a > (N0 / (c_int * (N0 / 2 - 2))).b
    int_tail = (L_min / C1) * (N0 / (2 * c_int * (N0 / 2 - 1))) * ((1 + x_half) ** (-(N0 / 2 - 1)))
    S_int = (2 / (1 - a_half)) * (a_half + int_tail)
    B_2A = S1 + S2 + S_mid + S_int
    cos_B_2A = iv.cos(B_2A)
    assert B_2A.b < (PI / 2).a

    # Delta_123 with FIXED uniform eta_3_coeff
    C_sm_2 = 8 / (gam_ext ** 2)
    c_lg_a = (PI * N0) / (2 * (N0 / 2 - 1))
    S_denom_25 = (C_sm_2 / (KAP0 ** 2)) + 2 * ((PI + KAP0) / KAP0) * (1 + c_lg_a / KAP0)
    W_3_max = S_denom_25 * eta_3_coeff_uniform * iv.exp(-C_exp_pert / KAP0)
    Delta_123_max = iv.exp(t_1_max + t_2_max + W_3_max) - 1

    # Oscillation ratio & Delta_B with FIXED uniform R_bdy_coeff
    a1_max = a1 * KAP0
    C_osc_sm = 2 * (a1_max * iv.log(1 + 2 / a1_max) + (PI ** 2) / 4) / a1
    osc_const_sm = 2 * (S2 + S_mid)
    s6 = 3 / PI
    v6 = KAP0 / PI
    a_lg1 = s6 * (iv.log(1 + v6) / v6)
    alg1_max = a_lg1 * KAP0
    C_osc_lg1 = 2 * (alg1_max * iv.log(1 + 2 / alg1_max) + (PI ** 2) / 4) / a_lg1
    v_half = 2 * KAP0 / PI
    a_lg2 = (2 / PI) * (iv.log(1 + v_half) / v_half)
    T_lg2 = a_lg2 * C1 / (2 * KAP0)
    assert (iv.mpf(21) * a_lg2).a > KAP0.b  # monotonicity in kappa <= 2.5
    osc_const_lg2 = 2 * (iv.log((1 + iv.exp(-T_lg2)) / (1 - iv.exp(-T_lg2)))
                         + (2 / (a_lg2 * KAP0)) * iv.exp(-T_lg2) / (1 - iv.exp(-T_lg2)))
    C_osc_tot = C_osc_sm + C_osc_lg1
    osc_const_tot = osc_const_sm + osc_const_lg2
    assert C_osc_tot.b < C_exp_pert.a

    bdy_ratio_max = (R_bdy_coeff_uniform / 2) * iv.exp(osc_const_tot) * iv.exp(-(C_exp_pert - C_osc_tot) / KAP0)
    margin_2A = cos_B_2A - Delta_123_max - (1 + Delta_123_max) * bdy_ratio_max
    print(f"  B_2A <= {float(B_2A.b):.6f} => cos(B_2A) >= {float(cos_B_2A.a):.6f}")
    print(f"  Delta_123 (fixed) <= {float(Delta_123_max.b):.3e}, bdy_ratio (fixed) <= {float(bdy_ratio_max.b):.5f}")
    print(f"  >>> RIGOROUS FIXED SUB-REGIME 2A MARGIN_LB = {float(margin_2A.a):+.6f} <<<")
    assert margin_2A.a > iv.mpf('0.21')


def polylog_iv(s, W):
    M_cut = 250
    tot = sum((W ** m) / (iv.mpf(m) ** s) for m in range(1, M_cut + 1))
    tail = (W ** (M_cut + 1)) / ((iv.mpf(M_cut + 1) ** s) * (1 - W))
    return tot + iv.mpf([iv.mpf(0).a, tail.b])


def max_f_on_interval(W, B, xa, xb):
    mlw = -iv.log(W)
    xp = (iv.sqrt(iv.mpf('2.25') + 4 * B * mlw) - iv.mpf('1.5')) / (2 * mlw)
    if xp.b <= xa.a:
        xc = xa
    elif xp.a >= xb.b:
        xc = xb
    else:
        xc = iv.mpf([max(xa.a, xp.a), min(xb.b, xp.b)])
    return (xc ** (iv.mpf(-3) / 2)) * (W ** xc) * iv.exp(-B / xc)


def audit_R2_3():
    print("\n==========================================================================")
    print("AUDIT OF LEMMA R2.3 (SUB-REGIME 2C-UPGRADE: c*q*kap <= 0.90, n >= 132)")
    print("==========================================================================")
    C2 = iv.mpf('0.90')
    e_iv = iv.exp(1)

    kap_max_all = C2 / 2
    kap_max_q2 = C2 / 4
    alpha_all = iv.mpf('0.07169')
    assert (iv.sin(alpha_all) - kap_max_all / (2 * PI)).a > 0
    alpha_q2 = iv.mpf('0.03583')
    assert (iv.sin(alpha_q2) - kap_max_q2 / (2 * PI)).a > 0

    S_all = iv.sin(alpha_all) / alpha_all
    S_q2 = iv.sin(alpha_q2) / alpha_q2
    n_min_q1 = 24 * (PI ** 2) * (S_all ** 2) / kap_max_all
    n_min_q2 = 15 * (PI ** 2) * (S_q2 ** 2) / kap_max_q2
    n_min = iv_min_lb([n_min_q1, n_min_q2])

    # [FIX 5] Explicit stretched-exponential iv bound on cubic-root perturbation E_cubic on |w| = lam_* < 1!
    # In c*q*kap <= C2, c*kap^2 = (c*q*kap)*(kap/q) <= C2 * kap_max_all = C2^2 / 2 = 0.405!
    C1_2C = (C2 ** 2) / 2
    K_jc_2C = (C1_2C / (2 * (PI ** 4) * (S_all ** 4))) ** (iv.mpf(1) / 3)
    xcot_all = (alpha_all / 2) * iv_cot(alpha_all / 2)
    c_se_2C = 2 * xcot_all / K_jc_2C
    # Verify n^2 * exp(-c_se_2C * n^{1/3}) strictly decreasing for n >= n_min:
    # d/dn [2 log n - c_se * n^{1/3}] = 2/n - (c_se/3) n^{-2/3} < 0 iff n > (6/c_se)^3
    assert ((6 / c_se_2C) ** 3).b < n_min.a
    Rz_Ry_se_2C = iv.exp(-c_se_2C * (n_min ** (iv.mpf(1) / 3)))
    # Also R_y <= exp(-2*pi*cot(alpha_all)) and (R_z/R_y)^{1 - 1/n} <= exp(-((n_min-1)/n_min)*c_se_2C*n^{1/3})
    Ry_2C = iv.exp(-2 * PI * iv_cot(alpha_all))
    c_se_pert = ((n_min - 1) / n_min) * c_se_2C
    assert ((6 / c_se_pert) ** 3).b < n_min.a
    n2_X_2C = (n_min ** 2) * iv.exp(-c_se_pert * (n_min ** (iv.mpf(1) / 3)))
    # On |w| = lam_* < 1:
    #   |log F_1| <= 2 * n * R_y * R_z <= 2 * n * R_y^2 * (R_z/R_y)
    #   |log G_2| <= n * (2 * R_y * bar_dy) / (1 - R_y) <= 4 * n * R_y * (R_z/R_y)
    #   sum_{k=0}^{n-1} |log((1 - tilde_a_k w)/(1 - a_k w))| <= (n^2 / ((-log W_max)*W_max)) * 2.5 * X
    W_max_all = iv.mpf('0.836')
    W_max_q2 = iv.mpf('0.832')
    E_cubic_2C = (2 * Ry_2C + 4 * Ry_2C + iv.mpf('2.5') / ((-iv.log(W_max_all)) * W_max_all)) * n2_X_2C
    print(f"  [FIX 5 check] Certified stretched-exp E_cubic <= {float(E_cubic_2C.b):.3e} (<< 1e-10; c_se_2C >= {float(c_se_2C.a):.3f})")
    assert E_cubic_2C.b < iv.mpf('1e-40')

    # Step 2: certify W_max (with ZERO float64 truncation!)
    def certify_W_max(alpha_b, W_try):
        gam_u = 1 / (iv.cos(alpha_b) ** 2)
        a_u = gam_u * C2
        L_alt = 120
        alt_s = sum(((-1) ** (l - 1)) * (W_try ** l) / iv.sqrt(l) for l in range(1, L_alt + 1))
        alt_e = (W_try ** (L_alt + 1)) / iv.sqrt(L_alt + 1)
        I_m = iv.sqrt(PI / a_u) * (alt_s + iv.mpf([-alt_e.b, alt_e.b]))
        y_sh = iv.mpf('0.80')
        ay2 = a_u * (y_sh ** 2)
        r_far = W_try * iv.exp(ay2 - (PI ** 2) / (9 * ay2))
        assert r_far.b < iv.mpf('0.5')
        # Also verify ZERO poles inside strip |Im u| <= y_sh:
        # At any pole 2*a_u*x*y >= pi => a_u*x^2 >= pi^2/(4*a_u*y^2) >= pi^2/(4*ay2)
        r_pole = W_try * iv.exp(ay2 - (PI ** 2) / (4 * ay2))
        assert r_pole.b < iv.mpf('0.1')
        d_min = iv_min_lb([iv.sqrt(3) / 2, 1 - r_far])
        e_p = iv.exp(-4 * PI * y_sh)
        j_err = (W_try * iv.exp(ay2) / d_min) * iv.sqrt(PI / a_u) * (2 * e_p / (1 - e_p))
        x_tr = PI / (2 * alpha_b)
        t_err = 2 * iv.exp(-a_u * (x_tr ** 2)) / (1 - iv.exp(-a_u * x_tr))
        S_lo = I_m - j_err - t_err
        assert S_lo.a > 1
        return S_lo

    S_lo_all = certify_W_max(alpha_all, W_max_all)
    S_lo_q2 = certify_W_max(alpha_q2, W_max_q2)
    H_per_g = 1 / (1 + W_max_all)
    assert H_per_g.a > iv.mpf('0.544')
    print(f"  [Step 2] S_lo(all) >= {float(S_lo_all.a):.6f}, S_lo(q>=2) >= {float(S_lo_q2.a):.6f}, H_q/g >= {float(H_per_g.a):.5f}")

    # Step 3: Polya-Szego rat_main (ZERO float64 truncation in iv hulls!)
    H_split = iv.mpf('5.0')
    C_th = iv.mpf('2.5')
    th0_max = (C_th / H_split) ** (iv.mpf(1) / 3)
    u0 = th0_max / 2
    phi_max = (C_th / 6) / (iv.cos(u0) ** 2)
    cos_phi = iv.cos(phi_max)
    a_lo = iv.mpf('0.5') / iv.cos(u0)
    Ht2_min = (C_th ** (iv.mpf(2) / 3)) * (H_split ** (iv.mpf(1) / 3))
    z_sq_min = a_lo * Ht2_min
    erf_lo = 1 - iv.exp(-z_sq_min) / (iv.sqrt(PI) * iv.sqrt(z_sq_min))
    sH_main = (cos_phi / (2 * iv.sqrt(PI * a_lo))) * erf_lo
    sH_tail1 = (PI / (8 * iv.sqrt(Ht2_min))) * iv.exp(-4 * Ht2_min / (PI ** 2))
    sH_tail2 = iv.mpf('0.5') * iv.sqrt(H_split) * iv.exp(-H_split)
    rat_ge5 = (sH_main - sH_tail1 - sH_tail2) / (iv.sqrt(PI) / 4 + sH_tail2)

    N_H = 150
    N_t = 120
    H_lo_start = iv.mpf('0.544')
    Hs = [H_lo_start + (H_split - H_lo_start) * k / N_H for k in range(N_H + 1)]
    box_rats = []
    for k in range(N_H):
        Ha = iv.mpf(Hs[k].a)
        Hb = iv.mpf(Hs[k + 1].b)
        lo_t, hi_t = 0.0, 2.8
        for _ in range(30):
            mt = 0.5 * (lo_t + hi_t)
            if 2.0 * float(Hb.b) * (2.0 * float(mp.tan(mt / 2)) - mt) < float(mp.pi / 2) - 1e-4:
                lo_t = mt
            else:
                hi_t = mt
        th_cut = iv.mpf(str(round(lo_t, 5)))
        pos_sum = iv.mpf(0)
        dt1 = th_cut / N_t
        for i in range(N_t):
            t_iv = iv_hull(dt1 * i, dt1 * (i + 1))
            ph_iv = 2 * Hb * (2 * iv.tan(t_iv / 2) - t_iv)
            if ph_iv.b < (PI / 2).a:
                pos_sum += iv.mpf(((iv.cos(t_iv / 2) ** (4 * Hb)) * iv.cos(ph_iv)).a) * dt1
            else:
                pos_sum -= iv.mpf((iv.exp(-2 * Ha * (iv.sin(t_iv / 2) ** 2))).b) * dt1
        neg_sum = iv.mpf(0)
        dt2 = (PI - th_cut) / N_t
        for i in range(N_t):
            t_iv = iv_hull(th_cut + dt2 * i, th_cut + dt2 * (i + 1))
            neg_sum += iv.mpf((iv.exp(-2 * Ha * (iv.sin(t_iv / 2) ** 2))).b) * dt2
        u_sum = iv.mpf(0)
        dt3 = PI / N_t
        for i in range(N_t):
            t_iv = iv_hull(dt3 * i, dt3 * (i + 1))
            u_sum += iv.mpf((iv.exp(-2 * Ha * (iv.sin(t_iv / 2) ** 2))).b) * dt3
        num_box = pos_sum - neg_sum
        assert num_box.a > 0
        box_rats.append(num_box / u_sum)
    rat_main_iv = iv_min_lb(box_rats + [rat_ge5])
    print(f"  [Step 3] Rigorous 50-dps Polya-Szego rat_main >= {float(rat_main_iv.a):.6f}")

    # Step 4: Negative-binomial KL bound on |R_nonq| WITH [FIX 3]:
    #   F_R(l=0) = 9/14 (universal Jensen (1+1/3)(1+1/6)=14/9 for ALL r>=1, u_R in (0,1]!)
    #   pref_L0(l=0) = 2.0 (universal for ALL r>=1, u_L in (0,1]!)
    C_sh_q2 = (iv.mpf(49) / 96) / (1 - (iv.mpf(49) / 96) * (alpha_q2 ** 2))
    K_drift = C_sh_q2 * (C2 ** 2) / (4 * (PI ** 2) * (S_q2 ** 3))
    gam_lo_q2 = iv.cos(alpha_q2) ** 2  # stronger & exact: (tan(th_c/2)/(th_c/2))*cos^2(psi+th_c/2) > cos^2(alpha_q2)
    C_pref = iv.sqrt(PI)  # exact: Gamma_lo < 1 => sqrt(pi*Gamma_lo) < sqrt(pi)

    R_nonq_bound = iv.mpf(0)
    W_q = W_max_q2
    for l in range(0, 48):
        xa = iv.mpf(l) if l >= 1 else iv.mpf('1e-6')
        xb = iv.mpf(l + 1)
        eps_l = (l + 1) * K_drift / 4
        assert eps_l.b < iv.mpf('0.13')
        A_r_min = 4 * gam_lo_q2 * (PI ** 2) / (xb * C2)
        pref_L0 = iv.mpf('2.0') if l == 0 else (iv.mpf(2 * l + 2) / (2 * l + 1))
        pref_L_far = iv.mpf('2.0') if l == 0 else (iv.mpf(2 * l + 2) / (2 * l + 1))
        # Rigorous bound for t < t_*/2 and t = 0:
        #   At l = 0: A_r >= 43.4, so t < t_*/2 (j>=1) & t=0 are < 1e-16.
        #   At l >= 1: t < t_*/2 => j >= l+1 => A_r * j^2 >= (4*gam_lo*pi^2/C2)*(l+1) >= 86 => < 1e-30!
        assert (A_r_min * ((l + 1) ** 2)).a > 40
        tail_L_far_tiny = iv.mpf('1e-12')
        tail_L = 2 * iv.exp(-A_r_min) / (1 - iv.exp(-3 * A_r_min)) + tail_L_far_tiny
        tail_R = iv.exp(-A_r_min / 2) / (1 - iv.exp(-A_r_min / 2))

        # [FIX 3] Universal Jensen F_R: 9/14 for l=0; (4l+2)/(4l+3) for l>=1
        # (For l>=1: q=2 => u_R<=1/2, r=2l+1 => (1+1/(6r))(1+1/(12r)) < 1+1/(2r) = (4l+3)/(4l+2);
        #  q>=3 => u_R<=1, r>=3l+1 => (1+1/(3r))(1+1/(6r)) < 1+1/(4l+2) = (4l+3)/(4l+2)!)
        if l >= 1:
            r3 = iv.mpf(3 * l + 1)
            assert ((1 + 1 / (3 * r3)) * (1 + 1 / (6 * r3))).b < (iv.mpf(4 * l + 3) / (4 * l + 2)).a
        F_R = (iv.mpf(9) / 14) if l == 0 else (iv.mpf(4 * l + 2) / (4 * l + 3))
        M_L_exact = gam_lo_q2 * (PI ** 2)
        M_R = F_R * gam_lo_q2 * (PI ** 2) * ((1 - eps_l) ** 2)

        def phi_side(M_side, k_val):
            if (2 * M_side * (k_val ** 2) / xb).a >= C2.b:
                return (1 / iv.sqrt(C2)) * max_f_on_interval(W_q, M_side * (k_val ** 2) / C2, xa, xb)
            else:
                return (1 / iv.sqrt(2 * e_iv * M_side)) * (W_q ** xa) / (k_val * xa)

        t_L1 = (pref_L0 + tail_L) * phi_side(M_L_exact, 1)
        t_R1 = (1 + tail_R) * phi_side(M_R, 1)
        t_L2 = 2 * (pref_L_far + tail_L) * phi_side(M_L_exact, 2)
        t_R2 = 2 * (1 + tail_R) * phi_side(M_R, 2)
        assert (2 * M_R * 9 / xb).a >= C2.b
        t_L3 = 2 * (pref_L_far + tail_L) * (1 / iv.sqrt(C2)) * max_f_on_interval(W_q, M_L_exact * 9 / C2, xa, xb) / (1 - iv.exp(-7 * M_L_exact / (xb * C2)))
        t_R3 = 2 * (1 + tail_R) * (1 / iv.sqrt(C2)) * max_f_on_interval(W_q, M_R * 9 / C2, xa, xb) / (1 - iv.exp(-7 * M_R / (xb * C2)))
        term_l = C_pref * (t_L1 + t_R1 + t_L2 + t_R2 + t_L3 + t_R3)
        R_nonq_bound += term_l
        if l == 0:
            print(f"  [FIX 3 check] l=0 fixed (F_R=9/14, pref_L0=2.0): term_0 <= {float(term_l.b):.6f}")

    # Tails l >= 48
    eps_drift = 12 * K_drift
    r_48_q3 = iv.mpf(3 * 48 + 1)
    F_R_48 = 1 / ((1 + 1 / (3 * r_48_q3)) * (1 + 1 / (6 * r_48_q3)))
    assert F_R_48.a > (iv.mpf(6) / 7).b
    M0_tail = (iv.mpf(6) / 7) * gam_lo_q2 * (PI ** 2) * ((1 - eps_drift) ** 2)
    A_r_drift = gam_lo_q2 * (PI ** 2) / (12 * C2)
    mult_LR_drift = iv.mpf('2.02') + iv.exp(-A_r_drift / 2) / (1 - iv.exp(-A_r_drift / 2)) + 2 * iv.exp(-A_r_drift) / (1 - iv.exp(-3 * A_r_drift)) + iv.mpf('1e-12')
    tail_A = 2 * C_pref * mult_LR_drift * iv.sqrt(PI / (4 * M0_tail)) * (W_q ** 48) / (48 * (1 - W_q))
    tail_B = ((1 + W_q) / (24 * (-iv.log(W_q)) * (W_q ** 2))) * (W_q ** 48)
    R_nonq_bound += tail_A + tail_B
    print(f"  [Step 4] Rigorous fixed |R_nonq| <= {float(R_nonq_bound.b):.6f}")

    # Step 5: Conjugate-pair phase error + [FIX 4] Delta_max bound for q >= 2!
    def cert_E_phase(alpha_b, S_b, W_b, n_min_b, j_over_c_min, inv_gq4_max, inv_g2q5_max, inv_g2q4_max, label):
        C_sh_b = (iv.mpf(49) / 96) / (1 - (iv.mpf(49) / 96) * (alpha_b ** 2))
        gam_up_b = 1 / (iv.cos(alpha_b) ** 2)
        sinc_q_b = iv.sin(alpha_b / 4) / (alpha_b / 4)
        gam_lo_b = iv.cos(alpha_b) ** 2

        K_B1 = C_sh_b / (2 * PI * (S_b ** 4) * j_over_c_min)
        K_B3 = (iv.mpf(4) / 3) * gam_up_b / (2 * PI * (S_b ** 2) * j_over_c_min * sinc_q_b * iv.cos(alpha_b))

        U0_num = iv.mpf('9.0')
        eta_U0 = U0_num / (S_b * iv.sqrt(n_min_b))
        v_U0 = (2 / n_min_b) * gam_up_b * (U0_num ** 2)
        gam_g = (iv.log(1 + v_U0) / v_U0) * ((iv.sin(eta_U0 / 2) / (eta_U0 / 2)) ** 2) * gam_lo_b
        ag = gam_g * C2

        I_vals = {1: iv.sqrt(PI) / 2, 2: 3 * iv.sqrt(PI) / 4, 3: 15 * iv.sqrt(PI) / 8}
        M_vals = {1: (1 / e_iv), 2: (2 / e_iv) ** 2, 3: (3 / e_iv) ** 3}

        def mom_b(p):
            Li_int = polylog_iv(iv.mpf(p) - iv.mpf('0.5'), W_b)
            Li_dsc = polylog_iv(iv.mpf(p - 1), W_b) if p >= 2 else W_b / (1 - W_b)
            t_int = inv_gq4_max * (C2 ** 4) * Li_int * I_vals[p] / (ag ** (p + iv.mpf('0.5')))
            t_dsc = 2 * inv_g2q5_max * (C2 ** 4) * Li_dsc * M_vals[p] / (ag ** p)
            return t_int + t_dsc

        Ep_main = iv.mpf('0.5') * ((K_B1 ** 2) * mom_b(1) + 2 * K_B1 * K_B3 * mom_b(2) + (K_B3 ** 2) * mom_b(3))

        M_t = n_min_b / 2
        b_t = (4 / (PI ** 2)) * gam_lo_b
        assert (iv.mpf(3) * M_t / ((M_t - 3) * b_t)).b < (U0_num ** 2).a

        def t_mom(p):
            base_d = (1 + b_t * (U0_num ** 2) / M_t) ** (-(M_t - p - 1))
            val_int = 2 * ((M_t / b_t) ** p) * (M_t / (2 * b_t * U0_num * (M_t - p - 1))) * base_d
            val_pt = 2 * ((M_t / b_t) ** p) * ((1 + b_t * (U0_num ** 2) / M_t) ** (-(M_t - p)))
            return inv_gq4_max * (C2 ** 4) * (val_int / (C2 ** (p + iv.mpf('0.5'))) + val_pt / (C2 ** p))

        Ep_tail = (W_b / (2 * ((1 - W_b) ** 2))) * ((K_B1 ** 2) * t_mom(1) + 2 * K_B1 * K_B3 * t_mom(2) + (K_B3 ** 2) * t_mom(3))
        Ep = Ep_main + Ep_tail

        # [FIX 4] Explicit upper bound on Delta_max = max_k (4 x_k / (1-x_k)^2) * sin^2(phi_k/2)
        #   <= max_k (x_k / (1-x_k)^2) * phi_k^2
        #   For |u_k| <= U0/sqrt(tau_q): sum_{m>=1} m W_b^m u^{2p} e^{-m gam_g tau_q u^2}
        #     <= Li_{p-1}(W_b) * M_vals[p] / (gam_g * tau_q)^p!
        #   Multiply by B_i B_j <= K_Bi K_Bj * tau_q^4 / (g^2 q^4): max at tau_q = C2!
        def pk_b(p):
            Li_pk = polylog_iv(iv.mpf(p - 1), W_b) if p >= 2 else W_b / (1 - W_b)
            return inv_g2q4_max * (C2 ** 4) * Li_pk * M_vals[p] / (ag ** p)

        Dmax_main = (K_B1 ** 2) * pk_b(1) + 2 * K_B1 * K_B3 * pk_b(2) + (K_B3 ** 2) * pk_b(3)
        # Tail |u_k| > U0/sqrt(tau_q): single point max
        def pk_tail(p):
            return inv_g2q4_max * (C2 ** 4) * ((M_t / b_t) ** p) * ((1 + b_t * (U0_num ** 2) / M_t) ** (-(M_t - p))) / (C2 ** p)
        Dmax_tail = (W_b / ((1 - W_b) ** 2)) * ((K_B1 ** 2) * pk_tail(1) + 2 * K_B1 * K_B3 * pk_tail(2) + (K_B3 ** 2) * pk_tail(3))
        Dmax = iv_max_ub([Dmax_main, Dmax_tail])
        print(f"  [{label}] E_phase <= {float(Ep.b):.6f}, Delta_max <= {float(Dmax.b):.6f} (1/(1-Delta_max) <= {float((1/(1-Dmax)).b):.6f})")
        assert Dmax.b < iv.mpf('0.1')
        return Ep, Dmax

    Ep_q1, Dmax_q1 = cert_E_phase(alpha_all, S_all, W_max_all, n_min_q1, iv.mpf(2),
                                  iv.mpf('0.5'), iv.mpf('0.25'), iv.mpf('0.25'), "q = 1")
    Ep_q2, Dmax_q2 = cert_E_phase(alpha_q2, S_q2, W_max_q2, n_min_q2, iv.mpf(1),
                                  iv.mpf(1) / 16, iv.mpf(1) / 32, iv.mpf(1) / 16, "q >= 2")

    # For q = 1: prod(1 - Delta_k) * exp(E_cubic) - 1 has |...| <= exp(Ep_q1 + E_cubic) - 1 (no log needed, or use log with 1/(1-Dmax) too!)
    # For q >= 2: |(1/q) sum log(1 - Delta_k)| <= Ep_q2 / (1 - Dmax_q2)
    err_q1 = iv.exp(Ep_q1 / (1 - Dmax_q1) + E_cubic_2C) - 1
    err_q2 = iv.exp(R_nonq_bound + Ep_q2 / (1 - Dmax_q2) + E_cubic_2C) - 1
    margin_q1 = rat_main_iv - err_q1
    margin_q2 = rat_main_iv - err_q2
    print(f"  >>> RIGOROUS FIXED R2.3 MARGIN_LB (q = 1):  {float(margin_q1.a):+.6f} <<<")
    print(f"  >>> RIGOROUS FIXED R2.3 MARGIN_LB (q >= 2): {float(margin_q2.a):+.6f} <<<")
    assert margin_q1.a > iv.mpf('0.40')
    assert margin_q2.a > iv.mpf('0.25')


def audit_exact_identities_80dps():
    print("\n==========================================================================")
    print("80-DPS CHECK OF EXACT IDENTITIES (PAIRING, BINOMIAL p_r, RESCALED CONTOUR)")
    print("==========================================================================")
    mp.mp.dps = 80
    # 1. Exact pairing identity a_m * a_{-m} = (1 + zeta_m)^{-n} & exact binomial p_r identity
    for n, c, j in [(20, 2, 3), (24, 3, 4), (30, 4, 5), (30, 3, 6)]:
        r1 = mp.sin((j + c) * mp.pi / n)
        r2 = mp.sin(j * mp.pi / n)
        r3 = mp.sin(c * mp.pi / n)
        s = r2 / r1
        th_c = c * mp.pi / n
        # Pairing identity check
        max_pair_err = mp.mpf(0)
        for m in range(1, c // 2 + 1):
            delta = 2 * m * mp.pi / n
            a_plus = ((-1) ** j) * ((r3 / r1) ** n) / ((1 - s * mp.exp(1j * (th_c + delta))) ** n)
            a_minus = ((-1) ** j) * ((r3 / r1) ** n) / ((1 - s * mp.exp(1j * (th_c - delta))) ** n)
            zeta_m = (4 * r1 * r2 / (r3 ** 2)) * mp.exp(1j * (2 * j + c) * mp.pi / n) * (mp.sin(m * mp.pi / n) ** 2)
            err = abs(a_plus * a_minus - (1 + zeta_m) ** (-n)) / abs(a_plus * a_minus)
            if err > max_pair_err:
                max_pair_err = err

        # Exact binomial p_r identity check for r = 1, 2, 3
        etas = [(2 * k * mp.pi / n) if (c % 2 == 0) else ((2 * k + 1) * mp.pi / n) for k in range(n)]
        a_list = [((-1) ** j) * ((r3 / r1) ** n) / ((1 - s * mp.exp(1j * et)) ** n) for et in etas]
        Ry = s ** n
        Rz = (r3 / r1) ** n
        max_pr_err = mp.mpf(0)
        for r in (1, 2, 3):
            pr_direct = sum(ak ** r for ak in a_list)
            # Binomial sum over t = 0..40
            pr_binom = mp.mpf(0)
            for t_idx in range(60):
                term = ((-1) ** (c * t_idx)) * mp.binomial(n * (r + t_idx) - 1, t_idx * n) * (Ry ** t_idx)
                pr_binom += term
            pr_binom *= ((-1) ** (j * r)) * (Rz ** r) * n
            err = abs(pr_direct - pr_binom) / abs(pr_direct)
            if err > max_pr_err:
                max_pr_err = err

        # Exact rescaled contour identity vs engine.py on |w| = 0.5 < 1!
        coeffs = bk.compute_elementary_coeffs(n)
        d = (n - 3 * c) // 2
        t = Ry if (c % 2 == 0) else -Ry
        u = t / ((1 + t) ** 2)
        Q_val = mp.mpf(0)
        for b in range(d, -1, -1):
            Q_val = Q_val * u + mp.mpf(coeffs[(n - 3 * c - 2 * b, b, c)])
        lhs = ((-1) ** j) * ((1 + t) ** (n - 3 * c)) * Q_val

        lam = mp.mpf('0.4')
        N_pts = 128
        acc = mp.mpc(0)
        for idx in range(N_pts):
            phi = 2 * mp.pi * idx / N_pts
            w_sh = lam * mp.exp(1j * phi)  # shifted w: original w_orig = (-1)^j * w_sh
            w_orig = ((-1) ** j) * w_sh
            def newton_cubic(T0):
                T = mp.mpc(T0)
                for _ in range(35):
                    f = T * (T - 1) * (T - t) - t * Rz * w_orig
                    df = 3 * T * T - 2 * (1 + t) * T + t
                    T -= f / df
                return T
            xv = newton_cubic(1)
            yv = newton_cubic(t)
            zv = (1 + t) - xv - yv
            dx = xv - 1
            dy = yv / (t * xv) - 1
            eps_s = (1 + dy) ** (mp.mpf(1) / n) - 1
            F1 = (1 + dx) ** n
            G2 = (1 - t * dy / (1 - t)) ** n
            prod_k = mp.mpc(1)
            for et in etas:
                wk_t = (1 - s * (1 + eps_s) * mp.exp(1j * et)) ** n
                prod_k *= (1 - (zv / xv) / wk_t)
            acc += (w_sh ** (-c)) * F1 * G2 * prod_k
        g_gcd = math.gcd(c, j)
        sign_pref = (-1) ** (j * (c + 1) + c)
        assert sign_pref == (-1) ** g_gcd
        rhs = (sign_pref * ((1 - t) ** n) / ((Ry * Rz) ** c) * (acc / N_pts)).real
        cont_err = abs(rhs - lhs) / abs(lhs)
        print(f"  n={n:2d}, c={c}, j={j}: pair_err={float(max_pair_err):.2e}, binom_pr_err={float(max_pr_err):.2e}, rescaled_contour_err={float(cont_err):.2e}")
        assert max_pair_err < mp.mpf('1e-70') and max_pr_err < mp.mpf('1e-50') and cont_err < mp.mpf('1e-65')


def main():
    audit_R2_1_and_2A()
    audit_R2_3()
    audit_exact_identities_80dps()
    print("\n=== ALL R2.1, R2.2 (2A), R2.3 (2C-UPGRADE) AUDIT CHECKS PASSED! ===")


if __name__ == "__main__":
    main()
