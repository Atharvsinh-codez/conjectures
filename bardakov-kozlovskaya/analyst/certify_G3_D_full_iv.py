"""
COMPLETE 50-DIGIT mpmath.iv CERTIFICATION OF PIECE G3 CASE G3-D (|D| = 1, tau in [1.05, 40.0])!

100% pure 50-digit mpf interval arithmetic (ZERO float64 truncation inside iv.mpf([...]),
ZERO complex-interval wrapping, 100% explicit finite-n perturbation bounds):
  1. Exact 2-step (1-w)-FREE Li_2/W_2 identity for g_{1,fin}(u_0) (ZERO u=0 log singularity!).
  2. Exact 1-line (1-w)-FREE W_1 identity (2nd power in P(u)!) for valley drop Re(Psi_fin(0) - Psi_fin(th)).
  3. Exact 2-step (1-w)-FREE W_2 identity (4th power in P(u)!) for peak curvature -Re Psi_fin''(th)
     and |Im Psi_fin''(th)| on |th| <= th_core.
  4. Exact per-box geometry: r_0*c*kap^2 <= min(42, 0.45*tau, tau^2/ct_lo) and kap <= min(0.45, tau/c_min),
     giving ultra-tight perturbations on BOTH moderate-ct fixed-u_0 tiers AND large-ct exact-saddle tiers!
  5. Exact Sub-Lemma G3.3 Poisson peak error E_pk(box) (optimal quadratic pole height!) and 2-pole safe
     valley factor V_bump <= 4.10 for ALL r_0 >= 1 and ALL sectors nu = 0..r_0-1!
"""
import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from math import pi as fpi, sin as fsin, cos as fcos, log as flog, exp as fexp, sqrt as fsqrt
from scipy.integrate import quad
import mpmath as mp
from mpmath import iv

t_start = time.time()
iv.dps = 50
mp.mp.dps = 18  # fast precision for Part 0 numerical quad self-tests (restored to 50 before Part 1!)
PI = iv.pi

# ==============================================================================
# PART 0: FAST MACHINE-PRECISION SELF-TEST OF ALL 5 EXACT (1-w)-FREE IDENTITIES
# ==============================================================================
def run_part0_selftests():
    print("==========================================================================")
    print("PART 0: SELF-TEST OF EXACT (1-w)-FREE g1(W1), g1(W2), W1, W2, Li3 (D=+-1)")
    print("==========================================================================")
    mp.mp.dps = 18
    def selftest_identities(c, j, r0, D, target_kap, u0=0.42):
        n_est = int(round(2.0 * fpi * fpi * j * (j + c) / (c * target_kap)))
        n = max(132, 2 * ((n_est + 1) // 2))
        kap = 2.0 * fpi * fsin(j * fpi / n) * fsin((j + c) * fpi / n) / fsin(c * fpi / n)
        tau = c * r0 * kap
        p0 = round(r0 * j / c)
        assert r0 * j - c * p0 == D
        rj = fsin(j * fpi / n) / fsin((j + c) * fpi / n)
        sj = fsin(c * fpi / n) / fsin((j + c) * fpi / n)
        th_c = c * fpi / n
        def zeta(u):
            eta = 2.0 * u * th_c
            return r0 * n * np.log(sj / (1.0 - rj * np.exp(1j * eta))) + 1j * fpi * r0 * j - 2j * fpi * c * p0 * (u + 0.5)
        def dzeta(u):
            w = rj * np.exp(2j * u * th_c)
            return (2j * fpi * c * ((r0 + p0) * w - p0)) / (1.0 - w)
        def d2zeta(u):
            w = rj * np.exp(2j * u * th_c)
            return -(4.0 * fpi * c * r0 * th_c * w) / ((1.0 - w) ** 2)
        def W1_free(u):
            w = rj * np.exp(2j * u * th_c)
            A = 2j * fpi * c * ((r0 + p0) * w - p0)
            B = 4.0 * fpi * c * r0 * th_c
            return (B * w) / (A ** 2)
        def W2_free(u):
            w = rj * np.exp(2j * u * th_c)
            A = 2j * fpi * c * ((r0 + p0) * w - p0)
            B = 4.0 * fpi * c * r0 * th_c
            return (3.0 * (B ** 2) * (w ** 2) + A * B * (2j * th_c) * w * (1.0 + w)) / (A ** 4)
        lam_W = fexp(-zeta(u0).real)
        L_u = min(n / (2.0 * c), 0.5 + fsqrt(35.0 / max(tau, 0.2)))
        def wrap(x): return (x + fpi) % (2.0 * fpi) - fpi
        def Psi_p(th):
            vr, _ = quad(lambda u: (lam_W * np.exp(zeta(u) + 1j * th) / (1.0 + lam_W * np.exp(zeta(u) + 1j * th))).real,
                         -L_u, L_u, points=[-u0, u0], limit=120, epsabs=1e-10, epsrel=1e-10)
            vi, _ = quad(lambda u: (lam_W * np.exp(zeta(u) + 1j * th) / (1.0 + lam_W * np.exp(zeta(u) + 1j * th))).imag,
                         -L_u, L_u, points=[-u0, u0], limit=120, epsabs=1e-10, epsrel=1e-10)
            return 1j * (-1.0 + vr + 1j * vi)
        def L1(u, th):
            z = lam_W * np.exp(zeta(u) + 1j * th)
            return np.log(1.0 + 1.0 / z) if abs(u) < u0 else np.log(1.0 + z)
        def ME(u, th):
            z = lam_W * np.exp(zeta(u) + 1j * th)
            return complex(mp.polylog(2, -1.0 / z)) if abs(u) < u0 else -complex(mp.polylog(2, -z))
        g1_num = (Psi_p(0.0) / 1j).real
        ap0 = wrap(zeta(u0).imag); am0 = wrap(zeta(-u0).imag)
        b1 = (2.0 * u0 - 1.0 - 1j * ap0 / dzeta(u0) + 1j * am0 / dzeta(-u0)).real
        i1_g1 = quad(lambda u: (-L1(u, 0.0) * W1_free(u)).real, -L_u, L_u, points=[-u0, u0], limit=120, epsabs=1e-10)[0]
        err_g1_W1 = abs(g1_num - (b1 + i1_g1))
        b2 = ((0.5 * ap0 * ap0 - fpi * fpi / 6.0) * (d2zeta(u0) / (dzeta(u0) ** 3) - d2zeta(-u0) / (dzeta(-u0) ** 3))).real
        i2_g1 = quad(lambda u: (ME(u, 0.0) * W2_free(u)).real, -L_u, L_u, points=[-u0, u0], limit=120, epsabs=1e-9)[0]
        err_g1_W2 = abs(g1_num - (b1 + b2 + i2_g1))
        def Re_Psi(th):
            def f(u):
                z = lam_W * np.exp(zeta(u) + 1j * th)
                az = abs(z)
                return (flog(az) + flog(max(abs(1.0 + 1.0 / z), 1e-300))) if az > 1.0 else flog(max(abs(1.0 + z), 1e-300))
            return -flog(lam_W) + quad(f, -L_u, L_u, points=[-u0, u0], limit=120, epsabs=1e-9)[0]
        dRe_num = Re_Psi(fpi) - Re_Psi(0.0)
        pts_t = [t for t in [fpi - abs(ap0), fpi - abs(am0)] if 0 < t < fpi]
        b_re = quad(lambda t: (wrap(zeta(u0).imag + t) / dzeta(u0) - wrap(zeta(-u0).imag + t) / dzeta(-u0)).real,
                    0.0, fpi, points=pts_t, limit=80, epsabs=1e-9)[0]
        w_re = quad(lambda u: ((ME(u, fpi) - ME(u, 0.0)) * (-W1_free(u))).real, -L_u, L_u, points=[-u0, u0], limit=120, epsabs=1e-9)[0]
        err_W1 = abs(dRe_num - (b_re + w_re))
        th = 0.8; eps = 1e-5
        P2_num = (Psi_p(th + eps) - Psi_p(th - eps)) / (2.0 * eps)
        ap = wrap(zeta(u0).imag + th); am = wrap(zeta(-u0).imag + th)
        H_term = 1.0 / dzeta(u0) - 1.0 / dzeta(-u0)
        B1_term = 1j * (ap * d2zeta(u0) / (dzeta(u0) ** 3) - am * d2zeta(-u0) / (dzeta(-u0) ** 3))
        i2r, _ = quad(lambda u: (L1(u, th) * W2_free(u)).real, -L_u, L_u, points=[-u0, u0], limit=120, epsabs=1e-9)
        i2i, _ = quad(lambda u: (L1(u, th) * W2_free(u)).imag, -L_u, L_u, points=[-u0, u0], limit=120, epsabs=1e-9)
        err_W2 = abs(P2_num - (H_term + B1_term - (i2r + 1j * i2i)))
        assert err_g1_W1 < 1e-8 and err_g1_W2 < 1e-8 and err_W1 < 1e-8 and err_W2 < 1e-6, (err_g1_W1, err_g1_W2, err_W1, err_W2)
        print(f"  [PASS] c={c:2d} j={j:2d} r0={r0} D={D:+d} tau={tau:5.2f} | err_g1(W1)={err_g1_W1:.1e} err_g1(W2)={err_g1_W2:.1e} | err_W1@pi={err_W1:.1e} | err_W2={err_W2:.1e}")
    selftest_identities(20, 21, 1, +1, 0.20, 0.42)
    selftest_identities(20, 39, 1, -1, 0.20, 0.42)
    selftest_identities(46, 61, 3, -1, 0.05, 0.45)
    selftest_identities(31, 47, 2, +1, 0.12, 0.46)
    _check_par_li3_id()
    mp.mp.dps = 50
    print("  [PASS] Exact parabola Li_3 / W_2 identity verified to < 1e-8 (th = 0.7, 2.3, pi)!!")

# ==============================================================================
# PART 1: EXACT G3.1 CONSTANTS & PURE 50-DIGIT iv PRIMITIVES (NO float() CASTS!)
# ==============================================================================
CT_MIN    = iv.mpf(10)
R0K_MAX   = iv.mpf('0.45')
CTAU2_MIN = iv.mpf(40)
CK2_MAX   = iv.mpf(42)
R3CK2_MAX = iv.mpf(60)

def box_geom(kap_ub):
    """Exact 50-digit iv geometry for any certified upper bound kap <= kap_ub <= 0.45."""
    # Since sin(x) >= x - x^3/6 for x in (0, 0.10), al_ub = (kap_ub/(2*pi)) / (1 - (kap_ub/(2*pi))^2/6):
    y = kap_ub / (2 * PI)
    al_b = iv.mpf((y / (1 - (y ** 2) / 5)).b)
    assert (iv.sin(al_b) - y).a > 0
    th_m = al_b / 2
    S_b = iv.sin(al_b) / al_b
    c_sh = iv.mpf(1) / (24 * (1 - (th_m ** 2) / 24))
    c_cos = iv.mpf(1) / (2 * (1 - (al_b ** 2) / 2))
    C_D0 = (c_cos - (iv.mpf(1) / 6) * (1 - (al_b ** 2) / 20) + c_cos + c_cos / 4 + c_sh + c_sh * c_cos * (al_b ** 2)) / (16 * (PI ** 2) * (S_b ** 4))
    gam_m = (kap_ub ** 2) / (2 * (PI ** 2) * (S_b ** 2) * (iv.cos(al_b) ** 3))
    xcot = (th_m / 2) * iv.cos(th_m / 2) / iv.sin(th_m / 2)
    te_lo = xcot
    te_hi = 1 / (iv.cos(al_b) ** 2)
    assert (1 - te_lo).b <= (te_hi - 1).a
    T_c5 = iv.mpf('0.05'); ut_c5 = iv.mpf('0.080')
    assert (iv.sin(ut_c5) - iv.sqrt(T_c5 / 8)).a > 0
    sinc_c5_sq = (iv.sin(ut_c5) / ut_c5) ** 2
    sinc_2c5 = iv.sin(2 * ut_c5) / (2 * ut_c5)
    c_re5 = te_lo * sinc_2c5 / (1 + T_c5)
    T_in = gam_m * iv.mpf('0.31'); ut_in = iv.mpf('0.55') * th_m
    assert T_in.b < T_c5.a
    sinc_in_sq = (iv.sin(ut_in) / ut_in) ** 2
    sinc_2in = iv.sin(2 * ut_in) / (2 * ut_in)
    c_rein = te_lo * sinc_2in / (1 + T_in)
    r_re0 = te_lo * (iv.sin(2 * iv.mpf('0.52') * th_m) / (2 * iv.mpf('0.52') * th_m)) / (1 + gam_m * iv.mpf('0.28'))
    S_half = iv.sin(iv.mpf('0.5') * th_m) / (iv.mpf('0.5') * th_m)
    K_D2 = 12 / ((S_half ** 3) * (iv.cos(al_b) ** 3) * (1 - (al_b ** 2) / 3))
    C_D2 = K_D2 * C_D0
    thc_hi = kap_ub / (4 * PI * (S_b ** 2))
    thc_1ms = kap_ub / (2 * PI * S_half * iv.cos(al_b))
    n_m = 2 * (PI ** 2) * (S_b ** 2) * 23 / kap_ub
    RzRy = iv.exp(-PI * iv.cos(al_b / 2) / iv.sin(al_b / 2))
    y_br = PI * iv.cos(al_b) / kap_ub
    r_br = iv.mpf('0.50') / y_br
    q_str = 1 - r_br
    te_str = te_hi * (iv.mpf(2) ** r_br) / (q_str ** 2)
    ut_50 = iv.mpf('0.255')
    assert (iv.sin(ut_50) - iv.mpf('0.501') / iv.sqrt(iv.mpf(4))).a > 0
    c_re50 = te_lo * ((iv.log(iv.mpf('1.5')) - iv.log(1 + T_in)) / (iv.mpf('0.5') - T_in)) * (iv.sin(2 * ut_50) / (2 * ut_50))
    return dict(al=al_b, th_m=th_m, S=S_b, S_half=S_half, C_D0=C_D0, K_D2=K_D2, C_D2=C_D2, gam=gam_m,
                te_lo=te_lo, te_hi=te_hi, te_str=te_str, T_c5=T_c5, sinc_c5_sq=sinc_c5_sq,
                c_re5=c_re5, c_re50=c_re50, c_rein=c_rein, r_re0=r_re0, T_in=T_in, sinc_in_sq=sinc_in_sq,
                thc_hi=thc_hi, thc_1ms=thc_1ms, n_min=n_m, RzRy=RzRy, y_br=y_br)

def iv_min(x, y): return iv.mpf(min(x.a, y.a))
def iv_max_hi(x, y): return iv.mpf(max(x.b, y.b))
def iv_min_hi(x, y): return iv.mpf(min(x.b, y.b))
def iv_abs_hi(x): return iv.mpf(max(abs(x.a), abs(x.b)))
def iv_sym(e): return iv.mpf([-e.b, e.b])

# Rigorous 50-digit iv enclosure of zeta(3) via partial sum M=2000 + monotone integral sandwich (ZERO hardcoded 1.21!):
_M_z3 = 2000
_z3_part = sum(iv.mpf(1) / (iv.mpf(m) ** 3) for m in range(1, _M_z3 + 1))
zeta3_D = _z3_part + iv.mpf([(iv.mpf(1) / (2 * (iv.mpf(_M_z3 + 1) ** 2))).a, (iv.mpf(1) / (2 * (iv.mpf(_M_z3) ** 2))).b])

def _F_mlog(x):
    """Exact antiderivative F(x) = (1+x)*log(1+x) - x*log(x) of log((1+x)/x) on x >= 0 (F(0) = 0)."""
    if x.b <= 0:
        return iv.mpf(0)
    return (1 + x) * iv.log(1 + x) - x * iv.log(x)

def int_mlog_ub(alpha, d1_lo, d2_hi, x_f_lo, r_f_ub, du):
    """Exact 50-digit iv upper bound on int_{ua}^{ub} (-log(1 - r_f(u))) du (#106 fix #3, ZERO hardcoded 14.0!)."""
    a_lo = iv.mpf(alpha.a)
    assert a_lo.a > 0, f"int_mlog_ub alpha <= 0: {alpha}"
    x1 = a_lo * iv.mpf(max(mp.mpf(0), d1_lo.a))
    x2 = a_lo * iv.mpf(d2_hi.b)
    res = (_F_mlog(x2) - _F_mlog(x1)) / a_lo
    assert res.b > 0
    if x_f_lo.a > mp.mpf('1e-6'):
        res = iv_min_hi(res, (-iv.log(1 - r_f_ub)) * du)
    return iv.mpf(res.b)

# Zero-hardcoded universal G3 Pigeonhole-Annulus E_cubic_G3_sub + shifted-line NB_err / Y_big_err:
G0_D = box_geom(R0K_MAX)
N_G3_D = iv.mpf(int(G0_D['n_min'].a))
assert N_G3_D.a >= 1000
K_jc_D = (CK2_MAX / (2 * (PI ** 4) * (G0_D['S'] ** 4))) ** (iv.mpf(1) / 3)
xcot_D = (G0_D['al'] / 2) * iv.cos(G0_D['al'] / 2) / iv.sin(G0_D['al'] / 2)
c_se_D = 2 * xcot_D / K_jc_D
c_se_D_eff = c_se_D * ((N_G3_D - 1) / N_G3_D)
N_sp_D = iv.mpf(50000)
assert ((12 / c_se_D_eff) ** 3).b < N_sp_D.a
M4_D = iv_max_hi((N_sp_D ** 4) * (G0_D['RzRy'] ** ((N_G3_D - 1) / N_G3_D)),
                 (N_sp_D ** 4) * iv.exp(-c_se_D_eff * (N_sp_D ** (iv.mpf(1) / 3))))
M3_D = iv_max_hi((N_sp_D ** 3) * (G0_D['RzRy'] ** ((N_G3_D - 1) / N_G3_D)),
                 (N_sp_D ** 3) * iv.exp(-c_se_D_eff * (N_sp_D ** (iv.mpf(1) / 3))))
eps_ph_D = iv.mpf('1e-5')
lam_max_D = iv.mpf('3.01'); lam_inv_D = iv.mpf('20.1')
B_pure_den_D = (1 / lam_max_D) * ((eps_ph_D / (1 + eps_ph_D / N_G3_D)) ** 2)
E_cubic_D = (-iv.log(1 - (25 * lam_inv_D * lam_max_D * M4_D) / B_pure_den_D)
             - iv.log(1 - (14 * lam_max_D / (iv.mpf('0.49') * eps_ph_D)) * M3_D)
             + 20 * lam_max_D * M3_D)
L1_shift_D = (1 - eps_ph_D / N_G3_D) ** (-N_G3_D)
E_cubic_G3_sub = (iv.exp(E_cubic_D) - 1) * L1_shift_D
assert E_cubic_G3_sub.b < mp.mpf('1e-5')

r02k_max_1 = iv.sqrt(R3CK2_MAX / CT_MIN)
NB_exp_D = (PI ** 2) * G0_D['c_rein'] * CT_MIN / r02k_max_1
NB_err_D = 2 * iv.exp(-NB_exp_D) / (1 - iv.exp(-3 * NB_exp_D))
Y_big_err_D = 100 * iv.exp(-9 * PI) / (1 - iv.exp(-9 * PI))
assert NB_exp_D.a > 39 and NB_err_D.b < mp.mpf('1e-16') and Y_big_err_D.b < mp.mpf('1e-10')
assert (PI * (1 + G0_D['C_D0'] * CK2_MAX) / N_G3_D + (iv.mpf('0.25') * G0_D['te_str'] / PI) * (G0_D['al'] / 2) * r02k_max_1).b < (iv.mpf('0.5') * iv.log(iv.mpf('1.5') / (1 + G0_D['T_in']))).a

def solve_out_neg_E_D1(G, u0_lo, u0_hi, Del0_hi, ct_over_tau_min, ct_tau_min, ct_min, Delta_iv):
    d_hi = 1 + Del0_hi
    c1 = d_hi / (G['c_re50'] * ct_tau_min)
    c2 = G['te_str'] / (4 * (PI ** 2) * G['c_re50'] * (ct_min ** 2))
    b1 = G['te_hi'] / (PI * ct_over_tau_min)
    b2 = G['te_str'] / (4 * (PI ** 2) * ct_min * ct_over_tau_min)
    def phi_ub(E):
        rhs = c1 * E + c2 * (E ** 2)
        d = rhs / (u0_lo + iv.sqrt(u0_lo ** 2 + rhs))
        return 2 * PI * d_hi * d + b1 * (u0_hi + d) * E + b2 * (E ** 2)
    lo = mp.mpf(0); hi = mp.mpf('1.0')
    while phi_ub(iv.mpf(hi)).b < Delta_iv.a:
        hi *= 2
    for _ in range(30):
        mid = (lo + hi) / 2
        if phi_ub(iv.mpf(mid)).b < Delta_iv.a:
            lo = mid
        else:
            hi = mid
    E_cert = iv.mpf(lo)
    assert E_cert.a > 0 and phi_ub(E_cert).b < Delta_iv.a
    return E_cert

def cert_E_higher_D1(G, u0_lo, u0_hi, Del0_hi, ct_over_tau_min, ct_tau_min, ct_min):
    d_hi = 1 + Del0_hi; d_lo = 1 - Del0_hi
    e_in_step = ((PI ** 2) / (G['te_hi'] * u0_hi)) * ct_over_tau_min
    assert e_in_step.a > mp.mpf('5.0'), f"e_in_step={e_in_step}"
    S_in_a = 2 * (-iv.log(1 - iv.exp(-e_in_step))) / (1 - iv.exp(-2 * e_in_step))
    S_in_b = iv.mpf(0)
    if u0_hi.b > (1 / (2 * d_hi)).a:
        d_m1 = 1 / (2 * d_hi)
        e_deep = G['c_rein'] * ct_tau_min * (2 * u0_lo * d_m1 - d_m1 ** 2) / (2 * d_hi)
        assert e_deep.a > mp.mpf('4.5'), f"e_deep={e_deep}"
        S_in_b = 2 * (-iv.log(1 - iv.exp(-e_deep))) + 2 * (-iv.log(1 - iv.exp(-2 * e_in_step))) / (1 - iv.exp(-2 * e_in_step))
    # y > 0 outer branch (k >= 2 on BOTH sides, #106 fix #5b!):
    d_lo_u0 = 1 - iv_max_hi(Del0_hi, (G['K_D2'] * (u0_hi ** 2) - 1) * Del0_hi)
    assert d_lo_u0.a > 0, f"d_lo_u0 <= 0: {d_lo_u0}"
    C_x = iv.sqrt(G['te_str'] / G['c_re50'])
    C_quad = (3 + d_hi / d_lo_u0) * G['te_str'] * C_x
    E_out_0 = (4 * (PI ** 2) * d_lo_u0 / G['te_str']) * ct_over_tau_min
    c_out = 2 * PI * iv.sqrt((2 * PI * ct_min * ct_over_tau_min) / C_quad)
    assert E_out_0.a > mp.mpf('8.0') and c_out.a > mp.mpf('8.0')
    S_out_pos = iv.mpf(0); K_cut = 16
    for kk in range(2, K_cut + 1):
        ek = iv.mpf(max(E_out_0.a, (c_out * iv.sqrt(kk - 1)).a))
        S_out_pos += 2 * (-iv.log(1 - iv.exp(-ek)))
    u_t = c_out * iv.sqrt(K_cut - 1)
    S_out_pos += 2 * iv.mpf('1.01') * (2 / (c_out ** 2)) * (u_t + 1) * iv.exp(-u_t)
    # y < 0 outer branch (k >= 2 on BOTH sides! #77 fix 2a!):
    S_out_neg = iv.mpf(0); E_16 = None
    for kk in range(2, K_cut + 1):
        ek = solve_out_neg_E_D1(G, u0_lo, u0_hi, Del0_hi, ct_over_tau_min, ct_tau_min, ct_min, 2 * PI * (kk - 1))
        S_out_neg += 2 * (-iv.log(1 - iv.exp(-ek)))
        if kk == K_cut: E_16 = ek
    c1 = d_hi / (G['c_re50'] * ct_tau_min)
    c2 = G['te_str'] / (4 * (PI ** 2) * G['c_re50'] * (ct_min ** 2))
    b1 = G['te_hi'] / (PI * ct_over_tau_min)
    b2 = G['te_str'] / (4 * (PI ** 2) * ct_min * ct_over_tau_min)
    s_d = iv.sqrt(c1 / E_16 + c2)
    A_tail = (2 * PI * d_hi * s_d + b1 * u0_hi) / E_16 + b1 * s_d + b2
    c_neg_tail = iv.sqrt(2 * PI / A_tail)
    assert c_neg_tail.a > mp.mpf('3.0')
    u_nt = c_neg_tail * iv.sqrt(K_cut - 1)
    S_out_neg += 2 * iv.mpf('1.01') * (2 / (c_neg_tail ** 2)) * (u_nt + 1) * iv.exp(-u_nt)
    return S_in_a + S_in_b + S_out_pos + S_out_neg

def inv_pow_re_im(a, b, n):
    r2 = a ** 2 + b ** 2
    if n == 1: return a / r2, -b / r2
    if n == 2:
        r4 = r2 ** 2
        return (a ** 2 - b ** 2) / r4, (-2 * a * b) / r4
    if n == 3:
        r6 = r2 ** 3
        return (a ** 3 - 3 * a * (b ** 2)) / r6, (-3 * (a ** 2) * b + b ** 3) / r6
    if n == 4:
        r8 = r2 ** 4
        return (a ** 4 - 6 * (a ** 2) * (b ** 2) + b ** 4) / r8, (-4 * (a ** 3) * b + 4 * a * (b ** 3)) / r8
    raise ValueError(n)

def li2_pos_ub(r_ub):
    """Exact upper bound on Li_2(r) = sum_{m>=1} r^m/m^2 for r in [0, 1]."""
    return r_ub + (r_ub ** 2) / 4 + (r_ub ** 3) / 9 + ((PI ** 2) / 6 - iv.mpf(49) / 36) * (r_ub ** 4)

def w1_par_mv(tau_lo, tau_hi, ua, ub):
    """Zero-dependency centered mean-value enclosure of w1(tau, u) = (tau/2)/(tau*u - i*pi)^2."""
    tm = iv.mpf(((tau_lo + tau_hi) / 2).a); um = iv.mpf(((ua + ub) / 2).a)
    dt = iv_max_hi(tau_hi - tm, tm - tau_lo); du = iv_max_hi(ub - um, um - ua)
    r0, i0 = inv_pow_re_im(tm * um, -PI, 2)
    r0 = (iv.mpf('0.5') * tm) * r0; i0 = (iv.mpf('0.5') * tm) * i0
    # |d/dtau w1| = 1 / (2*(tau^2*u^2 + pi^2)), |d/du w1| = tau^2 / (tau^2*u^2 + pi^2)^(3/2)
    dt_ub = iv.mpf('0.5') / ((tau_lo * ua) ** 2 + PI ** 2)
    # Max of t^2 / (t^2*ua^2 + pi^2)^(3/2) on [tau_lo, tau_hi] (unique peak at t* = sqrt(2)*pi/ua):
    du_ub = iv_max_hi((tau_lo ** 2) / (((tau_lo * ua) ** 2 + PI ** 2) ** iv.mpf('1.5')),
                      (tau_hi ** 2) / (((tau_hi * ua) ** 2 + PI ** 2) ** iv.mpf('1.5')))
    if ua.a > 0:
        t_st = iv.sqrt(2) * PI / ua
        if tau_lo.a <= t_st.b and t_st.a <= tau_hi.b:
            du_ub = iv_max_hi(du_ub, 2 / (3 * iv.sqrt(3) * PI * (ua ** 2)))
    rad = dt_ub * dt + du_ub * du
    return r0 + iv_sym(rad), i0 + iv_sym(rad)

def w2_par_mv(tau_lo, tau_hi, ua, ub):
    """Zero-dependency centered mean-value enclosure of w2(tau, u) = (3/4)*tau^2/(tau*u - i*pi)^4."""
    tm = iv.mpf(((tau_lo + tau_hi) / 2).a); um = iv.mpf(((ua + ub) / 2).a)
    dt = iv_max_hi(tau_hi - tm, tm - tau_lo); du = iv_max_hi(ub - um, um - ua)
    r0, i0 = inv_pow_re_im(tm * um, -PI, 4)
    r0 = (iv.mpf('0.75') * (tm ** 2)) * r0; i0 = (iv.mpf('0.75') * (tm ** 2)) * i0
    # |d/dtau w2| = (3/2)*tau / (tau^2*u^2 + pi^2)^2 (unique peak at t* = pi/(sqrt(3)*ua)):
    dt_ub = iv.mpf('1.5') * iv_max_hi(tau_lo / (((tau_lo * ua) ** 2 + PI ** 2) ** 2),
                                      tau_hi / (((tau_hi * ua) ** 2 + PI ** 2) ** 2))
    if ua.a > 0:
        t_st = PI / (iv.sqrt(3) * ua)
        if tau_lo.a <= t_st.b and t_st.a <= tau_hi.b:
            dt_ub = iv_max_hi(dt_ub, iv.mpf('1.5') * (3 * iv.sqrt(3)) / (16 * (PI ** 3) * ua))
    # |d/du w2| = 3*tau^3 / (tau^2*u^2 + pi^2)^(5/2) (unique peak at t* = sqrt(3/2)*pi/ua):
    du_ub = 3 * iv_max_hi((tau_lo ** 3) / (((tau_lo * ua) ** 2 + PI ** 2) ** iv.mpf('2.5')),
                          (tau_hi ** 3) / (((tau_hi * ua) ** 2 + PI ** 2) ** iv.mpf('2.5')))
    if ua.a > 0:
        t_st2 = iv.sqrt(iv.mpf('1.5')) * PI / ua
        if tau_lo.a <= t_st2.b and t_st2.a <= tau_hi.b:
            du_ub = iv_max_hi(du_ub, 3 * (iv.mpf('1.5') ** iv.mpf('1.5')) / ((iv.mpf('2.5') ** iv.mpf('2.5')) * (PI ** 2) * (ua ** 3)))
    rad = dt_ub * dt + du_ub * du
    return r0 + iv_sym(rad), i0 + iv_sym(rad)

def tau_c3_mv(tau_lo, tau_hi, u0_lo, u0_hi):
    """Zero-dependency centered mean-value enclosure of tau / (tau*u0 - i*pi)^3."""
    tm = iv.mpf(((tau_lo + tau_hi) / 2).a); um = iv.mpf(((u0_lo + u0_hi) / 2).a)
    dt = iv_max_hi(tau_hi - tm, tm - tau_lo); du = iv_max_hi(u0_hi - um, um - u0_lo)
    r0, i0 = inv_pow_re_im(tm * um, -PI, 3)
    r0 = tm * r0; i0 = tm * i0
    # |d/dtau (tau/(tau*u - i*pi)^3)| = sqrt(4*tau^2*u^2 + pi^2) / (tau^2*u^2 + pi^2)^2 <= 2 / (tau^2*u^2 + pi^2)^(3/2)
    dt_ub = 2 / (((tau_lo * u0_lo) ** 2 + PI ** 2) ** iv.mpf('1.5'))
    # |d/du (tau/(tau*u - i*pi)^3)| = 3*tau^2 / (tau^2*u^2 + pi^2)^2 (unique peak at t* = pi/u0_lo):
    du_ub = 3 * iv_max_hi((tau_lo ** 2) / (((tau_lo * u0_lo) ** 2 + PI ** 2) ** 2),
                          (tau_hi ** 2) / (((tau_hi * u0_lo) ** 2 + PI ** 2) ** 2))
    if u0_lo.a > 0:
        t_st = PI / u0_lo
        if tau_lo.a <= t_st.b and t_st.a <= tau_hi.b:
            du_ub = iv_max_hi(du_ub, 3 / (4 * (PI ** 2) * (u0_lo ** 2)))
    rad = dt_ub * dt + du_ub * du
    return r0 + iv_sym(rad), i0 + iv_sym(rad)

def cmul(ar, ai, br, bi): return ar * br - ai * bi, ar * bi + ai * br

def li2_neg_polar(r_iv, ph_iv, r_ub):
    ph_w = float(ph_iv.b) - float(ph_iv.a)
    rv = float(r_ub.b)
    M_max = 10 if rv < 0.45 else (18 if rv < 0.72 else (30 if rv < 0.88 else 45))
    M_use = min(M_max, max(4, int(1.6 / max(ph_w, 1e-6))))
    sr = iv.mpf(0); si = iv.mpf(0); rp = r_iv
    for m in range(1, M_use + 1):
        c = ((-1 if (m % 2 == 1) else 1) * rp) / (m * m)
        sr += c * iv.cos(m * ph_iv); si += c * iv.sin(m * ph_iv)
        rp = rp * r_iv
    t_int = (r_ub ** (M_use + 1)) / M_use
    if rv < 0.99:
        t_geo = (r_ub ** (M_use + 1)) / (((M_use + 1) ** 2) * (1 - r_ub))
        tb = iv_min_hi(t_int, t_geo)
    else:
        tb = t_int
    return sr + iv_sym(tb), si + iv_sym(tb)

def li3_even_cub_polar(r_iv, ph0_iv, th_iv, r_ub, sgn_th, need_cub=True):
    """
    Exact factored polar series for:
      EVEN: 0.5*(Li_3(-r e^{i(ph0+th)}) + Li_3(-r e^{i(ph0-th)})) - Li_3(-r e^{i ph0})
            = sum_{m>=1} ((-1)^m r^m / m^3) * (-2 * sin^2(m*th/2)) * e^{i m ph0}    [O(th^2)!]
      CUB:  0.5*(Li_3(-r e^{i(ph0+s*th)}) - Li_3(-r e^{i(ph0-s*th)})) - i*s*th*Li_2(-r e^{i ph0})
            = sum_{m>=1} ((-1)^m r^m / m^3) * (i * s * (sin(m*th) - m*th)) * e^{i m ph0}  [O(th^3)!]
    """
    ph_w = float(ph0_iv.b) - float(ph0_iv.a)
    th_w = float(th_iv.b) - float(th_iv.a)
    rv = float(r_ub.b)
    if need_cub:
        M_max = 8 if rv < 0.45 else (14 if rv < 0.75 else 24)
    else:
        M_max = 6 if rv < 0.45 else (10 if rv < 0.75 else 18)
    M_use = min(M_max, max(5, int(1.8 / max(ph_w + 0.5 * th_w, 1e-6))))
    er = iv.mpf(0); ei = iv.mpf(0); cr = iv.mpf(0); ci = iv.mpf(0)
    rp = r_iv
    for m in range(1, M_use + 1):
        base = ((-1 if (m % 2 == 1) else 1) * rp) / (m * m * m)
        cm = iv.cos(m * ph0_iv); sm = iv.sin(m * ph0_iv)
        fe = -2 * (iv.sin(m * th_iv / 2) ** 2)
        er += (base * fe) * cm; ei += (base * fe) * sm
        if need_cub:
            fc = sgn_th * (iv.sin(m * th_iv) - m * th_iv)
            cr += (base * fc) * (-sm); ci += (base * fc) * cm
        rp = rp * r_iv
    t_int = (r_ub ** (M_use + 1)) / (2 * (M_use ** 2))
    if rv < 0.99:
        t_geo = (r_ub ** (M_use + 1)) / (((M_use + 1) ** 3) * (1 - r_ub))
        t_base = iv_min_hi(t_int, t_geo)
    else:
        t_base = t_int
    t_even = 2 * t_base
    if rv < 0.99:
        t_even = iv_min_hi(t_even, iv.mpf('0.5') * (iv.mpf(th_iv.b) ** 2) * (r_ub ** (M_use + 1)) / ((M_use + 1) * (1 - r_ub)))
    if not need_cub:
        return er + iv_sym(t_even), ei + iv_sym(t_even), iv.mpf(0), iv.mpf(0)
    t_cub = t_base + iv.mpf(th_iv.b) * ((r_ub ** (M_use + 1)) / M_use)
    if rv < 0.99:
        t_cub = iv_min_hi(t_cub, ((iv.mpf(th_iv.b) ** 3) / 6) * (r_ub ** (M_use + 1)) / (1 - r_ub))
    return er + iv_sym(t_even), ei + iv_sym(t_even), cr + iv_sym(t_cub), ci + iv_sym(t_cub)

# Exact 1-line verification of parabola Li_3 / W_2 identity for -int (ME(th) - ME(0)) W_1 du:
def _check_par_li3_id():
    tau_t = 3.0; u0_t = 0.35; del_t = 0.5 - u0_t
    def wrap_t(x): return (x + fpi) % (2.0 * fpi) - fpi
    for th_t in [0.7, 2.3, fpi]:
        def ME_p(u, th):
            x = tau_t * (u0_t * u0_t - u * u)
            z = np.exp(x + 1j * (fpi + 2 * fpi * u + th))
            return complex(mp.polylog(2, -1.0 / z)) if abs(u) < u0_t else -complex(mp.polylog(2, -z))
        def NE_p(u, th):
            x = tau_t * (u0_t * u0_t - u * u)
            z = np.exp(x + 1j * (fpi + 2 * fpi * u + th))
            return complex(mp.polylog(3, -1.0 / z)) if abs(u) < u0_t else complex(mp.polylog(3, -z))
        def W1_p(u): return (tau_t / 2.0) / ((tau_t * u - 1j * fpi) ** 2)
        def W2_p(u): return (0.75 * tau_t * tau_t) / ((tau_t * u - 1j * fpi) ** 4)
        L_t = 3.8
        lhs_r = quad(lambda u: (-(ME_p(u, th_t) - ME_p(u, 0.0)) * W1_p(u)).real, -L_t, L_t, points=[-u0_t, u0_t], limit=200)[0]
        lhs_i = quad(lambda u: (-(ME_p(u, th_t) - ME_p(u, 0.0)) * W1_p(u)).imag, -L_t, L_t, points=[-u0_t, u0_t], limit=200)[0]
        def bdy_n(th):
            ap = wrap_t(-2 * fpi * del_t + th); am = wrap_t(2 * fpi * del_t + th)
            return (-1j * tau_t / 24.0) * (ap * (fpi**2 - ap**2) / ((tau_t * u0_t - 1j * fpi) ** 3) - am * (fpi**2 - am**2) / ((-tau_t * u0_t - 1j * fpi) ** 3))
        i_r = quad(lambda u: (-(NE_p(u, th_t) - NE_p(u, 0.0)) * W2_p(u)).real, -L_t, L_t, points=[-u0_t, u0_t], limit=200)[0]
        i_i = quad(lambda u: (-(NE_p(u, th_t) - NE_p(u, 0.0)) * W2_p(u)).imag, -L_t, L_t, points=[-u0_t, u0_t], limit=200)[0]
        rhs = (bdy_n(th_t) - bdy_n(0.0)) + (i_r + 1j * i_i)
        assert abs((lhs_r + 1j * lhs_i) - rhs) < 1e-8, (th_t, lhs_r + 1j * lhs_i, rhs)
mp.mp.dps = 50
def log1p_polar(r_iv, ph_iv, r_ub):
    rv = float(r_ub.b)
    if rv < 0.05:
        M = 4
        sr = iv.mpf(0); si = iv.mpf(0); rp = r_iv
        for m in range(1, M + 1):
            c = ((1 if (m % 2 == 1) else -1) * rp) / m
            sr += c * iv.cos(m * ph_iv); si += c * iv.sin(m * ph_iv)
            rp = rp * r_iv
        tb = (r_ub ** (M + 1)) / ((M + 1) * (1 - r_ub))
        return sr + iv_sym(tb), si + iv_sym(tb)
    c_p = iv.cos(ph_iv); s_p = iv.sin(ph_iv)
    xr = 1 + r_iv * c_p; yi = r_iv * s_p
    m2_a = xr ** 2 + yi ** 2
    m2_b = 1 + r_iv * (r_iv + 2 * c_p)
    m2 = iv.mpf([max(m2_a.a, m2_b.a), min(m2_a.b, m2_b.b)])
    assert m2.a > mp.mpf('1e-6'), f"Near -1 singularity: m2={m2}"
    re_part = iv.mpf('0.5') * iv.log(m2)
    if xr.a > 0 or yi.a > 0 or yi.b < 0:
        im_part = iv.atan2(yi, xr)
    else:
        im_part = iv_sym(PI)
    # Zero-dependency centered mean-value intersection (for r <= 1, Re(1+r*e^{i*ph}) > 0 on entire box since m2.a > 0!):
    if r_iv.b <= 1:
        r_a = iv.mpf(r_iv.a); r_b = iv.mpf(r_iv.b)
        p_a = iv.mpf(ph_iv.a); p_b = iv.mpf(ph_iv.b)
        rm = iv.mpf(((r_a + r_b) / 2).a); pm = iv.mpf(((p_a + p_b) / 2).a)
        dr = iv_max_hi(r_b - rm, rm - r_a); dp = iv_max_hi(p_b - pm, pm - p_a)
        xm = 1 + rm * iv.cos(pm); ym = rm * iv.sin(pm)
        re_m = iv.mpf('0.5') * iv.log(xm ** 2 + ym ** 2)
        im_m = iv.atan2(ym, xm)
        rad = (dr + r_b * dp) / iv.sqrt(iv.mpf(m2.a))
        re_mv = re_m + iv_sym(rad); im_mv = im_m + iv_sym(rad)
        re_part = iv.mpf([max(re_part.a, re_mv.a), min(re_part.b, re_mv.b)])
        im_part = iv.mpf([max(im_part.a, im_mv.a), min(im_part.b, im_mv.b)])
        assert re_part.a <= re_part.b and im_part.a <= im_part.b
    return re_part, im_part

def inv_quad_ub(tau_lo, tau_hi, au_lo, pwr):
    u2 = au_lo ** 2
    v1 = 2 * tau_lo * u2 + 2 * (PI ** 2) / tau_lo
    v2 = 2 * tau_hi * u2 + 2 * (PI ** 2) / tau_hi
    vmin = iv_min(v1, v2)
    if au_lo.b > 0:
        t_star = PI / au_lo
        if tau_lo.a <= t_star.b and t_star.a <= tau_hi.b:
            vmin = iv_min(vmin, 4 * PI * au_lo)
    return 1 / (vmin ** pwr)

def eps_P_bound(G, tau_lo, tau_hi, au_lo, au_hi, c_min_iv, Del0_hi):
    rat_re = (2 * tau_hi * au_hi) / iv.sqrt((2 * tau_hi * au_hi) ** 2 + (2 * PI) ** 2)
    b_ub = (1 + Del0_hi) / (2 * c_min_iv)
    dP_re_rel = ((G['te_hi'] - 1) + G['te_hi'] * b_ub + G['te_hi'] * (1 + b_ub) * ((2 * au_hi * G['thc_hi']) ** 2) / 6) * rat_re
    P_p_lo = iv.sqrt((2 * tau_lo * au_lo) ** 2 + (2 * PI) ** 2)
    im_t1 = 2 * PI * Del0_hi / P_p_lo
    im_t2 = G['te_hi'] * (1 + b_ub) * au_hi * G['thc_hi'] * rat_re
    dP_im_rel = iv_max_hi(im_t1, im_t2)
    ep = iv.sqrt(dP_re_rel ** 2 + dP_im_rel ** 2)
    assert ep.b < mp.mpf('0.30'), f"eps_P too large: {float(ep.b)}"
    return ep

def W1_W2_fin_tail_ub(G, tau_lo, tau_hi, L_iv, c_min_iv, Del0_hi):
    """Exact uniform upper bounds on |W1_fin|, |W2_fin| for u in [L_iv, u_0.05] AND on all of [-L_n, L_n]
    (for T(u) > 0.05), PLUS r_T5 & far_r5 >= 2*L_n*sup_{T(u)>0.05} r_f(u) (#106 fix #5a!)."""
    assert (G['gam'] * (L_iv ** 2)).b <= mp.mpf('0.05'), f"gam*L^2 > 0.05: {G['gam']*L_iv**2}"
    assert (L_iv * G['thc_hi']).b <= mp.mpf('0.080'), f"L*thc_hi > 0.080: {L_iv*G['thc_hi']}"
    d_lo = 1 - Del0_hi
    f_p_lo = (1 - 1 / c_min_iv) * (1 - 1 / (2 * c_min_iv))
    mob_out = f_p_lo * (G['S'] ** 4) * G['te_lo'] * G['sinc_c5_sq']
    u2 = L_iv ** 2
    v1 = 2 * mob_out * tau_lo * u2 + 2 * ((d_lo * PI) ** 2) / tau_lo
    v2 = 2 * mob_out * tau_hi * u2 + 2 * ((d_lo * PI) ** 2) / tau_hi
    vmin_out = iv_min(v1, v2)
    t_st = (d_lo * PI) / (iv.sqrt(mob_out) * L_iv)
    if tau_lo.a <= t_st.b and t_st.a <= tau_hi.b:
        vmin_out = iv_min(vmin_out, 4 * iv.sqrt(mob_out) * d_lo * PI * L_iv)
    W1_fin_out = G['te_hi'] / vmin_out
    P_lo_out = 2 * iv.sqrt(mob_out * ((tau_lo * L_iv) ** 2) + (d_lo * PI) ** 2)
    W2_fin_out = 3 * ((G['te_hi'] / vmin_out) ** 2) + 4 * G['te_hi'] * G['thc_1ms'] / (vmin_out * P_lo_out)
    P_lo_all = 2 * PI * d_lo
    W1_fin_all = 2 * G['te_hi'] * tau_hi / (P_lo_all ** 2)
    W2_fin_all = 3 * ((2 * G['te_hi'] * tau_hi) ** 2) / (P_lo_all ** 4) + 8 * G['te_hi'] * tau_hi * G['thc_1ms'] / (P_lo_all ** 3)
    n_lo = iv.mpf(max(N_G3_D.a, G['n_min'].a))
    assert (2 / iv.log(iv.mpf('1.05') / (1 + G['T_in']))).b < n_lo.a
    r_T5 = ((1 + G['T_in']) / iv.mpf('1.05')) ** (n_lo / 2)
    far_r5 = n_lo * r_T5
    assert far_r5.b < mp.mpf('1e-6')
    return W1_fin_out, W2_fin_out, W1_fin_all, W2_fin_all, r_T5, far_r5

def dzeta_err_bound(G, u_lo, u_hi, u0_lo, u0_hi, tau_hi, Del0_hi, r0ck2_hi, inside):
    if inside:
        dist_hi = u0_hi - u_lo
        dre = iv_max_hi(G['te_hi'] - 1, 1 - G['c_rein']) * tau_hi * (u0_hi ** 2 - u_lo ** 2)
        dim_1 = 2 * PI * (Del0_hi + G['C_D2'] * r0ck2_hi * (u0_hi ** 2)) * dist_hi
        dim_2 = 2 * PI * Del0_hi * iv.mpf('0.39')
        dim = iv_min_hi(dim_1, dim_2) if u0_hi.b <= mp.mpf('0.515') else dim_1
    else:
        dist_hi = u_hi - u0_lo
        # Exact local lower bound on -Re zeta'(t)/(2*tau*t) for t in [u0_lo, u_hi]:
        #   >= te_lo * (1 - (2/3)*(u_hi*thc_hi)^2) / (1 + gam * u_hi^2) (NO artificial T_c5 = 0.05 floor!!):
        c_re_loc = G['te_lo'] * (1 - (iv.mpf(2) / 3) * ((u_hi * G['thc_hi']) ** 2)) / (1 + G['gam'] * (u_hi ** 2))
        dre_coeff = iv_max_hi(G['te_hi'] - 1, 1 - c_re_loc)
        dre = dre_coeff * tau_hi * (u_hi ** 2 - u0_lo ** 2)
        dim = 2 * PI * (Del0_hi + G['C_D2'] * r0ck2_hi * (u_hi ** 2)) * dist_hi
    return iv.sqrt(dre ** 2 + dim ** 2)

def delta_fin_iv(G, u0_iv, Del0_hi, r0ck2_hi):
    del_par = iv.mpf('0.5') - u0_iv
    ad_hi = iv_abs_hi(del_par)
    u_top = iv.mpf(max(mp.mpf('0.5'), u0_iv.b))
    e1 = (Del0_hi + G['C_D2'] * r0ck2_hi * (u_top ** 2)) * ad_hi
    e2 = iv.mpf('0.39') * Del0_hi if u0_iv.b <= mp.mpf('0.505') else e1
    err = iv_min_hi(e1, e2)
    return del_par + iv_sym(err)

def H_fin0_iv(G, tau_iv, u0_iv, c_min_iv, Del0_hi):
    assert u0_iv.b <= mp.mpf('0.52'), f"H_fin0_iv u0 > 0.52: {u0_iv}"
    f_p_hi = (1 + 1 / c_min_iv) * (1 + 1 / (2 * c_min_iv))
    f_p_lo = (1 - 1 / c_min_iv) * (1 - 1 / (2 * c_min_iv))
    M_hi = f_p_hi * G['te_hi'] / (G['S'] ** 4)
    mob_in = f_p_lo * (G['S'] ** 4) * G['te_lo'] * G['sinc_in_sq']
    x_iv = tau_iv * u0_iv
    x_lo = iv.mpf(x_iv.a); x_hi = iv.mpf(x_iv.b)
    A_hi = (PI ** 2) * ((1 + Del0_hi) ** 2); B_hi = M_hi
    A_lo = (PI ** 2) * ((1 - Del0_hi) ** 2) / (1 + G['T_in']); B_lo = mob_in / (1 + G['T_in'])
    # Exact unimodal range of x / (A + B*x^2) on [x_lo, x_hi]: min is at endpoint, max at endpoint or sqrt(A/B)!
    h_lo = G['r_re0'] * iv_min(x_lo / (A_hi + B_hi * (x_lo ** 2)), x_hi / (A_hi + B_hi * (x_hi ** 2)))
    h_hi = G['te_hi'] * iv_max_hi(x_lo / (A_lo + B_lo * (x_lo ** 2)), x_hi / (A_lo + B_lo * (x_hi ** 2)))
    x_st = iv.sqrt(A_lo / B_lo)
    if x_lo.a <= x_st.b and x_st.a <= x_hi.b:
        h_hi = iv_max_hi(h_hi, G['te_hi'] / (2 * iv.sqrt(A_lo * B_lo)))
    return iv.mpf([h_lo.a, h_hi.b])

def make_params(ta, tb, ct_lo, ct_hi=None):
    # Exact integer lower bound r_0 >= ceil(tau_a / (0.45 * ct_hi)) since tau = ct * r_0 * (r_0*kap) <= 0.45 * ct_hi * r_0!
    r0_min = iv.mpf(1)
    if ct_hi is not None:
        r_lb = ta / (R0K_MAX * ct_hi)
        r0_int = 1
        while r_lb.a > iv.mpf(r0_int).b:
            r0_int += 1
        if r0_int > 1:
            assert r_lb.a > iv.mpf(r0_int - 1).b
        r0_min = iv.mpf(r0_int)
    c_min_iv = iv.mpf(max((r0_min * ct_lo).a, (ta / R0K_MAX).a))
    kap_hi = iv_min_hi(R0K_MAX / r0_min, tb / c_min_iv)
    # r_0*c*kap^2 <= 42 (for r_0=1: c*kap^2 < 42; for r_0>=2: (r_0^3*c*kap^2)/r_0^2 <= 60/4 = 15 < 42!)
    # and r_0*c*kap^2 = ct*(r_0*kap)^2 = tau^2/(ct*r_0^2) <= min(42, 0.45*tb, tb^2/(ct_lo*r0_min^2), 0.2025*ct_hi)!
    r0ck2_hi = iv_min_hi(iv_min_hi(CK2_MAX, R0K_MAX * tb), (tb ** 2) / (ct_lo * (r0_min ** 2)))
    if ct_hi is not None:
        r0ck2_hi = iv_min_hi(r0ck2_hi, (R0K_MAX ** 2) * ct_hi)
    G = box_geom(kap_hi)
    Del0_hi = G['C_D0'] * r0ck2_hi
    return G, c_min_iv, r0ck2_hi, Del0_hi

# ==============================================================================
# EVALUATOR A: g_{1,fin}(u_0) — INTERSECTION OF:
#   (1) Closed-form 1-step log1p_polar * W_1 parabola-diff (tight for small tau)
#   (2) ZERO-dz_e direct 2-step M_E * W_2 finite-n bound (tight for tau >= 5!)
# ==============================================================================
def _eval_g1_fin_single_iv(ta, tb, u0_iv, ct_lo, n_bins=56, ct_hi=None):
    assert u0_iv.b <= mp.mpf('0.52'), f"g1 u0 > 0.52: {u0_iv}"
    G, c_min_iv, r0ck2_hi, Del0_hi = make_params(ta, tb, ct_lo, ct_hi)
    tau_iv = iv.mpf([ta.a, tb.b])
    tau_lo = iv.mpf(ta.a); tau_hi = iv.mpf(tb.b)
    u0_lo = iv.mpf(u0_iv.a); u0_hi = iv.mpf(u0_iv.b)
    del_f = delta_fin_iv(G, u0_iv, Del0_hi, r0ck2_hi)

    # 1. Exact T_1 with cancellation of -2*del_par + 2*del_par AND 1-line linear drift bound (#77 fix 2c!):
    d_im_u0 = iv_max_hi(Del0_hi, G['C_D2'] * r0ck2_hi * (u0_hi ** 2))
    a_z = iv.mpf([(2 * G['r_re0'] * tau_lo * u0_lo).a, (2 * G['te_hi'] * tau_hi * u0_hi).b])
    b_z = 2 * PI * iv.mpf([(1 - d_im_u0).a, (1 + d_im_u0).b])
    del_par = iv.mpf('0.5') - u0_iv
    u_max_t1 = iv.mpf(max(mp.mpf('0.5'), u0_hi.b))
    conc_err = 2 * PI * (G['C_D2'] * r0ck2_hi * (u_max_t1 ** 2)) * iv_abs_hi(del_par)
    # Zero-dependency monotone form: a_z^2 / (a_z^2 + b_z^2) = 1 / (1 + (b_z/a_z)^2)!
    frac_a2 = 1 / (1 + ((b_z / a_z) ** 2))
    b_over_d = iv.mpf(b_z.b) / (iv.mpf(a_z.a) ** 2 + iv.mpf(b_z.a) ** 2)
    t1_iv = -2 * del_par * frac_a2 + iv_sym(2 * b_over_d * conc_err)

    # 2. Direct 2-step B_2 term at +-u_0 (for ZERO-dz_e direct W_2 enclosure):
    tc3r, _ = tau_c3_mv(tau_lo, tau_hi, u0_lo, u0_hi)
    fac_b2 = 2 * (PI ** 2) * (del_f ** 2) - (PI ** 2) / 6
    t2_par = fac_b2 * iv.mpf('0.5') * tc3r
    ep0 = eps_P_bound(G, tau_lo, tau_hi, u0_lo, u0_hi, c_min_iv, Del0_hi)
    inv2_0 = 1 / ((1 - ep0) ** 2)
    rel_W1_0 = ((G['te_hi'] - 1) + 2 * u0_hi * G['thc_hi']) * inv2_0 + (inv2_0 - 1)
    dz_p_lo = iv.sqrt((2 * tau_lo * u0_lo) ** 2 + (2 * PI) ** 2)
    rat_re0 = (2 * tau_hi * u0_hi) / iv.sqrt((2 * tau_hi * u0_hi) ** 2 + (2 * PI) ** 2)
    ep_dz0 = iv.sqrt((((G['te_hi'] - G['r_re0']) * rat_re0) ** 2) + ((2 * PI * d_im_u0 / dz_p_lo) ** 2))
    rel_t2 = (1 + rel_W1_0) / (1 - ep_dz0) - 1
    err_t2 = iv_abs_hi(fac_b2) * (tau_hi / 2) / (((tau_lo * u0_lo) ** 2 + PI ** 2) ** iv.mpf('1.5')) * rel_t2

    L_val = float(iv.sqrt(u0_hi ** 2 + iv.mpf('11.0') / (G['c_re5'] * tau_lo)).b)
    L_iv = iv.mpf(str(round(L_val, 2)))
    W1_fin_out, W2_fin_out, W1_fin_all, W2_fin_all, r_T5, far_r5 = W1_W2_fin_tail_ub(G, tau_lo, tau_hi, L_iv, c_min_iv, Del0_hi)
    u_lo_f = float(u0_lo.a); u_hi_f = float(u0_hi.b); L_f = float(L_iv.b)
    u_act = min(L_f, fsqrt(u_hi_f * u_hi_f + 4.5 / max(float(tau_lo.a), 0.5)))
    # Exact pure-mpf breakpoints at u0_lo and u0_hi (ZERO float64 slop around u0!):
    pts_iv = [iv.mpf(0)]
    for k in range(1, n_bins):
        p = u0_lo * k / n_bins
        if p.b < u0_lo.a: pts_iv.append(iv.mpf(p.a))
    pts_iv.append(u0_lo)
    if u0_hi.b > u0_lo.a: pts_iv.append(u0_hi)
    u_act_iv = iv.mpf(str(round(u_act, 6)))
    if u_act_iv.a <= u0_hi.b: u_act_iv = L_iv
    for k in range(1, n_bins + 1):
        p = u0_hi + (u_act_iv - u0_hi) * (iv.mpf(k) / n_bins) ** iv.mpf('1.25')
        if p.a > u0_hi.b and p.b < L_iv.a: pts_iv.append(iv.mpf(p.a))
    for k in range(1, 8):
        p = u_act_iv + (L_iv - u_act_iv) * k / 8
        if p.a > pts_iv[-1].b and p.b < L_iv.a: pts_iv.append(iv.mpf(p.a))
    if L_iv.a > pts_iv[-1].b: pts_iv.append(L_iv)

    i1_par = iv.mpf(0); err_i1 = iv.mpf(0); dir_w2 = iv.mpf(0)
    for i in range(len(pts_iv) - 1):
        ua = pts_iv[i]; ub = pts_iv[i + 1]; du = ub - ua
        if du.b <= mp.mpf('1e-14'): continue
        u_iv = iv.mpf([ua.a, ub.b])
        w1r, w1i = w1_par_mv(tau_lo, tau_hi, ua, ub)
        W1_abs = inv_quad_ub(tau_lo, tau_hi, ua, 1)
        W2_abs = 3 * inv_quad_ub(tau_lo, tau_hi, ua, 2)
        ep = eps_P_bound(G, tau_lo, tau_hi, ua, ub, c_min_iv, Del0_hi)
        inv2 = 1 / ((1 - ep) ** 2); inv4 = inv2 ** 2
        rel_W1 = ((G['te_hi'] - 1) + 2 * ub * G['thc_hi']) * inv2 + (inv2 - 1)
        P_over_2tau = iv.sqrt(ub ** 2 + (PI / tau_lo) ** 2) * (1 + ep)
        # Exact 1/(1-s_j) factor via thc_1ms (#106 fix #1!):
        rel_W2 = ((G['te_hi'] ** 2 - 1) + 4 * ub * G['thc_hi']) * inv4 + (inv4 - 1) + inv4 * (iv.mpf(4) / 3) * G['te_hi'] * G['thc_1ms'] * P_over_2tau

        if ub.b <= u0_lo.a:
            x_lo = tau_lo * (u0_lo ** 2 - ub ** 2)
            r_ub = iv.exp(-x_lo); r_iv = iv.exp(-tau_iv * (u0_iv ** 2 - u_iv ** 2))
            lr, li = log1p_polar(r_iv, PI - 2 * PI * u_iv, r_ub)
            pr, _ = cmul(-lr, -li, w1r, w1i)
            i1_par += pr * du
            x_f_lo = G['c_rein'] * tau_lo * (u0_lo ** 2 - ub ** 2)
            r_f_ub = iv.exp(-x_f_lo)
            dir_w2 += li2_pos_ub(r_f_ub) * (1 + rel_W2) * W2_abs * du
            dz_e = dzeta_err_bound(G, ua, ub, u0_lo, u0_hi, tau_hi, Del0_hi, r0ck2_hi, True)
            ph_lim = iv_max_hi(abs(PI - 2 * PI * ua), abs(PI - 2 * PI * ub)) + dz_e
            d_ang = iv.mpf(1) if ph_lim.b <= (PI / 2).a else (iv.sin(PI - ph_lim) if ph_lim.b < PI.a - mp.mpf('0.03') else iv.mpf(0))
            d_m = iv.mpf(max((1 - r_f_ub).a, d_ang.a))
            # Exact d_m & closed-form int_mlog_ub (#106 fix #3, ZERO hardcoded 0.05 / 14.0!):
            assert d_m.a > mp.mpf('1e-4'), f"g1 inside d_m={d_m} at u=[{ua},{ub}]"
            alpha_in = G['c_rein'] * tau_lo * (u0_lo + ua)
            imlog_du = int_mlog_ub(alpha_in, u0_lo - ub, u0_lo - ua, x_f_lo, r_f_ub, du)
            dl_ub = (r_f_ub / d_m) * dz_e
            l_ub_du = iv_min_hi((r_f_ub / d_m) * du, imlog_du)
            err_i1 += (dl_ub * du + l_ub_du * rel_W1) * W1_abs
        elif ua.a >= u0_hi.b:
            x_lo = tau_lo * (ua ** 2 - u0_hi ** 2)
            r_ub = iv.exp(-x_lo); r_iv = iv.exp(-tau_iv * (u_iv ** 2 - u0_iv ** 2))
            lr, li = log1p_polar(r_iv, PI + 2 * PI * u_iv, r_ub)
            pr, _ = cmul(-lr, -li, w1r, w1i)
            i1_par += pr * du
            x_f_lo = G['c_re5'] * tau_lo * (ua ** 2 - u0_hi ** 2)
            r_f_ub = iv.exp(-x_f_lo)
            dir_w2 += li2_pos_ub(r_f_ub) * (1 + rel_W2) * W2_abs * du
            dz_e = dzeta_err_bound(G, ua, ub, u0_lo, u0_hi, tau_hi, Del0_hi, r0ck2_hi, False)
            d_ang = iv.mpf(0)
            if ub.b <= mp.mpf('0.90'):
                ph_lim = iv_max_hi(abs(2 * PI * ua - PI), abs(2 * PI * ub - PI)) + dz_e
                d_ang = iv.mpf(1) if ph_lim.b <= (PI / 2).a else (iv.sin(PI - ph_lim) if ph_lim.b < PI.a - mp.mpf('0.03') else iv.mpf(0))
            d_m = iv.mpf(max((1 - r_f_ub).a, d_ang.a))
            assert d_m.a > mp.mpf('1e-4'), f"g1 outside d_m={d_m} at u=[{ua},{ub}]"
            alpha_out = G['c_re5'] * tau_lo * (ua + u0_hi)
            imlog_du = int_mlog_ub(alpha_out, ua - u0_hi, ub - u0_hi, x_f_lo, r_f_ub, du)
            dl_ub = (r_f_ub / d_m) * dz_e
            l_ub_du = iv_min_hi((r_f_ub / d_m) * du, imlog_du)
            err_i1 += (dl_ub * du + l_ub_du * rel_W1) * W1_abs
        else:
            ph_st0 = 2 * PI * iv_abs_hi(del_f)
            d_st0 = iv.mpf(1) if ph_st0.b <= (PI / 2).a else iv.sin(PI - ph_st0)
            assert d_st0.a > mp.mpf('0.05'), f"g1 straddle d_st0={d_st0}"
            err_i1 += (1 / d_st0) * W1_abs * (2 + rel_W1) * du
            dir_w2 += ((PI ** 2) / 6) * (1 + rel_W2) * W2_abs * du

    # Exact [L_iv, u_0.05] + [u_0.05, L_n] tails (#106 fix #5a, ZERO magic factors!):
    r_L = iv.exp(-G['c_re5'] * tau_lo * (L_iv ** 2 - u0_hi ** 2))
    I_tail_5 = r_L / (2 * G['c_re5'] * tau_lo * L_iv)
    r_L_par = iv.exp(-tau_lo * (L_iv ** 2 - u0_hi ** 2))
    I_tail_par = r_L_par / (2 * tau_lo * L_iv)
    W1_par_out = inv_quad_ub(tau_lo, tau_hi, L_iv, 1)
    tail_1 = (2 * (1 / (1 - r_L)) * I_tail_5 * W1_fin_out
              + 2 * (1 / (1 - r_L_par)) * I_tail_par * W1_par_out
              + (1 / (1 - r_T5)) * far_r5 * W1_fin_all)
    tail_2 = (2 * ((PI ** 2) / 6) * I_tail_5 * W2_fin_out
              + ((PI ** 2) / 6) * far_r5 * W2_fin_all)
    g1_par_diff = t1_iv + 2 * i1_par + iv_sym(2 * err_i1 + tail_1)
    g1_dir_W2   = t1_iv + t2_par + iv_sym(err_t2 + 2 * dir_w2 + tail_2)
    lo_int = max(g1_par_diff.a, g1_dir_W2.a)
    hi_int = min(g1_par_diff.b, g1_dir_W2.b)
    assert lo_int <= hi_int, f"g1 intersection empty: {g1_par_diff} vs {g1_dir_W2}"
    return iv.mpf([lo_int, hi_int])

def eval_g1_fin_iv(ta, tb, u0_iv, ct_lo, n_bins=56, ct_hi=None):
    """Exact hull of g_{1,fin}(u_0) over tau in [ta, tb], subdividing wide tau boxes to kill t1/i1 dependency."""
    rel_w = (float(tb.b) - float(ta.a)) / max(float(ta.a), 0.5)
    n_t = max(1, min(4, int(rel_w / 0.055 + 0.999999)))
    if n_t == 1:
        return _eval_g1_fin_single_iv(ta, tb, u0_iv, ct_lo, n_bins=n_bins, ct_hi=ct_hi)
    lo_all = None; hi_all = None
    for k in range(n_t):
        t_a = ta + (tb - ta) * k / n_t
        t_b = ta + (tb - ta) * (k + 1) / n_t
        g_k = _eval_g1_fin_single_iv(iv.mpf(t_a.a), iv.mpf(t_b.b), u0_iv, ct_lo, n_bins=n_bins, ct_hi=ct_hi)
        lo_all = g_k.a if lo_all is None else min(lo_all, g_k.a)
        hi_all = g_k.b if hi_all is None else max(hi_all, g_k.b)
    return iv.mpf([lo_all, hi_all])

# Exact continuous 2*pi-periodic cubic f_cub(a) = wrap(a)*(pi^2 - wrap(a)^2) on any interval a_iv!
def f_cub_periodic_iv(a_iv):
    k = int(round(float(a_iv.a) / (2.0 * fpi)))
    s = a_iv - 2 * k * PI
    while s.a < -PI.a:
        s = s + 2 * PI
    while s.a > PI.b:
        s = s - 2 * PI
    def eval_on(x_iv):
        c_pt = PI / iv.sqrt(3)
        vals = [iv.mpf(x_iv.a) * (PI ** 2 - iv.mpf(x_iv.a) ** 2), iv.mpf(x_iv.b) * (PI ** 2 - iv.mpf(x_iv.b) ** 2)]
        if x_iv.a <= c_pt.b and c_pt.a <= x_iv.b:
            vals.append(c_pt * (PI ** 2 - c_pt ** 2))
        if x_iv.a <= -c_pt.a and -c_pt.b <= x_iv.b:
            vals.append((-c_pt) * (PI ** 2 - c_pt ** 2))
        return iv.mpf([min(v.a for v in vals), max(v.b for v in vals)])
    if s.b <= PI.b:
        return eval_on(iv.mpf([max(s.a, -PI.b), min(s.b, PI.b)]))
    v1 = eval_on(iv.mpf([max(s.a, -PI.b), PI.b]))
    v2 = eval_on(iv.mpf([-PI.b, min((s - 2 * PI).b, PI.b)]))
    return iv.mpf([min(v1.a, v2.a), max(v1.b, v2.b)])

# ==============================================================================
# EVALUATOR B: EXACT PARABOLA Li_3 / W_2 VALUE + 2ND-POWER W_1 FINITE-n PERTURBATION!
#   - Parabola integral -int (ME(th) - ME(0)) W_1 du is evaluated EXACTLY via its
#     1-line Li_3 / W_2 identity (verified to 1e-9 above!): O(1/M^2) tail => width < 0.004!
#   - Finite-n difference |(Psi_fin - Psi_par)| uses ONLY 2nd-power W_1 (err <= 0.034)!
# ==============================================================================
def _prepare_Psi_W1_cache(ta, tb, u0_iv, ct_lo, ct_hi=None, allow_cub=False):
    """Lossless precomputation of all theta-independent u-grid quantities for eval_Psi_diff_W1_iv (exact same 50-dps iv math!)."""
    assert u0_iv.b <= mp.mpf('0.52'), f"W1 u0 > 0.52: {u0_iv}"
    G, c_min_iv, r0ck2_hi, Del0_hi = make_params(ta, tb, ct_lo, ct_hi)
    tau_iv = iv.mpf([ta.a, tb.b]); tau_lo = iv.mpf(ta.a); tau_hi = iv.mpf(tb.b)
    u0_lo = iv.mpf(u0_iv.a); u0_hi = iv.mpf(u0_iv.b)
    del_par = iv.mpf('0.5') - u0_iv
    del_f = delta_fin_iv(G, u0_iv, Del0_hi, r0ck2_hi)
    ad_hi = iv_abs_hi(del_f)
    ad_lo = iv.mpf(abs(del_f).a)
    th_hit_lo = PI * (1 - 2 * ad_hi)
    th_hit_hi = PI * (1 - 2 * ad_lo)
    H_f = H_fin0_iv(G, tau_iv, u0_iv, c_min_iv, Del0_hi)
    tc3r, tc3i = tau_c3_mv(tau_lo, tau_hi, u0_lo, u0_hi)
    fp_0 = f_cub_periodic_iv(-2 * PI * del_par)
    fm_0 = f_cub_periodic_iv( 2 * PI * del_par)

    L_val = float(iv.sqrt(u0_hi ** 2 + iv.mpf('11.0') / (G['c_re5'] * tau_lo)).b)
    L_iv = iv.mpf(str(round(L_val, 2)))
    W1_fin_out, _, W1_fin_all, _, r_T5, far_r5 = W1_W2_fin_tail_ub(G, tau_lo, tau_hi, L_iv, c_min_iv, Del0_hi)
    u_lo_f = float(u0_lo.a); u_hi_f = float(u0_hi.b); L_f = float(L_iv.b)
    u_act = min(L_f, fsqrt(u_hi_f * u_hi_f + 4.0 / max(float(tau_lo.a), 0.5)))
    in_raw = [0.0]
    n_in = max(10, int(round(u_lo_f / 0.032)))
    for k in range(1, n_in): in_raw.append(u_lo_f * k / n_in)
    if u_lo_f > 0.02:
        in_raw.append(u_lo_f - 0.008); in_raw.append(u_lo_f - 0.002)
    in_pts = [iv.mpf(str(round(x, 6))) for x in sorted(set(in_raw)) if 0.0 <= x < u_lo_f - 1e-7]
    in_pts = [p for p in in_pts if p.b < u0_lo.a]
    out_raw = [min(u_act, u_hi_f + 0.002), min(u_act, u_hi_f + 0.008)]
    n_out1 = max(10, int(round((u_act - u_hi_f) / 0.036)))
    for k in range(1, n_out1 + 1): out_raw.append(u_hi_f + (u_act - u_hi_f) * k / n_out1)
    for k in range(1, 6): out_raw.append(u_act + (L_f - u_act) * k / 6.0)
    out_pts = [iv.mpf(str(round(x, 6))) for x in sorted(set(out_raw)) if u_hi_f + 1e-7 < x < L_f - 1e-7]
    out_pts = [p for p in out_pts if p.a > u0_hi.b and p.b < L_iv.a]
    pos_pts_iv = in_pts + [u0_lo] + ([u0_hi] if u0_hi.b > u0_lo.a else []) + out_pts + [L_iv]

    bins = []
    for i in range(len(pos_pts_iv) - 1):
        ua = pos_pts_iv[i]; ub = pos_pts_iv[i + 1]; du = ub - ua
        if du.b <= mp.mpf('1e-14'): continue
        u_iv = iv.mpf([ua.a, ub.b])
        w2r, w2i = w2_par_mv(tau_lo, tau_hi, ua, ub)
        W1_abs = inv_quad_ub(tau_lo, tau_hi, ua, 1)
        ep = eps_P_bound(G, tau_lo, tau_hi, ua, ub, c_min_iv, Del0_hi)
        inv2 = 1 / ((1 - ep) ** 2)
        rel_W1 = ((G['te_hi'] - 1) + 2 * ub * G['thc_hi']) * inv2 + (inv2 - 1)
        if ub.b <= u0_lo.a or ua.a >= u0_hi.b:
            inside = (ub.b <= u0_lo.a)
            if inside:
                x_lo = tau_lo * (u0_lo ** 2 - ub ** 2)
                r_ub = iv.exp(-x_lo); r_iv = iv.exp(-tau_iv * (u0_iv ** 2 - u_iv ** 2))
                ph0_iv = PI - 2 * PI * u_iv
                x_f_lo = G['c_rein'] * tau_lo * (u0_lo ** 2 - ub ** 2)
                r_f_ub = iv.exp(-x_f_lo)
                alpha_loc = G['c_rein'] * tau_lo * (u0_lo + ua)
                imlog_du = int_mlog_ub(alpha_loc, u0_lo - ub, u0_lo - ua, x_f_lo, r_f_ub, du)
                dz_e = dzeta_err_bound(G, ua, ub, u0_lo, u0_hi, tau_hi, Del0_hi, r0ck2_hi, True)
                ph_abs_u = iv_max_hi(abs(PI - 2 * PI * ua), abs(PI - 2 * PI * ub))
                sgn_th = -1
            else:
                x_lo = tau_lo * (ua ** 2 - u0_hi ** 2)
                r_ub = iv.exp(-x_lo); r_iv = iv.exp(-tau_iv * (u_iv ** 2 - u0_iv ** 2))
                ph0_iv = PI + 2 * PI * u_iv
                x_f_lo = G['c_re5'] * tau_lo * (ua ** 2 - u0_hi ** 2)
                r_f_ub = iv.exp(-x_f_lo)
                alpha_loc = G['c_re5'] * tau_lo * (ua + u0_hi)
                imlog_du = int_mlog_ub(alpha_loc, ua - u0_hi, ub - u0_hi, x_f_lo, r_f_ub, du)
                dz_e = dzeta_err_bound(G, ua, ub, u0_lo, u0_hi, tau_hi, Del0_hi, r0ck2_hi, False)
                ph_abs_u = iv_max_hi(abs(2 * PI * ua - PI), abs(2 * PI * ub - PI)) if ub.b <= mp.mpf('0.90') else None
                sgn_th = +1
            ph_w = float(ph0_iv.b) - float(ph0_iv.a)
            rv = float(r_ub.b)
            M_cap_cub = 8 if rv < 0.45 else (14 if rv < 0.75 else 24)
            M_cap_ev  = 6 if rv < 0.45 else (10 if rv < 0.75 else 18)
            M_cap = M_cap_cub if allow_cub else M_cap_ev
            M_lim = min(M_cap, max(5, int(1.8 / max(ph_w, 1e-6))))
            harm = []
            rp = r_iv
            for m in range(1, M_lim + 1):
                base = ((-1 if (m % 2 == 1) else 1) * rp) / (m * m * m)
                harm.append((base, iv.cos(m * ph0_iv), iv.sin(m * ph0_iv)))
                rp = rp * r_iv
            m_re_triv_du = 2 * li2_pos_ub(r_f_ub) * du
            l_re_triv_du = 2 * imlog_du
            d_m0 = None
            if allow_cub:
                if ph_abs_u is not None:
                    ph0_lim = ph_abs_u + dz_e
                    d0_ang = iv.mpf(1) if ph0_lim.b <= (PI / 2).a else (iv.sin(PI - ph0_lim) if ph0_lim.b < PI.a - mp.mpf('0.03') else iv.mpf(0))
                else:
                    d0_ang = iv.mpf(0)
                d_m0 = iv.mpf(max((1 - r_f_ub).a, d0_ang.a))
                assert d_m0.a > mp.mpf('1e-4'), f"W1 d_m0={d_m0} at u=[{ua},{ub}]"
            bins.append(('side', du, w2r, w2i, W1_abs, rel_W1, r_ub, rv, ph_w, M_cap_ev, M_cap_cub,
                         harm, sgn_th, r_f_ub, imlog_du, dz_e, m_re_triv_du, l_re_triv_du, ph_abs_u, d_m0))
        else:
            W2_abs = 3 * inv_quad_ub(tau_lo, tau_hi, ua, 2)
            bins.append(('straddle', du, W1_abs, W2_abs, rel_W1))

    r_L = iv.exp(-G['c_re5'] * tau_lo * (L_iv ** 2 - u0_hi ** 2))
    I_tail_5 = r_L / (2 * G['c_re5'] * tau_lo * L_iv)
    r_L_par = iv.exp(-tau_lo * (L_iv ** 2 - u0_hi ** 2))
    I_tail_par = r_L_par / (2 * tau_lo * L_iv)
    W1_par_out = inv_quad_ub(tau_lo, tau_hi, L_iv, 1)
    W2_par_out = 3 * inv_quad_ub(tau_lo, tau_hi, L_iv, 2)
    tail_1_fin_re = 4 * ((PI ** 2) / 6) * I_tail_5 * W1_fin_out + 2 * ((PI ** 2) / 6) * far_r5 * W1_fin_all
    tail_1_par_re = 4 * ((PI ** 2) / 6) * I_tail_par * W1_par_out
    tail_ne_re    = 4 * zeta3_D * I_tail_par * W2_par_out
    return dict(allow_cub=allow_cub, del_par=del_par, del_f=del_f, th_hit_lo=th_hit_lo, th_hit_hi=th_hit_hi,
                H_f=H_f, tc3r=tc3r, tc3i=tc3i, fp_0=fp_0, fm_0=fm_0, bins=bins,
                r_L=r_L, I_tail_5=I_tail_5, r_L_par=r_L_par, I_tail_par=I_tail_par,
                W1_fin_out=W1_fin_out, W1_fin_all=W1_fin_all, r_T5=r_T5, far_r5=far_r5,
                W1_par_out=W1_par_out, W2_par_out=W2_par_out,
                tail_1_fin_re=tail_1_fin_re, tail_1_par_re=tail_1_par_re, tail_ne_re=tail_ne_re)

def eval_Psi_diff_W1_iv(ta, tb, u0_iv, th_iv, ct_lo, need_im=False, g1_iv=None, ct_hi=None, cache=None):
    if cache is None or (need_im and not cache['allow_cub']):
        cache = _prepare_Psi_W1_cache(ta, tb, u0_iv, ct_lo, ct_hi=ct_hi, allow_cub=need_im)
    del_par = cache['del_par']; del_f = cache['del_f']
    th_hit_lo = cache['th_hit_lo']; th_hit_hi = cache['th_hit_hi']
    H_f = cache['H_f']; tc3r = cache['tc3r']; tc3i = cache['tc3i']

    th_a = iv.mpf(th_iv.a); th_b = iv.mpf(th_iv.b)
    th_w = float(th_iv.b) - float(th_iv.a)
    on_peak = (th_b.b < th_hit_lo.a - mp.mpf('0.02')) and ((2 * PI * abs(del_par) + th_b).b < PI.a - mp.mpf('0.01'))

    # Direct 1-step W_1 boundary drop (ZERO parabola comparison, ZERO bdy_n_re subtraction!):
    if on_peak:
        main_1s_lo = iv.mpf('0.5') * iv.mpf(H_f.a) * (th_a ** 2)
        main_1s_hi = iv.mpf('0.5') * iv.mpf(H_f.b) * (th_b ** 2)
        coeff_re2 = iv.mpf('0.5') * H_f - iv.mpf('0.5') * PI * del_par * tc3i
        main_re_lo = iv.mpf(coeff_re2.a) * (th_a ** 2) if coeff_re2.a >= 0 else iv.mpf(coeff_re2.a) * (th_b ** 2)
        main_re_hi = iv.mpf(coeff_re2.b) * (th_b ** 2) if coeff_re2.b >= 0 else iv.mpf(coeff_re2.b) * (th_a ** 2)
        bdy_cub_im = (tc3r / 12) * (th_iv ** 3)
    else:
        f_a_lo = iv.mpf('0.5') * (th_a ** 2) - PI * iv.mpf(max(mp.mpf(0), (th_a - th_hit_lo).b))
        f_b_lo = iv.mpf('0.5') * (th_b ** 2) - PI * iv.mpf(max(mp.mpf(0), (th_b - th_hit_lo).b))
        f_min = iv_min(f_a_lo, f_b_lo)
        f_a_hi = iv.mpf('0.5') * (th_a ** 2) - PI * iv.mpf(max(mp.mpf(0), (th_a - th_hit_hi).a))
        f_b_hi = iv.mpf('0.5') * (th_b ** 2) - PI * iv.mpf(max(mp.mpf(0), (th_b - th_hit_hi).a))
        f_max = iv_max_hi(f_a_hi, f_b_hi)
        if th_a.a <= th_hit_hi.b and th_hit_hi.a <= th_b.b:
            f_max = iv_max_hi(f_max, iv.mpf('0.5') * (th_hit_hi ** 2))
        main_1s_lo = iv.mpf(H_f.a) * f_min if f_min.a >= 0 else iv.mpf(H_f.b) * f_min
        main_1s_hi = iv.mpf(H_f.b) * f_max if f_max.b >= 0 else iv.mpf(H_f.a) * f_max
        fp_th = f_cub_periodic_iv(-2 * PI * del_par + th_iv)
        fm_th = f_cub_periodic_iv( 2 * PI * del_par + th_iv)
        dfp = fp_th - cache['fp_0']; dfm = fm_th - cache['fm_0']
        bdy_n_re = ((dfp - dfm) * tc3i) / 24
        main_re_lo = iv.mpf((main_1s_lo - iv.mpf(bdy_n_re.b)).a)
        main_re_hi = iv.mpf((main_1s_hi - iv.mpf(bdy_n_re.a)).b)
        bdy_cub_im = iv.mpf(0)

    # Precompute theta-only harmonic factors ONCE for this theta-bin (24 trig calls instead of 24*65!):
    M_top = 24 if need_im else 18
    fe_arr = [None] + [-2 * (iv.sin(m * th_iv / 2) ** 2) for m in range(1, M_top + 1)]
    fc_arr = ([None] + [iv.sin(m * th_iv) - m * th_iv for m in range(1, M_top + 1)]) if need_im else None

    par_ne_re = iv.mpf(0); par_ne_cub_im = iv.mpf(0); err_re_tot = iv.mpf(0); err_im_cub_tot = iv.mpf(0)
    dir_re_tot = iv.mpf(0); dir_im_cub_tot = iv.mpf(0)
    for b in cache['bins']:
        if b[0] == 'side':
            (_, du, w2r, w2i, W1_abs, rel_W1, r_ub, rv, ph_w, M_cap_ev, M_cap_cub,
             harm, sgn_th, r_f_ub, imlog_du, dz_e, m_re_du, l_re_du, ph_abs_u, d_m0) = b
            M_max = M_cap_cub if need_im else M_cap_ev
            M_use = min(M_max, max(5, int(1.8 / max(ph_w + 0.5 * th_w, 1e-6))))
            er = iv.mpf(0); ei = iv.mpf(0); cr = iv.mpf(0); ci = iv.mpf(0)
            for m in range(1, M_use + 1):
                base, cm, sm = harm[m - 1]
                bf = base * fe_arr[m]
                er += bf * cm; ei += bf * sm
                if need_im:
                    bfc = base * (sgn_th * fc_arr[m])
                    cr += bfc * (-sm); ci += bfc * cm
            r_pow = r_ub ** (M_use + 1)
            t_int = r_pow / (2 * (M_use ** 2))
            if rv < 0.99:
                t_geo = r_pow / (((M_use + 1) ** 3) * (1 - r_ub))
                t_base = iv_min_hi(t_int, t_geo)
            else:
                t_base = t_int
            t_even = 2 * t_base
            if rv < 0.99:
                t_even = iv_min_hi(t_even, iv.mpf('0.5') * (th_b ** 2) * r_pow / ((M_use + 1) * (1 - r_ub)))
            er = er + iv_sym(t_even); ei = ei + iv_sym(t_even)
            pr, _ = cmul(-er, -ei, w2r, w2i)
            par_ne_re += 2 * pr * du
            if need_im:
                t_cub = t_base + th_b * (r_pow / M_use)
                if rv < 0.99:
                    t_cub = iv_min_hi(t_cub, ((th_b ** 3) / 6) * r_pow / (1 - r_ub))
                cr = cr + iv_sym(t_cub); ci = ci + iv_sym(t_cub)
                _, pi_ = cmul(-cr, -ci, w2r, w2i)
                par_ne_cub_im += 2 * pi_ * du
            if need_im or on_peak:
                if ph_abs_u is not None:
                    ph_lim = ph_abs_u + th_b + dz_e
                    d_lo_ang = iv.mpf(1) if ph_lim.b <= (PI / 2).a else (iv.sin(PI - ph_lim) if ph_lim.b < PI.a - mp.mpf('0.05') else iv.mpf(0))
                else:
                    d_lo_ang = iv.mpf(0)
                d_m = iv.mpf(max((1 - r_f_ub).a, d_lo_ang.a))
                if d_m.a > mp.mpf('0.08'):
                    m_re_du = iv_min_hi(m_re_du, (iv.mpf('0.5') * (th_b ** 2) * r_f_ub / d_m) * du)
                    l_re_du = iv_min_hi(l_re_du, (iv.mpf('0.5') * (th_b ** 2) * r_f_ub / (d_m ** 2)) * du)
                if need_im:
                    m_im_triv_du = ((PI ** 2) / 3) * r_f_ub * du + th_b * iv_min_hi((r_f_ub / d_m0) * du, imlog_du)
                    l_im_triv_du = 2 * imlog_du + th_b * (r_f_ub / (d_m0 ** 2)) * du
                    if d_m.a > mp.mpf('0.08'):
                        m_im_du = iv_min_hi(m_im_triv_du, (((th_b ** 3) / 6) * r_f_ub / (d_m ** 2)) * du)
                        l_im_du = iv_min_hi(l_im_triv_du, (((th_b ** 3) / 3) * r_f_ub / (d_m ** 3)) * du)
                    else:
                        m_im_du = m_im_triv_du; l_im_du = l_im_triv_du
                    err_im_cub_tot += 2 * (l_im_du * dz_e + m_im_du * rel_W1) * W1_abs
                    dir_im_cub_tot += 2 * m_im_du * (1 + rel_W1) * W1_abs
            err_re_tot += 2 * (l_re_du * dz_e + m_re_du * rel_W1) * W1_abs
            dir_re_tot += 2 * m_re_du * (1 + rel_W1) * W1_abs
        else:
            _, du, W1_abs, W2_abs, rel_W1 = b
            ph_st = 2 * PI * iv_abs_hi(del_f) + th_b
            d_st = iv.mpf(1) if ph_st.b <= (PI / 2).a else (iv.sin(PI - ph_st) if ph_st.b < PI.a - mp.mpf('0.05') else iv.mpf(0))
            if on_peak and d_st.a > mp.mpf('0.10'):
                m_re_st = iv_min_hi(2 * ((PI ** 2) / 6), iv.mpf('0.5') * (th_b ** 2) / d_st)
                n_re_st = iv_min_hi(2 * zeta3_D, iv.mpf('0.5') * (th_b ** 2) / d_st)
                dir_re_tot += 2 * m_re_st * W1_abs * (1 + rel_W1) * du
                err_re_tot += (2 * m_re_st * W1_abs * (2 + rel_W1) + 2 * n_re_st * W2_abs) * du
                if need_im:
                    ph_st0 = 2 * PI * iv_abs_hi(del_f)
                    d_st0 = iv.mpf(1) if ph_st0.b <= (PI / 2).a else iv.sin(PI - ph_st0)
                    assert d_st0.a > mp.mpf('0.05'), f"W1 straddle d_st0={d_st0}"
                    m_im_u = 2 * ((PI ** 2) / 6) + th_b / d_st0
                    m_im_st = iv_min_hi(m_im_u, ((th_b ** 3) / 6) / (d_st ** 2))
                    n_im_st = iv_min_hi(zeta3_D + th_b * (PI ** 2) / 6, ((th_b ** 3) / 6) / d_st)
                    dir_im_cub_tot += 2 * m_im_st * W1_abs * (1 + rel_W1) * du
                    err_im_cub_tot += (2 * m_im_st * W1_abs * (2 + rel_W1) + 2 * n_im_st * W2_abs) * du
            else:
                dir_re_tot += 4 * ((PI ** 2) / 6) * W1_abs * (1 + rel_W1) * du
                err_re_tot += (4 * ((PI ** 2) / 6) * W1_abs * (2 + rel_W1) + 4 * zeta3_D * W2_abs) * du
                if need_im:
                    ph_st0 = 2 * PI * iv_abs_hi(del_f)
                    d_st0 = iv.mpf(1) if ph_st0.b <= (PI / 2).a else iv.sin(PI - ph_st0)
                    assert d_st0.a > mp.mpf('0.05'), f"W1 straddle d_st0={d_st0}"
                    m_im_u = 2 * ((PI ** 2) / 6) + th_b / d_st0
                    n_im_u = zeta3_D + th_b * (PI ** 2) / 6
                    dir_im_cub_tot += 2 * m_im_u * W1_abs * (1 + rel_W1) * du
                    err_im_cub_tot += (2 * m_im_u * W1_abs * (2 + rel_W1) + 2 * n_im_u * W2_abs) * du

    err_re_tot += cache['tail_1_fin_re'] + cache['tail_1_par_re'] + cache['tail_ne_re']
    dir_re_tot += cache['tail_1_fin_re']
    if need_im:
        M_im_5   = 2 * ((PI ** 2) / 6) + th_b / (1 - cache['r_L'])
        M_im_far = 2 * ((PI ** 2) / 6) + th_b / (1 - cache['r_T5'])
        M_im_par = 2 * ((PI ** 2) / 6) + th_b / (1 - cache['r_L_par'])
        tail_1_fin_im = 2 * M_im_5 * cache['I_tail_5'] * cache['W1_fin_out'] + M_im_far * cache['far_r5'] * cache['W1_fin_all']
        tail_1_par_im = 2 * M_im_par * cache['I_tail_par'] * cache['W1_par_out']
        tail_ne_im    = 2 * (zeta3_D + th_b * (PI ** 2) / 6) * cache['I_tail_par'] * cache['W2_par_out']
        err_im_cub_tot += tail_1_fin_im + tail_1_par_im + tail_ne_im
        dir_im_cub_tot += tail_1_fin_im
    # INTERSECTION of parabola-diff W_1 bound AND ZERO-dz_e direct 1-step W_1 bound!!
    drop_lo_par = iv.mpf((main_re_lo - iv.mpf(par_ne_re.b) - iv.mpf(err_re_tot.b)).a)
    drop_hi_par = iv.mpf((main_re_hi - iv.mpf(par_ne_re.a) + iv.mpf(err_re_tot.b)).b)
    drop_lo_dir = iv.mpf((main_1s_lo - iv.mpf(dir_re_tot.b)).a)
    drop_hi_dir = iv.mpf((main_1s_hi + iv.mpf(dir_re_tot.b)).b)
    drop_lo = iv.mpf(max(drop_lo_par.a, drop_lo_dir.a))
    drop_hi = iv.mpf(min(drop_hi_par.b, drop_hi_dir.b))
    assert drop_lo.a <= drop_hi.b, f"drop intersection empty: par=[{drop_lo_par},{drop_hi_par}] vs dir=[{drop_lo_dir},{drop_hi_dir}]"
    if not need_im:
        return drop_lo
    if g1_iv is None:
        g1_iv = eval_g1_fin_iv(ta, tb, u0_iv, ct_lo, ct_hi=ct_hi)
    im_par = g1_iv * th_iv + bdy_cub_im + par_ne_cub_im + iv_sym(err_im_cub_tot)
    im_dir = g1_iv * th_iv + iv_sym(dir_im_cub_tot) if on_peak else im_par
    im_hi = iv_min_hi(iv_abs_hi(im_par), iv_abs_hi(im_dir))
    return drop_lo, drop_hi, im_hi, iv_min_hi(err_re_tot, dir_re_tot)

def eval_valley_drop_iv(ta, tb, u0_iv, th_iv, ct_lo, ct_hi=None, cache=None):
    return eval_Psi_diff_W1_iv(ta, tb, u0_iv, th_iv, ct_lo, need_im=False, ct_hi=ct_hi, cache=cache)

# ==============================================================================
# EVALUATOR C: PEAK CURVATURE -Re Psi_fin''(th) AND |Im Psi_fin''(th)| VIA W_2!
# ==============================================================================
def _eval_peak_curv_single_iv(ta, tb, u0_iv, th_iv, ct_lo, ct_hi=None):
    assert u0_iv.b <= mp.mpf('0.52'), f"curv u0 > 0.52: {u0_iv}"
    G, c_min_iv, r0ck2_hi, Del0_hi = make_params(ta, tb, ct_lo, ct_hi)
    tau_iv = iv.mpf([ta.a, tb.b]); tau_lo = iv.mpf(ta.a); tau_hi = iv.mpf(tb.b)
    u0_lo = iv.mpf(u0_iv.a); u0_hi = iv.mpf(u0_iv.b)
    del_f = delta_fin_iv(G, u0_iv, Del0_hi, r0ck2_hi)
    H_f = H_fin0_iv(G, tau_iv, u0_iv, c_min_iv, Del0_hi)

    ap = -2 * PI * del_f + th_iv; am = 2 * PI * del_f + th_iv
    assert ap.a > -PI.a + mp.mpf('0.12') and ap.b < PI.a - mp.mpf('0.12')
    assert am.a > -PI.a + mp.mpf('0.12') and am.b < PI.a - mp.mpf('0.12')
    # Exact zero-cancellation identity using tau_c3_mv:
    #   (tau/4)*(ap*c1i - am*c2i) = +PI * del_f * tc3i,   (tau/4)*(-ap*c1r + am*c2r) = +0.5 * th_iv * tc3r!
    tc3r, tc3i = tau_c3_mv(tau_lo, tau_hi, u0_lo, u0_hi)
    b1_re = PI * del_f * tc3i
    b1_im = iv.mpf('0.5') * th_iv * tc3r

    ep0 = eps_P_bound(G, tau_lo, tau_hi, u0_lo, u0_hi, c_min_iv, Del0_hi)
    inv2_0 = 1 / ((1 - ep0) ** 2)
    rel_W1_0 = ((G['te_hi'] - 1) + 2 * u0_hi * G['thc_hi']) * inv2_0 + (inv2_0 - 1)
    d_im_u0 = iv_max_hi(Del0_hi, G['C_D2'] * r0ck2_hi * (u0_hi ** 2))
    dz_p_lo = iv.sqrt((2 * tau_lo * u0_lo) ** 2 + (2 * PI) ** 2)
    rat_re0 = (2 * tau_hi * u0_hi) / iv.sqrt((2 * tau_hi * u0_hi) ** 2 + (2 * PI) ** 2)
    ep_dz0 = iv.sqrt((((G['te_hi'] - G['r_re0']) * rat_re0) ** 2) + ((2 * PI * d_im_u0 / dz_p_lo) ** 2))
    rel_B1 = (1 + rel_W1_0) / (1 - ep_dz0) - 1
    err_B1 = (tau_hi / 4) * (iv_abs_hi(ap) + iv_abs_hi(am)) / (((tau_lo * u0_lo) ** 2 + PI ** 2) ** iv.mpf('1.5')) * rel_B1
    # Exact O(theta) bound on |Im Psi_fin''(th)| <= M_3 * th (since Im Psi_fin''(0) == 0 identically by u -> -u symmetry!):
    b1_im_coeff = iv.mpf('0.5') * iv_abs_hi(tc3r) + (tau_hi / 2) / (((tau_lo * u0_lo) ** 2 + PI ** 2) ** iv.mpf('1.5')) * rel_B1

    L_val = float(iv.sqrt(u0_hi ** 2 + iv.mpf('11.0') / (G['c_re5'] * tau_lo)).b)
    L_iv = iv.mpf(str(round(L_val, 2)))
    _, W2_fin_out, _, W2_fin_all, r_T5, far_r5 = W1_W2_fin_tail_ub(G, tau_lo, tau_hi, L_iv, c_min_iv, Del0_hi)
    # Exact pure-mpf breakpoints at +-u0_lo and +-u0_hi:
    pos_iv = [iv.mpf(0)]
    for k in range(1, 16):
        p = u0_lo * k / 16
        if p.b < u0_lo.a: pos_iv.append(iv.mpf(p.a))
    pos_iv.append(u0_lo)
    if u0_hi.b > u0_lo.a: pos_iv.append(u0_hi)
    for k in range(1, 12):
        p = u0_hi + (L_iv - u0_hi) * (iv.mpf(k) / 12) ** iv.mpf('1.5')
        if p.a > u0_hi.b and p.b < L_iv.a: pos_iv.append(iv.mpf(p.a))
    pos_iv.append(L_iv)
    pts_iv = [-x for x in reversed(pos_iv[1:])] + pos_iv

    i2_re = iv.mpf(0); eL2 = iv.mpf(0); eW2 = iv.mpf(0); m3_int = iv.mpf(0)
    for i in range(len(pts_iv) - 1):
        ua = pts_iv[i]; ub = pts_iv[i + 1]; du = ub - ua
        if du.b <= mp.mpf('1e-14'): continue
        u_iv = iv.mpf([ua.a, ub.b])
        au_lo = iv.mpf(0) if (ua.a <= 0 <= ub.b) else iv.mpf(min(abs(ua.a), abs(ub.b)))
        au_hi = iv.mpf(max(abs(ua.a), abs(ub.b)))
        assert ua.a >= 0 or ub.b <= 0
        w2r, w2i = w2_par_mv(tau_lo, tau_hi, au_lo, au_hi)
        if ub.b <= 0:
            w2i = -w2i
        W2_abs = 3 * inv_quad_ub(tau_lo, tau_hi, au_lo, 2)
        ep = eps_P_bound(G, tau_lo, tau_hi, au_lo, au_hi, c_min_iv, Del0_hi)
        inv4 = 1 / ((1 - ep) ** 4)
        P_over_2tau = iv.sqrt(au_hi ** 2 + (PI / tau_lo) ** 2) * (1 + ep)
        # Exact 1/(1-s_j) factor via thc_1ms (#106 fix #1!):
        rel_W2 = ((G['te_hi'] ** 2 - 1) + 4 * au_hi * G['thc_hi']) * inv4 + (inv4 - 1) + inv4 * (iv.mpf(4) / 3) * G['te_hi'] * G['thc_1ms'] * P_over_2tau

        if au_hi.b <= u0_lo.a:
            x_lo = tau_lo * (u0_lo ** 2 - au_hi ** 2)
            r_ub = iv.exp(-x_lo); r_iv = iv.exp(-tau_iv * (u0_iv ** 2 - u_iv ** 2))
            ph_iv = PI - 2 * PI * u_iv - th_iv
            lr, li = log1p_polar(r_iv, ph_iv, r_ub)
            pr, _ = cmul(lr, li, w2r, w2i)
            i2_re += pr * du
            x_f_lo = G['c_rein'] * tau_lo * (u0_lo ** 2 - au_hi ** 2)
            r_f_ub = iv.exp(-x_f_lo)
            dz_e = dzeta_err_bound(G, au_lo, au_hi, u0_lo, u0_hi, tau_hi, Del0_hi, r0ck2_hi, True)
            # Distance to odd*pi over ALL t in [0, th_b] (valid for both Re at [th_a, th_b] AND M_3 FTC on [0, th_b]!):
            ph_a = (PI - 2 * PI * ub) - iv.mpf(th_iv.b) - iv.mpf(dz_e.b)
            ph_b = (PI - 2 * PI * ua) + iv.mpf(dz_e.b)
            k_min = int(round((0.5 * (float(ph_a.a) + float(ph_b.b)) / fpi - 1.0) / 2.0))
            dist_odd = PI
            for kk in (k_min - 1, k_min, k_min + 1):
                target = (2 * kk + 1) * PI
                if ph_a.a <= target.b and target.a <= ph_b.b:
                    dist_odd = iv.mpf(0)
                else:
                    dist_odd = iv_min(dist_odd, iv_min(abs(ph_a - target), abs(ph_b - target)))
            d_ang = iv.mpf(1) if dist_odd.a >= (PI / 2).b else (iv.sin(dist_odd) if dist_odd.a > mp.mpf('0.03') else iv.mpf(0))
            d_min = iv.mpf(max((1 - r_f_ub).a, d_ang.a))
            assert d_min.a > mp.mpf('0.005'), f"W2 inside d_min={d_min} at u={ua}..{ub}"
            eL2 += (r_f_ub / d_min) * dz_e * W2_abs * du
            eW2 += (r_f_ub / d_min) * rel_W2 * W2_abs * du
            m3_int += (r_f_ub / d_min) * (1 + rel_W2) * W2_abs * du
        elif au_lo.a >= u0_hi.b:
            x_lo = tau_lo * (au_lo ** 2 - u0_hi ** 2)
            r_ub = iv.exp(-x_lo); r_iv = iv.exp(-tau_iv * (u_iv ** 2 - u0_iv ** 2))
            ph_iv = PI + 2 * PI * u_iv + th_iv
            lr, li = log1p_polar(r_iv, ph_iv, r_ub)
            pr, _ = cmul(lr, li, w2r, w2i)
            i2_re += pr * du
            x_f_lo = G['c_re5'] * tau_lo * (au_lo ** 2 - u0_hi ** 2)
            r_f_ub = iv.exp(-x_f_lo)
            dz_e = dzeta_err_bound(G, au_lo, au_hi, u0_lo, u0_hi, tau_hi, Del0_hi, r0ck2_hi, False)
            ph_a = (PI + 2 * PI * ua) - iv.mpf(dz_e.b)
            ph_b = (PI + 2 * PI * ub) + iv.mpf(th_iv.b) + iv.mpf(dz_e.b)
            k_min = int(round((0.5 * (float(ph_a.a) + float(ph_b.b)) / fpi - 1.0) / 2.0))
            dist_odd = PI
            for kk in (k_min - 1, k_min, k_min + 1):
                target = (2 * kk + 1) * PI
                if ph_a.a <= target.b and target.a <= ph_b.b:
                    dist_odd = iv.mpf(0)
                else:
                    dist_odd = iv_min(dist_odd, iv_min(abs(ph_a - target), abs(ph_b - target)))
            d_ang = iv.mpf(1) if dist_odd.a >= (PI / 2).b else (iv.sin(dist_odd) if dist_odd.a > mp.mpf('0.03') else iv.mpf(0))
            d_min = iv.mpf(max((1 - r_f_ub).a, d_ang.a))
            assert d_min.a > mp.mpf('0.005'), f"W2 outside d_min={d_min} at u={ua}..{ub}"
            eL2 += (r_f_ub / d_min) * dz_e * W2_abs * du
            eW2 += (r_f_ub / d_min) * rel_W2 * W2_abs * du
            m3_int += (r_f_ub / d_min) * (1 + rel_W2) * W2_abs * du
        else:
            dist_s = PI - iv_max_hi(iv_abs_hi(ap), iv_abs_hi(am)) - iv.mpf('0.08')
            assert dist_s.a > mp.mpf('0.03')
            d_min = iv.sin(iv_min(dist_s, iv.mpf('1.5')))
            # Include BOTH finite-n and skipped parabola L1*W2 in eW2 ((2 + rel_W2), #106 fix #4!):
            eW2 += (1 / d_min) * W2_abs * (2 + rel_W2) * du
            m3_int += (1 / d_min) * W2_abs * (1 + rel_W2) * du

    r_L = iv.exp(-G['c_re5'] * tau_lo * (L_iv ** 2 - u0_hi ** 2))
    I_tail_5 = r_L / (2 * G['c_re5'] * tau_lo * L_iv)
    r_L_par = iv.exp(-tau_lo * (L_iv ** 2 - u0_hi ** 2))
    I_tail_par = r_L_par / (2 * tau_lo * L_iv)
    W2_par_out = 3 * inv_quad_ub(tau_lo, tau_hi, L_iv, 2)
    tail_2_fin = 2 * (1 / (1 - r_L)) * I_tail_5 * W2_fin_out + (1 / (1 - r_T5)) * far_r5 * W2_fin_all
    tail_2_par = 2 * (1 / (1 - r_L_par)) * I_tail_par * W2_par_out
    err_all = err_B1 + eL2 + eW2 + tail_2_fin + tail_2_par
    neg_Re_P2 = H_f - b1_re + i2_re + iv_sym(err_all)
    M3_ub = b1_im_coeff + m3_int + tail_2_fin
    return neg_Re_P2, M3_ub, H_f, (err_all / iv.mpf(H_f.a))

def eval_peak_curv_iv(ta, tb, u0_iv, th_iv, ct_lo, ct_hi=None):
    """Exact hull of peak curvature over tau in [ta, tb], splitting wide tau boxes (>10% rel width) to kill dependency."""
    rel_w = (float(tb.b) - float(ta.a)) / max(float(ta.a), 0.5)
    n_t = 1 if rel_w <= 0.105 else (2 if rel_w <= 0.22 else 3)
    if n_t == 1:
        return _eval_peak_curv_single_iv(ta, tb, u0_iv, th_iv, ct_lo, ct_hi=ct_hi)
    p_lo = None; p_hi = None; m3_hi = None; h_lo = None; h_hi = None; pr_hi = None
    for k in range(n_t):
        t_a = iv.mpf((ta + (tb - ta) * k / n_t).a)
        t_b = iv.mpf((ta + (tb - ta) * (k + 1) / n_t).b)
        np2, m3, hf, pr = _eval_peak_curv_single_iv(t_a, t_b, u0_iv, th_iv, ct_lo, ct_hi=ct_hi)
        p_lo = np2.a if p_lo is None else min(p_lo, np2.a)
        p_hi = np2.b if p_hi is None else max(p_hi, np2.b)
        m3_hi = m3.b if m3_hi is None else max(m3_hi, m3.b)
        h_lo = hf.a if h_lo is None else min(h_lo, hf.a)
        h_hi = hf.b if h_hi is None else max(h_hi, hf.b)
        pr_hi = pr.b if pr_hi is None else max(pr_hi, pr.b)
    return iv.mpf([p_lo, p_hi]), iv.mpf(m3_hi), iv.mpf([h_lo, h_hi]), iv.mpf(pr_hi)

# ==============================================================================
# SUB-LEMMA G3.3: OPTIMAL QUADRATIC POLE HEIGHT & POISSON PEAK ERROR E_pk
# ==============================================================================
def certify_E_pk_optimal(ta, tb, u0_iv, ct_lo, th_core, ct_hi=None):
    G, _, r0ck2_hi, Del0_hi = make_params(ta, tb, ct_lo, ct_hi)
    u0_lo = iv.mpf(u0_iv.a); u0_hi = iv.mpf(u0_iv.b)
    # Certify lam_* in [0.05, 3.00] on this box (justifying Pigeonhole-Annulus Lemma!):
    log_inv_lam_hi = G['te_hi'] * tb * iv.mpf(max(mp.mpf(0), (iv.mpf('0.25') - u0_lo ** 2).b))
    log_lam_hi     = G['te_hi'] * tb * iv.mpf(max(mp.mpf(0), (u0_hi ** 2 - iv.mpf('0.25')).b))
    assert log_inv_lam_hi.b < iv.log(iv.mpf(20)).a, f"lam_* too small: exp(-{log_inv_lam_hi})"
    assert log_lam_hi.b < iv.log(iv.mpf(3)).a, f"lam_* too large: exp(+{log_lam_hi})"

    E_higher = cert_E_higher_D1(G, u0_lo, u0_hi, Del0_hi, ct_lo / tb, ct_lo * ta, ct_lo)
    E_line_val = E_higher + NB_err_D + Y_big_err_D
    del_f = delta_fin_iv(G, u0_iv, Del0_hi, r0ck2_hi)
    ad_hi = iv_abs_hi(del_f)
    ad_lo = iv.mpf(abs(del_f).a)
    th_hit_lo = PI * (1 - 2 * ad_hi)
    Del_th = th_hit_lo - th_core
    assert Del_th.a > mp.mpf('0.12'), f"Peak too close to hit: Del_th={Del_th}"
    A = 2 * PI * (1 + Del0_hi)
    tau_iv = iv.mpf([ta.a, tb.b])
    K = G['c_rein'] * tau_iv
    C = 2 * G['te_hi'] * tau_iv * u0_iv
    KC = K * C
    B_q = 2 * KC * u0_iv + A ** 2
    def exp_in_for_dist(dth):
        disc_hi = iv.mpf(max(mp.mpf(0), (iv.mpf(B_q.b) ** 2 - 4 * iv.mpf(KC.a) * iv.mpf(A.a) * iv.mpf(dth.a)).b))
        d_lo = iv.mpf((2 * iv.mpf(A.a) * iv.mpf(dth.a) / (iv.mpf(B_q.b) + iv.sqrt(disc_hi))).a)
        d_c = iv.mpf(min(d_lo.a, u0_iv.a))
        assert d_c.a > 0, f"d_c <= 0: {d_c}"
        y_p_lo = iv.mpf((iv.mpf(K.a) * (2 * iv.mpf(u0_iv.a) * d_c - d_c ** 2) / iv.mpf(A.b)).a)
        return 2 * PI * ct_lo * y_p_lo
    def exp_out_for_dist(dth):
        return solve_out_neg_E_D1(G, u0_lo, u0_hi, Del0_hi, ct_lo / tb, ct_lo * ta, ct_lo, dth)
    def exp_for_dist(dth):
        return iv_min(exp_in_for_dist(dth), exp_out_for_dist(dth))
    exp_pole = exp_for_dist(Del_th)
    assert exp_pole.a > mp.mpf('1.5'), f"exp_pole too small: {float(exp_pole.a)}"
    # All 4 m = 1 boundary poles (2 inner y>0 + 2 outer y<0) with distances >= (th_hit-t, th_hit, pi-t, pi)! (#77 fix 2b!)
    def E_pk_at_th(th_hi, ct_val):
        sc = ct_val / ct_lo
        e1 = exp_for_dist(th_hit_lo - th_hi) * sc
        e2 = exp_for_dist(th_hit_lo) * sc
        e3 = exp_for_dist(PI - th_hi) * sc
        e4 = exp_for_dist(PI) * sc
        return (-iv.log(1 - iv.exp(-e1)) - iv.log(1 - iv.exp(-e2))
                - iv.log(1 - iv.exp(-e3)) - iv.log(1 - iv.exp(-e4)) + E_line_val)
    # Certified 1-line monotonicity 2-pole + 2-opposite-half-plane-pole valley factor (ZERO hardcoded 4.10!!):
    e_opp_0 = exp_for_dist(PI)
    opp_fac_0 = (1 + iv.exp(-e_opp_0)) ** 2
    sep = iv_min(4 * PI * ad_lo, 2 * PI - 4 * PI * ad_hi)
    V_4safe = 4 * opp_fac_0 * iv.exp(E_line_val)
    if sep.a >= mp.mpf('0.25'):
        e_half_0 = exp_for_dist(sep / 2)
        V_safe_box = iv_min_hi(2 * (1 + iv.exp(-e_half_0)) * opp_fac_0 * iv.exp(E_line_val), V_4safe)
    else:
        V_safe_box = V_4safe
    def V_bin_at(tj_a, tj_b, ct_val):
        sc = ct_val / ct_lo
        opp_f = (1 + iv.exp(-e_opp_0 * sc)) ** 2
        V_cap = 4 * opp_f * iv.exp(E_line_val)
        if sep.a >= mp.mpf('0.25'):
            e_h = exp_for_dist(sep / 2) * sc
            V_cap = iv_min_hi(2 * (1 + iv.exp(-e_h)) * opp_f * iv.exp(E_line_val), V_cap)
        if tj_b.b < th_hit_lo.a - mp.mpf('0.08'):
            d1 = th_hit_lo - tj_b
            d2 = iv_min(th_hit_lo + tj_a, PI - tj_b)
            e1 = exp_for_dist(d1) * sc
            e2 = exp_for_dist(d2) * sc
            return iv_min_hi((1 + iv.exp(-e1)) * (1 + iv.exp(-e2)) * opp_f * iv.exp(E_line_val), V_cap)
        return V_cap
    E_edge = E_pk_at_th(th_core, ct_lo)
    return E_edge, exp_pole, E_pk_at_th, V_bin_at, V_safe_box

# ==============================================================================
# MASTER BOX CERTIFIER FOR CASE G3-D
# ==============================================================================
def certify_G3D_box(ta, tb, u0_iv, ct_lo, ct_hi, g1_max_iv, th_core_val, label, exact_saddle=False):
    assert u0_iv.b <= mp.mpf('0.52'), f"box u0 > 0.52: {u0_iv}"
    th_core = iv.mpf(str(th_core_val))

    vl_drops = []
    N_vl = 10 if ct_hi is None else 12
    d_vl = (PI - th_core) / N_vl
    w1_cache = _prepare_Psi_W1_cache(ta, tb, u0_iv, ct_lo, ct_hi=ct_hi, allow_cub=(ct_hi is not None))
    for j in range(N_vl):
        tja = th_core + d_vl * j; tjb = th_core + d_vl * (j + 1)
        tj = iv.mpf([tja.a, tjb.b])
        dr_lo = eval_valley_drop_iv(ta, tb, u0_iv, tj, ct_lo, ct_hi=ct_hi, cache=w1_cache)
        vl_drops.append((tja, tjb, d_vl, dr_lo))
    # Lazy exact v in [0, u0] valley drop floor for small u0 (tau <= 2.85, ct >= 40) ONLY when W1 drop < 0.015:
    if tb.b <= mp.mpf('2.85') and u0_iv.b <= mp.mpf('0.30') and ct_lo.a >= mp.mpf('40.0'):
        if min(dr.a for _, _, _, dr in vl_drops) < mp.mpf('0.015'):
            G_vl, _, r0ck2_vl, Del0_vl = make_params(ta, tb, ct_lo, ct_hi)
            del_f_vl = delta_fin_iv(G_vl, u0_iv, Del0_vl, r0ck2_vl)
            th_hit_hi_vl = PI * (1 - 2 * iv.mpf(abs(del_f_vl).a))
            v0_drop_floor, v0_curv_lo = v0_valley_drop_floor_iv(ta, tb, iv.mpf(u0_iv.a), iv.mpf(u0_iv.b), ct_lo)
            if v0_drop_floor.a > 0:
                new_vl = []
                for tja, tjb, dth, dr_lo in vl_drops:
                    if tja.a >= th_hit_hi_vl.b:
                        dr_lo = iv.mpf(max(dr_lo.a, v0_drop_floor.a))
                    elif v0_curv_lo.a > 0:
                        dr_lo = iv.mpf(max(dr_lo.a, min((iv.mpf('0.5') * v0_curv_lo * (tja ** 2)).a, v0_drop_floor.a)))
                    new_vl.append((tja, tjb, dth, dr_lo))
                vl_drops = new_vl

    E_pk, exp_pole, E_pk_at, V_bin_at, V_safe_box = certify_E_pk_optimal(ta, tb, u0_iv, ct_lo, th_core, ct_hi=ct_hi)

    if ct_hi is not None:
        g1_box = iv.mpf(0) if exact_saddle else eval_g1_fin_iv(ta, tb, u0_iv, ct_lo, ct_hi=ct_hi)
        g1_abs_hi = iv_abs_hi(g1_box)
        N_tot = 15; h = th_core / N_tot
        D_lo_arr = []; D_hi_arr = []; Phi_arr = []; max_pert = mp.mpf(0)
        for k in range(N_tot):
            tk = iv.mpf([(h * k).a, (h * (k + 1)).b])
            d_lo, d_hi, im_hi, err_w = eval_Psi_diff_W1_iv(ta, tb, u0_iv, tk, ct_lo, need_im=True, g1_iv=g1_box, ct_hi=ct_hi, cache=w1_cache)
            D_lo_arr.append(d_lo); D_hi_arr.append(d_hi); Phi_arr.append(im_hi)
            if err_w.b > max_pert: max_pert = err_w.b
        # Also intersect on every valid prefix bin with exact 2nd-derivative curvature + O(th^3) M_3 phase integral!
        # (Exact conjugation symmetry Psi(-th) = conj(Psi(th)) gives Re Psi'(0) == 0 and Im Psi''(0) == 0 for EVERY u_0!)
        N_pk_a = 3; dt_pka = th_core / N_pk_a
        cl_a = []; ch_a = []; m3_a = []
        for j in range(N_pk_a):
            tj = iv.mpf([(dt_pka * j).a, (dt_pka * (j + 1)).b])
            try:
                nReP2, M3_j, _, _ = eval_peak_curv_iv(ta, tb, u0_iv, tj, ct_lo, ct_hi=ct_hi)
            except AssertionError:
                break
            if nReP2.a <= 0:
                break
            cl_a.append(iv.mpf(nReP2.a)); ch_a.append(iv.mpf(nReP2.b)); m3_a.append(iv.mpf(M3_j.b))
        n_ok = len(cl_a)
        if n_ok > 0:
            p_lo = iv.mpf(0); p_hi = iv.mpf(0); dp_lo = iv.mpf(0); dp_hi = iv.mpf(0)
            im_c = iv.mpf(0); dim_c = iv.mpf(0)
            for k in range(N_tot):
                jb = (k * N_pk_a) // N_tot
                if jb >= n_ok:
                    break
                cl = cl_a[jb]; ch = ch_a[jb]; m3 = m3_a[jb]
                tk_a = h * k; tk_b = h * (k + 1)
                p_lo_next = p_lo + dp_lo * h + iv.mpf('0.5') * cl * (h ** 2)
                dp_lo = dp_lo + cl * h
                p_hi_next = p_hi + dp_hi * h + iv.mpf('0.5') * ch * (h ** 2)
                dp_hi = dp_hi + ch * h
                im_next = im_c + dim_c * h + m3 * ((tk_b ** 3 - tk_a ** 3) / 6 - iv.mpf('0.5') * (tk_a ** 2) * h)
                dim_c = dim_c + iv.mpf('0.5') * m3 * (tk_b ** 2 - tk_a ** 2)
                # Intersect rigorously with W1 bounds on bin k:
                D_lo_arr[k] = iv.mpf(max(D_lo_arr[k].a, p_lo.a))
                D_hi_arr[k] = iv.mpf(min(D_hi_arr[k].b, p_hi_next.b))
                assert D_lo_arr[k].a <= D_hi_arr[k].b, f"curv/W1 drop intersection empty on bin {k}: [{D_lo_arr[k].a}, {D_hi_arr[k].b}]"
                Phi_arr[k] = iv_min_hi(Phi_arr[k], g1_abs_hi * tk_b + im_next)
                p_lo = p_lo_next; p_hi = p_hi_next; im_c = im_next
        min_rat = mp.mpf('1e9'); worst_ct = float(ct_lo.a); worst_diag = None
        N_ct = 24
        r_geo = iv.mpf((ct_hi / ct_lo) ** (iv.mpf(1) / N_ct))
        for m in range(N_ct):
            ca = ct_lo * (r_geo ** m)
            cb = ct_lo * (r_geo ** (m + 1)) if m < N_ct - 1 else ct_hi
            pk_re = iv.mpf(0); pk_abs = iv.mpf(0)
            for k in range(N_tot):
                assert D_lo_arr[k].a >= 0, f"{label} negative D_lo_arr[{k}]={D_lo_arr[k]}"
                E_k = E_pk_at(h * (k + 1), ca)
                phi_k = E_k + cb * Phi_arr[k]
                e_hi = iv.exp(E_k) * iv.exp(-ca * D_lo_arr[k])
                pk_abs += iv.mpf(e_hi.b) * h
                if phi_k.b < (PI / 2).a:
                    e_lo = iv.exp(-E_k) * iv.exp(-cb * D_hi_arr[k]) * iv.cos(phi_k)
                    pk_re += iv.mpf(e_lo.a) * h
                elif phi_k.b < PI.a:
                    # Exact lower bound cos(phi) >= cos(phi_max) in (-1, 0] when phi in [pi/2, pi):
                    pk_re += iv.mpf((e_hi * iv.cos(iv.mpf(phi_k.b))).a) * h
                else:
                    pk_re -= iv.mpf(e_hi.b) * h
            vl_abs = iv.mpf(0)
            for tja, tjb, dth, dr_lo in vl_drops:
                assert dr_lo.a >= 0, f"{label} negative valley drop dr_lo={dr_lo}"
                vl_abs += V_bin_at(tja, tjb, ca) * dth * iv.exp(-ca * dr_lo)
            rat_m = (pk_re - vl_abs) / (pk_abs + vl_abs) - E_cubic_G3_sub
            if rat_m.a < min_rat:
                min_rat = rat_m.a; worst_ct = float(ca.a)
                worst_diag = (float(pk_re.a), float(pk_abs.b), float(vl_abs.b), float(g1_box.a), float(g1_box.b))
        min_dr = min(float(dr.a) for _, _, _, dr in vl_drops)
        print(f"  [{label:30s}] ct=[{float(ct_lo.a):4.0f},{float(ct_hi.b):4.0f}] | errW1<={float(max_pert):.4f} | E_pk<={float(E_pk.b):.4f} (ep={float(exp_pole.a):.2f}) V<={float(V_safe_box.b):.2f} | min_vl={min_dr:.4f} | rat>={float(min_rat):+.4f} (worst ct={worst_ct:.1f})", flush=True)
        if float(min_rat) <= 0.05:
            print(f"    DIAG: pk_re={worst_diag[0]:+.4f} pk_abs={worst_diag[1]:.4f} vl_abs={worst_diag[2]:.4f} g1=[{worst_diag[3]:+.4f},{worst_diag[4]:+.4f}]", flush=True)
        assert float(min_rat) > 0.05, f"FAILED {label}: min_rat = {float(min_rat)}"
        return float(min_rat)

    # UNBOUNDED TIER B (ct in [ct_lo, inf), exact saddle g1_max_iv == 0, pure O(th^3) phase!):
    N_pk = 3; dt_pk = th_core / N_pk
    c_lo_list = []; c_hi_list = []; m3_hi_list = []; H_f_ref = None; max_pert = mp.mpf(0)
    for j in range(N_pk):
        tj = iv.mpf([(dt_pk * j).a, (dt_pk * (j + 1)).b])
        nReP2, M3_j, H_f, pert_rel = eval_peak_curv_iv(ta, tb, u0_iv, tj, ct_lo)
        if H_f_ref is None: H_f_ref = H_f
        assert nReP2.a > 0, f"{label} non-positive peak curv on bin {j}: {nReP2}"
        c_lo_list.append(iv.mpf(nReP2.a))
        c_hi_list.append(iv.mpf(nReP2.b))
        m3_hi_list.append(iv.mpf(M3_j.b))
        if pert_rel.b > max_pert: max_pert = pert_rel.b

    N_sub = 8; N_tot = N_pk * N_sub; h = th_core / N_tot
    D_lo_arr = []; D_hi_arr = []; Phi_arr = []
    p_lo = iv.mpf(0); p_hi = iv.mpf(0); dp_lo = iv.mpf(0); dp_hi = iv.mpf(0)
    im_hi = iv.mpf(0); dim_hi = iv.mpf(0)
    for k in range(N_tot):
        jb = k // N_sub
        cl = c_lo_list[jb]; ch = c_hi_list[jb]; m3 = m3_hi_list[jb]
        tk_a = h * k; tk_b = h * (k + 1)
        D_lo_arr.append(p_lo)
        p_lo_next = p_lo + dp_lo * h + iv.mpf('0.5') * cl * (h ** 2)
        dp_lo = dp_lo + cl * h
        p_hi_next = p_hi + dp_hi * h + iv.mpf('0.5') * ch * (h ** 2)
        dp_hi = dp_hi + ch * h
        D_hi_arr.append(p_hi_next)
        im_next = im_hi + dim_hi * h + m3 * ((tk_b ** 3 - tk_a ** 3) / 6 - iv.mpf('0.5') * (tk_a ** 2) * h)
        dim_hi = dim_hi + iv.mpf('0.5') * m3 * (tk_b ** 2 - tk_a ** 2)
        Phi_arr.append(im_next)
        p_lo = p_lo_next; p_hi = p_hi_next; im_hi = im_next

    min_rat = mp.mpf('1e9'); worst_ct = float(ct_lo.a)
    dr_vl_min = iv.mpf(min(dr.a for _, _, _, dr in vl_drops))
    assert dr_vl_min.a > 0, f"{label} non-positive valley drop: {dr_vl_min}"
    ct_cap = iv.mpf(max((ct_lo * iv.mpf('1.1')).b, (iv.mpf('30.0') / iv.mpf(H_f_ref.a)).b, (iv.mpf('1.0') / dr_vl_min).b))
    if ct_lo.a < ct_cap.b:
        N_ct = 32
        r_geo = iv.mpf((ct_cap / ct_lo) ** (iv.mpf(1) / N_ct))
        for m in range(N_ct):
            ca = ct_lo * (r_geo ** m)
            cb = ct_lo * (r_geo ** (m + 1)) if m < N_ct - 1 else ct_cap
            pk_re = iv.mpf(0); pk_abs = iv.mpf(0)
            for k in range(N_tot):
                t_hi = h * (k + 1)
                E_k = E_pk_at(t_hi, ca)
                phi_k = E_k + cb * (Phi_arr[k] + g1_max_iv * t_hi)
                e_hi = iv.exp(E_k) * iv.exp(-ca * D_lo_arr[k])
                pk_abs += iv.mpf(e_hi.b) * h
                if phi_k.b < (PI / 2).a:
                    e_lo = iv.exp(-E_k) * iv.exp(-cb * D_hi_arr[k]) * iv.cos(phi_k)
                    pk_re += iv.mpf(e_lo.a) * h
                elif phi_k.b < PI.a:
                    # Exact lower bound cos(phi) >= cos(phi_max) in (-1, 0] when phi in [pi/2, pi):
                    pk_re += iv.mpf((e_hi * iv.cos(iv.mpf(phi_k.b))).a) * h
                else:
                    pk_re -= iv.mpf(e_hi.b) * h
            vl_abs = iv.mpf(0)
            for tja, tjb, dth, dr_lo in vl_drops:
                vl_abs += V_bin_at(tja, tjb, ca) * dth * iv.exp(-ca * dr_lo)
            rat_m = (pk_re - vl_abs) / (pk_abs + vl_abs) - E_cubic_G3_sub
            if rat_m.a < min_rat:
                min_rat = rat_m.a; worst_ct = float(ca.a)

    assert g1_max_iv.b == 0
    c_lo_all = iv.mpf(min(x.a for x in c_lo_list))
    c_hi_all = iv.mpf(max(x.b for x in c_hi_list))
    m3_all   = iv.mpf(max(x.b for x in m3_hi_list))
    Z_c = iv.mpf('1.55') / iv.sqrt(c_lo_all)
    assert (Z_c / iv.sqrt(ct_cap)).b <= th_core.a
    N_z = 40; dz = Z_c / N_z; L_pk = iv.mpf(0)
    for i in range(N_z):
        z_iv = iv.mpf([(dz * i).a, (dz * (i + 1)).b])
        phi_z = E_pk + (m3_all / (6 * iv.sqrt(ct_cap))) * (z_iv ** 3)
        if phi_z.b < (PI / 2).a:
            L_pk += iv.exp(-E_pk) * iv.mpf((iv.exp(-iv.mpf('0.5') * c_hi_all * (z_iv ** 2)) * iv.cos(phi_z)).a) * dz
        else:
            L_pk -= iv.exp(E_pk) * iv.mpf((iv.exp(-iv.mpf('0.5') * c_lo_all * (z_iv ** 2))).b) * dz
    U_pt = iv.exp(E_pk) * iv.exp(-iv.mpf('0.5') * c_lo_all * (Z_c ** 2)) / (c_lo_all * Z_c)
    assert ct_cap.a > (1 / (2 * dr_vl_min)).b
    U_vl = V_safe_box * (PI - th_core) * iv.sqrt(ct_cap) * iv.exp(-ct_cap * dr_vl_min)
    U_tot = iv.exp(E_pk) * iv.sqrt(PI / (2 * c_lo_all)) + U_vl
    rat_inf = (L_pk - U_pt - U_vl) / U_tot - E_cubic_G3_sub
    if rat_inf.a < min_rat:
        min_rat = rat_inf.a; worst_ct = float(ct_cap.a)

    min_dr = min(float(dr.a) for _, _, _, dr in vl_drops)
    print(f"  [{label:30s}] ct>={float(ct_lo.a):5.1f} | pert<={float(max_pert):.3f} | E_pk<={float(E_pk.b):.4f} (ep={float(exp_pole.a):.2f}) V<={float(V_safe_box.b):.2f} | min_vl={min_dr:.4f} | rat>={float(min_rat):+.4f} (worst ct={worst_ct:.1f})", flush=True)
    assert float(min_rat) > 0.05, f"FAILED {label}: min_rat = {float(min_rat)}"
    return float(min_rat)

# ==============================================================================
# HARMONICITY / CAUCHY HORIZONTAL-SADDLE HELPERS FOR SMALL-tau BANDS (D1, D2)
# ==============================================================================
def CW3_free_D1(G, c_min, Del0_u, T_loc):
    # Exact 1/(1-s_j) factor (#106 fix #1!):
    a_p = (1 + Del0_u) / (G['S'] * iv.cos(G['al']) * c_min)
    b_p = iv.sqrt((1 + 1 / c_min) * (2 + 1 / c_min)) / (G['S'] * iv.cos(G['al']))
    P_th = iv.sqrt(a_p ** 2 + (b_p ** 2) * T_loc)
    return (15 * ((2 * G['te_hi']) ** 3) + 40 * ((2 * G['te_hi']) ** 2) * P_th + 48 * G['te_hi'] * (P_th ** 2)) / 120

def far_T5_helpers_D1(G, c_min, Del0_u):
    """Exact T(u) > 0.05 far-tail helpers for D = +-1 (#106 fix #5a!)."""
    n_lo = iv.mpf(max(N_G3_D.a, G['n_min'].a))
    assert (2 / iv.log(iv.mpf('1.05') / (1 + G['T_in']))).b < n_lo.a
    far_r5 = n_lo * (((1 + G['T_in']) / iv.mpf('1.05')) ** (n_lo / 2))
    assert far_r5.b < mp.mpf('1e-6')
    Q_T5 = G['thc_1ms'] / (iv.sqrt(iv.mpf('0.05')) * (1 - 1 / c_min) * iv.cos(G['th_m'] / 2))
    C_W3_T5 = (15 * ((2 * G['te_hi']) ** 3) * (Q_T5 ** 2)
               + 40 * ((2 * G['te_hi']) ** 2) * Q_T5 * G['thc_1ms']
               + 24 * (2 * G['te_hi']) * (G['thc_1ms'] ** 2))
    return far_r5, C_W3_T5

def cert_B1_B2_pert_D1(G, c_min, u0_hi, Del0_u, r0ck2_u):
    b_ub = (1 + Del0_u) / (2 * c_min)
    d_re = (G['te_hi'] - 1) + G['te_hi'] * b_ub + G['te_hi'] * (1 + b_ub) * ((2 * u0_hi * G['thc_hi']) ** 2) / 6
    d_im = iv_max_hi(Del0_u, G['te_hi'] * (1 + b_ub) * u0_hi * G['thc_hi'])
    ep_P = iv.sqrt(d_re ** 2 + d_im ** 2)
    inv2 = 1 / ((1 - ep_P) ** 2)
    rel_W1 = ((G['te_hi'] - 1) + 2 * u0_hi * G['thc_hi']) * inv2 + (inv2 - 1)
    d_im_z = iv_max_hi(Del0_u, G['C_D2'] * r0ck2_u * (u0_hi ** 2))
    ep_dz = iv.sqrt(((G['te_hi'] - G['r_re0']) ** 2) + (d_im_z ** 2))
    rel_B1 = (1 + rel_W1) / (1 - ep_dz) - 1
    a_p0 = (1 + Del0_u) / (G['S'] * iv.cos(G['al']) * c_min)
    Pth_2t = iv.sqrt(a_p0 ** 2 + (u0_hi * G['thc_1ms']) ** 2) * (1 + ep_P)
    mult_W2 = ((G['te_hi'] ** 2) + (iv.mpf(4) / 3) * G['te_hi'] * Pth_2t) / ((1 - ep_P) ** 4)
    return rel_B1, mult_W2 / (1 - ep_dz)

def exact_R3_W3_D1(G, ta, tb, u_max, ct_min, Del0):
    """Exact 50-digit iv upper bound on |R_{3,fin}| / (tau / (2*pi^2)) for ALL u0 in [0, u_max], |th| < th_hit (#106 fixes #1, #2, #5a!)."""
    CW3_in  = CW3_free_D1(G, ct_min, Del0, G['T_in'])
    CW3_out = CW3_free_D1(G, ct_min, Del0, iv.mpf('0.05'))
    f_p_lo = (1 - 1 / ct_min) * (1 - 1 / (2 * ct_min))
    mob_in  = f_p_lo * (G['S'] ** 4) * G['te_lo'] * G['sinc_in_sq']
    mob_out = f_p_lo * (G['S'] ** 4) * G['te_lo'] * G['sinc_c5_sq']
    d_lo = 1 - Del0
    far_r5, C_W3_T5 = far_T5_helpers_D1(G, ct_min, Del0)
    far_R3 = ((PI ** 2) / 6) * far_r5 * C_W3_T5 * (2 * (PI ** 2)) / ((2 * PI * d_lo) ** 4)
    L = iv.mpf('4.0')
    assert (G['gam'] * (L ** 2)).b <= mp.mpf('0.05') and (L * G['thc_hi']).b <= mp.mpf('0.080')
    mv = mp.mpf(0)
    for m in range(4):
        t1 = ta + (tb - ta) * m / 4; t2 = ta + (tb - ta) * (m + 1) / 4
        N = 60
        s_in = iv.mpf(0); du1 = u_max / 12
        for k in range(12):
            ua = du1 * k
            s_in += CW3_in * ((PI ** 2) / 6) / ((mob_in * ((t1 * ua) ** 2) + (d_lo * PI) ** 2) ** 3) * du1
        s_out = iv.mpf(0); du2 = (L - u_max) / N
        for k in range(N):
            ua = u_max + du2 * k
            rf = iv.exp(-G['c_re5'] * t1 * (ua ** 2 - u_max ** 2))
            s_out += CW3_out * li2_pos_ub(rf) / ((mob_out * ((t1 * ua) ** 2) + (d_lo * PI) ** 2) ** 3) * du2
        tail = CW3_out * ((PI ** 2) / 6) * iv.exp(-G['c_re5'] * t1 * (L ** 2 - u_max ** 2)) / (2 * G['c_re5'] * t1 * L) / ((mob_out * ((t1 * L) ** 2) + (d_lo * PI) ** 2) ** 3)
        val = (iv.mpf(15) / 2) * (PI ** 2) * (t2 ** 2) * (s_in + s_out + tail) + far_R3
        if val.b > mv: mv = val.b
    return iv.mpf(mv)

def curv_W3_small_u_iv(ta, tb, u_max, ct_lo, th_factor_pi=iv.mpf(2), R3_pre=None):
    """50-digit iv enclosure [c_lo, c_hi] of -Re Psi_fin''(th) / (tau/(2*pi^2)) on u0 in [0, u_max], |th| <= th_factor_pi*pi*u0."""
    assert u_max.b <= mp.mpf('0.52'), f"curv_W3 u_max > 0.52: {u_max}"
    G, c_min, r0ck2, Del0 = make_params(ta, tb, ct_lo)
    R3 = exact_R3_W3_D1(G, ta, tb, u_max, ct_lo, Del0) if R3_pre is None else R3_pre
    p_lo = mp.mpf('10'); p_hi = mp.mpf('0')
    for i in range(20):
        ua = u_max * i / 20; ub = u_max * (i + 1) / 20
        u_i = iv.mpf([ua.a, ub.b]); s_i = iv.mpf([mp.mpf(0), (tb * ub / PI).b])
        val_i = (2 * u_i) / (1 + s_i ** 2) + (1 - 2 * u_i) * (1 - 3 * (s_i ** 2)) / ((1 + s_i ** 2) ** 3)
        if val_i.a < p_lo: p_lo = val_i.a
        if val_i.b > p_hi: p_hi = val_i.b
    rel_B1, mult_B2 = cert_B1_B2_pert_D1(G, c_min, u_max, Del0, r0ck2)
    err_HB1 = (1 - G['r_re0']) + (G['te_hi'] - 1) + rel_B1 + 2 * Del0
    if th_factor_pi.b <= mp.mpf('1.5'):
        # For |th| <= th_factor_pi*pi*u_0, |ap|, |am| <= pi*(1 - (2 - th_factor_pi)*u_0 + 2*Del0):
        b2_fac = mp.mpf(0)
        rem_f = 2 - th_factor_pi
        for i in range(20):
            ua = u_max * i / 20; ub = u_max * (i + 1) / 20
            ui = iv.mpf([ua.a, ub.b])
            pos_f = iv.mpf('1.5') * ((1 - rem_f * ui + 2 * Del0) ** 2) - iv.mpf('0.5')
            f_hi = iv_max_hi(iv.mpf('0.5'), iv_abs_hi(pos_f))
            v = ub * f_hi
            if v.b > b2_fac: b2_fac = v.b
        ReB2 = mult_B2 * (iv.mpf(5) / (2 * (PI ** 2))) * (tb ** 2) * iv.mpf(b2_fac)
    else:
        ReB2 = mult_B2 * (iv.mpf(5) / (2 * (PI ** 2))) * (tb ** 2) * u_max
    c_lo = iv.mpf(p_lo) - err_HB1 - ReB2 - R3
    c_hi = iv.mpf(p_hi) + err_HB1 + ReB2 + R3
    K3 = (1 + rel_B1) * (3 * tb * u_max / (PI ** 2)) + mult_B2 * (3 * tb / (2 * (PI ** 2)))
    return G, Del0, R3, rel_B1, mult_B2, c_lo, c_hi, K3

def cauchy_horiz_saddle_iv(G, ta, tb, ct_lo, Del0):
    """Exact 50-digit iv Cauchy shift to Im u = y_s (zero side-lines by 2L_n-periodicity, zero interior poles!):"""
    dw = (1 + Del0) / (2 * ct_lo)
    a_s_lo = G['te_lo'] * (1 - dw) / ((1 + dw) ** 2)
    a_s_hi = G['te_hi'] * (1 + dw) / ((1 - dw) ** 2)
    exp_lo = (PI ** 2) * ((1 - Del0) ** 2) / (a_s_hi * tb)
    exp_hi = (PI ** 2) * ((1 + Del0) ** 2) / (a_s_lo * ta)
    assert (PI * (1 + Del0) / (a_s_lo * ta)).b < (iv.mpf('0.5') * G['y_br']).a
    X_s = iv.mpf('2.0')
    kap_hi = tb / ct_lo
    phi_s_max = (iv.mpf(4) / 3) * a_s_hi * ((kap_hi / iv.sqrt(ta)) / (2 * PI * G['S'] * iv.cos(G['al']) * (1 - dw))) * (X_s ** 3)
    assert phi_s_max.b < mp.mpf('0.35')
    I1_core = iv.mpf(0); dt_s = X_s / 40
    for j in range(40): I1_core += 2 * iv.exp(-a_s_hi * ((dt_s * (j + 1)) ** 2)) * dt_s
    c_re5_s = a_s_lo * (G['c_re5'] / G['te_lo'])
    # Rigorous far-x tail where T_s(x) > 0.05 (|x| <= L_n = n/(2c) <= n/2): bounded by n_min * (1.05)^(-n_min/2)!
    n_lo = iv.mpf(max(N_G3_D.a, G['n_min'].a))
    far_s_err = n_lo * (iv.mpf('1.05') ** (-n_lo / 2))
    assert far_s_err.b < mp.mpf('1e-8')
    I1_tail = 2 * iv.exp(-c_re5_s * (X_s ** 2)) / (2 * c_re5_s * X_s) + far_s_err
    c_I1_lo = (iv.cos(phi_s_max) * I1_core - I1_tail) / iv.sqrt(PI)
    assert c_I1_lo.a > mp.mpf('0.70')
    I1_lo = c_I1_lo * iv.sqrt(PI / tb) * iv.exp(-exp_hi)
    I1_hi = (iv.sqrt(PI / (c_re5_s * ta)) + far_s_err) * iv.exp(-exp_lo)
    rho_I = (1 / (iv.sqrt(2 * c_re5_s) * c_I1_lo) + far_s_err) * iv.exp(-exp_lo) / (1 - iv.exp(-exp_lo))
    assert rho_I.b < mp.mpf('0.5'), f"cauchy rho_I too large: {rho_I}"
    c_v0_ratio = 2 * (1 + iv.mpf('0.5') * rho_I) / (1 - rho_I)
    rem_hi = (iv.sqrt(PI / (2 * c_re5_s * ta)) + far_s_err) * iv.exp(-2 * exp_lo) / (1 - iv.exp(-exp_lo))
    assert (I1_lo - rem_hi).a > 0, f"cauchy I1_lo <= rem_hi: {I1_lo} vs {rem_hi}"
    res_err = Del0 / (1 - Del0)
    neg_g1_lo = I1_lo - rem_hi - res_err
    neg_g1_hi = I1_hi + rem_hi + res_err
    return neg_g1_lo, neg_g1_hi, c_v0_ratio, exp_lo

def v0_valley_drop_floor_iv(ta, tb, u0_lo, u0_hi, ct_lo):
    """Certified exact v in [0, u0] valley drop floor >= c_drop * ta * u0_lo^2 (for th >= th_hit) and universal |th| < th_hit curvature."""
    G, Del0, R3_v0, rel_B1, mult_B2, c_lo_all, c_hi_all, _ = curv_W3_small_u_iv(ta, tb, u0_hi, ct_lo, th_factor_pi=iv.mpf(2))
    _, _, _, _, _, _, c_hi_0, _ = curv_W3_small_u_iv(ta, tb, u0_hi, ct_lo, th_factor_pi=iv.mpf('0.5'), R3_pre=R3_v0)
    _, _, c_v0_ratio, _ = cauchy_horiz_saddle_iv(G, ta, tb, ct_lo, Del0)
    s_max = tb * u0_hi / PI
    d_lo = 1 - Del0
    far_r5, _ = far_T5_helpers_D1(G, ct_lo, Del0)
    far_W2 = ((PI ** 2) / 6) * far_r5 * (3 * ((2 * G['te_hi'] * tb) ** 2) / ((2 * PI * d_lo) ** 4) + 8 * G['te_hi'] * tb * G['thc_1ms'] / ((2 * PI * d_lo) ** 3))
    d_lam_th_ub = -(1 - (s_max ** 2) - 2 * Del0) + (1 + rel_B1) * (tb / (6 * PI)) + mult_B2 * ((tb ** 2) / (8 * (PI ** 2))) * (2 * u0_hi + iv.sqrt(PI / (G['c_re5'] * ta))) + far_W2
    # Include exact D = -1 residue term at u_0 = 0+ (#96 fix #3!):
    res_drop = c_v0_ratio * (Del0 / (1 - Del0)) / (ta * (u0_lo ** 2)) if u0_lo.a > 0 else iv.mpf('1e9')
    c_drop = (-d_lam_th_ub) * G['r_re0'] - c_v0_ratio * c_hi_0 * G['te_hi'] * tb / (2 * (PI ** 2)) - res_drop
    dr_fl = c_drop * ta * (u0_lo ** 2) if c_drop.a > 0 else iv.mpf(0)
    curv_lo = c_lo_all * ta / (2 * (PI ** 2)) if c_lo_all.a > 0 else iv.mpf(0)
    return dr_fl, curv_lo

def harmonic_u_star_bracket_iv(ta, tb, ct_lo, n_sub_t=4):
    """Certified exact [u_L, u_R] bracket of u_0^* via Cauchy horizontal-saddle g_1(0+) + sub-sliced Harmonicity integral."""
    uL_all = mp.mpf('10'); uR_all = mp.mpf('0')
    for k in range(n_sub_t):
        t1 = iv.mpf((ta + (tb - ta) * k / n_sub_t).a)
        t2 = iv.mpf((ta + (tb - ta) * (k + 1) / n_sub_t).b)
        G, _, _, Del0 = make_params(t1, t2, ct_lo)
        ng_lo, ng_hi, _, _ = cauchy_horiz_saddle_iv(G, t1, t2, ct_lo, Del0)
        assert ng_lo.a > 0, f"cauchy ng_lo <= 0 on [{t1},{t2}]"
        dv = iv.mpf('0.02')
        int_lo = iv.mpf(0); int_hi = iv.mpf(0)
        uL_k = None; uR_k = None
        for i in range(24):
            va = dv * i; vb = dv * (i + 1)
            if i < 4:
                _, _, _, _, _, cl_n, ch_n, _ = curv_W3_small_u_iv(t1, t2, vb, ct_lo, th_factor_pi=iv.mpf('0.05'))
                curv_lo = cl_n * t1 / (2 * (PI ** 2)); curv_hi = ch_n * t2 / (2 * (PI ** 2))
            else:
                nP2, _, _, _ = eval_peak_curv_iv(t1, t2, iv.mpf([va.a, vb.b]), iv.mpf(0), ct_lo)
                curv_lo = iv.mpf(nP2.a); curv_hi = iv.mpf(nP2.b)
            assert curv_lo.a > 0, f"non-positive curvature on v=[{va},{vb}]"
            dz_lo = G['r_re0'] * t1 * (vb ** 2 - va ** 2)
            dz_hi = G['te_hi'] * t2 * (vb ** 2 - va ** 2)
            d_lo = curv_lo * dz_lo; d_hi = curv_hi * dz_hi
            if uL_k is None and (int_hi + d_hi).b >= ng_lo.a:
                rem = iv.mpf(max(mp.mpf(0), (ng_lo - int_hi).a))
                uL_k = iv.sqrt(va ** 2 + rem / (curv_hi * G['te_hi'] * t2))
            int_lo += d_lo; int_hi += d_hi
            if uR_k is None and int_lo.a > ng_hi.b:
                rem = iv.mpf(max(mp.mpf(0), (ng_hi - (int_lo - d_lo)).b))
                uR_k = iv.sqrt(va ** 2 + rem / (curv_lo * G['r_re0'] * t1))
                break
        assert uL_k is not None and uR_k is not None, f"harmonic bracket failed on [{t1},{t2}]"
        if uL_k.a < uL_all: uL_all = uL_k.a
        if uR_k.b > uR_all: uR_all = uR_k.b
    return iv.mpf(uL_all), iv.mpf(uR_all)

def certify_G3D_harmonic_tail(ta, tb, ct_min, u_max, c_u, label):
    """
    100% self-contained 50-digit iv harmonicity + v=0 saddle tail certifier for ALL ct >= ct_min, tau in [ta, tb]!
    Chooses u_0 = max(u_{0,fin}^*, min(u_max, c_u / (ct * tau^3)^(1/4))) with f_core = 1.25 (so |th| <= 1.25*pi*u_0).
    Certifies u_0^* <= u_max via exact Cauchy shift, exact 4-pole E_pk, exact [th_use, PI] buffer, and exact D=-1 residue!
    """
    f_core = iv.mpf('1.25')
    ct_tau = ct_min * ta
    G, Del0, R3, rel_B1, mult_B2, c_lo, c_hi, K3 = curv_W3_small_u_iv(ta, tb, u_max, ct_min, th_factor_pi=f_core)
    # Separate universal curvature c_lo_buf valid on ALL |th| < th_hit (for the [th_use, th_hit] buffer!):
    _, _, _, _, _, c_lo_buf, _, _ = curv_W3_small_u_iv(ta, tb, u_max, ct_min, th_factor_pi=iv.mpf(2))
    assert c_lo.a > mp.mpf('0.35') and c_lo_buf.a > mp.mpf('0.20'), f"{label} c_lo={c_lo}, c_lo_buf={c_lo_buf}"
    # Certify lam_* in [0.05, 3.00] (u_0 in [0, u_max] < 0.5 => lam_* <= 1, 1/lam_* <= exp(te_hi*tb/4)):
    assert (G['te_hi'] * tb / 4).b < iv.log(iv.mpf(20)).a
    E_higher = cert_E_higher_D1(G, iv.mpf('0.001'), u_max, Del0, ct_min / tb, ct_tau, ct_min)
    V_BUMP = 4 * iv.exp(E_higher + NB_err_D + Y_big_err_D)
    ep_opp = solve_out_neg_E_D1(G, iv.mpf('0.0'), u_max, Del0, ct_min / tb, ct_tau, ct_min, PI)
    V_full = V_BUMP * ((1 + iv.exp(-ep_opp)) ** 2)
    neg_g1_lo, neg_g1_hi, c_v0_ratio, exp_s_lo = cauchy_horiz_saddle_iv(G, ta, tb, ct_min, Del0)
    # Certify u_0^* <= u_max (hence g_{1,fin}(u_max) > 0 and g_{1,fin}(u_0) >= 0 ALWAYS, #96 fix #3!):
    u_star_ub = iv.sqrt(2 * (PI ** 2) * neg_g1_hi / (c_lo * G['r_re0'] * (ta ** 2)))
    assert u_star_ub.b <= u_max.a, f"{label} u_star_ub={float(u_star_ub.b):.4f} > u_max={float(u_max.a):.4f}"
    # Exact harmonicity misalignment (including exact D = -1 residue +Del0/(1-Del0) at u_0 = 0+!):
    sq_g1 = (c_hi * G['te_hi'] * (c_u ** 2)) / (iv.sqrt(2) * PI) + (iv.sqrt(2) * PI * G['C_D0'] * (tb ** iv.mpf('1.5'))) / ((1 - Del0) * iv.sqrt(ct_min))
    g1_mis = sq_g1 ** 2
    # Exact v in [0, u_0] valley drop on [th_hit, pi] (including exact D = -1 residue term & T(u) > 0.05 far tail!):
    s_max = tb * u_max / PI
    d_lo_h = 1 - Del0
    far_r5_h, _ = far_T5_helpers_D1(G, ct_min, Del0)
    far_W2_h = ((PI ** 2) / 6) * far_r5_h * (3 * ((2 * G['te_hi'] * tb) ** 2) / ((2 * PI * d_lo_h) ** 4) + 8 * G['te_hi'] * tb * G['thc_1ms'] / ((2 * PI * d_lo_h) ** 3))
    d_lam_th_ub = -(1 - (s_max ** 2) - 2 * Del0) + (1 + rel_B1) * (tb / (6 * PI)) + mult_B2 * ((tb ** 2) / (8 * (PI ** 2))) * (2 * u_max + iv.sqrt(PI / (G['c_re5'] * ta))) + far_W2_h
    res_drop_D = c_v0_ratio * (G['C_D0'] * (tb ** 2) / (1 - Del0)) / iv_min(ct_tau * (u_max ** 2), (c_u ** 2) * iv.sqrt(ct_min * (ta ** 2)) / (tb ** iv.mpf('1.5')))
    c_drop = (-d_lam_th_ub) * G['r_re0'] - c_v0_ratio * c_hi * G['te_hi'] * tb / (2 * (PI ** 2)) - res_drop_D
    assert c_drop.a > mp.mpf('0.35'), f"{label} c_drop too small: {c_drop}"
    drop_floor = c_drop * ct_tau * (u_max ** 2)
    sqrt_drop  = c_drop * (c_u ** 2) * PI * iv.sqrt(2) / tb
    # All 4 m = 1 boundary poles at |th| <= f_core*pi*u_0 with exact angle-dependent ratio gamma = t/th_use:
    th_core_max = f_core * PI * u_max
    def ep_in_for_f(f_eff):
        th_eff = f_eff * PI * u_max
        d_st1 = (2 * PI * u_max * (1 - Del0) - th_eff) / (2 * PI * (1 + Del0) + 2 * G['te_hi'] * tb * u_max)
        ep1 = G['c_rein'] * ct_tau * (2 * u_max * d_st1 - d_st1 ** 2) / (1 + Del0)
        rho = ((2 - f_eff) * PI - 2 * PI * Del0) / (2 * PI * (1 + Del0) + 2 * G['te_hi'] * tb * u_max)
        assert rho.a > 0 and rho.b < 1
        ep2 = G['c_rein'] * (2 * rho - rho ** 2) * (c_u ** 2) * iv.sqrt(ct_min * (ta ** 2)) / ((tb ** iv.mpf('1.5')) * (1 + Del0))
        return iv_min(ep1, ep2)
    ep1_out = solve_out_neg_E_D1(G, u_max, u_max, Del0, ct_min / tb, ct_tau, ct_min, PI - th_core_max)
    assert (2 * PI * (1 + Del0) * (iv.sqrt(iv.mpf('1.65')) - 1) + iv.mpf('0.35')).b < ((1 / u_max - f_core) * PI).a
    ep2_out = G['c_re5'] * iv.mpf('0.65') * (c_u ** 2) * iv.sqrt(ct_min * (ta ** 2)) / ((tb ** iv.mpf('1.5')) * iv.mpf('1.02') * (1 + Del0))
    # Explicit 50-digit iv certification of ALL 3 constants (T_c5, 1.02, 0.35) + u*th_c <= 0.080 in ep2_out (#96 fix #4 + #106 fix #5!):
    assert (G['gam'] * (u_max ** 2) * iv.mpf('1.65')).b < iv.mpf('0.05')
    assert (iv.sqrt(iv.mpf('1.65')) * u_max * G['thc_hi']).b <= mp.mpf('0.080')
    assert (1 + G['te_str'] * ep2_out / (4 * (PI ** 2) * (1 + Del0) * (ct_min / tb))).b <= mp.mpf('1.02')
    assert (G['te_hi'] * iv.sqrt(iv.mpf('1.65')) * ep2_out / (PI * (ct_min / tb))
            + G['te_str'] * (ep2_out ** 2) * (tb ** iv.mpf('1.75')) / (4 * (PI ** 2) * c_u * (ct_min ** iv.mpf('1.75')))).b <= mp.mpf('0.35')
    ep_in_0 = ep_in_for_f(iv.mpf(0))
    ep_out = iv_min(ep1_out, ep2_out)
    def E_pk_at_gamma(g_iv):
        g_hi = iv_min_hi(iv.mpf(1), iv.mpf(g_iv.b))
        e_near = ep_in_for_f(g_hi * f_core)
        return -iv.log(1 - iv.exp(-e_near)) - iv.log(1 - iv.exp(-ep_in_0)) - 2 * iv.log(1 - iv.exp(-ep_out)) + E_higher + NB_err_D + Y_big_err_D
    E_pk = E_pk_at_gamma(iv.mpf(1))
    # Single-shot positivity for ALL X >= X_min (exact Z_core(X) = z_core_coeff * X^(1/4), universal c_lo_buf on buffer!):
    X_min = ct_tau / (2 * (PI ** 2))
    z_core_coeff = f_core * (c_u / iv.sqrt(2)) * ((2 * (PI ** 2)) ** iv.mpf('0.25')) / iv.sqrt(tb)
    X_big = iv.mpf(max(iv.mpf('25.0').a, X_min.a))
    Z_c = iv.mpf('1.55')
    Z_core_big = iv_min(th_core_max * iv.sqrt(X_big), z_core_coeff * (X_big ** iv.mpf('0.25')))
    assert Z_c.b <= Z_core_big.a, f"{label} Z_c > Z_core_big={Z_core_big}"
    E_pk_big = E_pk_at_gamma(Z_c / Z_core_big)
    N_z = 60; dz = Z_c / N_z; L_pk_big = iv.mpf(0)
    for i in range(N_z):
        z_iv = iv.mpf([(dz * i).a, (dz * (i + 1)).b])
        phi_z = E_pk_big + iv.sqrt(g1_mis) * z_iv + iv.mpf('0.5') * R3 * (z_iv ** 2) + (K3 / (6 * iv.sqrt(X_big))) * (z_iv ** 3)
        if phi_z.b < (PI / 2).a:
            L_pk_big += iv.exp(-E_pk_big) * iv.mpf((iv.exp(-iv.mpf('0.5') * c_hi * (z_iv ** 2)) * iv.cos(phi_z)).a) * dz
        else:
            L_pk_big -= iv.exp(E_pk_big) * iv.mpf((iv.exp(-iv.mpf('0.5') * c_lo * (z_iv ** 2))).b) * dz
    buf_tail_big = V_full * iv.exp(-iv.mpf('0.5') * c_lo_buf * (Z_core_big ** 2)) / (c_lo_buf * Z_core_big)
    U_pt_big = iv.exp(E_pk) * iv.exp(-iv.mpf('0.5') * c_lo * (Z_c ** 2)) / (c_lo * Z_c) + buf_tail_big
    assert iv.sqrt(X_big).a > (1 / sqrt_drop).b
    U_vl_big = V_full * PI * iv.sqrt(X_big) * iv.exp(-sqrt_drop * iv.sqrt(X_big))
    U_tot_big = iv.exp(E_pk) * iv.sqrt(PI / (2 * c_lo)) + buf_tail_big + U_vl_big
    min_rat = ((L_pk_big - U_pt_big - U_vl_big) / U_tot_big - E_cubic_G3_sub).a
    worst_X = float(X_big.a)
    if X_min.a < X_big.b:
        N_X = 40; N_t = 40
        for k in range(N_X):
            Xa = X_min + (X_big - X_min) * k / N_X; Xb = X_min + (X_big - X_min) * (k + 1) / N_X
            th_use = iv_min(th_core_max, z_core_coeff * (Xa ** iv.mpf('0.25')) / iv.sqrt(Xb))
            dt = th_use / N_t; pk_sum = iv.mpf(0); pk_abs = iv.mpf(0)
            for i in range(N_t):
                t_iv = iv.mpf([(dt * i).a, (dt * (i + 1)).b])
                E_t = E_pk_at_gamma(iv.mpf(t_iv.b) / th_use)
                phi_t = E_t + iv.sqrt(g1_mis * Xb) * t_iv + iv.mpf('0.5') * R3 * Xb * (t_iv ** 2) + (K3 / 6) * Xb * (t_iv ** 3)
                e_hi = iv.exp(E_t) * iv.exp(-iv.mpf('0.5') * c_lo * Xa * (t_iv ** 2))
                pk_abs += iv.mpf(e_hi.b) * dt
                if phi_t.b < (PI / 2).a:
                    pk_sum += iv.mpf((iv.exp(-E_t) * iv.exp(-iv.mpf('0.5') * c_hi * Xb * (t_iv ** 2)) * iv.cos(phi_t)).a) * dt
                else:
                    pk_sum -= iv.mpf(e_hi.b) * dt
            # Exact conservative buffer [th_use, PI] using universal c_lo_buf valid on ALL |th| < th_hit:
            mid_abs = iv.mpf(0); N_m = 20; dm = (PI - th_use) / N_m
            for i in range(N_m):
                t_iv = iv.mpf([(th_use + dm * i).a, (th_use + dm * (i + 1)).b])
                mid_abs += V_full * iv.mpf((iv.exp(-iv.mpf('0.5') * c_lo_buf * Xa * (t_iv ** 2))).b) * dm
            d_box = iv.mpf(max(drop_floor.a, (sqrt_drop * iv.sqrt(Xa)).a))
            vl_abs = V_full * (PI - th_use) * iv.exp(-d_box)
            rat_b = (pk_sum - mid_abs - vl_abs) / (pk_abs + mid_abs + vl_abs) - E_cubic_G3_sub
            if rat_b.a < min_rat:
                min_rat = rat_b.a; worst_X = float(Xa.a)
    print(f"  [{label:30s}] ct>={float(ct_min.a):5.1f} | c=[{float(c_lo.a):.3f},{float(c_hi.b):.3f}] (buf>={float(c_lo_buf.a):.3f}) u*<={float(u_star_ub.b):.3f} | E_pk<={float(E_pk.b):.4f} V<={float(V_full.b):.2f} | rat>={float(min_rat):+.4f} (worst X={worst_X:.1f})", flush=True)
    assert float(min_rat) > 0.05, f"FAILED {label}: min_rat = {float(min_rat)}"
    return float(min_rat)

def find_ivt_bracket(ta, tb, ct_lo, u_guess, ct_hi=None, use_cauchy_uL=False):
    uL_iv = None; uR_iv = None
    if use_cauchy_uL and tb.b <= mp.mpf('2.85'):
        uL_iv, uR_iv = harmonic_u_star_bracket_iv(ta, tb, ct_lo, n_sub_t=4)
        # Also tighten via eval_g1_fin_iv when possible (e.g., tau >= 2.20)!
        if ta.a >= mp.mpf('2.19'):
            nb_use = 90
            uL_try = min(u_guess, float(uR_iv.b))
            gL = eval_g1_fin_iv(ta, tb, iv.mpf(str(round(uL_try, 3))), ct_lo, n_bins=nb_use, ct_hi=ct_hi)
            for _ in range(14):
                if gL.b < 0: break
                uL_try -= 0.015
                if uL_try <= float(uL_iv.a) + 0.005: break
                gL = eval_g1_fin_iv(ta, tb, iv.mpf(str(round(uL_try, 3))), ct_lo, n_bins=nb_use, ct_hi=ct_hi)
            if gL.b < 0 and iv.mpf(str(round(uL_try, 3))).a > uL_iv.a:
                uL_iv = iv.mpf(str(round(uL_try, 3)))
        nb_R = 90
        uR_try = min(float(uR_iv.b), max(u_guess, float(uL_iv.a) + 0.03))
        gR = eval_g1_fin_iv(ta, tb, iv.mpf(str(round(uR_try, 3))), ct_lo, n_bins=nb_R, ct_hi=ct_hi)
        for _ in range(14):
            if gR.a > 0: break
            uR_try += 0.015
            if uR_try >= float(uR_iv.b): break
            gR = eval_g1_fin_iv(ta, tb, iv.mpf(str(round(uR_try, 3))), ct_lo, n_bins=nb_R, ct_hi=ct_hi)
        if gR.a > 0 and iv.mpf(str(round(uR_try, 3))).b < uR_iv.b:
            uR_iv = iv.mpf(str(round(uR_try, 3)))
        assert uL_iv.a < uR_iv.b
        return iv.mpf([uL_iv.a, uR_iv.b]), float(uL_iv.a), float(uR_iv.b)
    uL = u_guess; uR = u_guess
    nb_use = 90 if tb.b <= mp.mpf('3.50') else (64 if tb.b <= mp.mpf('6.00') else 40)
    gL = eval_g1_fin_iv(ta, tb, iv.mpf(str(round(uL, 3))), ct_lo, n_bins=nb_use, ct_hi=ct_hi)
    for _ in range(20):
        if gL.b < 0: break
        uL -= 0.015
        if uL < 0.05: break
        gL = eval_g1_fin_iv(ta, tb, iv.mpf(str(round(uL, 3))), ct_lo, n_bins=nb_use, ct_hi=ct_hi)
    if gL.b >= 0:
        gL = eval_g1_fin_iv(ta, tb, iv.mpf(str(round(uL, 3))), ct_lo, n_bins=120, ct_hi=ct_hi)
    nb_R = 90 if tb.b <= mp.mpf('3.50') else (64 if tb.b <= mp.mpf('6.00') else 40)
    gR = eval_g1_fin_iv(ta, tb, iv.mpf(str(round(uR, 3))), ct_lo, n_bins=nb_R, ct_hi=ct_hi)
    for _ in range(20):
        if gR.a > 0: break
        uR += 0.015
        if uR > 0.52: uR = 0.52; break
        gR = eval_g1_fin_iv(ta, tb, iv.mpf(str(round(uR, 3))), ct_lo, n_bins=nb_R, ct_hi=ct_hi)
    if gR.a <= 0:
        gR = eval_g1_fin_iv(ta, tb, iv.mpf(str(round(uR, 3))), ct_lo, n_bins=120, ct_hi=ct_hi)
    assert gL.b < 0 and gR.a > 0, f"IVT failed on [{float(ta.a)},{float(tb.b)}] ct>={float(ct_lo.a)}: uL={uL:.3f} gL=[{float(gL.a):+.4f},{float(gL.b):+.4f}], uR={uR:.3f} gR=[{float(gR.a):+.4f},{float(gR.b):+.4f}]"
    return iv.mpf([iv.mpf(str(round(uL, 3))).a, iv.mpf(str(round(uR, 3))).b]), round(uL, 3), round(uR, 3)

def certify_G3D_exact_subdiv(ta, tb, ct_lo, ct_hi, u_g, th_c_cap, label_prefix, max_du=0.018, use_cauchy_uL=False):
    """Certify exact finite-n saddle u0* in [uL, uR] by subdividing [uL, uR] into thin sub-intervals (width <= max_du)."""
    u0_all, uL, uR = find_ivt_bracket(ta, tb, ct_lo, u_g, ct_hi=ct_hi, use_cauchy_uL=use_cauchy_uL)
    n_sub = max(1, int((uR - uL) / max_du + 0.999999))
    uL_iv = iv.mpf(u0_all.a); uR_iv = iv.mpf(u0_all.b)
    rats = []
    for k in range(n_sub):
        ua = uL_iv + (uR_iv - uL_iv) * k / n_sub
        ub = uL_iv + (uR_iv - uL_iv) * (k + 1) / n_sub
        u_sub = iv.mpf([ua.a, ub.b])
        uaf = float(ua.a); ubf = float(ub.b)
        th_hit_sub = fpi * (1 - 2 * max(abs(0.5 - uaf), abs(0.5 - ubf)))
        if tb.b <= mp.mpf('2.85'):
            th_c = min(th_c_cap, round(min(th_hit_sub - 0.18, 0.78 * th_hit_sub), 2))
        else:
            th_c = min(th_c_cap, round(min(th_hit_sub - 0.35, 0.68 * th_hit_sub), 2))
        th_c = max(th_c, 0.25)
        lbl = f"{label_prefix} [{uaf:.3f},{ubf:.3f}]" if n_sub > 1 else f"{label_prefix} u*=[{uL:.2f},{uR:.2f}]"
        rats.append(certify_G3D_box(ta, tb, u_sub, ct_lo, ct_hi, iv.mpf(0), th_c, lbl, exact_saddle=True))
    return min(rats)

# ==============================================================================
# PART 2: PARALLEL WORKER FUNCTIONS + FULL 50-DIGIT iv CERTIFICATION OF CASE G3-D!
# ==============================================================================
def _safe_run(fn, *a, **kw):
    try:
        return fn(*a, **kw), None
    except Exception as e:
        print(f"  [!!FAIL!!] {e}", flush=True)
        return -99.0, str(e)

def _run_D3_slice(args):
    t_a, t_b, u_g = args
    ta = iv.mpf(str(t_a)); tb = iv.mpf(str(t_b))
    ct_lo = iv.mpf(max(CT_MIN.a, ((ta ** 2) / R3CK2_MAX).a))
    rats = []; fails = []
    def add(res):
        r, err = res
        rats.append(r)
        if err is not None: fails.append(err)
    if t_a < 4.00:
        ct_1 = iv.mpf(14); ct_2 = iv.mpf(24); ct_3 = iv.mpf(48)
        ct_sw = iv.mpf(90) if t_a < 3.18 else (iv.mpf(75) if t_a < 3.68 else iv.mpf(60))
        u_1 = 0.42 if t_a < 3.15 else (0.395 if t_a < 3.45 else 0.400)
        u_2 = round(0.45 * u_1 + 0.55 * u_g, 3)
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf(str(u_1)), ct_lo, ct_1, iv.mpf(0), 1.45, f"D3A1 [{t_a:4.2f},{t_b:4.2f}] u={u_1:.2f}"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf(str(u_2)), ct_1, ct_2, iv.mpf(0), 1.38, f"D3A2 [{t_a:4.2f},{t_b:4.2f}] u={u_2:.2f}"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf(str(u_g)), ct_2, ct_3, iv.mpf(0), 1.20, f"D3A3 [{t_a:4.2f},{t_b:4.2f}] u={u_g:.2f}"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf(str(u_g)), ct_3, ct_sw, iv.mpf(0), 1.15, f"D3A4 [{t_a:4.2f},{t_b:4.2f}] u={u_g:.2f}"))
        add(_safe_run(certify_G3D_exact_subdiv, ta, tb, ct_sw, None, u_g, 1.10, f"D3B  [{t_a:4.2f},{t_b:4.2f}]", 0.010 if t_a < 2.98 else 0.015))
    elif t_a < 9.50:
        ct_lo_f = float(ct_lo.a)
        ct_sw_f = 60.0 if t_a < 4.82 else 55.0
        ct_r2 = int(t_a / 0.45) - 0.5
        raw_cuts = [12.5, 16.0, 24.0, 36.0] if t_a < 4.82 else [14.0, 22.0, 34.0, 45.0]
        if t_a >= 4.95 and ct_lo_f + 1.0 < ct_r2 < ct_sw_f - 1.5:
            raw_cuts.append(ct_r2)
        cuts = [ct_lo_f]
        for x in sorted(raw_cuts):
            if ct_lo_f + 1.0 < x < ct_sw_f - 1.5 and x - cuts[-1] >= 1.5:
                cuts.append(x)
        cuts.append(ct_sw_f)
        th_c_base = 1.45 if t_a < 4.82 else (1.55 if t_a < 6.05 else (1.70 if t_a < 7.3 else 1.85))
        for idx in range(len(cuts) - 1):
            ca = ct_lo if idx == 0 else iv.mpf(str(cuts[idx]))
            cb = iv.mpf(str(cuts[idx + 1]))
            if t_a >= 5.68 and cuts[idx] >= 13.5:
                u_use = round(u_g - 0.004, 3)
                th_use = 1.30
            else:
                u_use = u_g
                th_use = th_c_base
            add(_safe_run(certify_G3D_box, ta, tb, iv.mpf(str(u_use)), ca, cb, iv.mpf(0), th_use, f"D3A{idx+1} [{t_a:4.2f},{t_b:4.2f}] u={u_use:.3f}"))
        add(_safe_run(certify_G3D_exact_subdiv, ta, tb, iv.mpf(str(ct_sw_f)), None, u_g, min(th_c_base, 1.45), f"D3B  [{t_a:4.2f},{t_b:4.2f}]", 0.015))
    else:
        ct_lo_f = float(ct_lo.a)
        ct_sw_f = max(ct_lo_f, 55.0 if t_a < 20.5 else 45.0)
        th_c_base = 1.95 if t_a < 16.5 else 2.15
        ct_r3 = int(t_a / 0.90) - 0.5
        ct_r2 = int(t_a / 0.45) - 0.5
        raw_cuts = [ct_r3, ct_r2, 28.0, 42.0] if t_a < 20.5 else [ct_r3, ct_r2, 28.0]
        cuts = [ct_lo_f]
        for x in sorted(raw_cuts):
            if ct_lo_f + 0.5 < x < ct_sw_f - 0.5:
                if x - cuts[-1] >= 1.5:
                    cuts.append(x)
        if ct_sw_f > cuts[-1] + 0.5:
            cuts.append(ct_sw_f)
        for idx in range(len(cuts) - 1):
            ca = ct_lo if idx == 0 else iv.mpf(str(cuts[idx]))
            cb = iv.mpf(str(cuts[idx + 1]))
            if t_a < 13.5 and cuts[idx] >= 18.0:
                u_use = round(u_g - 0.003, 3)
                th_use = 1.30
            else:
                u_use = u_g
                th_use = th_c_base
            add(_safe_run(certify_G3D_box, ta, tb, iv.mpf(str(u_use)), ca, cb, iv.mpf(0), th_use, f"D3A{idx+1} [{t_a:4.1f},{t_b:4.1f}] u={u_use:.3f}"))
        ct_B = ct_lo if ct_lo_f >= ct_sw_f - 0.5 else iv.mpf(str(ct_sw_f))
        add(_safe_run(certify_G3D_exact_subdiv, ta, tb, ct_B, None, u_g, 1.55 if t_a < 20.5 else 1.65, f"D3B  [{t_a:4.1f},{t_b:4.1f}]", 0.015))
    return min(rats), fails

def _run_D2_slice(args):
    t_a, t_b, u_g = args
    ta = iv.mpf(str(t_a)); tb = iv.mpf(str(t_b))
    ct_min_box = iv.mpf(max(CT_MIN.a, (CTAU2_MIN / (tb ** 2)).a))
    rats = []; fails = []
    def add(res):
        r, err = res
        rats.append(r)
        if err is not None: fails.append(err)
    if t_a < 1.65:
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.43'), ct_min_box, iv.mpf(24), iv.mpf(0), 1.40, f"D2-T1 [{t_a:.2f},{t_b:.2f}] u=0.43"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.35'), iv.mpf(24), iv.mpf(52), iv.mpf(0), 1.22, f"D2-T2 [{t_a:.2f},{t_b:.2f}] u=0.35"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.28'), iv.mpf(52), iv.mpf(88), iv.mpf(0), 1.00, f"D2-T3 [{t_a:.2f},{t_b:.2f}] u=0.28"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.23'), iv.mpf(88), iv.mpf(135), iv.mpf(0), 0.80, f"D2-T4 [{t_a:.2f},{t_b:.2f}] u=0.23"))
        add(_safe_run(certify_G3D_harmonic_tail, ta, tb, iv.mpf(135), iv.mpf('0.26'), iv.mpf('1.08'), f"D2-H  [{t_a:.2f},{t_b:.2f}]"))
    elif t_a < 1.82:
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.45'), ct_min_box, iv.mpf(20), iv.mpf(0), 1.45, f"D2-T1 [{t_a:.2f},{t_b:.2f}] u=0.45"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.36'), iv.mpf(20), iv.mpf(44), iv.mpf(0), 1.25, f"D2-T2 [{t_a:.2f},{t_b:.2f}] u=0.36"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.28'), iv.mpf(44), iv.mpf(95), iv.mpf(0), 1.00, f"D2-T3 [{t_a:.2f},{t_b:.2f}] u=0.28"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.22'), iv.mpf(95), iv.mpf(185), iv.mpf(0), 0.78, f"D2-T4 [{t_a:.2f},{t_b:.2f}] u=0.22"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.185'), iv.mpf(185), iv.mpf(360), iv.mpf(0), 0.62, f"D2-T5 [{t_a:.2f},{t_b:.2f}] u=0.185"))
        add(_safe_run(certify_G3D_exact_subdiv, ta, tb, iv.mpf(360), None, u_g, 0.90, f"D2-B  [{t_a:.2f},{t_b:.2f}]", 0.025, True))
    elif t_a < 2.00:
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.46'), ct_min_box, iv.mpf(17), iv.mpf(0), 1.45, f"D2-T1 [{t_a:.2f},{t_b:.2f}] u=0.46"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.38'), iv.mpf(17), iv.mpf(36), iv.mpf(0), 1.30, f"D2-T2 [{t_a:.2f},{t_b:.2f}] u=0.38"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.30'), iv.mpf(36), iv.mpf(78), iv.mpf(0), 1.06, f"D2-T3 [{t_a:.2f},{t_b:.2f}] u=0.30"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.24'), iv.mpf(78), iv.mpf(155), iv.mpf(0), 0.84, f"D2-T4 [{t_a:.2f},{t_b:.2f}] u=0.24"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.205'), iv.mpf(155), iv.mpf(280), iv.mpf(0), 0.70, f"D2-T5 [{t_a:.2f},{t_b:.2f}] u=0.205"))
        add(_safe_run(certify_G3D_exact_subdiv, ta, tb, iv.mpf(280), None, u_g, 0.95, f"D2-B  [{t_a:.2f},{t_b:.2f}]", 0.012, True))
    elif t_a < 2.20:
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.46'), ct_min_box, iv.mpf(17), iv.mpf(0), 1.45, f"D2-T1 [{t_a:.2f},{t_b:.2f}] u=0.46"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.38'), iv.mpf(17), iv.mpf(36), iv.mpf(0), 1.30, f"D2-T2 [{t_a:.2f},{t_b:.2f}] u=0.38"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.30'), iv.mpf(36), iv.mpf(78), iv.mpf(0), 1.06, f"D2-T3 [{t_a:.2f},{t_b:.2f}] u=0.30"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.25'), iv.mpf(78), iv.mpf(145), iv.mpf(0), 0.88, f"D2-T4 [{t_a:.2f},{t_b:.2f}] u=0.25"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.225'), iv.mpf(145), iv.mpf(220), iv.mpf(0), 0.78, f"D2-T5 [{t_a:.2f},{t_b:.2f}] u=0.225"))
        add(_safe_run(certify_G3D_exact_subdiv, ta, tb, iv.mpf(220), None, u_g, 1.00, f"D2-B  [{t_a:.2f},{t_b:.2f}]", 0.010, True))
    elif t_a < 2.40:
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.45'), ct_min_box, iv.mpf(17), iv.mpf(0), 1.45, f"D2-T1 [{t_a:.2f},{t_b:.2f}] u=0.45"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.37'), iv.mpf(17), iv.mpf(38), iv.mpf(0), 1.28, f"D2-T2 [{t_a:.2f},{t_b:.2f}] u=0.37"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.30'), iv.mpf(38), iv.mpf(85), iv.mpf(0), 1.05, f"D2-T3 [{t_a:.2f},{t_b:.2f}] u=0.30"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.255'), iv.mpf(85), iv.mpf(150), iv.mpf(0), 0.88, f"D2-T4 [{t_a:.2f},{t_b:.2f}] u=0.255"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.245'), iv.mpf(150), iv.mpf(220), iv.mpf(0), 0.84, f"D2-T5 [{t_a:.2f},{t_b:.2f}] u=0.245"))
        add(_safe_run(certify_G3D_exact_subdiv, ta, tb, iv.mpf(220), None, u_g, 1.05, f"D2-B  [{t_a:.2f},{t_b:.2f}]", 0.010, True))
    elif t_a < 2.60:
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.44'), ct_min_box, iv.mpf(17), iv.mpf(0), 1.45, f"D2-T1 [{t_a:.2f},{t_b:.2f}] u=0.44"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.37'), iv.mpf(17), iv.mpf(38), iv.mpf(0), 1.28, f"D2-T2 [{t_a:.2f},{t_b:.2f}] u=0.37"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.30'), iv.mpf(38), iv.mpf(85), iv.mpf(0), 1.05, f"D2-T3 [{t_a:.2f},{t_b:.2f}] u=0.30"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.275'), iv.mpf(85), iv.mpf(180), iv.mpf(0), 0.95, f"D2-T4 [{t_a:.2f},{t_b:.2f}] u=0.275"))
        add(_safe_run(certify_G3D_exact_subdiv, ta, tb, iv.mpf(180), None, u_g, 1.10, f"D2-B  [{t_a:.2f},{t_b:.2f}]", 0.008, True))
    else:
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.43'), ct_min_box, iv.mpf(17), iv.mpf(0), 1.45, f"D2-T1 [{t_a:.2f},{t_b:.2f}] u=0.43"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.36'), iv.mpf(17), iv.mpf(38), iv.mpf(0), 1.28, f"D2-T2 [{t_a:.2f},{t_b:.2f}] u=0.36"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.31'), iv.mpf(38), iv.mpf(80), iv.mpf(0), 1.08, f"D2-T3 [{t_a:.2f},{t_b:.2f}] u=0.31"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.295'), iv.mpf(80), iv.mpf(135), iv.mpf(0), 1.02, f"D2-T4 [{t_a:.2f},{t_b:.2f}] u=0.295"))
        add(_safe_run(certify_G3D_box, ta, tb, iv.mpf('0.290'), iv.mpf(135), iv.mpf(185), iv.mpf(0), 1.00, f"D2-T5 [{t_a:.2f},{t_b:.2f}] u=0.290"))
        add(_safe_run(certify_G3D_exact_subdiv, ta, tb, iv.mpf(185), None, u_g, 1.10, f"D2-B  [{t_a:.2f},{t_b:.2f}]", 0.008, True))
    return min(rats), fails

def _run_D1_slice(args):
    t_a, t_b, u_g = args
    ta = iv.mpf(str(t_a)); tb = iv.mpf(str(t_b))
    ct_min_box = iv.mpf(max(CT_MIN.a, (CTAU2_MIN / (tb ** 2)).a))
    if t_a < 1.25:
        ct1 = iv.mpf(48); ct2 = iv.mpf(80); ct3 = iv.mpf(115)
    elif t_a < 1.37:
        ct1 = iv.mpf(42); ct2 = iv.mpf(72); ct3 = iv.mpf(105)
    else:
        ct1 = iv.mpf(42); ct2 = iv.mpf(78); ct3 = iv.mpf(120)
    rats = []; fails = []
    for res in [
        _safe_run(certify_G3D_box, ta, tb, iv.mpf('0.42'), ct_min_box, ct1, iv.mpf(0), 1.25, f"D1-T1 [{t_a:.2f},{t_b:.2f}] u=0.42"),
        _safe_run(certify_G3D_box, ta, tb, iv.mpf('0.33'), ct1, ct2, iv.mpf(0), 1.02, f"D1-T2 [{t_a:.2f},{t_b:.2f}] u=0.33"),
        _safe_run(certify_G3D_box, ta, tb, iv.mpf('0.26'), ct2, ct3, iv.mpf(0), 0.82, f"D1-T3 [{t_a:.2f},{t_b:.2f}] u=0.26"),
        _safe_run(certify_G3D_harmonic_tail, ta, tb, ct3, iv.mpf('0.24'), iv.mpf('1.08'), f"D1-H  [{t_a:.2f},{t_b:.2f}]"),
    ]:
        rats.append(res[0])
        if res[1] is not None: fails.append(res[1])
    return min(rats), fails

def _run_any_slice(item):
    band, args = item
    if band == 'D3': r, fs = _run_D3_slice(args)
    elif band == 'D2': r, fs = _run_D2_slice(args)
    else: r, fs = _run_D1_slice(args)
    return band, r, fs

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--part", choices=["all", "D1", "D2", "D3a", "D3b", "D3a1", "D3a2", "D3b1", "D3b2"], default="all")
  parser.add_argument("--workers", type=int, default=3)
  args = parser.parse_args()
  run_part0_selftests()
  n_workers = max(1, min(args.workers, os.cpu_count() or 2))
  print("\n==========================================================================")
  print(f"PART 2: FULL 50-DIGIT mpmath.iv CERTIFICATION OF CASE G3-D (part={args.part}, workers={n_workers})")
  print("==========================================================================")

  taus_D3 = [2.80, 2.98, 3.18, 3.40, 3.54, 3.68, 3.84, 4.00,
             4.18, 4.38, 4.60, 4.82, 5.08, 5.36, 5.68, 6.05, 6.50, 7.10, 7.80, 8.60, 9.50,
             11.2, 13.5, 16.5, 20.5, 25.5, 32.0, 40.0]
  u_par_D3 = [0.320, 0.335, 0.351, 0.363, 0.374, 0.384, 0.394,
              0.402, 0.411, 0.420, 0.429, 0.437, 0.446, 0.454, 0.461, 0.469, 0.476, 0.483, 0.489, 0.494,
              0.499, 0.502, 0.504, 0.504, 0.503, 0.502, 0.501]
  taus_D2 = [1.50, 1.65, 1.82, 2.00, 2.20, 2.40, 2.60, 2.80]
  u_par_D2 = [0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.31]
  taus_D1 = [1.05, 1.25, 1.37, 1.50]
  u_par_D1 = [0.12, 0.15, 0.18]

  jobs_D1 = [('D1', (taus_D1[k], taus_D1[k + 1], u_par_D1[k])) for k in range(len(taus_D1) - 1)]
  jobs_D2 = [('D2', (taus_D2[k], taus_D2[k + 1], u_par_D2[k])) for k in range(len(taus_D2) - 1)]
  jobs_D3 = [('D3', (taus_D3[k], taus_D3[k + 1], u_par_D3[k])) for k in range(len(taus_D3) - 1)]
  if args.part == "D1": jobs = jobs_D1
  elif args.part == "D2": jobs = jobs_D2
  elif args.part == "D3a": jobs = jobs_D3[:14]
  elif args.part == "D3b": jobs = jobs_D3[14:]
  elif args.part == "D3a1": jobs = jobs_D3[:7]
  elif args.part == "D3a2": jobs = jobs_D3[7:14]
  elif args.part == "D3b1": jobs = jobs_D3[14:21]
  elif args.part == "D3b2": jobs = jobs_D3[21:]
  else: jobs = jobs_D2 + jobs_D1 + jobs_D3

  all_fails = []
  band_rats = {'D1': [], 'D2': [], 'D3': []}
  if n_workers == 1:
    for item in jobs:
      band, r, fs = _run_any_slice(item)
      band_rats[band].append(r)
      all_fails.extend(fs)
  else:
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
      for band, r, fs in pool.map(_run_any_slice, jobs):
        band_rats[band].append(r)
        all_fails.extend(fs)

  if all_fails:
    print(f"\n[FAIL] {len(all_fails)} box(es) failed:")
    for msg in all_fails:
      print("  -", msg)
    sys.exit(1)

  summary_parts = []
  for bname in ['D1', 'D2', 'D3']:
    if band_rats[bname]:
      m = min(band_rats[bname])
      assert m > 0.05, f"Band {bname} min_rat={m} <= 0.05"
      summary_parts.append(f"{bname} >= {m:+.4f}")
  print("\n==========================================================================")
  print(f"CASE G3-D (part={args.part}) 100% 50-DIGIT iv CERTIFIED!!")
  print(f"  Global min rat_G3: {', '.join(summary_parts)}")
  print(f"  Elapsed: {time.time() - t_start:.2f}s")
  print("==========================================================================")
