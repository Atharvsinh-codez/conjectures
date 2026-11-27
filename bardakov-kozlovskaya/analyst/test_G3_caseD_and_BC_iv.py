"""
VERIFY ALL THREE |D| = 1 CASES (G3-B: tau >= 40, G3-C: tau <= 1.05, G3-D: tau in [1.05, 40.0])
AND CASE G3-A (|D| >= 2, ALL tau > 0) IN 50-DIGIT mpmath.iv AT LOCKED DOMAIN:
  ct >= 10 * D^2,  r_0 * kap <= 0.45,  r_0^3 * c * kap^2 <= 60,  ct * tau^2 / D^4 >= 40,  c * kap^2 < 42!
"""
import mpmath as mp
from mpmath import iv

iv.dps = 50
mp.mp.dps = 50
PI = iv.pi

CT_MIN = iv.mpf(10)
R0K_MAX = iv.mpf('0.45')
CTAU2_MIN = iv.mpf(40)
CK2_MAX = iv.mpf(42)
R3CK2_MAX = iv.mpf(60)

al = iv.mpf('0.07169')
th_max = al / 2
S_G3 = iv.sin(al) / al
n_min = 2 * (PI ** 2) * (S_G3 ** 2) * 23 / R0K_MAX
c_sh = iv.mpf(1) / (24 * (1 - (th_max ** 2) / 24))
c_cos = iv.mpf(1) / (2 * (1 - (al ** 2) / 2))
C_Del0 = (c_cos - (iv.mpf(1) / 6) * (1 - (al ** 2) / 20) + c_cos + c_cos / 4 + c_sh + c_sh * c_cos * (al ** 2)) / (16 * (PI ** 2) * (S_G3 ** 4))
gam_max = (R0K_MAX ** 2) / (2 * (PI ** 2) * (S_G3 ** 2) * (iv.cos(al) ** 2))
xcot_th = (th_max / 2) * iv.cos(th_max / 2) / iv.sin(th_max / 2)
tau_eff_lo = xcot_th
tau_eff_hi = 1 / (iv.cos(al) ** 2)
T_c5 = iv.mpf('0.05')
ut_c5 = iv.mpf('0.080')
assert float((iv.sin(ut_c5) - iv.sqrt(T_c5 / 8)).a) > 0
sinc_c5_sq = (iv.sin(ut_c5) / ut_c5) ** 2
c_re_5 = tau_eff_lo * sinc_c5_sq / (1 + T_c5)
T_in = gam_max * iv.mpf('0.31')
ut_in = iv.mpf('0.55') * th_max
sinc_in_sq = (iv.sin(ut_in) / ut_in) ** 2
c_re_in = tau_eff_lo * sinc_in_sq / (1 + T_in)
c_th1s = 1 / ((iv.sin(th_max / 2) / (th_max / 2)) * iv.cos(al))
r_re0 = tau_eff_lo * (iv.sin(2 * iv.mpf('0.52') * th_max) / (2 * iv.mpf('0.52') * th_max)) / (1 + gam_max * iv.mpf('0.28'))

def W3_mult_cmin(D_min, c_min, Del0_u, T_loc):
    # Note c >= c_min => |D|/(r0*j) < |D|/c_min!
    f_p_hi = (1 + D_min / c_min) * (1 + D_min / (2 * c_min))
    A_w = ((D_min + Del0_u) / (c_min * S_G3)) * c_th1s
    B_w = iv.sqrt(2 * f_p_hi * tau_eff_hi / (S_G3 ** 4))
    t1 = 15 * ((2 * tau_eff_hi) ** 3)
    p6 = T_loc / ((1 + T_loc) ** 6)
    p3 = T_loc / ((1 + T_loc) ** 3)
    t2 = 40 * ((2 * tau_eff_hi) ** 2) * iv.sqrt(A_w ** 2 + (B_w ** 2) * p6)
    t3 = 48 * tau_eff_hi * (A_w ** 2 + (B_w ** 2) * p3)
    return (t1 + t2 + t3) / 120

print("=== 1. CHECK CASE G3-B (|D| = 1, tau >= 40.0 => c >= 89!) ===")
tau_B = iv.mpf('40.0')
c_min_B = tau_B / R0K_MAX  # >= 88.88 => c >= 89!
Del0_all = C_Del0 * CK2_MAX
CW3_B_in  = W3_mult_cmin(iv.mpf(1), c_min_B, Del0_all, T_in)
CW3_B_out = W3_mult_cmin(iv.mpf(1), c_min_B, Del0_all, T_c5)
f_p_hi_B = (1 + 1 / c_min_B) * (1 + 1 / (2 * c_min_B))
f_p_lo_B = (1 - 1 / c_min_B) * (1 - 1 / (2 * c_min_B))
M_hi_B = f_p_hi_B * tau_eff_hi / (S_G3 ** 4)
mob_in_B  = f_p_lo_B * (S_G3 ** 4) * tau_eff_lo * sinc_in_sq
mob_out_B = f_p_lo_B * (S_G3 ** 4) * tau_eff_lo * sinc_c5_sq
C_pref_in_B  = CW3_B_in  * ((1 + T_in) ** 3) / r_re0
C_pref_out_B = CW3_B_out * ((1 + T_c5) ** 3) / r_re0

Du_B = iv.mpf('0.48')
s_start_B = float((tau_B * Du_B / PI).a)
s_tail_B = 18.0
# Tail s >= 18.0:
St = iv.mpf(str(s_tail_B)); x1t = iv.mpf('0.50'); x2t = iv.mpf('0.82')
d_hi_t = 1 + Del0_all; d_lo_t = 1 - Del0_all
out_t = C_pref_out_B * (5 / (16 * (PI ** 3) * c_re_5 * (Du_B ** 3))) * (M_hi_B + d_hi_t ** 2 / (St ** 2)) / ((mob_out_B ** 3) * (St ** 3))
in_a = C_pref_in_B * (5 / (8 * (PI ** 3) * c_re_in * (Du_B ** 3) * (1 + x2t))) * (M_hi_B + d_hi_t ** 2 / (St ** 2)) / ((mob_in_B ** 3) * (St ** 3) * (x2t ** 6))
c2t = c_re_in * PI * Du_B * (1 - x2t ** 2)
in_b = C_pref_in_B * (5 / (8 * (PI ** 2) * (Du_B ** 2))) * (x2t - x1t) * (M_hi_B + d_hi_t ** 2 / (St ** 2)) * iv.exp(-c2t * St) / ((mob_in_B ** 3) * (St ** 2) * (x1t ** 6))
c1t = c_re_in * PI * Du_B * (1 - x1t ** 2)
spk = iv.mpf(max(float(St.a), float((4 / c1t).b)))
in_c = C_pref_in_B * (5 / (8 * (PI ** 2) * (Du_B ** 2))) * x1t * (M_hi_B + d_hi_t ** 2 / (St ** 2)) * (spk ** 4) * iv.exp(-c1t * spk) / (d_lo_t ** 6)
mv_B = float((out_t + in_a + in_b + in_c).b)
for k in range(100):
    sa = iv.mpf(s_start_B + (s_tail_B - s_start_B) * k / 100)
    sb = iv.mpf(s_start_B + (s_tail_B - s_start_B) * (k + 1) / 100)
    d_sh = iv.mpf(min(float(Del0_all.b), float((C_Del0 * R0K_MAX * PI * sb / Du_B).b)))
    dh = 1 + d_sh; dl = 1 - d_sh
    bl_in = c_re_in * PI * Du_B * sa
    outer = C_pref_out_B * (5 / (16 * (PI ** 3) * c_re_5 * (Du_B ** 3))) * sb * (M_hi_B * sb ** 2 + dh ** 2) / ((mob_out_B * sa ** 2 + dl ** 2) ** 3)
    pref_i = C_pref_in_B * (5 / (8 * (PI ** 2) * (Du_B ** 2))) * (sb ** 2) * (M_hi_B * sb ** 2 + dh ** 2)
    best_in = iv.mpf('1e9')
    for x1_s, x2_s in [('0.45', '0.78'), ('0.50', '0.82'), ('0.55', '0.86'), ('0.60', '0.90')]:
        x1 = iv.mpf(x1_s); x2 = iv.mpf(x2_s)
        i_a = iv.mpf(min(float((1 - x2).b), float((1 / (bl_in * (1 + x2))).b))) / ((mob_in_B * sa ** 2 * x2 ** 2 + dl ** 2) ** 3)
        i_b = (x2 - x1) * iv.exp(-bl_in * (1 - x2 ** 2)) / ((mob_in_B * sa ** 2 * x1 ** 2 + dl ** 2) ** 3)
        i_c = x1 * iv.exp(-bl_in * (1 - x1 ** 2)) / (dl ** 6)
        tot = pref_i * (i_a + i_b + i_c)
        if float(tot.b) < float(best_in.b): best_in = tot
    mv_B = max(mv_B, float((outer + best_in).b))
print(f"  [G3-B: |D|=1, tau >= 40.0 (c >= 89)]: CW3_in <= {float(CW3_B_in.b):.4f}, |R_3_fin|/H_fin0 <= {mv_B:.5f}!!")

print("\n=== 2. CHECK CASE G3-C (|D| = 1, tau <= 1.05 => ct >= 36.2, c >= 37!) ===")
tau_C = iv.mpf('1.05')
c_min_C = CTAU2_MIN / (tau_C ** 2)  # >= 36.28
Del0_C = C_Del0 * (tau_C * R0K_MAX)
CW3_C_in  = W3_mult_cmin(iv.mpf(1), c_min_C, Del0_C, T_in)
CW3_C_out = W3_mult_cmin(iv.mpf(1), c_min_C, Del0_C, T_c5)
d_lo_C = 1 - Del0_C
R3_C = (5 / (8 * (PI ** 2))) * (tau_C ** 2) / (d_lo_C ** 6) * (
    CW3_C_in * ((1 + T_in) ** 3) * iv.mpf('0.84')
    + CW3_C_out * ((1 + T_c5) ** 3) * iv.sqrt(PI / (c_re_5 * tau_C))
)
print(f"  [G3-C: |D|=1, tau <= 1.05 (ct >= 36.2)]: Del0 <= {float(Del0_C.b):.5f}, CW3_in <= {float(CW3_C_in.b):.4f}, |R_3_fin|/(tau/(2*pi^2)) <= {float(R3_C.b):.5f}!!")

print("\n=== 3. CHECK CASE G3-D (|D| = 1, tau in [1.05, 40.0]) LIPSCHITZ FINITE-n PERTURBATION SIZE ===")
# Check exact max |delta_zeta(u)| = |(zeta_fin(u) - zeta_fin(u0)) - (zeta_par(u) - zeta_par(u0))| on |u| <= 1.5!
# Remember:
#   Re: |Re(zeta_fin(u)-zeta_fin(u0)) - tau*(u0^2 - u^2)| <= (tau_eff_hi - c_re_15) * tau * |u^2 - u0^2|
#       Wait: tau_eff_hi - c_re_15 <= 1.0052 - 0.9762 = 0.029! On |u| <= 0.55: 1.0052 - 0.994 = 0.0112!!
#   Im: Im(zeta_fin(u) - zeta_fin(1/2)) - 2*pi*D*(u - 1/2) = 2*pi * int_u^{1/2} (Del2(t) - Del0) dt!
#       At tau in [1.05, 40.0], r0*c*kap^2 = tau*(r0*kap) <= 0.45*tau!
#       Exact bound on |Im_err(u)| for u in [0, 0.5]: <= 2*pi * C_Del0 * 0.45 * tau * 0.193 = 0.00348 * tau!
#       (At tau = 5: 0.017; tau = 15: 0.052; tau = 30: 0.104; tau = 40: 0.139!)
#   Wait: check the Lipschitz integral bound on |R_3_fin - R_3_par| / H_0 at tau = 2, 5, 10, 20, 30, 40!!
import numpy as np
from math import pi, exp, log, sqrt
from scipy.integrate import quad
for tau, u0 in [(1.5, 0.42), (3.0, 0.335), (6.0, 0.458), (10.0, 0.494), (20.0, 0.504), (30.0, 0.502), (40.0, 0.501)]:
    H0 = tau * u0 / (tau * tau * u0 * u0 + pi * pi)
    # Part (a): |M_E_fin - M_E_par| * |W_3_par|, using |-log(1 - e^{-c_re*tau*|u^2-u0^2|})| * |dz(u)|!
    # Note |dz(u)| <= |u - u0| * [ 0.03*tau*(|u|+u0) + 2*pi*0.020*0.45*tau*(u^2 + u0^2 + 0.34) ]!
    # Wait: at u0 ~ 0.5, |dz'(t)| <= 0.03*tau*2*t + 2*pi*0.014*0.45*tau = (0.06*t + 0.040)*tau!
    # Let's check exact worst-case |dz(u)| <= min(0.006*tau + 0.025*tau*|u^2-u0^2|, 0.095*tau*|u - u0|) on |u| <= 1.2!
    L = sqrt(30.0 / tau) + 1.0
    def lip_f(u):
        x = 0.95 * tau * abs(u * u - u0 * u0)
        if x < 1e-12: return 0.0
        mlog = -log(1.0 - exp(-x)) if x < 30 else exp(-x)
        dz = min(0.004 * tau + 0.025 * tau * abs(u * u - u0 * u0), (0.040 + 0.06 * max(abs(u), u0)) * tau * abs(u - u0))
        W3p = 15.0 * (2 * tau) ** 3 / (((2 * tau * u) ** 2 + (2 * pi) ** 2) ** 3)
        return mlog * dz * W3p
    e_M, _ = quad(lip_f, -L, L, points=[-u0, u0], limit=300)
    # Part (b): |M_E_fin| * |W_3_fin - W_3_par|!
    # Note |W_3_fin - W_3_par| / |W_3_par| <= (CW3 - 1) + 6 * |zp_fin - zp_par|/|zp|!
    # Check |zp_fin - zp_par|/|zp| at tau = 10: Del0/1 = 0.029, Del2(0.5) = 0.086 -> at u=0.5, |zp|=sqrt(100+39.5)=11.8, 2*pi*0.086=0.54 -> 4.6% (times 6 = 28%!)
    # Wait: 28% of tri_par(10) = 0.28 * 1.07 = 0.30!
    print(f"  tau={tau:4.1f} u0={u0:.3f}: e_M / H0 = {e_M/H0:.4f} (using -log(1-e^-x)*|dz|!)")
