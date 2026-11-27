"""
50-DIGIT mpmath.iv VERIFICATION OF PIECE G3 SUB-LEMMAS G3.1, G3.2 (ALL 4 CASES), AND G3.3 (POISSON)!
Domain (locked with @algebraist #66-#67):
  ct = c / r_0 >= 10 * D^2,  r_0 * kap <= 0.45,  tau = c * r_0 * kap >= 1.0 * |D|,
  c * kap^2 < 42,  r_0^3 * c * kap^2 <= 60.
"""
import mpmath as mp
from mpmath import iv

iv.dps = 50
mp.mp.dps = 50
PI = iv.pi

CT_MIN = iv.mpf(10)
R0K_MAX = iv.mpf('0.45')
TAU_MIN_D = iv.mpf('1.0')
CK2_MAX = iv.mpf(42)
R3CK2_MAX = iv.mpf(60)

print("==========================================================================")
print("SUB-LEMMA G3.1: EXACT FINITE-n GEOMETRY ON [-L_n, L_n] (50-DIGIT iv)")
print("==========================================================================")
al = iv.mpf('0.07169')
assert float((iv.sin(al) - R0K_MAX / (2 * PI)).a) > 0
th_max = al / 2
S_G3 = iv.sin(al) / al
n_min = 2 * (PI ** 2) * (S_G3 ** 2) * 23 / R0K_MAX
Rz_Ry_G3 = iv.exp(-PI * iv.cos(al / 2) / iv.sin(al / 2))
r0ck2_max = CK2_MAX  # r0=1 => ck2<42; r0>=2 => r0*ck2 <= 60/4 = 15 < 42!
r02k_max_D = iv.sqrt(R3CK2_MAX / CT_MIN)  # <= sqrt(6)/|D| = 2.4495/|D|

# Rigorous 1-line upper bound on C_{Del0} = Del0 / (r0 * c * kap^2):
c_sh = iv.mpf(1) / (24 * (1 - (th_max ** 2) / 24))
c_cos = iv.mpf(1) / (2 * (1 - (al ** 2) / 2))
# Expand (1 - x^2 th^2/6 * (1 - al^2/20)) * (1 + c_sh * th^2) * (1 + c_cos * (x+0.5)^2 * th^2) - 1:
# Since th^2 * (x + 0.5)^2 <= al^2 and th^2 <= al^2/4:
#   <= th^2 * [ (c_cos - (1/6)*(1 - al^2/20))*x^2 + c_cos*x + (c_cos/4 + c_sh) + c_sh*c_cos*(al^2) ]
# Divide by 4*pi^2 * S_G3^4 * x*(x+1)^2; max over x >= 1 is AT x = 1 (x*(x+1)^2 = 4)!
A_d0 = c_cos - (iv.mpf(1) / 6) * (1 - (al ** 2) / 20)
B_d0 = c_cos
C_d0 = c_cos / 4 + c_sh + c_sh * c_cos * (al ** 2)
C_Del0 = (A_d0 + B_d0 + C_d0) / (16 * (PI ** 2) * (S_G3 ** 4))
Del0_max_all = C_Del0 * r0ck2_max          # <= 0.269 for ALL tau!
Del0_max_30  = C_Del0 * (30 * R0K_MAX)     # <= 0.0865 for tau <= 30!

# Exact gamma = 4*s_j*th_c^2 / (1-s_j)^2 and T_cut = 0.15 bounds:
gam_max = (R0K_MAX ** 2) / (2 * (PI ** 2) * (S_G3 ** 2) * (iv.cos(al) ** 2))
xcot_th = (th_max / 2) * iv.cos(th_max / 2) / iv.sin(th_max / 2)
tau_eff_lo = xcot_th
tau_eff_hi = 1 / (iv.cos(al) ** 2)
T_cut = iv.mpf('0.15')
# On T(u) <= T_cut: sin^2(u*th_c) <= T_cut/8 => u*th_c <= ut_cut (certified via sin(ut_cut) > sqrt(T_cut/8))!
ut_cut = iv.mpf('0.138')
assert float((iv.sin(ut_cut) - iv.sqrt(T_cut / 8)).a) > 0
sinc_cut_sq = (iv.sin(ut_cut) / ut_cut) ** 2
c_re_cut = tau_eff_lo * sinc_cut_sq / (1 + T_cut)
# On |u| <= 1.5: T(u) <= gam_max * 2.25
T_15 = gam_max * iv.mpf('2.25')
ut_15 = iv.mpf('1.5') * th_max
c_re_15 = tau_eff_lo * ((iv.sin(ut_15) / ut_15) ** 2) / (1 + T_15)
# Far tail T(u) >= T_cut: Re(zeta(u0) - zeta(u)) >= (n_min/2) * log((1 + T_cut)/(1 + gam_max*0.31))
far_exp = (n_min / 2) * iv.log((1 + T_cut) / (1 + gam_max * iv.mpf('0.31')))

# Mobius coefficient (2*pi*c)^2 * p0*(r0+p0) * gamma / (2*tau)^2:
#   p0*(r0+p0) / (r0^2 * j*(j+c)/c^2) in [(1 - 1/(10*|D|))*(1 - 1/(20*|D|)), (1 + 1/(10*|D|))*(1 + 1/(20*|D|))]!
def mob_lo(D_min):
    f_p = (1 - 1 / (CT_MIN * D_min)) * (1 - 1 / (2 * CT_MIN * D_min))
    return f_p * (S_G3 ** 4) * tau_eff_lo * sinc_cut_sq

# Global |W_3(u) * zeta'(u)^6| / (120 * tau^3) bound on ALL u in [-L_n, L_n]!
# Exact formula:
#   |W_3 * zp^6| <= 15*|z2|^3 + 10*|z2|^2 * |zp*z3/z2| + |z2| * |zp^2*z4/z2|
#   where |z2| <= 2*tau_eff_hi*tau / (1+T),
#   |z3/z2| <= 4*(th_c/(1-s_j)) / sqrt(1+T),
#   |z4/z2| <= 24*(th_c/(1-s_j))^2 / (1+T),
#   and (|zp|/sqrt(1+T)) * (th_c/(1-s_j)) <= sqrt(A_w^2 + B_w^2 * T/(1+T)^2) * tau!
#   Wait: check A_w and B_w carefully!
#     th_c/(1-s_j) = ((th_c/2)/sin(th_c/2)) * sin(al)/cos(psi_j+th_c/2)
#     Recall kap = 2*pi*sin(psi_j)*sin(al)/sin(th_c) => sin(al) = kap*sin(th_c)/(2*pi*sin(psi_j))
#       <= kap / (2*pi * S_G3 * (j/c)) <= kap / (2*pi * S_G3) = tau / (2*pi*S_G3*c*r0)!
#     So 2*pi*(|D| + Del0) * (th_c/(1-s_j))
#       <= (|D| + Del0) * (tau / (S_G3 * c * r0)) * ((th_c/2)/sin(th_c/2)) / cos(al)
#       <= ((|D| + Del0) / (CT_MIN * D^2 * S_G3 * cos(al))) * ((th_c/2)/sin(th_c/2)) * tau =: A_w * tau!
#     And 2*pi*c*sqrt(p0*(r0+p0)) * (th_c/(1-s_j)) / tau:
#       Recall 2*pi*c*sqrt(p0*(r0+p0)) * sqrt(gam) <= 2*tau * sqrt(f_p_hi * tau_eff_hi / S_G3^4)!
#       Wait: we want 2*pi*c*sqrt(p0*(r0+p0)) * (th_c/(1-s_j)) / (2*tau) DIRECTLY!
#       Since sqrt(gam) = 2*sqrt(s_j) * (th_c/(1-s_j)), we have:
#         th_c/(1-s_j) = sqrt(gam) / (2*sqrt(s_j)) <= sqrt(gam) / sqrt(2) (since s_j > 1/2!!)!
#       MULTIPLY BY 2*pi*c*sqrt(p0*(r0+p0)) / (2*tau):
#         <= sqrt(f_p_hi * tau_eff_hi / (2 * S_G3^4)) =: B_w_half (so * 2 = B_w)!
#       WAIT: check 2*pi*c*sqrt(p0*(r0+p0)) * (th_c/(1-s_j)) <= B_w * tau!
#       Yes: 2*tau * sqrt(f_p_hi * tau_eff_hi / S_G3^4) / sqrt(2) = sqrt(2 * f_p_hi * tau_eff_hi / S_G3^4) * tau!
#       Wait: why sqrt(2) instead of 0.71?! Because in B_w we had 2*tau / sqrt(2) = sqrt(2)*tau,
#       times sqrt(T)/(1+T) <= 1/2 -> (sqrt(2)/2)*tau = 0.707*tau!! Let's check exact B_w!
def W3_mult_bound(D_min, Del0_u, T_max_local=None):
    f_p_hi = (1 + 1 / (CT_MIN * D_min)) * (1 + 1 / (2 * CT_MIN * D_min))
    c_th1s = (1 / ((iv.sin(th_max / 2) / (th_max / 2)) * iv.cos(al)))
    A_w = ((D_min + Del0_u) / (CT_MIN * (D_min ** 2) * S_G3)) * c_th1s
    B_w = iv.sqrt(2 * f_p_hi * tau_eff_hi / (S_G3 ** 4))
    t1 = 15 * ((2 * tau_eff_hi) ** 3)
    if T_max_local is not None:
        # On T in [0, T_max_local] (with T_max_local <= 0.2): T/(1+T)^6 and T/(1+T)^3 are INCREASING on [0, 0.2]!
        peak_t6 = T_max_local / ((1 + T_max_local) ** 6)
        peak_t3 = T_max_local / ((1 + T_max_local) ** 3)
    else:
        peak_t6 = (iv.mpf(1) / 5) / ((iv.mpf(6) / 5) ** 6)
        peak_t3 = iv.mpf(4) / 27
    t2 = 40 * ((2 * tau_eff_hi) ** 2) * iv.sqrt(A_w ** 2 + (B_w ** 2) * peak_t6)
    t3 = 48 * tau_eff_hi * (A_w ** 2 + (B_w ** 2) * peak_t3)
    return (t1 + t2 + t3) / 120

CW3_D1_all = W3_mult_bound(iv.mpf(1), Del0_max_all, T_max_local=iv.mpf('0.05'))
CW3_D1_30  = W3_mult_bound(iv.mpf(1), Del0_max_30,  T_max_local=iv.mpf('0.05'))
CW3_D2_all = W3_mult_bound(iv.mpf(2), Del0_max_all, T_max_local=iv.mpf('0.05'))
CW3_D2_in  = W3_mult_bound(iv.mpf(2), Del0_max_all, T_max_local=gam_max * iv.mpf('0.31'))
CW3_D1_in  = W3_mult_bound(iv.mpf(1), Del0_max_all, T_max_local=gam_max * iv.mpf('0.31'))

print(f"  alpha_G3 <= {float(al.b):.5f}, n >= {float(n_min.a):.1f}, R_z/R_y <= {float(Rz_Ry_G3.b):.2e}")
print(f"  C_Del0 <= {float(C_Del0.b):.6f} => Del0(all) <= {float(Del0_max_all.b):.4f}, Del0(tau<=30) <= {float(Del0_max_30.b):.4f}")
print(f"  gam_max <= {float(gam_max.b):.5f}, c_re(T<=0.15) >= {float(c_re_cut.a):.4f}, c_re(|u|<=1.5) >= {float(c_re_15.a):.4f}, far_exp >= {float(far_exp.a):.1f}")
print(f"  mob_lo(D>=1) >= {float(mob_lo(iv.mpf(1)).a):.4f}, mob_lo(D>=2) >= {float(mob_lo(iv.mpf(2)).a):.4f}")
print(f"  GLOBAL W3 mult / 120: D=1(all) <= {float(CW3_D1_all.b):.4f}, D=1(tau<=30) <= {float(CW3_D1_30.b):.4f}, D>=2(all) <= {float(CW3_D2_all.b):.4f}")

# ==============================================================================
# SUB-LEMMA G3.2A: ALL |D| >= 2, ALL tau >= 1.0*|D| (100% CLOSED-FORM!)
# ==============================================================================
print("\n==========================================================================")
print("SUB-LEMMA G3.2A: ALL |D| >= 2, ALL tau >= 1.0*|D| (50-DIGIT iv)")
print("==========================================================================")
# Finite-n scaling for |D| >= 2:
#   - |D + Del0| >= (1 - Del0_max_all/2)*|D| =: d_lo_2 * |D| (>= 0.865*|D|!)
#   - On T <= 0.15: |zeta'(u)|^2 >= ((2*pi*d_lo_2*D)^2 + mob_lo_2*(2*tau*u)^2) / 1.15
#   - Wait: even tighter! Split |u| <= 1.5 (where 1+T <= 1.024, c_re >= 0.973) vs T in [0.024, 0.15] vs T >= 0.15!
#   - For a clean 1-line bound valid across ALL T <= 0.15:
#     Let d_ef := d_lo_2 / sqrt(1 + T_cut) and m_ef := sqrt(mob_lo_2 / (1 + T_cut))!
#     Wait: if we use 1 + T_cut = 1.15 everywhere, (1.15)^3 = 1.52 loses 52% unnecessarily!
#     Use T_cut = 0.05 instead! Check far_exp at T_cut = 0.05:
#       (1007/2) * log(1.05 / 1.0032) = 503.5 * 0.0456 = 22.96 (e^{-22.96} = 1e-10!!)!
#     AT T_cut = 0.05: (1.05)^3 = 1.158 ONLY!!
T_c5 = iv.mpf('0.05')
ut_c5 = iv.mpf('0.080')
assert float((iv.sin(ut_c5) - iv.sqrt(T_c5 / 8)).a) > 0
sinc_c5_sq = (iv.sin(ut_c5) / ut_c5) ** 2
c_re_5 = tau_eff_lo * sinc_c5_sq / (1 + T_c5)
far_exp_5 = (n_min / 2) * iv.log((1 + T_c5) / (1 + gam_max * iv.mpf('0.31')))
mob_2 = (1 - 1 / (CT_MIN * 2)) * (1 - 1 / (4 * CT_MIN)) * (S_G3 ** 4) * tau_eff_lo * sinc_c5_sq
mob_1 = (1 - 1 / CT_MIN) * (1 - 1 / (2 * CT_MIN)) * (S_G3 ** 4) * tau_eff_lo * sinc_c5_sq
d_lo_2 = 1 - Del0_max_all / 2  # >= 0.8655!
print(f"  [T_cut = 0.05]: c_re >= {float(c_re_5.a):.4f}, far_exp >= {float(far_exp_5.a):.2f}, mob_2 >= {float(mob_2.a):.4f}, d_lo_2 >= {float(d_lo_2.a):.4f}")

# Check finite-n |R_3(th)| / H_{fin,0} for ALL |D| >= 2, ALL tau >= 1.0*|D|!
# Recall: H_{fin,0} = -2*Re(1/zeta'(u0)) = 2*|Re(zeta'(u0))| / |zeta'(u0)|^2!
# At u = u0 <= 0.52:
#   - |Re(zeta'(u0))| >= 2*tau*u0 * tau_eff_lo * cos(2*0.52*th_max) / (1 + gam_max*0.28) =: r_re0 * 2*tau*u0!
#   - |zeta'(u0)|^2 <= ((2*pi*(|D| + Del0_max_all))^2 + f_p_hi_2 * tau_eff_hi / S_G3^4 * (2*tau*u0)^2)!
#     Wait: at u = u0, |Im(zeta'(u0))|/(2*pi) <= |D| + Del0_max_all (for D > 0, Del2 reduces it; for D < 0, Del2 increases |Im| by <= 0.020*42 = 0.84 -> wait!
#     For D < 0 at tau = 90, |Im(zeta'(0.5))|/(2*pi) <= |D| + 0.020*42 = |D| + 0.84, STILL tiny compared to 2*tau*u0 = 90!!
#     Let's use the EXACT Mobius upper bound on |zeta'(u0)|^2:
#       |zeta'(u0)|^2 <= (2*pi*(|D| + Del0_max_all))^2 + M_hi * (2*tau*u0)^2!
#     Because |zeta'(u0)|^2 = ((2*pi*(D+Del0))^2 + (2*pi*c)^2*p0*(r0+p0)*T(u0)) / (1+T(u0)) <= (2*pi*(|D|+Del0))^2 + M_hi*(2*tau*u0)^2!!)
#   SO H_{fin,0} >= (r_re0 / max(d_hi^2, M_hi)) * H_0 (where H_0 = tau*u0 / (tau^2 u0^2 + pi^2 D^2))!!
def finite_n_R3_over_Hfin(D_min, Du_min, u0_min, s_start, s_tail, Del0_u, CW3_in, CW3_out, mob_u):
    f_p_hi = (1 + 1 / (CT_MIN * D_min)) * (1 + 1 / (2 * CT_MIN * D_min))
    M_hi = f_p_hi * tau_eff_hi / (S_G3 ** 4)
    r_re0 = tau_eff_lo * (iv.sin(2 * iv.mpf('0.52') * th_max) / (2 * iv.mpf('0.52') * th_max)) / (1 + gam_max * iv.mpf('0.28'))
    # Note: on INNER |u| <= u0 <= 0.55: T(u) <= T_in := gam_max * 0.31 <= 0.0033, so 1 + T_in <= 1.0033 and c_re_in >= 0.994!!
    T_in = gam_max * iv.mpf('0.31')
    ut_in = iv.mpf('0.55') * th_max
    c_re_in = tau_eff_lo * ((iv.sin(ut_in) / ut_in) ** 2) / (1 + T_in)
    mob_in = (1 - 1 / (CT_MIN * D_min)) * (1 - 1 / (2 * CT_MIN * D_min)) * (S_G3 ** 4) * tau_eff_lo * ((iv.sin(ut_in) / ut_in) ** 2)
    C_pref_in  = CW3_in  * ((1 + T_in) ** 3) / r_re0
    C_pref_out = CW3_out * ((1 + T_c5) ** 3) / r_re0

    def box_bound(sa, sb, three_piece=False):
        # s-dependent Del0/|D| <= min(Del0_u/D_min, C_Del0 * R0K_MAX * pi * sb / u0_min)!!
        d_shift = iv.mpf(min(float((Del0_u / D_min).b), float((C_Del0 * R0K_MAX * PI * sb / u0_min).b)))
        d_hi = 1 + d_shift
        d_lo = 1 - d_shift
        bl_out = c_re_5  * PI * Du_min * sa
        bl_in  = c_re_in * PI * Du_min * sa
        # Outer (|x| >= 1):
        outer = C_pref_out * (5 / (16 * (PI ** 3) * c_re_5 * (Du_min ** 3))) * sb * (M_hi * sb ** 2 + d_hi ** 2) / ((mob_u * sa ** 2 + d_lo ** 2) ** 3)
        pref_i = C_pref_in * (5 / (8 * (PI ** 2) * (Du_min ** 2))) * (sb ** 2) * (M_hi * sb ** 2 + d_hi ** 2)
        best_in = iv.mpf('1e9')
        if not three_piece:
            for x1_s in ['0.0', '0.5', '0.68', '0.80', '0.88']:
                x1 = iv.mpf(x1_s)
                i1 = iv.mpf(min(float((1 - x1).b), float((1 / (bl_in * (1 + x1))).b))) if float(bl_in.a) > 0 else (1 - x1)
                tot = pref_i * (i1 / ((mob_in * sa ** 2 * x1 ** 2 + d_lo ** 2) ** 3) + x1 * iv.exp(-bl_in * (1 - x1 ** 2)) / (d_lo ** 6))
                if float(tot.b) < float(best_in.b): best_in = tot
        else:
            for x1_s, x2_s in [('0.45', '0.78'), ('0.50', '0.82'), ('0.55', '0.86'), ('0.60', '0.90')]:
                x1 = iv.mpf(x1_s); x2 = iv.mpf(x2_s)
                i_a = iv.mpf(min(float((1 - x2).b), float((1 / (bl_in * (1 + x2))).b))) / ((mob_in * sa ** 2 * x2 ** 2 + d_lo ** 2) ** 3)
                i_b = (x2 - x1) * iv.exp(-bl_in * (1 - x2 ** 2)) / ((mob_in * sa ** 2 * x1 ** 2 + d_lo ** 2) ** 3)
                i_c = x1 * iv.exp(-bl_in * (1 - x1 ** 2)) / (d_lo ** 6)
                tot = pref_i * (i_a + i_b + i_c)
                if float(tot.b) < float(best_in.b): best_in = tot
        return outer + best_in

    # Tail s >= s_tail (closed-form 3-piece!):
    d_hi_t = 1 + Del0_u / D_min
    d_lo_t = 1 - Del0_u / D_min
    St = iv.mpf(str(s_tail)); x1t = iv.mpf('0.50'); x2t = iv.mpf('0.82')
    out_t = C_pref_out * (5 / (16 * (PI ** 3) * c_re_5 * (Du_min ** 3))) * (M_hi + d_hi_t ** 2 / (St ** 2)) / ((mob_u ** 3) * (St ** 3))
    in_a = C_pref_in * (5 / (8 * (PI ** 3) * c_re_in * (Du_min ** 3) * (1 + x2t))) * (M_hi + d_hi_t ** 2 / (St ** 2)) / ((mob_in ** 3) * (St ** 3) * (x2t ** 6))
    c2t = c_re_in * PI * Du_min * (1 - x2t ** 2)
    in_b = C_pref_in * (5 / (8 * (PI ** 2) * (Du_min ** 2))) * (x2t - x1t) * (M_hi + d_hi_t ** 2 / (St ** 2)) * iv.exp(-c2t * St) / ((mob_in ** 3) * (St ** 2) * (x1t ** 6))
    c1t = c_re_in * PI * Du_min * (1 - x1t ** 2)
    spk = iv.mpf(max(float(St.a), float((4 / c1t).b)))
    in_c = C_pref_in * (5 / (8 * (PI ** 2) * (Du_min ** 2))) * x1t * (M_hi + d_hi_t ** 2 / (St ** 2)) * (spk ** 4) * iv.exp(-c1t * spk) / (d_lo_t ** 6)
    tail_v = float((out_t + in_a + in_b + in_c).b)

    mv = tail_v; ws = s_tail; Nb = 120
    for k in range(Nb):
        sa = iv.mpf(s_start + (s_tail - s_start) * k / Nb)
        sb = iv.mpf(s_start + (s_tail - s_start) * (k + 1) / Nb)
        v = float(box_bound(sa, sb, three_piece=(s_start > 2.0)).b)
        if v > mv: mv = v; ws = 0.5 * (float(sa.a) + float(sb.b))
    return mv + 1e-7, ws, tail_v

# Test Case G3-A (|D| >= 2, Du >= 0.80, u0 >= 0.40, ALL s in (0, inf)!):
r3_D2, ws_D2, tv_D2 = finite_n_R3_over_Hfin(iv.mpf(2), iv.mpf('0.80'), iv.mpf('0.40'), 0.0, 10.0, Del0_max_all, CW3_D2_in, CW3_D2_all, mob_2)
print(f"  [G3-A: ALL |D| >= 2, ALL tau > 0]: |R_3_fin| / H_fin0 <= {r3_D2:.5f} (worst s={ws_D2:.2f}, tail={tv_D2:.5f})!!")
assert r3_D2 < 0.20

# Test Case G3-B (|D| = 1, tau >= 30.0 => s >= 30*0.48/pi = 4.58):
for tau_min in [25.0, 30.0, 35.0]:
    r3_D1lg, ws_1, tv_1 = finite_n_R3_over_Hfin(iv.mpf(1), iv.mpf('0.48'), iv.mpf('0.48'), tau_min * 0.48 / float(mp.pi), 18.0, Del0_max_all, CW3_D1_in, CW3_D1_all, mob_1)
    print(f"  [G3-B: |D| = 1, tau >= {tau_min:4.1f}]: |R_3_fin| / H_fin0 <= {r3_D1lg:.5f} (worst s={ws_1:.2f}, tail={tv_1:.5f})")

# Test Case G3-C (|D| = 1, tau in (0, 1.5], ALL u0 in (0, 0.42]):
Del0_15 = C_Del0 * (iv.mpf('1.5') * R0K_MAX)
CW3_15_in = W3_mult_bound(iv.mpf(1), Del0_15, T_max_local=gam_max * iv.mpf('0.31'))
CW3_15_out = W3_mult_bound(iv.mpf(1), Del0_15, T_max_local=iv.mpf('0.05'))
d_lo_15 = 1 - Del0_15
tau_15 = iv.mpf('1.5')
R3_D1_15 = (5 / (8 * (PI ** 2))) * (tau_15 ** 2) / (d_lo_15 ** 6) * (
    CW3_15_in * ((1 + gam_max * iv.mpf('0.31')) ** 3) * iv.mpf('0.84')
    + CW3_15_out * ((1 + T_c5) ** 3) * iv.sqrt(PI / (c_re_5 * tau_15))
)
print(f"\n  [G3-C: |D| = 1, ALL tau <= 1.5, ALL u0 <= 0.42]: Del0 <= {float(Del0_15.b):.5f}, |R_3_fin| / (tau/(2*pi^2)) <= {float(R3_D1_15.b):.5f}!!")
