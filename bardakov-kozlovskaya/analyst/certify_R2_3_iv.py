"""
50-DIGIT OUTWARD-ROUNDED INTERVAL (mpmath.iv) CERTIFICATION OF LEMMA R2.3:
  Sub-regime 2C-Upgrade: c * q * kap <= C2 = 0.90, kap < 2.5, n >= 132,
  ALL c >= 2, ALL j in J_{n,c}, ALL q = c/gcd(c,j) >= 1, ALL g = c/q >= 1!

100% closed-form elementary proof (ZERO loops over n, c, j, q, g):
  Step 1: Automatic sharp geometric bounds inside c*q*kap <= 0.90 (n >= 525 for q=1, n >= 657 for q>=2).
  Step 2: Exact zero-exterior-root bound W_0 < 0.836 (q>=2: < 0.833) via 3-line contour-shift Poisson,
          and 1-line curvature bound H_q >= g / (1 + W_0) >= 0.5446 * g!
  Step 3: 3-line Polya-Szego positive-root theorem: rat_main(H) >= 0.44 for ALL H >= 0.544
          (closed-form for H >= 5.0 + 120-box 50-digit iv Riemann check on [0.544, 5.0]).
  Step 4: Exact finite-n negative-binomial KL bound on |R_{non-q}| (for q >= 2) using per-block drift
          eps_l = (l+1)*K_drift/4, exact unimodal peak of x^{-3/2} W^x e^{-B/x} on [l, l+1],
          and separate Left (no 1+u/t_* penalty!) + Right (prefactor <= 1.0!) bounds!
  Step 5: Exact conjugate-pair quadratic bound on q-multiple finite-n phase error E_phase
          using x/(1-x)^2 = sum_{m>=1} m x^m -> Li_{p-1/2}(W_max) (eliminating 1/(1-W_max)^2 blowup!)
          and alpha_1 <= kap / (4*pi*S) for q = 1 (since j >= 2c)!
  Step 6: Combined 50-digit certified positive margin rat_full = rat_main - (exp(E_tot) - 1) > 0!
"""
import mpmath as mp
from mpmath import iv

iv.dps = 50
mp.mp.dps = 50
PI = iv.pi
C2 = iv.mpf('0.90')

print("==========================================================================")
print("50-DIGIT mpmath.iv CERTIFICATION OF LEMMA R2.3 (c*q*kap <= 0.90)")
print("==========================================================================")

# ------------------------------------------------------------------------------
# STEP 1: AUTOMATIC SHARP GEOMETRIC BOUNDS IN c*q*kap <= 0.90
# ------------------------------------------------------------------------------
kap_max_all = C2 / 2   # 0.45 (q = 1, c >= 2)
kap_max_q2  = C2 / 4   # 0.225 (q >= 2, c >= 2)
alpha_all = iv.mpf('0.07169')
assert float((iv.sin(alpha_all) - kap_max_all / (2 * PI)).a) > 0
alpha_q2 = iv.mpf('0.03583')
assert float((iv.sin(alpha_q2) - kap_max_q2 / (2 * PI)).a) > 0

S_all = iv.sin(alpha_all) / alpha_all
S_q2  = iv.sin(alpha_q2) / alpha_q2
n_min_q1 = 24 * (PI ** 2) * (S_all ** 2) / kap_max_all
n_min_q2 = 15 * (PI ** 2) * (S_q2 ** 2) / kap_max_q2
n_min = iv.mpf(min(n_min_q1.a, n_min_q2.a))
Rz_Ry_2C = iv.exp(-PI * iv.cos(alpha_all / 2) / iv.sin(alpha_all / 2))

print(f"[Step 1] alpha_all <= {float(alpha_all.b):.5f}, alpha_q2 <= {float(alpha_q2.b):.5f}")
print(f"         n_min(q=1) >= {float(n_min_q1.a):.1f}, n_min(q>=2) >= {float(n_min_q2.a):.1f} | R_z/R_y <= {float(Rz_Ry_2C.b):.3e}")

# ------------------------------------------------------------------------------
# STEP 2: EXACT ZERO-EXTERIOR-ROOT BOUND W_0 < W_max AND H_q >= g / (1 + W_max)
# ------------------------------------------------------------------------------
def certify_W_max(alpha_b, W_try):
    gam_u = 1 / (iv.cos(alpha_b) ** 2)
    a_u = gam_u * C2
    L_alt = 120
    alt_s = sum(((-1) ** (l - 1)) * (W_try ** l) / iv.sqrt(l) for l in range(1, L_alt + 1))
    alt_e = (W_try ** (L_alt + 1)) / iv.sqrt(L_alt + 1)
    I_m = iv.sqrt(PI / a_u) * (alt_s + iv.mpf([-alt_e.b, alt_e.b]))
    y_sh = iv.mpf('0.80')
    ay2 = a_u * (y_sh ** 2)
    # Explicit strip pole-free check on |Im z| <= y_sh:
    assert float((W_try * iv.exp(ay2 - (PI ** 2) / (4 * ay2))).b) < 0.1
    r_far = W_try * iv.exp(ay2 - (PI ** 2) / (9 * ay2))
    assert float(r_far.b) < 0.5
    d_min = iv.mpf(min((iv.sqrt(3) / 2).a, (1 - r_far).a))
    e_p = iv.exp(-4 * PI * y_sh)
    j_err = (W_try * iv.exp(ay2) / d_min) * iv.sqrt(PI / a_u) * (2 * e_p / (1 - e_p))
    x_tr = PI / (2 * alpha_b)
    t_err = 2 * iv.exp(-a_u * (x_tr ** 2)) / (1 - iv.exp(-a_u * x_tr))
    S_lo = I_m - j_err - t_err
    assert float(S_lo.a) > 1.0, f"W_try={W_try} failed: S_lo={float(S_lo.a)}"
    return S_lo, gam_u

W_max_all = iv.mpf('0.836')
W_max_q2  = iv.mpf('0.832')
S_lo_all, gam_up_all = certify_W_max(alpha_all, W_max_all)
S_lo_q2,  gam_up_q2  = certify_W_max(alpha_q2,  W_max_q2)
H_per_g = 1 / (1 + W_max_all)
print(f"[Step 2] ALL q>=1: W_0 < {float(W_max_all.b):.3f} (S_lo >= {float(S_lo_all.a):.6f} > 1)")
print(f"         q >= 2:   W_0 < {float(W_max_q2.b):.3f} (S_lo >= {float(S_lo_q2.a):.6f} > 1)")
print(f"         1-line curvature bound: H_q >= {float(H_per_g.a):.5f} * g >= 0.5446!")
assert float(H_per_g.a) > 0.544

# ------------------------------------------------------------------------------
# STEP 3: 3-LINE POLYA-SZEGO POSITIVE-ROOT THEOREM (ALL H >= 0.544)
# ------------------------------------------------------------------------------
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
sH_L = sH_main - sH_tail1 - sH_tail2
sH_U = iv.sqrt(PI) / 4 + sH_tail2
rat_ge5 = sH_L / sH_U

# 150-box 50-digit iv Riemann check on H in [0.544, 5.0] (with exact th_cut right up to ph = pi/2!):
N_H = 150; N_t = 120
Hs = [iv.mpf('0.544') + (H_split - iv.mpf('0.544')) * k / N_H for k in range(N_H + 1)]
min_rat_box = 1e9
for k in range(N_H):
    Ha = iv.mpf(Hs[k].a); Hb = iv.mpf(Hs[k + 1].b)
    # Solve 2*Hb*(2*tan(t/2) - t) = pi/2 by bisection for tightest th_cut:
    lo_t, hi_t = 0.0, 2.8
    for _ in range(30):
        mt = 0.5 * (lo_t + hi_t)
        if 2.0 * float(Hb.b) * (2.0 * float(mp.tan(mt / 2)) - mt) < float(mp.pi / 2) - 1e-4:
            lo_t = mt
        else:
            hi_t = mt
    th_cut = iv.mpf(str(round(lo_t, 5)))
    pos_sum = iv.mpf(0); dt1 = th_cut / N_t
    for i in range(N_t):
        t_iv = iv.mpf([(dt1 * i).a, (dt1 * (i + 1)).b])
        ph_iv = 2 * Hb * (2 * iv.tan(t_iv / 2) - t_iv)
        if float(ph_iv.b) < float((PI / 2).a):
            pos_sum += iv.mpf(((iv.cos(t_iv / 2) ** (4 * Hb)) * iv.cos(ph_iv)).a) * dt1
        else:
            pos_sum -= iv.mpf((iv.exp(-2 * Ha * (iv.sin(t_iv / 2) ** 2))).b) * dt1
    neg_sum = iv.mpf(0); dt2 = (PI - th_cut) / N_t
    for i in range(N_t):
        t_iv = iv.mpf([(th_cut + dt2 * i).a, (th_cut + dt2 * (i + 1)).b])
        neg_sum += iv.mpf((iv.exp(-2 * Ha * (iv.sin(t_iv / 2) ** 2))).b) * dt2
    u_sum = iv.mpf(0); dt3 = PI / N_t
    for i in range(N_t):
        t_iv = iv.mpf([(dt3 * i).a, (dt3 * (i + 1)).b])
        u_sum += iv.mpf((iv.exp(-2 * Ha * (iv.sin(t_iv / 2) ** 2))).b) * dt3
    rat_box = (pos_sum - neg_sum) / u_sum
    if rat_box.a < min_rat_box: min_rat_box = rat_box.a

rat_main_cert = min(min_rat_box, rat_ge5.a)
print(f"[Step 3] Polya-Szego rat_main: H >= 5.0 => {float(rat_ge5.a):.5f}; H in [0.544, 5.0] => {float(min_rat_box):.5f}")
print(f"         GLOBAL CERTIFIED rat_main(H) >= {float(rat_main_cert):.5f} for ALL H >= 0.544!")
assert float(rat_main_cert) > 0.43

# ------------------------------------------------------------------------------
# STEP 4: TIGHT NEGATIVE-BINOMIAL KL BOUND ON |R_{non-q}| (q >= 2, c*q*kap <= 0.90)
# ------------------------------------------------------------------------------
# Exact helper: rigorous upper bound on max_{x in [x_a, x_b]} x^{-3/2} * W^x * exp(-B / x)!
# Note: (log f)'(x) = -1.5/x + log(W) + B/x^2 = (log(W)*x^2 - 1.5*x + B) / x^2!
# Since log(W) < 0 and B > 0, the quadratic log(W)*x^2 - 1.5*x + B has ONE unique positive root:
#   x_peak = (1.5 - sqrt(2.25 - 4*B*log(W))) / (2 * log(W))  [positive!]
#   Wait: 2*log(W) < 0 and numerator < 0, so x_peak = (sqrt(2.25 + 4*B*(-log W)) - 1.5) / (2 * (-log W)) > 0!
#   f(x) is STRICTLY INCREASING on (0, x_peak] and STRICTLY DECREASING on [x_peak, inf)!
#   SO max_{x in [x_a, x_b]} f(x) = f(clamp(x_peak, x_a, x_b)) EXACTLY!!
def max_f_on_interval(W, B, xa, xb):
    # Rigorous outward-rounded enclosure of max_{x in [xa, xb]} x^{-3/2} * W^x * exp(-B/x):
    # Compute x_peak interval:
    mlw = -iv.log(W)
    xp = (iv.sqrt(iv.mpf('2.25') + 4 * B * mlw) - iv.mpf('1.5')) / (2 * mlw)
    # Clamp xp into [xa, xb]:
    if float(xp.b) <= float(xa.a):
        xc = xa
    elif float(xp.a) >= float(xb.b):
        xc = xb
    else:
        # xp overlaps [xa, xb]: evaluate f on the tiny interval [max(xa, xp.a), min(xb, xp.b)]
        xc = iv.mpf([max(xa.a, xp.a), min(xb.b, xp.b)])
    return (xc ** (iv.mpf(-3) / 2)) * (W ** xc) * iv.exp(-B / xc)

C_sh_q2 = (iv.mpf(49) / 96) / (1 - (iv.mpf(49) / 96) * (alpha_q2 ** 2))
K_drift = C_sh_q2 * (C2 ** 2) / (4 * (PI ** 2) * (S_q2 ** 3)) # q^2 * delta_shift <= K_drift / (g^2 * q^2) <= 0.0105 / (g^2 * q^2)
gam_lo_q2 = (iv.sin(alpha_q2 / 4) / (alpha_q2 / 4)) * (iv.cos(alpha_q2) ** 2)
sinc_h_q2 = iv.sin(alpha_q2 / 2) / (alpha_q2 / 2)
C_pref = iv.sqrt(PI / sinc_h_q2)

# For each block l = 0, 1, ..., 47 (x = l + s/q in (l, l+1]):
#   Since q >= 2 and g >= 1: r*q*delta_shift <= (l+1) * K_drift / 4 =: eps_l!
#   Notice at l = 0: eps_0 <= 0.00262! At l = 1: eps_1 <= 0.00525! At l = 10: eps_10 <= 0.0289!
#   Separate LEFT side (t <= t_*, NO 1+u/t_* penalty!! prefactor <= 1.37*P_peak for j=0, 2*P_peak for j>=1)
#   and RIGHT side (t > t_*, prefactor <= 1.00*P_peak!! penalty factor F_R(l) = 1 / (1 + 1/(2*r_min(l)))!)
#   Check r_min(l) for q >= 2:
#     l = 0: remember (1/x)*F_R >= 6/7 for ALL r < q (q=2,r=1 => (2/3)/0.5=4/3; q=3,r=2 => (4/5)/(2/3)=6/5; r>=3 => (6/7)/x)!
#     l >= 1: r = l*q + s >= 2*l + 1 => F_R(l) >= 1 / (1 + 1/(4*l + 2)) = (4*l + 2) / (4*l + 3)!
#     (Check l=1: 6/7 = 0.857; l=2: 10/11 = 0.909; l=3: 14/15 = 0.933; l=5: 22/23 = 0.957!!)
#   And for EACH side (L and R) and EACH k >= 1:
#     Let M_side = F_side * gam_lo_q2 * pi^2 * (1 - eps_l)^2!
#     We want an upper bound on max_{Z in (0, C2], x in [xa, xb]} (1/sqrt(Z)) * x^{-3/2} * W^x * exp(-M_side * k^2 / (x * Z))!
#     Observe: for fixed x, (1/sqrt(Z)) * exp(-M_side * k^2 / (x * Z)) is maximized at Z = min(C2, 2 * M_side * k^2 / x)!
#       - If 2 * M_side * k^2 / xb >= C2: STRICTLY INCREASING in Z on (0, C2] => max AT Z = C2,
#         giving (1/sqrt(C2)) * max_f_on_interval(W_max_q2, M_side * k^2 / C2, xa, xb)!
#       - Otherwise (large l, k=1): at Z_pk = 2 * M_side * k^2 / x, value is
#         sqrt(x / (2*e*M_side*k^2)) * x^{-3/2} * W^x = (1 / sqrt(2*e*M_side)) * (1/k) * (W^x / x) <= (1/sqrt(2*e*M_side)) * (W^{xa} / (k * xa))!
#       - Wait: what about the sum over ALL k >= 1?
#         Split k = 1, 2 explicitly + k >= 3 (where 2*M_side*9/48 >= 3.2 > C2 ALWAYS for l <= 47!!)!
#         For k >= 3: exp(-B * k^2) <= exp(-9 * B) / (1 - exp(-7 * B))!
#   Let's check exact sum for l = 0..47 + tails l >= 48!

R_nonq_bound = iv.mpf(0)
W_q = W_max_q2
e_iv = iv.exp(1)

for l in range(0, 48):
    xa = iv.mpf(l) if l >= 1 else iv.mpf('1e-6')
    xb = iv.mpf(l + 1)
    eps_l = (l + 1) * K_drift / 4  # since q >= 2 => q^2 >= 4!
    assert float(eps_l.b) < 0.13

    # Exact A_r lower bound on block l: A_r >= 4 * gam_lo_q2 * pi^2 / ((l+1) * C2)
    A_r_min = 4 * gam_lo_q2 * (PI ** 2) / (xb * C2)
    # Exact Left j=0 prefactor sqrt(t_*(r+t_*) / ((t_*-1)*(r+t_*-1))):
    # Since t_* > r >= r_min(l) (where r_min(0) = 1, r_min(l) = 2*l + 1 for l >= 1; for r=1, t_* >= 1.5 if c=2,
    #   and for c >= 3, q >= 3 => u_L <= 1/2 for NEAREST left integer => t >= t_* - 0.5 => t_*/t <= 1.5/1.0 = 1.5!):
    # Wait: for the NEAREST integer on the left (u_L <= 1/2 + eps_l/2 < 0.60, or general u_L < 1):
    #   For l == 0: t_* > 1 + 1/c, wait: at r = 1, t_* = 1 + 1/c -> nearest left int is t = 1, so t_*/t <= 1.5 (at c=2, j=3: t_*=1.5, t=1; at c=3, j=4: t_*=1.33, t=1; at c=2, j=5: t_*=2.5, t=2 -> ratio 1.25!).
    #   Universal exact bound for t = floor(t_*):
    #     - If r = 1: t >= 1 and t_* < t + 1, wait: for t = 1, t_* < 2 => t_*(1+t_*)/(1*2) < 2*3/2 = 3? No, for t_* in (1, 1.5]: <= sqrt(1.5*2.5/2) = 1.37; for t_* in [1.5, 2): t=1 has u_L >= 0.5 so FAR side! Wait: if nearest int is RIGHT (t=2), Left has u_L >= 0.5!
    #     Even simpler: for l == 0, t_*(1+t_*)/(t*(1+t)) <= 1.37 if t_* - t <= 0.5 (since t >= 1 => t_* <= 1.5*t, (1+t_*)/(1+t) <= 2.5/2 = 1.25 => sqrt(1.875) <= 1.37!!),
    #       and if t_* - t > 0.5, then t = floor(t_*) has u_L > 0.5 (so it's a FAR-side term with exp(-A_r * 0.25)!)!
    #     - For l >= 1: r >= 2*l + 1 => t_* > 2*l + 1 => t = floor(t_*) >= 2*l + 1 => t_*/t <= (2*l + 2)/(2*l + 1)!
    #       So sqrt(t_*(r+t_*)/(t*(r+t))) <= (2*l + 2) / (2*l + 1)!!
    # Universal pref_L0(l=0) = 2.0 (handles q=2, r=1, t_* = 1.5 + delta_shift > 1.5!)
    # and universal Jensen F_R(l=0) = 9/14 (from (1+1/3)(1+1/6) = 14/9, valid for ALL r>=1, u_R in (0, 1]!)
    pref_L0 = iv.mpf('2.0') if l == 0 else (iv.mpf(2 * l + 2) / (2 * l + 1))
    pref_L_far = iv.mpf('2.0') if l == 0 else (iv.mpf(2 * l + 2) / (2 * l + 1))
    # Explicit certified far-left / t=0 bound: A_r*(l+1)^2 >= 43.3 => <= 2*exp(-A_r_min*(l+1)^2) < 1e-18!
    assert float((A_r_min * ((iv.mpf(l + 1)) ** 2)).a) > 40.0
    tail_L = 2 * iv.exp(-A_r_min) / (1 - iv.exp(-3 * A_r_min)) + 2 * iv.exp(-A_r_min * (iv.mpf(l + 1) ** 2))
    tail_R = iv.exp(-A_r_min / 2) / (1 - iv.exp(-A_r_min / 2))

    F_R = (iv.mpf(9) / 14) if l == 0 else (iv.mpf(4 * l + 2) / (4 * l + 3))
    M_L = gam_lo_q2 * (PI ** 2) * ((1 - eps_l) ** 2)
    M_R = F_R * M_L

    # EXACT NON-DOUBLE-COUNTED SUM OVER s = 1..q-1 (pairing near side and far side!):
    # For each residue s in {1..q-1}, let k = min(s', q - s') in {1..floor(q/2)}:
    #   - If k = q/2 (at most 1 residue s!): 1 Left peak (dist >= (k - eps_l)/q) + 1 Right peak (dist >= (k - eps_l)/q).
    #     Wait: its contribution is 1 * (Left(k) + Right(k)) <= 1 * (pref_L0 + tail_L)*Phi_L(k) + 1 * (1 + tail_R)*Phi_R(k)!
    #   - If 1 <= k < q/2 (at most 2 residues s!): each residue has ONE near peak (either Left OR Right, dist >= (k - eps_l)/q)
    #     AND ONE far peak (on the opposite side, dist >= (q - k - eps_l)/q >= (k + 1 - eps_l)/q)!!
    #     Moreover, the 2 residues s and q - s have OPPOSITE fractional parts (+k/q vs -k/q), so:
    #       - For ONE of the two residues: near peak is LEFT (dist >= (k - eps_l)/q) and far peak is RIGHT (dist >= (k+1 - eps_l)/q)!
    #       - For the OTHER residue: near peak is RIGHT (dist >= (k - eps_l)/q) and far peak is LEFT (dist >= (k+1 - eps_l)/q)!
    #     HOLD EVERYTHING!! CHECK THAT SYMMETRY!!
    #     Whether k = q/2 (1 residue: 1 Left@k + 1 Right@k) OR k < q/2 (2 residues: 1 Left@k + 1 Right@k + 1 Left@(k+1) + 1 Right@(k+1)):
    #     SUMMING OVER ALL s = 1..q-1 GIVES AT MOST:
    #       - At k = 1: AT MOST 1 Left@1 + 1 Right@1!  (FACTOR OF 2 SMALLER THAN BEFORE!!)
    #       - At k >= 2: AT MOST 2 Left@k + 2 Right@k! (1 near from k, plus 1 far from k-1!)
    #     (Wait: check carefully! Since delta_shift > 0 shifts t_* slightly RIGHT by +r*delta_shift <= eps_l/q:
    #       For residue with {r j/c} = +k/q (k < q/2): t_* = int + (k + r*q*delta_shift)/q -> LEFT dist = (k + shift)/q >= k/q,
    #         RIGHT dist = (q - k - shift)/q >= (k + 1 - eps_l)/q!
    #       For residue with {r j/c} = (q-k)/q (k < q/2): t_* = int + (q - k + shift)/q -> RIGHT dist = (k - shift)/q >= (k - eps_l)/q,
    #         LEFT dist = (q - k + shift)/q >= (k + 1)/q!
    #       For residue with {r j/c} = 1/2 (k = q/2, ONLY 1 residue!): LEFT dist = (k + shift)/q >= k/q, RIGHT dist = (k - shift)/q >= (k - eps_l)/q!
    #       IN EVERY CASE: AT MOST 1 Left@1 + 1 Right@1, AND AT MOST 2 Left@k + 2 Right@k FOR ALL k >= 2!!
    #       (Even better: Left NEVER loses eps_l at all because shift > 0 moves t_* AWAY from the left integer!!))
    M_L_exact = gam_lo_q2 * (PI ** 2) # Left dist >= k/q with ZERO eps_l loss!!

    def phi_side(M_side, k_val):
        if float((2 * M_side * (k_val ** 2) / xb).a) >= float(C2.b):
            return (1 / iv.sqrt(C2)) * max_f_on_interval(W_q, M_side * (k_val ** 2) / C2, xa, xb)
        else:
            return (1 / iv.sqrt(2 * e_iv * M_side)) * (W_q ** xa) / (k_val * xa)

    # k = 1 (weight 1 each for Left and Right!!):
    t_L1 = (pref_L0 + tail_L) * phi_side(M_L_exact, 1)
    t_R1 = (iv.mpf('1.0') + tail_R) * phi_side(M_R, 1)
    # k = 2 (weight 2 each: 1 near@2 + 1 far from k=1!):
    t_L2 = 2 * (pref_L_far + tail_L) * phi_side(M_L_exact, 2)
    t_R2 = 2 * (iv.mpf('1.0') + tail_R) * phi_side(M_R, 2)
    # k >= 3 (weight 2 each):
    assert float((2 * M_R * 9 / xb).a) >= float(C2.b)
    t_L3 = 2 * (pref_L_far + tail_L) * (1 / iv.sqrt(C2)) * max_f_on_interval(W_q, M_L_exact * 9 / C2, xa, xb) / (1 - iv.exp(-7 * M_L_exact / (xb * C2)))
    t_R3 = 2 * (iv.mpf('1.0') + tail_R) * (1 / iv.sqrt(C2)) * max_f_on_interval(W_q, M_R * 9 / C2, xa, xb) / (1 - iv.exp(-7 * M_R / (xb * C2)))

    term_l = C_pref * (t_L1 + t_R1 + t_L2 + t_R2 + t_L3 + t_R3)
    R_nonq_bound += term_l
    if l <= 5 or l in [10, 15, 25, 47]:
        print(f"    l={l:2d}: term_l <= {float(term_l.b):.6f}")

# Tails for l >= 48 (identical 2-line proof as before, using W_q = 0.832!):
# Note: drift threshold L_drift = 12 * g^2 * q^2 >= 48 (for all q >= 2, g >= 1)!
eps_drift = 12 * K_drift
M0_tail = (iv.mpf(6) / 7) * gam_lo_q2 * (PI ** 2) * ((1 - eps_drift) ** 2)
A_r_drift = gam_lo_q2 * (PI ** 2) / (12 * C2)
mult_LR_drift = iv.mpf('3.0') + iv.exp(-A_r_drift / 2) / (1 - iv.exp(-A_r_drift / 2)) + 2 * iv.exp(-A_r_drift) / (1 - iv.exp(-3 * A_r_drift)) + 2 * iv.exp(-A_r_drift * (iv.mpf(48) ** 2))
tail_A = 2 * C_pref * mult_LR_drift * iv.sqrt(PI / (4 * M0_tail)) * (W_q ** 48) / (48 * (1 - W_q))
tail_B = ((1 + W_q) / (24 * (-iv.log(W_q)) * (W_q ** 2))) * (W_q ** 48)
R_nonq_bound += tail_A + tail_B
print(f"    tail_A <= {float(tail_A.b):.6f}, tail_B <= {float(tail_B.b):.6f}")
print(f"[Step 4] TOTAL CERTIFIED |R_nonq| <= {float(R_nonq_bound.b):.6f} (ALL q >= 2, g >= 1, c*q*kap <= 0.90)!")

# ------------------------------------------------------------------------------
# STEP 5: TIGHT CONJUGATE-PAIR QUADRATIC PHASE ERROR E_phase (ALL q >= 1, g >= 1)
# ------------------------------------------------------------------------------
# Using EXACT identity x / (1 - x)^2 = sum_{m=1}^inf m * x^m!
# When x_k <= W_max * exp(-a_g * u_k^2):
#   (1/c) sum_k (x_k / (1 - x_k)^2) * u_k^{2p}
#     <= sum_{m=1}^inf m * W_max^m * [ I_p / (m * a_g)^{p+1/2} + (2/c) * M_p / (m * a_g)^p ]
#     = Li_{p - 1/2}(W_max) * I_p / a_g^{p+1/2} + (2/c) * Li_{p - 1}(W_max) * M_p / a_g^p !!
# Compute Li_s(W) rigorously in 50-digit mpmath.iv for s in {0, 0.5, 1, 1.5, 2, 2.5}:
def polylog_iv(s, W):
    # Sum m = 1..250 + tail bound W^{251} / (251^s * (1 - W))
    M_cut = 250
    tot = sum((W ** m) / (iv.mpf(m) ** s) for m in range(1, M_cut + 1))
    tail = (W ** (M_cut + 1)) / ((iv.mpf(M_cut + 1) ** s) * (1 - W))
    return tot + iv.mpf([0, tail.b])

# Evaluate BOTH cases (q = 1 with g = c >= 2 and j >= 2c -> alpha_1 <= kap/(4*pi*S_all)!!
#                      AND q >= 2 with g >= 1 -> 1/(g*q^4) <= 1/16 and alpha_1 <= alpha_q2!)
def cert_E_phase_case(alpha_b, S_b, W_b, n_min_b, j_over_c_min, inv_gq4_max, inv_g2q5_max, inv_g2q4_max, label):
    C_sh_b = (iv.mpf(49) / 96) / (1 - (iv.mpf(49) / 96) * (alpha_b ** 2))
    gam_up_b = 1 / (iv.cos(alpha_b) ** 2)
    sinc_q_b = iv.sin(alpha_b / 4) / (alpha_b / 4)
    gam_lo_b = sinc_q_b * (iv.cos(alpha_b) ** 2)

    # Recall kap >= 2*pi*S_b * (j/c) * (c/j?) No: kap = 2*pi * sin(psi_j)*sin(alpha_1)/sin(theta_c)
    #   >= 2*pi * S_b * (j/c) * sin(alpha_1) >= 2*pi * S_b^2 * j_over_c_min * alpha_1!
    # So alpha_1 <= kap / (2*pi * S_b^2 * j_over_c_min)!
    # Wait: check! sin(psi_j)/sin(theta_c) >= S_b * (j/c) (since sin(psi_j) >= S_b*psi_j and sin(theta_c) <= theta_c!),
    #   and sin(alpha_1) >= S_b * alpha_1!
    #   SO kap = 2*pi * (sin(psi_j)/sin(theta_c)) * sin(alpha_1) >= 2*pi * S_b^2 * (j/c) * alpha_1!!
    #   Check for q = 1 (j/c >= 2): alpha_1 <= kap / (4*pi*S_b^2)!!
    #   Check for q >= 2 (j/c >= 1): alpha_1 <= kap / (2*pi*S_b^2)!!
    # Also (j/c)*alpha_1^2 = [(j/c)*alpha_1] * alpha_1 <= (kap / (2*pi*S_b^2)) * (kap / (2*pi*S_b^2*j_over_c_min))!
    # Therefore:
    #   B_1 = 2*pi*c*q*delta_shift <= C_sh_b * (c*q*kap)^2 / (2*pi * S_b^4 * j_over_c_min * g * q^2)
    #   B_3 <= (4/3)*gam_up_b * (c*q*kap)^2 / (2*pi * S_b^2 * j_over_c_min * sinc_q_b * cos(alpha_b) * g * q^2)
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

    # Factor 1/2 from (2 x / (1-x)^2) * sin^2(phi/2) <= (1/2) * (x / (1-x)^2) * phi^2!
    # Plus exact log-factor 1/(1 - Delta_max) where Delta_max = max_k (4 x_k/(1-x_k)^2) sin^2(phi_k/2) <= max_k (x_k/(1-x_k)^2) * phi_k^2!
    def pt_mom_b(p):
        Li_dsc = polylog_iv(iv.mpf(p - 1), W_b) if p >= 2 else W_b / (1 - W_b)
        return inv_g2q4_max * (C2 ** 4) * Li_dsc * M_vals[p] / (ag ** p)
    Delta_max = (K_B1 ** 2) * pt_mom_b(1) + 2 * K_B1 * K_B3 * pt_mom_b(2) + (K_B3 ** 2) * pt_mom_b(3)
    assert float(Delta_max.b) < 0.05
    log_fac = 1 / (1 - Delta_max)
    Ep_main = iv.mpf('0.5') * (
        (K_B1 ** 2) * mom_b(1) + 2 * K_B1 * K_B3 * mom_b(2) + (K_B3 ** 2) * mom_b(3)
    )
    # Tail |u_k| > U0_num / sqrt(tau_q):
    M_t = n_min_b / 2
    b_t = (4 / (PI ** 2)) * gam_lo_b
    def t_mom(p):
        base_d = (1 + b_t * (U0_num ** 2) / M_t) ** (-(M_t - p - 1))
        val_int = 2 * ((M_t / b_t) ** p) * (M_t / (2 * b_t * U0_num * (M_t - p - 1))) * base_d
        val_pt = 2 * ((M_t / b_t) ** p) * ((1 + b_t * (U0_num ** 2) / M_t) ** (-(M_t - p)))
        return inv_gq4_max * (C2 ** 4) * (val_int / (C2 ** (p + iv.mpf('0.5'))) + val_pt / (C2 ** p))
    Ep_tail = (W_b / (2 * ((1 - W_b) ** 2))) * (
        (K_B1 ** 2) * t_mom(1) + 2 * K_B1 * K_B3 * t_mom(2) + (K_B3 ** 2) * t_mom(3)
    )
    Ep = log_fac * (Ep_main + Ep_tail)
    print(f"  [{label}] K_B1<={float(K_B1.b):.4f}, K_B3<={float(K_B3.b):.4f}, Delta_max<={float(Delta_max.b):.5f} => E_phase <= {float(Ep.b):.6f} (tail={float(Ep_tail.b):.2e})")
    return Ep

print("\n[Step 5] Conjugate-pair quadratic phase error:")
# Case q = 1 (g = c >= 2, j >= 2c => j/c >= 2, 1/(g*q^4) <= 1/2, 1/(g^2*q^5) <= 1/4, 1/(g^2*q^4) <= 1/4, R_nonq == 0!):
E_phase_q1 = cert_E_phase_case(alpha_all, S_all, W_max_all, n_min_q1, iv.mpf(2), iv.mpf('0.5'), iv.mpf('0.25'), iv.mpf('0.25'), "q = 1")
# Case q >= 2 (g >= 1, j >= c+1 => j/c >= 1, 1/(g*q^4) <= 1/16, 1/(g^2*q^5) <= 1/32, 1/(g^2*q^4) <= 1/16):
E_phase_q2 = cert_E_phase_case(alpha_q2, S_q2, W_max_q2, n_min_q2, iv.mpf(1), iv.mpf(1)/16, iv.mpf(1)/32, iv.mpf(1)/16, "q >= 2")

# ------------------------------------------------------------------------------
# STEP 6: EXPLICIT CUBIC PERTURBATION E_cubic ON |w| = lam_* < 1 & FINAL MARGINS!
# ------------------------------------------------------------------------------
# Inside c*q*kap <= 0.90: c*kap^2 <= C2^2/2 = 0.405, n >= n_min >= 525, alpha_1 <= alpha_all <= 0.07169!
# Stretched-exponential bound on R_z/R_y inside c*kap^2 <= 0.405:
#   K_2C = (0.405 / (2*pi^4*S_all^4))^{1/3}, c_se_2C = 2*((alpha_all/2)*cot(alpha_all/2)) / K_2C >= 15.6!
#   Since 1 - W_0^{1/q} >= (-log W_0)*W_0/q >= (-log W_0)*W_0/n, summing over all n roots yields factor n^2:
#   For ALL n >= n_min: n^2 * (R_z/R_y)^{1 - 1/n} is strictly decreasing for n > (6 / c_se_eff)^3!
C1_2C = (C2 ** 2) / 2
K_2C = (C1_2C / (2 * (PI ** 4) * (S_all ** 4))) ** (iv.mpf(1) / 3)
xcot_2C = (alpha_all / 2) * iv.cos(alpha_all / 2) / iv.sin(alpha_all / 2)
c_se_2C = 2 * xcot_2C / K_2C
c_se_2C_eff = c_se_2C * ((n_min - 1) / n_min)
assert float(((6 / c_se_2C_eff) ** 3).b) < float(n_min.a)
RzRy_se_2C = iv.exp(-c_se_2C * (n_min ** (iv.mpf(1) / 3)))
E_cubic_2C = 20 * (n_min ** 2) * (RzRy_se_2C ** ((n_min - 1) / n_min)) / ((-iv.log(W_max_all)) * W_max_all)
assert float(E_cubic_2C.b) < 1e-40

err_q1 = iv.exp(E_phase_q1 + E_cubic_2C) - 1
err_q2 = iv.exp(R_nonq_bound + E_phase_q2 + E_cubic_2C) - 1
margin_q1 = iv.mpf(rat_main_cert) - err_q1
margin_q2 = iv.mpf(rat_main_cert) - err_q2

print("\n==========================================================================")
print(f"[FINAL LEMMA R2.3 CERTIFICATION (ALL c*q*kap <= 0.90, n >= 132, c >= 2)]:")
print(f"  Explicit cubic perturbation: c_se_2C >= {float(c_se_2C.a):.4f} => E_cubic <= {float(E_cubic_2C.b):.3e}")
print(f"  rat_main (3-line Polya-Szego, H_q >= 0.544) >= {float(rat_main_cert):.6f}")
print(f"  Case q = 1 (R_nonq = 0):  E_phase <= {float(E_phase_q1.b):.6f} => exp(E)-1 <= {float(err_q1.b):.6f} => MARGIN >= {float(margin_q1.a):+.6f}!!")
print(f"  Case q >= 2: |R_nonq| <= {float(R_nonq_bound.b):.6f}, E_phase <= {float(E_phase_q2.b):.6f} => exp(E)-1 <= {float(err_q2.b):.6f} => MARGIN >= {float(margin_q2.a):+.6f}!!")
print("==========================================================================")
assert float(margin_q1.a) > 0.25
assert float(margin_q2.a) > 0.10
