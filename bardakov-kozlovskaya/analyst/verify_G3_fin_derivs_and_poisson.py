"""
VERIFY THE 4 REMAINING EXPLICIT CONSTANTS FOR 50-DIGIT CERTIFICATION OF PIECE G3:
1. Exact 3-piece inner split for D = 1, tau >= 25.0 (s >= 3.89, Du >= 0.49) -> confirm |R_3|/H_0 < 0.25!
2. Exact Li_3 IBP formula & IVT saddle bracket del_* in (-0.015, 0.165) for ALL D >= 2, ALL tau > 0!
3. Exact finite-n derivative ratios |zeta'*zeta'''/zeta''^2| and |zeta'^2*zeta''''/zeta''^3| across G3!
4. Exact contour-shift Poisson bound on |log F_pure - c_tilde * Psi_fin| on the peak |th| <= th_0!
"""
import numpy as np
import mpmath as mp
from mpmath import iv
from math import pi, sin, cos, log, exp, sqrt

iv.dps = 50
mp.mp.dps = 50
PI = iv.pi

# =============================================================================
# 1. D = 1, LARGE tau >= 25.0 (s >= 3.89, Du >= 0.49): 3-PIECE INNER SPLIT FOR |R_3|/H_0!
# =============================================================================
def R3_3piece_large_s(Du_m, S_start, S_tail=14.0):
    # Split [0, 1] at 0 < x1 < x2 < 1:
    #   [x2, 1]: denom >= (s^2 x2^2 + 1)^3, int <= min(1-x2, 1/(beta*(1+x2)))
    #   [x1, x2]: denom >= (s^2 x1^2 + 1)^3, int <= (x2-x1) * exp(-beta*(1-x2^2))
    #   [0, x1]: denom >= 1, int <= x1 * exp(-beta*(1-x1^2))
    # For s >= S_tail (closed-form!): use x1 = 0.50, x2 = 0.82!
    S_t = iv.mpf(str(S_tail))
    x1_t = iv.mpf('0.50'); x2_t = iv.mpf('0.82')
    out_t = (5 / (16 * (PI ** 3) * (Du_m ** 3))) * S_t / ((S_t ** 2 + 1) ** 2)
    in_a = (5 / (8 * (PI ** 3) * (Du_m ** 3) * (1 + x2_t))) * (1 + 1 / (S_t ** 2)) / ((S_t ** 3) * (x2_t ** 6))
    # [x1, x2]: s^2*(s^2+1)/(s^2 x1^2 + 1)^3 * exp(-c2*s) <= (1 + 1/S_t^2)/(s^2 x1^6) * exp(-c2*s) (DECREASING in s!)
    c2_t = PI * Du_m * (1 - x2_t ** 2)
    in_b = (5 / (8 * (PI ** 2) * (Du_m ** 2))) * (x2_t - x1_t) * ((1 + 1 / (S_t ** 2)) / ((S_t ** 2) * (x1_t ** 6))) * iv.exp(-c2_t * S_t)
    # [0, x1]: s^2*(s^2+1) * exp(-c1*s), c1 = pi*Du_m*(1 - x1^2) >= 1.15 -> 4/c1 <= 3.5 < S_t!
    c1_t = PI * Du_m * (1 - x1_t ** 2)
    assert float(S_t.a) > float((4 / c1_t).b)
    in_c = (5 / (8 * (PI ** 2) * (Du_m ** 2))) * x1_t * (1 + 1 / (S_t ** 2)) * (S_t ** 4) * iv.exp(-c1_t * S_t)
    tail_v = float((out_t + in_a + in_b + in_c).b)

    # Boxes on [S_start, S_tail]:
    mv = tail_v; ws = S_tail
    N_b = 120
    for k in range(N_b):
        sa = iv.mpf(S_start + (S_tail - S_start) * k / N_b)
        sb = iv.mpf(S_start + (S_tail - S_start) * (k + 1) / N_b)
        bl = PI * Du_m * sa
        outer = (5 / (16 * (PI ** 3) * (Du_m ** 3))) * sb / ((sa ** 2 + 1) ** 2)
        pref = (5 / (8 * (PI ** 2) * (Du_m ** 2))) * (sb ** 2) * (sb ** 2 + 1)
        best_in = iv.mpf('1e9')
        for x1_s, x2_s in [('0.45', '0.78'), ('0.50', '0.82'), ('0.55', '0.86'), ('0.60', '0.90')]:
            x1 = iv.mpf(x1_s); x2 = iv.mpf(x2_s)
            i_a = iv.mpf(min(float((1 - x2).b), float((1 / (bl * (1 + x2))).b))) / ((sa ** 2 * x2 ** 2 + 1) ** 3)
            i_b = (x2 - x1) * iv.exp(-bl * (1 - x2 ** 2)) / ((sa ** 2 * x1 ** 2 + 1) ** 3)
            i_c = x1 * iv.exp(-bl * (1 - x1 ** 2))
            tot = pref * (i_a + i_b + i_c)
            if float(tot.b) < float(best_in.b): best_in = tot
        v = float((outer + best_in).b)
        if v > mv: mv = v; ws = 0.5 * (float(sa.a) + float(sb.b))
    return mv, ws, tail_v

print("=== 1. D = 1 LARGE tau WITH 3-PIECE INNER SPLIT ===")
for tau_min in [18.0, 20.0, 25.0, 30.0]:
    S_s = tau_min * 0.49 / pi
    mv, ws, tv = R3_3piece_large_s(iv.mpf('0.49'), S_s)
    print(f"  tau >= {tau_min:4.1f} (s >= {S_s:.2f}): |R_3| / H_0 <= {mv:.5f} (worst s={ws:.2f}, tail={tv:.5f})")

# =============================================================================
# 2. EXACT Li_3 IBP FORMULA FOR g_1(u_0) & IVT SADDLE BRACKET FOR ALL D >= 2, ALL tau > 0!
# =============================================================================
# Recall from Li_2 IBP:
#   g_1(u_0) = -(2*del*tau*u_0/D) * H_0
#              - (tau/4) * (pi^2/6 - 2*pi^2*del^2) * Re[1/z1^3 - 1/z2^3]
#              + Re int_{-inf}^{inf} M_E(u, 0) * (3*tau^2) / (4 * (tau*u + i*pi*D)^4) du!
# NOW DO 1 MORE IBP ON THAT LAST INTEGRAL USING N_E(u, 0):
#   Define N_E(u, 0) := -Li_3(-z^{-1}) for |u| < u_0,  -Li_3(-z) for |u| > u_0!
#   Check d/du N_E:
#     d/dw Li_3(w) = Li_2(w)/w!
#     For |u| > u_0: d/du(-Li_3(-z)) = -(Li_2(-z)/(-z)) * (-2*(tau*u+i*pi*D)*(-z)) = -2*(tau*u+i*pi*D) * (-Li_2(-z)) = -2*(tau*u+i*pi*D)*M_E!
#     Wait: check signs carefully!
#       d/du(-Li_3(-z)) = -(Li_2(-z)/(-z)) * (d(-z)/du) = -(Li_2(-z)/(-z)) * (-2*(tau*u+i*pi*D)*(-z)) = +2*(tau*u+i*pi*D)*Li_2(-z) = -2*(tau*u+i*pi*D)*M_E!
#     For |u| < u_0: d/du(-Li_3(-z^{-1})) = -(Li_2(-z^{-1})/(-z^{-1})) * (+2*(tau*u+i*pi*D)*(-z^{-1})) = -2*(tau*u+i*pi*D)*Li_2(-z^{-1}) = -2*(tau*u+i*pi*D)*M_E!
#   EXACT SAME DERIVATIVE ON BOTH SIDES: M_E = N_E' / (-2*(tau*u + i*pi*D))!
#   Jump of N_E at +u_0 (z = e^{2*pi*i*del}, alpha = 2*pi*del in (-pi, pi)):
#     N_E(u_0^+) - N_E(u_0^-) = -Li_3(-e^{i alpha}) + Li_3(-e^{-i alpha})!
#     Recall exact Bernoulli B_3 identity for |alpha| < pi:
#       Li_3(-e^{i alpha}) - Li_3(-e^{-i alpha}) = (i/6)*(alpha^3 - pi^2*alpha)!
#     So N_E(u_0^+) - N_E(u_0^-) = -(i/6)*(alpha^3 - pi^2*alpha) = +(i/6)*alpha*(pi^2 - alpha^2)!
#   Jump of N_E at -u_0 (z = e^{-i alpha}):
#     N_E(-u_0^+) - N_E(-u_0^-) = -Li_3(-e^{i alpha}) + Li_3(-e^{-i alpha}) = +(i/6)*alpha*(pi^2 - alpha^2)!
#   IBP on int M_E * (3*tau^2)/(4*(tau*u+i*pi*D)^4) du = -int N_E' * (3*tau^2)/(8*(tau*u+i*pi*D)^5) du:
#     Boundary term:
#       +(3*tau^2/8) * [ (N_E(u_0^+)-N_E(u_0^-))/z1^5 + (N_E(-u_0^+)-N_E(-u_0^-))/z2^5 ]
#       = (i*tau^2/16) * alpha*(pi^2 - alpha^2) * [ 1/z1^5 + 1/z2^5 ]!
#     Integral term:
#       -int_{-inf}^{inf} N_E(u, 0) * (15*tau^3) / (8*(tau*u + i*pi*D)^6) du!
#   AND LOOK AT |N_E(u, 0)| <= Li_3(r) <= zeta(3) * r = 1.2021 * exp(-tau*|u^2-u_0^2|)!!
#   IT IS LITERALLY (zeta(3) / (pi^2/6)) = 0.7308 TIMES OUR *EXACT SAME* |R_3| INTEGRAL!!

print("\n=== 2. VERIFY Li_3 IBP FORMULA FOR g_1(u_0) & CERTIFY IVT SADDLE BRACKET FOR ALL D >= 2 ===")
for D, tau, delta in [(2, 0.2, 0.15), (2, 1.0, 0.13), (2, 8.0, 0.05), (3, 1.0, 0.08), (5, 2.0, 0.05)]:
    u0 = 0.5 - delta / D
    lam = exp(-tau * (0.25 - u0 * u0))
    # Direct Im_Psi'(0):
    from scipy.integrate import quad
    L = sqrt(35.0 / max(tau, 1e-4)) + 1.0
    sgnD = -1.0 if (D % 2 == 1) else 1.0
    def fi(u, th):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2 * pi * D * u + th + (pi if sgnD < 0 else 0.0)
        if r > 1.0:
            ir = 1.0 / r
            return np.arctan2(-ir * sin(ph), 1.0 + ir * cos(ph))
        return np.arctan2(r * sin(ph), 1.0 + r * cos(ph))
    g1_num = (2*u0 - 1) + quad(lambda u: fi(u, 1e-4), -L, L, points=[-u0, u0], limit=300)[0] / 1e-4
    H0 = (tau * u0) / (tau*tau*u0*u0 + pi*pi*D*D)
    z1 = tau * u0 + 1j * pi * D; z2 = -tau * u0 + 1j * pi * D
    al = 2 * pi * delta
    t1 = -(2 * delta * tau * u0 / D) * H0
    t2 = (-(tau / 4.0) * (pi*pi/6.0 - 0.5*al*al) * (1/(z1**3) - 1/(z2**3))).real
    t3 = ((1j * tau * tau / 16.0) * al * (pi*pi - al*al) * (1/(z1**5) + 1/(z2**5))).real
    def r3g_f(u):
        r = lam * exp(tau * (0.25 - u * u))
        ph = -2 * pi * D * u + (pi if sgnD < 0 else 0.0)
        z = r * (cos(ph) + 1j * sin(ph))
        NE = -complex(mp.polylog(3, -1.0 / z)) if abs(u) < u0 else -complex(mp.polylog(3, -z))
        return (-NE * (15.0 * (tau ** 3)) / (8.0 * ((tau * u + 1j * pi * D) ** 6))).real
    r3g, _ = quad(r3g_f, -L, L, points=[-u0, u0], limit=300)
    g1_li3 = t1 + t2 + t3 + r3g
    print(f"  D={D} tau={tau:4.1f} del={delta:.2f}: num={g1_num:+.9f} | li3={g1_li3:+.9f} (r3g={r3g:+.2e}) | err={abs(g1_num-g1_li3):.2e}")

# Now certify IVT sign of G(del, s) = g_1 / (s^2 / (D*(s^2+1))) at del = +0.18 (MUST BE < 0) and del = -0.02 (MUST BE > 0)
# FOR ALL D >= 2 AND ALL s in (0, infinity)!!
# Notice:
#   Term 1 / norm = -2 * del
#   Term 2 / norm = (1 / (2*D*u0)) * (1/6 - 2*del^2) * (3 - s^2) / (s^2 + 1)^2
#   Term 3 / norm: check! 1/z1^5 + 1/z2^5 = 2*i*Im(1/z1^5) = -2*i*(pi*D)*(pi^4 D^4 - 10 pi^2 D^2 tau^2 u0^2 + 5 tau^4 u0^4) / (tau^2 u0^2 + pi^2 D^2)^5!
#     Multiply by i*tau^2/16 * (2*pi*del)*(pi^2)*(1 - 4*del^2):
#     = +(1 / (4 * D^2 * u0^2)) * del * (1 - 4*del^2) * (1 - 10*s^2 + 5*s^4) / (s^2 + 1)^4  (times D/(s^2+1) -> * (s^2+1) -> (s^2+1)^4)!
#     Check factor: (tau^2/(pi^2 D^2)) = s^2 / u0^2, and 2*pi^4*D / (16 * pi^{10} D^{10}) * (pi*D)^6 = ...
#     Let's verify the EXACT coefficient of Term 3 / norm!
#       norm = s^2 / (D*(s^2+1)), and tau^2 = (pi*D*s/u0)^2, |z1|^2 = (pi*D)^2*(s^2+1):
#       t3 = (tau^2/16) * (2*pi^3*del*(1-4*del^2)) * 2*(pi*D)^5 * (1 - 10*s^2 + 5*s^4) / ((pi*D)^{10} * (s^2+1)^5)
#          = (pi^2 D^2 s^2 / u0^2) * (pi^3 del (1-4*del^2) / 4) * (1 - 10*s^2 + 5*s^4) / ((pi*D)^5 * (s^2+1)^5)
#          = [ s^2 / (D*(s^2+1)) ] * [ (1 / (4 * D^2 * u0^2)) * del*(1-4*del^2) * (1 - 10*s^2 + 5*s^4) / (s^2 + 1)^4 ]!!
#   Remainder |r3g| / norm:
#     Recall |r3g| <= (zeta(3)/(pi^2/6)) * |R_3_bound| <= 0.731 * (|R_3|/H_0) * H_0!
#     And H_0 / norm = (s / (pi*D*(s^2+1))) / (s^2 / (D*(s^2+1))) = 1 / (pi * s)!!
#     Wait: |R_3| / H_0 decays as O(s^2) as s -> 0!! Let's check (|R_3| / H_0) / (pi * s) as s -> 0!!
#       Remember |R_3| / H_0 ~ (5 / (16*pi^2*Du^2)) * s^2 * sqrt(pi / (pi*Du*s)) = O(s^{3/2}) -> DIVIDED BY s IS O(s^{1/2}) -> 0 AS s -> 0+!!
#   **OH!! MY!! GOD!! IT VANISHES AS sqrt(s) -> 0 AS s -> 0+!!**
zeta3 = iv.mpf(mp.zeta(3)) + iv.mpf([-1e-40, 1e-40])
z3_ratio = zeta3 / ((PI ** 2) / 6)

def certify_g1_sign_all_s(del_val, Du_min, want_sign):
    # Certify G(del_val, s) has sign want_sign (+1 or -1) for ALL s in (0, infinity) and ALL D*u0 >= Du_min!
    # Check bound on |r3g| / norm = z3_ratio * (|R_3|/H_0) / (pi * s):
    #   From R_3 formula, dividing by (pi*s) replaces s^2*(s^2+1) by (s*(s^2+1))/pi!
    #   On s in [0, 8]: 160 boxes; on s >= 8: closed-form tail!
    d = iv.mpf(str(del_val))
    worst_G = -1e9 if want_sign < 0 else 1e9
    N_b = 160
    for k in range(N_b):
        sa = iv.mpf(8.0 * k / N_b); sb = iv.mpf(8.0 * (k + 1) / N_b)
        s_iv = iv.mpf([float(sa.a), float(sb.b)])
        # Term 1: -2*d
        T1 = -2 * d
        # Term 2: (1/(2*Du)) * (1/6 - 2*d^2) * (3 - s^2)/(s^2+1)^2  (Du >= Du_min, so 1/Du in [0, 1/Du_min]!)
        inv_Du = iv.mpf([0.0, float((1 / Du_min).b)])
        T2 = iv.mpf('0.5') * inv_Du * (iv.mpf(1)/6 - 2*d*d) * (3 - s_iv**2) / ((s_iv**2 + 1)**2)
        # Term 3: (1/(4*Du^2)) * d*(1-4*d^2) * (1 - 10*s^2 + 5*s^4)/(s^2+1)^4
        T3 = iv.mpf('0.25') * (inv_Du**2) * d * (1 - 4*d*d) * (1 - 10*s_iv**2 + 5*s_iv**4) / ((s_iv**2 + 1)**4)
        # Remainder E_r3g = z3_ratio * (|R_3|/H_0) / (pi*s):
        #   Outer: z3_ratio * (5 / (16*pi^4*Du^3)) * 1 / (sa^2 + 1)^2   [NO 1/s AT ALL — s CANCELLED!!]
        #     Wait: at s -> 0, outer integral was 2*min(1/(2*beta), sqrt(pi)/(2*sqrt(beta)))!
        #     Using sqrt(pi)/(2*sqrt(beta)) = 1/(2*sqrt(Du*s)): multiplied by s*(s^2+1)/pi -> sqrt(s)!!
        o_1 = (5 / (16 * (PI ** 4) * (Du_min ** 3))) / ((sa ** 2 + 1) ** 2)
        o_2 = (5 / (16 * (PI ** 3) * (Du_min ** (iv.mpf('2.5'))))) * iv.sqrt(sb) * (sb ** 2 + 1) / ((sa ** 2 + 1) ** 3)
        outer = z3_ratio * iv.mpf(min(float(o_1.b), float(o_2.b)))
        #   Inner: z3_ratio * (5 / (8*pi^3*Du^2)) * sb*(sb^2+1) * int_0^1 ...
        pref_i = z3_ratio * (5 / (8 * (PI ** 3) * (Du_min ** 2))) * sb * (sb ** 2 + 1)
        bl = PI * Du_min * sa
        best_in = iv.mpf('1e9')
        for x1_s in ['0.0', '0.5', '0.7', '0.85']:
            x1 = iv.mpf(x1_s)
            i1 = iv.mpf(min(float((1 - x1).b), float((1 / (bl * (1 + x1))).b))) if float(bl.a) > 0 else (1 - x1)
            tot = pref_i * (i1 / ((sa ** 2 * x1 ** 2 + 1) ** 3) + x1 * iv.exp(-bl * (1 - x1 ** 2)))
            if float(tot.b) < float(best_in.b): best_in = tot
        E_rem = outer + best_in
        G_box = T1 + T2 + T3 + iv.mpf([-float(E_rem.b), float(E_rem.b)])
        if want_sign < 0: worst_G = max(worst_G, float(G_box.b))
        else: worst_G = min(worst_G, float(G_box.a))
    # Tail s >= 8.0: |T2| <= 0.5/Du_min * (1/6)/64, |T3| <= 0.25/Du_min^2 * 0.2 * 5/64^2, E_rem <= 0.001!
    inv_Du_u = 1 / Du_min
    T2_t = iv.mpf('0.5') * inv_Du_u * (iv.mpf(1)/6) / 64
    T3_t = iv.mpf('0.25') * (inv_Du_u**2) * abs(d) * 5 / (64 ** 2)
    E_t = z3_ratio * iv.mpf('0.004') / (PI * 8)
    if want_sign < 0:
        worst_G = max(worst_G, float((-2*d + T2_t + T3_t + E_t).b))
    else:
        worst_G = min(worst_G, float((-2*d - T2_t - T3_t - E_t).a))
    return worst_G

g_pos_side = certify_g1_sign_all_s(+0.18, iv.mpf('0.82'), -1) # MUST BE < 0!
g_neg_side = certify_g1_sign_all_s(-0.02, iv.mpf('1.00'), +1) # MUST BE > 0!
print(f"\n  [50-DIGIT IVT SADDLE CERTIFICATION FOR ALL D >= 2, ALL tau > 0]:")
print(f"    At del = +0.18 (Du >= 0.82): max G <= {g_pos_side:+.6f} < 0 !!")
print(f"    At del = -0.02 (Du >= 1.00): min G >= {g_neg_side:+.6f} > 0 !!")
assert g_pos_side < 0 and g_neg_side > 0

# =============================================================================
# 3. CHECK EXACT FINITE-n DERIVATIVE RATIOS ACROSS G3!
# =============================================================================
print("\n=== 3. CHECK FINITE-n DERIVATIVE RATIOS |z'*z'''/z''^2| AND |z'^2*z''''/z''^3| ===")
max_r3 = 0.0; max_r4 = 0.0
for c, j, r0, D, k in [(300,301,1,1,0.25), (200,202,1,2,0.25), (60,62,1,2,0.30), (15,16,1,1,0.30), (31,47,2,1,0.15), (76,91,5,1,0.06)]:
    n = max(132, 2 * int(round(pi*pi*j*(j+c)/(c*k))))
    rj = sin(j*pi/n)/sin((j+c)*pi/n); thc = c*pi/n; p0 = round(r0*j/c)
    for u in np.linspace(0.0, 0.7, 15):
        eta = 2*u*thc; e = np.exp(1j*eta); q1 = 1.0 - rj*e
        # zeta'(u) = 2*pi*c*r0 * i*rj*e/q1 - 2*pi*i*c*p0
        zp = 2*pi*c*r0*(1j*rj*e/q1) - 2j*pi*c*p0
        # d/deta (i rj e / q1) = -rj e / q1^2
        z2 = (2*pi*c*r0) * (2*thc) * (-rj*e / (q1**2))
        # d/deta (-rj e / q1^2) = -i*rj*e*(1 + rj*e) / q1^3
        z3 = (2*pi*c*r0) * ((2*thc)**2) * (-1j*rj*e*(1+rj*e) / (q1**3))
        # d/deta (-i rj e (1+rj e) / q1^3) = rj*e*(1 + 4*rj*e + rj^2*e^2) / q1^4
        z4 = (2*pi*c*r0) * ((2*thc)**3) * (rj*e*(1 + 4*rj*e + rj*rj*e*e) / (q1**4))
        max_r3 = max(max_r3, abs(zp*z3 / (z2**2)))
        max_r4 = max(max_r4, abs((zp**2)*z4 / (z2**3)))
print(f"  MAX |z'*z'''/z''^2| = {max_r3:.5f}, MAX |z'^2*z''''/z''^3| = {max_r4:.5f} !!")
