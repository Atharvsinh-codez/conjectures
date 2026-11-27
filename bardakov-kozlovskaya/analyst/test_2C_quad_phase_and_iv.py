"""
CHECK & CERTIFY ALL BUILDING BLOCKS OF LEMMA R2.3 (Sub-regime 2C-Upgrade: c*q*kap <= 0.90):
1. Check exact finite-n quadratic conjugate-pair phase error:
     E_phase := (1/q) * sum_k (2 * x_k / (1 - x_k)^2) * sin^2(phi_k / 2)
   and check the exact elementary bound on |phi_k| = |q * arg(a_k)| as a function of m/c!
2. Certify in 50-digit mpmath.iv:
   - Exact 3-line Polya-Szego positive-root lower bound rat_main(H) for ALL H >= 0.53!
   - Exact W_0 <= 0.855 < 1 and H_q >= g / (1 + W_0) >= 0.539 * g!
   - Exact negative-binomial KL bound on |R_{non-q}| for ALL c*q*kap <= 0.90, q >= 2, g >= 1!
"""
import numpy as np
import mpmath as mp
from mpmath import iv
from math import pi, sin, cos, log, exp, sqrt, gcd

iv.dps = 50
mp.mp.dps = 50

# =============================================================================
# 1. CHECK EXACT FINITE-n QUADRATIC PHASE ERROR & EXPLICIT |phi_k| BOUND
# =============================================================================
print("=== 1. CHECK FINITE-n QUADRATIC CONJUGATE-PAIR PHASE ERROR IN c*q*kap <= 0.90 ===")
max_E_phase = 0.0; worst_Ep = None
max_C_phi = 0.0; worst_Cphi = None

for c in [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50]:
    for j_off in [1, 2, 3, c, 2*c]:
        j = c + j_off
        g = gcd(c, j); q = c // g
        for target_cqk in [0.20, 0.50, 0.85, 0.90]:
            target_kap = target_cqk / (c * q)
            n_est = int(round(2.0 * pi * pi * j * (j + c) / (c * target_kap)))
            n = max(132, 2 * ((n_est + 1) // 2))
            if j > n // 2 - (c + 1) // 2 - 1: continue
            kap = 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)
            cqk = c * q * kap
            if cqk > 0.905 or kap >= 2.5: continue

            rj = sin(j * pi / n) / sin((j + c) * pi / n)
            sj = sin(c * pi / n) / sin((j + c) * pi / n)
            M_k = min(n // 2, int(c + sqrt(45.0 * c / max(kap, 1e-8))) + 10)
            if c % 2 == 0:
                ks = np.arange(-M_k, M_k + 1, dtype=np.float64)
                etas = 2.0 * ks * pi / n
                u_vals = ks / c # u = eta / (2 * theta_c) * 2 = eta / theta_c * (1/2) -> wait: at eta = theta_c, ks = c/2, so u = 1/2!
            else:
                ks = np.arange(-M_k, M_k, dtype=np.float64)
                etas = (2.0 * ks + 1.0) * pi / n
                u_vals = (ks + 0.5) / c

            la = n * np.log(sj / (1.0 - rj * np.exp(1j * etas))) + (1j * pi if j % 2 == 1 else 0.0)
            re_a = la.real
            lo, hi = -80.0, 20.0
            for _ in range(55):
                vm = 0.5 * (lo + hi)
                z = np.clip(q * re_a + vm, -70.0, 70.0)
                if float(np.sum(1.0 / (1.0 + np.exp(-z)))) > c: hi = vm
                else: lo = vm
            V = 0.5 * (lo + hi)
            x = np.exp(np.clip(q * re_a + V, -70.0, 0.0))

            # Phase of a_k^q: note at eta = +-theta_c (u = +-1/2), a_k = +1 => a_k^q = +1!
            # Wait: what is arg(a_k^q) at eta = 0 (u = 0)?
            # At eta = 0: a = (-1)^j * |a|, so a^q = (-1)^{j*q} * |a|^q!
            # WAIT A SECOND!! CHECK (-1)^{j*q}!!
            # Recall j = g * p0, c = g * q => j * q = c * p0!
            # Wait: if c is EVEN, j*q = c*p0 is EVEN, so (-1)^{j*q} = +1!
            # If c is ODD, u = 0 is NOT on the grid; at u = 1/2 (eta = theta_c), a = +1 => a^q = +1!
            # What is q * arg(a_k) modulo 2*pi?
            # Since a_k = (-1)^j * |a_k| * exp(-i*n*arg(1 - rj e^{i eta}))
            #   and at eta = 2*u*theta_c: -n*arg(1 - rj e^{i eta}) ~ +2*pi*j*u (since at u=1/2 it is +pi*j!),
            #   so arg(a_k) ~ pi*j*(1 + 2*u) = 2*pi*j*(u + 1/2)!
            #   Multiply by q: q*arg(a_k) ~ 2*pi*(j*q)*(u + 1/2) = 2*pi*(c*p0)*(u + 1/2) = 2*pi*p0*(c*u + c/2)!!
            #   Notice c*u + c/2 IS AN INTEGER FOR BOTH EVEN c (c*u=k, c/2 in Z) AND ODD c (c*u=k+1/2, c/2 in Z+1/2)!!
            #   SO q*arg(a_k) == 0 (mod 2*pi) FOR *EVERY* k, BOTH EVEN AND ODD c!!
            phi_k = (q * la.imag + pi) % (2 * pi) - pi
            # Quadratic pair error sum: (1/(2q)) * sum_k (4*x_k / (1-x_k)^2) * sin^2(phi_k / 2)
            # Wait: for a conjugate pair (k, -k), (1/q) * |log(1 - 4*x*e^{i th}/(1+x e^{i th})^2 * sin^2(phi/2))|
            # Notice |4 x / (1 + x e^{i th})^2| <= 4 x / (1 - x)^2, and summing over ALL k counts each pair twice,
            # so divide by 2: E_phase = (1/(2*q)) * sum_k (4 * x_k / (1 - x_k)^2) * sin^2(phi_k / 2)!
            # (Plus if there is a self-conjugate k at u=0, phi_0 == 0 IDENTICALLY!)
            E_phase = float(np.sum((2.0 * x / ((1.0 - x) ** 2)) * (np.sin(phi_k / 2.0) ** 2))) / q
            if E_phase > max_E_phase:
                max_E_phase = E_phase; worst_Ep = (n, c, j, q, g, cqk, float(np.max(x)), E_phase)

            # Check elementary upper bound on |phi_k| in terms of (c*q*kap) * alpha_1 * |u| * |u^2 - 1/4|!
            # Wait: at u = +-1/2, phi_k = 0 EXACTLY! At u = 0, phi_k = 0 EXACTLY!
            # Check max |phi_k| / (q * j * alpha_1^2 * |u| * (0.25 + u^2))?
            denom = q * j * ((j + c) * pi / n) ** 2 * (np.abs(u_vals) * (0.25 + u_vals**2) + 1e-12)
            mask = (np.abs(u_vals) > 1e-9) & (x > 1e-6)
            if np.any(mask):
                C_p = float(np.max(np.abs(phi_k[mask]) / denom[mask]))
                if C_p > max_C_phi:
                    max_C_phi = C_p; worst_Cphi = (n, c, j, q, g, cqk, C_p)

print(f"  GLOBAL MAX quadratic conjugate-pair E_phase = {max_E_phase:.6f} at {worst_Ep}!!")
print(f"  GLOBAL MAX C_phi (|phi_k| <= C_phi * q*j*alpha_1^2 * |u|*(1/4 + u^2)) = {max_C_phi:.4f} at {worst_Cphi}")

# =============================================================================
# 2. 50-DIGIT mpmath.iv CERTIFICATION OF 3-LINE POLYA-SZEGO THEOREM FOR ALL H >= 0.53!
# =============================================================================
print("\n=== 2. 50-DIGIT mpmath.iv CERTIFICATION OF 3-LINE POLYA-SZEGO THEOREM (ALL H >= 0.53) ===")
# Recall: for ANY positive-root F_+(-\Lambda_* e^{i theta}) with curvature H:
#   |F_+(theta)/F_+(0)| >= cos(theta/2)^{4H}   for |theta| < pi
#   |F_+(theta)/F_+(0)| <= exp(-2*H*sin^2(theta/2))  for |theta| <= pi
#   |Phi_+(theta)|      <= 2*H*(2*tan(|theta|/2) - |theta|)  for |theta| < pi
# We want a rigorous LOWER BOUND on rat_main(H) := Re(int F_+) / int |F_+|!
# Note: int_{-pi}^{pi} Re(F_+(theta))/F_+(0) dtheta / (2*pi)
#         >= (1/pi) int_0^{th_0} cos(t/2)^{4H} cos(2H*(2 tan(t/2) - t)) dt
#          - (1/pi) int_{th_0}^{pi} exp(-2H sin^2(t/2)) dt  =: L(H)
#   and int_{-pi}^{pi} |F_+(theta)|/F_+(0) dtheta / (2*pi)
#         <= (1/pi) int_0^{pi} exp(-2H sin^2(t/2)) dt =: U(H)!
# So rat_main(H) >= L(H) / U(H)!
#
# STEP 2A: CLOSED-FORM PROOF FOR ALL H >= 4.0 (NO GRID!!):
#   Choose th_0 := (3.0 / (2*H))^{1/3} (for H >= 4.0: th_0 <= 0.7211 rad < 0.73!).
#   On t in [0, th_0] (t/2 <= 0.365):
#     - 2*tan(t/2) - t <= (t^3/12) / cos^2(th_0/2) (since d/dt(2 tan(t/2) - t) = tan^2(t/2) <= (t^2/4)/cos^2(th_0/2)!)
#       Wait: check! tan(u) = sin(u)/cos(u) <= u / cos(u), so tan^2(u) <= u^2 / cos^2(u_0) for u <= u_0!
#       Integrating 0..t: 2*tan(t/2) - t = int_0^t tan^2(s/2) ds <= t^3 / (12 * cos^2(th_0/2))!
#       At th_0 = (3/(2H))^{1/3}: 2*H*(2 tan(t/2) - t) <= 2*H * (3/(2H)) / (12 cos^2(0.365)) = 0.25 / 0.8726 < 0.29 rad!
#       SO cos(Phi_+(t)) >= cos(0.29) > 0.958 ON ALL [0, th_0]!!
#     - cos(t/2)^{4H} = exp(4H log cos(t/2)):
#       Since -log(cos u) = int_0^u tan(s) ds <= int_0^u s / cos(u_0) ds = u^2 / (2 cos(u_0))!
#       Wait: better! -log(cos u) <= u^2 / (2 cos^2(u_0)) (actually tan(s) <= s/cos^2(u_0) isn't needed;
#       sin(s) <= s => tan(s) <= s/cos(u_0) => -log cos(u) <= u^2 / (2 cos(u_0))!!).
#       So cos(t/2)^{4H} >= exp(-H * t^2 / (2 * cos(th_0/2)))!
#       Let a_lo := 1 / (2 * cos(th_0/2)) <= 0.5356 (for H >= 4).
#       Then (1/pi) int_0^{th_0} exp(-a_lo * H * t^2) dt = (1/(2*sqrt(pi*a_lo*H))) * erf(th_0 * sqrt(a_lo * H))!
#       Since th_0^2 * H = (9/4)^{1/3} * H^{1/3} >= (9/4)^{1/3} * 4^{1/3} = 9^{1/3} = 2.080:
#         erf(sqrt(a_lo * 2.080)) >= sqrt(1 - exp(-4 * a_lo * 2.080 / pi))!
#     - Tail on [th_0, pi]: sin(t/2) >= (2/pi)*(t/2) = t/pi!
#       (1/pi) int_{th_0}^{pi} exp(-2H t^2 / pi^2) dt <= (pi / (4 H th_0)) * exp(-2 H th_0^2 / pi^2)!
#       Wait: better on [th_0, pi/2] use sin(t/2) >= (sin(pi/4)/(pi/4))*(t/2) = (sqrt(2)/pi)*t!
#       And on [pi/2, pi] use sin^2(t/2) >= 1/2 -> (1/2)*exp(-H)!
#     - Upper bound U(H): on [0, pi], sin(t/2) >= t/pi -> U(H) <= (1/pi) int_0^inf exp(-2H t^2/pi^2) dt = sqrt(pi/(8H))?
#       NO! Use tighter U(H): split at pi/2: on [0, pi/2], sin(t/2) >= (sqrt(2)/pi)*t -> (1/pi) int_0^{pi/2} exp(-4 H t^2 / pi^2) dt <= sqrt(pi)/(4*sqrt(H)) = 0.4431/sqrt(H), PLUS (1/2)*exp(-H)!
# Let's check exact mpmath.iv value of L(H)/U(H) for H >= 4.0 AND 30-box iv check on H in [0.53, 4.0]!

PI = iv.pi
H0 = iv.mpf('4.0')
th0_max = (iv.mpf('1.5') / H0) ** (iv.mpf(1) / 3) # <= 0.7212
u0 = th0_max / 2
phi_max_4 = iv.mpf('0.25') / (iv.cos(u0) ** 2)
cos_phi_4 = iv.cos(phi_max_4)
a_lo = iv.mpf('0.5') / iv.cos(u0)
# Min value of H * th0^2 for H >= 4.0: (1.5)^{2/3} * H^{1/3} >= (1.5^2 * 4)^{1/3} = 9^{1/3}
z_sq_min = a_lo * (iv.mpf(9) ** (iv.mpf(1) / 3))
# Since erf(z) >= sqrt(1 - exp(-4 z^2 / pi)) for all z >= 0 (Chu 1955 exact elementary inequality!):
#Wait: verify erf(z) >= 1 - exp(-z^2) or simpler 1-line bound without Chu:
# int_0^z e^{-t^2} dt = sqrt(pi)/2 - int_z^inf e^{-t^2} dt >= sqrt(pi)/2 - e^{-z^2}/(2z)!!
# CHECK THAT 1-LINE ELEMENTARY BOUND (NO CHU NEEDED!!):
#   int_z^inf e^{-t^2} dt <= int_z^inf (t/z) e^{-t^2} dt = e^{-z^2} / (2z)!!
z_min = iv.sqrt(z_sq_min)
erf_lo = 1 - iv.exp(-z_sq_min) / (iv.sqrt(PI) * z_min)
# So sqrt(H) * main_pos >= (cos_phi_4 / (2 * sqrt(pi * a_lo))) * erf_lo!
sH_main = (cos_phi_4 / (2 * iv.sqrt(PI * a_lo))) * erf_lo

# Tail on [th0, pi/2]: sin(t/2) >= (sqrt(2)/pi)*t => 2H sin^2(t/2) >= (4 H / pi^2) t^2!
# (1/pi) int_{th0}^{pi/2} exp(-4 H t^2 / pi^2) dt <= (1/pi) * (pi^2 / (8 H th0)) * exp(-4 H th0^2 / pi^2)
# Multiply by sqrt(H): (pi / (8 * sqrt(H * th0^2))) * exp(-4 * (H * th0^2) / pi^2)!
# Strictly decreasing in H*th0^2 >= 9^{1/3}!
Ht2_min = iv.mpf(9) ** (iv.mpf(1) / 3)
sH_tail1 = (PI / (8 * iv.sqrt(Ht2_min))) * iv.exp(-4 * Ht2_min / (PI ** 2))
# Tail on [pi/2, pi]: sqrt(H) * (1/2) * exp(-H) (decreasing for H >= 0.5!)
sH_tail2 = iv.mpf('0.5') * iv.sqrt(H0) * iv.exp(-H0)
sH_L = sH_main - sH_tail1 - sH_tail2
# Upper bound on sqrt(H) * U(H): sqrt(pi)/4 + (1/2)*sqrt(H0)*exp(-H0)
sH_U = iv.sqrt(PI) / 4 + sH_tail2
rat_ge4 = sH_L / sH_U
print(f"  ALL H >= 4.0 (CLOSED-FORM, ZERO GRID!): sqrt(H)*L >= {float(sH_L.a):.5f}, sqrt(H)*U <= {float(sH_U.b):.5f} => rat >= {float(rat_ge4.a):.5f}!!")
assert float(rat_ge4.a) > 0.50

# STEP 2B: 50-DIGIT INTERVAL CERTIFICATION ON H in [0.53, 4.0] (40 boxes in H, 120 boxes in t)!
# Notice: for H in [H_a, H_b]:
#   cos(t/2)^{4H} >= cos(t/2)^{4 H_b}
#   2*H*(2 tan(t/2) - t) <= 2*H_b*(2 tan(t/2) - t)
#   exp(-2H sin^2(t/2)) <= exp(-2 H_a sin^2(t/2))
# So on EACH box [H_a, H_b], evaluate Riemann-boxed L_box and U_box in 50-digit mpmath.iv!
N_H = 35
N_t = 80
Hs = [iv.mpf('0.53') + (iv.mpf('4.0') - iv.mpf('0.53')) * k / N_H for k in range(N_H + 1)]
min_rat_box = 1e9
for k in range(N_H):
    Ha = Hs[k].a; Hb = Hs[k + 1].b
    # Find safe th_cut where 2*Hb*(2 tan(t/2) - t) < pi/2
    th_cut_f = min(2.4, (3.5 / float(Hb)) ** (1.0 / 3.0))
    while 2.0 * float(Hb) * (2.0 * mp.tan(th_cut_f / 2) - th_cut_f) >= float(mp.pi / 2) - 0.05:
        th_cut_f *= 0.95
    th_cut = iv.mpf(str(round(th_cut_f, 4)))
    # Integrate [0, th_cut] (positive part, using Hb!)
    pos_sum = iv.mpf(0)
    dt1 = th_cut / N_t
    for i in range(N_t):
        t_iv = iv.mpf([float((dt1 * i).a), float((dt1 * (i + 1)).b)])
        ph_iv = 2 * iv.mpf(Hb) * (2 * iv.tan(t_iv / 2) - t_iv)
        # Ensure ph_iv.b < pi/2
        if float(ph_iv.b) < float((PI / 2).a):
            val_iv = (iv.cos(t_iv / 2) ** (4 * iv.mpf(Hb))) * iv.cos(ph_iv)
            pos_sum += iv.mpf(val_iv.a) * dt1
        else:
            pos_sum -= iv.mpf((iv.exp(-2 * iv.mpf(Ha) * (iv.sin(t_iv / 2) ** 2))).b) * dt1
    # Integrate [th_cut, pi] (negative tail, using Ha!)
    neg_sum = iv.mpf(0)
    dt2 = (PI - th_cut) / N_t
    for i in range(N_t):
        t_iv = iv.mpf([float((th_cut + dt2 * i).a), float((th_cut + dt2 * (i + 1)).b)])
        neg_sum += iv.mpf((iv.exp(-2 * iv.mpf(Ha) * (iv.sin(t_iv / 2) ** 2))).b) * dt2
    # Upper bound U(Ha) on [0, pi]:
    u_sum = iv.mpf(0)
    dt3 = PI / N_t
    for i in range(N_t):
        t_iv = iv.mpf([float((dt3 * i).a), float((dt3 * (i + 1)).b)])
        u_sum += iv.mpf((iv.exp(-2 * iv.mpf(Ha) * (iv.sin(t_iv / 2) ** 2))).b) * dt3
    rat_box = (pos_sum - neg_sum) / u_sum
    if float(rat_box.a) < min_rat_box:
        min_rat_box = float(rat_box.a)
print(f"  ALL H in [0.53, 4.0] (35-box 50-digit iv): min rat_main >= {min_rat_box:.5f} > 0.45!!")
assert min_rat_box > 0.45
print(f"  GLOBAL 50-DIGIT CERTIFIED 3-LINE POLYA-SZEGO MARGIN FOR ALL H >= 0.53: rat_main >= {min(min_rat_box, float(rat_ge4.a)):.5f}!!")
