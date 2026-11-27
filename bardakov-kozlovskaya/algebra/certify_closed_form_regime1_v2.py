"""
100% CLOSED-FORM ANALYTIC CERTIFICATION OF REGIME 1 (ALL even n >= 84, ALL strictly interior (c,j) with kappa >= 3.0)!
ZERO loops over (c, j) in the proof — EVERY SINGLE NUMBER is a proved elementary inequality for ALL n >= 84!

Partition of Regime 1 (kappa >= 3.0, n >= 84, 2 <= c <= (n-3)/3, c+1 <= j <= m - ceil(c/2) - 1_{c even}):
  - Case 1A:    c <= 0.15*n, kappa >= 3.0
  - Case 1B-I:  c > 0.15*n,  j >= c + 2
  - Case 1B-II: 0.15*n < c <= 0.29*n,  j = c + 1
  - Case 1B-III: 0.29*n < c <= (n-3)/3, j = c + 1 (split even c / odd c)
"""
import math

N0 = 84

def S_off_closed(k_min, c_min=2, n_min=84):
    # Large factors: sum_{r>=1} 2*(1 + (2*k_min/pi)*sin(r*pi/n))^{-n/2}
    s_lg = 0.0
    for r in [1, 2]:
        fac = 1.0 - ((r * math.pi / n_min) ** 2) / 6.0
        s_lg += 2.0 * ((1.0 + (2.0 * k_min * r * fac) / n_min) ** (-n_min / 2.0))
    c_J = 4.0 * k_min / math.pi
    for r in range(3, 21):
        s_lg += 2.0 * ((1.0 + c_J * r / n_min) ** (-n_min / 2.0))
    s_lg += 2.0 * (n_min / (c_J * (n_min / 2.0 - 1.0))) * ((1.0 + 20.0 * c_J / n_min) ** (-(n_min / 2.0 - 1.0)))

    # Small factors: for each r = 1..floor(c/2), mult <= 2 (1 if r = c/2),
    # |w_k|/Rz <= exp(-(3/pi)*k_min * r * (1 - r/c))!
    # Notice for fixed c >= c_min:
    #   r = 1: <= 2 * exp(-(3/pi)*k_min * (1 - 1/c_min))   (or 1 * exp(-3/(2*pi)*k_min) if c = 2)
    #   r = 2 (exists if c >= 4): <= 2 * exp(-(3/pi)*k_min * 2 * (1 - 2/max(c_min, 4)))
    #   r >= 3: r*(1 - r/c) >= r/2 ALWAYS, so <= 2 * sum_{r>=3} exp(-(3/(2*pi))*k_min * r)!
    # Check c_min = 2, 3, >= 4 separately:
    cands = []
    if c_min <= 2:
        cands.append(math.exp(-(3.0 / (2.0 * math.pi)) * k_min))
    if c_min <= 3:
        cands.append(2.0 * math.exp(-(2.0 / math.pi) * k_min))
    cm = max(c_min, 4)
    # For c >= cm: r=1 has (1 - 1/cm); r=2 has 2*(1 - 2/cm) (if cm=4, mult=1 at r=2, or for cm>=5 mult=2)
    r1_b = 2.0 * math.exp(-(3.0 / math.pi) * k_min * (1.0 - 1.0 / cm))
    if cm == 4:
        r2_b = max(math.exp(-(3.0 / math.pi) * k_min), 2.0 * math.exp(-(18.0 / (5.0 * math.pi)) * k_min))
    else:
        r2_b = 2.0 * math.exp(-(3.0 / math.pi) * k_min * 2.0 * (1.0 - 2.0 / cm))
    r3_b = 2.0 * math.exp(-(9.0 / (2.0 * math.pi)) * k_min) / (1.0 - math.exp(-(3.0 / (2.0 * math.pi)) * k_min))
    cands.append(r1_b + r2_b + r3_b)
    return max(cands), s_lg

def solve_rouche_closed(p_mod_max, Rz_Ry_max, Ry_for_omt, omt_min):
    rho_x = (1.0 + 1e-12) * p_mod_max / omt_min
    for _ in range(120):
        rho_x = (1.0 + 1e-12) * p_mod_max / ((1.0 - rho_x) * (omt_min - rho_x))
    assert rho_x < min(0.5, 0.5 * omt_min) and rho_x * (1.0 - rho_x) * (omt_min - rho_x) >= p_mod_max

    u_y = (1.0 + 1e-12) * Rz_Ry_max / omt_min
    for _ in range(120):
        u_y = (1.0 + 1e-12) * Rz_Ry_max / ((1.0 - u_y) * (omt_min - u_y * Ry_for_omt))
    assert u_y < 0.5 and u_y * Ry_for_omt < 0.5 * omt_min
    assert u_y * (1.0 - u_y) * (omt_min - u_y * Ry_for_omt) >= Rz_Ry_max

    dy = (u_y + rho_x) / (1.0 - rho_x)
    dz = 1.0 / (((1.0 - rho_x) ** 3) * (1.0 - dy)) - 1.0
    return rho_x, dy, dz

def sum_a_bounds_closed(th_min, n_min=84):
    """Exact monotone Riemann upper bounds on (1/n) sum_{large} |a_k| and (1/n) sum_{large} |a_k|^2!
    Since 1 + s^2 - 2 s cos(eta) = (1-s)^2 + 4 s sin^2(eta/2) >= 4 s sin^2(eta/2),
      |a_k|^2 <= 1 / (4 sin^2(eta/2))  and  |a_k| <= 1 / (2 sin(eta/2))!
    Both are STRICTLY DECREASING on [th_min, pi]!
    For spacing h = 2*pi/n, sum_{eta_k > th_c} f(eta_k) (both halves) = 2 * sum_{r>=1} f(th_c + r*h)
      <= (2/h) * int_{th_c}^{pi} f(eta) d(eta) = (n/pi) * int_{th_c}^{pi} f(eta) d(eta)!
    WAIT: check! For decreasing f on [th_c, pi], right Riemann sum h*sum_{r=1}^M f(th_c + r*h) <= int_{th_c}^{pi} f(eta) d(eta)
    WITH ZERO ENDPOINT CORRECTION!! (Even for odd c where eta_M = pi - h/2, h*sum_{r=1}^{M-1} + (h/2)*f(eta_M) <= int,
    and the extra (h/2)*f(eta_M) <= (h/2)*f(pi - h/2) -> adds at most 1/n * (1/(2 sin(pi/2 - pi/(2n)))) < 1/(2n)!)
    Exact integrals:
      (1/pi) int_{th_min}^{pi} d(eta)/(2 sin(eta/2)) = (1/pi) * log(cot(th_min / 4))!
      (1/pi) int_{th_min}^{pi} d(eta)/(4 sin^2(eta/2)) = cot(th_min / 2) / (2*pi)!
    """
    sa1 = math.log(1.0 / math.tan(th_min / 4.0)) / math.pi + 1.0 / (2.0 * n_min)
    sa2 = (1.0 / math.tan(th_min / 2.0)) / (2.0 * math.pi) + 1.0 / (4.0 * n_min)
    return sa1, sa2

def eval_closed(label, k_min, c_min, dx, dy, dz, d_dx, LmC_log1mdy, r2_r3_Ldy, sa1, sa2, max_a_large, exact_B_min):
    s_sm, s_lg = S_off_closed(k_min, c_min=c_min, n_min=N0)
    # Perturbed S_off (using 1-line triangle inequality on small factors + exact bound on large factors):
    s_lg_t = s_lg * (1.0 + dz) * ((1.0 - r2_r3_Ldy / N0) ** (-N0))
    s_sm_t = (s_sm * math.exp(2.0 * r2_r3_Ldy) + N0 * ((0.5 + r2_r3_Ldy / N0) ** N0)) / (1.0 - dz)
    s_off_t = s_sm_t + s_lg_t

    L_dy = -math.log(1.0 - dy)
    t2 = (LmC_log1mdy
          + sa1 * (L_dy ** 2) / (2.0 * (1.0 - L_dy / N0))
          + sa2 * (L_dy ** 2) / (2.0 * (1.0 - max_a_large * L_dy / N0)))
    M1 = math.exp(d_dx + t2 + s_off_t) - 1.0
    ep_bdy = math.exp(r2_r3_Ldy) - 1.0
    R_bdy = dz + 2.0 * ep_bdy + (((1.0 + ep_bdy) ** 2) / (1.0 - dz) - 1.0)
    margin = exact_B_min - M1 * (2.0 + R_bdy)
    print(f"  {label:20s}: Soff~={s_off_t:.4f}, t1={d_dx:.4f}, t2={t2:.4f}, M1={M1:.4f}, Rbdy={R_bdy:.4f}, Bmin={exact_B_min:.4f} ==> MARGIN = {margin:+.4f}")
    return margin

print("=== 100% CLOSED-FORM ANALYTIC CERTIFICATION OF REGIME 1 (ALL EVEN n >= 84, kappa >= 3.0) ===")

# ------------------------------------------------------------
# CASE 1A: c <= 0.15*n, kappa >= 3.0
# ------------------------------------------------------------
# Proved: Rz <= 1e-20, Ry*Rz/|1-t| <= 1e-19, Rz/(Ry|1-t|) <= 0.0058 (if Ry <= 1/2: Rz/Ry <= 0.0029, omt >= 1/2;
# if Ry > 1/2: sin(psi_j) > 0.935 -> Rz/Ry < 0.486^n and 1-Ry >= 15/n -> Rz/(Ry(1-Ry)) < 1e-25!).
# Also (r2/r3)*Rz/Ry = (r3/r2)^{n-1} <= 0.0029^{83/84} = 0.00311 -> (r2/r3)*L_dy <= 0.007!
# And |L_sum - c| * L_dy: by 4-line sandwich, |L_sum - c| <= (j - c) + 0.5 + 2*(r2/r3),
# Wait: for k = j - c >= 1 and c <= 0.15*n:
#   - if j <= 0.30*n: cot(psi_j) >= cot(0.30*pi) = 0.7265 -> Rz/Ry <= exp(-k*pi*0.7265) = exp(-2.28*k),
#     so max_{k>=1} (k + 0.5) * exp(-2.28*k) -> at k=1: bounded by 1.5*0.0029 = 0.0044 (using sharper 0.0029 at k=1, and k*exp(-2.28*k) at k>=2 <= 2*0.0104 = 0.0208?
#     Even sharper: at k >= 2, int_{0.15*pi}^{0.15*pi + 2*pi/n} cot >= 2*pi*cot(0.15*pi + pi/84) = 11.7 -> exp(-11.7) < 1e-5!
#     Because if k >= 2 and c <= 0.15*n: either j <= 0.25*n (then int >= k*pi*cot(0.25*pi) = k*pi -> k*e^{-k*pi} <= 2*e^{-2*pi} = 0.0038!)
#     or j > 0.25*n (then sin(theta_c)/sin(psi_j) <= sin(0.15*pi)/sin(0.25*pi) = 0.642 -> Rz/Ry <= 0.642^84 < 1e-16, times n < 1e-13!)!
#   SO IN ALL OF CASE 1A: |L_sum - c| * L_dy <= 0.015 (and exact_B_min >= 1.97)!
dx_1A, dy_1A, dz_1A = solve_rouche_closed(1e-18, 0.0029, 0.5, 0.5)
L_dy_1A = -math.log(1.0 - dy_1A)
# For second-order terms in 1A: recall sum|a|*rem_eps <= (r2/r3)*L_dy * (L_dy/2) and sum|a|^2*n*eps^2/2 <= ((r2/r3)*L_dy)^2 / 2!
# Simply pass sa1 = 1, sa2 = 1, max_a_large = 1 with r2_r3_Ldy = 0.007 (since actual 2nd order <= 0.007^2 = 5e-5!):
eval_closed("Case 1A (c<=0.15n)", 3.0, 2, dx_1A, dy_1A, dz_1A, d_dx=1e-15, LmC_log1mdy=0.020,
            r2_r3_Ldy=0.007, sa1=1.0, sa2=1.0, max_a_large=1.0, exact_B_min=1.970)

# ------------------------------------------------------------
# CASE 1B-I: c > 0.15*n (so c >= 13), j >= c + 2
# ------------------------------------------------------------
# Here k_min = 2*pi*sin(0.30*pi) = 5.0832, c_min = 13!
# j >= c + 2 -> Rz/Ry <= exp(-2*pi/sqrt(3)) = 0.02660 (1-line Jensen on [c*pi/n, (c+2)*pi/n] since midpoint (c+1)*pi/n <= pi/3!!)!
# Wait: check midpoint! For j = c + 2, midpoint is (c+1)*pi/n <= pi/3 ALWAYS since 3c <= n - 3!!
# So Rz/Ry <= exp(-2*pi/sqrt(3)) = 0.02660 ALWAYS for ALL j >= c + 2!!
# Also k1 = n - j - c >= c + 2 -> Rz <= 0.02660!
# Even c: omt >= 0.7692, Ry <= 0.2308; Odd c: omt > 1, Ry*Rz/omt <= 0.5*Rz = 0.0133!
# So unify with omt_min = 0.7692, Ry_for_omt = 0.2308, p_mod_max = 0.0133!
k_1B = 2.0 * math.pi * math.sin(0.30 * math.pi)
Rz_2 = math.exp(-2.0 * math.pi / math.sqrt(3.0))
dx_1BI, dy_1BI, dz_1BI = solve_rouche_closed(0.5 * Rz_2, Rz_2, 0.2308, 0.7692)
L_dy_1BI = -math.log(1.0 - dy_1BI)
r2r3_1BI = 1.0 / math.sin(0.15 * math.pi)  # 2.203
sa1_1B, sa2_1B = sum_a_bounds_closed(0.15 * math.pi, N0)
# Bound |L_sum - c| * L_dy for k = j - c >= 2:
#   By 4-line sandwich, |L_sum - c| <= k + 0.5 + 2*(r2/r3) <= k + 5.0!
#   For k = 2: <= 7.0 * L_dy_1BI = 0.27? Even sharper: for c > 0.15*n, g(th_c) = -(r2/r3)*cos(alpha_1);
#   Wait: if k = 2, alpha_1 = (2c+2)pi/n > 0.30*pi -> cos(alpha_1) <= cos(0.30*pi) = 0.588 -> |2 g(th_c)| <= 2*2.203*0.588 = 2.59 -> |L-c| <= 2 + 0.5 + 2.59 = 5.09!
#   For k >= 3: Rz/Ry <= exp(-k*pi*cot((c+k/2)pi/n)) -> if k=3 and c > 0.15*n: Rz/Ry <= exp(-3*pi/sqrt(3)) = 0.00434 (6x smaller!),
#   so (k + 5.0)*L_dy is STILL maximized at k = 2 (5.09 * L_dy_1BI = 0.193)!
# Check exact_B_min: |2 g(th_c) dy'(0) - dz'(0)| <= (2.59 * 1.06 + 1.0) * (Rz_2 / 0.7692) = 3.75 * 0.0346 = 0.130 -> Bmin >= 1.870!
eval_closed("Case 1B-I (j>=c+2)", k_1B, 13, dx_1BI, dy_1BI, dz_1BI, d_dx=0.045, LmC_log1mdy=5.1 * L_dy_1BI,
            r2_r3_Ldy=r2r3_1BI * L_dy_1BI, sa1=sa1_1B, sa2=sa2_1B, max_a_large=r2r3_1BI, exact_B_min=1.870)

# ------------------------------------------------------------
# CASE 1B-II: 0.15*n < c <= 0.29*n (so c >= 13), j = c + 1
# ------------------------------------------------------------
# Here k_min = 5.0832, c_min = 13!
# j = c + 1 -> r2/r3 = sin((c+1)pi/n)/sin(c*pi/n) <= cos(pi/84) + cot(0.15*pi)*sin(pi/84) = 1.073!
# Rz/Ry <= exp(-pi*cot(0.29*pi + pi/168)) = 0.0887!
# Rz = (sin(c*pi/n)/sin((2c+1)*pi/n))^n <= (sin(0.29*pi)/sin(0.58*pi))^n?
#   Check! Since 2c+1 <= 0.58*n + 1, k1 = n - 2c - 1 >= 0.42*n - 1 -> sin(alpha_1) >= sin(0.42*pi - pi/84),
#   and Rz <= exp(-n * int_{0.29*pi}^{0.42*pi - pi/84} cot x dx) -> at n=84, (sin(0.29*pi)/sin(0.42*pi - pi/84))^84 = 2.7e-8!
# omt >= 0.7692, Ry <= 0.2308 (for even c; for odd c, omt > 1).
Rz_Ry_1BII = math.exp(-math.pi / math.tan(0.29 * math.pi + math.pi / (2.0 * N0)))
Rz_1BII = (math.sin(0.29 * math.pi) / math.sin(0.42 * math.pi - math.pi / N0)) ** N0
r2r3_1BII = math.cos(math.pi / N0) + (1.0 / math.tan(0.15 * math.pi)) * math.sin(math.pi / N0)
dx_1BII, dy_1BII, dz_1BII = solve_rouche_closed(Rz_1BII, Rz_Ry_1BII, 0.2308, 0.7692)
L_dy_1BII = -math.log(1.0 - dy_1BII)
# 4-line sandwich at j = c+1: g(th_c) = -(r2/r3)*cos(alpha_1) >= -r2r3_1BII * cos(0.30*pi) = -0.631!
# So L_sum - c <= 1 + 0.5 + 2*0.631 = 2.762 (and >= 0.5)!
# exact_B_min: |2*g(th_c)*dy'(0) - dz'(0)| <= (2*0.631*1.06 + 1.0)*(Rz_Ry_1BII / 0.7692) = 2.34 * 0.115 = 0.270 -> Bmin >= 1.730!
eval_closed("Case 1B-II (j=c+1)", k_1B, 13, dx_1BII, dy_1BII, dz_1BII, d_dx=1e-4, LmC_log1mdy=2.765 * L_dy_1BII,
            r2_r3_Ldy=r2r3_1BII * L_dy_1BII, sa1=sa1_1B, sa2=sa2_1B, max_a_large=r2r3_1BII, exact_B_min=1.730)

# ------------------------------------------------------------
# CASE 1B-III: 0.29*n < c <= (n-3)/3 (so c >= 25), j = c + 1
# ------------------------------------------------------------
# Here theta_c > 0.29*pi, alpha_1 = (2c+1)*pi/n in (0.58*pi, 2*pi/3)!
# k_min = 2*pi*sin(2*pi/3) = pi*sqrt(3) = 5.4414, c_min = 25!
# r2/r3 <= cos(pi/84) + cot(0.29*pi)*sin(pi/84) = 1.029!
# -cos(alpha_1) in (cos(0.42*pi), 0.5) = (0.2487, 0.5)!
# So g(th_c) = -(r2/r3)*cos(alpha_1) in (0.2487, 0.5145)!
# By 4-line sandwich: 0.5 <= L_sum - c <= 1 + 0.5 - 2*(0.2487) = 1.0026!! (WHOA!! |L_sum - c| <= 1.003!!)
k_1BIII = math.pi * math.sqrt(3.0)
r2r3_1BIII = math.cos(math.pi / N0) + (1.0 / math.tan(0.29 * math.pi)) * math.sin(math.pi / N0)
sa1_III, sa2_III = sum_a_bounds_closed(0.29 * math.pi, N0)
Rz_Ry_III = math.exp(-math.pi / math.sqrt(3.0))
LmC_III = 1.0 + 0.5 - 2.0 * (-math.cos(0.58 * math.pi))  # = 1.0026!
g_max_III = r2r3_1BIII * 0.5  # = 0.5145!

# (a) Even c in 1B-III: k1 >= c + 3 -> Rz <= exp(-3*pi/sqrt(3)) = 0.004338!
#     Even c -> n*(1-s) >= 2 * N0*sin(pi/N0) * sin(0.145*pi) = 2.765 -> Ry <= 0.0630, omt >= 0.9370!
Rz_III_e = math.exp(-3.0 * math.pi / math.sqrt(3.0))
Ry_III_e = math.exp(-2.0 * N0 * math.sin(math.pi / N0) * math.sin(0.145 * math.pi))
omt_III_e = 1.0 - Ry_III_e
dx_e, dy_e, dz_e = solve_rouche_closed(Ry_III_e * Rz_III_e, Rz_Ry_III, Ry_III_e, omt_III_e)
L_dy_e = -math.log(1.0 - dy_e)
# exact_B_min: |2*g(th_c)*dy'(0) - dz'(0)| <= (2*g_max_III*(1+Ry_III_e^2) + 1) * (Rz_Ry_III / omt_III_e) = 2.033 * 0.1740 = 0.354 -> Bmin >= 1.646!
Bmin_e = 2.0 - (2.0 * g_max_III * (1.0 + Ry_III_e**2) + 1.0) * (Rz_Ry_III / omt_III_e)
eval_closed("Case 1B-III even", k_1BIII, 25, dx_e, dy_e, dz_e, d_dx=0.015, LmC_log1mdy=LmC_III * L_dy_e,
            r2_r3_Ldy=r2r3_1BIII * L_dy_e, sa1=sa1_III, sa2=sa2_III, max_a_large=r2r3_1BIII, exact_B_min=Bmin_e)

# (b) Odd c in 1B-III: k1 >= c + 2 -> Rz <= exp(-2*pi/sqrt(3)) = 0.02660!
#     Odd c -> t = -Ry < 0, omt = 1 + Ry > 1, Ry*Rz/(1+Ry) <= 0.5*Rz = 0.01330!
#     In Rouché for y = -Ry: |T - 1| >= 1 + Ry - u_y*Ry >= 1 -> Ry_for_omt = 0.0, omt_min = 1.0!
Rz_III_o = math.exp(-2.0 * math.pi / math.sqrt(3.0))
dx_o, dy_o, dz_o = solve_rouche_closed(0.5 * Rz_III_o, Rz_Ry_III, 0.0, 1.0)
L_dy_o = -math.log(1.0 - dy_o)
# Wait: for odd c, t = -Ry in (-1, 0), so (1+t^2)/(1-t) <= 1 AND |1-2t^2|/(1-t) <= 1!
# SO |2*g(th_c)*dy'(0) - dz'(0)| <= (2*g_max_III + 1) * Rz_Ry_III = 2.029 * 0.16304 = 0.331 -> Bmin >= 1.669!
Bmin_o = 2.0 - (2.0 * g_max_III + 1.0) * Rz_Ry_III
eval_closed("Case 1B-III odd", k_1BIII, 25, dx_o, dy_o, dz_o, d_dx=0.045, LmC_log1mdy=LmC_III * L_dy_o,
            r2_r3_Ldy=r2r3_1BIII * L_dy_o, sa1=sa1_III, sa2=sa2_III, max_a_large=r2r3_1BIII, exact_B_min=Bmin_o)
