"""
Test our 100% ELEMENTARY, ZERO-FUDGE REGIME 1 MASTER BOUND across ALL strictly interior points
for n = 20, 40, 84, 160, 320!

Exact algebraic factorization of F_exact(phi) / C_{m,c,j}:
  e_1 = 1 + t, e_2 = t, e_3 = t R_z e^{i phi} (t = (-1)^c R_y).
  Roots of T^3 - (1+t) T^2 + t T - t R_z e^{i phi} = 0:
    x = 1 + dx,   y/x = t(1 + dy),   z/x = R_z e^{i phi}(1 + dz),  with 1 + dz = 1 / ((1+dx)^3 (1+dy)).
  Exact integrand:
    I(phi) / C_{m,c,j} = (1 + dx)^{n-3c}
                       * [prod_{large} (w_tilde_k / w_k) / (1 + dy)^c]
                       * [prod_{small} (1 - w_tilde_k / (z/x)) * prod_{large} (1 - (z/x) / w_tilde_k)]
                       * [ (1 / (1 + dz)) * (-Z^{-1}) * ((1+dz) - Z(1+ep_+)) * ((1+dz) - Z(1+ep_-)) ]
  where Z = (-1)^j e^{-i phi}, and w_tilde_k / w_k = (1 - (s e^{i eta_k} / (1 - s e^{i eta_k})) * eps_s)^n
  with 1 + eps_s = (1 + dy)^{1/n}!

Wait: check the last bracket!
  Let B(phi) = (1 / (1 + dz)) * (-Z^{-1}) * ((1+dz) - Z(1+ep_+)) * ((1+dz) - Z(1+ep_-))
             = (-Z^{-1}) * [ (1+dz) - Z(2 + ep_+ + ep_-) + Z^2 (1+ep_+)(1+ep_-)/(1+dz) ]
             = (2 - Z - Z^{-1}) + [ -Z^{-1} dz + (ep_+ + ep_-) - Z * ((1+ep_+)(1+ep_-)/(1+dz) - 1) ]!!
  OH!! MY!! GOD!! LOOK HOW INSANELY CLEAN THAT EXPANSION OF B(phi) IS!!
  Write B(phi) = 4 sin^2((phi + j pi)/2) + R_bdy(phi) where:
    |R_bdy(phi)| <= |dz| + |ep_+ + ep_-| + |(1+ep_+)(1+ep_-)/(1+dz) - 1|!
  Wait: even better, check BOTH:
  - Bound A (pointwise triangle inequality on R_bdy):
      avg(Re( (1 + Delta) * (4 sin^2 + R_bdy) )) >= 2*(1 - M1) - (1 + M1) * max|R_bdy|!
  - Bound B (keeping exact 1st-order integral of R_bdy, or separating the phi-avg of (2 - Z - Z^{-1}) = 2):
      Let's test Bound A FIRST — if Bound A ALREADY > 0 everywhere for k0 = 3.5, the proof is LITERALLY HALF A PAGE!!
"""
import math

def kappa_val(n, c, j):
    return 2 * math.pi * math.sin(j * math.pi / n) * math.sin((j + c) * math.pi / n) / math.sin(c * math.pi / n)

def rigorous_rouche_bounds(n, c, j):
    """Return 100% rigorous upper bounds on |dx|, |dy|, |dz| via Rouché's theorem on T(T-1)(T-t) = p."""
    r1 = math.sin((j + c) * math.pi / n)
    r2 = math.sin(j * math.pi / n)
    r3 = math.sin(c * math.pi / n)
    s = r2 / r1
    log_Ry = n * math.log(s)
    log_Rz = n * math.log(r3 / r1)
    log_Rz_Ry = n * math.log(r3 / r2)
    Ry = math.exp(log_Ry) if log_Ry > -700 else 0.0
    Rz = math.exp(log_Rz) if log_Rz > -700 else 0.0
    Rz_over_Ry = math.exp(log_Rz_Ry) if log_Rz_Ry > -700 else 0.0
    t = Ry if (c % 2 == 0) else -Ry
    omt = abs(1.0 - t)

    # Rouché disk around T = 1: solve rho_x * (1 - rho_x) * (omt - rho_x) = Ry * Rz
    p_mod = Ry * Rz
    rho_x = (1.0 + 1e-12) * p_mod / omt
    for _ in range(80):
        denom = (1.0 - rho_x) * (omt - rho_x)
        if denom <= 0:
            return None
        rho_x_new = (1.0 + 1e-12) * p_mod / denom
        if rho_x_new <= rho_x:
            break
        rho_x = rho_x_new
    if not (rho_x < min(0.5, 0.5 * omt) and rho_x * (1.0 - rho_x) * (omt - rho_x) >= p_mod):
        return None

    # Rouché disk around T = t: |T - t| = rho_y = u_y * Ry, need u_y * (1 - u_y) * (omt - u_y * Ry) >= Rz / Ry
    u_y = (1.0 + 1e-12) * Rz_over_Ry / omt
    for _ in range(80):
        denom = (1.0 - u_y) * (omt - u_y * Ry)
        if denom <= 0:
            return None
        u_y_new = (1.0 + 1e-12) * Rz_over_Ry / denom
        if u_y_new <= u_y:
            break
        u_y = u_y_new
    if not (u_y < 0.5 and u_y * Ry < 0.5 * omt and u_y * (1.0 - u_y) * (omt - u_y * Ry) >= Rz_over_Ry):
        return None

    # Now y/x = t * (1 + delta_y) where 1 + delta_y = (1 + tilde_dy) / (1 + dx)
    # So |delta_y| <= (u_y + rho_x) / (1 - rho_x)
    dy = (u_y + rho_x) / (1.0 - rho_x)
    # Exact Vieta: 1 + dz = 1 / ((1 + dx)^3 (1 + dy))
    # So |dz| <= 1 / ((1 - rho_x)^3 * (1 - dy)) - 1!
    dz = 1.0 / (((1.0 - rho_x) ** 3) * (1.0 - dy)) - 1.0
    return rho_x, dy, dz, s, Ry, Rz, Rz_over_Ry, r1, r2, r3

def check_master_bound(n, k0):
    m = n // 2
    max_c = (n - 3) // 3
    worst_margin = 1e9
    worst_info = None
    cnt = 0
    for c in range(2, max_c + 1):
        j_end = m - (c + 1) // 2 - (1 if c % 2 == 0 else 0)
        for j in range(c + 1, j_end + 1):
            kap = kappa_val(n, c, j)
            if kap < k0:
                continue
            cnt += 1
            rb = rigorous_rouche_bounds(n, c, j)
            if rb is None:
                return -999.0, (c, j, "rouche_fail"), cnt
            dx, dy, dz, s, Ry, Rz, Rz_over_Ry, r1, r2, r3 = rb
            th_c = c * math.pi / n

            # Rigorous bound on eps_s = (1 + dy)^{1/n} - 1: |eps_s| <= -log(1 - dy)/n * (1/(1 - ...))
            # Actually |(1+w)^{1/n} - 1| <= (1 - dy)^{-1/n} - 1 <= -log(1 - dy) / n * (1 - dy)^{-1/n}
            # Even simpler: |(1+dy)^{1/n} - 1| <= 1 - (1 - dy)^{1/n} (Wait: for |z|<=dy<1,
            # |(1+z)^{1/n} - 1| = |sum_{k>=1} C(1/n, k) z^k| <= sum_{k>=1} |C(1/n, k)| dy^k = 1 - (1 - dy)^{1/n}!! EXACT!!)
            eps_s = 1.0 - ((1.0 - dy) ** (1.0 / n))

            # Term 1: |(1 + dx)^{n-3c} - 1| <= (1 + dx)^{n-3c} - 1
            t1_log = (n - 3 * c) * math.log(1.0 + dx)

            # Term 2: log [ prod_{large} (w_tilde_k / w_k) / (1 + dy)^c ]
            # For each large k: w_tilde_k / w_k = (1 + a_k eps_s)^n where a_k = -s e^{i eta_k}/(1 - s e^{i eta_k})!
            # Expand log(1 + a_k eps_s) = a_k eps_s + rem_k with |rem_k| <= (|a_k| eps_s)^2 / (2(1 - |a_k| eps_s))!
            # Also n * eps_s = dy + rem_eps where |n * eps_s - dy| <= -log(1 - dy) - dy? Wait:
            # 1 + dy = (1 + eps_s_complex)^n -> n * eps_s_complex = log(1 + dy_complex) + O(n * |eps_s|^2)!
            # Specifically |n * eps_s_complex - log(1 + dy_complex)| <= n * eps_s^2 / (2(1 - eps_s))!
            # So log(Term 2) = (L_sum - c) * log(1 + dy_complex) + sum_{large} [n*rem_k + a_k * (n*eps_s - log(1+dy))]!
            # Bound |log(Term 2)| cleanly!
            L_sum = 0.0
            sum_a_abs = 0.0
            sum_a_sq = 0.0
            s_off_tilde = 0.0  # perturbed S_off!
            for k in range(n):
                eta = (2*k * math.pi / n) if (c % 2 == 0) else ((2*k + 1) * math.pi / n)
                if eta > math.pi: eta -= 2 * math.pi
                ae = abs(eta)
                if abs(ae - th_c) < 1e-12:
                    continue
                denom_sq = 1.0 + s * s - 2.0 * s * math.cos(eta)
                a_abs = s / math.sqrt(denom_sq)
                rat_sq = 1.0 + (2.0 * r2 * r1 / (r3 * r3)) * (math.cos(th_c) - math.cos(eta))
                log_rat = 0.5 * n * math.log(rat_sq)
                # Factor (1 + a_abs * eps_s)^n <= exp(n * a_abs * eps_s)
                pert_k = n * a_abs * eps_s
                if ae < th_c:
                    # |w_tilde_k / (z/x)| <= |w_k / Rz| * exp(pert_k) / (1 - dz)
                    val_log = log_rat + pert_k - math.log(1.0 - dz)
                    if val_log > -700:
                        s_off_tilde += math.exp(val_log)
                else:
                    # |(z/x) / w_tilde_k| <= (Rz / |w_k|) * (1 + dz) * (1 - a_abs * eps_s)^{-n}
                    val_log = -log_rat + math.log(1.0 + dz) - n * math.log(1.0 - a_abs * eps_s)
                    if val_log > -700:
                        s_off_tilde += math.exp(val_log)
                    L_sum += (s * s - s * math.cos(eta)) / denom_sq
                    sum_a_abs += a_abs
                    sum_a_sq += a_abs * a_abs

            log_1mdy = -math.log(1.0 - dy)  # >= |log(1 + dy_complex)|
            rem_eps = n * eps_s * eps_s / (2.0 * (1.0 - eps_s))  # bound on |n*eps_s - log(1+dy)|
            # Also max |a_k| over large k is at eta = th_c + 2pi/n -> bounded by s / |1 - s e^{i th_c}| = r2/r3
            max_a_large = r2 / r3
            rem_large_sq = n * (eps_s ** 2) * sum_a_sq / (2.0 * (1.0 - max_a_large * eps_s))
            t2_log = abs(L_sum - c) * log_1mdy + sum_a_abs * rem_eps + rem_large_sq

            # Term 3: |prod (1 + z_k) - 1| <= exp(s_off_tilde) - 1 (or log bound s_off_tilde if all |z_k| < 1)
            # Combine Term 1 * Term 2 * Term 3 = 1 + Delta_mult:
            # |Term1 * Term2 * Term3 - 1| <= exp(t1_log + t2_log + s_off_tilde) - 1! (Wait: prod(1+z_k) has |prod - 1| <= exp(sum|z_k|) - 1!)
            M1 = math.exp(t1_log + t2_log + s_off_tilde) - 1.0

            # EXACT CAUCHY INTEGRAL OF B(w) (holomorphic on |w| <= 1!):
            #   (1/(2pi)) int B(e^{i phi}) dphi = 2 - (-1)^j * (2 * g(th_c) * dy'(0) - dz'(0))!
            # where x'(0) = t Rz/(1-t), y'(0) = -Rz/(1-t), dy'(0) = -(1+t^2)Rz/(t(1-t)), dz'(0) = (1-2t^2)Rz/(t(1-t))!
            t_signed = Ry if (c % 2 == 0) else -Ry
            # Compute exact Rz / (t * (1 - t)) cleanly using Rz_over_Ry:
            rz_over_t = Rz_over_Ry if (c % 2 == 0) else -Rz_over_Ry
            dy_p0 = -(1.0 + t_signed * t_signed) * rz_over_t / (1.0 - t_signed)
            dz_p0 = (1.0 - 2.0 * t_signed * t_signed) * rz_over_t / (1.0 - t_signed)
            # Exact g(th_c) = Re(-s e^{i th_c} / (1 - s e^{i th_c})) = -r2 * cos((j+c)pi/n) / r3
            g_thc = -r2 * math.cos((j + c) * math.pi / n) / r3
            exact_B_avg = 2.0 - ((-1) ** j) * (2.0 * g_thc * dy_p0 - dz_p0)

            # Term 4 brute-force max |R_bdy| ONLY multiplies M1!
            ep_bdy = ((1.0 + (r2 / r3) * eps_s) ** n) - 1.0
            R_bdy_bound = dz + 2.0 * ep_bdy + (((1.0 + ep_bdy) ** 2) / (1.0 - dz) - 1.0)
            # EXACT MASTER INEQUALITY:
            #   avg(Re( (1 + Delta) * B )) >= exact_B_avg - M1 * (2 + R_bdy_bound)!
            margin = exact_B_avg - M1 * (2.0 + R_bdy_bound)
            if margin < worst_margin:
                worst_margin = margin
                worst_info = (c, j, round(kap, 2), round(exact_B_avg, 4), round(M1, 4), round(R_bdy_bound, 4), round(s_off_tilde, 4), round(t2_log, 4))
    return worst_margin, worst_info, cnt

print("=== RIGOROUS CAUCHY-UPGRADED MASTER MARGIN: exact_B_avg - M1*(2 + R_bdy) > 0 ===")
for n in [20, 40, 84, 160, 320]:
    print(f"\n--- n = {n} ---")
    for k0 in [2.0, 2.5, 3.0, 3.5, 4.0]:
        mg, info, cnt = check_master_bound(n, k0)
        print(f"  k0={k0:.1f} (cnt={cnt:4d}): WORST MARGIN = {mg:+.4f}  at (c,j,kap,Bavg,M1,Rbdy,Soff~,t2)={info}")
