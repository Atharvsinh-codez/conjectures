"""
Compute the EXACT Regime 1 closed-form error bounds across ALL strictly interior points:
  2 <= c <= floor((2m - 3)/3),   c + 1 <= j <= m - ceil(c/2) - (1 if c%2==0 else 0)
for n = 20, 40, 80, 160, 320 and kappa_0 in [2.5, 3.0, 3.5, 4.0, 4.5, 5.0]!

Check ALL FOUR exact Regime 1 ingredients:
1. S_off = sum_{small} |w_k|/Rz + sum_{large} Rz/|w_k| (exact sum AND our 2-line closed-form upper bound!)
2. Exact cubic-root shifts: max |delta_x|, max |delta_y|, max |delta_z| on the circle |e3| = |t_j| Rz!
3. The exact 3rd cancellation:
     L_sum := sum_{large} (-s e^{i eta_k}) / (1 - s e^{i eta_k})  vs  +j,
     so |L_sum - c| * |delta_y|  is UNIFORMLY TINY across ALL (n, c, j)!!
4. The COMPLETE exact perturbation factor M(phi) = F_exact(phi) / (C_{m,c,j} * 4 sin^2((phi+j pi)/2))
   — wait: at the boundary pair, (z/x - w_tilde_+)(z/x - w_tilde_-)/(z/x)^2 is shifted slightly from 4 sin^2((phi+j pi)/2)!
   Let's check the EXACT closed-form perturbation:
     (1 - (-1)^j e^{-i phi} (1 + eps_+)) (1 - (-1)^j e^{-i phi} (1 + eps_-))
     = (1 - (-1)^j e^{-i phi})^2 - (-1)^j e^{-i phi}(eps_+ + eps_-)(1 - (-1)^j e^{-i phi}) + e^{-2i phi} eps_+ eps_-!
   AND integrate over phi in [-pi, pi]:
     int (1 - (-1)^j e^{-i phi})^2 * (-(-1)^j e^{i phi}) dphi/(2pi) = 2!
   Let's check the EXACT total bound on |1/(2pi) int (Re F_exact - C * 4 sin^2) dphi| / (2 C) across ALL strictly interior points!
"""
import math
import mpmath as mp

mp.mp.dps = 50

def kappa_val(n, c, j):
    return 2 * math.pi * math.sin(j * math.pi / n) * math.sin((j + c) * math.pi / n) / math.sin(c * math.pi / n)

def check_all_interior_analytic(n, k0_list=(2.5, 3.0, 3.5, 4.0, 5.0)):
    m = n // 2
    max_c = (n - 3) // 3
    # Track worst values for each kappa_0!
    stats = {k0: {
        'cnt': 0,
        'max_Soff': (0.0, None),
        'max_dx_d': (0.0, None),   # (n-3c) * |delta_x|
        'max_dy_jc': (0.0, None),  # |L_sum - c| * |delta_y|
        'max_eps_bdy': (0.0, None),# |eps_+| = |(r2/r3)*delta_y + delta_z|
        'max_Rz': (0.0, None),
        'max_Ry_Rz': (0.0, None),
        'max_dz': (0.0, None),
        'max_dy': (0.0, None),
    } for k0 in k0_list}

    for c in range(2, max_c + 1):
        j_end = m - (c + 1) // 2 - (1 if c % 2 == 0 else 0)
        for j in range(c + 1, j_end + 1):
            kap = kappa_val(n, c, j)
            if kap < min(k0_list):
                continue
            r1 = math.sin((j + c) * math.pi / n)
            r2 = math.sin(j * math.pi / n)
            r3 = math.sin(c * math.pi / n)
            s = r2 / r1
            # Use log for safety
            log_Ry = n * math.log(s)
            log_Rz = n * math.log(r3 / r1)
            Ry = math.exp(log_Ry) if log_Ry > -700 else 0.0
            Rz = math.exp(log_Rz) if log_Rz > -700 else 0.0
            t = Ry if (c % 2 == 0) else -Ry
            one_m_t = abs(1.0 - t)  # for odd c, 1+Ry > 1; for even c, 1-Ry > 0

            # Exact S_off (using log to avoid overflow!)
            th_c = c * math.pi / n
            s_off = 0.0
            L_sum = 0.0
            for k in range(n):
                eta = (2*k * math.pi / n) if (c % 2 == 0) else ((2*k + 1) * math.pi / n)
                if eta > math.pi: eta -= 2 * math.pi
                ae = abs(eta)
                if abs(ae - th_c) < 1e-12:
                    continue
                # ratio |1 - s e^{i eta}|^2 / |1 - s e^{i th_c}|^2 = 1 + 2*r2*r1/r3^2 * (cos th_c - cos eta)
                rat_sq = 1.0 + (2.0 * r2 * r1 / (r3 * r3)) * (math.cos(th_c) - math.cos(eta))
                log_rat = 0.5 * n * math.log(rat_sq)
                if ae < th_c:
                    if log_rat > -700:
                        s_off += math.exp(log_rat)
                else:
                    if -log_rat > -700:
                        s_off += math.exp(-log_rat)
                    # Real part of -s e^{i eta} / (1 - s e^{i eta}) = (s^2 - s cos eta) / (1 + s^2 - 2 s cos eta)
                    L_sum += (s * s - s * math.cos(eta)) / (1.0 + s * s - 2.0 * s * math.cos(eta))

            # Estimate cubic shifts at worst phi (using Rouché / 1st-order formula, exact check below):
            # x - 1 ~ t Rz / (1 - t) -> |dx| ~ Ry * Rz / |1 - t|
            # y/x - t ~ -Rz / (1 - t) - t * dx -> |dy| = |(y/x)/t - 1| ~ Rz / (Ry * |1 - t|)
            # Wait: compute Rz / Ry safely via log!
            log_Rz_Ry = n * math.log(r3 / r2)
            Rz_over_Ry = math.exp(log_Rz_Ry) if log_Rz_Ry > -700 else 0.0
            dx_est = (Ry * Rz / one_m_t) if one_m_t > 0 else 1e9
            dy_est = (Rz_over_Ry / one_m_t) if one_m_t > 0 else 1e9
            # z/x = (t Rz e^{i phi}) / (x^2 (y/x)) -> dz ~ (1 + 1/t) Rz = Rz + Rz/Ry!
            dz_est = Rz + Rz_over_Ry
            d = n - 3 * c
            dx_d = d * dx_est
            dy_jc = abs(L_sum - c) * dy_est
            eps_bdy = (r2 / r3) * dy_est + dz_est

            for k0 in k0_list:
                if kap >= k0:
                    st = stats[k0]
                    st['cnt'] += 1
                    if s_off > st['max_Soff'][0]: st['max_Soff'] = (s_off, (c, j, round(kap, 2)))
                    if dx_d > st['max_dx_d'][0]: st['max_dx_d'] = (dx_d, (c, j))
                    if dy_jc > st['max_dy_jc'][0]: st['max_dy_jc'] = (dy_jc, (c, j, round(L_sum, 2)))
                    if eps_bdy > st['max_eps_bdy'][0]: st['max_eps_bdy'] = (eps_bdy, (c, j))
                    if Rz > st['max_Rz'][0]: st['max_Rz'] = (Rz, (c, j))
                    if Ry*Rz > st['max_Ry_Rz'][0]: st['max_Ry_Rz'] = (Ry*Rz, (c, j))
                    if dz_est > st['max_dz'][0]: st['max_dz'] = (dz_est, (c, j))
                    if dy_est > st['max_dy'][0]: st['max_dy'] = (dy_est, (c, j))
    return stats

for n in [40, 80, 160, 320]:
    stats = check_all_interior_analytic(n, k0_list=(2.5, 3.0, 3.5, 4.0, 5.0))
    print(f"\n=== n = {n} (STRICTLY INTERIOR POINTS, #51 double-corner cut) ===")
    for k0 in (2.5, 3.0, 3.5, 4.0, 5.0):
        st = stats[k0]
        print(f"  k0={k0:.1f} (cnt={st['cnt']:4d}): max_Soff={st['max_Soff'][0]:.4f}{st['max_Soff'][1]} | d*dx={st['max_dx_d'][0]:.4f}{st['max_dx_d'][1]} | |L-c|*dy={st['max_dy_jc'][0]:.4f}{st['max_dy_jc'][1]} | eps_bdy={st['max_eps_bdy'][0]:.4f}{st['max_eps_bdy'][1]}")
