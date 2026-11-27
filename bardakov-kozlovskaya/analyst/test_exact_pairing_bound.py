"""
TEST OUR 4-LINE EXACT ELEMENTARY PAIRING BOUND ON |Phi_exact(psi)| AT lambda = 1!!
Exact pairing bound (100% rigorous whenever 2j + c <= n/2, which holds for ALL kappa < 3.5):
  a_m * a_{-m} = (1 + zeta_m)^{-n},
  zeta_m = [4 * sin(j*pi/n) * sin((j+c)*pi/n) / sin^2(c*pi/n)] * exp(-i*(2j+c)*pi/n) * sin^2(m*pi/n).
Since Re(zeta_m) >= 0, |1 + zeta_m| >= 1, so for EVERY r >= 1:
  |a_m^r - a_{-m}^{-r}| / r <= n * |zeta_m| * |a_{-m}^{-1}|^r.
Summing over r >= 1:
  sum_{r>=1} |a_m^r - a_{-m}^{-r}| / r <= n * |zeta_m| * |a_{-m}^{-1}| / (1 - |a_{-m}^{-1}|).

Total bound on sup_psi |Phi_exact(psi)| <= sum_{r>=1} |B_r - A_r| / r:
  B_pair = 2 * sum_{m=1}^{floor((c-1)/2)} n * |zeta_m| * |a_{-m}^{-1}| / (1 - |a_{-m}^{-1}|)
         + (1_{c even}) * log(1 / (1 - |a_{-c/2}^{-1}|))
         + 2 * sum_{m = floor((c-1)/2)+1}^{floor((n-c)/2)} log(1 / (1 - |a_m|))   [Wait: careful with m=(n-c)/2 single term!]

Let's check:
1. Exact value of B_pair vs max|Phi_exact| across ALL Regime-2 points (kap < 3.5) at n = 80, 160, 320, 1000!
2. Exact value of B_pair on extreme corners at n = 10000 and n = 100000!
3. What is the EXACT threshold C_0 such that c * kap^2 >= C_0 AND kap < 3.5 => B_pair < pi/2 (with margin!)?!
   Also check: what about kap >= 3.5 (even if 2j+c > n/2, naive bound works)?
"""
import numpy as np
from math import gcd

def kappa_val(n, c, j):
    return 2.0 * np.pi * np.sin(j * np.pi / n) * np.sin((j + c) * np.pi / n) / np.sin(c * np.pi / n)

def compute_B_pair(n, c, j):
    rj = np.sin(j * np.pi / n) / np.sin((j + c) * np.pi / n)
    sj = np.sin(c * np.pi / n) / np.sin((j + c) * np.pi / n)
    th_star = np.pi * (n - c) / n
    # Coefficient in |zeta_m|:
    Z_coef = 4.0 * np.sin(j * np.pi / n) * np.sin((j + c) * np.pi / n) / (np.sin(c * np.pi / n) ** 2)

    m_half = (c - 1) // 2
    # 1. Paired terms m = 1 .. m_half (counted twice for symmetry m <-> -c-m)
    pair_sum = 0.0
    if m_half >= 1:
        ms = np.arange(1, m_half + 1, dtype=np.float64)
        # theta_{-m} = th_star + 2*m*pi/n
        th_ext = th_star + 2.0 * ms * np.pi / n
        # |a_{-m}^{-1}| = (|1 + rj e^{i th_ext}| / sj)^n
        mod_ext = np.exp(n * np.log(np.abs(1.0 + rj * np.exp(1j * th_ext)) / sj))
        zeta_abs = Z_coef * (np.sin(ms * np.pi / n) ** 2)
        # Note: exact bound using min(n*zeta_abs, 2.0) or n*zeta_abs / |1+zeta_m| (we use even simpler n*zeta_abs!)
        # Wait: if n*zeta_abs is big, we can ALSO bound |a_m^r - a_{-m}^{-r}|/r by log(1/(1-|a_m|)) + log(1/(1-|a_{-m}^{-1}|))!
        # Let's check BOTH pure n*zeta_abs and min(pair, individual)!
        term_pair = n * zeta_abs * mod_ext / (1.0 - mod_ext)
        pair_sum = 2.0 * float(np.sum(term_pair))

    # 2. Unpaired exterior midpoint m = c/2 when c is even:
    mid_ext = 0.0
    if c % 2 == 0:
        # theta = pi -> |1 - rj| / sj
        mod_mid = np.exp(n * np.log(abs(1.0 - rj) / sj))
        mid_ext = float(-np.log(1.0 - mod_mid))

    # 3. Unpaired interior terms: m = m_half + 1 .. (n - c) // 2
    # Carefully count multiplicity: total interior roots = n - c - 1.
    # For m = 1 .. floor((n-c-1)/2), there are TWO roots (theta_star - 2m*pi/n and -(theta_star - 2m*pi/n)).
    # If (n - c) is even, there is ONE extra root at theta = 0 (m = (n-c)/2).
    m_int_max = (n - c - 1) // 2
    tail_int = 0.0
    if m_half + 1 <= m_int_max:
        # Only need up to where mod_int > 1e-16
        m_hi = min(m_int_max, m_half + 2000)
        ms_i = np.arange(m_half + 1, m_hi + 1, dtype=np.float64)
        th_int = th_star - 2.0 * ms_i * np.pi / n
        lm_int = n * np.log(sj / np.abs(1.0 + rj * np.exp(1j * th_int)))
        keep = lm_int > -40.0
        if np.any(keep):
            mod_int = np.exp(lm_int[keep])
            tail_int = 2.0 * float(np.sum(-np.log(1.0 - mod_int)))
    if (n - c) % 2 == 0:
        m_zero = (n - c) // 2
        if m_zero > m_half:
            lm_0 = n * np.log(sj / (1.0 + rj))
            if lm_0 > -40.0:
                tail_int += float(-np.log(1.0 - np.exp(lm_0)))

    return pair_sum + mid_ext + tail_int, pair_sum, mid_ext, tail_int

print("=== TEST 4-LINE EXACT B_pair BOUND ACROSS n = 80, 160, 320, 1000 (kap < 3.5) ===")
for n in [80, 160, 320, 1000]:
    tot = 0; ok_all = 0
    max_B_ck2_12 = 0.0; worst_12 = None
    max_B_ck2_10 = 0.0; worst_10 = None
    max_ck2_fail = 0.0; worst_fail = None
    for c in range(2, int(0.11 * n) + 2):
        for j in range(c + 1, int(0.20 * n) + 2):
            if j > (n - c) // 2: continue
            kap = kappa_val(n, c, j)
            if kap >= 3.5: continue
            tot += 1
            ck2 = c * kap * kap
            B, ps, me, ti = compute_B_pair(n, c, j)
            if B < np.pi / 2:
                ok_all += 1
            else:
                if ck2 > max_ck2_fail:
                    max_ck2_fail = ck2; worst_fail = (c, j, kap, ck2, B, ps, me, ti)
            if ck2 >= 12.0 and B > max_B_ck2_12:
                max_B_ck2_12 = B; worst_12 = (c, j, kap, ck2, B, ps, me, ti)
            if ck2 >= 10.0 and B > max_B_ck2_10:
                max_B_ck2_10 = B; worst_10 = (c, j, kap, ck2, B, ps, me, ti)
    print(f"n={n:4d} ({tot:4d} R2 pts): B_pair<pi/2: {ok_all}/{tot} | MAX ck2 among fails={max_ck2_fail:.3f} {worst_fail[:5] if worst_fail else None}")
    print(f"   max B_pair for ck2>=10: {max_B_ck2_10:.4f} at {worst_10[:5] if worst_10 else None}")
    print(f"   max B_pair for ck2>=12: {max_B_ck2_12:.4f} at {worst_12[:5] if worst_12 else None}")

print("\n=== CHECK EXTREME LARGE-n POINTS WITH ck2 >= 12 (n = 10^5, 10^6) ===")
for n, c, j in [
    (100000, 100, 200), (100000, 200, 201), (100000, 500, 501),
    (100000, 50, 100), (100000, 20, 80), (100000, 10, 60),
    (1000000, 1000, 1001), (1000000, 500, 1000), (1000000, 200, 600),
]:
    kap = kappa_val(n, c, j)
    ck2 = c * kap * kap
    B, ps, me, ti = compute_B_pair(n, c, j)
    print(f"  n={n:7d} c={c:4d} j={j:4d}: kap={kap:.4f} ck2={ck2:7.2f} -> B_pair={B:.4f} (pair={ps:.4f}, mid={me:.4f}, tail={ti:.4f})")
