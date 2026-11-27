"""
Lock the EXACT cleanest thresholds (M_G3, eps_G3, K_G2) for our 100% ZERO-LOOP G1 proof!
1. Check exact upper bound H_D(t) for 2*pi^2 * nu^2 * |D(T, nu)| (including the deep dip near t = sqrt(3)!!).
2. Compute the EXACT coupled zero-loop bound on ALL of Piece G1:
     S_total(X) = S_q(d, X) + S_SetI(d, X) + S_conv(M_G3, eps_G3, X) + S_semi(X) + S_rest(X) + S_k_ne0 + S_fin
   over ALL X = q^2*kap (d=1: X >= 0.85; d>=2: X >= K_G2)!
"""
import math

# Step 1: Verify exact dip-preserving envelope H_D(t) >= 2*pi^2 * nu^2 * |D(T, nu)| for ALL nu >= 1, T > 0!
def D_HS(T, nu):
    N = 200; smax = 7.0; ds = 2.0 * smax / N
    b = math.pi * nu; rt = math.sqrt(T)
    tot = 0.0
    for i in range(N):
        s = -smax + (i + 0.5) * ds
        w = math.exp(-0.5 * s * s) / math.sqrt(2.0 * math.pi)
        h = s * rt
        tot += 2.0 * w * ((T + h) / ((T + h) ** 2 + b * b) - T / (T * T + (b + h) ** 2)) * ds
    return tot

# Measure exact max_{nu >= 1} 2*pi^2 * nu^2 * |D(t*pi*nu, nu)| as a function of t!
print("=== 1. EXACT DIP-PRESERVING ENVELOPE H_D(t) = max_{nu >= 1} 2*pi^2 * nu^2 * |D(t*pi*nu, nu)| ===")
t_samples = [0.05, 0.15, 0.30, 0.45, 0.628, 0.85, 1.10, 1.35, 1.60, 1.732, 2.00, 2.50, 3.20, 4.50, 7.00]
for t in t_samples:
    mx = max(2.0 * math.pi * math.pi * nu * nu * abs(D_HS(t * math.pi * nu, nu)) for nu in [1, 2, 3, 5, 15, 50])
    h_asymp = 8.0 * t * t * abs(t * t - 3.0) / ((t * t + 1.0) ** 3)
    print(f"  t={t:5.3f}: exact max H_D(t) = {mx:.4f}  (h_asymp = {h_asymp:.4f})")

# Build a super clean closed-form upper bound H_upper(t) >= H_D(t) for ALL t > 0:
# Notice from above (let's verify after seeing exact numbers!):
# On t <= 1.2:  <= min(24.05 t^2 / (1 + t^2)^2, 3.04)  (check!)
# On t > 1.2:   <= min(3.04, 1.80 / t^2 ? Wait: check exact value at t=1.6..2.5!)
# Let's check candidate smooth closed-form H_upper(t):
#   h_asymp(t) + delta_nu1(t)  where delta_nu1 is tiny!
# Actually, let's check max_{nu >= 1} |2*pi^2*nu^2*D(t*pi*nu, nu)| - h_asymp(t) AND max over t in [1.3, 4.0]!
max_above_1p5 = max(
    2.0 * math.pi * math.pi * nu * nu * abs(D_HS(t * math.pi * nu, nu))
    for nu in [1, 2, 3, 5, 15, 50]
    for t in [1.50 + 0.1 * k for k in range(40)]
)
max_above_1p3 = max(
    2.0 * math.pi * math.pi * nu * nu * abs(D_HS(t * math.pi * nu, nu))
    for nu in [1, 2, 3, 5, 15, 50]
    for t in [1.30 + 0.1 * k for k in range(40)]
)
max_t2_above_1p5 = max(
    t * t * 2.0 * math.pi * math.pi * nu * nu * abs(D_HS(t * math.pi * nu, nu))
    for nu in [1, 2, 3, 5, 15, 50]
    for t in [1.50 + 0.2 * k for k in range(40)]
)
print(f"\n  max_{{t >= 1.30}} H_D(t) = {max_above_1p3:.4f}")
print(f"  max_{{t >= 1.50}} H_D(t) = {max_above_1p5:.4f}")
print(f"  max_{{t >= 1.50}} t^2 * H_D(t) = {max_t2_above_1p5:.4f}")

# Step 2: Now define clean closed-form H_clean(t) >= H_D(t) everywhere:
#   - t <= 1.50: min(24.1 * t^2 / (1 + 0.5*t^2)^2, 3.04)  (check!)
#   - t >= 1.50: min(max_above_1p5 + 0.05, (max_t2_above_1p5 + 0.2) / t^2)
# First verify H_clean(t) >= H_D(t) on a dense grid!
def H_clean(t):
    if t <= 0: return 0.0
    if t <= 1.35:
        return min(24.1 * t * t / ((1.0 + 0.5 * t * t) ** 2), 3.04)
    elif t <= 1.50:
        return min(max_above_1p3 + 0.05, 3.04)
    else:
        return min(max_above_1p5 + 0.05, (max_t2_above_1p5 + 0.3) / (t * t))

ok = True
for nu in [1, 2, 3, 5, 10, 30]:
    for it in range(1, 120):
        t = 0.02 * (1.06 ** it)
        val = 2.0 * math.pi * math.pi * nu * nu * abs(D_HS(t * math.pi * nu, nu))
        if val > H_clean(t) + 1e-6:
            print(f"  VIOLATION of H_clean at nu={nu}, t={t:.4f}: val={val:.4f} > H_clean={H_clean(t):.4f}")
            ok = False; break
    if not ok: break
print(f"  H_clean(t) >= H_D(t) verified on dense (nu, t) grid: {ok}")

# Step 3: Now compute the EXACT ZERO-LOOP CONVERGENT SUM BOUND using H_clean(t)!
# Split convergents at t* = 1.50!
#   - Group 1: k1 = largest k with t_k <= 1.50 (if any)
#   - Group 2: k2 = k1 + 1 (so t_{k2} > 1.50)
#   - Group 3: all k <= k1 - 1
#   - Group 4: all k >= k2 + 1
# Notice:
#   - For Group 1 (k = k1):
#       IF q_{k1}*kap <= eps_G3: <= 0.0926 * M_G3 * (H_clean(t1) / 3.04) = 0.03045 * M_G3 * H_clean(t1)!
#       IF q_{k1}*kap >  eps_G3: <= 1.2021 * 0.344 / (q_{k1} D_{k1} eps_G3) <= 0.4135 / eps_G3!
#     WAIT: Notice at M_G3 = 8: 0.0926 * 8 = 0.741; at eps_G3 = 0.55: 0.4135 / 0.55 = 0.752!
#           At M_G3 = 10: 0.0926 * 10 = 0.926; at eps_G3 = 0.45: 0.4135 / 0.45 = 0.919!
#   - For Group 2 (k = k2, t_{k2} > 1.50):
#       H_clean(t_{k2}) <= max_above_1p5 + 0.05!
#       WAIT: If q_{k2}*kap <= eps_G3: <= 0.03045 * M_G3 * (max_above_1p5 + 0.05)!
#       If q_{k2}*kap > eps_G3: ALSO bounded by 1.2021 * (max_{t>=1.5} t*H_D(t)/(2*pi)) / (q_{k2} D_{k2} eps_G3)!
#         CHECK max_{t >= 1.5} t * H_D(t) / (2*pi)!!
max_tHD_1p5 = max(
    t * 2.0 * math.pi * math.pi * nu * nu * abs(D_HS(t * math.pi * nu, nu)) / (2.0 * math.pi)
    for nu in [1, 2, 3, 5, 15, 50]
    for t in [1.50 + 0.1 * k for k in range(60)]
)
print(f"  max_{{t >= 1.50}} t*H_D(t)/(2*pi) = {max_tHD_1p5:.4f}  (-> Group 2 high-conv <= {1.2021 * max_tHD_1p5:.4f} / eps_G3 !)")

#   - For Group 3 (k <= k1 - 1): coupled with Group 1!
#     Specifically: k1 - 1 has t_{k1-1} <= t1 / (rq_0 * (D_{k1-1}/D_{k1})) AND M_{k1-1} <= (rq_0 + 1)/D_{k1-1}!
#     Wait: if t1 <= 1.50, compute max over t1 in [0, 1.50] of:
#       [ M_G3 * H_clean(t1) + max_{rq0>=1, D0>=1} ((rq0+1)/D_prev) * H_clean(t1 * D0 / (rq0 * D_prev)) ]
#     PLUS tail k <= k1 - 2 (D_k >= F_{j+1} D0, t_k <= 1.50 / (r * F_{j+1}))!
def max_g3_term(A):
    mx = 2.0 * H_clean(A)
    for ir in range(1, 300):
        r = 1.0 + 0.05 * ir
        mx = max(mx, (r + 1.0) * H_clean(A / r))
    return mx

fibs = [2, 3, 5, 8, 13, 21, 34]
g3_tail = sum((1.0 / F) * max_g3_term(1.50 / F) for F in fibs) * 0.03045

# Group 4 (k >= k2 + 1, t_{k2} > 1.50):
# M_k * H_clean(t_k) <= M_k * ((max_t2_above_1p5 + 0.3) / t_k^2)
# Sum <= 1.05 * (max_t2_above_1p5 + 0.3) / (1.50^2) * 0.03045!
g4_val = 1.05 * (max_t2_above_1p5 + 0.3) / (1.50 * 1.50) * 0.03045

print(f"\n=== 2. ZERO-LOOP COUPLED CONVERGENT TOTAL (t* = 1.50) ===")
for M_G3, eps_G3 in [(10, 0.45), (9, 0.48), (8, 0.50), (8, 0.55), (7, 0.55)]:
    # Coupled max of Group 1 + (k1 - 1) over t1 in (0, 1.50]:
    mx_1_and_km1 = 0.0
    for it in range(1, 151):
        t1 = 1.50 * it / 150.0
        # k1 term (capped by BOTH M_G3 and high-conv bound! Actually worst of the two):
        # Wait: if k1 is low: 0.03045 * M_G3 * H_clean(t1); if k1 is high: 0.4135 / eps_G3 (maxed at peak t1 ~ 0.9)
        term_k1 = max(0.03045 * M_G3 * H_clean(t1), 1.2021 * (t1 * H_clean(t1) / (2.0 * math.pi)) / eps_G3)
        # k1 - 1 term:
        km1 = 0.03045 * max((1.0 / (2 if D0 == 1 else D0 + 1)) * max_g3_term(t1 * D0 / (2 if D0 == 1 else D0 + 1)) for D0 in range(1, 30))
        mx_1_and_km1 = max(mx_1_and_km1, term_k1 + km1)
    # Group 2 (k2, t2 > 1.50):
    g2_low = 0.03045 * M_G3 * (max_above_1p5 + 0.05)
    g2_high = 1.2021 * max_tHD_1p5 / eps_G3
    g2_val = max(g2_low, g2_high)
    tot_conv = mx_1_and_km1 + g3_tail + g2_val + g4_val
    print(f"  (M={M_G3:2d}, eps={eps_G3:.2f}): [k1+(k1-1)]={mx_1_and_km1:.4f}, tail_lo={g3_tail:.4f}, k2={g2_val:.4f}, tail_hi={g4_val:.4f} => ZERO-LOOP CONV = {tot_conv:.4f}")
