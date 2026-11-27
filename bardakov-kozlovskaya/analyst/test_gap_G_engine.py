"""
TEST ALL THREE ENGINES FOR THE GAP G = {c*kap^2 < 45, c*q*kap > 0.18}:
1. Positive-root 5-line saddle theorem (for D = 0, ALL g >= 1):
   Check minimum H_min for which our elementary global decay + cubic phase lower bound > 0!
   Also check actual H(lam_*) for ALL q=1 cases with c >= 2 and c*kap > 0.18!
2. Exact 2-step & 3-step Abel summation on |B_r - A_r|/r for Minor Arcs at lam = 1:
   Verify our closed-form bound on |B_r - A_r|/r against exact finite-n values across c = 20..500!
3. Check the EXACT small-c / small-tau overlap between 2C (c*q*kap <= C_2), lam=1, and Major Arcs!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd

# =============================================================================
# 1. POSITIVE-ROOT 5-LINE SADDLE LOWER BOUND AS A FUNCTION OF H
# =============================================================================
# Recall: for F(-lam_* e^{i theta}) = F(-lam_*) * prod_k (1 + x_k e^{i theta})/(1 + x_k) * e^{-i c theta},
# with ALL x_k > 0 and sum x_k/(1+x_k) = c, H = sum x_k/(1+x_k)^2:
#   - Upper bound on |F/F0| for ALL theta in [-pi, pi]:
#       |F/F0| <= exp(-2 * H * sin^2(theta/2))
#   - Lower bound on |F/F0| for |theta| <= pi/2:
#       Since log(1 - u) >= -u / (1 - u) for u in [0, 1) and 4*x/(1+x)^2 * sin^2(theta/2) <= sin^2(theta/2) <= 1/2:
#       (1/2) log(1 - u) >= -(u/2) * (1/(1 - sin^2(theta/2))) = -(2x/(1+x)^2) * tan^2(theta/2) * cos^2?
#       Wait: u = (4x/(1+x)^2)*sin^2(theta/2) <= sin^2(theta/2).
#       Since -log(1 - u)/u is increasing in u, for u <= s^2 = sin^2(theta/2):
#         (1/2) log(1 - u) >= (u / (2 s^2)) * log(1 - s^2) = (2x/(1+x)^2) * log(cos(theta/2))!
#       SUM OVER ALL k (EXACT!!):
#         |F/F0| >= (cos(theta/2))^{4 H}   FOR ALL |theta| < pi!!
#       OH!! MY!! GOD!! CHECK THAT 1-LINE EXACT LOWER BOUND ON |F/F0|!!
#         At x = 1: |1 + e^{i theta}| / 2 = cos(theta/2) = (cos(theta/2))^{4 * (1/4)}  (100% SHARP AT x = 1!!)
#         For ANY x > 0: let p = 4x/(1+x)^2 in (0, 1]. By convexity of t -> (1 - t*s^2) / (1 - s^2)^t?
#         By Bernoulli / concavity of log(1 - t*s^2): since log(1 - t*s^2) is CONCAVE in t on [0, 1],
#         it lies ABOVE the secant from t=0 (val 0) to t=1 (val log(1 - s^2) = 2*log(cos(theta/2))):
#           log(1 - p*s^2) >= p * log(1 - s^2) = 2 * p * log(cos(theta/2))!!
#         Multiply by 1/2 and sum over k (sum p_k = 4 H):
#           log|F(-lam_* e^{i theta}) / F(-lam_*)| >= 4 * H * log(cos(theta/2))  FOR ALL theta in (-pi, pi)!!
#   - Phase bound for |theta| <= pi/2:
#       |Phi(theta)| <= K_ph * H * |theta|^3, where K_ph = 1/3 universally (and smaller if |1-x|/(1+x) bounded)!
#       Even sharper! Integrate 4*sin^2(t/2) / (1 - sin^2(t/2))? Wait:
#         1 + 2x cos t + x^2 = (1+x)^2 * (1 - p sin^2(t/2)) >= (1+x)^2 * cos^2(t/2)!
#       So |d/dt (arg(1 + x e^{it}) - (x/(1+x))t)| <= (2x/(1+x)^2) * (1 - cos t) / (2 cos^2(t/2))
#                                                   = (2x/(1+x)^2) * tan^2(t/2)!
#       Integrate tan^2(t/2) from 0 to theta: int_0^theta tan^2(t/2) dt = 2*tan(theta/2) - theta!!
#       SUM OVER ALL k:
#         |Phi(theta)| <= 2 * H * (2 * tan(|theta|/2) - |theta|)   FOR ALL |theta| < pi!!
#       (Notice as theta -> 0: 2*(2 tan(t/2) - t) ~ 2*(t^3/12) = t^3/6 — FACTOR OF 2 BETTER THAN t^3/3,
#        AND VALID ON ALL (-pi, pi) WITH ZERO RESTRICTION TO |theta| <= pi/2!!)

print("=== 1. TEST EXACT 3-LINE POSITIVE-ROOT LOWER BOUND FOR ALL H > 0 ===")
ths = np.linspace(-pi, pi, 8000, endpoint=False)
ath = np.abs(ths)
c2 = np.cos(ath / 2.0)
s2 = np.sin(ath / 2.0)
# Exact phase upper bound: 2 * H * (2*tan(ath/2) - ath)  [where ath < pi]
t_fn = np.where(ath < pi - 1e-6, 2.0 * np.tan(ath / 2.0) - ath, 1e9)
for H in [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0, 50.0]:
    ph = 2.0 * H * t_fn
    lb = np.where(ph < pi/2, (c2 ** (4.0 * H)) * np.cos(ph), -np.exp(-2.0 * H * (s2 ** 2)))
    val = float(np.mean(lb))
    # Also compare with exact I / avg|F| lower bound:
    ub_L1 = float(np.mean(np.exp(-2.0 * H * (s2 ** 2))))
    print(f"  H={H:5.1f}: RIGOROUS I/F0 >= {val:+.6f} | avg|F/F0| <= {ub_L1:.6f} | rat >= {val/ub_L1:+.5f} {'PASS!' if val > 0 else 'fail'}")

# Also check actual H(lam_*) for q = 1 (D = 0) across c = 2..100 and c*kap in [0.15, 20]!
print("\n=== CHECK ACTUAL H(lam_*) FOR q = 1 (D = 0) AT SMALL c & SMALL c*kap ===")
for c in [2, 3, 4, 6, 8, 10, 20]:
    min_H = 1e9; worst_t = None
    for tau in [0.18, 0.30, 0.50, 1.0, 2.0, 5.0, 10.0]:
        kap = tau / c
        if kap >= 2.5: continue
        alpha = kap / c
        M = int(np.ceil(np.sqrt(50.0 / alpha))) + 10
        ms = np.arange(-(c + M), M + 1, dtype=np.float64)
        ex = -alpha * ms * (c + ms)
        # Solve sum 1/(1 + exp(-(ex + v))) = c
        lo, hi = -20.0, 20.0
        for _ in range(50):
            vm = 0.5 * (lo + hi)
            z = np.clip(ex + vm, -60.0, 60.0)
            if float(np.sum(1.0 / (1.0 + np.exp(-z)))) > c: hi = vm
            else: lo = vm
        v = 0.5 * (lo + hi)
        z = np.clip(ex + v, -60.0, 60.0)
        p = 1.0 / (1.0 + np.exp(-z))
        H = float(np.sum(p * (1.0 - p)))
        if H < min_H: min_H = H; worst_t = (tau, kap, exp(v), H)
    print(f"  c={c:2d}: min H(lam_*) over tau>=0.18, kap<2.5 is {min_H:.4f} at tau={worst_t[0]}, kap={worst_t[1]:.3f}, lam*={worst_t[2]:.3f}")
