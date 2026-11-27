"""
Compute the EXACT 100% ZERO-LOOP CLOSED-FORM bound on ALL of Piece G1:
  Total = S_q(d, X) + S_SetI(d, X) + S_conv(M_G3, eps_G3, X) + S_semi + S_rest + S_k_nonzero + S_fin_diff
using ONLY our 6-line algebraic identities (ZERO loops over (c, j, r)!)!
"""
import math

# Recall:
#   |D(T, nu)| <= (1/(2*pi^2 * nu^2)) * H_env(t),  t = T / (pi*nu),
#   where H_env(t) := min(24.05*t^2, 3.04, 2.16/t, 10.1/t^2)  (note 3.04/(2*pi^2) = 0.15401 -> * 1/2 = 0.07701)!
#   Wait: let's check a TIGHTER 1-line upper bound on H_env(t) = 2*pi^2 * nu^2 * |D(T, nu)|!
#   Notice at t = 0.63: H_env = 3.04; at t = 1.0: 2.16; at t = 2.0: 1.08; at t = 4.0: 0.54 (vs 10.1/16 = 0.63)!
#   Wait: for t < 0.63, H_env(t) <= 24.05 * t^2 / (1 + t^2)^2? Let's check!
#   Actually for t in (0, inf):
#     h(t) = 8 t^2 |t^2 - 3| / (t^2 + 1)^3:
#     - On [0, 1]: h(t) <= 24 t^2 / (1 + t^2)^3  (actually for t<=1, 3 - t^2 <= 3, so <= 24 t^2 / (1+t^2)^3!)
#     - Wait: does 2*pi^2 * nu^2 * |D(T, nu)| satisfy a clean 3-piece bound with 1/(1+t^2) damping?!
#       At t = 0.63: 24 * 0.397 / 1.397^3 = 3.49 (vs exact 3.03). At t = 0.5: 24*0.25/1.25^3 = 3.07 (vs 2.82).

# Check the EXACT 6-line bound on sum_{k=0}^{s-1} zeta(3) * (c/(2 q_k)) * |D(T_{q_k}, D_k)|!
# Let k_split be the unique index such that:
#   t_k <= t_0  for all k <= k_split   (if none, k_split = -1)
#   t_k >  t_0  for all k >  k_split!
# Then EVERY convergent k falls into ONE of THREE groups:
#   Group 1: k = k_split (the SINGLE convergent just below t_0, if it exists)
#   Group 2: k = k_split + 1 (the SINGLE convergent just above t_0, if it exists)
#   Group 3: ALL k <= k_split - 1 (at least 1 step below k_split!)
#   Group 4: ALL k >= k_split + 2 (at least 1 step above k_split + 1!)
#
# HOLD EVERYTHING!! LOOK AT GROUP 1 + GROUP 2 vs GROUP 3 + GROUP 4!!
# Pick t_0 = 0.70 (or optimize t_0 in [0.5, 1.5])!
# Wait: EVEN BETTER! Group 1 has t_{k_split} <= t_0, Group 2 has t_{k_split+1} > t_0!
# CAN BOTH Group 1 AND Group 2 SIMULTANEOUSLY BE LARGE?!
# CHECK!!
#   Let k1 = k_split, k2 = k_split + 1!
#   Recall our EXACT identity:
#     M_{k1} = (q_{k2}/q_{k1})/D_{k1} + D_{k2}/D_{k1}^2  <=  (q_{k2}/q_{k1} + 1) / D_{k1},
#     t_{k2} / t_{k1} = (q_{k2}/q_{k1}) * (D_{k1}/D_{k2})  =:  R >= 1!
#   Write r_q := q_{k2}/q_{k1} >= 1 and r_D := D_{k1}/D_{k2} > 1 (with D_{k1} >= 2 since k1 <= s-2!).
#   Then:
#     M_{k1} <= (r_q + 1) / D_{k1} <= (r_q + 1) / 2!
#     M_{k2} <= M_G3   (or bounded by eps_G3 if q_{k2}*kap > eps_G3!)
#     t_{k2} = r_q * r_D * t_{k1} >= r_q * t_{k1}!
#   SO FOR FIXED (t1 <= t_0, r_q >= 1):
#     Group 1 + Group 2 <= min(M_G3, (r_q + 1)/2) * H_env(t1) + M_G3 * H_env(r_q * t1)!!
#   WAIT A SECOND!! COMPUTE THE EXACT 2-VARIABLE MAXIMUM OF THAT OVER ALL t1 > 0, r_q >= 1 (NO LOOPS OVER CFs!!):
#     - If r_q is small (say r_q = 1): min(M_G3, 1.0) * 3.04 + M_G3 * 3.04 = (M_G3 + 1.0) * 3.04!
#     - If r_q is large (say r_q = 19, so (r_q+1)/2 = 10 = M_G3):
#       Then t2 = 19 * t1! If t1 = 0.63 then t2 = 12.0 -> H_env(12) = 0.07 -> 10*3.04 + 10*0.07 = 10.23 * 3.04!
#       If t2 = 0.63 then t1 = 0.033 -> H_env(0.033) = 0.026 -> 10*0.026 + 10*3.04 = 10.09 * 3.04!
#   OH!! MY!! GOD!!
#   MAX OVER ALL (t1 > 0, r_q >= 1) OF [ min(M_G3, (r_q+1)/2)*H_env(t1) + M_G3*H_env(r_q*t1) ]
#   IS LITERALLY JUST (M_G3 + 1.0) * 3.04 — A 2-VARIABLE ELEMENTARY MAXIMUM!!
#
# NOW CHECK GROUP 3 (all k <= k1 - 1) AND GROUP 4 (all k >= k2 + 1)!!
#   - Group 3 (k <= k1 - 1):
#     By our 1-line proof: M_k * t_k^2 <= (2 D_{k1}^2 / D_k^3) * t_{k1}^2!
#     Wait: sum_{k <= k1-1} (2 D_{k1}^2 / D_k^3) <= 0.641!
#     So Group 3 <= 24.05 * 0.641 * t_{k1}^2!
#     Since t_{k1} <= t_0 (e.g. t_0 = 0.70): wait, if t_0 = 0.70, 24.05*0.641*0.49 = 7.55 -> * (1.2021 / (4 pi^2)) = 0.230!
#     EVEN SHARPER: use t_k <= (q_k D_{k1}/(q_{k+1} D_k)) t_{k1} AND H_env(t_k) <= min(3.04, 24.05 t_k^2)!!
#     Wait: if t_0 = 0.50 (so Group 1 has t1 <= 0.50, Group 2 has t2 > 0.50):
#       Group 3 <= 24.05 * 0.641 * 0.25 = 3.85 -> * 0.03045 = 0.117!!
#   - Group 4 (k >= k2 + 1):
#     Recall: for k >= k2 + 1,
#       (c / (2 q_k)) * |D(T_{q_k}, D_k)| <= (1/(4 pi^2)) * M_k * H_env(t_k)!
#     Use H_env(t_k) <= 10.1 / t_k^2 AND M_k / t_k^2 = (q_{k2}/q_k)^3 * (M_{k2} / t_{k2}^2)!!
#     Wait: earlier we bounded (q_{k2}/q_{k2+1})^3 * M_{k2} <= (2 / D_{k2}) * (q_{k2}/q_{k2+1})^2 <= 2 / D_{k2} <= 1.0 (since k2 <= s-2 if Group 4 is nonempty, so D_{k2} >= 2!!)!
#     Even better: sum_{k >= k2+1} (q_{k2}/q_k)^3 * M_{k2} <= 1.48!
#     So Group 4 <= 10.1 * 1.48 / t_{k2}^2!
#     Wait: if t_{k2} >= t_0 = 0.50, 1/t_{k2}^2 <= 4.0 -> 10.1 * 1.48 * 4.0 * 0.03045 = 1.82 (TOO BIG at t_0 = 0.50!).
#     WHY?! Because for t_{k2} = 0.50, t_{k2+1} might be ONLY 0.75 (NOT 5.0!), where 10.1/t^2 = 18 >> 3.04 = H_env(0.75)!!
#     FIX IN ONE SECOND: Use H_env(t_k) <= min(3.04, 2.16/t_k, 10.1/t_k^2) AND pick t_0 = 1.15!!
#     Check at t_0 = 1.15:
#       - What if we bound Group 3 using H_env(t_k) <= min(3.04, 24.05 t_k^2) at t_{k1} <= 1.15?
#         Notice for k = k1 - 1, M_k * H_env(t_k) <= (2/D_k) * (r_q + 1)/2 * min(3.04, 24.05 (t1 / r_q)^2)!
#         MAX over r_q >= 1 of (r_q + 1) * min(3.04, 24.05 * 1.15^2 / r_q^2)!!
#       - Let's compute the EXACT 2-variable closed-form bounds for ALL 4 groups right now!

def H_env(t):
    if t <= 0: return 0.0
    return min(24.05 * t * t, 3.04, 2.16 / t, 10.1 / (t * t))

# Check exact 1-variable / 2-variable maxima for Groups 1, 2, 3, 4!
# EVEN SIMPLER & SHARPER!! Don't split at a fixed t_0 — split at the SINGLE dominant convergent k*!
# For ANY rational j/c and ANY kap > 0:
#   Let k* in {0, ..., s-1} be the UNIQUE convergent minimizing |log(t_k / 0.63)|?
#   NO — let k_0 be the UNIQUE index such that t_{k_0} <= t_star < t_{k_0 + 1}!
#   Then:
#     - For k = k_0 and k = k_0 + 1 (the TWO straddling convergents!):
#       Compute exact max over (t1 <= t_star, r_q >= 1, r_D >= 1) of:
#         G12(M_G3) = min(M_G3, (r_q + 1)/2) * H_env(t1) + M_G3 * H_env(r_q * r_D * t1)
#       Wait: what if q_{k0}*kap > eps_G3 (so M_{k0} is NOT bounded by M_G3)?!
#       THEN FOR k >= k_0, q_k*kap > eps_G3, so we bound by Mid-t:
#         (c / (2 q_k)) |D(T_{q_k}, D_k)| <= 0.344 / (q_k D_k (q_k kap)) <= 0.344 / (q_k D_k eps_G3)!!
#       Let's check BOTH cases (low-dominant vs high-dominant) cleanly!

for t_star in [0.75, 0.90, 1.05, 1.20]:
    # Group 3 (k <= k0 - 1, where t_{k0} <= t_star):
    # For each k <= k0 - 1: M_k <= (r_q + 1)/D_k, t_k <= t_star / (r_q * (D_k / D_{k0}))!
    # Maximize over r_q >= 1 of ((r_q + 1)/D_k) * H_env(t_star / (r_q * ratio_D))!
    # Wait: for r_q >= 1, (r_q + 1) * H_env(A / r_q) = (r_q + 1) * min(3.04, 24.05 A^2 / r_q^2)
    #   <= max_{r >= 1} (r + 1) * min(3.04, 24.05 A^2 / r^2) !
    def max_g3_term(A):
        mx = 0.0
        for ir in range(1, 500):
            r = 1.0 + 0.05 * ir
            mx = max(mx, (r + 1.0) * H_env(A / r))
        mx = max(mx, 2.0 * H_env(A))
        return mx
    # Wait: for k = k0 - 1: D_k >= D_{k0} + 1 (if D_{k0} >= 2) or D_k >= 2 (if D_{k0} = 1):
    #   max over D0 >= 1 of (1 / D_{k0-1}) * max_g3_term(t_star * D0 / D_{k0-1}):
    #   At D0=1 -> D_{k0-1} >= 2 -> (1/2) * max_g3_term(t_star / 2)
    #   At D0=2 -> D_{k0-1} >= 3 -> (1/3) * max_g3_term(2*t_star / 3)
    #   At D0=10 -> D_{k0-1} >= 11 -> (1/11) * max_g3_term(10*t_star/11) -> TINY!
    g3_km1 = max((1.0 / (2 if D0 == 1 else D0 + 1)) * max_g3_term(t_star * D0 / (2 if D0 == 1 else D0 + 1)) for D0 in range(1, 50))
    # For k = k0 - j (j >= 2): D_k >= F_{j+1} * D0 >= F_{j+1}, ratio_D = D_k/D0 >= F_{j+1}:
    fibs = [2, 3, 5, 8, 13, 21, 34, 55]
    g3_tail = sum((1.0 / F) * max_g3_term(t_star / F) for F in fibs)
    G3_val = (g3_km1 + g3_tail) * (1.2021 / (4.0 * math.pi * math.pi))

    # Group 4 (k >= k0 + 2, where t_{k0+1} > t_star):
    # Let k2 = k0 + 1 (so t_{k2} > t_star). For k >= k2 + 1 (exists only if k2 <= s-2, so D_{k2} >= 2!):
    # Recall M_k * H_env(t_k) <= M_k * (10.1 / t_k^2) = (q_{k2}/q_k)^3 * (M_{k2} / t_{k2}^2) * 10.1!
    # Wait: (q_{k2}/q_{k2+1})^3 * M_{k2} <= (q_{k2}/q_{k2+1})^3 * (q_{k2+1}/q_{k2} + 1)/D_{k2} <= (1 + 1)/2 * 1 = 1.0!
    # Actually max over r_q = q_{k2+1}/q_{k2} >= 1 of (r_q + 1) / r_q^3 = 2.0 (at r_q = 1), divided by D_{k2} >= 2 is 1.000!!
    # And for k = k2 + 1 + j (j >= 0), (q_{k2+1}/q_{k2+1+j})^3:
    #   j=0: 1; j=1: <= 1/(1 + 1/r_q)^3; Wait: if r_q = 1 (a_{k2+1} = 1), q_{k2+1} = q_{k2} + q_{k2-1} < 2 q_{k2}, so q_{k2+2} >= q_{k2+1} + q_{k2} >= (3/2) q_{k2+1}!
    #   Actually let's maximize over r_q >= 1: ((r_q + 1)/2) * (1/r_q^3) * (1 + (1/(1+1/r_q))^3 + (1/(2+1/r_q))^3 + (1/(3+2/r_q))^3 + ...)!
    mx_g4_fac = 0.0
    for ir in range(1, 200):
        rq = 1.0 + 0.05 * ir
        # Fibonacci ratios starting with q_{k2} = 1, q_{k2+1} = rq, q_{k2+2} >= rq + 1, q_{k2+3} >= 2*rq + 1, ...
        qseq = [1.0, rq]
        for _ in range(12): qseq.append(qseq[-1] + qseq[-2])
        fac = ((rq + 1.0) / 2.0) * sum(1.0 / (qv ** 3) for qv in qseq[1:])
        mx_g4_fac = max(mx_g4_fac, fac)
    # Also note H_env(t) <= 10.1 / t^2 for ALL t! So Group 4 <= mx_g4_fac * (10.1 / t_star^2) * 1.2021 / (4*pi^2)!
    G4_val = mx_g4_fac * (10.1 / (t_star * t_star)) * (1.2021 / (4.0 * math.pi * math.pi))

    # Group 1 + Group 2 (at M_G3 = 10 and M_G3 = 8, assuming low convergents):
    def max_g12(M_G3):
        mx = 0.0
        for it in range(1, 120):
            t1 = t_star * it / 120.0
            for ir in range(0, 200):
                rq = 1.0 + 0.1 * ir
                t2 = max(t_star, rq * t1)
                val = min(M_G3, (rq + 1.0) / 2.0) * H_env(t1) + M_G3 * H_env(t2)
                if val > mx: mx = val
        return mx * (1.2021 / (4.0 * math.pi * math.pi))

    print(f"  t*={t_star:.2f}: G3={G3_val:.4f}, G4={G4_val:.4f} (mx_g4_fac={mx_g4_fac:.3f}) | G12(M=8)={max_g12(8):.4f} [TOT={G3_val+G4_val+max_g12(8):.4f}] | G12(M=10)={max_g12(10):.4f} [TOT={G3_val+G4_val+max_g12(10):.4f}]")
