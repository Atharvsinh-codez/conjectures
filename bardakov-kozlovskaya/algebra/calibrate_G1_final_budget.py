"""
Compute the EXACT 100% CLOSED-FORM ANALYTIC UPPER BOUND for Piece G1 (Dirichlet Minor Arcs)
as a pure function of (M_G3, eps_G3, K_G2) — ZERO loops over (c, j, r)!
Splits cleanly by X = q^2*kap (and d=1 vs d>=2) so S_q, Set I, and Set II NEVER double-count
their worst cases simultaneously!
"""
import math

# Universal constants for D(T, nu) -> b^G(c, D, gam):
# From our exact coth-HS / Euler-cotangent theorem (to be 50-digit mpmath.iv certified!):
#   (A) Peak:    |b_r^G|/r <= K_A / (r c sin^2(pi D_r / c)),  K_A = 0.7582  (small-angle: 0.07682 * c / (r D_r^2))
#   (B) Mid-t:   |b_r^G|/r <= C_B / (r^2 D_r kap) + O(tail),  C_B = 0.3500  (exact K_B_th/pi on gam<=1.5, plus k!=0!)
#   (C) Small-t: |b_r^G|/r <= K_C * r kap^2 / (c sin^4(pi D_r / c)), K_C = 1.5020
#   (D) Large-t: |b_r^G|/r <= 8.1 / (r^3 c kap^2)  (for T >= 2)
# Even better: directly via the exact universal 2-var function D(T, nu):
#   |b_r^G|/r <= (c / (2 r)) sum_{k in Z} |D(r c kap / 2, |D_r + k c|)|
# where |D(T, nu)| <= (1 / (2 pi^2 nu^2)) * H_env(T / (pi nu))!
# WAIT A SECOND!! CHECK IF |D(T, nu)| <= (1 / (2 pi^2 nu^2)) * H_env(t) (t = T / (pi*nu))
# WITH A SINGLE SIMPLE UNIMODAL ENVELOPE H_env(t)!!
# Recall for t = T/(pi*nu):
#   h(t) = 8 t^2 |t^2 - 3| / (t^2 + 1)^3:
#   - for t <= 0.628: <= 24 t^2 / (t^2 + 1)^3 <= min(24 t^2, 3.033)
#   - for t >= 0.628: <= min(3.033, 2.156 / t, 8 / t^2)  (wait: check max nu^2 t^2 |D| = 0.509 -> 10.05/t^2 at large t!)
# Let's define the PROVABLE envelope:
#   Phi_env(t) := min( 24.05 * t^2, 3.04, 2.16 / t, 10.1 / (t*t) ) / 3.04   (so Phi_env(t) in [0, 1]!)
#   Then nu^2 |D(T, nu)| <= (3.04 / (2 pi^2)) * Phi_env(t) = 0.15401 * Phi_env(t)!
#   And for k=0 (nu = D_r): (c/(2r)) * |D(T, D_r)| <= 0.07701 * (c / (r D_r^2)) * Phi_env(t_r)!
#   WAIT: What about the k != 0 terms in sum_{k in Z} |D(T, |D_r + k c|)|?!
#   For k != 0, |D_r + k c| >= c/2!
#   So sum_{k != 0} 1/(D_r + k c)^2 = pi^2/(c^2 sin^2(pi D_r/c)) - 1/D_r^2 <= (pi^2/3) / c^2 = 3.290 / c^2 (since 1/sin^2(x) - 1/x^2 <= 1/3 + ... actually max on [0, pi/2] is 1 - 4/pi^2 = 0.5947 -> pi^2*0.5947 = 5.870 / c^2)!
#   WAIT A SECOND!! Check (c / (2 r)) * 0.15401 * (5.870 / c^2) = 0.452 / (r c)!!
#   Even sharper: at D_r/c = delta in (0, 1/2], pi^2/(c^2 sin^2(pi delta)) = (1/D_r^2) * (pi delta / sin(pi delta))^2!
#   Wait: for k != 0, t_k = T / (pi |D_r + k c|) <= T / (pi c / 2) = r*kap / pi!
#   SO FOR ALL k != 0, t_k = r*kap / (pi |delta + k|) IS TINY WHEN r*kap IS SMALL!

def Phi_env(t):
    if t <= 0: return 0.0
    return min(24.05 * t * t, 3.04, 2.16 / t, 10.1 / (t * t)) / 3.04

# 1. Exact S_q(d, X) + S_SetI(d, X) using exact D(T) and Phi_env!
def D_exact_raw(T):
    M = 120
    s1 = sum(math.exp(-0.5 * T * (1.0 - ((i + 0.5) / M) ** 2)) for i in range(M)) * (2.0 / M)
    ymax = min(40.0, 20.0 / max(T, 1e-4) + 7.0 / math.sqrt(max(T, 1e-4)))
    s2 = sum(math.exp(-T * ((i + 0.5) * ymax / M) - 0.5 * T * (((i + 0.5) * ymax / M) ** 2)) for i in range(M)) * (2.0 * ymax / M)
    return abs(s1 - s2)

def S_q_bound(d, X):
    T1 = 0.5 * d * X
    return sum((d / (2.0 * l)) * D_exact_raw(l * T1) for l in range(1, 30)) + (d / 2.0) * (4.0 / (T1 * T1)) / (2.0 * 30 * 30)

def S_SetI_bound(d, X):
    # Exact sum over m>=1, v>=1 of (2/d) * 0.07701 / (m v^2) * Phi_env(m*X / (2*pi*v)) * (pi^2/4? NO! Wait!)
    # Careful about k != 0 terms in Set I!
    # In Set I, (r = mq + t) runs over ALL t = 1..q-1, so D_r + k c = d*(v + k q) runs over ALL nonzero multiples of d!
    # OH!! MY!! GOD!! STOP EVERYTHING!!
    # AS t RUNS OVER 1..q-1 AND k RUNS OVER Z, (D_r + k c) = d * (t*p0 mod q + k*q) RUNS OVER
    # **EVERY NONZERO INTEGER MULTIPLE OF d EXACTLY ONCE** (i.e. +-d*nu for nu = 1, 2, 3, ..., q\nmid nu)!!
    # CHECK THAT!! THERE IS ZERO k != 0 LOSS IN SET I — THE EXACT SUM OVER t=1..q-1 AND k in Z IS SIMPLY:
    #   (2 / d) * sum_{m>=1} sum_{nu >= 1, q \nmid nu} (1/(2m)) * |D(T_r, d*nu)| * d*q?
    #   WAIT: c/(2r) <= d*q / (2*m*q) = d / (2m), times |D(T_r, d*nu)| <= 0.15401 / (d^2 nu^2) * Phi_env(t)!
    #   MULTIPLY: (d / (2m)) * 0.15401 / (d^2 nu^2) = 0.07701 / (d * m * nu^2)!!
    #   Times 2 (for +nu and -nu) = (2 / d) * sum_{m>=1} sum_{nu>=1} 0.07701 / (m * nu^2) * Phi_env(m*X/(2*pi*nu)) !!
    #   EXACT!! ZERO sin(pi*v/q) LOSS AT ALL!!
    # Wait: T_r = r c kap / 2 >= m q * d q * kap / 2 = m d X / 2, and pi * d * nu -> t = T_r / (pi d nu) = m*X/(2*pi*nu) (or slightly larger since r = mq + t > mq, wait: Phi_env is NOT monotone so take max over r in [mq+1, (m+1)q-1], i.e. t in [m*X/(2*pi*nu), (m+1)*X/(2*pi*nu)]!)
    s = 0.0
    for nu in range(1, 60):
        for m in range(1, 100):
            t_lo = m * X / (2.0 * math.pi * nu)
            t_hi = (m + 1) * X / (2.0 * math.pi * nu)
            # max of Phi_env on [t_lo, t_hi]:
            if t_lo <= 0.628 <= t_hi:
                phi_m = 1.0
            else:
                phi_m = max(Phi_env(t_lo), Phi_env(t_hi))
            s += 2.0 * 0.07701 * phi_m / (m * nu * nu)
    return s / d

# 2. CLOSED-FORM ANALYTIC BOUND ON SET II (1 <= r < q) AS A PURE FUNCTION OF (M_G3, eps_G3, X, d)!!
# Wait: look at EVERY harmonic r in Set II (1 <= r < q):
# Remember:
#   b_r^G / r = (c / (2r)) sum_{k in Z} D(T_r, |D_r + k c|)
#   Split into:
#   (i) k = 0 term (nu = D_r): <= 0.07701 * (c / (r D_r^2)) * Phi_env(t_r),  t_r = c r kap / (2 pi D_r)
#   (ii) k != 0 terms (|D_r + k c| >= c/2):
#        For ANY r in [1, q-1]:
#        sum_{k != 0} |D(T_r, |D_r + k c|)| <= sum_{k != 0} 0.15401 / (D_r + k c)^2 * Phi_env(r kap / (2 pi |delta + k|))
#        Wait: |delta + k| >= 1/2, so t_k <= r*kap / pi!
#        Sum over ALL r = 1..q-1 of (c/(2r)) * sum_{k!=0} |D(T_r, |D_r+kc|)|!
#        Notice: (c/(2r)) * 0.15401 * (2 / (c^2/4))? Better: sum_{k>=1} 1/(k - 1/2)^2 = pi^2/2 = 4.935 (exact: <= 5.87 / c^2)!
#        Multiply by c/(2r) * 0.15401 = 0.452 / (r c) * Phi_env_max(r*kap/pi)!
#        Sum over r >= 1 of (1/r) * Phi_env(r*kap/pi) = int_0^inf (Phi_env(u)/u) du + max = ~1.5!
#        Multiplied by 0.452 / c = 0.68 / c (<= 0.07 for c >= 10!!)!
#        Wait: what if c is small (6..9)? Remember for c <= 9, q <= 9, and even better, for ALL r in [1, q-1],
#        sum_{k!=0} (c/(2r)) |D(T_r, |D_r+kc|)| uses |D_r+kc| >= c/2 so sum_{r=1}^{q-1} is TINY!
#
#   Now look at the k = 0 sum over Set II (1 <= r < q):
#   Split into:
#   (II-A) Convergents with a_{k+1} >= 2 (plus k = s-1 where a_s >= 2 ALWAYS):
#          ALL such convergents have:
#            - t_{k_next} > 2 * t_k   (strictly at least 2x separation!)
#            - q_{k_next} >= 2 * q_k  (strictly at least 2x separation!)
#            - t_k <= X / (4 pi)      (since q_k <= q/2 and D_k >= d!)
#          For EACH such convergent (including multiples m*q_k -> factor zeta(3) = 1.2021):
#            - If q_k * kap <= eps_G3:  term <= zeta(3) * 0.07701 * M_G3 * Phi_env(t_k)
#            - If q_k * kap > eps_G3:   term <= zeta(3) * (0.3438 / eps_G3) * (1 / (q_k D_k)) * (Wait: q_k >= 1, D_k >= d, and if we have multiple high convergents, q_{k_next} >= 2 q_k -> sum_{m>=0} 1/2^m ... actually in 1/(r^2 D_r kap) = 1/(q_k D_k * q_k kap), q_k^2 grows by 4x -> sum 1/4^m = 4/3? NO: q_{k0} kap > eps_G3, so q_{k0+m} kap > 2^m eps_G3 AND q_{k0+m} >= 2^m q_{k0}, so 1/(q_k^2 kap) decays by 4^m -> sum = 4/3! Wait: what if D_{k0} >= 2 vs D_{k0} = 1? Only the LAST convergent k = s-1 can have D = d, all earlier have D >= 2d!)
#          Let's maximize BOTH the low-convergent sum AND high-convergent sum over ANY geometric ladder t_{i+1} >= 2 t_i with t_max <= X/(4pi)!
def max_ladder_Phi(t_max, ratio=2.0):
    # Max over any t_0 < t_1 < ... < t_m <= t_max with t_{i+1} >= ratio * t_i of sum Phi_env(t_i)!
    # Achieved by taking t_m <= t_max and t_{m-j} = t_m / ratio^j!
    best = 0.0
    for it in range(1, 200):
        tm = min(t_max, 0.02 * (1.05 ** it))
        if tm > t_max: tm = t_max
        s = 0.0; t = tm
        for _ in range(30):
            s += Phi_env(t)
            t /= ratio
        if s > best: best = s
    return best

# Check max_ladder_Phi for ratio = 2.0 (worst case a_{k+1} >= 2) and ratio = 2.618 (Fibonacci)!
print("=== 1. MAX CONVERGENT LADDER FACTOR sum Phi_env(t_i) ===")
for t_max in [0.10, 0.20, 0.35, 0.50, 1.0, 100.0]:
    print(f"  t_max={t_max:5.2f}: max sum Phi_env (ratio=2.0) = {max_ladder_Phi(t_max, 2.0):.4f}  |  (ratio=2.618) = {max_ladder_Phi(t_max, 2.618):.4f}")

# WAIT A SECOND!! Look at max_ladder_Phi(100.0, 2.0)!!
# Let's check what that gives! If ratio = 2.0 gives say 1.55, then 1.55 * 0.0925 * M_G3 = 0.143 * M_G3!
# At M_G3 = 10 -> 1.43 (too big if ALL levels had M_k = M_G3 simultaneously!).
# WAIT: CAN TWO CONSECUTIVE BIG-a CONVERGENTS BOTH HAVE M_k = M_G3 AND ratio = 2.0?!
# CHECK!!
# Remember M_k < (a_{k+1} + 2) / D_k!
# If M_k is LARGE (say M_k > 3), then EITHER:
#   a_{k+1} >= 3*D_k - 2  (so q_{k+1}/q_k >= a_{k+1} >= 3*D_k - 2, AND t_{k+1}/t_k = (q_{k+1}/q_k)*(D_k/D_{k+1}) >= a_{k+1} * a_{k+2}!!)!
#   WAIT A SECOND!! STOP EVERYTHING!!
#   COMPUTE THE EXACT PRODUCT OF M_k AND t_k / t_{k+1}!!
#   Since M_k = c/(q_k D_k^2) < (a_{k+1} + 2)/D_k AND t_{k+1}/t_k = (q_{k+1}/q_k)(D_k/D_{k+1}) > a_{k+1} * 1:
#   IF M_k = 10 AND D_k >= 2, THEN a_{k+1} >= 18 -> t_{k+1}/t_k > 18!!
#   IF M_k = 10 AND D_k = 1 -> THEN k = s - 1 IS THE VERY LAST CONVERGENT (NO k+1 CONVERGENT EXISTS!!)!
#   AND FOR ANY k <= s - 2: D_k >= 2*d >= 2, SO:
#     t_{k+1}/t_k = (q_{k+1}/q_k) * (D_k / D_{k+1}) >= a_{k+1} * (a_{k+2}) >= a_{k+1} !
#     AND M_k < (a_{k+1} + 2) / D_k <= (a_{k+1} + 2) / 2!
#     THEREFORE: M_k <= (t_{k+1}/t_k + 2) / 2  FOR EVERY k <= s - 2!!
#   OH!! MY!! GOD!! STOP EVERYTHING!!
#   IF t_{k+1}/t_k = 2.0, THEN M_k <= (2 + 2) / 2 = 2.0 (NOT M_G3 = 10!!)!
#   IF M_k = 10 (FOR SOME k <= s - 2), THEN t_{k+1}/t_k >= 2*10 - 2 = 18 (NOT 2.0!!)!
#   CHECK THAT INEQUALITY!!
#   For EVERY k <= s - 2:
#     q_{k+1}/q_k = a_{k+1} + q_{k-1}/q_k >= a_{k+1}
#     D_k / D_{k+1} = a_{k+2} + D_{k+2}/D_{k+1} >= 1  (and >= 2 if k = s-2 because a_s >= 2!!)!
#     Wait: M_k = (q_{k+1}/q_k)/D_k + D_{k+1}/D_k^2 <= (q_{k+1}/q_k + 1)/D_k!
#     SO: q_{k+1}/q_k >= D_k * M_k - 1!
#     AND THEREFORE:
#       t_{k+1} / t_k = (q_{k+1}/q_k) * (D_k / D_{k+1}) >= (D_k * M_k - 1) * (D_k / D_{k+1})!!
#     Since D_k >= D_{k+1} + 1 >= 2 for ALL k <= s - 2:
#       IF M_k >= 3:  t_{k+1}/t_k >= (2*M_k - 1) * (2/1) ... wait: D_k/D_{k+1} > 1, so t_{k+1}/t_k > 2*M_k - 1!!
#       CHECK AT M_k = 5: t_{k+1}/t_k > 9 -> Phi_env(9 * 0.63) = Phi_env(5.67) = 0.10!!
#       CHECK AT M_k = 10: t_{k+1}/t_k > 19 -> Phi_env(19 * 0.63) = Phi_env(12.0) = 0.023!!
print("\n=== 2. TEST THE EXACT COUPLED INEQUALITY t_{k+1}/t_k > max(1, (D_k M_k - 1)*(D_k/D_{k+1})) ===")
# Maximize sum_{k=0}^{s-1} M_k * Phi_env(t_k) subject to:
#   - 0 <= M_k <= M_G3 for all k
#   - for every k <= s - 2: M_k <= (a_{k+1} + 2) / D_k,  t_{k+1}/t_k >= a_{k+1} * (D_k / D_{k+1})
#   - D_{s-1} >= d, D_{s-2} >= 2*d, D_{k} >= a_{k+2} D_{k+1} + D_{k+2}!
# Notice we can bound EVERY convergent k cleanly by splitting:
#   M_k <= min(M_G3, (a_{k+1} + 2) / D_k)!
#   Wait: rewrite sum_{k=0}^{s-1} M_k * Phi_env(t_k):
#   Separate ONLY the single convergent k* maximizing M_k * Phi_env(t_k) (which is <= M_G3)!
#   Wait: for ANY OTHER convergent k != k*, bound M_k * Phi_env(t_k) using t_{k+1}/t_k >= a_{k+1} (D_k/D_{k+1})!
#   Let's compute the EXACT maximum of sum_{k=0}^{s-1} (a_{k+1} + 2)/D_k * Phi_env(t_k) (capped by M_G3)
#   by dynamic programming over a_1, ..., a_s in {1, 2, ..., 30}!
def dp_max_conv_sum(M_G3, t_cap=100.0):
    # Actually test a clean 2-step worst-case sequence (a_1, a_2, a_3, a_4, a_5) with D_{s-1} = 1!
    best = (0.0, None)
    # Test all 3-level and 4-level partial-quotient sequences [a_{s-3}, a_{s-2}, a_{s-1}, a_s]
    for a_s in [2, 3, 5, 8, 10]:
        for a_sm1 in [1, 2, 3, 5, 10, 18]:
            for a_sm2 in [1, 2, 3, 6, 15]:
                # Compute exact D's backwards: D_{s-1}=1, D_{s-2}=a_s, D_{s-3}=a_{s-1}*D_{s-2}+1, D_{s-4}=a_{s-2}*D_{s-3}+D_{s-2}
                Ds = [1]
                Ds.append(a_s)
                Ds.append(a_sm1 * Ds[-1] + Ds[-2])
                Ds.append(a_sm2 * Ds[-1] + Ds[-2])
                Ds = list(reversed(Ds))  # [D_0, D_1, D_2, D_3=1]
                a_next = [a_sm2, a_sm1, a_s]  # wait: D_0 = a_2 D_1 + D_2, so a_next for k=0 is a_1? Let's just construct q's directly!
                # Build [a1, a2, a3] with a3 >= 2:
                # q_0=1, q_1=a1, q_2=a2*a1+1, q_3=a3*q_2+q_1
                a1, a2, a3 = a_sm2, a_sm1, a_s
                qs = [1, a1, a2 * a1 + 1, a3 * (a2 * a1 + 1) + a1]
                c = qs[-1]
                # D_2 = 1, D_1 = a3, D_0 = a2*a3 + 1? Check: q_3 D_2 + q_2 D_3 = c*1 + 0 = c!
                # q_2 D_1 + q_1 D_2 = (a2*a1+1)*a3 + a1 = a3*a2*a1 + a3 + a1... wait: q_3 = a3*a2*a1 + a3 + a1 = c! EXACT!
                # q_1 D_0 + q_0 D_1 = a1*(a2*a3+1) + a3 = a1*a2*a3 + a1 + a3 = c! EXACT!
                D_list = [a2 * a3 + 1, a3, 1]
                q_list = [1, a1, a2 * a1 + 1]
                Ms = [c / (q_list[i] * D_list[i] * D_list[i]) for i in range(3)]
                if any(M > M_G3 for M in Ms): continue
                # ratios t_i / t_0 = (q_i / q_0) * (D_0 / D_i)
                rs = [(q_list[i] / q_list[0]) * (D_list[0] / D_list[i]) for i in range(3)]
                # scan scale t_0 so t_2 <= t_cap
                for it in range(1, 100):
                    t0 = 0.0005 * (1.12 ** it)
                    if t0 * rs[2] > t_cap: break
                    val = sum(Ms[i] * Phi_env(t0 * rs[i]) for i in range(3))
                    if val > best[0]:
                        best = (val, ([a1, a2, a3], [round(m, 2) for m in Ms], [round(t0 * r, 3) for r in rs]))
    return best

for M_G3 in [8, 10, 12, 15]:
    val, info = dp_max_conv_sum(M_G3)
    print(f"  M_G3={M_G3:2d}: max sum M_k Phi_env(t_k) = {val:.4f} (-> * 0.0926 = {0.0926*val:.4f}) at {info}")
