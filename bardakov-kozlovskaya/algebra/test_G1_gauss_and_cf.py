"""
1. Check exact worst-case constants of b_r^G (Gaussian) and b_r^exact (finite-n) for
   K_A (peak), K_B (mid-t), K_C (small-t), AND check max sum |b_r^exact - b_r^G| / r across G1!
2. Compute the EXACT 100% CLOSED-FORM CF UPPER BOUND (pure math formula of M_G3, eps_G3, K_G2, ZERO loops over (c,j,r)!)
   and find the exact cleanest (M_G3, eps_G3, K_G2) that gives < 1.35 < pi/2!
"""
import math

def b_r_gauss(r, c, j, kap):
    th = 2.0 * math.pi * ((r * j) % c) / c
    B = sum(math.exp(-r * kap * m * (1.0 - m / c)) * math.cos(m * th) for m in range(1, c) if r * kap * m * (1.0 - m / c) < 40)
    A = 0.0; m = 1
    while True:
        ex = r * kap * m * (1.0 + m / c)
        if ex > 40: break
        A += 2.0 * math.exp(-ex) * math.cos(m * th); m += 1
    return B - A

print("=== 1. EXACT GLOBAL CONSTANTS FOR b_r^G (GAUSSIAN) vs b_r^main ===")
# Check max of:
#   gA = |b_r^G|/r * (r c sin^2(th/2))
#   gB = |b_r^G|/r * (r c gam sin(th/2))   [for gam <= 1.5]
#   gC = |b_r^G|/r * (c sin^4(th/2)) / (r kap^2)
#   gD = |b_r^G|/r * (r^3 c kap^2)         [for T = r*c*kap/2 >= 2]
gA = (0.0, None); gB = (0.0, None); gC = (0.0, None)
for c in [4, 6, 10, 15, 25, 50, 100]:
    for D in range(1, c // 2 + 1):
        th = 2.0 * math.pi * D / c
        s = math.sin(0.5 * th)
        for gam in [0.005, 0.02, 0.06, 0.15, 0.30, 0.50, 0.80, 1.20, 1.50, 2.50]:
            # Only depends on (c, D, gam) because r*kap = gam! Wait: b_r^G / r = b_1^G(c, D, gam)!!
            # CHECK THAT! In b_r^G, r ONLY enters as r*kap = gam and (r*j mod c) = D!
            # So |b_r^G| = |b_1^G(c, D, gam)| -> |b_r^G|/r * (r c ...) = |b_1^G(c, D, gam)| * (c ...)!!
            val = abs(b_r_gauss(1, c, D, gam))
            vA = val * c * s * s
            if vA > gA[0]: gA = (vA, (c, D, gam))
            if gam <= 1.5:
                vB = val * c * gam * s
                if vB > gB[0]: gB = (vB, (c, D, gam))
            vC = val * c * (s ** 4) / (gam * gam)
            if vC > gC[0]: gC = (vC, (c, D, gam))
print(f"  GAUSSIAN EXACT CONSTANTS (compare main: 0.7581, 1.0970, 1.5001):")
print(f"    gA = max |b_1^G| * c * sin^2(th/2)        = {gA[0]:.5f} at {gA[1]}")
print(f"    gB = max |b_1^G| * c * gam * sin(th/2)    = {gB[0]:.5f} at {gB[1]}")
print(f"    gC = max |b_1^G| * c * sin^4(th/2)/gam^2  = {gC[0]:.5f} at {gC[1]}")

print("\n=== 2. CHECK EXACT CLOSED-FORM SET I + S_q USING EXACT F(gam, th) (NO 1.74x BLOAT!) ===")
# Exact D(T)
def D_exact_raw(T):
    M = 150
    s1 = sum(math.exp(-0.5 * T * (1.0 - ((i + 0.5) / M) ** 2)) for i in range(M)) * (2.0 / M)
    ymax = min(40.0, 20.0 / max(T, 1e-4) + 7.0 / math.sqrt(max(T, 1e-4)))
    s2 = sum(math.exp(-T * ((i + 0.5) * ymax / M) - 0.5 * T * (((i + 0.5) * ymax / M) ** 2)) for i in range(M)) * (2.0 * ymax / M)
    return abs(s1 - s2)

def S_q_exact(d, X):
    T1 = 0.5 * d * X
    return sum((d / (2.0 * l)) * D_exact_raw(l * T1) for l in range(1, 30)) + (d / 2.0) * (4.0 / (T1 * T1)) / (2.0 * 30 * 30)

# Exact 1-variable function h_ex(t) := 8 t^2 |t^2 - 3| / (t^2 + 1)^3!
# Remember for small th = 2*pi*v/q, F(gam, th) / (r c) = |h_ex(t)| / (4*pi^2 * d * m * v^2) where t = m * X / (2*pi*v)!
# Check exact worst-case Set I (exact F(m*q*kap, 2*pi*v/q) summed over m>=1, 1<=v<=q/2) vs asymptotic h_ex(m*X/(2*pi*v))!
def S_SetI_exact_F(d, q, X):
    kap = X / (q * q)
    s = 0.0
    mmax = min(400, int(18.0 / (q * kap)) + 15)
    for m in range(1, mmax):
        gam = m * q * kap
        if gam > 35: break
        eg = math.exp(-gam)
        for v in range(1, q):
            th = 2.0 * math.pi * v / q
            z = eg * complex(math.cos(th), math.sin(th))
            F = 4.0 * gam * abs((z * (1.0 + z) / ((1.0 - z) ** 3)).real)
            # r >= m*q -> 1/(r*c) <= 1/(m*q * d*q) = 1/(m * d * q^2)
            s += F / (m * d * q * q)
    return s

def S_SetI_asymp(d, X):
    # (2/d) * sum_{m>=1, v>=1} |h_ex(m*X/(2*pi*v))| / (4*pi^2 * m * v^2)
    s = 0.0
    for v in range(1, 80):
        for m in range(1, 120):
            t = m * X / (2.0 * math.pi * v)
            ht = 8.0 * (t * t) * abs(t * t - 3.0) / ((t * t + 1.0) ** 3)
            s += 2.0 * ht / (4.0 * math.pi * math.pi * m * v * v)
    return s / d

for d, X in [(1, 0.85), (1, 0.90), (1, 1.0), (1, 1.5), (1, 2.5), (1, 4.0), (1, 8.0),
             (2, 2.0), (2, 2.2), (2, 2.5), (2, 4.0), (5, 2.0), (5, 2.5)]:
    sq = S_q_exact(d, X)
    sa = S_SetI_asymp(d, X)
    # Compare with exact finite q = 5, 20:
    sf5 = S_SetI_exact_F(d, 5, X)
    sf20 = S_SetI_exact_F(d, 20, X)
    print(f"  d={d} X={X:4.2f}: S_q={sq:.4f} | SetI_asymp={sa:.4f} (q=5:{sf5:.4f}, q=20:{sf20:.4f}) | SUM={sq+sa:.4f}")

print("\n=== 3. EXACT UNIVERSAL BOUND ON SET II BY CONVERGENT LADDER (PURE CLOSED-FORM FORMULA!) ===")
# Let's check the EXACT worst possible value of Set II = sum_{1 <= r < q} |b_r^main|/r
# over ALL continued fraction sequences (q_0, q_1, ..., q_s = q) satisfying:
#   for every k < s:  c/(q_k D_k^2) < M_G3   OR   q_k * kap > eps_G3 !
# Wait: FIRST verify our analytical bound on Set II for ANY CF ladder!
# Look at what makes Set II easy to bound analytically:
# For EVERY level k in {0, ..., s-1}:
#   1) Convergent multiples m*q_k (1 <= m <= a_{k+1}):
#      Sum <= zeta(3) * B(q_k, D_k), where B(q_k, D_k) = min(gA/(q_k c s_k^2), gB/(q_k^2 c kap s_k), gC q_k kap^2/(c s_k^4), bE)!
#      WAIT: Notice t_k = (q_k kap) / (2 sin(pi D_k / c))!
#      Even better: write s_k = sin(pi D_k / c), gamma_k = q_k kap!
#      Check the ratio t_{k+1} / t_k!
#      For k <= s-2: q_{k+1} >= q_k + q_{k-1} and D_k >= D_{k+1} + D_{k+2} -> t_{k+1} / t_k = (q_{k+1}/q_k) * (s_k/s_{k+1})!
#      Wait: can we bound sum_{k=0}^{s-1} zeta(3) * B(q_k, D_k) PURELY AS A FUNCTION OF (M_G3, eps_G3, X)?!
#      CHECK! Let k* be the convergent maximizing B(q_k, D_k)!
#      Wait: at k = k*, EITHER q_{k*} kap <= eps_G3 (so c/(q_{k*} D_{k*}^2) < M_G3 -> B <= 0.0768 * M_G3!)
#                      OR q_{k*} kap > eps_G3 (so B <= 0.349 / (q_{k*} D_{k*} kap) <= 0.349 / (D_{k*} eps_G3) <= 0.349 / eps_G3!)
#      HOLD ON! What about the OTHER convergents k != k*?!
#      Look! If k = k* - 1 has t_{k*-1} <= t_{k*} / 2.618, COULD k* - 1 ALSO have M_{k*-1} = M_G3 AND t_{k*-1} near the peak?!
#      YES!! Wait: if t_{k*} = 1.1 and t_{k*-1} = 0.42 (Fibonacci ratio 2.618!), BOTH t = 0.42 and t = 1.1 are near the peak t = 0.63!
#      TEST THAT EXACT WORST-CASE CF!! Can TWO consecutive convergents BOTH have c/(q_k D_k^2) ~ M_G3?!
#      CHECK THE EXACT CF IDENTITY q_{k+1} D_k + q_k D_{k+1} = c!!
#      Multiply (c / (q_k D_k^2)) * (c / (q_{k+1} D_{k+1}^2))?! NO, look at M_k = c / (q_k D_k^2):
#        Since q_{k+1} D_k < c -> D_k < c / q_{k+1},
#        AND q_k D_k < c * q_k / q_{k+1} = c / a_{k+1}!
#        Check: M_{k-1} = c / (q_{k-1} D_{k-1}^2) <= c / (q_{k-1} * (c / (2 q_k))^2) = 4 q_k^2 / (c q_{k-1})!
#        Meanwhile M_k = c / (q_k D_k^2) >= c / (q_k * (c / q_{k+1})^2) = q_{k+1}^2 / (c q_k)!
#        WAIT A SECOND!! LOOK AT THE PRODUCT OF M_{k-1} AND M_k? No, LOOK AT M_{k-1} / M_k!!
#          D_{k-1} = a_{k+1} D_k + D_{k+1} >= D_k + D_{k+1}, AND q_k = a_k q_{k-1} + q_{k-2}!
#          Compute (q_{k-1} D_{k-1}^2) / (q_k D_k^2) = (q_{k-1}/q_k) * (D_{k-1}/D_k)^2!
#          In Fibonacci (a_i = 1): q_{k-1}/q_k = 1/phi, D_{k-1}/D_k = phi -> ratio = phi = 1.618 -> M_{k-1} = M_k / phi!!
#          IN GENERAL: D_{k-1} >= a_{k+1} D_k and q_k <= (a_k + 1) q_{k-1} -> M_{k-1} / M_k = (q_k / q_{k-1}) / (D_{k-1} / D_k)^2!
#          HOLD ON! If a_k is HUGE (say a_k = 10) and a_{k+1} = 1, then q_k / q_{k-1} = 10 while D_{k-1}/D_k = 1.5, so M_{k-1} > M_k!
#          SO WHAT ABOUT THE SUM sum_{k=0}^{s-1} c/(q_k D_k^2) * |h(t_k)|?!
#          Let's maximize sum_{k=0}^{s-1} min(M_k, 0.349/...) * |h(t_k)| over ALL CF partial-quotient sequences [a_1, a_2, ..., a_s]!

def test_worst_CF_setII(M_G3, eps_G3):
    # Maximize Set II over a huge search of CF partial quotients [a1, a2, ..., a_s] and kap!
    worst = (0.0, None)
    # Test 2-level, 3-level, 4-level, 5-level CFs (including Fibonacci [1,1,1,1,1,2] and double-spike [1, M, 1, M, 2]!)
    cf_list = []
    for a1 in [1, 2, 3, 5, 10, 20]:
        for a2 in [1, 2, 3, 5, 10, 20]:
            cf_list.append([a1, a2 + 1])
            for a3 in [1, 2, 4, 8, 15]:
                cf_list.append([a1, a2, a3 + 1])
    cf_list += [[1,1,1,1,1,2], [1,1,1,1,1,1,2], [2,2,2,2,2], [1,5,1,5,2], [1,8,1,8,2], [1,10,1,10,2], [2,8,2,8,2]]
    for cf in cf_list:
        # Build p/q from cf = [a1, a2, ..., as] (a_s >= 2)
        p, q = 1, cf[-1]
        for a in reversed(cf[:-1]):
            p, q = q, a * q + p
        # Try d = 1, 2, 3, 5, 10 (only if c = d*q <= 5000)
        for d in [1, 2, 5]:
            c = d * q
            if c < 6 or c > 2000: continue
            j = c + d * p
            # Build convergents q_0 = 1, q_1, ..., q_s = q
            qs = [1]; p0, p1, q0, q1 = 0, 1, 1, 0
            a_tmp, b_tmp = p, q
            while b_tmp != 0:
                ai = a_tmp // b_tmp
                p2, q2 = ai * p1 + p0, ai * q1 + q0
                if q2 > qs[-1]: qs.append(q2)
                p0, p1, q0, q1 = p1, p2, q1, q2
                a_tmp, b_tmp = b_tmp, a_tmp % b_tmp
            Dks = [min((qk * j) % c, c - ((qk * j) % c)) for qk in qs]
            # Candidate kaps: peak of each convergent + eps_G3/qk
            kaps = set()
            for qk, Dk in zip(qs[:-1], Dks[:-1]):
                if Dk > 0:
                    kaps.add(2.0 * math.pi * Dk / (c * qk) * 0.63)
                    kaps.add(1.01 * eps_G3 / qk)
                    kaps.add(0.99 * eps_G3 / qk)
            for kap in kaps:
                if kap >= 2.5 or c * kap * kap >= 42.0 or c * q * kap <= 0.85: continue
                if d >= 2 and q * q * kap <= 2.2: continue
                # Check G3!
                if any(qk * kap <= eps_G3 and Dk > 0 and c >= M_G3 * qk * Dk * Dk for qk, Dk in zip(qs[:-1], Dks[:-1])):
                    continue
                # Compute EXACT sum_{k=0}^{s-1} zeta(3) * |b_{qk}^main|/qk AND full sum_{1<=r<q} |b_r^main|/r!
                S_conv_ex = 0.0
                for qk, Dk in zip(qs[:-1], Dks[:-1]):
                    gam = qk * kap; th = 2.0 * math.pi * Dk / c
                    z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
                    S_conv_ex += 1.20206 * (4.0 * kap / c) * abs((z * (1.0 + z) / ((1.0 - z) ** 3)).real)
                if S_conv_ex > worst[0]:
                    worst = (S_conv_ex, (cf, d, c, round(kap, 4), [(qk, Dk, round(c/(qk*Dk*Dk), 1)) for qk, Dk in zip(qs[:-1], Dks[:-1])]))
    print(f"  Worst CONVERGENT-LADDER sum (M={M_G3}, eps={eps_G3}): {worst[0]:.4f} at {worst[1]}")

for M_G3, eps_G3 in [(12, 0.35), (10, 0.40), (10, 0.45), (8, 0.45), (8, 0.50)]:
    test_worst_CF_setII(M_G3, eps_G3)
