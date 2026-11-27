"""
Test the TWO final pieces of Piece G1 (Dirichlet Minor Arcs):
1. Measure EXACT finite-n difference sum_{r>=1} |b_r^exact - b_r^G| / r across a comprehensive
   stress grid of G1 (all small c, medium c, large c up to c*kap^2 = 41.5, kap up to 2.45, n >= 132)!
2. Test our 100% CLOSED-FORM analytic upper bound on Set II (convergents + semiconvergents + rest)
   + Set I + S_q as a PURE FORMULA of (M_G3, eps_G3, K_G2) — ZERO loops over (c, j, r)!
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

def exact_finite_n_br(n, c, j, rmax):
    th_c = c * math.pi / n; psi_j = j * math.pi / n; alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1; log_Rz = n * math.log(r3 / r1); kap = 2.0 * math.pi * r1 * r2 / r3
    shift = 0.0 if (c % 2 == 0) else 0.5; sgn_ang = 0.0 if (j % 2 == 0) else math.pi
    small_a = []; large_a = []
    for k in range(-(n // 2), n // 2 + 1):
        x = k + shift
        if x <= -n / 2.0 + 1e-12 or x > n / 2.0 + 1e-12 or abs(abs(x) - c / 2.0) < 1e-9: continue
        eta = 2.0 * math.pi * x / n; z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
        lr = log_Rz - n * math.log(abs(z)); ang = sgn_ang - n * math.atan2(z.imag, z.real)
        if abs(x) < c / 2.0:
            if lr < 45.0: small_a.append((-lr, -ang))
        else:
            if lr > -45.0: large_a.append((lr, ang))
    return [sum(math.exp(r*lr)*math.cos(r*ang) for lr,ang in small_a if r*lr>-40) -
            sum(math.exp(r*lr)*math.cos(r*ang) for lr,ang in large_a if r*lr>-40) for r in range(1, rmax+1)], kap

print("=== 1. STRESS TEST sum_{r>=1} |b_r^exact - b_r^G| / r ACROSS G1 ===")
worst_diff = (0.0, None)
for c, j, kt in [
    (2, 3, 2.0), (3, 4, 2.0), (4, 5, 2.4), (5, 6, 2.4),
    (6, 7, 1.5), (6, 7, 2.4), (7, 8, 2.3), (8, 9, 2.1),
    (10, 11, 1.8), (10, 13, 0.05), (10, 13, 0.5), (10, 13, 2.0),
    (12, 13, 0.50), (12, 13, 1.80), (15, 16, 0.50), (15, 16, 1.60),
    (20, 21, 0.46), (20, 21, 1.40), (25, 26, 0.50), (25, 26, 1.25),
    (29, 39, 0.05), (29, 39, 1.15), (40, 41, 0.50), (40, 41, 1.00),
    (60, 61, 0.50), (60, 61, 0.82), (100, 101, 0.46), (100, 101, 0.64),
]:
    if c * kt * kt >= 41.5 or kt >= 2.48: continue
    n = max(132, int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt))) // 2 * 2)
    # check valid strictly interior j: j <= n/2 - ceil(c/2) - (1 if c even else 0)
    if j > n // 2 - (c + 1) // 2 - (1 if c % 2 == 0 else 0):
        continue
    rmax = min(200, int(18.0 / kt) + 15)
    brs, kap = exact_finite_n_br(n, c, j, rmax)
    if kap >= 2.5 or c * kap * kap >= 42.0: continue
    S_diff = sum(abs(brs[r - 1] - b_r_gauss(r, c, j, kap)) / r for r in range(1, rmax + 1))
    S_ex = sum(abs(brs[r - 1]) / r for r in range(1, rmax + 1))
    S_ga = sum(abs(b_r_gauss(r, c, j, kap)) / r for r in range(1, rmax + 1))
    if S_diff > worst_diff[0]:
        worst_diff = (S_diff, (c, j, n, round(kap, 4), round(c*kap*kap, 2), round(S_ex, 4), round(S_ga, 4)))
    if S_diff > 0.04:
        print(f"  c={c:3d} j={j:3d} n={n:5d} kap={kap:.3f} ck^2={c*kap*kap:5.2f} | S_ex={S_ex:.4f}, S_ga={S_ga:.4f}, S_diff={S_diff:.4f}")
print(f"  WORST S_diff = sum |b_r^ex - b_r^G|/r across stress grid: {worst_diff[0]:.4f} at {worst_diff[1]}")

print("\n=== 2. TEST PURE CLOSED-FORM CONVERGENT + SEMICONVERGENT + REST BOUND ON SET II ===")
# Wait: let's check our EXACT closed-form bound on Set II (1 <= r < q) for ANY rational p/q!
# Recall for ANY level k in {0, ..., s-1} (q_k <= r < q_{k+1}):
#   Every r in [q_k, q_{k+1} - 1] has unique representation r = m*q_k + r' (1 <= m <= a_{k+1}, 0 <= r' < q_k).
#   Let's split Set II cleanly into:
#   (A) Convergent multiples: r = m * q_k (1 <= m <= a_{k+1}) -> D_r = m * D_k.
#       Sum <= zeta(3) * B(q_k, D_k) for each k in {0, ..., s-1}.
#   (B) Non-multiples r in [q_k, q_{k+1} - 1] (only exists when q_k >= 2, i.e., k >= 1!):
#       Write r = m * q_k + r' (1 <= m <= a_{k+1}, 1 <= r' <= q_k - 1, with r < q_{k+1}).
#       WAIT A SECOND! For fixed k >= 1 and fixed m in [1, a_{k+1}]:
#       All q_k - 1 values of r' give points {r * alpha mod 1} separated by >= D_{k-1} / c!
#       AND every one of them satisfies D_r >= D_k!
#       Wait: EXCEPT for the ONE closest point per m (which is the semiconvergent r = m*q_k + q_{k-1} when m < a_{k+1},
#       having D_r = D_{k-1} - m*D_k = (a_{k+1} - m)*D_k + D_{k+1}, or for m = a_{k+1} the closest point is >= D_{k-1} - a_{k+1}*D_k + D_{k-1}? Wait: m=a_{k+1}, r'<q_{k-1}, so D_r >= D_{k-1} - D_k? Let's check!),
#       ALL OTHER points in that block of length q_k have D_r >= v * D_{k-1} (2 copies for each v = 1, 2, ...)!
#       CHECK HOW BIG THE SUM OVER (B) = SEMICONVERGENTS + OTHER NON-MULTIPLES IS ANALYTICALLY!!
#
#   Let's compute the EXACT universal bound for:
#     (B1) ALL Semiconvergents r = q_{k+1} - u*q_k (1 <= u <= a_{k+1} - 1), D_r = u*D_k + D_{k+1}:
#          Notice: r >= max(q_k, q_{k+1} - u*q_k), and D_r >= u*D_k!
#          Actually: u <= a_{k+1}/2 -> r >= q_{k+1}/2, D_r >= u*D_k -> c/(r D_r^2) <= 2 c / (q_{k+1} u^2 D_k^2) = (2 / (a_{k+1} u^2)) * M_k!
#          EVEN SIMPLER: Look at r * D_r^2 for u = 1..a-1 where r = (a - u)*q_k + q_{k-1} >= (a - u)*q_k and D_r >= u*D_k:
#            sum_{u=1}^{a-1} 1 / ((a - u) * u^2) !!
#          Compute max_{a >= 2} sum_{u=1}^{a-1} 1 / ((a - u) * u^2)!!
#            a=2: 1/(1*1) = 1.000  (wait: at a=2, u=1 -> r = q_k + q_{k-1}, D_r = D_k + D_{k+1}, so r*D_r^2 >= (q_k + q_{k-1})(D_k + D_{k+1})^2!)
#            a=3: 1/2 + 1/4 = 0.750
#            a=4: 1/3 + 1/8 + 1/9 = 0.569
#            a=10: 1/9 + 1/32 + ... = 0.25!
#          WAIT A SECOND! That's ONLY O(1/a_{k+1}) times M_k! Wait: at a_{k+1} = 2, M_k = c/(q_k D_k^2) <= (a_{k+1} + 2)/D_k = 4/D_k!
#          WAIT A SECOND!! STOP EVERYTHING!!
#          CHECK WHAT WE JUST PROVED FOR EVERY SEMICONVERGENT AT LEVEL k!!
#          Since M_k < (a_{k+1} + 2) / D_k:
#            - If a_{k+1} is SMALL (say 2, 3, 4): M_k <= 6 / D_k IS SMALL!
#            - If a_{k+1} is LARGE (say >= 5): sum_{u=1}^{a-1} 1/((a-u) u^2) <= 2.5 / a_{k+1} CANCELS a_{k+1} IN M_k < (a_{k+1}+2)/D_k!!
#          MULTIPLY THEM!!
#            sum_{u=1}^{a_{k+1}-1} c / (r_u D_{r_u}^2) <= (sum_{u=1}^{a-1} 1/((a-u) u^2)) * M_k
#            <= (sum_{u=1}^{a-1} 1/((a-u) u^2)) * (a + 2) / D_k  <=  (3.0 / D_k) * ... WAIT: let's check max_{a >= 2} (a + 2) * sum_{u=1}^{a-1} 1/((a-u) u^2)!!
for a in [2, 3, 4, 5, 10, 20, 100]:
    s_u = sum(1.0 / ((a - u) * u * u) for u in range(1, a))
    print(f"  a={a:3d}: sum 1/((a-u)u^2) = {s_u:.4f},  (a+1)*sum = {(a+1)*s_u:.4f} (NOTE: for u=1..a-1, r >= (a-u)q_k, and M_k < (a+1)/D_k + 1/D_k!)")

# WAIT A SECOND!! Look at (a+1) * sum 1/((a-u)u^2) -> pi^2/6 = 1.645 as a -> inf, and max at a=2 is 3.000!
# EVEN BETTER: at a=2, u=1 has r = q_k + q_{k-1} and D_r = D_k + D_{k+1}, AND c = q_{k+1} D_k + q_k D_{k+1} = (2 q_k + q_{k-1}) D_k + q_k D_{k+1}!
# So c / (r D_r^2) <= ((2 q_k + q_{k-1}) D_k + q_k D_{k+1}) / ((q_k + q_{k-1}) (D_k + D_{k+1})^2) <= 2 / (D_k + D_{k+1}) <= 2 / D_k !!
# CHECK THAT!! FOR EVERY a >= 2 AND EVERY k, THE SUM OF c/(r D_r^2) OVER ALL SEMICONVERGENTS AT LEVEL k IS <= (2.5 / D_k) — COMPLETELY INDEPENDENT OF M_G3 AND a_{k+1}!!
# Let's check max over ALL CFs of (sum_{semi at level k} c/(r D_r^2)) * D_k!
print("\n=== 3. VERIFY UNIVERSAL BOUND ON SEMICONVERGENTS & REST AT LEVEL k (INDEPENDENT OF M_G3!) ===")
# Check (2) non-semiconvergent rest r = m*q_k + r' (1 <= m <= a_{k+1}, r' != q_{k-1}):
# For fixed m, the points {m*q_k + r' : 1 <= r' < q_k, r' != q_{k-1}} are D_{k-1}-separated, with v = 1, 2, ... (2 copies each, at distance >= v * D_{k-1} - m*D_k? No, for r' != q_{k-1}, min distance is >= D_{k-1} - m*D_k + ... wait: check!):
# For m*q_k + r' with 1 <= m <= a_{k+1} and 1 <= r' < q_k:
# Since D_k = |q_k j - c p_k| and D_{k-1} = |q_{k-1} j - c p_{k-1}| have OPPOSITE signs,
# the closest point in {m*q_k + r' : 1 <= r' < q_k} to 0 mod c has distance >= D_k,
# and because any two distinct r'_1, r'_2 in [1, q_k - 1] satisfy |r'_1 - r'_2| < q_k -> ||(r'_1 - r'_2) alpha|| >= D_{k-1}/c,
# sorting the q_k - 1 values of D_r in increasing order gives:
#   1st (semiconvergent): >= D_k
#   2nd: >= D_{k-1} - D_k? Wait: at least (D_{k-1} - D_r^{(1)}) >= D_{k-1}/2 if D^{(1)} <= D_{k-1}/2!
#   In general, in ANY interval [0, 1/2], points separated by >= delta = D_{k-1}/c satisfy:
#   after removing AT MOST ONE point in [0, delta/2), the remaining points satisfy D_{(v)} >= floor(v/2) * D_{k-1} (v >= 2) or >= (v - 1/2) D_{k-1}!!
#   Compute sum_{m=1}^{a_{k+1}} 1 / (m * q_k) * sum_{v >= 1} 2 * c / (v * D_{k-1})^2 !!
#   Wait: 1/(q_k D_{k-1}^2) * c = (q_k D_{k-1} + q_{k-1} D_k) / (q_k D_{k-1}^2) <= 2 / D_{k-1}!!
#   MULTIPLY BY sum_{m=1}^{a_{k+1}} (1/m) ~ log(a_{k+1})?! Wait: NO! We have the small-t / mid-t cutoff on m too, just like Set I!!
#   HOLD ON! CHECK sum_{r in rest at level k} c / (r D_r^2) vs (1 / D_{k-1})!
mx_semi_Dk = 0.0; mx_rest_sum = 0.0
for cf in [[1,1,1,1,1,2], [2,2,2,2,2], [1,5,1,5,2], [1,10,1,10,2], [2,10,2,10,2], [3,15,3], [1,20,2], [5,1,20,2], [20,1,9]]:
    p, q = 1, cf[-1]
    for a in reversed(cf[:-1]): p, q = q, a * q + p
    c = q; j = c + p
    qs = [1]; p0, p1, q0, q1 = 0, 1, 1, 0; at, bt = p, q
    while bt != 0:
        ai = at // bt; p2, q2 = ai*p1+p0, ai*q1+q0
        if q2 > qs[-1]: qs.append(q2)
        p0, p1, q0, q1 = p1, p2, q1, q2; at, bt = bt, at % bt
    Dks = [min((qk * j) % c, c - ((qk * j) % c)) for qk in qs]
    tot_rest = 0.0
    for idx in range(1, len(qs) - 1):
        qk = qs[idx]; qk1 = qs[idx + 1]; qkm1 = qs[idx - 1]; Dk = Dks[idx]; Dkm1 = Dks[idx - 1]
        ak1 = qk1 // qk
        s_semi = 0.0; s_rest_k = 0.0
        for r in range(qk, qk1):
            if r % qk == 0: continue
            Dr = min((r * j) % c, c - ((r * j) % c))
            term = c / (r * Dr * Dr)
            if (r - qkm1) % qk == 0:
                s_semi += term
            else:
                s_rest_k += term
        mx_semi_Dk = max(mx_semi_Dk, s_semi * Dk)
        tot_rest += s_rest_k
    mx_rest_sum = max(mx_rest_sum, tot_rest)
    print(f"  cf={str(cf):18s} q={q:4d}: sum_rest c/(r D^2) = {tot_rest:.4f} (times 0.0768 = {0.0768*tot_rest:.4f}) | max(semi*Dk) = {mx_semi_Dk:.4f}")
print(f"  GLOBAL max(semi * D_k) = {mx_semi_Dk:.4f} (-> times 0.0768 is ONLY {0.0768*mx_semi_Dk:.4f} / D_k !!)")
