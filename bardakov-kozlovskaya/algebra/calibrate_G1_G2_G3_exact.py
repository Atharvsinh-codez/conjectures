"""
Calibrate the exact G1 <-> G2 <-> G3 thresholds and test a 100% CLOSED-FORM
Continued-Fraction upper bound on sum_{r>=1} |b_r| / r!
Also check whether G3 (D != 0 Major Arcs) stays super easy at looser (M_G3, eps_G3)!
"""
import math

def b_r_gauss(r, c, j, kap):
    th = 2.0 * math.pi * ((r * j) % c) / c
    B = 0.0
    for m in range(1, c):
        ex = r * kap * m * (1.0 - m / c)
        if ex < 40.0:
            B += math.exp(-ex) * math.cos(m * th)
    A = 0.0
    m = 1
    while True:
        ex = r * kap * m * (1.0 + m / c)
        if ex > 40.0:
            break
        A += 2.0 * math.exp(-ex) * math.cos(m * th)
        m += 1
    return B - A

def b_r_main(r, c, j, kap):
    gam = r * kap
    if gam > 40.0:
        return 0.0
    th = 2.0 * math.pi * ((r * j) % c) / c
    z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
    return 4.0 * r * (kap / c) * (z * (1.0 + z) / ((1.0 - z) ** 3)).real

print("=== 1. CHECK EXACT WORST S = sum |b_r|/r AT r0=1, D=1, r0*kap = 0.25..0.55 ===")
for eps in [0.25, 0.30, 0.35, 0.40, 0.45, 0.50]:
    worst = (0.0, None)
    for c in [12, 15, 18, 20, 22, 25, 28, 32, 40, 60, 100]:
        j = c + 1
        kap = eps + 1e-4
        if c * kap * kap >= 42.0:
            continue
        S = sum(abs(b_r_gauss(r, c, j, kap)) / r for r in range(1, int(25.0 / kap) + 5))
        # Also check max |Phi(phi)| directly!
        brs = [b_r_gauss(r, c, j, kap) / r for r in range(1, min(200, int(25.0 / kap) + 5))]
        mx_phi = max(abs(sum(brs[r - 1] * math.sin(r * phi) for r in range(1, len(brs) + 1)))
                     for phi in [math.pi * k / 200.0 for k in range(1, 200)])
        if S > worst[0]:
            worst = (S, (c, kap, mx_phi))
    print(f"  eps_G3={eps:.2f}: worst sum|b_r|/r = {worst[0]:.4f} (max|Phi|={worst[1][2]:.4f}) at c={worst[1][0]}, kap={worst[1][1]:.4f}")

print("\n=== 2. CHECK EXACT 4-REGIME CLOSED-FORM BOUND ON |b_r^main|/r ===")
# Check exact constants C_A, C_B, C_C, C_D for |b_r^main|/r = (4*gam/(r*c)) * |Re z(1+z)/(1-z)^3|
# over ALL (gam > 0, th in (0, pi])!
# Define phi_th := 2*sin(th/2) in (0, 2], or delta = th/(2pi) in (0, 1/2], so D_r = c*delta!
# Let F(gam, th) := 4*gam * |Re z(1+z)/(1-z)^3| (so |b_r^main|/r = F(gam, th) / (r*c)).
# Check:
#   (A) Peak:      F(gam, th) * sin^2(th/2)  <=  K_A  (theory small-th: 3.03224 / 4 = 0.75806)
#   (B) Mid-t:     F(gam, th) * gam * sin(th/2) <= K_B (theory small-th: 2.1558 / 2 = 1.0779)
#   (C) Small-t:   F(gam, th) * sin^4(th/2) / gam^2 <= K_C (theory small-th: 24 / 16 = 1.500)
#   (D) Large-t:   F(gam, th) * gam^2 <= K_D (theory small-gam: 8.000)
#   (E) Exp-decay: F(gam, th) * sin^2(th/2) * exp(gam/2) <= K_E? Or exact F(gam, th) * (1 - e^{-gam})^3 / (4*gam*e^{-gam}*(1+e^{-gam})) <= 1!
K_A = K_B = K_C = K_D = 0.0
for ig in range(1, 400):
    gam = 0.002 * (1.025 ** ig)
    for it in range(1, 200):
        th = math.pi * it / 200.0
        z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
        F = 4.0 * gam * abs((z * (1.0 + z) / ((1.0 - z) ** 3)).real)
        s = math.sin(th / 2.0)
        K_A = max(K_A, F * s * s)
        K_B = max(K_B, F * gam * s)
        K_C = max(K_C, F * (s ** 4) / (gam * gam))
        K_D = max(K_D, F * gam * gam)
print(f"  EXACT GLOBAL CONSTANTS for F(gam,th) = 4*gam*|Re z(1+z)/(1-z)^3|:")
print(f"    K_A (F * sin^2)       = {K_A:.5f}  (-> C_A = K_A / (4*pi^2 * (2/pi)^2)? Wait: sin(th/2) >= 2*delta = 2 D_r/c -> K_A/4 = {K_A/4:.5f})")
print(f"    K_B (F * gam * sin)   = {K_B:.5f}  (-> K_B/2 = {K_B/2:.5f})")
print(f"    K_C (F * sin^4/gam^2) = {K_C:.5f}  (-> K_C/16 = {K_C/16:.5f})")
print(f"    K_D (F * gam^2)       = {K_D:.5f}")

# Wait! Check F(gam, th) vs sin(th/2) >= (2/pi)*(th/2) = 2*delta!
# WAIT: for small th, sin(th/2) = pi*delta, NOT 2*delta! Losing pi/2 = 1.57x costs 2.47x in sin^2!
# Instead check F(gam, th) directly in terms of delta = th/(2*pi) in (0, 1/2]!
C_A = C_B = C_C = 0.0
for ig in range(1, 400):
    gam = 0.002 * (1.025 ** ig)
    for it in range(1, 200):
        delta = 0.5 * it / 200.0
        th = 2.0 * math.pi * delta
        z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
        F = 4.0 * gam * abs((z * (1.0 + z) / ((1.0 - z) ** 3)).real)
        C_A = max(C_A, F * delta * delta)
        C_B = max(C_B, F * gam * delta)
        C_C = max(C_C, F * (delta ** 4) / (gam * gam))
print(f"  EXACT GLOBAL CONSTANTS IN delta = D_r/c in (0, 1/2]:")
print(f"    C_A (F * delta^2)       = {C_A:.5f}  (|b_r^main|/r <= C_A * c / (r D_r^2))")
print(f"    C_B (F * gam * delta)   = {C_B:.5f}  (|b_r^main|/r <= C_B / (r^2 D_r kap))")
print(f"    C_C (F * delta^4/gam^2) = {C_C:.5f}  (|b_r^main|/r <= C_C * c^3 r kap^2 / D_r^4)")

print("\n=== 3. TEST CLOSED-FORM CONTINUED-FRACTION UPPER BOUND ON sum_{q nmid r} |b_r^main|/r ===")
# Wait: for ANY (c, j, kap) in G1, let 1 = q_0 <= q_1 < ... < q_s = q be the CF convergents of j/c!
# Let's check how a PURELY ANALYTIC CF SUM BOUND behaves!
# Idea: For EVERY r >= 1 with q \nmid r, let k in {0, ..., s-1} be the unique index with q_k <= r < q_{k+1}
# (or if r >= q, reduce r mod q: every block [m*q, (m+1)*q - 1] has {r*j/c} = {1/q, 2/q, ..., (q-1)/q} EXACTLY!).
# WAIT A SECOND!! STOP!!
# For r >= q: EVERY block of q consecutive integers has EXACTLY TWO points at distance v/q (1 <= v <= q/2)!
# And for r < q: split into blocks [q_k, q_{k+1} - 1] for k = 0, ..., s-1!
# Wait: inside [q_k, q_{k+1} - 1], write r = m*q_k + r' (1 <= m <= a_{k+1}, 0 <= r' < q_k):
#   - r' = 0 (multiples m*q_k): D_r = m * D_k!
#   - r' > 0: any two distinct r in a sub-block of length q_k have |r1 - r2| < q_k -> ||(r1-r2)alpha|| >= delta_{k-1}
#     AND every r in [q_k, q_{k+1}-1] has ||r*alpha|| >= delta_k!
# Even simpler: evaluate BOTH exact sum_{q \nmid r} |b_r^main|/r AND our analytic CF bound across 5000 pairs!

def cf_convergents(num, den):
    # Returns list of (p_k, q_k) with q_0=1, ..., q_s = den/gcd(num,den)
    g = math.gcd(num, den)
    a, b = num // g, den // g
    # fractional part a % b / b
    a = a % b
    if a == 0:
        return [(0, 1)]
    convs = [(0, 1)]
    p_prev, q_prev = 1, 0
    p_curr, q_curr = 0, 1
    x, y = b, a
    while y != 0:
        q_part = x // y
        p_next = q_part * p_curr + p_prev
        q_next = q_part * q_curr + q_prev
        convs.append((p_next, q_next))
        p_prev, q_prev = p_curr, q_curr
        p_curr, q_curr = p_next, q_next
        x, y = y, x % y
    # Wait: p/q above are convergents of a/b starting with p0/q0=0/1, p1/q1=1/a1, ...
    # Actually let's compute standard CF of a/b properly:
    res = []
    p0, p1 = 0, 1
    q0, q1 = 1, 0
    num0, den0 = a, b
    while den0 != 0:
        ai = num0 // den0
        p2 = ai * p1 + p0
        q2 = ai * q1 + q0
        if q2 >= 1:
            res.append((p2, q2))
        p0, p1 = p1, p2
        q0, q1 = q1, q2
        num0, den0 = den0, num0 % den0
    return res

# Test the CF identity q_{k+1} D_k + q_k D_{k+1} = c and sum of 1/(r D_r^2) vs peak!
# Wait: for a given convergent ladder q_0 < q_1 < ... < q_s = q, how big can
#   S_CF_exact := sum_{k=0}^{s-1} (sum_{m=1}^{a_{k+1}} term(m*q_k, m*D_k)) + non-multiples
# be when ALL convergents satisfy: (c/(q_k D_k^2) < M_G3 OR q_k*kap > eps_G3)?!
def exact_term_bound(r, D, c, kap):
    # Exact 4-way min for |b_r^main|/r at (r, D)
    gam = r * kap
    delta = D / c
    if delta == 0:
        return 0.0
    bA = C_A * c / (r * D * D)
    bB = C_B / (r * r * D * kap)
    bC = C_C * (c ** 3) * r * kap * kap / (D ** 4)
    bD = 8.0 / (r ** 3 * c * kap * kap)
    bE = (4.0 * kap / c) * math.exp(-gam) * (1.0 + math.exp(-gam)) / ((1.0 - math.exp(-gam)) ** 3) if gam < 40 else 0.0
    return min(bA, bB, bC, bD, bE)

def check_4way_sum(M_G3, eps_G3, K_G2):
    worst_nonq = (0.0, None)
    worst_tot = (0.0, None)
    # Test tricky rationals: Fibonacci convergents, small c, medium c, large c with 2 simultaneous near-G3 convergents!
    test_cases = []
    # 1) All small c = 6..35, various j
    for c in range(6, 36):
        for j in range(c + 1, c + c // 2 + 2):
            test_cases.append((c, j))
    # 2) Larger c with crafted CF [0; a1, a2, a3]
    for a1 in [2, 3, 4, 5, 8]:
        for a2 in [1, 2, 3, 5, 10]:
            for a3 in [1, 2, 3, 5]:
                # rational p/q = 1/(a1 + 1/(a2 + 1/a3))
                q = a1 * (a2 * a3 + 1) + a3
                p = a2 * a3 + 1
                for d in [1, 2, 3, 5]:
                    test_cases.append((d * q, d * (q + p)))
    for c, j in test_cases:
        d = math.gcd(c, j)
        q = c // d
        # Test kap values around every convergent scale 1/q_k and around eps_G3 / q_k!
        convs = cf_convergents(j, c)
        kaps = set([0.005, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 0.80, 1.20])
        for _, qk in convs:
            for mult in [0.5, 1.01, 1.5, 2.5]:
                kaps.add(mult * eps_G3 / qk)
            Dk = min((qk * j) % c, c - ((qk * j) % c))
            if Dk > 0:
                kaps.add(2.0 * math.pi * Dk / (c * qk) * 0.63)
        if d >= 2:
            kaps.add((K_G2 + 0.01) / (q * q))
        for kap in kaps:
            if kap >= 2.5 or kap <= 1e-6 or c * kap * kap >= 42.0 or c * q * kap <= 0.85:
                continue
            if d >= 2 and q * q * kap <= K_G2:
                continue
            # Check G3:
            in_G3 = False
            for _, qk in convs:
                if qk == q:
                    continue
                Dk = min((qk * j) % c, c - ((qk * j) % c))
                if Dk > 0 and qk * kap <= eps_G3 and c >= M_G3 * qk * Dk * Dk:
                    in_G3 = True
                    break
            if in_G3:
                continue
            # Also check ALL r0 <= eps_G3/kap just in case a non-convergent triggers G3 (convergents suffice for c/(r0 D^2) >= 2!)
            S_nonq = 0.0
            rmax = min(3000, int(18.0 / kap) + 20)
            for r in range(1, rmax):
                if r % q == 0:
                    continue
                Dr = min((r * j) % c, c - ((r * j) % c))
                S_nonq += exact_term_bound(r, Dr, c, kap)
            # q|r exact bound (outside G2):
            T1 = 0.5 * c * q * kap
            if d == 1:
                S_q = 0.6138
            else:
                S_q = 1.1477 / (q * q * kap)
            if S_nonq > worst_nonq[0]:
                worst_nonq = (S_nonq, (c, j, d, q, kap, S_q))
            if S_nonq + S_q > worst_tot[0]:
                worst_tot = (S_nonq + S_q, (c, j, d, q, kap, S_nonq, S_q))
    print(f"  (M_G3={M_G3:2d}, eps_G3={eps_G3:.2f}, K_G2={K_G2:.1f}): worst 4-way S_nonq = {worst_nonq[0]:.4f}, worst S_nonq+S_q = {worst_tot[0]:.4f} at {worst_tot[1]}")

for M_G3, eps_G3, K_G2 in [
    (15, 0.25, 2.0),
    (12, 0.35, 2.0),
    (10, 0.40, 2.0),
    (10, 0.45, 2.2),
    (8,  0.50, 2.5),
]:
    check_4way_sum(M_G3, eps_G3, K_G2)
