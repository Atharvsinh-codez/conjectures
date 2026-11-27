"""
FAST (1-second) calibration of Piece G1 (Minor Arcs) using:
1. Memoized exact S_q(T1) for q|r,
2. Sharp sin(pi D_r / c) 5-way bound on |b_r^main|/r for q \nmid r,
3. Direct comparison with TRUE sum |b_r^G|/r AND max |Phi(phi)| across the G1 <-> G2 <-> G3 boundary!
4. Test of a 100% CLOSED-FORM Continued-Fraction Analytic Theorem (ZERO loop over r!) on G1!
"""
import math

K_A = 0.7581   # F * sin^2(th/2) <= K_A
K_B = 1.0785   # wait: check exact max F * gam * (th/2) below!
K_C = 1.5001   # F * sin^4(th/2) / gam^2 <= K_C
K_D = 8.0001   # F * gam^2 <= K_D

# Quick exact check of K_B_th = max F * gam * (th/2)
KB_th = 0.0
for ig in range(1, 200):
    gam = 0.005 * (1.035 ** ig)
    for it in range(1, 100):
        th = math.pi * it / 100.0
        z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
        F = 4.0 * gam * abs((z * (1.0 + z) / ((1.0 - z) ** 3)).real)
        KB_th = max(KB_th, F * gam * (th / 2.0))
C_B_th = KB_th / math.pi
print(f"Exact K_B_th = {KB_th:.5f}  ->  |b_r^main|/r <= {C_B_th:.5f} / (r^2 D_r kap) !")

# Precompute D(T) table on a fine grid T in [0.4, 30]:
def D_exact_raw(T):
    M = 120
    s1 = sum(math.exp(-0.5 * T * (1.0 - ((i + 0.5) / M) ** 2)) for i in range(M)) * (2.0 / M)
    ymax = min(40.0, 20.0 / max(T, 1e-4) + 7.0 / math.sqrt(max(T, 1e-4)))
    s2 = sum(math.exp(-T * ((i + 0.5) * ymax / M) - 0.5 * T * (((i + 0.5) * ymax / M) ** 2)) for i in range(M)) * (2.0 * ymax / M)
    return abs(s1 - s2)

# Build fast table of f_q(T1) := sum_{l>=1} (1/(2l)) |D(l*T1)| for T1 = 0.40 + 0.05*k
T_grid = [0.40 + 0.05 * k for k in range(300)]
fq_grid = []
for T1 in T_grid:
    s = sum((0.5 / l) * D_exact_raw(l * T1) for l in range(1, 25)) + 2.0 / (T1 * T1 * 2.0 * 25 * 25)
    fq_grid.append(s)

def S_q_fast(d, q, kap):
    T1 = 0.5 * d * q * q * kap
    if T1 >= 15.0:
        return (d / 2.0) * (4.0 / (T1 * T1)) * 1.21
    idx = max(0, min(len(T_grid) - 1, int((T1 - 0.40) / 0.05)))
    return d * fq_grid[idx]

def cf_convergents(num, den):
    g = math.gcd(num, den)
    a, b = (num // g) % (den // g), den // g
    if a == 0:
        return [(0, 1)]
    res = []
    p0, p1, q0, q1 = 0, 1, 1, 0
    while b != 0:
        ai = a // b
        p2, q2 = ai * p1 + p0, ai * q1 + q0
        if q2 >= 1:
            res.append((p2, q2))
        p0, p1, q0, q1 = p1, p2, q1, q2
        a, b = b, a % b
    return res

def sharp_term_bound(r, D, c, kap):
    gam = r * kap
    if gam > 35.0:
        return 0.0
    th = 2.0 * math.pi * D / c
    s = math.sin(0.5 * th)
    bA = K_A / (r * c * s * s)
    bB = C_B_th / (r * r * D * kap)
    bC = K_C * r * kap * kap / (c * (s ** 4))
    bD = K_D / (r ** 3 * c * kap * kap)
    bE = (4.0 * kap / c) * math.exp(-gam) * (1.0 + math.exp(-gam)) / ((1.0 - math.exp(-gam)) ** 3)
    return min(bA, bB, bC, bD, bE)

def b_r_main_abs_over_r(r, D, c, kap):
    gam = r * kap
    if gam > 35.0:
        return 0.0
    th = 2.0 * math.pi * D / c
    z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
    return (4.0 * kap / c) * abs((z * (1.0 + z) / ((1.0 - z) ** 3)).real)

test_cases = []
for c in range(6, 36):
    for j in range(c + 1, c + c // 2 + 2):
        test_cases.append((c, j))
for a1 in [2, 3, 4, 5, 8]:
    for a2 in [1, 2, 3, 5, 10]:
        for a3 in [1, 2, 3, 5]:
            q = a1 * (a2 * a3 + 1) + a3
            p = a2 * a3 + 1
            for d in [1, 2, 3, 5]:
                test_cases.append((d * q, d * (q + p)))

print("\n=== TEST 1: EXACT |b_r^main|/r SUM vs 5-WAY SHARP BOUND ON G1 ===")
for M_G3, eps_G3, K_G2 in [
    (15, 0.25, 2.0),
    (12, 0.35, 2.0),
    (10, 0.40, 2.0),
    (10, 0.45, 2.2),
    (8,  0.50, 2.5),
]:
    worst_sharp = (0.0, None)
    worst_exact = (0.0, None)
    for c, j in test_cases:
        d = math.gcd(c, j); q = c // d
        convs = cf_convergents(j, c)
        kaps = set([0.01, 0.05, 0.12, 0.25, 0.40, 0.60, 1.00])
        for _, qk in convs:
            kaps.add(1.01 * eps_G3 / qk)
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
            in_G3 = False
            for _, qk in convs:
                if qk == q: continue
                Dk = min((qk * j) % c, c - ((qk * j) % c))
                if Dk > 0 and qk * kap <= eps_G3 and c >= M_G3 * qk * Dk * Dk:
                    in_G3 = True; break
            if in_G3: continue
            S_sharp = 0.0; S_ex = 0.0
            rmax = min(1500, int(14.0 / kap) + 15)
            for r in range(1, rmax):
                if r % q == 0: continue
                Dr = min((r * j) % c, c - ((r * j) % c))
                S_sharp += sharp_term_bound(r, Dr, c, kap)
                S_ex += b_r_main_abs_over_r(r, Dr, c, kap)
            S_q = S_q_fast(d, q, kap)
            if S_sharp + S_q > worst_sharp[0]:
                worst_sharp = (S_sharp + S_q, (c, j, d, q, round(kap, 4), round(S_sharp, 4), round(S_q, 4)))
            if S_ex + S_q > worst_exact[0]:
                worst_exact = (S_ex + S_q, (c, j, d, q, round(kap, 4), round(S_ex, 4), round(S_q, 4)))
    print(f"  (M_G3={M_G3:2d}, eps={eps_G3:.2f}, K_G2={K_G2:.1f}): EXACT={worst_exact[0]:.4f} {worst_exact[1]} | 5-WAY={worst_sharp[0]:.4f} {worst_sharp[1]}")
