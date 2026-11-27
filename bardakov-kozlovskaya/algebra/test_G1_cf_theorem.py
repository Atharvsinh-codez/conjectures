"""
1. Re-run the sharp 5-way bound on G1 using the CORRECT mid-t constant K_B = 1.0970 (on gam <= 1.5)
   and exact exponential bound on gam > 1.5!
2. Test our 100% CLOSED-FORM Continued-Fraction bound (ZERO loop over r!) on G1!
3. Check exact finite-n Hayman health of G3 (D != 0 Major Arcs) at the looser corner
   c/(r0 D^2) = 10..15, r0*kap = 0.25..0.45 so we know what is effortless for @analyst!
"""
import math

K_A = 0.7581    # F * sin^2(th/2) <= K_A (all gam > 0)
K_B = 1.0970    # F * gam * sin(th/2) <= K_B (all 0 < gam <= 1.5)
K_C = 1.5001    # F * sin^4(th/2) / gam^2 <= K_C (all gam > 0)
K_D = 8.0001    # F * gam^2 <= K_D (all gam > 0)

# Fast D(T) lookup
def D_exact_raw(T):
    M = 120
    s1 = sum(math.exp(-0.5 * T * (1.0 - ((i + 0.5) / M) ** 2)) for i in range(M)) * (2.0 / M)
    ymax = min(40.0, 20.0 / max(T, 1e-4) + 7.0 / math.sqrt(max(T, 1e-4)))
    s2 = sum(math.exp(-T * ((i + 0.5) * ymax / M) - 0.5 * T * (((i + 0.5) * ymax / M) ** 2)) for i in range(M)) * (2.0 * ymax / M)
    return abs(s1 - s2)

T_grid = [0.40 + 0.05 * k for k in range(300)]
fq_grid = [sum((0.5 / l) * D_exact_raw(l * T1) for l in range(1, 25)) + 2.0 / (T1 * T1 * 2.0 * 25 * 25) for T1 in T_grid]

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

def sharp_term_bound_fixed(r, D, c, kap):
    gam = r * kap
    if gam > 35.0:
        return 0.0
    th = 2.0 * math.pi * D / c
    s = math.sin(0.5 * th)
    bA = K_A / (r * c * s * s)
    bB = (K_B / (r * c * gam * s)) if gam <= 1.5 else (1.5040 / (r * c * gam * s))
    bC = K_C * r * kap * kap / (c * (s ** 4))
    bD = K_D / (r ** 3 * c * kap * kap)
    bE = (4.0 * kap / c) * math.exp(-gam) * (1.0 + math.exp(-gam)) / ((1.0 - math.exp(-gam)) ** 3)
    return min(bA, bB, bC, bD, bE)

# Also test the 1-LINE EXACT |Re| <= |cosh/sinh^3| bound + small-t factor!
def exact_envelope_bound(r, D, c, kap):
    gam = r * kap
    if gam > 35.0:
        return 0.0
    th = 2.0 * math.pi * D / c
    u = 0.5 * gam
    s = math.sin(0.5 * th)
    sh = math.sinh(u)
    # 1-line |Re| <= |cosh(w)/sinh^3(w)|:
    b1 = (gam * math.cosh(u) / ((sh * sh + s * s) ** 1.5)) / (r * c)
    # 1-line small-u bound (since at u=0, Re=0!):
    b2 = K_C * r * kap * kap / (c * (s ** 4))
    return min(b1, b2)

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

print("=== 1. FIXED 5-WAY BOUND AND 2-LINE ENVELOPE BOUND ON G1 ===")
for M_G3, eps_G3, K_G2 in [
    (15, 0.25, 2.0),
    (12, 0.35, 2.0),
    (10, 0.40, 2.0),
    (10, 0.45, 2.2),
    (8,  0.45, 2.2),
    (8,  0.50, 2.5),
]:
    w5 = (0.0, None); w2 = (0.0, None)
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
            S5 = 0.0; S2 = 0.0
            rmax = min(1500, int(14.0 / kap) + 15)
            for r in range(1, rmax):
                if r % q == 0: continue
                Dr = min((r * j) % c, c - ((r * j) % c))
                S5 += sharp_term_bound_fixed(r, Dr, c, kap)
                S2 += exact_envelope_bound(r, Dr, c, kap)
            S_q = S_q_fast(d, q, kap)
            if S5 + S_q > w5[0]: w5 = (S5 + S_q, (c, j, d, q, round(kap, 4), round(S5, 4), round(S_q, 4)))
            if S2 + S_q > w2[0]: w2 = (S2 + S_q, (c, j, d, q, round(kap, 4), round(S2, 4), round(S_q, 4)))
    print(f"  (M={M_G3:2d}, eps={eps_G3:.2f}, K_G2={K_G2:.1f}): 5-WAY={w5[0]:.4f} {w5[1]} | 2-LINE={w2[0]:.4f} {w2[1]}")

print("\n=== 2. CHECK EXACT FINITE-n G3 HEALTH AT LOOSER CORNERS (M=8..15, r0*kap=0.25..0.50) ===")
def exact_finite_n_roots_complex(n, c, j):
    th_c = c * math.pi / n; psi_j = j * math.pi / n; alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1; log_Rz = n * math.log(r3 / r1)
    roots = []; shift = 0.0 if (c % 2 == 0) else 0.5; sgn_ang = 0.0 if (j % 2 == 0) else math.pi
    k = 0
    while True:
        added = False
        klist = [0] if (k == 0 and c % 2 == 0) else ([0, -1] if (k == 0 and c % 2 == 1) else [k, -k - (1 if c % 2 == 1 else 0)])
        for kk in klist:
            eta = 2.0 * (kk + shift) * math.pi / n
            if abs(eta) > math.pi + 1e-12: continue
            z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
            log_r = log_Rz - n * math.log(abs(z))
            ang = sgn_ang - n * math.atan2(z.imag, z.real)
            if log_r > -50.0:
                roots.append((log_r, ang)); added = True
        if not added and k > c: break
        k += 1
        if k > n // 2 + 2: break
    return roots, 2.0 * math.pi * r1 * r2 / r3

import cmath
def check_G3_rat(n, c, j, r0):
    roots, kap = exact_finite_n_roots_complex(n, c, j)
    sgn_pref = 1.0 if ((j * (c + 1) + c) % 2 == 0) else -1.0
    # Scan lam in [0.7, 1.4] for best rat_L1!
    best_rat = -1.0; best_lam = 1.0
    for lam in [0.85, 1.0, 1.15, 1.30]:
        ll = math.log(lam)
        M = max(128, 32 * r0)
        vals = []
        max_lr = -1e9
        for k in range(M):
            phi = -math.pi + 2.0 * math.pi * (k + 0.5) / M
            la = -c * ll; ph = -c * phi
            for lr, ang in roots:
                lm = lr + ll; th = ang + phi
                if lm < -40.0: continue
                elif lm > 40.0: la += lm; ph += th + math.pi
                else:
                    r = math.exp(lm); z = complex(1 - r * math.cos(th), -r * math.sin(th))
                    la += math.log(abs(z)); ph += math.atan2(z.imag, z.real)
            if la > max_lr: max_lr = la
            vals.append((la, ph))
        num = sum(sgn_pref * math.exp(la - max_lr) * math.cos(ph) for la, ph in vals)
        den = sum(math.exp(la - max_lr) for la, ph in vals)
        if num / den > best_rat:
            best_rat = num / den; best_lam = lam
    return best_rat, best_lam, kap

for c, j, r0, D, r0k in [
    (10, 11, 1, 1, 0.45), (12, 13, 1, 1, 0.45), (15, 16, 1, 1, 0.35), (20, 21, 1, 1, 0.30),
    (20, 21, 1, 1, 0.45), (20, 31, 2, 2, 0.40), (25, 38, 2, 1, 0.45), (30, 41, 3, 3, 0.40),
    (12, 13, 1, 1, 0.50), (15, 16, 1, 1, 0.50),
]:
    kap = r0k / r0
    n = max(132, int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kap))) // 2 * 2)
    rat, lam, kap_ex = check_G3_rat(n, c, j, r0)
    print(f"  c={c:2d} j={j:2d} r0={r0} D={D} c/(r0 D^2)={c/(r0*D*D):5.1f} r0*kap={r0*kap_ex:.3f} (n={n}): G3 rat_L1 = {rat:+.4f} (at lam={lam:.2f})")
