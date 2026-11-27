"""
Test the clean overlap between G1 (Minor Arcs), G2 (D=0 Major Arcs), and G3 (D!=0 Major Arcs):
1. Test G2 (D=0) across ALL d >= 1 (including d = 1, 2, 3, 4!) for d*(q^2*kap)^2 <= 16 (and/or q^2*kap <= 2.0)!
2. Test G1 (Minor Arcs) using exact finite-interval Gaussian b_r^G(c, j, kap, r) AND our analytic bound
   min(b0_exact(T), b2_exact(theta), b4_exact(theta)) across ALL Minor Arcs outside G2 and G3!
"""
import math
import cmath

# First: check exact 1-variable function D(T) = |int_{-1}^1 exp(-T/2*(1-u^2)) du - 2 int_1^inf exp(-T/2*(u^2-1)) du|
# for T = r*c*kap/2 > 0!
def D_exact(T):
    # Numerical quadrature of exact Gaussian continuum integral:
    M = 400
    s1 = sum(math.exp(-0.5 * T * (1.0 - ((i + 0.5) / M) ** 2)) for i in range(M)) * (2.0 / M)
    # 2 * int_0^inf exp(-T*y - 0.5*T*y^2) dy:
    ymax = min(50.0, 25.0 / max(T, 1e-4) + 10.0 / math.sqrt(max(T, 1e-4)))
    s2 = sum(math.exp(-T * ((i + 0.5) * ymax / M) - 0.5 * T * (((i + 0.5) * ymax / M) ** 2)) for i in range(M)) * (2.0 * ymax / M)
    return abs(s1 - s2), s1 - s2

print("=== 1. EXACT 1-VARIABLE FUNCTION D(T) (Order 0 bound: |B_r - A_r| <= (c/2) * D(T)) ===")
for T in [0.425, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0]:
    ad, sd = D_exact(T)
    # Compare with 4/T^2 (which corresponds to 8/(r^2 c kap^2)):
    print(f"  T={T:6.3f}: s1-s2={sd:+.5f}, |D(T)|={ad:.5f}, |D(T)|*T^2={ad*T*T:.4f} (asymptotic = 4.000)")

# Check exact finite-n G2 for d = 1, 2, 3, 4, 5 at q^2*kap in [0.85/d, 2.0]:
import sys
sys.path.insert(0, r'd:\solving-conjectures\work\bardakov-kozlovskaya\algebra')

def exact_finite_n_roots(n, c, j):
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
            if log_r > -55.0:
                roots.append((log_r, ang)); added = True
        if not added and k > c: break
        k += 1
        if k > n // 2 + 2: break
    return roots, 2.0 * math.pi * r1 * r2 / r3

def solve_Gq_saddle(roots, c, q):
    d = c // q; x = 0.0
    for _ in range(40):
        val = -float(d); der = 0.0
        for lr, ang in roots:
            lmq = q * (lr + x); angq = q * ang
            if lmq < -45.0: continue
            elif lmq > 45.0: val += 1.0 / q
            else:
                u = math.exp(lmq) * complex(math.cos(angq), math.sin(angq))
                t = u / (1.0 + u); val += t.real / q; der += (t / (1.0 + u)).real
        if abs(val) < 1e-12 or der <= 0: break
        x -= max(-0.2, min(0.2, val / der))
    H_psi = 0.0
    for lr, ang in roots:
        lmq = q * (lr + x); angq = q * ang
        if abs(lmq) > 45.0: continue
        u = math.exp(lmq) * complex(math.cos(angq), math.sin(angq))
        H_psi += (u / ((1.0 + u) ** 2)).real / q
    return math.exp(x), H_psi

print("\n=== 2. TEST G2 (D=0) FOR SMALL d = 1, 2, 3, 4 (EXACT FINITE-n) ===")
for d, q, p0, kt in [
    (1, 2, 3, 0.30), (1, 3, 4, 0.15), (1, 4, 5, 0.08), (1, 5, 6, 0.04), (1, 10, 11, 0.015),
    (2, 1, 2, 0.60), (2, 2, 3, 0.25), (2, 3, 4, 0.10), (2, 5, 6, 0.04), (2, 10, 11, 0.01),
    (3, 1, 2, 0.30), (3, 1, 2, 1.00), (3, 2, 3, 0.20), (3, 4, 5, 0.05),
    (4, 1, 2, 0.30), (4, 1, 2, 1.20), (4, 2, 3, 0.20), (4, 5, 6, 0.02),
]:
    c = d * q; j = d * p0
    n = int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt)))
    n = max(132, n + (n % 2))
    roots, kap = exact_finite_n_roots(n, c, j)
    lam, H_psi = solve_Gq_saddle(roots, c, q)
    # Check exact rat at lam:
    M = max(256, 64 * q)
    sgn_pref = 1.0 if ((j * (c + 1) + c) % 2 == 0) else -1.0
    log_lam = math.log(lam)
    vals = []
    mx = -1e300
    for i in range(M):
        phi = 2.0 * math.pi * (i + 0.5) / M
        la = -c * log_lam; ph = -c * phi
        for lr, ang in roots:
            lm = lr + log_lam; th = ang + phi
            if lm < -45.0: continue
            elif lm > 45.0: la += lm; ph += th + math.pi
            else:
                r = math.exp(lm); z = complex(1.0 - r * math.cos(th), -r * math.sin(th))
                la += math.log(abs(z)); ph += math.atan2(z.imag, z.real)
        vals.append((la, ph))
        if la > mx: mx = la
    s_re = sum(math.exp(la - mx) * math.cos(ph + (0.0 if sgn_pref > 0 else math.pi)) for la, ph in vals if la - mx > -700)
    s_abs = sum(math.exp(la - mx) for la, ph in vals if la - mx > -700)
    print(f"  d={d} q={q:2d} c={c:2d} kap={kap:.4f} dq^2k={d*q*q*kap:5.2f} q^2k={q*q*kap:5.2f} | lam={lam:.4f} H_psi={H_psi:6.2f} | rat={s_re/s_abs:+.5f}")

print("\n=== 3. TEST EXACT G1 MINOR-ARC BOUND WITH FINITE-T ORDER 0 + EXACT b_r^main ===")
# Exact finite-T Order 0 bound: for T = r*c*kap/2:
#   |b_r| / r <= (1/r) * min(9.5 / (r^2 * c * kap^2),  c * (0.45 / T^0.5 if T < 1.5 else 4.5 / T^2)) ...
# Wait: from D(T) above, let's check exact Gaussian |b_r^G| / r on ALL Minor Arcs outside:
#   G2: d >= 1 and d*(q^2*kap)^2 <= 18.0  (Wait: or q^2*kap <= 2.0 AND d*(q^2*kap)^2 <= 25.0!)
#   G3: exists r0 >= 1 coprime to p0 with D != 0, c/(r0*D^2) >= 15.0, r0*kap <= 0.25!
def br_gauss_exact(c, j, kap, r):
    """Exact Gaussian b_r = B_r(1) - A_r(1) (true finite-c sum!)."""
    alpha = kap / c
    shift = 0.0 if (c % 2 == 0) else 0.5
    # B_r(1) = sum_{|x| < c/2} exp(r*alpha*(x^2 - c^2/4)) * cos(2*pi*r*j*(x - c/2)/c)
    # Wait: exact phase is exp(-2*pi*i*r*j*(x - c/2)/c) -> real part is cos(2*pi*r*j*(x - c/2)/c)!
    Br = 0.0
    for m in range(1, c):
        x = 0.5 * c - m
        Br += math.exp(r * alpha * (x * x - 0.25 * c * c)) * math.cos(2.0 * math.pi * r * j * m / c)
    Ar = 0.0
    m = 1
    while True:
        x = 0.5 * c + m
        ex = -r * alpha * (x * x - 0.25 * c * c)
        if ex < -35.0: break
        Ar += 2.0 * math.exp(ex) * math.cos(2.0 * math.pi * r * j * m / c)
        m += 1
        if m > 5000: break
    return Br - Ar

worst_sum_G1 = (0.0, None)
over_14 = []
for c in [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 50, 100]:
    for j in range(c + 1, min(2 * c + 1, c + 35)):
        d = math.gcd(c, j); q = c // d
        for kap in [0.005, 0.01, 0.02, 0.04, 0.08, 0.15, 0.30, 0.60, 1.20, 2.00]:
            if c * kap * kap >= 42.0 or c * q * kap <= 0.85: continue
            # Exclude G2: d*(q^2*kap)^2 <= 18.0:
            if d * ((q * q * kap) ** 2) <= 18.0: continue
            # Exclude G3: exists r0 coprime to p0, D != 0, c/(r0*D^2) >= 15.0, r0*kap <= 0.25:
            in_G3 = False
            for r0 in range(1, min(c, int(0.25 / kap) + 1) + 1):
                if r0 * kap > 0.25 + 1e-12: break
                p0 = round(r0 * j / c); D = r0 * j - c * p0
                if D != 0 and math.gcd(r0, p0) == 1 and c / (r0 * D * D) >= 15.0:
                    in_G3 = True; break
            if in_G3: continue
            # Sum true |b_r^G| / r over r = 1..150:
            tot = 0.0
            for r in range(1, 151):
                tot += abs(br_gauss_exact(c, j, kap, r)) / r
            if tot > worst_sum_G1[0]:
                worst_sum_G1 = (tot, (c, j, d, q, kap, round(d*(q*q*kap)**2, 2)))
            if tot > 1.40:
                over_14.append((tot, c, j, d, q, kap))
print(f"  TRUE sum |b_r^G|/r across ALL Minor Arcs (outside G2: dq^4k^2<=18, G3: c/(r0 D^2)>=15, r0*k<=0.25):")
print(f"    WORST = {worst_sum_G1[0]:.4f} at {worst_sum_G1[1]}; hits > 1.40: {len(over_14)}")
for h in sorted(over_14, reverse=True)[:10]:
    print(f"    tot={h[0]:.4f}: c={h[1]}, j={h[2]}, d={h[3]}, q={h[4]}, kap={h[5]}")
