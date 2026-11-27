"""
1. Verify exact finite-n 4th-order paired Abel identity to machine precision (< 1e-13)!
2. Measure worst-case exact finite-n constants K4_paired and K0 across ALL of G* = {kap < 2.5, c*kap^2 < 42, c*q*kap > 0.85, n >= 132}.
3. Check the EXACT Dirichlet Minor-Arc sum bound:
     Phi_bound = sum_{r >= 1} min( K0 / (r^3 * c * kap^2),  K4 * r * kap^2 / (16 * c * sin^4(pi*r*j/c)) )
   to determine the EXACT thresholds (d_min, dq^4k^2_max for G2; c/(r0*D^2)_min, r0*kap_max for G3) where Phi_bound <= 1.45 < pi/2!
"""
import math
import cmath

def exact_UV(n, c, j, r):
    th_c = c * math.pi / n
    psi_j = j * math.pi / n
    alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1
    log_Rz = n * math.log(r3 / r1)
    kap = 2.0 * math.pi * r1 * r2 / r3
    w0 = cmath.exp(2j * math.pi * r * j / c)
    sgn_j = 0.0 if (j % 2 == 0) else math.pi

    def a_val(eta):
        z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
        lr = log_Rz - n * math.log(abs(z))
        ang = sgn_j - n * math.atan2(z.imag, z.real)
        return lr, ang

    U = []
    for m in range(c + 1):
        lr, ang = a_val(th_c - 2.0 * m * math.pi / n)
        rlr = -r * lr
        val = (math.exp(rlr) if rlr > -700 else 0.0) * cmath.exp(1j * (-r * ang - 2.0 * math.pi * r * j * m / c))
        U.append(val)

    V = []
    M_max = (n - c) // 2 + 2
    for m in range(M_max + 1):
        eta = th_c + 2.0 * m * math.pi / n
        if eta > math.pi + 1e-12:
            break
        lr, ang = a_val(eta)
        rlr = r * lr
        val = (math.exp(rlr) if rlr > -700 else 0.0) * cmath.exp(1j * (r * ang - 2.0 * math.pi * r * j * m / c))
        V.append(val)
        if m > c + 15 and abs(val) < 1e-18:
            break

    Br = sum((U[m] * (w0 ** m)).real for m in range(1, c))
    Ar = 0.0
    for m in range(1, len(V)):
        eta = th_c + 2.0 * m * math.pi / n
        wt = 1.0 if abs(eta - math.pi) < 1e-12 else 2.0
        Ar += wt * (V[m] * (w0 ** m)).real
    return U, V, Br - Ar, kap, w0

def test_abel4_exact(n, c, j, r):
    U, V, br, kap, w0 = exact_UV(n, c, j, r)
    s4 = 16.0 * (math.sin(math.pi * r * j / c) ** 4)
    if s4 < 1e-12:
        return None
    def Uget(m):
        m = m % c
        return U[m]
    def Vget(m):
        if m < 0: return V[-m].conjugate()
        if m < len(V): return V[m]
        return 0j
    def D2(f, m): return f(m + 1) - 2.0 * f(m) + f(m - 1)
    def D4(f, m): return D2(f, m + 1) - 2.0 * D2(f, m) + D2(f, m - 1)

    bdy = 2.0 * ((Uget(2) - Vget(2)) - 4.0 * (Uget(1) - Vget(1))).real
    sum_D4_U = sum((D4(Uget, m) * (w0 ** m)).real for m in range(1, c))
    sum_D4_V = 2.0 * sum((D4(Vget, m) * (w0 ** m)).real for m in range(1, len(V) - 2))
    rhs = bdy + sum_D4_U - sum_D4_V
    lhs = s4 * br
    err = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-18)

    # Paired L1 norm:
    half = c // 2
    paired_L1 = abs(bdy)
    for m in range(1, half):
        paired_L1 += 2.0 * abs(D4(Uget, m) - D4(Vget, m))
    if c % 2 == 0:
        paired_L1 += abs(D4(Uget, half) - 2.0 * D4(Vget, half))
    else:
        paired_L1 += 2.0 * abs(D4(Uget, half) - D4(Vget, half))
    for m in range(half + 1, len(V) - 2):
        paired_L1 += 2.0 * abs(D4(Vget, m))

    scale = (r * r * kap * kap) / c
    # Also check paired L1 bound on |br| directly:
    L1_br = sum(abs(U[m] - V[m]) for m in range(1, min(c, len(V)))) + sum(abs(V[m]) for m in range(c, len(V)))  # crude
    return err, abs(lhs) / scale, paired_L1 / scale, abs(br) * (r * r * c * kap * kap), kap

print("=== 1. EXACT IDENTITY VERIFICATION & WORST CONSTANTS IN G* ===")
def mk_n(c, j, kt):
    n = int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt)))
    return max(132, n + (n % 2))

max_err = 0.0
worst_pK4 = (0.0, None)
worst_tK4 = (0.0, None)
worst_K0 = (0.0, None)
for c, j, kt in [
    (2, 3, 2.0), (3, 4, 1.5), (4, 5, 1.0), (5, 7, 1.2), (6, 7, 0.8),
    (7, 11, 0.5), (10, 13, 0.4), (10, 13, 1.8), (15, 22, 1.5),
    (20, 27, 0.25), (30, 47, 0.8), (50, 81, 0.3), (100, 162, 0.15),
    (200, 203, 0.04), (200, 201, 0.01), (200, 201, 0.005),
    (30, 31, 0.05), (80, 107, 0.015), (50, 81, 0.02),
]:
    d = math.gcd(c, j); q = c // d
    n = mk_n(c, j, kt)
    if c * (kt ** 2) >= 42.0 or c * q * kt <= 0.85:
        continue
    for r in [1, 2, 3, 5]:
        out = test_abel4_exact(n, c, j, r)
        if out is None: continue
        err, tK4, pK4, K0, kap = out
        max_err = max(max_err, err)
        if pK4 > worst_pK4[0]: worst_pK4 = (pK4, (n, c, j, r, round(kap, 4), round(c*kap, 3)))
        if tK4 > worst_tK4[0]: worst_tK4 = (tK4, (n, c, j, r, round(kap, 4), round(c*kap, 3)))
        if K0 > worst_K0[0]: worst_K0 = (K0, (n, c, j, r, round(kap, 4), round(c*kap, 3)))

print(f"  MAX identity error: {max_err:.2e}  (SHOULD BE < 1e-12!)")
print(f"  WORST true_K4:   {worst_tK4[0]:.3f} (K4/16 = {worst_tK4[0]/16:.3f}) at {worst_tK4[1]}")
print(f"  WORST paired_K4: {worst_pK4[0]:.3f} (K4/16 = {worst_pK4[0]/16:.3f}) at {worst_pK4[1]}")
print(f"  WORST K0:        {worst_K0[0]:.3f} at {worst_K0[1]}")

print("\n=== 2. CHECK DIRICHLET MINOR-ARC SUM BOUND WITH (K0, K4/16) ===")
# For a given (c, j, kap) in G*, compute:
#   S_bound = sum_{r=1}^{500} min( K0 / (r^3 * c * kap^2),  (K4/16) * r * kap^2 / (c * sin^4(pi*r*j/c)) )
# Wait: check exact K0 = 12, K4/16 = 3.5 (conservative!) vs K0 = 9, K4/16 = 3.1!
def sum_bound_G1(c, j, kap, K0=12.0, C4=3.5):
    tot = 0.0
    # Sum r = 1..500 + tail r > 500 bounded by K0/(c*kap^2) * 1/(2*500^2)
    for r in range(1, 501):
        b_l1 = K0 / ((r ** 3) * c * kap * kap)
        s = abs(math.sin(math.pi * r * j / c))
        if s < 1e-12:
            tot += b_l1
        else:
            b_a4 = C4 * r * kap * kap / (c * (s ** 4))
            tot += min(b_l1, b_a4)
    tot += K0 / (c * kap * kap * 2.0 * 500.0 * 500.0)
    return tot

# Test across a HUGE sweep of (c, j, kap) in G* to see what points have S_bound > 1.45!
def best_dirichlet(c, j, max_r0=30):
    best = None
    for r0 in range(1, min(max_r0, c) + 1):
        p0 = round(r0 * j / c)
        if math.gcd(r0, p0) != 1: continue
        D = r0 * j - c * p0
        score = (0, r0) if D == 0 else (1, r0 * D * D)
        if best is None or score < best[0]:
            best = (score, r0, p0, D)
    return best[1], best[2], best[3]

for K0, C4 in [(10.0, 3.2), (12.0, 3.5), (15.0, 4.0)]:
    hits_d0 = []
    hits_dn = []
    for c in [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 50, 100, 200]:
        for j in range(c + 1, min(2 * c + 1, c + 40)):
            d = math.gcd(c, j); q = c // d
            for kap in [0.005, 0.01, 0.02, 0.04, 0.08, 0.15, 0.3, 0.6, 1.2, 2.0]:
                if c * kap * kap >= 42.0 or c * q * kap <= 0.85: continue
                sb = sum_bound_G1(c, j, kap, K0=K0, C4=C4)
                if sb > 1.45:
                    r0, p0, D = best_dirichlet(c, j, max_r0=30)
                    if D == 0:
                        hits_d0.append((c / r0, r0 * kap, (c / r0) * ((r0 * r0 * kap) ** 2), c, j, kap, r0, sb))
                    else:
                        hits_dn.append((c / (r0 * D * D), r0 * kap, c, j, kap, r0, D, sb))
    print(f"\n  With K0={K0}, C4={C4}:")
    if hits_d0:
        print(f"    D=0  hits={len(hits_d0):4d}: min(d)={min(h[0] for h in hits_d0):.1f}, max(q*kap)={max(h[1] for h in hits_d0):.3f}, max(d*q^4*kap^2)={max(h[2] for h in hits_d0):.2f}")
    if hits_dn:
        print(f"    D!=0 hits={len(hits_dn):4d}: min(c/(r0*D^2))={min(h[0] for h in hits_dn):.2f}, max(r0*kap)={max(h[1] for h in hits_dn):.3f}")
        for h in sorted(hits_dn, key=lambda x: x[0])[:5]:
            print(f"      worst c/(r0*D^2): c/(r0*D^2)={h[0]:.2f}, r0*kap={h[1]:.3f}, c={h[2]}, j={h[3]}, kap={h[4]}, r0={h[5]}, D={h[6]}, sb={h[7]:.2f}")
