"""
Test the THREE-WAY PAIRED ABEL BOUND (Order 0 with gamma_3, Order 2, Order 4) for G1 (Minor Arcs):
  |b_r| / r <= min(
      K0 * gamma_3(r*c*kap/2) / (r^3 * c * kap^2),      [Order 0: direct paired L1]
      (K2 / 4) / (r * c * sin^2(pi*r*j/c)),             [Order 2: 2nd-order paired Abel]
      (K4 / 16) * r * kap^2 / (c * sin^4(pi*r*j/c)),    [Order 4: 4th-order paired Abel]
  )
1. Measure EXACT finite-n worst-case K0, K2/4, K4/16 across G*!
2. Compute S_bound = sum_{r>=1} |b_r|/r across a huge grid of (c, j, kap) in G* = {kap < 2.5, c*kap^2 < 42, c*q*kap > 0.85},
   and properly classify every S_bound > 1.45 hit into G2 (D = 0) vs G3 (D != 0) by its TRUE dominant approximant!
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
        return log_Rz - n * math.log(abs(z)), sgn_j - n * math.atan2(z.imag, z.real)
    U = []
    for m in range(c + 1):
        lr, ang = a_val(th_c - 2.0 * m * math.pi / n)
        rlr = -r * lr
        U.append((math.exp(rlr) if rlr > -700 else 0.0) * cmath.exp(1j * (-r * ang - 2.0 * math.pi * r * j * m / c)))
    V = []
    M_max = (n - c) // 2 + 2
    for m in range(M_max + 1):
        eta = th_c + 2.0 * m * math.pi / n
        if eta > math.pi + 1e-12: break
        lr, ang = a_val(eta)
        rlr = r * lr
        val = (math.exp(rlr) if rlr > -700 else 0.0) * cmath.exp(1j * (r * ang - 2.0 * math.pi * r * j * m / c))
        V.append(val)
        if m > c + 20 and abs(val) < 1e-18: break
    Br = sum((U[m] * (w0 ** m)).real for m in range(1, c))
    Ar = sum((1.0 if abs(th_c + 2.0 * m * math.pi / n - math.pi) < 1e-12 else 2.0) * (V[m] * (w0 ** m)).real for m in range(1, len(V)))
    return U, V, Br - Ar, kap, w0

def measure_exact_K0_K2_K4(n, c, j, r):
    U, V, br, kap, w0 = exact_UV(n, c, j, r)
    def Uget(m): return U[m % c]
    def Vget(m):
        if m < 0: return V[-m].conjugate()
        return V[m] if m < len(V) else 0j
    def D2(f, m): return f(m + 1) - 2.0 * f(m) + f(m - 1)
    def D4(f, m): return D2(f, m + 1) - 2.0 * D2(f, m) + D2(f, m - 1)

    half = c // 2
    # Order 0 paired L1:
    L1_0 = 0.0
    for m in range(1, half):
        L1_0 += 2.0 * abs(Uget(m) - Vget(m))
    if c % 2 == 0:
        L1_0 += abs(Uget(half) - 2.0 * Vget(half))
    else:
        L1_0 += 2.0 * abs(Uget(half) - Vget(half))
    for m in range(half + 1, len(V)):
        L1_0 += 2.0 * abs(Vget(m))

    # Order 2 paired L1 (exact identity: 4 sin^2(th/2) * (Br - Ar) = 2 Re[V(1) - U(1)] + Re sum_{m=1}^{c-1} D2 U w0^m - 2 Re sum_{m>=1} D2 V w0^m):
    bdy2 = 2.0 * (Vget(1) - Uget(1)).real
    L1_2 = abs(bdy2)
    for m in range(1, half):
        L1_2 += 2.0 * abs(D2(Uget, m) - D2(Vget, m))
    if c % 2 == 0:
        L1_2 += abs(D2(Uget, half) - 2.0 * D2(Vget, half))
    else:
        L1_2 += 2.0 * abs(D2(Uget, half) - D2(Vget, half))
    for m in range(half + 1, len(V) - 1):
        L1_2 += 2.0 * abs(D2(Vget, m))

    # Check exact Order 2 identity error:
    s2 = 4.0 * (math.sin(math.pi * r * j / c) ** 2)
    rhs2 = bdy2 + sum((D2(Uget, m) * (w0 ** m)).real for m in range(1, c)) - 2.0 * sum((D2(Vget, m) * (w0 ** m)).real for m in range(1, len(V) - 1))
    # Wait: check sign of Order 2!
    err2 = min(abs(s2 * br - rhs2), abs(s2 * br + rhs2)) / (abs(s2 * br) + abs(rhs2) + 1e-18)

    # Order 4 paired L1:
    bdy4 = 2.0 * ((Uget(2) - Vget(2)) - 4.0 * (Uget(1) - Vget(1))).real
    L1_4 = abs(bdy4)
    for m in range(1, half):
        L1_4 += 2.0 * abs(D4(Uget, m) - D4(Vget, m))
    if c % 2 == 0:
        L1_4 += abs(D4(Uget, half) - 2.0 * D4(Vget, half))
    else:
        L1_4 += 2.0 * abs(D4(Uget, half) - D4(Vget, half))
    for m in range(half + 1, len(V) - 2):
        L1_4 += 2.0 * abs(D4(Vget, m))

    T = 0.5 * r * c * kap
    gam3 = 1.0 - math.exp(-T) * (1.0 + T + 0.5 * T * T) if T < 50 else 1.0
    K0_val = L1_0 / (gam3 / (r * r * c * kap * kap))
    K2_val = L1_2 / (1.0 / c)
    K4_val = L1_4 / ((r * r * kap * kap) / c)
    return err2, K0_val, K2_val / 4.0, K4_val / 16.0, kap

def mk_n(c, j, kt):
    n = int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt)))
    return max(132, n + (n % 2))

print("=== 1. EXACT FINITE-n PAIRED L1 CONSTANTS (K0, K2/4, K4/16) ===")
wK0, wK2, wK4 = (0, None), (0, None), (0, None)
for c, j, kt in [
    (4, 5, 1.0), (5, 7, 1.2), (6, 7, 0.8), (7, 11, 0.5),
    (10, 13, 0.4), (10, 13, 1.8), (15, 22, 1.5), (20, 27, 0.25),
    (30, 47, 0.8), (50, 69, 0.005), (50, 69, 0.02), (50, 81, 0.3),
    (100, 107, 0.3), (100, 162, 0.15), (200, 203, 0.04), (200, 201, 0.005),
]:
    d = math.gcd(c, j); q = c // d
    if c * kt * kt >= 42.0 or c * q * kt <= 0.85: continue
    n = mk_n(c, j, kt)
    for r in [1, 2, 3, 5]:
        err2, k0v, k2v, k4v, kap = measure_exact_K0_K2_K4(n, c, j, r)
        if k0v > wK0[0]: wK0 = (k0v, (c, j, r, round(kap, 4), round(r*c*kap, 2)))
        if k2v > wK2[0]: wK2 = (k2v, (c, j, r, round(kap, 4), round(r*c*kap, 2)))
        if k4v > wK4[0]: wK4 = (k4v, (c, j, r, round(kap, 4), round(r*c*kap, 2)))
print(f"  WORST K0 (paired L1 / gam3): {wK0[0]:.3f} at {wK0[1]}")
print(f"  WORST K2/4 (paired L1):      {wK2[0]:.3f} at {wK2[1]}")
print(f"  WORST K4/16 (paired L1):     {wK4[0]:.3f} at {wK4[1]}")

print("\n=== 2. THREE-WAY MIN SUM BOUND ON ALL OF G* (CLASSIFYING BY DOMINANT TERM!) ===")
def term_r(c, j, kap, r, K0, C2, C4):
    T = 0.5 * r * c * kap
    gam3 = 1.0 - math.exp(-T) * (1.0 + T + 0.5 * T * T) if T < 50 else 1.0
    b0 = K0 * gam3 / ((r ** 3) * c * kap * kap)
    s = abs(math.sin(math.pi * r * j / c))
    if s < 1e-12:
        return b0
    b2 = C2 / (r * c * s * s)
    b4 = C4 * r * kap * kap / (c * (s ** 4))
    return min(b0, b2, b4)

def sum_3way(c, j, kap, K0=10.0, C2=1.5, C4=3.6):
    tot = 0.0
    best_r = (0.0, 1)
    for r in range(1, 401):
        tr = term_r(c, j, kap, r, K0, C2, C4)
        tot += tr
        if tr > best_r[0]:
            best_r = (tr, r)
    tot += K0 / (c * kap * kap * 2.0 * 400.0 * 400.0)
    return tot, best_r[1]

# Now test ALL (c, j, kap) in G* with TRUE dominant r0 classification:
for K0, C2, C4 in [(10.0, 1.5, 3.6), (12.0, 1.8, 4.0)]:
    hits_d0 = []
    hits_dn = []
    for c in [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 50, 100, 200]:
        for j in range(c + 1, min(2 * c + 1, c + 45)):
            d = math.gcd(c, j); q = c // d
            for kap in [0.005, 0.01, 0.02, 0.04, 0.08, 0.15, 0.3, 0.6, 1.2, 2.0]:
                if c * kap * kap >= 42.0 or c * q * kap <= 0.85: continue
                sb, r_dom = sum_3way(c, j, kap, K0=K0, C2=C2, C4=C4)
                if sb > 1.45:
                    # Reduce r_dom to primitive r0 (if D = 0, r0 = q; else r0 = r_dom / gcd(r_dom, p_dom))!
                    p_dom = round(r_dom * j / c)
                    D_dom = r_dom * j - c * p_dom
                    if D_dom == 0:
                        r0 = q; D = 0
                        hits_d0.append((d, q * kap, d * ((q * q * kap) ** 2), c, j, kap, r0, sb))
                    else:
                        g = math.gcd(r_dom, p_dom)
                        r0 = r_dom // g; D = D_dom // g
                        hits_dn.append((c / (r0 * D * D), r0 * kap, c, j, kap, r0, D, sb))
    print(f"\n  With K0={K0}, C2={C2}, C4={C4}:")
    if hits_d0:
        print(f"    D=0  hits={len(hits_d0):4d}: min(d)={min(h[0] for h in hits_d0)}, max(q*kap)={max(h[1] for h in hits_d0):.3f}, max(d*q^4*kap^2)={max(h[2] for h in hits_d0):.2f}")
        for h in sorted(hits_d0, key=lambda x: x[0])[:5]:
            print(f"      smallest d: d={h[0]}, q*kap={h[1]:.3f}, dq^4k^2={h[2]:.2f}, c={h[3]}, j={h[4]}, kap={h[5]}, q={h[6]}, sb={h[7]:.2f}")
    if hits_dn:
        print(f"    D!=0 hits={len(hits_dn):4d}: min(c/(r0*D^2))={min(h[0] for h in hits_dn):.2f}, max(r0*kap)={max(h[1] for h in hits_dn):.3f}")
        for h in sorted(hits_dn, key=lambda x: x[0])[:6]:
            print(f"      smallest c/(r0*D^2): {h[0]:5.2f}, r0*kap={h[1]:.3f}, c={h[2]}, j={h[3]}, kap={h[4]}, r0={h[5]}, D={h[6]}, sb={h[7]:.2f}")
