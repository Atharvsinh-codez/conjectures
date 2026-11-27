"""
Test the 100% EXACT FINITE-n 4th-order paired Abel identity for G1 (Minor Arcs):
  b_r := B_r(1) - A_r(1) = Re sum_{m=1}^{c-1} U_r(m) w0^m - Re sum_{m=1}^{n-c-1} V_r(m) w0^m,
where w0 = exp(2*pi*i * r * j / c), U_r(m) = w0^{-m} a(th_c - 2m*pi/n)^{-r}, V_r(m) = w0^{-m} a(th_c + 2m*pi/n)^r.

1. Verify the EXACT algebraic identity (zero approximation!) for 16*sin^4(pi*r*j/c) * b_r!
2. Measure the EXACT finite-n ratio:
     K4_eff = (16 * sin^4(pi*r*j/c) * |b_r|) / (r^2 * kap^2 / c)
   AND the exact L1 ratio:
     K0_eff = |b_r| / (1 / (r^2 * c * kap^2))
   across diverse (n, c, j) in G* = {kap < 2.5, c*kap^2 < 42, c*q*kap > 0.85, n >= 132}!
3. Compute the exact sum sum_{r>=1} min(K0/(r^3 c kap^2), K4*r*kap^2 / (16 c sin^4(pi*r*j/c)))
   across ALL minor arcs to determine exact G1/G2 and G1/G3 overlap thresholds!
"""
import math
import cmath

def exact_UV(n, c, j, r):
    """Compute exact finite-n U_r(m) for m=0..c and V_r(m) for m=0..M_lg (where |V_r| > 1e-18)."""
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

    # U_r(m) = w0^{-m} * a(th_c - 2m*pi/n)^{-r} for m = 0..c:
    U = []
    for m in range(c + 1):
        lr, ang = a_val(th_c - 2.0 * m * math.pi / n)
        rlr = -r * lr
        val = (math.exp(rlr) if rlr > -700 else 0.0) * cmath.exp(1j * (-r * ang - 2.0 * math.pi * r * j * m / c))
        U.append(val)

    # V_r(m) = w0^{-m} * a(th_c + 2m*pi/n)^r for m = 0..M_lg (stop when |V| < 1e-16 and m >= c):
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
        if m > c + 10 and abs(val) < 1e-18:
            break

    # Exact B_r(1) and A_r(1):
    # B_r(1) = Re sum_{m=1}^{c-1} U[m] * w0^m
    Br = sum((U[m] * (w0 ** m)).real for m in range(1, c))
    # A_r(1) = 2 * Re sum_{m=1}^{len(V)-1} V[m] * w0^m (taking care if last point is eta = pi for even c!)
    Ar = 0.0
    for m in range(1, len(V)):
        eta = th_c + 2.0 * m * math.pi / n
        wt = 1.0 if abs(eta - math.pi) < 1e-12 else 2.0
        Ar += wt * (V[m] * (w0 ** m)).real
    return U, V, Br - Ar, kap, w0

def check_exact_abel4(n, c, j, r):
    """Verify exact 4th-order Abel identity and compute L1 norm of Delta^4(U) and Delta^4(V) (paired!)."""
    U, V, br, kap, w0 = exact_UV(n, c, j, r)
    s4 = 16.0 * (math.sin(math.pi * r * j / c) ** 4)
    if s4 < 1e-12:
        return None
    # Extend U by U[-m] = conj(U[m]), U[c+m] = conj(U[c-m]) = U[m], and V by V[-m] = conj(V[m]), V[m>=len(V)] = 0:
    def Uget(m):
        if 0 <= m <= c: return U[m]
        if m < 0: return U[-m].conjugate()
        return U[m - c]  # since U[c-k] = conj(U[k]) -> U[c+k] = conj(U[c-k]) = U[k]!
    def Vget(m):
        if m < 0: return V[-m].conjugate()
        if m < len(V): return V[m]
        return 0j

    # Wait: check exact identity!
    # Let S_U = Re sum_{m=1}^{c-1} U[m] w0^m, S_V = 2 Re sum_{m=1}^{M} V[m] w0^m.
    # For any sequence f on Z with f(0)=1, f(-m)=conj(f(m)):
    #   (2 - w0 - w0^{-1}) * (2 Re sum_{m=1}^M f(m) w0^m) = 2 Re[f(0) w0 - f(1)] - 2 Re sum_{m=1}^M D2f(m) w0^m ... wait:
    # Let's check exact signs numerically!
    def D2(f, m): return f(m + 1) - 2.0 * f(m) + f(m - 1)
    def D4(f, m): return D2(f, m + 1) - 2.0 * D2(f, m) + D2(f, m - 1)

    bdy = 2.0 * ((Vget(2) - Uget(2)) - 4.0 * (Vget(1) - Uget(1))).real
    sum_D4_V = 2.0 * sum((D4(Vget, m) * (w0 ** m)).real for m in range(1, len(V) - 2))
    sum_D4_U = sum((D4(Uget, m) * (w0 ** m)).real for m in range(1, c))
    rhs = bdy + sum_D4_V - sum_D4_U
    lhs = s4 * br
    err = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-18)

    # Compute PAIRED L1 norm: split sum_{m=1}^{c-1} D4 U(m) w0^m into 2 Re sum_{m=1}^{half} (handling center m=c/2 cleanly)
    # and pair with 2 Re sum_{m=1}^{half} D4 V(m) w0^m!
    half = c // 2
    paired_L1 = abs(bdy)
    for m in range(1, half):
        paired_L1 += 2.0 * abs(D4(Vget, m) - D4(Uget, m))
    if c % 2 == 0:
        paired_L1 += abs(2.0 * D4(Vget, half) - D4(Uget, half))
    else:
        paired_L1 += 2.0 * abs(D4(Vget, half) - D4(Uget, half))
    for m in range(half + 1, len(V) - 2):
        paired_L1 += 2.0 * abs(D4(Vget, m))

    # Also compute UNPAIRED sum of |D4 V| and |D4 U|:
    unpaired_L1 = abs(bdy) + 2.0 * sum(abs(D4(Vget, m)) for m in range(1, len(V) - 2)) + sum(abs(D4(Uget, m)) for m in range(1, c))

    scale = (r * r * kap * kap) / c
    return err, abs(lhs) / scale, paired_L1 / scale, unpaired_L1 / scale, abs(br) * (r * r * c * kap * kap), kap

print("=== TEST 1: EXACT FINITE-n 4TH-ORDER PAIRED ABEL IDENTITY & CONSTANTS ===")
test_pts = [
    (200, 10, 13, 1), (200, 10, 13, 2), (200, 10, 13, 3),
    (500, 15, 22, 1), (500, 15, 22, 2),
    (2000, 30, 47, 1), (2000, 30, 47, 2), (2000, 30, 47, 5),
    (10000, 50, 81, 1), (10000, 50, 81, 3),
    (50000, 100, 162, 1), (50000, 100, 162, 2), (50000, 100, 162, 5),
    (200000, 200, 203, 1), (200000, 200, 203, 2),
    (5000, 20, 27, 1), (5000, 20, 27, 2),
    (150, 4, 7, 1), (150, 5, 8, 1), (150, 6, 9, 1), (150, 6, 10, 1),
    (1000, 7, 11, 1), (1000, 7, 11, 2), (1000, 8, 13, 1),
]
max_err = 0.0
max_paired_K4 = 0.0
max_true_K4 = 0.0
max_K0 = 0.0
for n, c, j, r in test_pts:
    out = check_exact_abel4(n, c, j, r)
    if out is None: continue
    err, true_K4, p_K4, u_K4, K0, kap = out
    max_err = max(max_err, err)
    max_true_K4 = max(max_true_K4, true_K4)
    max_paired_K4 = max(max_paired_K4, p_K4)
    max_K0 = max(max_K0, K0)
    print(f"n={n:6d} c={c:3d} j={j:3d} r={r} kap={kap:.3f} ck^2={c*kap*kap:5.2f} | "
          f"id_err={err:.1e} | true_K4={true_K4:6.2f} paired_K4={p_K4:7.2f} unp_K4={u_K4:8.2f} | K0={K0:6.2f}")

print(f"\nMAX id_err={max_err:.2e}, max_true_K4={max_true_K4:.2f}, max_paired_K4={max_paired_K4:.2f}, max_K0={max_K0:.2f}")
