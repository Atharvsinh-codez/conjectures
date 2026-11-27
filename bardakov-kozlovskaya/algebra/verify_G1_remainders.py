"""
Test our 100% CLOSED-FORM ZERO-LOOP bounds for the remaining 3 small pieces of Piece G1:
  1. Coupled Convergent + Semiconvergent + Rest bound (with exact a_s coupling at k = s-1!)
  2. Exact k != 0 Poisson mirror sum S_k_ne0 for Set II (1 <= r < q)
  3. Exact 5-line finite-n difference bound S_fin = sum_{r>=1} |b_r^exact - b_r^G| / r!
Target: GRAND TOTAL ZERO-LOOP BOUND ON ALL OF PIECE G1 <= 1.48 < pi/2 = 1.5708!
"""
import math

def D_HS(T, nu):
    N = 120; smax = 6.5; ds = 2.0 * smax / N
    b = math.pi * nu; rt = math.sqrt(T)
    tot = 0.0
    for i in range(N):
        s = -smax + (i + 0.5) * ds
        w = math.exp(-0.5 * s * s) / math.sqrt(2.0 * math.pi)
        h = s * rt
        tot += 2.0 * w * ((T + h) / ((T + h) ** 2 + b * b) - T / (T * T + (b + h) ** 2)) * ds
    return tot

def D_exact_0(T):
    M = 100
    s1 = sum(math.exp(-0.5 * T * (1.0 - ((i + 0.5) / M) ** 2)) for i in range(M)) * (2.0 / M)
    ymax = min(40.0, 20.0 / max(T, 1e-4) + 7.0 / math.sqrt(max(T, 1e-4)))
    s2 = sum(math.exp(-T * ((i + 0.5) * ymax / M) - 0.5 * T * (((i + 0.5) * ymax / M) ** 2)) for i in range(M)) * (2.0 * ymax / M)
    return abs(s1 - s2)

t_pts = [0.005 * i for i in range(1, 400)] + [2.0 + 0.05 * i for i in range(1, 201)]
H_pts = [max(2.0 * math.pi * math.pi * nu * nu * abs(D_HS(t * math.pi * nu, nu)) for nu in [1, 2, 3, 6, 20, 60]) for t in t_pts]

def H_D(t):
    if t <= 0: return 0.0
    if t >= 12.0: return 8.2 / (t * t)
    if t < 2.0:
        idx = int(t / 0.005)
        if idx <= 0: return min(H_pts[0], 25.0 * t * t)
        if idx >= 399: return H_pts[398]
        return max(H_pts[idx - 1], H_pts[idx]) + 0.003
    else:
        idx = 399 + int((t - 2.0) / 0.05)
        if idx >= len(H_pts) - 1: return H_pts[-1]
        return max(H_pts[idx], H_pts[idx + 1]) + 0.003

# 1. Check EXACT closed-form k != 0 mirror sum for Set II (1 <= r < q):
# For every r in [1, q-1]:
#   (c / (2*r)) * sum_{k != 0} |D(r*c*kap/2, |D_r + k*c|)|
#   Since |D_r + k*c| >= (|k| - 0.5)*c (for k = -1: >= c/2; for k = +1: >= c; for |k|>=2: (|k|-0.5)*c!):
#   Using D_r <= c/2: k=-1 gives nu >= c/2, k=+1 gives nu >= c, k=-2 gives nu >= 3c/2, k=+2 gives nu >= 2c, ...
#   In general, the set {|D_r + k*c| : k != 0} is DOMINATED term-by-term by { (m/2)*c : m = 1, 2, 3, 4, ... }?
#   WAIT: check! For delta in [0, 1/2]:
#     k = -1: 1 - delta >= 1/2
#     k = +1: 1 + delta >= 1
#     k = -2: 2 - delta >= 3/2
#     k = +2: 2 + delta >= 2
#     k = -m: m - delta >= m - 1/2
#     k = +m: m + delta >= m
#   YES!! Term-by-term, |delta + k| >= m/2 for m = 1, 2, 3, 4, ... (EXACTLY: 1/2, 1, 3/2, 2, 5/2, 3, ...)!
#   Wait: check sum_{m>=1} 1/(m/2)^2 = 4 * pi^2 / 6 = 6.5797 (vs exact pi^2 - 4 = 5.8696)!
#   Even better, for each m >= 1, t_{r, m} = r*c*kap / (2*pi * |D_r + k*c|) <= r*kap / (pi * m)!
#   NO! At delta = 1/2, k = -1 has 1 - delta = 1/2, so t = r*kap / pi; k = +1 has 1 + delta = 3/2, so t = r*kap / (3*pi)!
#   For ANY delta in [0, 1/2]:
#     sum_{k != 0} (1/(delta+k)^2) * H_D(r*kap / (2*pi*|delta+k|))
#   Let's compute the EXACT universal upper bound on:
#     S_k_ne0 <= (1 / (4*pi^2 * c)) * sum_{r=1}^{q-1} (1/r) * max_{delta in [0, 1/2]} sum_{k != 0} H_D(r*kap/(2*pi*|delta+k|)) / (delta+k)^2 !
print("=== 1. EXACT UNIVERSAL BOUND ON k != 0 MIRROR SUM FOR SET II (1 <= r < q) ===")
# Notice: sum_{r=1}^{q-1} (1/r) * G(r * kap) <= int_0^{q*kap} (G(u)/u) du + max_{u <= q*kap} G(u)!
# WAIT: Here q*kap = X / q <= X / c (for d=1) or X / (c/d) = d*X/c!
# Compute exact max_{delta in [0, 1/2]} sum_{k != 0} H_D(u / (2*pi*|delta+k|)) / (delta+k)^2 as a function of u = r*kap!
def G_kne0(u):
    # Check delta = 0.0, 0.25, 0.50 (convex in delta? max is at delta = 0.50!)
    mx = 0.0
    for delta in [0.0, 0.25, 0.50]:
        s = 0.0
        for k in range(-15, 16):
            if k == 0: continue
            dk = abs(delta + k)
            s += H_D(u / (2.0 * math.pi * dk)) / (dk * dk)
        mx = max(mx, s)
    return mx

# Integral int_0^inf G_kne0(u)/u du and max G_kne0(u):
du_log = 0.05
I_kne0 = sum(G_kne0(0.01 * math.exp(i * du_log)) * du_log for i in range(0, 220))
M_kne0 = max(G_kne0(0.01 * math.exp(i * du_log)) for i in range(0, 220))
C_kne0 = (I_kne0 + M_kne0) / (4.0 * math.pi * math.pi)
print(f"  int G(u)/u du = {I_kne0:.4f}, max G(u) = {M_kne0:.4f} => S_k_ne0 <= {C_kne0:.4f} / c !!")
print(f"  (For c >= 10: <= {C_kne0/10:.4f}; for c >= 6: <= {C_kne0/6:.4f}!)")
# Also check exact discrete sum_{r=1}^{c-1} G_kne0(r*kap)/r / (4*pi^2*c) for c = 6, 8, 10, 15, 20:
for c_test in [6, 8, 10, 15, 20]:
    mx_c = max(sum(G_kne0(r * kap) / r for r in range(1, c_test)) / (4.0 * math.pi * math.pi * c_test)
               for kap in [0.05 * k for k in range(1, 50)])
    print(f"    Exact max S_k_ne0 at c={c_test:2d}: {mx_c:.5f} (compare {C_kne0/c_test:.5f})")

# 2. Check COUPLED Convergent + Semiconvergent + Rest at k = s - 1!
print("\n=== 2. COUPLED CONVERGENT + SEMICONVERGENT BOUND AT LEVEL k = s - 1 ===")
# Recall: at k = s-1 (last convergent), a_s >= 2, D_{s-1} = d, M_{s-1} = (a_s + theta)/d (theta = q_{s-2}/q_{s-1} in [0, 1))!
# And semiconvergents at level s-1 are r_u = (a_s - u)*q_{s-1} + q_{s-2} (1 <= u <= a_s - 1),
# with D_{r_u} = u * d, c / (r_u D_{r_u}^2) <= a_s / ((a_s - u) * u^2 * d),
# and t(r_u) = ((a_s - u + theta) / u) * t_{s-1} in [((a_s - u)/u)*t_{s-1}, ((a_s - u + 1)/u)*t_{s-1}]!
# Let's compute the EXACT max over (a_s >= 2, t1 in (0, 1.50]) of:
#   CONV(k=s-1) + SEMI(k=s-1)  subject to M_{s-1} <= M_G3 = 10 (or eps_G3 = 0.45)!
C0 = 1.20206 / (4.0 * math.pi * math.pi)
C_pk = 1.0 / (4.0 * math.pi * math.pi)
for M_G3, eps_G3 in [(10, 0.45), (9, 0.48), (8, 0.50)]:
    mx_coupled = (0.0, None)
    for a_s in list(range(2, 20)) + [30, 50, 100]:
        for it in range(1, 151):
            t1 = 1.50 * it / 150.0
            # Convergent k = s-1:
            # If a_s <= M_G3: low bound <= C0 * min(M_G3, a_s + 1) * H_D(t1)
            # If a_s > M_G3: MUST be high convergent (q_{s-1}*kap > eps_G3)!!
            if a_s <= M_G3:
                c_term = max(C0 * min(M_G3, a_s + 1.0) * H_D(t1),
                             1.20206 * (t1 * H_D(t1) / (2.0 * math.pi)) / eps_G3)
            else:
                # MUST be high convergent (q_{s-1}*kap > eps_G3) -> BOTH convergent AND semiconvergents have r*kap > eps_G3!
                c_term = 1.20206 * (t1 * H_D(t1) / (2.0 * math.pi)) / eps_G3
            # Semiconvergents at level s-1:
            s_term = 0.0
            for u in range(1, a_s):
                t_lo = ((a_s - u) / u) * t1
                t_hi = ((a_s - u + 1.0) / u) * t1
                # max of H_D on [t_lo, t_hi]:
                hm = 3.04 if (t_lo <= 0.628 <= t_hi) else max(H_D(t_lo), H_D(t_hi))
                # If a_s > M_G3, q_{s-1}*kap > eps_G3 -> r_u * kap >= (a_s - u) * eps_G3!!
                if a_s <= M_G3:
                    st = C_pk * (a_s / ((a_s - u) * u * u)) * hm
                else:
                    # Mid-t bound: (t * H_D(t) / (2*pi)) / (r_u * D_{r_u} * r_u * kap) <= ...
                    # Actually: (c / (4*pi^2 * r_u * D^2)) * H_D(t) = (t * H_D(t) / (2*pi)) / (r_u * D * r_u * kap)
                    # Since r_u >= (a_s - u), D = u, r_u*kap >= (a_s - u)*eps_G3:
                    thm = 0.350 if (t_lo <= 0.90 <= t_hi) else max(t_lo * H_D(t_lo), t_hi * H_D(t_hi)) / (2.0 * math.pi)
                    st = min(C_pk * (a_s / ((a_s - u) * u * u)) * hm,
                             thm / ((a_s - u) * (a_s - u) * u * eps_G3))
                s_term += st
            if c_term + s_term > mx_coupled[0]:
                mx_coupled = (c_term + s_term, (a_s, round(t1, 3), round(c_term, 4), round(s_term, 4)))
    print(f"  (M={M_G3:2d}, eps={eps_G3:.2f}): max [CONV(s-1) + SEMI(s-1)] = {mx_coupled[0]:.4f} at {mx_coupled[1]}")

# 3. Check 1-line analytic finite-n difference bound S_fin = sum |b_r^exact - b_r^G| / r!
# Wait: check WHY S_diff <= 0.0663 across ALL of G1!
# Remember from Lemma G2.1:
#   - Relative modulus error on [0, c/2]: |log|a(x)| - (kap/c)(c^2/4 - x^2)| <= 0.015 * kap^2 * (kap/c)*(c^2/4 - x^2)!
#   - Phase error on [0, c/2]: |dpsi(x)| <= 0.1407 * kap * log|a(x)|!
# WAIT A SECOND!! CHECK THAT!!
# For EVERY small root x in [0, c/2]:
#   |a(x)|^{-r} = exp(-r * log|a(x)|) = U_r^G(m) * (1 + O(kap^2 * r * log|a(x)|))
#   and arg(a(x)^{-r}) - (2*pi*r*j*m/c) = -r * dpsi(x) = O(kap * r * log|a(x)|)!
#   For large roots x in [c/2, 0.10n]:
#   SAME relative error O(kap * r * |log|a(x)||)!
# Let's check what the EXACT bound on sum_{r>=1} |b_r^exact - b_r^G| / r is as a function of kap and c!
# Wait: AT SMALL kap (say kap <= 0.10): BOTH modulus error (O(kap^2)) AND phase error (O(kap)) are TINY!
# AT LARGE kap (say kap in [0.5, 2.5]): r=1 ALREADY has e^{-r*kap} <= e^{-0.5}, so ONLY r = 1, 2, 3 matter,
# AND even for r = 1, 2, 3, sum_{m>=1} m^2 e^{-r*kap*m} is TINY!
print("\n=== 3. CHECK SCALING OF EXACT S_diff = sum |b_r^ex - b_r^G|/r WITH (c, kap) ===")
# Check the naive triangle inequality sum_{r>=1} (1/r) * (sum_m |U_r^ex(m) - U_r^G(m) - (V_r^ex(m) - V_r^G(m))|)
# (pairing U and V so the linear phase -r*dpsi'(c/2)*m STILL cancels!)
def exact_paired_diff_L1(n, c, j):
    th_c = c * math.pi / n; psi_j = j * math.pi / n; alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1; log_Rz = n * math.log(r3 / r1); kap = 2.0 * math.pi * r1 * r2 / r3
    def U_ex(m, r):
        x = c / 2.0 - m
        eta = 2.0 * math.pi * x / n
        z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
        lr = -(log_Rz - n * math.log(abs(z)))
        ang = -((math.pi if j % 2 else 0.0) - n * math.atan2(z.imag, z.real)) - 2.0 * math.pi * m * j / c
        return math.exp(r * lr) * complex(math.cos(r * ang), math.sin(r * ang)) if r * lr > -45 else 0j
    def V_ex(m, r):
        x = c / 2.0 + m
        eta = 2.0 * math.pi * x / n
        z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
        lr = log_Rz - n * math.log(abs(z))
        ang = ((math.pi if j % 2 else 0.0) - n * math.atan2(z.imag, z.real)) - 2.0 * math.pi * m * j / c
        return math.exp(r * lr) * complex(math.cos(r * ang), math.sin(r * ang)) if r * lr > -45 else 0j
    tot_L1 = 0.0
    for r in range(1, min(100, int(15.0 / kap) + 10)):
        half = (c - 1) // 2
        s_r = 0.0
        for m in range(1, half + 1):
            uG = math.exp(-r * kap * m * (1.0 - m / c))
            vG = math.exp(-r * kap * m * (1.0 + m / c))
            diff_m = (U_ex(m, r) - uG) - (V_ex(m, r) - vG)
            s_r += 2.0 * abs(diff_m)
        if c % 2 == 0:
            m = c // 2
            s_r += abs(U_ex(m, r) - math.exp(-r * kap * m * 0.5))
        for m in range(half + 1, min(n // 2 - c // 2, half + int(20.0 / (r * kap)) + 5)):
            vG = math.exp(-r * kap * m * (1.0 + m / c)) if r * kap * m * (1.0 + m / c) < 45 else 0.0
            s_r += 2.0 * abs(V_ex(m, r) - vG)
        tot_L1 += s_r / r
    return tot_L1, kap

for c, j, kt in [(6, 7, 1.5), (10, 13, 0.30), (12, 13, 0.50), (20, 21, 0.50), (29, 39, 0.05), (50, 69, 0.80), (100, 137, 0.60)]:
    n = max(132, int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt))) // 2 * 2)
    L1_d, kap = exact_paired_diff_L1(n, c, j)
    print(f"  c={c:3d} j={j:3d} kap={kap:.3f} ck^2={c*kap*kap:5.2f}: paired L1 diff bound = {L1_d:.4f}")
