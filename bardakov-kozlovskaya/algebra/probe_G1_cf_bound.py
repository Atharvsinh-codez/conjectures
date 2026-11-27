"""
Two fast precision checks in one script:
1. G1 (Minor Arcs) boundary vs G3 (c/(r0 D^2) >= M_G3, r0*kap <= eps_G3):
   Check WORST exact sum |b_r|/r on (c, j, kap) outside G2 (d>=2, q^2*kap <= 2.0) and G3 (M_G3=15, eps_G3=0.25),
   specifically probing r0*kap just above 0.25 (0.25..0.60) and c/(r0 D^2) = 8..20!
   Also test a 100% CLOSED-FORM continued-fraction bound for sum_{q \nmid r} |b_r|/r!
2. G2 (D=0 Major Arcs) large-d corner check:
   Test d = 50, 200, 1000 with q = 1, 2, 3, 5, 10 at q^2*kap = 2.0 (where c*kap^2 < 42)
   to verify exact finite-n delta_psi impact on H_q, lam_*, and strip-Poisson R_{non-q}!
"""
import math
import cmath

def b_r_gauss(r, c, j, kap):
    # Exact Gaussian model B_r - A_r = sum_{m=1}^{c-1} exp(-r*kap*m*(1-m/c)) cos(2pi r j m / c)
    #                                - 2 sum_{m>=1} exp(-r*kap*m*(1+m/c)) cos(2pi r j m / c)
    th = 2.0 * math.pi * ((r * j) % c) / c
    # Small roots m=1..c-1
    B = 0.0
    for m in range(1, c):
        ex = r * kap * m * (1.0 - m / c)
        if ex < 40.0:
            B += math.exp(-ex) * math.cos(m * th)
    # Large roots m>=1
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
    gamma = r * kap
    if gamma > 40.0:
        return 0.0
    th = 2.0 * math.pi * ((r * j) % c) / c
    z = math.exp(-gamma) * complex(math.cos(th), math.sin(th))
    val = z * (1.0 + z) / ((1.0 - z) ** 3)
    return 4.0 * r * (kap / c) * val.real

print("=== CHECK 1A: STRESS TEST G3 BOUNDARY (r0*kap just > 0.25, c/(r0 D^2) large!) ===")
worst_G3_bdy = (0.0, None)
for r0 in [1, 2, 3, 5]:
    for D in [1, 2]:
        for c in [12, 16, 20, 30, 50, 80, 120, 200, 400, 600]:
            # pick j coprime-ish with r0*j - c*p0 = D
            # try p0 = r0 + 1 so j = (c*p0 + D) / r0 if divisible
            found_j = None
            for p0 in range(r0 + 1, 2 * r0 + 2):
                if math.gcd(p0, r0) == 1 and (c * p0 + D) % r0 == 0:
                    jj = (c * p0 + D) // r0
                    if jj > c:
                        found_j = jj
                        break
            if found_j is None:
                continue
            j = found_j
            d = math.gcd(c, j)
            q = c // d
            for r0k in [0.251, 0.28, 0.32, 0.40, 0.55]:
                kap = r0k / r0
                if kap >= 2.5 or c * kap * kap >= 42.0 or c * q * kap <= 0.85:
                    continue
                if d >= 2 and q * q * kap <= 2.0:
                    continue
                # Compute sum_{r>=1} |b_r|/r
                S = 0.0
                rmax = max(30, int(25.0 / kap) + 5)
                for r in range(1, min(rmax, 4000)):
                    S += abs(b_r_gauss(r, c, j, kap)) / r
                if S > worst_G3_bdy[0]:
                    worst_G3_bdy = (S, (c, j, r0, D, kap, r0k, c / (r0 * D * D), c * kap * kap))
                if S > 1.35:
                    print(f"  HIGH at G3 bdy: c={c} j={j} r0={r0} D={D} r0*kap={r0k:.3f} c/(r0 D^2)={c/(r0*D*D):.1f} ck^2={c*kap*kap:.2f} -> S={S:.4f}")
print(f"  WORST at G3 r0*kap>0.25 bdy: S={worst_G3_bdy[0]:.4f} at {worst_G3_bdy[1]}")

print("\n=== CHECK 1B: EXACT MAX OF 4*gamma*|Re z(1+z)/(1-z)^3| FOR gamma >= eps_G3 ===")
# If r0*kap > eps, then for EVERY r = l*r0, gamma = l*r0*kap > l*eps!
# Notice |b_r^main|/r = (4/c) * kap * |Re z(1+z)/(1-z)^3|
# Wait: multiply by c*kap (since c*kap^2 < 42 -> kap/c = (c*kap^2)/(c^2*kap), wait NO!
# If r0*kap >= eps, what is max_{c, theta} (4*gamma/(r*c)) |Re z(1+z)/(1-z)^3| subject to c*kap^2 < 42?!
# Notice 4*gamma/(r*c) = 4*kap/c! And c*kap^2 < 42 -> kap/c = (c*kap^2)/c^2 (small for large c!)
# Wait: for SMALL c, e.g. c=20, kap=0.26 -> c*kap^2 = 1.35, 4*kap/c = 0.052!
# Even better: (4*gamma/(r*c)) * f(gamma, theta):
# For ANY gamma > 0 and theta:
#   f(gamma, theta) = |Re z(1+z)/(1-z)^3| <= min( 0.07681 / (gamma * (theta/(2pi))^2) * ..., e^{-gamma}(1+e^{-gamma})/(1-e^{-gamma})^3 )!
# Compute exact max over theta in [0, pi] of |Re z(1+z)/(1-z)^3| vs theta=0!
for gam in [0.25, 0.30, 0.35, 0.40, 0.50, 0.70, 1.0]:
    mx = 0.0; th_mx = 0.0
    for k in range(2000):
        th = math.pi * k / 2000.0
        z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
        v = abs((z * (1.0 + z) / ((1.0 - z) ** 3)).real)
        if v > mx:
            mx = v; th_mx = th
    v0 = math.exp(-gam) * (1.0 + math.exp(-gam)) / ((1.0 - math.exp(-gam)) ** 3)
    # Also max over c with c*kap^2 < 42 and r=1 -> c < 42/gam^2, wait c >= 2!
    # Wait: if r0*kap >= gam and c/(r0 D^2) >= 15, then c >= 15 r0 D^2 >= 15 r0!
    # So 4*kap/c = 4*(r0*kap)/(r0*c) <= 4*gam / (15 * r0^2) <= 4*gam / 15!!
    print(f"  gam={gam:.2f}: max_th |Re|={mx:8.3f} (at th={th_mx:.3f}, th=0 is {v0:8.3f}) | if c>=15 r0: (4*gam/15)*mx = {(4*gam/15)*mx:.4f}")

print("\n=== CHECK 1C: STRESS TEST c/(r0 D^2) < M_G3 (e.g. c=6..14, r0=1, D=1) ===")
worst_small_M = (0.0, None)
for c in range(6, 20):
    for j in range(c + 1, 2 * c):
        d = math.gcd(c, j)
        q = c // d
        for kap in [0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.35, 0.50, 0.80, 1.20]:
            if c * kap * kap >= 42.0 or c * q * kap <= 0.85:
                continue
            if d >= 2 and q * q * kap <= 2.0:
                continue
            # Check if inside G3 (M_G3=15, eps_G3=0.25):
            in_G3 = False
            for r0 in range(1, c):
                if r0 * kap > 0.25:
                    break
                Dr = min((r0 * j) % c, c - ((r0 * j) % c))
                if Dr > 0 and c >= 15.0 * r0 * Dr * Dr:
                    in_G3 = True
                    break
            if in_G3:
                continue
            S = 0.0
            rmax = max(30, int(25.0 / kap) + 5)
            for r in range(1, min(rmax, 2000)):
                S += abs(b_r_gauss(r, c, j, kap)) / r
            if S > worst_small_M[0]:
                worst_small_M = (S, (c, j, d, q, kap))
print(f"  WORST small-c (6..19) outside G2(q^2k<=2) U G3(15, 0.25): S={worst_small_M[0]:.4f} at {worst_small_M[1]}")

print("\n=== CHECK 2: BOX G2 LARGE-d CORNER (EXACT FINITE-n!) ===")
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
    H_q = sum((math.exp(q*(lr+x))*complex(math.cos(q*ang), math.sin(q*ang)) / ((1 + math.exp(q*(lr+x))*complex(math.cos(q*ang), math.sin(q*ang)))**2)).real / q
              for lr, ang in roots if abs(q*(lr+x)) < 45.0)
    # Also check max |q*delta_psi| among active roots |q*(lr+x)| <= 10!
    max_qdpsi_active = 0.0
    max_qdpsi_all = 0.0
    for lr, ang in roots:
        # wrap q*ang to (-pi, pi]
        qa = ((q * ang + math.pi) % (2.0 * math.pi)) - math.pi
        if abs(q * (lr + x)) <= 10.0:
            max_qdpsi_active = max(max_qdpsi_active, abs(qa))
        if abs(q * (lr + x)) <= 45.0:
            max_qdpsi_all = max(max_qdpsi_all, abs(qa))
    return math.exp(x), H_q, max_qdpsi_active, max_qdpsi_all

for d, q, p0, q2k_target in [
    (50, 1, 2, 0.80), (200, 1, 2, 0.45), (500, 1, 2, 0.28),
    (50, 2, 3, 2.00), (200, 2, 3, 1.20),
    (50, 3, 4, 2.00), (200, 3, 4, 2.00),
    (50, 5, 6, 2.00), (200, 5, 6, 2.00), (1000, 5, 6, 2.00),
    (50, 10, 11, 2.00), (500, 10, 11, 2.00),
]:
    c = d * q; j = d * p0; kt = q2k_target / (q * q)
    if c * kt * kt >= 41.5 or kt >= 2.45:
        continue
    n = max(132, int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt))) // 2 * 2)
    roots, kap = exact_finite_n_roots(n, c, j)
    lam, H_q, qd_act, qd_all = solve_Gq_saddle(roots, c, q)
    print(f"  d={d:4d} q={q:2d} c*k^2={c*kap*kap:5.2f} q^2k={q*q*kap:.3f} lam={lam:.4f} H_q={H_q:6.2f} | max|q*dpsi|(active)={qd_act:.4f} (all)={qd_all:.4f}")
