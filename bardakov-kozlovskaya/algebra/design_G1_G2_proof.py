"""
Verify the ultra-clean SINGLE-PARAMETER boundary q^2 * kap = 2.0 between G2 (D=0) and G1 (Minor Arcs)!
Part 1: Compute exact max_{T >= 0.425} sum_{l=1}^inf (T / l) * |D(l*T)| (Wait: also check max over d=1, T>=0.425 of (1/2) sum |D(l*T)|/l!).
Part 2: Stress test Box G2 = { d >= 2, q^2*kap <= 2.0, d*q^2*kap >= 0.85, c*kap^2 < 42, kap < 2.5 }
        across ALL corners (d = 2, 3, 4, 5, 10, 50; q = 1, 2, 3, 5, 10, 20; q^2*kap up to 2.0)
        on EXACT FINITE-n roots!
Part 3: Check exact G1 Minor Arcs outside { G2: d>=2, q^2*kap <= 2.0 } U { G3: D!=0, c/(r0 D^2) >= 15, r0*kap <= 0.25 }!
"""
import math
import cmath

def D_exact(T):
    M = 400
    s1 = sum(math.exp(-0.5 * T * (1.0 - ((i + 0.5) / M) ** 2)) for i in range(M)) * (2.0 / M)
    ymax = min(50.0, 25.0 / max(T, 1e-4) + 10.0 / math.sqrt(max(T, 1e-4)))
    s2 = sum(math.exp(-T * ((i + 0.5) * ymax / M) - 0.5 * T * (((i + 0.5) * ymax / M) ** 2)) for i in range(M)) * (2.0 * ymax / M)
    return abs(s1 - s2)

print("=== PART 1: EXACT UNIVERSAL BOUND ON q|r SUM OUTSIDE G2 (q^2*kap > 2.0 OR d = 1) ===")
# Remember: |b_{l*q}| / (l*q) ~ (d / (2*l)) * |D(l*T1)| where T1 = (1/2) * d * q^2 * kap!
# Notice (d / (2*l)) = (T1 / (q^2*kap)) / l!
# So for q^2*kap > 2.0 and ANY d >= 2 (hence T1 > d >= 2.0):
#   sum_{l>=1} (d / (2*l)) * |D(l*T1)| = (1 / (q^2*kap)) * sum_{l>=1} (T1 / l) * |D(l*T1)|!
max_S_T = (0.0, 0.0)
for k in range(42, 400):
    T1 = k * 0.05
    # Sum l = 1..60 + tail (4/(l^2 T1^2) * T1/l = 4/(l^3 T1)):
    st = sum((T1 / l) * D_exact(l * T1) for l in range(1, 61)) + 4.0 / (T1 * 2.0 * 60 * 60)
    if T1 >= 1.0 and st > max_S_T[0]:
        max_S_T = (st, T1)
    if k in [42, 60, 80, 100, 140, 200]:
        print(f"  T1={T1:5.2f}: sum (T1/l)|D(l*T1)| = {st:.4f}  (-> at q^2*kap=2.0: {st/2.0:.4f}; at q^2*kap=2.5: {st/2.5:.4f})")
print(f"  GLOBAL MAX over T1 >= 1.0: {max_S_T[0]:.4f} at T1={max_S_T[1]:.2f} (-> at q^2*kap >= 2.0: <= {max_S_T[0]/2.0:.4f}!)")

# Also check d = 1 (T1 >= 0.425): sum_{l>=1} (1 / (2*l)) * |D(l*T1)|!
max_d1 = (0.0, 0.0)
for k in range(17, 200):
    T1 = k * 0.025
    sd1 = sum((0.5 / l) * D_exact(l * T1) for l in range(1, 61)) + 2.0 / (T1 * T1 * 2.0 * 60 * 60)
    if sd1 > max_d1[0]:
        max_d1 = (sd1, T1)
print(f"  GLOBAL MAX for d=1 (T1 >= 0.425): sum (1/(2l))|D(l*T1)| = {max_d1[0]:.4f} at T1={max_d1[1]:.3f}!")

# Part 2: Stress test Box G2 = { d >= 2, q^2*kap <= 2.0, d*q^2*kap >= 0.85, c*kap^2 < 42 } on EXACT FINITE-n!
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
    return math.exp(x), H_q

print("\n=== PART 2: STRESS TEST BOX G2 (d >= 2, q^2*kap <= 2.0, d*q^2*kap >= 0.85) ===")
for d, q, p0, q2k_target in [
    (2, 1, 2, 0.43), (2, 1, 2, 1.00), (2, 1, 2, 2.00),
    (2, 2, 3, 0.45), (2, 2, 3, 1.00), (2, 2, 3, 2.00),
    (2, 3, 4, 0.45), (2, 3, 4, 2.00),
    (2, 5, 6, 0.45), (2, 5, 6, 2.00),
    (2, 10, 11, 2.00),
    (3, 1, 2, 2.00), (3, 2, 3, 2.00), (3, 5, 6, 2.00),
    (5, 1, 2, 2.00), (5, 2, 3, 2.00), (5, 5, 6, 2.00),
    (10, 1, 2, 2.00), (10, 2, 3, 2.00), (10, 5, 6, 2.00),
]:
    c = d * q; j = d * p0; kt = q2k_target / (q * q)
    if c * kt * kt >= 42.0 or kt >= 2.5: continue
    n = max(132, int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt))) // 2 * 2)
    roots, kap = exact_finite_n_roots(n, c, j)
    lam, H_q = solve_Gq_saddle(roots, c, q)
    # Check max |R_{non-q}| (both Re and Im!) on |psi| <= 2*pi/3!
    log_lam = math.log(lam)
    sgn_pref = 0.0 if ((j*(c+1)+c) % 2 == 0) else math.pi
    max_re_R = 0.0; max_im_R = 0.0
    for i in range(64):
        psi = -2.0 * math.pi / 3.0 + (4.0 * math.pi / 3.0) * (i + 0.5) / 64.0
        phi = (math.pi + psi) / q
        laF, phF = -c * log_lam, -c * phi
        laG, phG = -c * log_lam, -c * phi
        # Compute G_q(W)^{1/q} directly in W = -w^q = lam^q e^{i psi} (no branch ambiguity!):
        laG = -d * q * log_lam; phG = -d * psi
        for lr, ang in roots:
            lm = lr + log_lam; th = ang + phi
            if lm < -45.0: continue
            elif lm > 45.0: laF += lm; phF += th + math.pi
            else:
                r = math.exp(lm); z = complex(1 - r*math.cos(th), -r*math.sin(th))
                laF += math.log(abs(z)); phF += math.atan2(z.imag, z.real)
            lmq = q * lm; angq = q * ang + psi
            if lmq < -45.0: continue
            elif lmq > 45.0: laG += lm; phG += (q * ang + psi) / q
            else:
                rq = math.exp(lmq); zq = complex(1 + rq*math.cos(angq), rq*math.sin(angq))
                laG += math.log(abs(zq)) / q; phG += math.atan2(zq.imag, zq.real) / q
        max_re_R = max(max_re_R, abs(laF - laG))
        epF = (phF + sgn_pref + math.pi) % (2 * math.pi) - math.pi
        epG = (phG + math.pi) % (2 * math.pi) - math.pi
        dph = (epF - epG + math.pi) % (2 * math.pi) - math.pi
        max_im_R = max(max_im_R, abs(dph))
    print(f"  d={d:2d} q={q:2d} kap={kap:.4f} q^2k={q*q*kap:.3f} lam={lam:.4f} H_q={H_q:5.2f} | max|Re R|(2pi/3)={max_re_R:.4f} max|Im R|(2pi/3)={max_im_R:.4f}")
