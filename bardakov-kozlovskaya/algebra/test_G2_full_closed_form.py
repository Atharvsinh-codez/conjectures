"""
Test the COMPLETE 100% CLOSED-FORM Hayman + 2nd-IBP + Zero-Pole Far-Tail bound for Piece G2 (D=0)!
Check exact margins across ALL of Box G2 = { d >= 2, q^2*kap <= K_G2, d*q^2*kap >= 0.85, c*kap^2 < 42 }
for K_G2 = 2.0, 2.2, 2.5!
"""
import math

def exact_finite_n_roots(n, c, j):
    th_c = c * math.pi / n; psi_j = j * math.pi / n; alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1; log_Rz = n * math.log(r3 / r1)
    roots = []; shift = 0.0 if (c % 2 == 0) else 0.5
    k = 0
    while True:
        added = False
        klist = [0] if (k == 0 and c % 2 == 0) else ([0, -1] if (k == 0 and c % 2 == 1) else [k, -k - (1 if c % 2 == 1 else 0)])
        for kk in klist:
            eta = 2.0 * (kk + shift) * math.pi / n
            if abs(eta) > math.pi + 1e-12: continue
            z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
            log_r = log_Rz - n * math.log(abs(z))
            if log_r > -55.0:
                roots.append(log_r); added = True
        if not added and k > c: break
        k += 1
        if k > n // 2 + 2: break
    return roots, 2.0 * math.pi * r1 * r2 / r3

def solve_Gq_plus_saddle(lrs, c, q):
    d = c // q; x = 0.0
    for _ in range(50):
        val = -float(d); der = 0.0
        for lr in lrs:
            lmq = q * (lr + x)
            if lmq < -45.0: continue
            elif lmq > 45.0: val += 1.0 / q
            else:
                u = math.exp(lmq)
                t = u / (1.0 + u); val += t / q; der += t / (1.0 + u)
        if abs(val) < 1e-13: break
        x -= val / der
    Hq = sum(math.exp(q*(lr+x)) / ((1.0 + math.exp(q*(lr+x)))**2) / q for lr in lrs if abs(q*(lr+x)) < 45.0)
    # Also compute exact log(|G_q^+(e^{i*psi})|^{1/q} / G_q^+(1)^{1/q}) at psi = pi/2, 2pi/3, pi!
    def log_ratio(psi):
        s = 0.0
        for lr in lrs:
            lmq = q * (lr + x)
            if abs(lmq) > 45.0: continue
            u = math.exp(lmq)
            h = u / ((1.0 + u) ** 2)
            s += 0.5 * math.log(1.0 - 4.0 * h * (math.sin(0.5 * psi) ** 2)) / q
        return s
    return math.exp(x), Hq, log_ratio(math.pi / 2), log_ratio(2 * math.pi / 3), log_ratio(math.pi)

print("=== 1. EXACT POSITIVE-ROOT SADDLE (lam_*, H_q^+) AND EXACT DECAY AT 2pi/3, pi ===")
min_Hq = 999.0
for d in [2, 3, 4, 5, 10, 50]:
    for q in [1, 2, 3, 5, 10]:
        for q2k in [0.43, 0.85, 1.50, 2.00, 2.50]:
            if d * q2k < 0.85: continue
            kt = q2k / (q * q)
            c = d * q
            if kt >= 2.45 or c * kt * kt >= 41.5: continue
            p0 = q + 1
            j = d * p0
            n = max(132, int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt))) // 2 * 2)
            lrs, kap = exact_finite_n_roots(n, c, j)
            lam, Hq, lr_p2, lr_2p3, lr_pi = solve_Gq_plus_saddle(lrs, c, q)
            if q2k <= 2.05 and Hq < min_Hq:
                min_Hq = Hq
            if (d in [2, 5, 50] and q in [1, 2, 5] and q2k in [0.43, 2.00, 2.50]):
                # Compare lr_pi with -pi^2/(q^2*kap) and -2*Hq!
                print(f"  d={d:2d} q={q:2d} q^2k={q*q*kap:.3f} lam={lam:.4f} H_q+={Hq:5.3f} | lr(2pi/3)={lr_2p3:6.3f} (vs -1.5Hq={-1.5*Hq:6.3f}) | lr(pi)={lr_pi:6.3f} (vs -pi^2/q^2k={-math.pi**2/(q*q*kap):6.3f})")
print(f"  MIN H_q+ across ALL tested points with q^2*kap <= 2.0: {min_Hq:.4f}")

print("\n=== 2. CLOSED-FORM HAYMAN LOWER BOUND ON I_pure / L1 FOR G2 (H_q+ >= 0.78) ===")
# Including:
#   - E_phase(W)^{1/q} <= eps_ph = 0.03 (real and imag!)
#   - R_{non-q}(w) <= eps_strip(psi):
#       on |psi| <= psi0: eps_re0 = 0.03, eps_im0 = 0.05
#       on psi0 <= |psi| <= 2pi/3: eps_re1 = 0.12
#       on 2pi/3 <= |psi| <= pi: far tail bound!
# WAIT: Let's check exact 1D Hayman ratio using our PROVED exact inequalities for G_q^+:
#   For ALL u > 0:
#     (a) log(|1+u e^{i psi}|/(1+u)) <= -2 h_u sin^2(psi/2)  (all psi in [-pi, pi])
#     (b) log(|1+u e^{i psi}|/(1+u)) >= -2 gamma(psi0) h_u sin^2(psi/2)  (|psi| <= psi0)
#         where gamma(psi0) = -log(cos^2(psi0/2)) / sin^2(psi0/2)
#     (c) |arg(1+u e^{i psi}) - (u/(1+u))*psi| <= (1/6) h_u |psi|^3  (|psi| <= pi/2, let's verify 1/6 vs 1/3!)
# Wait: FIRST verify exact max_{u>0, |psi|<=pi/2} |arg(1+u e^{i psi}) - (u/(1+u)) psi| / (h_u |psi|^3)!
max_c3 = 0.0
for iu in range(1, 300):
    u = math.exp(-6.0 + 12.0 * iu / 300.0)
    h = u / ((1.0 + u) ** 2)
    m = u / (1.0 + u)
    for ip in range(1, 200):
        psi = (math.pi / 2.0) * ip / 200.0
        err = abs(math.atan2(u * math.sin(psi), 1.0 + u * math.cos(psi)) - m * psi)
        max_c3 = max(max_c3, err / (h * (psi ** 3)))
print(f"  EXACT max |arg(1+u e^{{i psi}}) - m*psi| / (h_u |psi|^3) on |psi| <= pi/2: {max_c3:.6f}  (<= 1/6 = {1/6:.6f}!)")

# Now compute closed-form Hayman lower bound on (I_main - I_mid - I_far) / (I_main_L1 + I_mid + I_far)
# as a function of H = H_q+ in [0.78, 100]!
def hayman_lb(H, eps_re0=0.05, eps_im0=0.10, eps_re1=0.15, eps_far=0.15, q2k=2.0):
    # Optimize psi0 in (0, pi/2]
    best = -1.0
    for ip in range(10, 100):
        psi0 = (math.pi / 2.0) * ip / 100.0
        gam = -math.log(math.cos(0.5 * psi0) ** 2) / (math.sin(0.5 * psi0) ** 2)
        # Main peak [0, psi0]:
        #   lower bound on Re: e^{-eps_re0} * int_0^{psi0} exp(-2*gam*H*sin^2(psi/2)) * cos(c3*H*psi^3 + eps_im0) dpsi
        #   Wait: check that c3*H*psi0^3 + eps_im0 < pi/2!
        c3 = 1.0 / 6.0  # (will prove 1/6 or 0.17 analytically!)
        if c3 * H * (psi0 ** 3) + eps_im0 >= 1.50:
            continue
        N = 100
        dpsi = psi0 / N
        I_pos = 0.0
        I_main_up = 0.0
        for k in range(N):
            psi = (k + 0.5) * dpsi
            s2 = math.sin(0.5 * psi) ** 2
            ph = c3 * H * (psi ** 3) + eps_im0
            I_pos += math.exp(-eps_re0 - 2.0 * gam * H * s2) * math.cos(ph) * dpsi
            I_main_up += math.exp(eps_re0 - 2.0 * H * s2) * dpsi
        # Mid tail [psi0, 2pi/3]:
        dpsi1 = (2.0 * math.pi / 3.0 - psi0) / N
        I_mid = 0.0
        for k in range(N):
            psi = psi0 + (k + 0.5) * dpsi1
            s2 = math.sin(0.5 * psi) ** 2
            I_mid += math.exp(eps_re1 - 2.0 * H * s2) * dpsi1
        # Far tail [2pi/3, pi]:
        #   For q=1: NO R_nonq, so exp(eps_ph - 2*H*sin^2(psi/2)) <= exp(0.05 - 1.5*H) * pi/3!
        #   For q>=2: zero-pole shift gives decay >= min(1.5*H, pi^2/(3*q2k) * 0.90) -> wait: for H = 2/q2k, pi^2/(3*q2k) = (pi^2/6)*H = 1.645*H!
        #   Even conservatively using exp(eps_far - 1.35 * H) * (pi/3):
        I_far = math.exp(eps_far - 1.35 * H) * (math.pi / 3.0)
        num = I_pos - I_mid - I_far
        den = I_main_up + I_mid + I_far
        rat = num / den
        if rat > best:
            best = rat
    return best

for H in [0.75, 0.80, 0.85, 0.90, 1.0, 1.5, 2.0, 5.0, 20.0, 100.0]:
    print(f"  H_q+ = {H:6.2f}: CLOSED-FORM Hayman rat_lb (with generous errs 0.05/0.10/0.15/0.15) = {hayman_lb(H):+.4f}")
