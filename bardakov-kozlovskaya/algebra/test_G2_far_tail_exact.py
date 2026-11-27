"""
Verify the EXACT pointwise & integral far-tail bound for Piece G2 (q >= 2, |psi| in [2pi/3, pi])!
Check the ratio:
  integral_{-n/2}^{n/2} (1/q) * log( P_mod(|b(x)|*lam_*, psi) / (1 + (|b(x)|*lam_*)^q) ) dx
vs -c_far * H_q+ across ALL corners of Box G2 (d >= 2, q >= 2, q^2*kap <= 2.0, d*q^2*kap >= 0.85)!
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
    return math.exp(x), Hq

# Check the EXACT universal 1-variable function:
# For ANY v in R, q >= 2, beta in [0, pi/3] (psi = pi - beta):
#   f(v; q, beta) := log( P_mod(e^v; pi-beta) / (1 + e^{qv}) )
#                  = 0.5*log(1 - 4*h(qv)*cos^2(beta/2)) + 0.5*log( (sinh^2(v/2) + sin^2(pi/(6q))) / (sinh^2(v/2) + sin^2(beta/(2q))) )
# Wait: compare f(v; q, beta) with -c0 * h(qv) = -c0 / (4*cosh^2(qv/2))!
# Wait: as v -> inf, f(v) > 0 (O(q^{-2} e^{-v})), whereas on the integral over x:
print("=== 1. CHECK EXACT FAR-TAIL INTEGRAL RATIO (-log_mod / H_q+) ON ALL BOX G2 CORNERS ===")
worst_ratio = (999.0, None)
for d in [2, 3, 5, 10, 50]:
    for q in [2, 3, 5, 10]:
        for q2k in [0.43, 0.50, 0.85, 1.40, 2.00]:
            if d * q2k < 0.85: continue
            kt = q2k / (q * q); c = d * q
            if kt >= 2.45 or c * kt * kt >= 41.5: continue
            j = d * (q + 1)
            n = max(132, int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt))) // 2 * 2)
            lrs, kap = exact_finite_n_roots(n, c, j)
            lam, Hq = solve_Gq_plus_saddle(lrs, c, q)
            log_lam = math.log(lam)
            # Compute sum_k (1/q) * f(lr_k + log_lam; q, beta) over beta in [0, pi/3]
            for beta in [0.0, math.pi / 6.0, math.pi / 3.0]:
                tot = 0.0
                for lr in lrs:
                    v = lr + log_lam
                    qv = q * v
                    if abs(qv) > 45.0:
                        # still add the small positive shift term 0.5*log((sinh^2(v/2)+sin^2(pi/6q))/(sinh^2(v/2)+sin^2(beta/2q)))!
                        if abs(v) < 30.0:
                            sh2 = math.sinh(0.5 * v) ** 2
                            s1 = math.sin(math.pi / (6.0 * q)) ** 2
                            s0 = math.sin(beta / (2.0 * q)) ** 2
                            tot += 0.5 * math.log((sh2 + s1) / (sh2 + s0)) / q
                        continue
                    # Exact: P_mod(e^v; pi-beta) / (1 + e^{qv})
                    # = |1 - e^{qv - i*beta}| / (1 + e^{qv}) * |1 - e^{v + i*pi/(3q)}| / |1 - e^{v + i*beta/q}|
                    # At beta = 0, use limit |1 - e^{qv}|/|1 - e^v| = sum_{k=0}^{q-1} e^{kv} = sinh(qv/2)/sinh(v/2) * e^{(q-1)v/2}!
                    sh_q = math.sinh(0.5 * qv)
                    ch_q = math.cosh(0.5 * qv)
                    sh_1 = math.sinh(0.5 * v)
                    num1 = sh_q * sh_q + math.sin(0.5 * beta) ** 2
                    den1 = sh_1 * sh_1 + math.sin(0.5 * beta / q) ** 2
                    if abs(v) < 1e-8 and beta == 0.0:
                        rat_sq = (q * math.sin(math.pi / (6.0 * q))) ** 2
                    else:
                        rat_sq = (num1 / den1) * (sh_1 * sh_1 + math.sin(math.pi / (6.0 * q)) ** 2) / (ch_q * ch_q)
                    tot += 0.5 * math.log(rat_sq) / q
                ratio = -tot / Hq
                if ratio < worst_ratio[0]:
                    worst_ratio = (ratio, (d, q, round(q*q*kap, 3), round(beta, 3), round(tot, 4), round(Hq, 4)))
print(f"  WORST FAR-TAIL DECAY RATIO (-log_mod / H_q+) across ALL Box G2: {worst_ratio[0]:.4f} at {worst_ratio[1]}")

# Also check: what is the MAXIMUM of log_mod(v; q, beta) itself, OR better:
# Wait: compare f(v; q, beta) vs f(v; q, pi/3) = 0.5*log(1 - 3*h(qv))!
# WAIT A SECOND!! CHECK IF f(v; q, beta) <= f(v; q, pi/3) + something TINY, OR check
# EXACT pointwise max over v, q >= 2, beta in [0, pi/3] of f(v; q, beta) + c_pt * h(qv)!!
print("\n=== 2. CHECK IF P_mod(e^v; pi-beta) / (1+e^{qv}) HAS AN EVEN SIMPLER MONOTONICITY IN beta! ===")
# Check for each v, is f(v; q, beta) <= 0.5*log(1 - 3*h(qv))? (NO, at large v, beta=0 has positive shift!)
# Instead check what c_far gives: sum_{v} f(v) <= -c_far * sum_{v} h(qv) + tiny!
# Better: check if we replace theta_{t0} by pi/(3q) ONLY FOR |v| <= v_cut (NO! That breaks analyticity for strip-Poisson!).
# Instead: check how big the positive tail (where f(v) > 0) can EVER be!
for q in [2, 3, 5, 10, 50]:
    # Integrate max(0, f(v; q, 0)) over v in [0, inf) (times 2 for both sides of v=0):
    # Wait: f(v; q, 0) = log( (sinh(qv/2)/(cosh(qv/2)*sinh(v/2))) * sqrt(sinh^2(v/2) + sin^2(pi/6q)) )
    # Notice u = q*v -> (1/q) int f(u/q) du/q!
    pos_int = 0.0
    tot_int = 0.0
    dv = 0.005 / q
    for k in range(1, 8000):
        v = (k - 0.5) * dv
        sh_q = math.sinh(0.5 * q * v)
        ch_q = math.cosh(0.5 * q * v)
        sh_1 = math.sinh(0.5 * v)
        val = math.log((sh_q / (ch_q * sh_1)) * math.sqrt(sh_1 * sh_1 + math.sin(math.pi / (6.0 * q)) ** 2))
        tot_int += 2.0 * val * dv
        if val > 0:
            pos_int += 2.0 * val * dv
    # Remember exact tot_int should be -(pi^2 / (3*q))? Wait: 2/q times? Let's check q * tot_int!
    print(f"  q={q:2d}: q * int_{{-inf}}^{{inf}} f(v;q,0) dv = {q*tot_int:+.5f} (theory -pi^2/6 = {-math.pi**2/6:+.5f}? Wait -pi^2/3={-math.pi**2/3:+.5f}),  q * pos_part = {q*pos_int:+.5f}")
