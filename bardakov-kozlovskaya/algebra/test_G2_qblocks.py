"""
Test Method B (Exact q-Block Product Comparison) on EXACT FINITE-n D = 0 Major Arcs (d >= 2, q >= 2):
  Group roots x = k + shift into consecutive blocks B_s of length q:
    P_s(w) = prod_{t=0}^{q-1} (1 - a_{q*s+t} w)   vs   Q_s(w) = 1 - A_s w^q,   A_s = prod_{t=0}^{q-1} a_{q*s+t}.
Check:
  1. For each block s, how small is |P_s(w) / Q_s(w) - 1| on the peak |psi| <= psi0?
  2. What is the SUM sum_s |P_s(w)/Q_s(w) - 1| across ALL blocks s?
  3. On the ENTIRE circle (including tail |psi| > psi0), how does |F(w)| compare to |prod_s Q_s(w)| (NO zero-order mismatch!)?
"""
import math
import cmath

def exact_roots_indexed(n, c, j):
    """Return dict x_idx -> (log|a|, arg(a)) for k in [-K, K], where x = k + shift."""
    th_c = c * math.pi / n; psi_j = j * math.pi / n; alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1; log_Rz = n * math.log(r3 / r1)
    kap = 2.0 * math.pi * r1 * r2 / r3
    shift = 0.0 if (c % 2 == 0) else 0.5
    sgn_ang = 0.0 if (j % 2 == 0) else math.pi
    K = int(math.sqrt(0.25 * c * c + 50.0 * c / kap)) + 20
    K = min(K, n // 2)
    out = {}
    for k in range(-K, K + 1):
        eta = 2.0 * (k + shift) * math.pi / n
        if abs(eta) > math.pi + 1e-12: continue
        z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
        lr = log_Rz - n * math.log(abs(z))
        ang = sgn_ang - n * math.atan2(z.imag, z.real)
        out[k] = (lr, ang)
    return out, kap, K

def test_q_blocks(d, q, p0, kt):
    c = d * q; j = d * p0
    n = int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt)))
    n = max(132, n + (n % 2))
    rmap, kap, K = exact_roots_indexed(n, c, j)
    # Partition k into blocks of length q starting at k0 = -((K // q) * q)
    # Better: align blocks symmetrically around k = -c/2 (or k = 0)!
    # If c even: k = -c/2 is the left boundary root! Starting blocks at k = -c/2 mod q aligns boundary roots cleanly!
    # Actually test centered offset off so block centers are symmetric around x = 0:
    # For even c (shift=0): x = k. Symmetric block of length q has center 0 (odd q: k in {-(q-1)/2..+(q-1)/2})
    # or center 0.5 (even q: k in {-q/2+1..q/2}).
    # For odd c (shift=0.5): q must be odd, x = k + 0.5, center 0.5 -> k in {-(q-1)/2..+(q-1)/2}!
    off = -(q // 2) if (q % 2 == 1) else (-(q // 2) + 1)
    s_min = (-K - off) // q + 1
    s_max = (K - off) // q - 1
    blocks = []
    for s_idx in range(s_min, s_max + 1):
        ks = [off + s_idx * q + t for t in range(q)]
        if all(k in rmap for k in ks):
            blocks.append([rmap[k] for k in ks])

    # Solve exact saddle Lam = lam^q for Q(W) = prod_s (1 + A_s W) (W = -w^q):
    x = 0.0
    for _ in range(40):
        val = -float(d); der = 0.0
        for blk in blocks:
            lr_sum = sum(lr for lr, _ in blk) + q * x
            ang_sum = sum(ang for _, ang in blk)
            if lr_sum < -45.0: continue
            elif lr_sum > 45.0: val += 1.0
            else:
                u = math.exp(lr_sum) * complex(math.cos(ang_sum), math.sin(ang_sum))
                t = u / (1.0 + u); val += t.real; der += (t / (1.0 + u)).real * q
        if abs(val) < 1e-12 or der <= 0: break
        x -= max(-0.2, min(0.2, val / der))
    lam = math.exp(x)
    # Compute H_Q (curvature of Q in psi = q*phi - pi):
    H_Q = 0.0
    for blk in blocks:
        lr_sum = sum(lr for lr, _ in blk) + q * x
        ang_sum = sum(ang for _, ang in blk)
        if abs(lr_sum) > 45.0: continue
        u = math.exp(lr_sum) * complex(math.cos(ang_sum), math.sin(ang_sum))
        H_Q += (u / ((1.0 + u) ** 2)).real

    # Evaluate ratio prod_s (P_s(w) / Q_s(w)) at the peak phi0 = pi/q (psi = 0) AND on |psi| <= pi/2!
    def ratio_at_phi(phi):
        prod_rat = 1.0 + 0j
        sum_dev = 0.0
        for blk in blocks:
            lr_sum = sum(lr for lr, _ in blk) + q * x
            if lr_sum < -45.0: continue
            if lr_sum > 45.0:
                # P_s / Q_s = prod(1 - 1/(a w)) / (1 - 1/(A w^q)) ~ 1
                ps = 1.0 + 0j
                for lr, ang in blk:
                    lm = lr + x
                    if lm < 45.0:
                        u = math.exp(-lm) * cmath.exp(-1j * (ang + phi))
                        ps *= (1.0 - u)
                # wait: careful sign! (-a w) * (1 - 1/(a w)) -> prod(-a w) = (-1)^q A w^q;
                # Q_s = 1 - A w^q = -A w^q (1 - 1/(A w^q))!
                # Check sign of prod(-a_t w) / (-A w^q) = (-1)^{q-1}!!
                # Wait: for even q, prod(1 - u_t) = 1 - ... + (prod u_t), while for odd q it is 1 - ... - (prod u_t)!
                # OH!! CATCH!! prod_{t=0}^{q-1} (1 - u_t) has top term (-1)^q * (prod u_t),
                # whereas prod_{t=0}^{q-1} (1 - R*zeta^t*w) = 1 - R^q w^q BECAUSE prod_{t=0}^{q-1} zeta^t = (-1)^{q-1}!!
                # SO A_s := (-1)^{q-1} * prod a_t ALWAYS HAS TOP TERM -A_s w^q!!
                continue
            ps = 1.0 + 0j
            for lr, ang in blk:
                u = math.exp(lr + x) * cmath.exp(1j * (ang + phi))
                ps *= (1.0 - u)
            ang_sum = sum(ang for _, ang in blk)
            # Exact sign of prod(a_t): remember prod_{t=0}^{q-1} e^{2pi i t p0 / q} = (-1)^{p0*(q-1)} = (-1)^{q-1} (if q even then p0 odd!)
            # So (-1)^q * prod(a_t) = - (-1)^{q-1} prod(a_t)!
            # Check: in Q_s = 1 + (-1)^q * (prod a_t) * w^q = 1 - ((-1)^{q-1} prod a_t) w^q!
            qs = 1.0 + ((-1) ** q) * math.exp(lr_sum) * cmath.exp(1j * (ang_sum + q * phi))
            r_val = ps / qs
            prod_rat *= r_val
            sum_dev += abs(r_val - 1.0)
        return prod_rat, sum_dev

    # Check exact sign of (-1)^{q-1} * prod(a_t) at psi = 0!
    phi0 = math.pi / q
    r0, sd0 = ratio_at_phi(phi0)
    max_dev_half = 0.0
    max_log_rat_all = -1e9
    for k_i in range(64):
        psi = -math.pi + 2.0 * math.pi * (k_i + 0.5) / 64
        phi = (math.pi + psi) / q
        rv, sd = ratio_at_phi(phi)
        if abs(psi) <= math.pi / 2.0:
            max_dev_half = max(max_dev_half, abs(rv - 1.0))
        max_log_rat_all = max(max_log_rat_all, math.log(abs(rv)))
    print(f"d={d:2d} q={q:2d} kap={kap:.4f} dq^2k={d*q*q*kap:5.2f} dq^4k^2={d*(q*q*kap)**2:5.2f} | "
          f"lam={lam:.4f} H_Q={H_Q:5.2f} | |F/Q-1|(pk)={abs(r0-1):.4f} sum_dev(pk)={sd0:.4f} "
          f"|F/Q-1|(|psi|<=pi/2)={max_dev_half:.4f} max_log|F/Q|_all={max_log_rat_all:+.4f}")

print("=== METHOD B: EXACT q-BLOCK COMPARISON ON D = 0 MAJOR ARCS (d >= 2, q >= 2) ===")
for args in [
    (2, 2, 3, 0.25), (2, 3, 4, 0.10), (2, 5, 6, 0.04), (2, 10, 11, 0.01),
    (3, 2, 3, 0.20), (3, 4, 5, 0.05),
    (5, 2, 3, 0.25), (5, 2, 3, 0.50), (5, 3, 4, 0.15), (5, 5, 6, 0.06), (5, 7, 9, 0.03), (5, 10, 11, 0.015),
    (10, 2, 3, 0.15), (10, 3, 4, 0.10), (15, 5, 7, 0.04), (30, 3, 4, 0.03), (50, 2, 3, 0.03),
]:
    test_q_blocks(*args)
