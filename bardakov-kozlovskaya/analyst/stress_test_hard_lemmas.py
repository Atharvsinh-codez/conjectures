"""
STRESS TEST HARD LEMMAS A, B, C BEFORE WRITING PROOFS:
1. Check Hard Lemma B (D = 0, q >= 2) for 0.14 < c*q*kap <= 2.0 and g = 1..10:
   Compare exact discrete F_pure(lam e^{i psi}) vs the positive-root block product prod_s (1 + B_s W)!
   Specifically: what is the EXACT error between log F_pure(lam e^{i psi}) and the q-periodic sum
   of positive real roots (via our exact Fourier split: q|r vs q!|r)?
2. Check Hard Lemma C (D = 1, small c_tilde = 30, 50, 100 and tau = 0.4..5.0):
   Specifically test c=30 j=31 tau=1.5, c=50 j=51 tau=1.5, c=50 j=51 tau=0.8, c=100 j=101 tau=0.8:
   Find the EXACT explicit radius lam(c_tilde, D, tau) where H > 0 AND rat_L1 > 0.80 AND lam stays
   safely separated from e^{-tau/4}!
3. Check the EXACT small-c gap 0.14 < c*q*kap <= 1.0 for c = 4, 5, 6, 8, 10!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd

def eval_exact_model_info(c, j, kap, lam, psi0, M_psi=512):
    alpha = kap / c
    M_cut = int(np.ceil(sqrt(45.0 / max(alpha, 1e-12)))) + 10
    ms = np.arange(-(c + M_cut), M_cut + 1, dtype=np.float64)
    re_a = -alpha * ms * (c + ms) + log(lam)
    im_a = -2.0 * pi * ms * j / c
    keep = re_a > -35.0
    re_a = re_a[keep]; im_a = im_a[keep]
    big = re_a > 0.0

    # Curvature H and g1 at psi0:
    g1 = -c + 0.0j; H = 0.0j
    if np.any(big):
        iz = np.exp(-re_a[big] - 1j*(im_a[big] + psi0))
        g1 += np.sum(1.0 + iz / (1.0 - iz))
        H += np.sum(-iz / ((1.0 - iz)**2))
    sm = ~big
    if np.any(sm):
        z = np.exp(re_a[sm] + 1j*(im_a[sm] + psi0))
        g1 += np.sum(-z / (1.0 - z))
        H += np.sum(-z / ((1.0 - z)**2))

    # Exact integral rat_L1 on |w| = lam:
    psis = np.linspace(-pi, pi, M_psi, endpoint=False)
    log_F = -c * (log(lam) + 1j * psis) + (1j * pi if (j*(c+1)+c)%2 == 1 else 0.0)
    if np.any(big):
        rb = re_a[big]; ib = im_a[big]
        log_F = log_F + np.sum(rb) + 1j * (np.sum(ib) + len(rb)*pi) + len(rb)*1j*psis
        iu = np.exp(-rb - 1j*ib)
        epm = np.exp(-1j*psis)
        for k0 in range(0, len(iu), 250):
            ch = iu[k0:k0+250]
            log_F += np.sum(np.log(1.0 - ch[None, :]*epm[:, None]), axis=1)
    if np.any(sm):
        rs = re_a[sm]; is_ = im_a[sm]
        u = np.exp(rs + 1j*is_)
        epp = np.exp(1j*psis)
        for k0 in range(0, len(u), 250):
            ch = u[k0:k0+250]
            near1 = np.abs(ch) > 1.0 - 1e-9
            if np.any(near1):
                p1 = np.prod(1.0 - ch[near1][None, :]*epp[:, None], axis=1)
                with np.errstate(divide='ignore'):
                    log_F += np.where(np.abs(p1) > 0, np.log(p1), -1e9)
            far = ~near1
            if np.any(far):
                log_F += np.sum(np.log(1.0 - ch[far][None, :]*epp[:, None]), axis=1)
    sh = float(np.max(log_F.real))
    F = np.exp(log_F - sh)
    rat = float(np.mean(F).real) / float(np.mean(np.abs(F)))
    active = (log_F.real - sh) > -4.6
    aa = float(np.max(np.abs(np.angle(F[active]))))
    return g1.real, H.real, rat, aa

print("=== 1. CHECK HARD LEMMA C: D = 1 AT SMALL c = 30, 50, 100 AND tau = 0.8..3.0 ===")
for c, tau in [(30, 1.5), (30, 3.0), (50, 1.5), (50, 2.0), (50, 3.0), (100, 0.8), (100, 1.5), (100, 2.0)]:
    j = c + 1; kap = tau / c; em = exp(-tau / 4.0)
    print(f"\n  c={c}, j={j}, tau={tau:.1f}, e^(-tau/4)={em:.4f}:")
    for mult in [1.05, 1.12, 1.20, 1.30, 1.40]:
        lam = min(1.0, em * mult)
        g1, H, rat, aa = eval_exact_model_info(c, j, kap, lam, pi)
        print(f"    lam={lam:.4f} (lam/em={lam/em:.3f}): g1={g1:+7.3f} H={H:+8.2f} | rat={rat:.5f} argAct={aa:.4f}")
    g1, H, rat, aa = eval_exact_model_info(c, j, kap, 1.0, pi)
    print(f"    lam=1.0000 (lam/em={1.0/em:.3f}): g1={g1:+7.3f} H={H:+8.2f} | rat={rat:.5f} argAct={aa:.4f}")

print("\n=== 2. CHECK SMALL-c GAP (0.14 < c*q*kap <= 1.0) FOR c = 4, 5, 6, 8 ===")
for c, j in [(4, 5), (4, 6), (5, 6), (6, 7), (6, 9), (8, 9), (8, 12)]:
    q = c // gcd(c, j)
    print(f"  c={c} j={j} q={q}:")
    for cqk in [0.15, 0.25, 0.40, 0.65, 1.00]:
        kap = cqk / (c * q)
        # Test at lam = 1.0 AND at block-saddle lam_block!
        g = c // q; alpha = kap / c
        # Block saddle:
        M = int(np.ceil(c/(2*q) + sqrt(50.0/max(q*alpha, 1e-9))/q)) + 5
        ss = np.arange(-M, M+1, dtype=np.float64)
        mid = q*ss + 0.5*(q - 1)
        ex = -alpha*q*(mid*(c + mid) + (q*q - 1.0)/12.0)
        lo, hi = -20.0, 10.0
        for _ in range(40):
            vm = 0.5*(lo + hi)
            z = np.clip(ex + q*vm, -60.0, 60.0)
            if float(np.sum(1.0/(1.0 + np.exp(-z)))) > g: hi = vm
            else: lo = vm
        lam_b = exp(0.5*(lo + hi))
        _, _, rat1, aa1 = eval_exact_model_info(c, j, kap, 1.0, pi/q, M_psi=256)
        g1b, Hb, ratb, aab = eval_exact_model_info(c, j, kap, lam_b, pi/q, M_psi=256)
        print(f"    cqk={cqk:.2f} (kap={kap:.4f}): @1 rat={rat1:.4f} aa={aa1:.3f} | @lam_b={lam_b:.4f} g1={g1b:+.2f} H={Hb:+.2f} rat={ratb:.4f} aa={aab:.3f}")
