"""
TEST THE EXACT r_0-MULTISECTION IDENTITY + 3-LINE POLYA-SZEGO POSITIVE-ROOT THEOREM!
For ANY r_0 >= 1 and ANY lam > 0 (avoiding roots):
  log F_pure(lam e^{i psi}) = (1/r_0) log F_main(lam^{r_0} e^{i r_0 psi}) + R_{non-r_0}(lam e^{i psi}),
where F_main(W) = (-1)^{...} W^{-c} prod_k (1 - a_k^{r_0} W) is EXACT (ZERO block-averaging error!).

Test:
1. D = 0 (r_0 = q >= 2): across ALL c in [4, 200], g = c/q in [1..10], c*q*kap in [0.14, 20]:
   - Compute exact H_q at exact saddle Lambda_* = lam_*^q!
   - Compute exact max|Re R_{non-q}| and max|Im R_{non-q}| on |w| = lam_*!
   - Check if our 3-line Polya-Szego theorem + R_{non-q} IMMEDIATELY proves positivity!
2. Also check g = 1 and g = 2 direct power-sum expansion at w=0 for 0.14 < c*q*kap <= 2.5!
3. Check the EXACT finite-n R_{non-q} and phase error on real (n, c, j) points!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd

def analyze_D0_multisection(c, j, kap, M_psi=512):
    g = gcd(c, j)
    q = c // g
    alpha = kap / c
    M_cut = int(np.ceil(sqrt(45.0 / max(alpha, 1e-12)))) + 10
    ms = np.arange(-(c + M_cut), M_cut + 1, dtype=np.float64)
    re_a = -alpha * ms * (c + ms) # log |a_m|
    im_a = -2.0 * pi * ms * j / c
    keep = re_a > -40.0
    re_a = re_a[keep]; im_a = im_a[keep]

    # Solve EXACT positive-root saddle for W = -Lambda_* = -lam_*^q:
    #   sum_m 1 / (1 + exp(-(q * re_a[m] + V))) = c,  where V = q * log(lam_*)!
    lo, hi = -80.0, 40.0
    for _ in range(60):
        vm = 0.5 * (lo + hi)
        z = np.clip(q * re_a + vm, -70.0, 70.0)
        if float(np.sum(1.0 / (1.0 + np.exp(-z)))) > c:
            hi = vm
        else:
            lo = vm
    V = 0.5 * (lo + hi)
    lam_star = exp(V / q)
    z = np.clip(q * re_a + V, -70.0, 70.0)
    p = 1.0 / (1.0 + np.exp(-z))
    H_q = float(np.sum(p * (1.0 - p))) / q

    # Check if lam_star hits any root: min |re_a + log(lam_star)|
    min_gap = float(np.min(np.abs(re_a + V / q)))

    # Evaluate exact log F_model(lam_* e^{i psi}) and exact (1/q) log F_main(lam_*^q e^{i q psi}):
    psis = np.linspace(-pi, pi, M_psi, endpoint=False)
    # Full log F_model (shifted to avoid overflow):
    def eval_logF(r_a, i_a, c_pow, L_val, angles, sign_par):
        lg = -c_pow * (log(L_val) + 1j * angles) + (1j * pi if sign_par % 2 == 1 else 0.0)
        rv = r_a + log(L_val)
        big = rv > 0.0
        sm = ~big
        if np.any(big):
            rb = rv[big]; ib = i_a[big]
            lg = lg + np.sum(rb) + 1j * (np.sum(ib) + len(rb)*pi) + len(rb)*1j*angles
            iu = np.exp(-rb - 1j*ib)
            epm = np.exp(-1j*angles)
            for k0 in range(0, len(iu), 250):
                ch = iu[k0:k0+250]
                lg += np.sum(np.log(1.0 - ch[None, :]*epm[:, None]), axis=1)
        if np.any(sm):
            rs = rv[sm]; is_ = i_a[sm]
            u = np.exp(rs + 1j*is_)
            epp = np.exp(1j*angles)
            for k0 in range(0, len(u), 250):
                ch = u[k0:k0+250]
                lg += np.sum(np.log(1.0 - ch[None, :]*epp[:, None]), axis=1)
        return lg

    logF_full = eval_logF(re_a, im_a, c, lam_star, psis, (j*(c+1)+c) % 2)
    # Exact main term: (1/q) log [ (-1)^? W^{-c} prod (1 - r_k^q W) ] at W = lam_*^q e^{i q psi}
    # Check sign of F_main at W = -Lambda_* (theta = pi, i.e., q*psi = pi):
    # Note: for r_0 = q, a_k^q = r_k^q > 0 for all k, and at psi = pi/q: w^q = -Lambda_*.
    # Let's extract the exact q-multiples part via discrete q-averaging over omega^q = 1!
    # Wait: (1/q) sum_{t=0}^{q-1} logF_full(psi + 2*pi*t/q) IS the exact q-multisection (up to 2pi i)!
    # Easier: compute directly from (1/q) eval_logF(q*re_a, 0, c, lam_*^q, q*psis, g_sign):
    logF_main_raw = (1.0 / q) * eval_logF(q * re_a, np.zeros_like(re_a), c, lam_star**q, q * psis, c % 2)
    # Wait: adjust constant imaginary offset by matching at psi = pi/q (where W = -Lambda_* is positive real!)
    diff = logF_full - logF_main_raw
    # Unwrap imaginary part modulo 2*pi around mean:
    d_re = diff.real
    d_im = np.angle(np.exp(1j * (diff.imag - np.mean(diff.imag))))
    # Check if mean(d_re) is ~0 (must be 0 by Parseval!)
    max_R_re = float(np.max(np.abs(d_re)))
    max_R_im = float(np.max(np.abs(d_im)))

    # Also compute exact rat_L1 of F_full on |w| = lam_*:
    sh = float(np.max(logF_full.real))
    F = np.exp(logF_full - sh)
    rat = float(np.mean(F).real) / float(np.mean(np.abs(F)))
    return q, g, lam_star, H_q, min_gap, max_R_re, max_R_im, rat

print("=== 1. EXACT q-MULTISECTION FOR D = 0 (r_0 = q >= 2) ===")
for c, j in [(4, 5), (4, 6), (5, 6), (6, 7), (6, 9), (8, 9), (8, 12), (10, 12), (12, 16), (20, 25), (30, 45), (50, 75)]:
    q = c // gcd(c, j); g = gcd(c, j)
    print(f"\nc={c:2d} j={j:2d} (q={q:2d}, g={g:2d}):")
    for cqk in [0.15, 0.30, 0.60, 1.00, 2.00, 5.00]:
        kap = cqk / (c * q)
        if kap >= 2.5 or c * kap * kap >= 42.0: continue
        _, _, lam, Hq, mg, mre, mim, rat = analyze_D0_multisection(c, j, kap)
        print(f"  cqk={cqk:4.2f} (kap={kap:.4f}): lam*={lam:.4f} H_q={Hq:6.3f} gap={mg:.3f} | max|Re R|={mre:.4f} max|Im R|={mim:.4f} | rat={rat:.5f}")
