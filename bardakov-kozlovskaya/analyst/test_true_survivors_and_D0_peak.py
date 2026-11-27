"""
1. CHECK TRUE PEAK-ONLY |Re R|, |Im R| AND GLOBAL UPPER BOUND max(Re R) FOR D = 0 (q >= 2)
   WHEN c*q*kap > 0.85! (No longer fooled by negative -inf root-hit dips in the tail!)
2. CHECK TRUE SURVIVING SET OF B_minor_3way >= 1.45 CLASSIFIED BY THE *ACTUAL* DOMINANT
   HARMONIC r_0 = argmax_r b_r (NOT truncated at sqrt(c))!
   Compare BOTH:
     (a) Triangle pair bound C_p2 = 0.60 everywhere
     (b) Exact signed pair bound for q | r (where cos(2*pi*r*m*j/c) = +1!)
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd

# =============================================================================
# 1. TRUE PEAK BEHAVIOR & GLOBAL UPPER BOUND OF R_{non-q} FOR D = 0 (q >= 2)
# =============================================================================
print("=== 1. TRUE PEAK & GLOBAL UPPER BOUND OF R_{non-q} FOR D = 0 (q >= 2, c*q*kap > 0.85) ===")
for c, j, cqk in [
    (6, 9, 1.2), (8, 12, 1.5), (12, 18, 1.5), (12, 18, 3.0), (12, 18, 6.0),
    (15, 20, 1.5), (15, 20, 3.0), (20, 30, 1.5), (20, 30, 3.0), (20, 25, 1.5),
    (20, 25, 3.0), (20, 25, 6.0), (30, 45, 1.5), (30, 45, 3.0), (30, 45, 6.0),
    (30, 40, 1.5), (30, 40, 3.0), (40, 50, 1.5), (40, 50, 3.0), (50, 75, 2.0),
    (50, 75, 5.0), (50, 60, 1.5), (50, 60, 3.0),
]:
    g = gcd(c, j); q = c // g; kap = cqk / (c * q)
    alpha = kap / c
    M_cut = int(np.ceil(sqrt(45.0 / max(alpha, 1e-12)))) + 10
    ms = np.arange(-(c + M_cut), M_cut + 1, dtype=np.float64)
    re_a = -alpha * ms * (c + ms)
    im_a = -2.0 * pi * ms * j / c
    keep = re_a > -40.0
    re_a = re_a[keep]; im_a = im_a[keep]

    # Exact positive-root saddle for W = -lam^q:
    lo, hi = -80.0, 40.0
    for _ in range(60):
        vm = 0.5 * (lo + hi)
        z = np.clip(q * re_a + vm, -70.0, 70.0)
        if float(np.sum(1.0 / (1.0 + np.exp(-z)))) > c: hi = vm
        else: lo = vm
    V = 0.5 * (lo + hi)
    lam = exp(V / q)
    rv = re_a + log(lam)
    if float(np.min(np.abs(rv))) < 1e-4:
        lam *= exp(0.002 / q)
        V = q * log(lam); rv = re_a + log(lam)

    z = np.clip(q * re_a + V, -70.0, 70.0)
    p = 1.0 / (1.0 + np.exp(-z))
    H_q = float(np.sum(p * (1.0 - p))) / q

    # Evaluate EXACT ratio F_full(lam e^{i psi}) / F_main(lam^q e^{i q psi})^{1/q} directly on ONE sector psi in [0, 2*pi/q)!
    # On this single sector, peak is at psi_0 = pi / q (theta = q*psi - pi in [-pi, pi])!
    # Wait: F_full(lam e^{i psi}) = (-1)^{j(c+1)+c} (lam e^{i psi})^{-c} prod_k (1 - a_k lam e^{i psi})
    #       F_main_q(lam e^{i psi}) := (lam^{-c} e^{-i c psi}) * [ prod_k (1 - |a_k|^q lam^q e^{i q psi}) ]^{1/q}
    # To avoid ANY branch cut in log, compute the exact complex ratio:
    #   Q(theta) := F_full(lam e^{i (pi + theta)/q}) / F_main_pos(theta)
    # where F_main_pos(theta) > 0 at theta = 0 (at theta = 0, w^q = -lam^q, so 1 - |a_k|^q w^q = 1 + |a_k|^q lam^q > 1!).
    thetas = np.linspace(-pi, pi, 512, endpoint=False)
    psis = (pi + thetas) / q
    # Compute log |F_full| and arg(F_full) safely (for big and small roots):
    big = rv > 0.0; sm = ~big
    def eval_F_log(r_arr, i_arr, c_p, angs, sgn_p):
        lg = -c_p * (log(lam) + 1j * angs) + (1j * pi if sgn_p % 2 == 1 else 0.0)
        rv_loc = r_arr + log(lam)
        b_loc = rv_loc > 0.0; s_loc = ~b_loc
        if np.any(b_loc):
            rb = rv_loc[b_loc]; ib = i_arr[b_loc]
            lg = lg + np.sum(rb) + 1j*(np.sum(ib) + len(rb)*pi) + len(rb)*1j*angs
            ch = np.exp(-rb - 1j*ib)
            ep = np.exp(-1j*angs)
            for k0 in range(0, len(ch), 200):
                lg += np.sum(np.log(1.0 - ch[k0:k0+200][None, :] * ep[:, None]), axis=1)
        if np.any(s_loc):
            rs = rv_loc[s_loc]; is_ = i_arr[s_loc]
            ch = np.exp(rs + 1j*is_)
            ep = np.exp(1j*angs)
            for k0 in range(0, len(ch), 200):
                lg += np.sum(np.log(1.0 - ch[k0:k0+200][None, :] * ep[:, None]), axis=1)
        return lg

    # F_full on sector psi = (pi + theta)/q:
    lg_full = eval_F_log(re_a, im_a, c, psis, (j*(c+1)+c) % 2)
    # Exact positive-root F_main(theta) = F_main(0) * exp(-i g theta) * prod_k ((1 + x_k e^{i theta})/(1 + x_k))^{1/q}
    # with x_k = |a_k|^q lam^q = exp(q * rv):
    # Note at theta = 0, F_main(0) = lam^{-c} * prod_k (1 + exp(q*rv))^{1/q} > 0!
    log_F0_main = -c * log(lam) + (1.0 / q) * float(np.sum(np.where(q*rv > 30, q*rv, np.log1p(np.exp(np.clip(q*rv, -80, 30))))))
    # Exact log(F_main(theta) / F_main(0)):
    lg_main_rel = -1j * (c / q) * thetas
    # Split q*rv > 0 vs <= 0 for exact numerical stability:
    if np.any(big):
        ix = np.exp(-q * rv[big])
        epm = np.exp(-1j * thetas)
        lg_main_rel = lg_main_rel + (len(ix) / q) * 1j * thetas
        for k0 in range(0, len(ix), 200):
            ch = ix[k0:k0+200]
            lg_main_rel += (1.0 / q) * np.sum(np.log((1.0 + ch[None, :]*epm[:, None]) / (1.0 + ch[None, :])), axis=1)
    if np.any(sm):
        x = np.exp(q * rv[sm])
        epp = np.exp(1j * thetas)
        for k0 in range(0, len(x), 200):
            ch = x[k0:k0+200]
            lg_main_rel += (1.0 / q) * np.sum(np.log((1.0 + ch[None, :]*epp[:, None]) / (1.0 + ch[None, :])), axis=1)
    lg_main = log_F0_main + lg_main_rel

    # Ratio Q(theta) = F_full / F_main = exp(lg_full - lg_main)!
    d_re = lg_full.real - lg_main.real
    d_im = np.angle(np.exp(1j * (lg_full.imag - lg_main.imag)))
    pk = np.abs(thetas) <= pi / 2
    idx0 = int(np.argmin(np.abs(thetas)))
    print(f"  c={c:2d} j={j:2d} (q={q},g={g:2d}) cqk={cqk:.1f} Hq={H_q:5.2f}: @0=({d_re[idx0]:+.4f},{d_im[idx0]:+.4f}) | PK: |Re|<={np.max(np.abs(d_re[pk])):.4f} |Im|<={np.max(np.abs(d_im[pk])):.4f} | ALL: maxRe={np.max(d_re):+.4f} minRe={np.min(d_re):+.4f}")

# =============================================================================
# 2. TRUE SURVIVING SET OF B_minor_3way >= 1.45 (CLASSIFIED BY DOMINANT r0!)
# =============================================================================
print("\n=== 2. TRUE SURVIVING SET OF B_minor_3way >= 1.45 (USING EXACT SIGNED q|r BOUND!) ===")
def B_minor_3way_detail(c, jm, kap, use_signed_q=True):
    nu = jm / c
    g = gcd(c, jm) if jm > 0 else c
    q = c // g
    tot = 0.0
    max_b = 0.0; dom_r = 1
    for r in range(1, 101):
        tau_r = r * c * kap
        if use_signed_q and (r % q == 0):
            # Exact signed bound for q|r (tau_r >= 0.85): |B_r - A_r| / c <= min(0.36, 9.5/tau_r^2)
            # Wait: check tau_r < 0.85? In G*, c*q*kap > 0.85, so for ALL r with q|r, tau_r = r*c*kap >= c*q*kap > 0.85!
            bp = (c / r) * min(0.36, 9.5 / (tau_r * tau_r))
        else:
            bp = (c / r) * min(4.8 / (tau_r * tau_r), 0.60 / sqrt(tau_r))
        s2 = sin(pi * r * nu) ** 2
        if s2 > 1e-12:
            b2 = (1.35 / s2) * min(1.0 / (r * c), sqrt(kap / (r * c)))
            b4 = (4.5 / (s2 * s2)) * (r * kap * kap / c)
            br = min(bp, b2, b4)
        else:
            br = bp
        tot += br
        if br > max_b:
            max_b = br; dom_r = r
    tot += 4.8 / (2.0 * 10000.0 * c * kap * kap)
    return tot, dom_r, max_b

kaps = np.geomspace(0.002, 2.49, 60)
for c in [2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50, 100]:
    surv_D0 = []
    surv_Dne0 = []
    for jm in range(0, c // 2 + 1):
        g = gcd(c, jm) if jm > 0 else c
        q = c // g
        for kap in kaps:
            if c * q * kap <= 0.85 or c * kap * kap >= 42.0: continue
            bm, r0, mb = B_minor_3way_detail(c, jm, kap, use_signed_q=True)
            if bm >= 1.45:
                p0 = round(r0 * jm / c)
                D = abs(r0 * jm - c * p0)
                if D == 0:
                    surv_D0.append((jm, q, g, kap, bm, r0))
                else:
                    surv_Dne0.append((jm, q, r0, D, c / r0, kap, c * r0 * kap, bm, mb))
    d0_qs = sorted(set((x[1], x[2]) for x in surv_D0))
    dne0_min = min((x[4] / (x[3]**2) for x in surv_Dne0), default=999.0)
    dne0_s = sorted(set((x[2], x[3], round(x[4], 1)) for x in surv_Dne0))[:6]
    print(f"  c={c:3d}: D=0 surv={len(surv_D0):3d} (q,g)={d0_qs[:6]} | D!=0 surv={len(surv_Dne0):3d} min(ct/D^2)={dne0_min:6.1f} (r0,D,ct)={dne0_s}")
