"""
CRITICAL DIAGNOSTIC: Test golden-ratio j/c (e.g. c=200, j=323; c=200, j=299; c=89, j=144)
at kap = 1.0, 0.3, 0.1, 0.03, 0.01 across a grid of radii lam!
Does a circle |w| = lam STILL achieve rat >= 0.5 even for irrational-like j/c?!
Also check the EXACT sign of (-1)^(j+c) [w^c] prod (1 - a_m w) (sanity check!).
"""
import numpy as np

def eval_model_at_lam(alpha, c, j, lam, M_psi=1024):
    M_cut = int(np.ceil(np.sqrt(60.0 / max(alpha, 1e-12)))) + 10
    ms = np.arange(-(c + M_cut), M_cut + 1, dtype=np.float64)
    re_a = -alpha * ms * (c + ms) + np.log(lam)
    im_a = -2.0 * np.pi * ms * j / c
    keep = re_a > -45.0
    ms = ms[keep]; re_a = re_a[keep]; im_a = im_a[keep]
    big = re_a > 0.0

    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    log_F = -c * (np.log(lam) + 1j * psis)
    if (j + c) % 2 == 1:
        log_F = log_F + 1j * np.pi

    if np.any(big):
        rb = re_a[big]; ib = im_a[big]
        log_F = log_F + np.sum(rb) + 1j * (np.sum(ib) + len(rb) * np.pi) + len(rb) * 1j * psis
        inv_u = np.exp(-rb - 1j * ib)
        ep_neg = np.exp(-1j * psis)
        for k0 in range(0, len(inv_u), 200):
            ch = inv_u[k0:k0+200]
            log_F += np.sum(np.log(1.0 - ch[None, :] * ep_neg[:, None]), axis=1)

    sm = ~big
    if np.any(sm):
        rs = re_a[sm]; is_ = im_a[sm]
        u = np.exp(rs + 1j * is_)
        ep_pos = np.exp(1j * psis)
        for k0 in range(0, len(u), 200):
            ch = u[k0:k0+200]
            near1 = np.abs(ch) > 1.0 - 1e-9
            if np.any(near1):
                prod_n1 = np.prod(1.0 - ch[near1][None, :] * ep_pos[:, None], axis=1)
                with np.errstate(divide='ignore'):
                    log_F += np.where(np.abs(prod_n1) > 0, np.log(prod_n1), -1e9)
            far = ~near1
            if np.any(far):
                log_F += np.sum(np.log(1.0 - ch[far][None, :] * ep_pos[:, None]), axis=1)

    sh = float(np.max(log_F.real))
    F = np.exp(log_F - sh)
    val = float(np.mean(F).real)
    rat = val / float(np.mean(np.abs(F)))
    idx_pk = int(np.argmax(log_F.real))
    return rat, psis[idx_pk], float(np.angle(F[idx_pk]))

lams = [0.5, 0.7, 0.85, 0.92, 0.96, 1.0, 1.04, 1.08, 1.15, 1.3, 1.6]
for c, j, label in [(89, 144, "Fibonacci 144/89 ~ phi"), (200, 299, "299/200 ~ 3/2"), (200, 323, "323/200 ~ phi"), (55, 89, "Fibonacci 89/55 ~ phi")]:
    print(f"\n=== {label}: c={c}, j={j} ===")
    for kap in [1.0, 0.3, 0.1, 0.03, 0.01, 0.003]:
        best_rat = -1e9; best_info = None
        rat_at_1 = None
        for lam in lams:
            rat, ppk, phpk = eval_model_at_lam(kap/c, c, j, lam, M_psi=1024)
            if lam == 1.0: rat_at_1 = rat
            if rat > best_rat:
                best_rat = rat; best_info = (lam, ppk, phpk)
        print(f"  kap={kap:6.3f} (c*k={c*kap:6.2f}): rat@1={rat_at_1:8.5f} | BEST rat={best_rat:8.5f} at lam={best_info[0]:.2f} (psi_pk={best_info[1]:6.3f}, arg={best_info[2]:6.3f})")
