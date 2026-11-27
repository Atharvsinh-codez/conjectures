"""
Inspect the exact behavior of F(lam * e^{i psi}) for the tricky coprime-large-c case
(c=200, j=201, kap=0.01 and c=50, j=51, kap=0.1) AND the q=1 case (c=500, j=1000, kap=0.01)!
Fix the (-1)^c sign factor so (-1)^(j+c) gives the exact positive integral!
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
    # Also find where |F| peaks and what log|F| range is:
    idx_pk = int(np.argmax(log_F.real))
    psi_pk = psis[idx_pk]
    phi_pk = np.angle(F[idx_pk])
    # Exact saddle equation sum_m (-a_m w)/(1 - a_m w) at w = lam*e^{i psi_pk}:
    w_pk = lam * np.exp(1j * psi_pk)
    # Using z = a_m w: -z/(1-z) = 1 - 1/(1-z)
    # Wait: compute log_F'(w) = d/d(log w) log F(w) numerically via finite diff in psi!
    d_re = (log_F[(idx_pk+1)%M_psi].real - log_F[(idx_pk-1)%M_psi].real) / (2 * 2*np.pi/M_psi)
    d_im = (np.angle(F[(idx_pk+1)%M_psi] / F[(idx_pk-1)%M_psi])) / (2 * 2*np.pi/M_psi)
    return rat, psi_pk, phi_pk, d_im, float(np.max(log_F.real) - np.min(log_F.real))

for c, j, kap in [(200, 201, 0.01), (50, 51, 0.1), (50, 51, 0.01), (200, 201, 0.002), (500, 1000, 0.01)]:
    print(f"\n=== c={c}, j={j}, kap={kap} (c*kap={c*kap:.2f}, alpha={kap/c:.2e}) ===")
    for lam in [0.4, 0.6, 0.75, 0.85, 0.92, 0.96, 1.0, 1.04, 1.08, 1.15, 1.3, 1.5]:
        rat, ppk, phpk, dim, rng = eval_model_at_lam(kap/c, c, j, lam, M_psi=1024)
        print(f"  lam={lam:5.2f}: rat={rat:9.5f} | psi_pk={ppk:7.3f} | arg(F_pk)={phpk:7.3f} | d(arg)/dpsi@pk={dim:8.3f} | log|F|_rng={rng:6.1f}")
