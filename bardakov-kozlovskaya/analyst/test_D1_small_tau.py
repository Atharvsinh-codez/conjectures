"""
CRITICAL STRESS TEST: c = 2000, j = 2001, tau = 0.40 (kap = 0.0002, D = 1)
and c = 5000, j = 5001, tau = 0.40 (kap = 0.00008, D = 1)!
Here c * g1@1 = 2000 * 0.002033 = 4.07 and 5000 * 0.002033 = 10.17 (so lam = 1 fails!),
AND lam_* in the continuum sits right at e^{-tau/4} = 0.9048374!
What happens to rat(lam), argAct(lam), H(lam), and psi_pk(lam) as lam varies around 0.90484..1.00?!
"""
import numpy as np

def eval_F_fast(alpha, c, j, lam, M_psi=512):
    M_cut = int(np.ceil(np.sqrt(40.0 / max(alpha, 1e-12)))) + 10
    ms = np.arange(-(c + M_cut), M_cut + 1, dtype=np.float64)
    re_a = -alpha * ms * (c + ms) + np.log(lam)
    im_a = -2.0 * np.pi * ms * j / c
    keep = re_a > -35.0
    ms = ms[keep]; re_a = re_a[keep]; im_a = im_a[keep]
    big = re_a > 0.0

    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    log_F = -c * (np.log(lam) + 1j * psis)
    if (j * (c + 1) + c) % 2 == 1:
        log_F = log_F + 1j * np.pi

    if np.any(big):
        rb = re_a[big]; ib = im_a[big]
        log_F = log_F + np.sum(rb) + 1j * (np.sum(ib) + len(rb) * np.pi) + len(rb) * 1j * psis
        inv_u = np.exp(-rb - 1j * ib)
        ep_neg = np.exp(-1j * psis)
        for k0 in range(0, len(inv_u), 400):
            ch = inv_u[k0:k0+400]
            log_F += np.sum(np.log(1.0 - ch[None, :] * ep_neg[:, None]), axis=1)

    sm = ~big
    if np.any(sm):
        rs = re_a[sm]; is_ = im_a[sm]
        u = np.exp(rs + 1j * is_)
        ep_pos = np.exp(1j * psis)
        for k0 in range(0, len(u), 400):
            ch = u[k0:k0+400]
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
    l1 = float(np.mean(np.abs(F)))
    idx = int(np.argmax(log_F.real))
    active = (log_F.real - sh) > -4.6
    aa = float(np.max(np.abs(np.angle(F[active]))))
    return val / l1, psis[idx], aa, float(np.max(log_F.real) - np.min(log_F.real)), int(np.sum(big))

for c in [2000, 5000]:
    j = c + 1; tau = 0.40; kap = tau / c; alpha = kap / c
    em = np.exp(-tau / 4.0)
    print(f"\n=== c={c}, j={j}, tau={tau}, e^(-tau/4)={em:.6f} ===")
    for lam in [em*0.99, em*1.0005, em*1.002, em*1.005, em*1.01, em*1.02, 0.94, 0.96, 0.98, 1.00]:
        rat, ppk, aa, rng, nb = eval_F_fast(alpha, c, j, lam, M_psi=512)
        print(f"  lam={lam:8.6f} (#big={nb:4d}): rat={rat:9.5f} | psi_pk={ppk:7.3f} | argAct={aa:7.4f} | rng={rng:6.1f}")
