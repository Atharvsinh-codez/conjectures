"""
Diagnose the exact structure of the gap G = {c*kap^2 < 16 AND c*q*kap > 1.0}:
1. Test c = 200 across a sweep of j/c fractions (201/200, 202/200, 205/200, 220/200, 241/200,
   250/200, 267/200, 299/200, 300/200, 323/200, 400/200) at kap = 0.01 and 0.03!
   Print: q, max|Phi@1|, rat@1, lam_L2, rat@lam_L2, #peaks, peak psi, log|F| range!
2. Also check the EXACT tightened Tool-1 threshold C_1 (using exact |1-(1+zeta_m)^{-n}| and
   min with individual log bound) across n = 80, 160, 320, 1000!
"""
import numpy as np
from math import gcd

def eval_F_model(alpha, c, j, lam, M_psi=1024):
    M_cut = int(np.ceil(np.sqrt(60.0 / max(alpha, 1e-12)))) + 10
    ms = np.arange(-(c + M_cut), M_cut + 1, dtype=np.float64)
    re_a = -alpha * ms * (c + ms) + np.log(lam)
    im_a = -2.0 * np.pi * ms * j / c
    keep = re_a > -45.0
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
    l1 = float(np.mean(np.abs(F)))
    l2 = float(np.sqrt(np.mean(np.abs(F)**2)))
    log_L2 = sh + np.log(l2)
    # Count local peaks within 1.0 of max log|F|
    reF = log_F.real - sh
    # Find dominant peak angle
    idx_pk = int(np.argmax(reF))
    # Max |arg(F)| where |F| > 0.01 * max|F| (the active peak region!)
    active = reF > -4.6
    max_arg_active = float(np.max(np.abs(np.angle(F[active]))))
    frac_active = float(np.mean(active))
    return val / l1, log_L2, psis[idx_pk], max_arg_active, frac_active, int(np.sum(big))

def find_lam_L2(alpha, c, j, M_psi=512):
    a, b = -2.0, 1.0
    gr = (np.sqrt(5.0) - 1.0) / 2.0
    c1 = b - gr * (b - a); c2 = a + gr * (b - a)
    f1 = eval_F_model(alpha, c, j, np.exp(c1), M_psi=M_psi)[1]
    f2 = eval_F_model(alpha, c, j, np.exp(c2), M_psi=M_psi)[1]
    for _ in range(25):
        if f1 < f2:
            b = c2; c2 = c1; f2 = f1; c1 = b - gr * (b - a)
            f1 = eval_F_model(alpha, c, j, np.exp(c1), M_psi=M_psi)[1]
        else:
            a = c1; c1 = c2; f1 = f2; c2 = a + gr * (b - a)
            f2 = eval_F_model(alpha, c, j, np.exp(c2), M_psi=M_psi)[1]
    return float(np.exp(0.5 * (a + b)))

print("=== SWEEP OF j/c AT c = 200, kap = 0.01 (c*k = 2.0, c*k^2 = 0.02) ===")
print(f"{'j':>4} {'j/c':>7} {'q':>4} | {'rat@1':>8} {'argAct@1':>8} | {'lam_L2':>7} {'#big':>4} {'rat@L2':>8} {'argAct@L2':>9} {'fracAct':>7} {'psi_pk':>7}")
c = 200; kap = 0.01; alpha = kap / c
for j in [201, 202, 205, 220, 241, 250, 267, 299, 300, 323, 400]:
    q = c // gcd(c, j)
    r1, _, _, aa1, _, _ = eval_F_model(alpha, c, j, 1.0, M_psi=1024)
    lam = find_lam_L2(alpha, c, j, M_psi=512)
    r2, _, ppk, aa2, fa2, nb = eval_F_model(alpha, c, j, lam, M_psi=1024)
    print(f"{j:4d} {j/c:7.3f} {q:4d} | {r1:8.5f} {aa1:8.4f} | {lam:7.4f} {nb:4d} {r2:8.5f} {aa2:9.4f} {fa2:7.3f} {ppk:7.3f}")

print("\n=== SWEEP OF j/c AT c = 200, kap = 0.03 (c*k = 6.0, c*k^2 = 0.18) ===")
kap = 0.03; alpha = kap / c
for j in [201, 205, 220, 241, 250, 267, 299, 300, 323, 400]:
    q = c // gcd(c, j)
    r1, _, _, aa1, _, _ = eval_F_model(alpha, c, j, 1.0, M_psi=1024)
    lam = find_lam_L2(alpha, c, j, M_psi=512)
    r2, _, ppk, aa2, fa2, nb = eval_F_model(alpha, c, j, lam, M_psi=1024)
    print(f"{j:4d} {j/c:7.3f} {q:4d} | {r1:8.5f} {aa1:8.4f} | {lam:7.4f} {nb:4d} {r2:8.5f} {aa2:9.4f} {fa2:7.3f} {ppk:7.3f}")
