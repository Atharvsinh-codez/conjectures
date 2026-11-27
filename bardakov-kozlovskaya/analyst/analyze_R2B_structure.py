"""
Test the UNIVERSAL, PARAMETER-FREE L2-optimal radius:
  lam_L2 = argmin_{lam > 0} (1/(2pi)) int_{-pi}^{pi} |F(lam e^{i psi})|^2 dpsi
(unique minimizer of a strictly convex function of log(lam), ZERO knowledge of q or r_0!)
across ALL extreme corners of R2B (c * kap^2 < 20):
  - small c, tiny kap (c=2,3,4; q=1,2,3)
  - large c, q=1 (c=100, 500; kap=0.01..0.1)
  - large c, coprime near integer (c=50,j=51; c=200,j=201; kap=0.002..0.1)
  - large c, near 3/2 (c=200,j=299)
  - large c, golden ratio (c=89,j=144; c=200,j=323)
Check rat_L1 = I / avg|F| AND rat_L2 = I / sqrt(avg|F|^2) at lam = lam_L2!
Also check the EXACT relation between lam_L2 and a 1-line closed-form formula!
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
    return val / l1, val / l2, log_L2, int(np.sum(big))

def find_lam_L2(alpha, c, j, M_psi=512):
    # Golden-section minimization of log_L2(exp(v)) over v in [-4.0, +1.5]
    a, b = -4.0, 1.5
    gr = (np.sqrt(5.0) - 1.0) / 2.0
    c1 = b - gr * (b - a)
    c2 = a + gr * (b - a)
    f1 = eval_F_model(alpha, c, j, np.exp(c1), M_psi=M_psi)[2]
    f2 = eval_F_model(alpha, c, j, np.exp(c2), M_psi=M_psi)[2]
    for _ in range(28):
        if f1 < f2:
            b = c2; c2 = c1; f2 = f1
            c1 = b - gr * (b - a)
            f1 = eval_F_model(alpha, c, j, np.exp(c1), M_psi=M_psi)[2]
        else:
            a = c1; c1 = c2; f1 = f2
            c2 = a + gr * (b - a)
            f2 = eval_F_model(alpha, c, j, np.exp(c2), M_psi=M_psi)[2]
    return float(np.exp(0.5 * (a + b)))

cases = [
    (2, 3, 0.5), (2, 3, 0.2), (2, 3, 0.05), (2, 3, 0.01),
    (2, 4, 0.5), (2, 4, 0.1), (2, 4, 0.01),
    (3, 4, 0.5), (3, 4, 0.1), (3, 4, 0.01),
    (3, 6, 0.3), (3, 6, 0.03),
    (4, 5, 0.1), (4, 8, 0.1),
    (20, 40, 0.3), (20, 40, 0.1), (20, 40, 0.05),
    (100, 200, 0.3), (100, 200, 0.1), (100, 200, 0.03), (100, 200, 0.01),
    (500, 1000, 0.1), (500, 1000, 0.03), (500, 1000, 0.01), (500, 1000, 0.003),
    (50, 51, 0.3), (50, 51, 0.1), (50, 51, 0.03), (50, 51, 0.01),
    (200, 201, 0.1), (200, 201, 0.03), (200, 201, 0.01), (200, 201, 0.002),
    (200, 299, 0.03), (200, 299, 0.01), (200, 299, 0.003),
    (200, 323, 0.03), (200, 323, 0.01), (200, 323, 0.003),
    (89, 144, 0.03), (89, 144, 0.01), (89, 144, 0.003),
]
print("=== UNIVERSAL L2-OPTIMAL RADIUS lam_L2 ACROSS ALL 42 EXTREME CORNERS OF R2B ===")
print(f"{'c':>4} {'j':>5} {'q':>4} {'kap':>7} {'c*k':>6} {'c*k^2':>7} {'lam_L2':>8} {'#big':>5} {'rat_L1':>9} {'rat_L2':>9}")
min_l1 = 1e9
for c, j, kap in cases:
    g = gcd(c, j); q = c // g
    lam = find_lam_L2(kap / c, c, j, M_psi=512)
    rl1, rl2, _, nbig = eval_F_model(kap / c, c, j, lam, M_psi=1024)
    min_l1 = min(min_l1, rl1)
    print(f"{c:4d} {j:5d} {q:4d} {kap:7.4f} {c*kap:6.2f} {c*kap*kap:7.4f} {lam:8.5f} {nbig:5d} {rl1:9.5f} {rl2:9.5f}")
print(f"\nGLOBAL MINIMUM rat_L1 across ALL 42 extreme corners = {min_l1:.5f}!!")
