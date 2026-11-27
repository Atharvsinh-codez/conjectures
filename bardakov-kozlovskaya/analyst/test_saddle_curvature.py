"""
TEST OUR 5-LINE GLOBAL-CURVATURE SADDLE ENGINE:
For ANY (c, j, kap) at lam = lam_L2 (and at explicit saddle lam_*):
Compute:
1. Exact first derivative: g1 = sum_k (a_k w_0)/(1 - a_k w_0) + c  (at w_0 = lam * e^{i psi_pk})
2. Exact second derivative (curvature): H = -d^2/dtheta^2 log|F(w_0 e^{i theta})|_{theta=0}
                                          = sum_k Re[ -a_k w_0 / (1 - a_k w_0)^2 ]
3. Exact third derivative (cubic phase): C3 = d^3/dtheta^3 arg F(w_0 e^{i theta})|_{theta=0}
4. Check the dimensionless saddle quality parameter:  Q_sad := |C3| / H^{3/2}  and  |g1| / sqrt(H)!!
   Remember: if H >> 1, |g1|/sqrt(H) << 1, and |C3|/H^{3/2} << 1, then a 5-line Laplace-method
   bound proves I(n,c,j) > 0 IMMEDIATELY!
"""
import numpy as np
from math import gcd

def saddle_diagnostics(c, j, kap, lam, psi0):
    alpha = kap / c
    M_cut = int(np.ceil(np.sqrt(60.0 / max(alpha, 1e-12)))) + 10
    ms = np.arange(-(c + M_cut), M_cut + 1, dtype=np.float64)
    re_a = -alpha * ms * (c + ms) + np.log(lam)
    im_a = -2.0 * np.pi * ms * j / c
    keep = re_a > -45.0
    ms = ms[keep]; re_a = re_a[keep]; im_a = im_a[keep]
    # z_m = a_m * w_0 = exp(re_a + i*(im_a + psi0))
    # For log(1 - z e^{i theta}):
    #   d/dtheta log(1 - z e^{i theta}) = -i z / (1 - z)
    #   d^2/dtheta^2 = z / (1 - z)^2
    #   d^3/dtheta^3 = i z (1 + z) / (1 - z)^3
    # Notice for large |z| >> 1:
    #   -z/(1-z) = 1 + 1/(z-1)
    #   z/(1-z)^2 = (1/z) / (1 - 1/z)^2
    #   z(1+z)/(1-z)^3 = -(1/z)(1+1/z)/(1-1/z)^3
    big = re_a > 0.0
    g1 = -c + 0.0j
    H = 0.0j
    C3 = 0.0j
    if np.any(big):
        inv_z = np.exp(-re_a[big] - 1j * (im_a[big] + psi0))
        g1 += np.sum(1.0 + inv_z / (1.0 - inv_z))
        H += np.sum(-inv_z / ((1.0 - inv_z) ** 2))
        C3 += np.sum(-inv_z * (1.0 + inv_z) / ((1.0 - inv_z) ** 3))
    sm = ~big
    if np.any(sm):
        z = np.exp(re_a[sm] + 1j * (im_a[sm] + psi0))
        g1 += np.sum(-z / (1.0 - z))
        H += np.sum(-z / ((1.0 - z) ** 2))
        C3 += np.sum(z * (1.0 + z) / ((1.0 - z) ** 3))
    # Note: d/dtheta log F = i * g1  => d(arg F)/dtheta = Re(g1)
    #       d^2/dtheta^2 log F = -H  => d^2(log|F|)/dtheta^2 = -Re(H)
    #       d^3/dtheta^3 log F = i*C3 => d^3(arg F)/dtheta^3 = Re(C3)
    return g1.real, H.real, C3.real, abs(g1.imag), abs(H.imag)

cases = [
    (500, 1000, 0.01, 1.32016, np.pi),
    (500, 1000, 0.003, 0.85608, np.pi),
    (100, 200, 0.10, 1.36081, np.pi),
    (100, 200, 0.01, 0.69931, np.pi),
    (200, 201, 0.03, 0.78634, np.pi),
    (200, 201, 0.01, 0.67049, np.pi),
    (200, 201, 0.002, 0.92715, np.pi),
    (200, 300, 0.01, 1.12020, np.pi/2),
    (200, 299, 0.01, 0.90707, np.pi/2),
    (200, 267, 0.01, 0.92380, np.pi),
    (200, 267, 0.01, 0.92380, np.pi/3),
    (200, 250, 0.01, 1.08280, np.pi/4),
    (50, 51, 0.03, 0.74436, np.pi),
    (2, 3, 0.05, 0.59684, np.pi/2),
    (2, 4, 0.01, 0.09469, np.pi),
]
print(f"{'c':>4} {'j':>5} {'kap':>6} {'lam':>7} {'psi0/pi':>7} | {'g1':>8} {'H':>9} {'C3':>10} | {'|g1|/sH':>8} {'|C3|/H^1.5':>10}")
for c, j, kap, lam, psi0 in cases:
    g1, H, C3, gi, Hi = saddle_diagnostics(c, j, kap, lam, psi0)
    sH = np.sqrt(max(H, 1e-12))
    print(f"{c:4d} {j:5d} {kap:6.3f} {lam:7.4f} {psi0/np.pi:7.3f} | {g1:8.3f} {H:9.2f} {C3:10.2f} | {abs(g1)/sH:8.4f} {abs(C3)/(sH**3):10.4f}")
