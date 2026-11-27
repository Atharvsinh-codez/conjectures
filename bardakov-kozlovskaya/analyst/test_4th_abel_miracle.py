"""
TEST THE 1-LINE 4TH-ORDER PAIRED COSINE ABEL MIRACLE IDENTITY!!
  16 * sin^4(th/2) * (B_r - A_r)
    = 2 * [ (F_i(2) - F_e(2)) - 4*(F_i(1) - F_e(1)) ]
      + 2 * (d4Fi(1) - d4Fe(1)) * cos(th)
      + 2 * sum_{m=2}^{inf} d4Fi(m) cos(m*th)
      - sum_{m=2}^{c-2} d4Fe(m) cos(m*th) !!
Verify:
1. Machine-precision (1e-14) identity check!
2. Tightness of the triangle upper bound ub4 vs exact |B_r - A_r| / r!
3. Closed-form scaling constant C4 in |B_r - A_r| / r <= (C4 / sin^4(pi*r*nu)) * (r * kap^2 / c) [or exact formula]!
4. Re-run the contrapositive map of B_minor >= 1.45 across c = 2..100 with BOTH tight b_r^pair AND b_r^{Abel,4}!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd

def test_4th_abel_miracle(c, jm, kap, r):
    nu = jm / c
    th = 2.0 * pi * r * nu
    s2 = 4.0 * (sin(th / 2.0) ** 2)
    if s2 < 1e-12: return None
    alpha = kap / c
    M = int(np.ceil(sqrt(45.0 / max(r * alpha, 1e-10)))) + 15
    ms_e = np.arange(0, c + 1, dtype=np.float64)
    Fe = np.exp(-r * alpha * ms_e * (c - ms_e))
    ms_i = np.arange(0, M + 1, dtype=np.float64)
    Fi = np.exp(-r * alpha * ms_i * (c + ms_i))
    exact = (2.0 * float(np.sum(Fi[1:] * np.cos(ms_i[1:] * th)))
             - float(np.sum(Fe[1:c] * np.cos(ms_e[1:c] * th))))

    # h(m) = F(m+1) - 2 F(m) + F(m-1) for m = 1..N-1, with h(0) = h(N) = 2*(F(1) - 1)
    he = np.concatenate([[2*(Fe[1] - 1.0)], Fe[2:] - 2*Fe[1:-1] + Fe[:-2], [2*(Fe[1] - 1.0)]])
    hi = np.concatenate([[2*(Fi[1] - 1.0)], Fi[2:] - 2*Fi[1:-1] + Fi[:-2], [0.0]])
    d4Fe = he[2:] - 2*he[1:-1] + he[:-2] # indices m = 1..c-1 (length c-1)
    d4Fi = hi[2:] - 2*hi[1:-1] + hi[:-2] # indices m = 1..M-1 (length M-1)

    bdy = 2.0 * ((Fi[2] - Fe[2]) - 4.0 * (Fi[1] - Fe[1]))
    # Note: d4Fe[0] == d4Fe[-1] (m=1 and m=c-1), and cos((c-1)*th) == cos(th)!
    # So pair m=1: 2*d4Fi[0] - (d4Fe[0] + d4Fe[-1]) = 2*(d4Fi[0] - d4Fe[0])!
    # Wait: if c == 2, len(d4Fe) == 1 (only m=1), wait: for c >= 3, m=1 and m=c-1 are distinct!
    rhs = (bdy + 2.0 * float(np.sum(d4Fi * np.cos(ms_i[1:-1] * th)))
           - float(np.sum(d4Fe * np.cos(ms_e[1:-1] * th)))) / (s2 * s2)
    err = abs(exact - rhs)

    # Triangle upper bound (pairing m=1!):
    if c >= 3:
        m1_term = 2.0 * abs(d4Fi[0] - d4Fe[0])
        rest_i = 2.0 * float(np.sum(np.abs(d4Fi[1:])))
        rest_e = float(np.sum(np.abs(d4Fe[1:-1])))
        # Even better: pair d4Fi[m] - d4Fe[m] for m = 1..floor((c-1)/2) or keep separate?
        # Check separate vs paired!
        ub4_sep = (abs(bdy) + m1_term + rest_i + rest_e) / (s2 * s2 * r)
        # Also paired up to mh = (c-1)//2:
        mh = (c - 1) // 2
        pair_d4 = 2.0 * float(np.sum(np.abs(d4Fi[:mh] - d4Fe[:mh])))
        mid_d4 = abs(d4Fe[mh]) if (c % 2 == 0) else 0.0
        tail_d4 = 2.0 * float(np.sum(np.abs(d4Fi[mh:])))
        ub4_pair = (abs(bdy) + pair_d4 + mid_d4 + tail_d4) / (s2 * s2 * r)
    else:
        ub4_sep = ub4_pair = (abs(bdy) + 2*float(np.sum(np.abs(d4Fi))) + float(np.sum(np.abs(d4Fe)))) / (s2 * s2 * r)

    # Compare num = (abs(bdy) + m1_term + rest_i + rest_e) / (16 * r) against scale!
    num_sep = (abs(bdy) + m1_term + rest_i + rest_e) / (16.0 * r) if c >= 3 else 0.0
    num_pair = (abs(bdy) + pair_d4 + mid_d4 + tail_d4) / (16.0 * r) if c >= 3 else 0.0
    return abs(exact) / r, ub4_sep, ub4_pair, err, num_sep, num_pair

print("=== 1. VERIFY 1-LINE 4TH-ORDER PAIRED ABEL MIRACLE ===")
for c, jm, kap, r in [
    (100, 1, 0.010, 1), (100, 2, 0.010, 1), (100, 3, 0.010, 1), (100, 5, 0.010, 1),
    (100, 3, 0.005, 1), (100, 3, 0.020, 1), (100, 9, 0.002, 11), (100, 9, 0.002, 1),
    (50, 1, 0.020, 1), (50, 2, 0.020, 1), (50, 5, 0.020, 1),
    (30, 1, 0.050, 1), (30, 2, 0.050, 1), (30, 4, 0.050, 1),
    (20, 1, 0.100, 1), (10, 1, 0.150, 1), (6, 1, 0.250, 1), (200, 123, 0.005, 1),
    (100, 1, 0.050, 1), (100, 1, 0.100, 1), (100, 1, 0.200, 1),
]:
    ex, u_sep, u_pair, err, ns, np_ = test_4th_abel_miracle(c, jm, kap, r)
    tau = r * c * kap
    # Theoretical scaling of num_sep / (r * kap^2 / c):
    # Wait: check BOTH r*kap^2/c AND (r*kap^2/c)*... for small/large tau!
    sc1 = r * (kap**2) / c
    # At large tau: int_0^inf |d^4/dm^4 e^{-r*kap*m}| dm ~ (r*kap)^3!
    #Wait: d4Fi(m) ~ (r*kap)^4 e^{-r*kap*m} => sum |d4Fi| ~ (r*kap)^3! Divided by r -> r^2 * kap^3!
    #Wait: THAT'S ONLY FOR UNPAIRED d4Fi! For PAIRED d4Fi - d4Fe:
    #  d4(e^{-r*kap*m} * 2*r*alpha*m^2) ~ (r*kap)^2 * (r*alpha) * ... => sum ~ r^2 * alpha * kap = r^2 * kap^2 / c!!
    #  Divided by r -> r * kap^2 / c!!
    print(f"  c={c:3d} jm={jm:3d} kap={kap:.3f} r={r:2d} (tau={tau:4.1f}): ex={ex:.5f} | u4_sep={u_sep:.5f} u4_pair={u_pair:.5f} | err={err:.1e} | num_p/(r*k^2/c)={np_/sc1:.3f}")
