"""
TEST THE TAUTOLOGICAL DICHOTOMY X := sum_{r>=1} b_r(c, j, kap) < 1.50  vs  >= 1.50!
where:
  b_r := min( pair_bound(r), abel_1step_bound(r) )
with:
  pair_bound(r) ~ (4 / (c * r^3 * kap^2))  [exact: sum_{m} |a_m^r - a_{-m}^{-r}|/r]
  abel_1step_bound(r) = (K_abel / sin^2(pi*r*nu)) * min(1/(r*c), sqrt(kap/(r*c)))  (K_abel = 1.5)
Check X on:
  - Golden ratio / Fibonacci points (c=89,j=144; c=200,j=323; c=500,j=809) at kap = 0.001..0.5!
  - Small c (c=2..15) with c*q*kap > 0.18!
  - Major arcs (c=200,j=201; c=200,j=299; c=200,j=267; c=50,j=51)!
"""
import numpy as np
from math import pi, sin, sqrt, gcd, exp, log

def compute_X(c, j, kap, K_abel=2.0, R_max=400):
    alpha = kap / c
    nu = j / c
    # Exact model pair_bound(r) = (2/r) * sum_{m=1}^{(c-1)//2} |e^{-r*alpha*m*(c+m)} - e^{-r*alpha*m*(c-m)}|
    #                            + (1_{c even}/r) * e^{-r*alpha*c^2/4}
    #                            + (2/r) * sum_{m > (c-1)//2} e^{-r*alpha*m*(c+m)}
    mh = (c - 1) // 2
    ms_p = np.arange(1, mh + 1, dtype=np.float64)
    M_cut = int(np.ceil(np.sqrt(40.0 / max(alpha, 1e-12)))) + 5
    ms_t = np.arange(mh + 1, mh + M_cut + 1, dtype=np.float64)
    tot = 0.0
    max_term = (0.0, 0)
    for r in range(1, R_max + 1):
        # Exact pair bound at r:
        ep_in = np.exp(-r * alpha * ms_p * (c + ms_p))
        ep_ex = np.exp(-r * alpha * ms_p * (c - ms_p))
        pb = (2.0 / r) * float(np.sum(np.abs(ep_ex - ep_in)))
        if c % 2 == 0:
            pb += exp(-r * alpha * c * c / 4.0) / r
        et = -r * alpha * ms_t * (c + ms_t)
        et = et[et > -40.0]
        if len(et) > 0:
            pb += (2.0 / r) * float(np.sum(np.exp(et)))
        # 1-step Abel bound at r:
        s2 = sin(pi * r * nu) ** 2
        if s2 > 1e-14:
            ab = (K_abel / s2) * min(1.0 / (r * c), sqrt(kap / (r * c)))
        else:
            ab = 1e9
        br = min(pb, ab)
        tot += br
        if br > max_term[0]: max_term = (br, r, pb, ab)
    # Add asymptotic tail for r > R_max: <= 2 / (c * R_max^2 * kap^2)
    tail = 2.0 / max(c * (R_max ** 2) * (kap ** 2), 1e-12)
    return tot + tail, tot, tail, max_term

print("=== TEST X = sum_{r>=1} min(pair_r, abel_r) WITH K_abel = 2.0 ===")
cases = [
    # Golden ratio / Fibonacci (should be < 1.50!)
    (89, 144, 0.1), (89, 144, 0.03), (89, 144, 0.01), (89, 144, 0.003),
    (200, 323, 0.1), (200, 323, 0.03), (200, 323, 0.01), (200, 323, 0.003), (200, 323, 0.001),
    (500, 809, 0.03), (500, 809, 0.01), (500, 809, 0.003), (500, 809, 0.001),
    # Small c with c*q*kap > 0.20
    (2, 3, 0.10), (2, 3, 0.30), (2, 3, 1.0),
    (3, 4, 0.05), (3, 4, 0.20), (3, 4, 0.80),
    (4, 5, 0.05), (4, 5, 0.20),
    (10, 11, 0.03), (10, 11, 0.10), (10, 11, 0.30),
    # Near major arcs
    (50, 51, 0.03), (50, 51, 0.10),
    (200, 205, 0.01), (200, 220, 0.01), (200, 241, 0.01), (200, 202, 0.01), (200, 201, 0.01),
    (200, 267, 0.01), (200, 299, 0.01),
]
for c, j, kap in cases:
    q = c // gcd(c, j)
    X, xs, tl, mt = compute_X(c, j, kap, K_abel=2.0, R_max=300)
    st = "PASS(<1.5)" if X < 1.50 else "MAJOR"
    print(f"c={c:3d} j={j:3d} q={q:3d} kap={kap:6.3f} (cqk={c*q*kap:7.2f}): X={X:7.3f} (tail={tl:.3f}, max_r={mt[1]:2d} val={mt[0]:.3f}) -> {st}")
