"""
VERIFY THE MIRACULOUS 4TH-ORDER PAIRED ABEL IDENTITY (EXACT CANCELLATION OF g''(0)!!):
For ANY sequence F_e(0..c) (symmetric: F_e(c-m) = F_e(m), F_e(0) = 1) and F_i(0..M) (F_i(0) = 1),
let G(m) := F_i(m) - F_e(m) for m = 0..c (G(0) = G(c) = 0? Wait: F_i(c) != 1, so do BOTH ends cleanly!):
Recall B_r - A_r = 2 * sum_{m=1}^{M} F_i(m) cos(m*th) - sum_{m=1}^{c-1} F_e(m) cos(m*th).
Write sum_{m=1}^{c-1} F_e(m) cos(m*th) with cos((c-m)*th) = cos(m*th):
  - Define g_full(m) := F_i(m) - F_e(m) for 1 <= m <= c-1, and F_i(m) for m >= c?
    No, F_e has factor 1 while F_i has factor 2!
  - Better: write 2*F_i(m) = F_i(m) + F_i(m)!
    Split at m = c/2? NO! Notice for ANY m in Z:
    sum_{m=1}^{c-1} F_e(m) cos(m th) = (1/2) sum_{m=1}^{c-1} [F_e(m) + F_e(c-m)] cos(m th).
    Even simpler! Apply the exact 1-line identity:
      16 sin^4(th/2) * sum_{m=1}^inf g(m) cos(m th) = (g(2) - 4 g(1)) + sum_{m=1}^inf d4g(m) cos(m th)
      (valid whenever g(0) = 0 and g(-1) := +g(1)!)
    DIRECTLY to:
      B_r - A_r = sum_{m=1}^inf [2*F_i(m) - 2*F_e_half(m)] cos(m th)  (Wait: at m = c/2, does F_e_half have a jump?!)
    Wait: check BOTH:
      (a) Full smooth F_i(m) - F_e(m) on 1..c-1 (using 2*F_i(m) = F_i(m) on m>=1 PLUS F_i(m) on m'=c+m, wait:
          cos((c+m)*th) = cos(m*th) because c*th = 2*pi*r*j in 2*pi*Z!!)
      OH!! MY!! GOD!!
      READ THAT AGAIN!!
      BECAUSE c * th = 2 * pi * r * j IS A MULTIPLE OF 2*pi:
      DEFINE A *SINGLE* FUNCTION H(m) FOR ALL m IN Z:
        - H(m) = +F_i(m) = +exp(-r*alpha*m*(c+m))   for m >= 1
        - H(0) = 0
        - H(m) = -F_e(-m) = -exp(-r*alpha*(-m)*(c+m))  for -c < m < 0
        - H(-c) = 0
        - H(m) = +F_i(-c-m) = +exp(-r*alpha*(-c-m)*(-m)) for m < -c!!
      HOLD EVERYTHING!!
      FOR ALL m IN Z, H(m) = sgn(m*(c+m)) * exp(-r*alpha*|m*(c+m)|)!!
      AND B_r - A_r = sum_{m in Z} H(m) cos(m*th)!!
      NOW JUST MULTIPLY BY (2 - 2*cos(th))^2 = 16 sin^4(th/2) ON ALL OF Z (ZERO BOUNDARY TERMS AT ALL!!):
        16 sin^4(th/2) * sum_{m in Z} H(m) cos(m*th) = sum_{m in Z} (Delta^4 H)(m) cos(m*th)!!
      CHECK (Delta^4 H)(0) AT THE CROSSING m = 0:
        H(m) near 0: H(0) = 0, H(1) = e^{-r*alpha*(c+1)}, H(-1) = -e^{-r*alpha*(c-1)},
        H(2) = e^{-2*r*alpha*(c+2)}, H(-2) = -e^{-2*r*alpha*(c-2)}!
        SO (Delta^4 H)(0) = (H(2) + H(-2)) - 4*(H(1) + H(-1)) + 6*0!!
        CHECK H(m) + H(-m) = e^{-r*kap*m} * (e^{-r*alpha*m^2} - e^{+r*alpha*m^2}) = -2 e^{-r*kap*m} sinh(r*alpha*m^2)!!
        SO (Delta^4 H)(0) = g(2) - 4 g(1) = O(r^2 * alpha * kap)!!
        AND (Delta^4 H)(m) = O(r^2 * alpha * kap) AT m = +-1 TOO, AND O((r*kap)^4) EVERYWHERE ELSE!!
        ZERO KINK AT m = c/2! ZERO BOUNDARY TERMS! ONE SINGLE LINE FOR ALL m IN Z!!

TEST THIS RIGHT NOW!!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd

def test_exact_Z_4th_abel(c, jm, kap, r):
    nu = jm / c
    th = 2.0 * pi * r * nu
    s2 = 4.0 * (sin(th / 2.0) ** 2)
    if s2 < 1e-12: return None
    alpha = kap / c
    M = int(np.ceil(sqrt(45.0 / max(r * alpha, 1e-10)))) + 15
    ms = np.arange(-(c + M), M + 1, dtype=np.float64)
    val = ms * (c + ms)
    H = np.where(val > 0, np.exp(-r * alpha * val),
        np.where(val < 0, -np.exp(r * alpha * val), 0.0))
    exact = float(np.sum(H * np.cos(ms * th)))
    # 1-step (2nd difference on Z, ZERO boundary terms!):
    d2H = H[2:] - 2.0 * H[1:-1] + H[:-2]
    rhs1 = -float(np.sum(d2H * np.cos(ms[1:-1] * th))) / s2
    ub1 = float(np.sum(np.abs(d2H))) / (s2 * r)
    # 2-step (4th difference on Z, ZERO boundary terms!):
    d4H = d2H[2:] - 2.0 * d2H[1:-1] + d2H[:-2]
    rhs2 = float(np.sum(d4H * np.cos(ms[2:-2] * th))) / (s2 * s2)
    ub2 = float(np.sum(np.abs(d4H))) / (s2 * s2 * r)
    # Check what makes up sum |d4H|: at m=0,-c vs m=+-1,-c+-1 vs smooth interior!
    return abs(exact) / r, ub1, ub2, abs(exact - rhs1), abs(exact - rhs2), float(np.sum(np.abs(d4H)))

print("=== 1. TEST 1-LINE Z-LATTICE 4TH-ORDER ABEL IDENTITY (ZERO BOUNDARY TERMS!) ===")
for c, jm, kap, r in [
    (100, 9, 0.002, 11), (100, 9, 0.002, 1), (100, 3, 0.005, 1), (100, 3, 0.02, 1),
    (50, 1, 0.02, 1), (50, 5, 0.02, 1), (30, 1, 0.05, 1), (30, 4, 0.05, 1),
    (20, 1, 0.10, 1), (10, 1, 0.15, 1), (6, 1, 0.25, 1), (200, 123, 0.005, 1),
    (100, 1, 0.01, 1), (100, 2, 0.01, 1), (100, 3, 0.01, 1), (100, 5, 0.01, 1),
]:
    ex, u1, u2, e1, e2, sd4 = test_exact_Z_4th_abel(c, jm, kap, r)
    # Compare sum|d4H| / r against theoretical scaling!
    tau = r * c * kap
    print(f"  c={c:3d} jm={jm:3d} kap={kap:.3f} r={r:2d} (tau={tau:4.2f}): ex={ex:.5f} | ub1={u1:.4f} | ub2={u2:.5f} | id_err=({e1:.1e},{e2:.1e}) | sd4={sd4:.3e}")
