"""
VERIFY THE TWO REMAINING ANALYTIC LEMMAS FOR GAP G = {c*kap^2 < 45, c*q*kap > 0.18}:
1. EXACT CLOSED-FORM 1-STEP PAIRED ABEL BOUND (Minor Arcs):
   Recall the EXACT identity (verified to 1e-15 in test_abel_exact.py!):
     4 * sin^2(pi * r * nu) * (B_r - A_r)
       = 2*Re(U_r(1) - V_r(1)) + Re sum_{m=1}^{c-1} Delta^2 U_r(m) w^m - Re sum_{m=1}^{n-c-1} Delta^2 V_r(m) w^m.
   Derive and test a 100% elementary closed-form upper bound on |B_r - A_r| / r of the form:
     |B_r - A_r| / r <= (K_abel / sin^2(pi*r*nu)) * min(1/(r*c), sqrt(kap/(r*c)))
   across ALL (n, c, j, r) in G!
2. EXACT BERNOULLI-B2 FOURIER COEFFICIENTS C_r(lam; D, tau) AND S_r(lam; D, tau) FOR MAJOR ARCS (D >= 1):
   Check how tight C_r vs C_r^main is across ALL D in {1, 2, 3, 5, 10} and tau in [0.2, 30.0]!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd
from scipy.integrate import quad

# =============================================================================
# 1. CHECK EXACT CONSTANT K_abel IN |B_r - A_r| / r <= (K_abel / sin^2(pi r nu)) * min(1/(r c), sqrt(kap/(r c)))
# =============================================================================
def check_K_abel(n, c, j, r):
    rj = sin(j * pi / n) / sin((j + c) * pi / n)
    sj = sin(c * pi / n) / sin((j + c) * pi / n)
    th_star = pi * (n - c) / n
    nu = sin(j * pi / n) * cos((j + c) * pi / n) / sin(c * pi / n)
    kap = 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)
    s2 = 4.0 * (sin(pi * r * nu) ** 2)
    if s2 < 1e-10:
        return None

    ms_e = np.arange(0, c + 1, dtype=np.float64)
    th_e = th_star + 2.0 * ms_e * pi / n
    log_ae_inv = n * np.log((1.0 + rj * np.exp(1j * th_e)) / (sj * np.exp(1j * j * pi / n)))
    U = np.exp(r * log_ae_inv + 2j * pi * r * nu * ms_e)

    M_i = min(n - c, int(c + np.ceil(np.sqrt(40.0 * c / max(r * kap, 1e-12)))) + 20)
    ms_i = np.arange(0, M_i + 1, dtype=np.float64)
    th_i = th_star - 2.0 * ms_i * pi / n
    log_ai = n * np.log((sj * np.exp(1j * j * pi / n)) / (1.0 + rj * np.exp(1j * th_i)))
    V = np.exp(r * log_ai + 2j * pi * r * nu * ms_i)

    d2U = U[2:] - 2.0 * U[1:-1] + U[:-2]
    d2V = V[2:] - 2.0 * V[1:-1] + V[:-2]
    bdy = 2.0 * abs(float((U[1] - V[1]).real))
    mh = (c - 1) // 2
    pair_d2 = 2.0 * float(np.sum(np.abs(d2U[:mh] - d2V[:mh])))
    mid_d2 = float(abs(d2U[mh])) if (c % 2 == 0) else 0.0
    tail_d2 = 2.0 * float(np.sum(np.abs(d2V[mh:])))
    num = bdy + pair_d2 + mid_d2 + tail_d2
    # Wait: scale is num / (4 * r) compared to min(1/(r*c), sqrt(kap/(r*c)))!
    # So (|B_r - A_r| / r) <= (num / (4*r)) / sin^2(pi*r*nu)!
    scale = min(1.0 / (r * c), sqrt(kap / (r * c)))
    K_eff = (num / (4.0 * r)) / scale
    return kap, num / (4.0 * r), scale, K_eff

print("=== 1. CHECK EFFECTIVE CONSTANT K_abel IN 1-STEP PAIRED ABEL BOUND ===")
max_K = 0.0
for n, c, j, r in [
    (1000, 10, 11, 1), (10000, 10, 11, 1), (10000, 10, 11, 2), (10000, 10, 11, 5),
    (10000, 50, 81, 1), (10000, 50, 81, 2), (10000, 50, 81, 5),
    (10000, 200, 323, 1), (10000, 200, 323, 2), (10000, 200, 323, 5),
    (100000, 200, 323, 1), (100000, 200, 323, 3), (100000, 200, 201, 1),
    (100000, 50, 81, 1), (100000, 50, 81, 5), (1000000, 200, 323, 1),
    (1000000, 500, 809, 1), (1000000, 500, 809, 3),
]:
    res = check_K_abel(n, c, j, r)
    if res is None: continue
    kap, val, sc, K = res
    max_K = max(max_K, K)
    print(f"  n={n:7d} c={c:3d} j={j:3d} r={r} (kap={kap:.4f}, r*c*k={r*c*kap:6.2f}): num/(4r)={val:.3e} scale={sc:.3e} -> K_eff={K:.3f}")
print(f"MAX K_eff across all 18 points = {max_K:.3f}!!")

# =============================================================================
# 2. CHECK BERNOULLI-B2 ACCURACY: C_r(lam; D, tau) vs C_r^main FOR D = 1, 2, 3, 5
# =============================================================================
print("\n=== 2. CHECK EXACT C_r VS BERNOULLI-B2 MAIN TERM C_r^main ===")
for D in [1, 2, 3, 5]:
    for tau in [0.5, 2.0, 5.0, 15.0]:
        # Choose delta = min(0.25, 1/(4*D)) for D>=2, or 0.25 for D=1
        delta = 0.25 if D == 1 else 1.0 / (4.0 * D)
        u0 = 0.5 - delta / D
        lam = exp(-tau * (0.25 - u0 * u0))
        ratios = []
        for r in [1, 2, 3]:
            # C_r = 2*int_0^{u0} exp(-2*r*tau*u0*t + r*tau*t^2) cos(2*pi*r*D*(u0 - t)) dt
            #     + 2*int_0^inf  exp(-2*r*tau*u0*t - r*tau*t^2) cos(2*pi*r*D*(u0 + t)) dt
            I1, _ = quad(lambda t: exp(-2*r*tau*u0*t + r*tau*t*t)*cos(2*pi*r*D*(u0 - t)), 0, u0, limit=200)
            L_out = sqrt(35.0 / (r * tau))
            I2, _ = quad(lambda t: exp(-2*r*tau*u0*t - r*tau*t*t)*cos(2*pi*r*D*(u0 + t)), 0, L_out, limit=200)
            C_r = 2.0 * (I1 + I2)
            C_main = ((-1)**(r*D)) * cos(2*pi*r*delta) * (2.0 * tau * u0) / (r * (tau*tau*u0*u0 + pi*pi*D*D))
            ratios.append((r, C_r, C_main, C_r / C_main if abs(C_main) > 1e-12 else 0.0))
        s_str = " | ".join(f"r={r}: C/Cm={rat:.4f} (Cm={cm:+.2e})" for r, cr, cm, rat in ratios)
        print(f"  D={D} tau={tau:4.1f} (del={delta:.3f}, u0={u0:.3f}): {s_str}")
