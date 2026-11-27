"""
TEST THE MIRACLE:
Compute F_n(phi2, phi3) := P_n(1, R2*e^{i*phi2}, R3*e^{i*phi3}) * e^{-i*b*phi2 - i*c*phi3} / (R2^b * R3^c * |c_{a,b,c}|)
on the ENTIRE torus (phi2, phi3) in [-pi, pi]^2 at the EXACT saddle point:
  theta_j = pi * k_j / n  (k1 = a+b+c, k2 = b+c, k3 = c),
  r2 = sin(theta_2)/sin(theta_1), r3 = sin(theta_3)/sin(theta_1),
  R2 = r2^n, R3 = (r2*r3)^n!

Check:
1. Does F_n(phi2, phi3) converge as n -> inf to a CLEAN O(1) function with NO rapid oscillations?!
2. What is its REAL PART vs IMAGINARY PART across [-pi, pi]^2?!
3. Does its exact average equal sign(c_{a,b,c}) = (-1)^{b(c+1)+c} with HUGE margin (NO cancellation)?!
4. ALSO check near the boundaries (small c, small b, small a) AND for ODD n!
"""
import mpmath as mp
import numpy as np
import sys
sys.path.insert(0, r'd:\solving-conjectures\work\bardakov-kozlovskaya\verify')
import engine

mp.mp.dps = 50

def eval_Pn_from_e(n, e2, e3):
    roots = mp.polyroots([mp.mpc(1), mp.mpc(-1), mp.mpc(e2), mp.mpc(-e3)], maxsteps=200, extraprec=50)
    x, y, z = roots
    alpha = x**(mp.mpf(1)/n)
    beta = y**(mp.mpf(1)/n)
    prod = mp.mpc(1)
    for k in range(n):
        prod *= (z - (alpha + mp.exp(2*mp.pi*1j*k/n)*beta)**n)
    return prod

# Test 1: Right triangle (k1, k2, k3) = (n/2, n/3, n/6) -> (a, b, c) = (n/6, n/6, n/6)
# Note for n=12: (2,2,2) [c even, b even -> +1]
# For n=18: (3,3,3) [c odd, b odd -> -1]
# For n=24: (4,4,4) [c even, b even -> +1]
# For n=30: (5,5,5) [c odd, b odd -> -1]
# Also test (a, b, c) with c even, b ODD -> sign = -1! E.g. n=18: (6, 3, 2) -> k=(11, 5, 2)
for n, a, b, c in [(12, 2, 2, 2), (18, 3, 3, 3), (24, 4, 4, 4), (30, 5, 5, 5), (18, 6, 3, 2), (24, 6, 5, 2)]:
    coeffs = engine.compute_elementary_coeffs(n)
    exact_c = coeffs[(a, b, c)]
    k1, k2, k3 = a + b + c, b + c, c
    r2 = mp.sin(mp.pi * k2 / n) / mp.sin(mp.pi * k1 / n)
    r3 = mp.sin(mp.pi * k3 / n) / mp.sin(mp.pi * k1 / n)
    R2 = r2**n
    R3 = (r2 * r3)**n
    scale = (R2**b) * (R3**c) * abs(exact_c)

    # Sample 16 x 16 grid on [-pi, pi)^2
    M = 16
    vals = []
    re_vals = []
    for j2 in range(M):
        phi2 = -mp.pi + 2*mp.pi*j2/M
        for j3 in range(M):
            phi3 = -mp.pi + 2*mp.pi*j3/M
            e2 = R2 * mp.exp(1j * phi2)
            e3 = R3 * mp.exp(1j * phi3)
            v = eval_Pn_from_e(n, e2, e3) * mp.exp(-1j * (b*phi2 + c*phi3)) / scale
            vals.append(v)
            re_vals.append(float(v.real))
    avg = sum(vals) / len(vals)
    abs_avg = sum(abs(v) for v in vals) / len(vals)
    min_re = min(re_vals)
    max_re = max(re_vals)
    sgn = 1 if exact_c > 0 else -1
    # Check how much cancellation there is: |avg| / avg(|F|)
    print(f"n={n:2d} (a,b,c)=({a},{b},{c}) [sgn={sgn:+d}]: avg={float(avg.real):+.6f}{float(avg.imag):+.1e}j | avg|F|={float(abs_avg):.4f} | Re(F) in [{min_re:+.3f}, {max_re:+.3f}] | ratio |avg|/avg|F|={float(abs(avg)/abs_avg):.4f}")
