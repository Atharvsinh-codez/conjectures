"""
Attack the next boundary strips algebraically (c = 2, a = 1, b = 1) AND inspect the e2 = -1 chart around (x, y, z) = (1, -1, 0):
1. c = 2 exact structure:
   We have Q_2(s, p) from z^2 coeff of p_{2m}(x, y, z)!
   Wait: check if Q_2(1, u) (all roots in (0, 1/4)) has a clean exact algebraic expression or sign-change proof!
2. a = 1 exact structure (odd c, ALL coefficients strictly negative!):
   Recall from Depressed-Cubic Resultant (or Lucas Resultant):
   At e1 = 0, a = 0 gave P_{2m}(0, e2, e3).
   Compute exact residue/derivative formula for a = 1: S_1(e2, e3) = (d/de1 P_{2m})|_{e1 = 0}!
3. e2 = -1 chart around (x, y, z) = (1, -1, 0):
   Write x = 1, y = -1 + delta, z = eps (so e1 = delta + eps, e2 = -1 + delta + delta*eps, e3 = (-1 + delta)*eps)!
   Check the exact Vandermonde-Cauchy integral on a torus around (delta, eps) = (0, 0) for small-a cases like (1, 10, 1), (2, 8, 2), (1, 7, 3) at n = 24!
"""
import math
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 50

# ============================================================
# 1. TEST THE e2 = -1 CHART AROUND (x, y, z) = (1, -1, 0) FOR SMALL-a CASES!
# ============================================================
# Let's derive the EXACT Cauchy integral in the e2 = 1 chart:
#   P_{2m}(e1, e2, e3) = sum c_{a,b,c} e1^a e2^b e3^c
#   Set e2 = -1 (just flips (-1)^b!):
#   P_{2m}(w1, -1, w3) = sum (-1)^b c_{a,b,c} w1^a w3^c!
#   So (-1)^b c_{a,b,c} = [w1^a w3^c] P_{2m}(w1, -1, w3)
#                       = (1/(2pi i)^2) oint oint P_{2m}(w1, -1, w3) / (w1^{a+1} w3^{c+1}) dw1 dw3!
# Now map (y, z) -> (w1, w3) via homogeneous scaling so e2 = -1!
# Wait: for x = 1, e1 = 1+y+z, e2 = y+z+yz, e3 = yz.
# Scaling (x, y, z) by lambda = (-e2)^{-1/2} gives:
#   w1 = e1 / (-e2)^{1/2},  w3 = e3 / (-e2)^{3/2},  P_{2m}(w1, -1, w3) = p_{2m}(1, y, z) / (-e2)^m!
# Compute Jacobian det d(w1, w3) / d(y, z) in exact SymPy!
y, z = sp.symbols('y z')
e1 = 1 + y + z
e2 = y + z + y * z
e3 = y * z
# Let s = (-e2)^{1/2}; w1 = e1 / s, w3 = e3 / s^3
# d(w1)/dy = (e1_y * (-e2) + (1/2) e1 e2_y) / s^3
# Let's compute directly in SymPy with s = sp.sqrt(-e2)!
s = sp.sqrt(-e2)
w1 = e1 / s
w3 = e3 / (s**3)
J = sp.simplify(sp.diff(w1, y) * sp.diff(w3, z) - sp.diff(w1, z) * sp.diff(w3, y))
print("=== 1. EXACT SYMPY JACOBIAN FOR e2 = -1 CHART (w1 = e1/(-e2)^{1/2}, w3 = e3/(-e2)^{3/2}) ===")
print("  det d(w1, w3)/d(y, z) =", sp.factor(J))
# Simplify J * (-e2)^3:
print("  J * (-e2)^3 =", sp.factor(sp.simplify(J * (-e2)**3)))

# Check the integrand in (y, z):
#   (-1)^b c_{a,b,c} = (1/(2pi i)^2) oint oint [p_{2m}(1, y, z) / (-e2)^m] / [ (e1/s)^{a+1} (e3/s^3)^{c+1} ] * J dy dz!
# Wait: exponent of s = (-e2)^{1/2}:
#   In denominator: 2m - (a+1) - 3(c+1) = 2m - a - 3c - 4 = 2b - 4!
#   And J has s^{-6} = (-e2)^{-3}! Wait: let's check exponent of (-e2) in J!
#   Total power of (-e2) in denominator: m - (a+1)/2 - 3(c+1)/2 + 3 = (2m - a - 3c - 4 + 6)/2 = b + 1!!
#   THE HALF-INTEGER POWERS OF (-e2) CANCEL IDENTICALLY!!
#   So the integrand is LITERALLY THE EXACT SAME 2-FORM:
#     (-1)^b c_{a,b,c} = (sign) * (1/(2pi i)^2) oint oint p_{2m}(1, y, z) * (1-y)(1-z)(y-z) / [ e1^{a+1} (-e2)^{b+1} e3^{c+1} ] dy dz!
#   Wait: check the sign and factor of 1/2! (Let's check numerically on n = 6 right now!)

# Let's test integrating on the torus y = -1 + r1*e^{i phi1} - r3*e^{i phi3}?
# Even simpler: directly parametrize (w1, w3) OR parametrize (u_1 = e1, u_3 = e3) around (y, z) = (-1, 0)!
# Notice near (y, z) = (-1, 0):
#   e1 = (1 + y) + z,   e3 = y * z!
#   If we pick |1+y| = R1, |z| = R3 with R3 < R1 << 1 (so |e1| ~ R1, |e3| ~ R3, |e2| ~ 1!),
#   the map (y, z) -> (e1, e3) has degree 1 around (y, z) = (-1, 0) (Jacobian at (-1,0) is -1 != 0)!
def eval_p2m(n, yv, zv):
    # prod_{k=0}^{n-1} (z - (1 + zeta^k y^{1/n})^n)
    r = abs(yv)**(1.0 / n)
    ang = np.angle(yv) / n
    prod = 1.0 + 0j
    for k in range(n):
        root = r * np.exp(1j * (ang + 2 * np.pi * k / n))
        prod *= (zv - (1.0 + root)**n)
    return prod

print("\n=== 2. NUMERICAL TEST OF e2 = -1 CHART ON n = 6 AND n = 24 SMALL-a CASES! ===")
# Check n = 6 first on ALL 7 coeffs with y = -1 + R1 e^{i t1}, z = R3 e^{i t3} (R3 < R1 << 1? Wait: we need |z| < |1+y| so y != z and e1 != 0!)
# Wait: dy dz = (i R1 e^{i t1})(i R3 e^{i t3}) dt1 dt3 = - R1 R3 e^{i(t1+t3)} dt1 dt3!
# So (1/(2pi i)^2) dy dz = (1/(2pi)^2) R1 R3 e^{i(t1+t3)} dt1 dt3!
def integral_e2_chart(n, a, b, c, R1, R3, N=64):
    t = np.linspace(0, 2*np.pi, N, endpoint=False)
    T1, T3 = np.meshgrid(t, t, indexing='ij')
    dy_fac = R1 * np.exp(1j * T1)
    dz_fac = R3 * np.exp(1j * T3)
    Y = -1.0 + dy_fac
    Z = dz_fac
    E1 = 1.0 + Y + Z
    E2 = Y + Z + Y * Z
    E3 = Y * Z
    # Vectorized p_n(1, Y, Z)
    rad = np.abs(Y)**(1.0 / n)
    ang = np.angle(Y) / n
    Pval = np.ones_like(Y, dtype=complex)
    for k in range(n):
        w = rad * np.exp(1j * (ang + 2 * np.pi * k / n))
        Pval *= (Z - (1.0 + w)**n)
    # Factor from J: check exact constant from J above!
    # We'll print ratio against exact coeff for n=6 first!
    V = (1.0 - Y) * (1.0 - Z) * (Y - Z)
    # Wait: dy dz / (2pi i)^2 = (R1 e^{iT1} R3 e^{iT3}) / (2pi)^2 dT1 dT3 (since i*i / i^2 = +1)!
    # Notice J had factor 1/2! Wait: why 1/2? Because w1 = e1/(-e2)^{1/2} is 2-to-1? No, check J * (-e2)^3!
    integrand = Pval * V * dy_fac * dz_fac / ((E1**(a + 1)) * (E2**(b + 1)) * (E3**(c + 1)))
    avg = np.mean(integrand)
    abs_avg = np.mean(np.abs(integrand))
    return avg, abs_avg

# Verify sign & constant at n = 6:
exact6 = {(6,0,0):1, (4,1,0):-12, (2,2,0):48, (0,3,0):-64, (3,0,1):-2754, (1,1,1):-12312, (0,0,2):185193}
for (a, b, c), ex in exact6.items():
    avg, abs_avg = integral_e2_chart(6, a, b, c, R1=0.25, R3=0.05, N=48)
    print(f"  n=6 ({a},{b},{c}): exact={ex:8d}, integral={avg.real:12.4f} (im={avg.imag:.1e}), ratio={avg.real/ex:.6f}, |avg|/avg|F|={abs(avg)/abs_avg:.4f}")

# Now test n = 24 small-a cases including (1, 10, 1), (2, 8, 2), (1, 7, 3), (0, 6, 4), (0, 9, 2), (1, 4, 5)!
print("\n  Testing n = 24 small-a cases in e2 = -1 chart:")
for (a, b, c) in [(1, 10, 1), (2, 8, 2), (1, 7, 3), (0, 9, 2), (0, 6, 4), (1, 4, 5), (2, 5, 4), (2, 2, 6)]:
    # Pick optimal radii R1, R3 by quick scan!
    best_r = (0.0, None, None, None)
    for R1 in [0.05, 0.15, 0.3, 0.5, 0.7]:
        for R3 in [1e-4, 1e-3, 0.01, 0.03, 0.08, 0.15]:
            if R3 >= 0.8 * R1: continue
            avg, abs_avg = integral_e2_chart(24, a, b, c, R1=R1, R3=R3, N=32)
            # Check if integral is clean (imag tiny)
            if abs_avg > 0 and abs(avg.imag) < 1e-3 * abs(avg.real):
                ratio = abs(avg) / abs_avg
                if ratio > best_r[0]:
                    best_r = (ratio, R1, R3, avg.real)
    print(f"  n=24 ({a},{b},{c}): best |avg|/avg|F| = {best_r[0]:.5f} at R1={best_r[1]}, R3={best_r[2]}, val={best_r[3]}")
