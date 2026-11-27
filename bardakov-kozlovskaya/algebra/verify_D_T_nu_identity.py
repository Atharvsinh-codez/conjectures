"""
Verify the EXACT contour-rotated Poisson identity for b_r^G:
  b^G(c, D, gam) = (c/2) * sum_{k in Z} D(c*gam/2, |D + k*c|)
where for all T > 0 and integer nu >= 1:
  D(T, nu) = 2 int_0^inf e^{-pi*nu*v - (T/2)v^2} sin(T*v) dv
           - 2 Re int_0^inf e^{-(T - i*pi*nu)*y - (T/2)y^2} dy
           = 2 int_0^inf e^{-pi*nu*v} (e^{-(T/2)v^2} - 1) sin(T*v) dv
           + 2 Re int_0^inf e^{-(T - i*pi*nu)*y} (1 - e^{-(T/2)y^2}) dy  !!

WAIT A SECOND!! Look at the SECOND integral too:
Rotate y -> e^{i*theta} u? No, simply integrate by parts ONCE or TWICE, OR evaluate
the EXACT elementary bound on |D(T, nu)| in terms of (T, nu)!
Let's check BOTH the exact identity AND the exact 3-line elementary bound on |D(T, nu)|!
"""
import math
import cmath

def I1_direct(T, nu):
    # (-1)^nu int_{-1}^1 exp(-0.5*T*(1-t^2)) cos(pi*nu*t) dt
    N = 1000
    dt = 2.0 / N
    s = 0.0
    for i in range(N):
        t = -1.0 + (i + 0.5) * dt
        s += math.exp(-0.5 * T * (1.0 - t * t)) * math.cos(math.pi * nu * t) * dt
    return ((-1) ** nu) * s

def I1_rotated(T, nu):
    # 2 int_0^inf exp(-pi*nu*v - 0.5*T*v^2) sin(T*v) dv
    N = 1000
    vmax = 30.0 / (math.pi * nu)
    dv = vmax / N
    s = 0.0
    for i in range(N):
        v = (i + 0.5) * dv
        s += 2.0 * math.exp(-math.pi * nu * v - 0.5 * T * v * v) * math.sin(T * v) * dv
    return s

def I2_direct(T, nu):
    # 2 int_0^inf exp(-T*y - 0.5*T*y^2) cos(pi*nu*y) dy
    N = 2000
    ymax = min(40.0, 25.0 / max(T, 1e-4) + 10.0 / math.sqrt(max(T, 1e-4)))
    dy = ymax / N
    s = 0.0
    for i in range(N):
        y = (i + 0.5) * dy
        s += 2.0 * math.exp(-T * y - 0.5 * T * y * y) * math.cos(math.pi * nu * y) * dy
    return s

def D_T_nu(T, nu):
    return I1_rotated(T, nu) - I2_direct(T, nu)

def b_r_gauss(c, D, gam):
    th = 2.0 * math.pi * D / c
    B = sum(math.exp(-gam * m * (1.0 - m / c)) * math.cos(m * th) for m in range(1, c) if gam * m * (1.0 - m / c) < 40)
    A = 0.0; m = 1
    while True:
        ex = gam * m * (1.0 + m / c)
        if ex > 40: break
        A += 2.0 * math.exp(-ex) * math.cos(m * th); m += 1
    return B - A

print("=== 1. VERIFY EXACT ROTATED IDENTITY FOR I1(T, nu) ===")
for T, nu in [(0.5, 1), (2.0, 1), (5.0, 1), (2.0, 3), (10.0, 2)]:
    a = I1_direct(T, nu); b = I1_rotated(T, nu)
    print(f"  T={T:4.1f} nu={nu}: I1_direct={a:+.10f}, I1_rotated={b:+.10f}, err={abs(a-b):.2e}")

print("\n=== 2. VERIFY EXACT POISSON IDENTITY b^G(c, D, gam) = (c/2) sum_k D(c*gam/2, |D+kc|) ===")
for c, D, gam in [(6, 1, 0.5), (10, 3, 0.2), (15, 1, 0.35), (20, 1, 0.5), (29, 3, 0.15)]:
    lhs = b_r_gauss(c, D, gam)
    T = 0.5 * c * gam
    rhs = 0.5 * c * sum(D_T_nu(T, abs(D + k * c)) for k in range(-25, 26))
    print(f"  c={c:2d} D={D} gam={gam:.2f} (T={T:4.2f}): b^G={lhs:+.8f}, Poisson={rhs:+.8f}, err={abs(lhs-rhs):.2e}")

print("\n=== 3. CHECK EXACT 2-VARIABLE STRUCTURE OF D(T, nu) ===")
# Wait: compare D(T, nu) with its leading small-(T/nu) / large-T asymptotic:
# Remember: in (c/2)*D(T, nu), the main term is:
#   4*gam/c * Re[ 2 / (gam - i*2*pi*nu/c)^3 ] = 8*T / (pi^3 nu^3) * ... Wait:
#   Let t := T / (pi*nu) = c*gam / (2*pi*nu)!
#   Check the EXACT relation between D(T, nu) and h(t) = 8 t^2 (t^2 - 3) / (t^2 + 1)^3!
#   Notice (c/2)*D(T, nu) ~ (c / (4*pi^2 * nu^2)) * h(t)  <=>  D(T, nu) ~ h(t) / (2*pi^2 * nu^2) !!
#   WAIT A SECOND!! CHECK THAT!!
#   Compute:
#     nu^2 * |D(T, nu)|        (theory max: 3.03224 / (2*pi^2) = 0.15362!)
#     nu^2 * t * |D(T, nu)|    (theory max: 2.1558 / (2*pi^2) = 0.10922!)
#     nu^2 / t^2 * |D(T, nu)|  (theory max: 24 / (2*pi^2) = 1.21585!)
#     nu^2 * t^2 * |D(T, nu)|  (theory max: 8 / (2*pi^2) = 0.40528!)
mx_0 = (0.0, None); mx_1 = (0.0, None); mx_small = (0.0, None); mx_large = (0.0, None)
for nu in [1, 2, 3, 5, 10]:
    for t in [0.05, 0.15, 0.30, 0.50, 0.63, 0.90, 1.30, 2.0, 3.5, 8.0]:
        T = t * math.pi * nu
        val = abs(D_T_nu(T, nu))
        v0 = nu * nu * val
        v1 = nu * nu * t * val
        vs = nu * nu * val / (t * t)
        vl = nu * nu * t * t * val
        if v0 > mx_0[0]: mx_0 = (v0, (nu, t, T))
        if v1 > mx_1[0]: mx_1 = (v1, (nu, t, T))
        if vs > mx_small[0]: mx_small = (vs, (nu, t, T))
        if vl > mx_large[0]: mx_large = (vl, (nu, t, T))
print(f"  EXACT GLOBAL CONSTANTS FOR UNIVERSAL 2-VAR FUNCTION D(T, nu) (t = T/(pi*nu)):")
print(f"    max nu^2 * |D(T, nu)|       = {mx_0[0]:.5f} at {mx_0[1]}  (theory h(t)/(2pi^2): 0.15362)")
print(f"    max nu^2 * t * |D(T, nu)|   = {mx_1[0]:.5f} at {mx_1[1]}  (theory t*h(t)/(2pi^2): 0.10922)")
print(f"    max nu^2/t^2 * |D(T, nu)|   = {mx_small[0]:.5f} at {mx_small[1]}  (theory 24/(2pi^2): 1.21585)")
print(f"    max nu^2 * t^2 * |D(T, nu)| = {mx_large[0]:.5f} at {mx_large[1]}  (theory 8/(2pi^2): 0.40528)")
