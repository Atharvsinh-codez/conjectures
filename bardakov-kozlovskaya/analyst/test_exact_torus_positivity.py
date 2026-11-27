"""
1. Adversarial Nelder-Mead & exact Hessian check of:
     F_mono(phi2, psi) = (-1)^{k1*k2+k1*k3+k2*k3} * Re( p_n(1, r2^n e^{i phi2}, r3^n e^{i psi}) / (y^{k2} z^{k3}) )
   at n = 4, 6, 8, 10, 12 (especially around the exact double point phi2 = pi*k3, psi = pi*k2)!
   Check if F_mono >= 0 is LITERALLY EXACT for ALL (phi2, psi) or if it dips slightly negative!
2. Also check the EXACT Hessian of F_mono at the double point for n=4 (2,1,1) in SymPy!
3. Also check: for n = 4, 6, 8, 10, 12, what is the EXACT minimum over ALL (a, b, c) with c>=1
   of the ratio |c_{a,b,c}| / (peak_integrand) and the EXACT closed-form asymptotic formula for c_{a,b,c}
   including the O(1) prefactor (from the Jacobi theta average!!)?!
"""
import sympy as sp
import numpy as np
from scipy.optimize import minimize
import mpmath as mp
import sys
sys.path.insert(0, r'd:\solving-conjectures\work\bardakov-kozlovskaya\verify')
import engine

# 1. Exact SymPy check for n = 4, (k1, k2, k3) = (2, 1, 1) around double point (-1/4, -1/4)!
u, v = sp.symbols('u v', real=True)
y = -sp.Rational(1, 4) * sp.exp(sp.I * u)
z = -sp.Rational(1, 4) * sp.exp(sp.I * v)
p4 = (1 + y**4 + z**4) - 4*(y + z + y**3 + z**3 + y**3*z + y*z**3) + 6*(y**2 + z**2 + y**2*z**2) - 124*(y*z + y**2*z + y*z**2)
expr = -p4 / (y * z)
# Take real part by expanding around u=0, v=0 to order 4!
re_series = sp.series(sp.series(sp.re(expr.rewrite(sp.cos)), u, 0, 5).removeO(), v, 0, 5).removeO()
print("=== 1. Exact Taylor expansion of -Re(p4(1, -e^{iu}/4, -e^{iv}/4) / (yz)) at double point ===")
# Faster: substitute cos/sin series directly!
y_s = -sp.Rational(1,4)*(1 + sp.I*u - u**2/2 - sp.I*u**3/6 + u**4/24)
z_s = -sp.Rational(1,4)*(1 + sp.I*v - v**2/2 - sp.I*v**3/6 + v**4/24)
def p4_fn(Y, Z):
    return (1 + Y**4 + Z**4) - 4*(Y + Z + Y**3 + Z**3 + Y**3*Z + Y*Z**3) + 6*(Y**2 + Z**2 + Y**2*Z**2) - 124*(Y*Z + Y**2*Z + Y*Z**2)
val_s = sp.expand(-p4_fn(y_s, z_s) * 16 * (1 - sp.I*u - u**2/2 + sp.I*u**3/6 + u**4/24) * (1 - sp.I*v - v**2/2 + sp.I*v**3/6 + v**4/24))
# Extract real part up to total degree 4
re_poly = sp.Poly(sp.re(val_s), u, v)
terms_le4 = sum(c * u**i * v**j for (i, j), c in re_poly.terms() if i + j <= 4)
print("  Re(-p4/(yz)) up to deg 4:", sp.factor(terms_le4))
# Also extract deg 2 specifically
deg2 = sum(c * u**i * v**j for (i, j), c in re_poly.terms() if i + j == 2)
print("  Quadratic form (deg 2):", sp.factor(deg2), "| expanded:", deg2)

# 2. Adversarial Nelder-Mead on F_mono(phi2, psi) for n = 4, 6, 8, 10, 12!
print("\n=== 2. Adversarial global minimization of F_mono(phi2, psi) for n = 4, 6, 8, 10, 12 ===")
for n in [4, 6, 8, 10, 12]:
    E = engine.compute_E_table(n)
    C = engine.compute_elementary_coeffs(n)
    min_val = 1e9
    worst_info = None
    for (a, b, c) in sorted(C.keys()):
        if c == 0: continue
        k1, k2, k3 = a+b+c, b+c, c
        M_val = ((-1)**(n - k3)) * E[n - k3][k2]
        sgn = 1 if M_val > 0 else -1
        r2 = np.sin(np.pi * k2 / n) / np.sin(np.pi * k1 / n)
        r3 = np.sin(np.pi * k3 / n) / np.sin(np.pi * k1 / n)
        Ry, Rz = r2**n, r3**n
        zetas = np.exp(2j * np.pi * np.arange(n) / n)
        def eval_f(ang):
            ph2, ps = ang[0], ang[1]
            y_v = Ry * np.exp(1j * ph2)
            z_v = Rz * np.exp(1j * ps)
            beta = r2 * np.exp(1j * ph2 / n)
            prod = np.prod(z_v - (1.0 + zetas * beta)**n)
            val = sgn * prod / (y_v**k2 * z_v**k3 * abs(M_val))
            return float(val.real)
        # Start near double point + 8 random starts
        dp = np.array([np.pi * (k3 % 2), np.pi * (k2 % 2)])
        starts = [dp + np.array([0.05, 0.03]), dp + np.array([0.05, -0.05]), dp + np.array([0.2, -0.1])]
        rng = np.random.default_rng(100*n + 10*b + c)
        for _ in range(6):
            starts.append(rng.uniform(-np.pi, np.pi, 2))
        for st in starts:
            res = minimize(eval_f, st, method='Nelder-Mead', options={'xatol':1e-9, 'fatol':1e-12, 'maxiter':400})
            if res.fun < min_val:
                min_val = res.fun
                worst_info = (a, b, c, float(res.x[0]), float(res.x[1]), float(res.fun))
    print(f"n={n:2d}: adversarial global min Re(sgn*Q/|M|) = {min_val:+.10e} at (a,b,c)={worst_info[:3]}, ang=({worst_info[3]:.4f},{worst_info[4]:.4f})")
