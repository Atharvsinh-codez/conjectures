"""
DEEP ADVERSARIAL CHECK OF CONJECTURED EXACT MONOMIAL TORUS POSITIVITY:
For even n = 2m and k1 >= k2 >= k3 >= 1, k1 + k2 + k3 = n, r_j = sin(pi*k_j/n):
  F_mono(phi2, psi) = (-1)^{k1*k2+k1*k3+k2*k3} * Re( p_n(1, r2^n e^{i phi2}, r3^n e^{i psi}) / (y^{k2} z^{k3}) ).

1. Check exact Hessian eigenvalues at the double point (phi2 = pi*k3, psi = pi*k2) for ALL (k1,k2,k3) at n=4..20!
2. Run dense grid + adversarial Nelder-Mead for n = 4, 6, 8, 10, 12, 14, 16, 18, 20!
   Wait: check BOTH (a) around the double point AND (b) across the whole [-pi, pi]^2!
3. ALSO test if the SAME inequality holds for ODD n (with (-1)^{k1*k2+k1*k3+k2*k3 + ...})!!
4. ALSO test what happens if we perturb r2 or r3 by 1%! (Does positivity ONLY hold at the exact sine radii r_j = sin(pi*k_j/n)?!)
"""
import numpy as np
from scipy.optimize import minimize
import sys
sys.path.insert(0, r'd:\solving-conjectures\work\bardakov-kozlovskaya\verify')
import engine

def eval_F_mono(n, k2, k3, r2, r3, ph2, ps, sgn, zetas, scale):
    y = (r2**n) * np.exp(1j * ph2)
    z = (r3**n) * np.exp(1j * ps)
    beta = r2 * np.exp(1j * ph2 / n)
    prod = np.prod(z - (1.0 + zetas * beta)**n)
    val = sgn * prod / (y**k2 * z**k3 * scale)
    return float(val.real)

print("=== 1. Hessian at double point + adversarial global min for even n = 4..20 ===")
for n in range(4, 22, 2):
    E = engine.compute_E_table(n)
    zetas = np.exp(2j * np.pi * np.arange(n) / n)
    min_global = 1e9
    min_hess_eig = 1e9
    worst_g = None
    worst_h = None
    n_cases = 0
    for k3 in range(1, n // 3 + 1):
        for k2 in range(k3, (n - k3) // 2 + 1):
            k1 = n - k2 - k3
            n_cases += 1
            M_val = abs(E[n - k3][k2])
            sgn = (-1)**(k1*k2 + k1*k3 + k2*k3)
            r2 = np.sin(np.pi * k2 / n) / np.sin(np.pi * k1 / n)
            r3 = np.sin(np.pi * k3 / n) / np.sin(np.pi * k1 / n)
            # Finite-difference Hessian at double point (dp2, dp3)
            dp2 = np.pi * (k3 % 2)
            dp3 = np.pi * (k2 % 2)
            h = 1e-4
            f0 = eval_F_mono(n, k2, k3, r2, r3, dp2, dp3, sgn, zetas, M_val)
            fxx = (eval_F_mono(n, k2, k3, r2, r3, dp2+h, dp3, sgn, zetas, M_val) - 2*f0 + eval_F_mono(n, k2, k3, r2, r3, dp2-h, dp3, sgn, zetas, M_val)) / h**2
            fyy = (eval_F_mono(n, k2, k3, r2, r3, dp2, dp3+h, sgn, zetas, M_val) - 2*f0 + eval_F_mono(n, k2, k3, r2, r3, dp2, dp3-h, sgn, zetas, M_val)) / h**2
            fxy = (eval_F_mono(n, k2, k3, r2, r3, dp2+h, dp3+h, sgn, zetas, M_val) - eval_F_mono(n, k2, k3, r2, r3, dp2+h, dp3-h, sgn, zetas, M_val) - eval_F_mono(n, k2, k3, r2, r3, dp2-h, dp3+h, sgn, zetas, M_val) + eval_F_mono(n, k2, k3, r2, r3, dp2-h, dp3-h, sgn, zetas, M_val)) / (4*h**2)
            eigs = np.linalg.eigvalsh(np.array([[fxx, fxy], [fxy, fyy]]))
            if eigs[0] < min_hess_eig:
                min_hess_eig = float(eigs[0])
                worst_h = (k1, k2, k3, float(eigs[0]), float(eigs[1]))
            # Dense 20x20 grid search + Nelder-Mead from worst 3 grid points!
            g2 = np.linspace(-np.pi, np.pi, 20, endpoint=False)
            g3 = np.linspace(-np.pi, np.pi, 20, endpoint=False)
            # Vectorized grid evaluation!
            P2, P3 = np.meshgrid(g2, g3, indexing='ij')
            p2f, p3f = P2.ravel(), P3.ravel()
            # Vectorized over grid points
            y_v = (r2**n) * np.exp(1j * p2f)
            z_v = (r3**n) * np.exp(1j * p3f)
            beta_v = r2 * np.exp(1j * p2f / n) # shape (400,)
            # (400, n)
            factors = z_v[:, None] - (1.0 + beta_v[:, None] * zetas[None, :])**n
            prods = np.prod(factors, axis=1)
            vals = (sgn * prods / (y_v**k2 * z_v**k3 * M_val)).real
            idx_min = int(np.argmin(vals))
            if vals[idx_min] < min_global:
                min_global = float(vals[idx_min])
                worst_g = (k1, k2, k3, float(p2f[idx_min]), float(p3f[idx_min]), float(vals[idx_min]))
            # Refine with Nelder-Mead from top 2 smallest grid points + slight nudge from dp
            for idx in np.argsort(vals)[:2]:
                res = minimize(lambda x: eval_F_mono(n, k2, k3, r2, r3, x[0], x[1], sgn, zetas, M_val),
                               [p2f[idx], p3f[idx]], method='Nelder-Mead', options={'maxiter':200})
                if res.fun < min_global:
                    min_global = float(res.fun)
                    worst_g = (k1, k2, k3, float(res.x[0]), float(res.x[1]), float(res.fun))
    print(f"n={n:2d} ({n_cases:2d} cases): min_Hess_eig={min_hess_eig:+.4e} at {worst_h[:3]} | min_global={min_global:+.6e} at k={worst_g[:3]} ang=({worst_g[3]:.3f},{worst_g[4]:.3f})")

print("\n=== 2. Check if min_global stays >= 0 when r2 is perturbed by +-2% (at n=6, k=(3,2,1)) ===")
n, k1, k2, k3 = 6, 3, 2, 1
E6 = engine.compute_E_table(6)
M_val = abs(E6[n - k3][k2])
sgn = (-1)**(k1*k2 + k1*k3 + k2*k3)
zetas = np.exp(2j * np.pi * np.arange(n) / n)
r2_0 = np.sin(np.pi * k2 / n) / np.sin(np.pi * k1 / n)
r3_0 = np.sin(np.pi * k3 / n) / np.sin(np.pi * k1 / n)
for mult in [0.95, 0.99, 1.00, 1.01, 1.05]:
    r2 = r2_0 * mult
    g2 = np.linspace(-np.pi, np.pi, 40, endpoint=False)
    P2, P3 = np.meshgrid(g2, g2, indexing='ij')
    p2f, p3f = P2.ravel(), P3.ravel()
    y_v = (r2**n)*np.exp(1j*p2f); z_v = (r3_0**n)*np.exp(1j*p3f)
    factors = z_v[:,None] - (1.0 + (r2*np.exp(1j*p2f/n))[:,None]*zetas[None,:])**n
    vals = (sgn * np.prod(factors, axis=1) / (y_v**k2 * z_v**k3 * M_val)).real
    print(f"  r2 mult={mult:.2f}: min Re = {vals.min():+.6e}")
