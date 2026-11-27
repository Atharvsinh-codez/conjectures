"""
ADVERSARIAL OPTIMIZATION OF STAGE 3 DROP AT d=2 (h2=1, h1=t+a, h0=(t+b1)(t+b2))
AND d=3:
Maximize N^2 * (1 - min_t gamma_{j+1}(t) / min_t gamma_j(t)) subject to min_t gamma_j(t) >= 1!
Using EXACT algebraic minimum of gamma_j(t) at d=2!
"""
import math, random
import numpy as np
from scipy.optimize import minimize

# At d=2:
# X(t) = (t + b1)(t + b2) = t^2 + s*t + p   (where 0 < b1 < b2 < 1, roots -b2 < -b1)
# M_j(t) = D_j(t + a) = (j+1) t^2 + (2 + a*j) t + a
# Y_j(t) = D_j^2(1) = j(j+1) t^2 + 3*j t + 1
# Wronskian of two quadratics (A2 t^2 + A1 t + A0) and (B2 t^2 + B1 t + B0):
# W = (A2 B1 - A1 B2) t^2 + 2(A2 B0 - A0 B2) t + (A1 B0 - A0 B1)  (a QUADRATIC in t!)
# Wait: W(X, M) = X M' - X' M = (A2 t^2 + A1 t + A0)(2 B2 t + B1) - (2 A2 t + A1)(B2 t^2 + B1 t + B0)
#       = (A1 B2 - A2 B1) t^2 + 2(A0 B2 - A2 B0) t + (A0 B1 - A1 B0) !!

def W_quad(A, B):
    # A = (A0, A1, A2), B = (B0, B1, B2) -> returns (C0, C1, C2) for A B' - A' B
    return (A[0]*B[1] - A[1]*B[0], 2*(A[0]*B[2] - A[2]*B[0]), A[1]*B[2] - A[2]*B[1])

def min_gamma_d2(a, b1, b2, j):
    X = (b1*b2, b1+b2, 1.0)
    M = (a, 2.0 + a*j, j + 1.0)
    Y = (1.0, 3.0*j, j*(j + 1.0))
    w1 = W_quad(X, M)
    w2 = W_quad(M, Y)
    w3 = W_quad(X, Y)
    # Check w1 > 0 and w2 > 0 on all of R:
    if w1[2] <= 0 or w1[1]**2 - 4*w1[0]*w1[2] >= 0: return -1.0, 0.0
    if w2[2] <= 0 or w2[1]**2 - 4*w2[0]*w2[2] >= 0: return -1.0, 0.0
    # Minimize 4*w1(t)*w2(t) / w3(t)^2 over ALL t in R!
    # Derivative of (w1 w2 / w3^2) = 0 <=> (w1' w2 + w1 w2') w3 - 2 w1 w2 w3' = 0 (deg 5 polynomial!)
    p1 = np.array([w1[2], w1[1], w1[0]])
    p2 = np.array([w2[2], w2[1], w2[0]])
    p3 = np.array([w3[2], w3[1], w3[0]])
    p12 = np.polymul(p1, p2)
    dp12 = np.polyder(p12)
    dp3 = np.polyder(p3)
    crit_poly = np.polysub(np.polymul(dp12, p3), 2.0 * np.polymul(p12, dp3))
    rts = np.roots(crit_poly)
    cands = [z.real for z in rts if abs(z.imag) < 1e-7] + [-1e6, 0.0, 1e6]
    mn = 1e100
    wt = 0.0
    for tv in cands:
        v1 = np.polyval(p1, tv); v2 = np.polyval(p2, tv); v3 = np.polyval(p3, tv)
        val = (4.0 * v1 * v2) / (v3 * v3) if abs(v3) > 1e-30 else 1e100
        if val < mn: mn = val; wt = tv
    # Also t -> inf:
    val_inf = (4.0 * w1[2] * w2[2]) / (w3[2]**2) if abs(w3[2]) > 1e-30 else 1e100
    if val_inf < mn: mn = val_inf; wt = float('inf')
    return mn, wt

print("=== ADVERSARIAL SEARCH AT d=2 FOR MAX (j+1)^2 * (1 - gamma_{j+1}/gamma_j) ===", flush=True)
random.seed(7)
for j in [6, 12, 24, 48, 96]:
    N = j + 1
    best_sc = 0.0
    best_info = None
    # 50000 random samples focusing on boundaries (near 0, near 1, near D_j roots, barely stable)
    for tr in range(40000):
        # Roots of D_j(1) = -1/j; so a must be in (1/j, 1) so p1 << D_j p2!
        # Wait: root of p1 is -a, root of D_j(1) is -1/j, so -a < -1/j <=> a > 1/j!
        mode = tr % 4
        if mode == 0:
            a = random.uniform(1.0/j + 1e-4, 0.999)
        elif mode == 1:
            a = 1.0/j + 10**random.uniform(-5, -0.5) * (1 - 1.0/j)
        elif mode == 2:
            a = 1.0 - 10**random.uniform(-5, -0.5) * (1 - 1.0/j)
        else:
            a = random.uniform(1.0/j + 1e-4, min(0.999, 5.0/j))
        # Roots of M_j(t) = (j+1)t^2 + (2+a*j)t + a:
        disc_M = (2 + a*j)**2 - 4*(j+1)*a
        if disc_M <= 0: continue
        rM1 = (-(2 + a*j) - math.sqrt(disc_M)) / (2*(j+1))
        rM2 = (-(2 + a*j) + math.sqrt(disc_M)) / (2*(j+1))
        # Need p0 << M_j: -1 < -b2 < rM1 < -a < -b1 < rM2 < 0
        # i.e. b2 in (-rM1, 1), b1 in (-rM2, a)
        if mode == 1:
            b2 = -rM1 + 10**random.uniform(-5, -0.1) * (1 + rM1)
            b1 = -rM2 + 10**random.uniform(-5, -0.1) * (a + rM2)
        else:
            b2 = random.uniform(-rM1 + 1e-6, 0.9999)
            b1 = random.uniform(-rM2 + 1e-6, a - 1e-6)
        if not (0 < b1 < b2 < 1): continue
        g1, wt1 = min_gamma_d2(a, b1, b2, j)
        if g1 < 1.0: continue
        g2, wt2 = min_gamma_d2(a, b1, b2, j + 1)
        sc = (N**2) * (1.0 - g2 / g1)
        if sc > best_sc:
            best_sc = sc
            best_info = (a, b1, b2, g1, g2, wt1, wt2)
    # Local Nelder-Mead refinement around best!
    if best_info:
        a0, b10, b20 = best_info[:3]
        def obj(z):
            av, b1v, b2v = z
            if not (1.0/j < av < 1.0 and 0 < b1v < av and b1v < b2v < 1.0): return 1e6
            g1, _ = min_gamma_d2(av, b1v, b2v, j)
            if g1 < 1.0: return 1e6 + 1e4*(1.0 - g1)
            g2, _ = min_gamma_d2(av, b1v, b2v, j + 1)
            return -(N**2) * (1.0 - g2 / g1)
        res = minimize(obj, [a0, b10, b20], method='Nelder-Mead', options={'maxiter': 1000})
        if -res.fun > best_sc:
            av, b1v, b2v = res.x
            g1, wt1 = min_gamma_d2(av, b1v, b2v, j)
            g2, wt2 = min_gamma_d2(av, b1v, b2v, j + 1)
            if g1 >= 1.0:
                best_sc = (N**2)*(1.0 - g2/g1)
                best_info = (av, b1v, b2v, g1, g2, wt1, wt2)
    print(f"  j={j:2d} (N={N:2d}): MAX N^2*(1 - g_{{j+1}}/g_j) = {best_sc:.6f}  | a={best_info[0]:.5f}, b1={best_info[1]:.5f}, b2={best_info[2]:.5f}, g1={best_info[3]:.5f}@{best_info[5]:.3g} -> g2={best_info[4]:.5f}@{best_info[6]:.3g}", flush=True)
