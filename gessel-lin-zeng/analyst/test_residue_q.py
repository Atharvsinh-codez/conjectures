"""CRITICAL ADVERSARIAL TEST of the residue inequality on the real zero curve Q_{k-1,q}(t0, y0) = 0:
  R(t0, y0) := t0*(1+t0) + y0 * Q_{k-1,q}(t0, q*y0) / Q_t(t0, y0) <= 0!
Remember:
  - Necessary for Im((D_{k+1} Q + y Q(t, qy)) / Q) >= 0 on H+ x H+!
  - Note: for EVERY real t0, Q_{k-1,q}(t0, y) has k-1 real roots y0 in R!
  - Check R(t0, y0) across ALL branches y0(t0) for t0 in [-1000, 1000], especially |t0| large AND t0 in (-1, 0)!
"""
import math
import numpy as np

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_Dtilde(j, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(poly):
        res[m] += c * (m + 1)
        res[m + 1] += c * (m + j + 1)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def build_htilde(kmax):
    p = [[None] * (k + 1) for k in range(kmax + 1)]
    p[1][0] = [1, 2]
    p[1][1] = [1]
    for k in range(2, kmax + 1):
        p[k][0] = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][0]))
        for i in range(1, k):
            p[k][i] = apply_Dtilde(k, poly_add(apply_Dtilde(k, p[k - 1][i]), p[k - 1][i - 1]))
        p[k][k] = apply_Dtilde(k, p[k - 1][k - 1])
    return p

H = build_htilde(10)

def eval_poly_and_deriv(coeffs, t0):
    val = 0.0; dval = 0.0
    for c in reversed(coeffs):
        dval = dval * t0 + val
        val = val * t0 + float(c)
    return val, dval

t_grid = [-1000.0, -100.0, -20.0, -5.0, -2.0, -1.2, -1.01, -0.99, -0.8, -0.5, -0.2, -0.05, -0.01,
          0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0, 1000.0]

for q in [0.5, 0.35, 0.25, 0.1, 0.01]:
    print(f"\n=== RESIDUE TEST R(t0, y0) = t0(1+t0) + y0*Q(t0,q*y0)/Q_t(t0,y0) FOR q = {q} ===")
    for k in range(2, 9):
        k1 = k - 1
        w = [q ** (i * (i - 1) // 2) for i in range(k1 + 1)]
        max_R = (-1e300, None)
        max_R_scaled = (-1e300, None)
        for t0 in t_grid:
            vals = [eval_poly_and_deriv(H[k1][i], t0) for i in range(k1 + 1)]
            ycoeffs = [w[i] * vals[i][0] for i in range(k1 + 1)]
            if abs(ycoeffs[-1]) < 1e-300: continue
            r_all = np.roots(list(reversed(ycoeffs)))
            if np.max(np.abs(r_all.imag)) > 1e-5 * (1 + np.max(np.abs(r_all))):
                print(f"  k={k}: NONREAL y-root at t0={t0}!")
                continue
            for y0 in sorted(r_all.real):
                Q_qy = sum(w[i] * vals[i][0] * ((q * y0) ** i) for i in range(k1 + 1))
                Qt = sum(w[i] * vals[i][1] * (y0 ** i) for i in range(k1 + 1))
                if abs(Qt) < 1e-300: continue
                R = t0 * (1.0 + t0) + y0 * Q_qy / Qt
                R_sc = R / (1.0 + t0 * t0)
                if R > max_R[0]:
                    max_R = (R, (t0, y0))
                if R_sc > max_R_scaled[0]:
                    max_R_scaled = (R_sc, (t0, y0))
        print(f"  k={k} (Q_{k1}): max R = {max_R[0]:+.4e} at t0={max_R[1][0]}, y0={max_R[1][1]:.3e} | "
              f"max R/(1+t^2) = {max_R_scaled[0]:+.4e} at t0={max_R_scaled[1][0]}")
