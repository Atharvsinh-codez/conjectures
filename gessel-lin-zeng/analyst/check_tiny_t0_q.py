"""Adversarial check of R(t0, y0) at TINY positive t0 (1e-8 .. 1e-1) for q = 0.5, 0.25, 0.1, 0.01!"""
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

H = build_htilde(9)

def eval_pd(coeffs, t0):
    val = 0.0; dval = 0.0
    for c in reversed(coeffs):
        dval = dval * t0 + val
        val = val * t0 + float(c)
    return val, dval

ts = [10**(-p) for p in np.linspace(1, 8, 29)] + [-1.0 - 10**(-p) for p in np.linspace(1, 8, 29)]
for q in [0.5, 0.25, 0.01]:
    print(f"\n=== TINY t0 > 0 and t0 < -1 CHECK FOR q = {q} ===")
    for k in range(2, 9):
        k1 = k - 1
        w = [q ** (i*(i-1)//2) for i in range(k1 + 1)]
        max_R = (-1e300, None)
        for t0 in ts:
            vals = [eval_pd(H[k1][i], t0) for i in range(k1 + 1)]
            ycoeffs = [w[i] * vals[i][0] for i in range(k1 + 1)]
            r_all = np.roots(list(reversed(ycoeffs)))
            for y0 in r_all.real:
                Q_qy = sum(w[i] * vals[i][0] * ((q * y0) ** i) for i in range(k1 + 1))
                Qt = sum(w[i] * vals[i][1] * (y0 ** i) for i in range(k1 + 1))
                R = t0 * (1.0 + t0) + y0 * Q_qy / Qt
                if R > max_R[0]:
                    max_R = (R, (t0, y0))
        print(f"  k={k}: max R = {max_R[0]:+.4e} at t0={max_R[1][0]:.3e}, y0={max_R[1][1]:.3e}")
