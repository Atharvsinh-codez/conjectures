import numpy as np

# === TEST 1: Find exact threshold of x = sqrt(ac) for Psi(x, u, v, xi) >= 0 for ALL u, v in [0, 1], xi >= 0! ===
print("=== TEST 1: Check Psi(x, u, v, xi) >= 0 for all u, v in [0, 1], xi >= 0 ===")
# Note: Psi(x, u, v, xi) = 4(1 + 2u xi + xi^2)(1 + 2v x xi + x^2 xi^2) - (2u + (1 + 4uvx - x^2)xi + 2vx xi^2)^2
us = np.linspace(0.0, 1.0, 21)
vs = np.linspace(0.0, 1.0, 21)
xis = np.concatenate([[0.0], 10.0**np.linspace(-3, 3, 200)])

U_g, V_g = np.meshgrid(us, vs, indexing='ij') # (21, 21)

for x in [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 6.5, 7.0]:
    min_val = 1e99
    worst = None
    for xi in xis:
        lhs = 4.0 * (1.0 + 2.0*U_g*xi + xi*xi) * (1.0 + 2.0*V_g*x*xi + x*x*xi*xi)
        rhs = (2.0*U_g + (1.0 + 4.0*U_g*V_g*x - x*x)*xi + 2.0*V_g*x*xi*xi)**2
        diff = lhs - rhs
        m = float(np.min(diff))
        if m < min_val:
            idx = np.unravel_index(np.argmin(diff), diff.shape)
            min_val = m
            worst = (float(us[idx[0]]), float(vs[idx[1]]), float(xi))
    print(f"  x=sqrt(ac)={x:4.1f} (ac={x*x:5.2f}): min Psi = {min_val:+.6e} at (u,v,xi)={worst}")

# === TEST 2: Check for GLZ what max(W(A,B)*W(C,E) / W(B,C)^2) is on ALL of R for k=4..11! ===
print("\n=== TEST 2: For GLZ (k=4..11), check max of W(A,B)*W(C,E) / W(B,C)^2 and W(A,C) > 0, W(B,E) > 0 ===")
def poly_trim(c):
    while len(c) > 1 and c[-1] == 0: c.pop()
    return c
def apply_D(j, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return poly_trim(out)
def add_poly(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return poly_trim(out)
def sub_poly(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] -= b[i]
    return poly_trim(out)
def mul_poly(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b): out[i + j] += ai * bj
    return poly_trim(out)
def deriv(p):
    return poly_trim([m * p[m] for m in range(1, len(p))]) if len(p) > 1 else [0]
def wronsk(f, g):
    return sub_poly(mul_poly(f, deriv(g)), mul_poly(deriv(f), g))
def eval_safe(p, t_arr):
    # Scale by max abs coeff to avoid overflow; return (vals, log_scale)
    mx = max(abs(c) for c in p)
    pf = [c / mx for c in p]
    acc = np.zeros_like(t_arr, dtype=float)
    for c in reversed(pf): acc = acc * t_arr + c
    return acc, np.log(float(mx))

ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 12):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

t_grid = np.concatenate([np.linspace(-5.0, -1.0, 80), np.linspace(-1.0, 0.0, 301), np.linspace(0.0, 5.0, 80)])
for k in range(4, 12):
    max_ac = 0.0
    min_wac = 1e99
    min_wbe = 1e99
    max_at = None
    for i in range(3, k):
        A = ht[(k-1, i-3)]
        B = apply_D(k+1, ht[(k-1, i-2)])
        C = apply_D(k+1, apply_D(k+1, ht[(k-1, i-1)]))
        E = apply_D(k+1, apply_D(k+1, apply_D(k+1, ht[(k-1, i)])))
        wab, wbc, wce = wronsk(A, B), wronsk(B, C), wronsk(C, E)
        wac, wbe = wronsk(A, C), wronsk(B, E)
        v_ab, l_ab = eval_safe(wab, t_grid)
        v_bc, l_bc = eval_safe(wbc, t_grid)
        v_ce, l_ce = eval_safe(wce, t_grid)
        v_ac, _ = eval_safe(wac, t_grid)
        v_be, _ = eval_safe(wbe, t_grid)
        min_wac = min(min_wac, float(np.min(v_ac)))
        min_wbe = min(min_wbe, float(np.min(v_be)))
        # ac = (wab * wce) / wbc^2
        rat = (v_ab * v_ce) / (v_bc * v_bc) * np.exp(l_ab + l_ce - 2*l_bc)
        # Also at t -> infinity: lc(wab)*lc(wce) / lc(wbc)^2!
        rat_inf = (wab[-1] * wce[-1]) / (wbc[-1] ** 2)
        m = max(float(np.max(rat)), float(rat_inf))
        if m > max_ac:
            max_ac = m; max_at = i
    print(f"  k={k:2d}: max(W_AB*W_CE / W_BC^2) = {max_ac:.6f} (at i={max_at}) | W_AC>0: {min_wac>0} | W_BE>0: {min_wbe>0}")
