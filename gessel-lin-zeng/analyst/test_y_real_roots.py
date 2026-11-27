"""CRITICAL TEST:
For FIXED t in (-1, 0) (GLZ / Jacobi-Stirling) or FIXED t in (-inf, 0) (Laguerre):
Does H_k(t, y) = sum_{i=0}^k htilde_{k,i}(t) y^i have ALL REAL ROOTS IN y?!
AND does the Atomic Product prod_{j=1}^N (D_{c_j, m_j} + lam_j y)(1) have ALL REAL ROOTS IN y for every t in (-1, 0)?!
AND does the Laguerre Product prod_{j=1}^N (L_{c_j, m_j} + lam_j y)(1) have ALL REAL ROOTS IN y for every t < 0?!

WAIT: If YES, then by Newton's inequalities, htilde_{k,i}(t)^2 >= htilde_{k,i-1}(t)*htilde_{k,i+1}(t) on (-1, 0)
follows IMMEDIATELY, which by our 5-line Miracle Interlacing Lemma PROVES GLZ Conjecture 15!!
"""
import numpy as np
import mpmath as mp

def apply_D(c, m, poly):
    """D_{c,m} = t(1+t)d/dt + c + m*t."""
    d = len(poly) - 1
    res = [0.0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * (r + m)
    return res

def apply_L(c, m, poly):
    """L_{c,m} = t d/dt + c + m*t."""
    d = len(poly) - 1
    res = [0.0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * m
    return res

def poly_add_scaled(a, b, s=1.0):
    n = max(len(a), len(b))
    res = [0.0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += s * v
    return res

def eval_poly(p, tv):
    acc = 0.0
    for c in reversed(p):
        acc = acc * tv + c
    return acc

def max_imag_y_roots(polys, ts):
    """Given polys[i](t) for i=0..k, evaluate sum polys[i](t) y^i at each t in ts and return max |Im(y)|."""
    worst = (0.0, None)
    for tv in ts:
        coeffs = [eval_poly(p, tv) for p in polys]
        # Strip tiny leading? no, keep exact leading
        if abs(coeffs[-1]) < 1e-30:
            continue
        rts = np.roots(coeffs[::-1])
        mi = float(np.max(np.abs(rts.imag)))
        if mi > worst[0]:
            worst = (mi, (tv, rts))
    return worst

print("=== 1. GLZ H_k(t, y) = sum htilde_{k,i}(t) y^i: y-ROOTS FOR FIXED t in (-1, 0) ===")
# Build exact integer htilde_{k,i}(t) up to k = 15, then evaluate with mpmath at dense t in (-1, 0)!
def apply_D_int(m, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + 1)
        res[r + 1] += coef * (r + m)
    return res

def poly_add_int(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

KMAX = 14
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D_int(m, H[k-1][0])
    for i in range(1, k):
        U[i] = poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1):
        H[k][i] = apply_D_int(m, U[i])

mp.mp.dps = 50
ts_test = [mp.mpf(s) for s in [
    "-1e-12", "-1e-8", "-1e-5", "-1e-3", "-0.01", "-0.05", "-0.1", "-0.2",
    "-0.3", "-0.4", "-0.5", "-0.6", "-0.7", "-0.8", "-0.9", "-0.95", "-0.99", "-0.999", "-0.9999"
]]
for k in range(2, KMAX+1):
    max_im = mp.mpf(0)
    worst_t = None
    nonreal_count = 0
    for tv in ts_test:
        coeffs = []
        for i in range(k+1):
            # Horner in mpmath
            val = mp.mpf(0)
            for c in reversed(H[k][i]):
                val = val * tv + c
            coeffs.append(val)
        if abs(coeffs[-1]) < mp.mpf('1e-40'):
            continue
        rts = mp.polyroots(coeffs[::-1], maxsteps=200, extraprec=30)
        mi = max(abs(mp.im(r)) for r in rts)
        if mi > mp.mpf('1e-25'):
            nonreal_count += 1
        if mi > max_im:
            max_im = mi
            worst_t = tv
    print(f"  k={k:2d}: nonreal t-pts = {nonreal_count}/{len(ts_test)}, max |Im(y)| = {mp.nstr(max_im, 6)} at t = {mp.nstr(worst_t, 5)}")

print("\n=== 2. GLZ H_k(t, y) AT FIXED t > 0 (just to see contrast!) ===")
for k in [2, 3, 4, 5]:
    tv = mp.mpf("0.01")
    coeffs = []
    for i in range(k+1):
        val = mp.mpf(0)
        for c in reversed(H[k][i]): val = val * tv + c
        coeffs.append(val)
    rts = mp.polyroots(coeffs[::-1], maxsteps=200, extraprec=30)
    mi = max(abs(mp.im(r)) for r in rts)
    print(f"  k={k}, t=+0.01: max |Im(y)| = {mp.nstr(mi, 6)}")

print("\n=== 3. GENERAL ATOMIC PRODUCT prod_{j=1}^N (D_{c_j, m_j} + lam_j y)(1): y-ROOTS FOR FIXED t in (-1, 0) ===")
rng = np.random.default_rng(42)
ts_float = np.concatenate([-10.0**np.linspace(-6, -1, 15), np.linspace(-0.15, -0.95, 25), [-0.99, -0.999]])
atomic_fails = 0
for trial in range(40):
    N = int(rng.integers(3, 10))
    polys = [[1.0]]
    for j in range(N):
        c = float(10**rng.uniform(-1, 1))
        m = c + float(10**rng.uniform(-1, 2))
        lam = 0.0 if rng.random() < 0.25 else float(10**rng.uniform(-1, 1))
        new_polys = [[0.0]] * (len(polys) + 1)
        for i, p in enumerate(polys):
            new_polys[i] = poly_add_scaled(new_polys[i], apply_D(c, m, p), 1.0)
            if lam != 0.0:
                new_polys[i+1] = poly_add_scaled(new_polys[i+1], p, lam)
        # strip trailing all-zero polys
        while len(new_polys) > 1 and all(v == 0.0 for v in new_polys[-1]):
            new_polys.pop()
        polys = new_polys
    mi, info = max_imag_y_roots(polys, ts_float)
    if mi > 1e-6:
        atomic_fails += 1
        print(f"  ATOMIC FAIL trial {trial} (N={N}): max |Im(y)| = {mi:.3e} at t={info[0]}")
print(f"  General Atomic D_{{c,m}} y-root fails across 40 wild trials: {atomic_fails}/40")

print("\n=== 4. GENERAL LAGUERRE PRODUCT prod_{j=1}^N (L_{c_j, m_j} + lam_j y)(1): y-ROOTS FOR FIXED t < 0 ===")
ts_lag = -10.0**np.linspace(-5, 3, 40)
lag_fails = 0
for trial in range(40):
    N = int(rng.integers(3, 10))
    polys = [[1.0]]
    for j in range(N):
        c = float(10**rng.uniform(-1, 1))
        m = float(10**rng.uniform(-1, 1))
        lam = 0.0 if rng.random() < 0.25 else float(10**rng.uniform(-1, 1))
        new_polys = [[0.0]] * (len(polys) + 1)
        for i, p in enumerate(polys):
            new_polys[i] = poly_add_scaled(new_polys[i], apply_L(c, m, p), 1.0)
            if lam != 0.0:
                new_polys[i+1] = poly_add_scaled(new_polys[i+1], p, lam)
        while len(new_polys) > 1 and all(v == 0.0 for v in new_polys[-1]):
            new_polys.pop()
        polys = new_polys
    mi, info = max_imag_y_roots(polys, ts_lag)
    if mi > 1e-6:
        lag_fails += 1
        print(f"  LAGUERRE FAIL trial {trial} (N={N}): max |Im(y)| = {mi:.3e} at t={info[0]}")
print(f"  General Laguerre L_{{c,m}} y-root fails across 40 wild trials: {lag_fails}/40")
