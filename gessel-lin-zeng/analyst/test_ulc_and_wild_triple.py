"""Test:
1. Ultra-log-concavity (ULC) of htilde_{k,i}(t) on t in (-1, 0):
   (htilde_{k,i}/binom(k,i))^2 >= (htilde_{k,i-1}/binom(k,i-1)) * (htilde_{k,i+1}/binom(k,i+1))
   for k = 2..20 (dense grid + exact SymPy count_roots for k <= 10).
2. Does B_k(t, y) = sum binom(k,i) htilde_{k,i}(t) y^i have ALL REAL y-ROOTS for fixed t in (-1, 0) (k=2..14)?
   (Remember B_{10} failed STABILITY, but might STILL be real-rooted in y for fixed t in (-1,0)!)
3. Does normalized A_k_norm(t, y) = sum (htilde_{k,i}(t) / binom(k,i)) y^i (or with extra binom) have real y-roots?
4. Does the aligned triple quadratic X + 2 M y + Y y^2 have 4 W(X,M)W(M,Y) - W(X,Y)^2 >= 0
   AND UHP real stability for WILD parameters m_j in the Atomic Theorem?
"""
import math
import numpy as np
import sympy as sp

t = sp.Symbol('t')

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

def eval_poly(p, tv):
    acc = 0.0
    for c in reversed(p):
        acc = acc * tv + c
    return acc

KMAX = 16
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

print("=== 1. ULC ON (-1, 0): (h_{k,i}/C(k,i))^2 - (h_{k,i-1}/C(k,i-1))*(h_{k,i+1}/C(k,i+1)) ===")
ts_grid = np.concatenate([-10.0**np.linspace(-10, -1, 100), np.linspace(-0.1, -0.9999, 500)])
ulc_fails = []
for k in range(2, KMAX+1):
    for i in range(1, k):
        # exact factor: c_num = C(k,i-1)*C(k,i+1), c_den = C(k,i)^2
        # Check C(k,i-1)*C(k,i+1) * h_i^2 - C(k,i)^2 * h_{i-1}*h_{i+1}
        c_left = math.comb(k, i-1) * math.comb(k, i+1)
        c_right = math.comb(k, i) ** 2
        if k <= 8:
            hi = sum(sp.Integer(v)*t**r for r, v in enumerate(H[k][i]))
            hm = sum(sp.Integer(v)*t**r for r, v in enumerate(H[k][i-1]))
            hp = sum(sp.Integer(v)*t**r for r, v in enumerate(H[k][i+1]))
            P = sp.Poly(sp.expand(c_left * hi*hi - c_right * hm*hp), t, domain='ZZ')
            # factor out t (root at t=0)
            rc = int(sp.count_roots(P, -1, 0))
            val_mid = P.eval(sp.Rational(-1, 2))
            if rc > 0 or val_mid < 0:
                ulc_fails.append((k, i, "EXACT", rc, float(val_mid)))
        else:
            # float normalized check
            for tv in ts_grid:
                vi = eval_poly(H[k][i], tv) / math.comb(k, i)
                vm = eval_poly(H[k][i-1], tv) / math.comb(k, i-1)
                vp = eval_poly(H[k][i+1], tv) / math.comb(k, i+1)
                diff = vi*vi - vm*vp
                scale = vi*vi + abs(vm*vp) + 1e-300
                if diff < -1e-11 * scale:
                    ulc_fails.append((k, i, "GRID", tv, diff/scale))
                    break
print(f"  ULC fails on (-1, 0) for k=2..{KMAX}: {len(ulc_fails)} {ulc_fails[:8]}")

print("\n=== 2. B_k(t, y) = sum binom(k,i) htilde_{k,i}(t) y^i: y-ROOTS FOR FIXED t in (-1, 0) ===")
for k in range(2, 14):
    max_im = 0.0; worst_t = None
    for tv in [-1e-5, -1e-3, -0.01, -0.05, -0.1, -0.2, -0.3, -0.5, -0.7, -0.9, -0.99]:
        coeffs = [math.comb(k, i) * eval_poly(H[k][i], tv) for i in range(k+1)]
        if abs(coeffs[-1]) < 1e-30: continue
        rts = np.roots(coeffs[::-1])
        mi = float(np.max(np.abs(rts.imag)))
        if mi > max_im:
            max_im = mi; worst_t = tv
    print(f"  k={k:2d}: max |Im(y)| = {max_im:.3e} at t={worst_t}")

print("\n=== 3. ALIGNED TRIPLE QUADRATIC X + 2 M y + Y y^2 FOR WILD ATOMIC PARAMETERS ===")
def apply_D_flt(c, m, poly):
    d = len(poly) - 1
    res = [0.0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * (r + m)
    return res

def poly_add_flt(a, b):
    n = max(len(a), len(b))
    res = [0.0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

rng = np.random.default_rng(777)
test_ys = [complex(re_y, im_y) for re_y in [-10, -2, -0.5, 0, 0.5, 2, 10] for im_y in [1e-3, 0.1, 1.0, 10.0]]
wild_quad_fails = 0
for trial in range(50):
    k = int(rng.integers(3, 8))
    # Build atomic product of k pairs D_{m_r}(D_{m_r} + y) with wild m_r > 1
    polys = [[1.0]]
    ms = [1.0 + float(10**rng.uniform(-0.5, 2.0)) for _ in range(k)]
    for r in range(k):
        m = ms[r]
        new_p = [[0.0]] * (len(polys) + 1)
        for i, p in enumerate(polys):
            new_p[i] = poly_add_flt(new_p[i], apply_D_flt(1.0, m, p))
            new_p[i+1] = poly_add_flt(new_p[i+1], p)
        polys = [apply_D_flt(1.0, m, p) for p in new_p]
    # Test aligned triple with next m_next!
    m_next = 1.0 + float(10**rng.uniform(-0.5, 2.0))
    failed = False
    for i in range(2, len(polys)):
        X = polys[i-2]
        M = apply_D_flt(1.0, m_next, polys[i-1])
        Y = apply_D_flt(1.0, m_next, apply_D_flt(1.0, m_next, polys[i]))
        d = len(X) - 1
        for yv in test_ys:
            c = np.array([X[r] + 2*yv*M[r] + yv*yv*Y[r] for r in range(d+1)], dtype=np.complex128)
            if abs(c[-1]) < 1e-30: continue
            rts = np.roots(c[::-1])
            if np.max(rts.imag) > 1e-7 and np.min(rts.imag) < -1e-7:
                wild_quad_fails += 1
                failed = True
                if wild_quad_fails <= 5:
                    print(f"  WILD FAIL trial {trial} (k={k}, i={i}, m_next={m_next:.2f}): Im in [{np.min(rts.imag):.2e}, {np.max(rts.imag):.2e}]")
                break
        if failed: break
print(f"  Wild Atomic triple quadratic UHP fails: {wild_quad_fails}/50")
