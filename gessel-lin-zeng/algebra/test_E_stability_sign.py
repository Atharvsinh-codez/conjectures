import math
import numpy as np

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

ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 14):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

def eval_c(p, z):
    acc = 0j
    for c in reversed(p): acc = acc * z + c
    return acc

test_ts = []
for re in [-20, -5, -1.5, -1.05, -0.95, -0.7, -0.5, -0.3, -0.1, -0.01, 0.0, 0.01, 0.1, 0.5, 1.0, 5.0, 20.0]:
    for im in [1e-5, 1e-3, 0.02, 0.2, 2.0]:
        test_ts.append(complex(re, im))

print("=== CHECK min(Im y) >= 0 (i.e. REAL STABILITY OF P(t, -y)) for binom-E, 1/i!-E, and Quad X - 2My + Yy^2 ===")
for k in range(3, 14):
    E = []
    for i in range(k):
        p = ht[(k-1, i)][:]
        for _ in range(i): p = apply_D(k + 1, p)
        E.append(p)
    # 1. binom E
    min_im_E = 1e9
    worst_E = None
    for tv in test_ts:
        yc = [math.comb(k-1, i) * eval_c(E[i], tv) for i in range(k-1, -1, -1)]
        mx = max(abs(c) for c in yc)
        yc = [c / mx for c in yc]
        mi = min(z.imag for z in np.roots(yc))
        if mi < min_im_E:
            min_im_E = mi; worst_E = tv
    # 2. Quad X - 2My + Yy^2 (test min Im y of X + 2My + Yy^2 = 0!)
    min_im_Q = 1e9
    worst_Q = None
    for i in range(2, k):
        X = ht[(k-1, i-2)]
        M = apply_D(k+1, ht[(k-1, i-1)])
        Y = apply_D(k+1, apply_D(k+1, ht[(k-1, i)]))
        for tv in test_ts:
            yc = [eval_c(Y, tv), 2*eval_c(M, tv), eval_c(X, tv)]
            mx = max(abs(c) for c in yc)
            yc = [c / mx for c in yc]
            mi = min(z.imag for z in np.roots(yc))
            if mi < min_im_Q:
                min_im_Q = mi; worst_Q = (i, tv)
    # 3. Also test 1/i! E (wait: fact(k-1)//fact(i) * E[i])
    min_im_F = 1e9
    for tv in test_ts:
        yc = [(math.factorial(k-1)//math.factorial(i)) * eval_c(E[i], tv) for i in range(k-1, -1, -1)]
        mx = max(abs(c) for c in yc)
        yc = [c / mx for c in yc]
        mi = min(z.imag for z in np.roots(yc))
        if mi < min_im_F: min_im_F = mi
    print(f"  k={k:2d}: binom-E minIm={min_im_E:+.3e} (at {worst_E}) | 1/i!-E minIm={min_im_F:+.3e} | Quad X-2My+Yy^2 minIm={min_im_Q:+.3e} (at {worst_Q})")
