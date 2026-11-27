import math
import numpy as np

def poly_trim(c):
    while len(c) > 1 and c[-1] == 0: c.pop()
    return c

def apply_S(L, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (L - m) * cm
    return poly_trim(out)

def add_poly(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return poly_trim(out)

def eval_c(p, z):
    acc = 0j
    for c in reversed(p): acc = acc * z + c
    return acc

A = {(0, 0): [1]}
for k in range(1, 17):
    for i in range(k + 1):
        t1 = apply_S(3*k - 3 - i, A[(k-1, i)]) if i <= k - 1 else [0]
        t2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S(3*k - 2 - i, add_poly(t1, t2))

# Test BOTH B_k(x, -x*y) = sum binom(k,i) A_{k,i}(x) (-x)^i y^i AND sum (1/i!) A_{k,i}(x) (-x)^i y^i
# AND test B_k(t, -t*y) = sum binom(k,i) ht_{k,i}(t) (-t)^i y^i!
# Remember: for real stability in (x, y), when Im(x) > 0 we need ALL y-roots to have Im(y) <= 0!
# Note: in sum w_i A_{k,i}(x) (-x)^i y^i = sum w_i A_{k,i}(x) (-x*y)^i:
# Let z = -x*y, so y = -z/x!
# We check max(Im y) <= 0 for x in UHP!
test_xs = []
for re in [-50, -10, -3, -1.0, -0.5, -0.1, -0.01, 0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0]:
    for im in [1e-6, 1e-4, 1e-2, 0.1, 1.0, 10.0]:
        test_xs.append(complex(re, im))

print("=== TEST 1: Real stability of B_k(x, -x*y) = sum binom(k,i) A_{k,i}(x) (-x)^i y^i (k=2..14) ===")
for k in range(2, 15):
    max_im = -1e9
    worst = None
    for xv in test_xs:
        yc = [math.comb(k, i) * eval_c(A[(k, i)], xv) * ((-xv)**i) for i in range(k, -1, -1)]
        mx = max(abs(c) for c in yc)
        if mx == 0: continue
        yc = [c / mx for c in yc]
        mi = max(r.imag for r in np.roots(yc))
        if mi > max_im:
            max_im = mi; worst = xv
    print(f"  k={k:2d}: B_k(x, -x*y) max Im(y) = {max_im:+.4e} at x={worst}")

print("\n=== TEST 2: Real stability of L_k(x, -x*y) = sum (k!/i!) A_{k,i}(x) (-x)^i y^i (k=2..14) ===")
for k in range(2, 15):
    max_im = -1e9
    worst = None
    for xv in test_xs:
        yc = [(math.factorial(k)//math.factorial(i)) * eval_c(A[(k, i)], xv) * ((-xv)**i) for i in range(k, -1, -1)]
        mx = max(abs(c) for c in yc)
        if mx == 0: continue
        yc = [c / mx for c in yc]
        mi = max(r.imag for r in np.roots(yc))
        if mi > max_im:
            max_im = mi; worst = xv
    print(f"  k={k:2d}: L_k(x, -x*y) max Im(y) = {max_im:+.4e} at x={worst}")

print("\n=== TEST 3: Real stability of B_k(t, -t*y) = sum binom(k,i) ht_{k,i}(t) (-t)^i y^i (k=2..14) ===")
def apply_D(j, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return poly_trim(out)
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 15):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

test_ts = [complex(re, im) for re in [-20, -5, -1.5, -1.05, -0.95, -0.5, -0.1, -0.01, 0.0, 0.01, 0.1, 0.5, 1.0, 5.0, 20.0] for im in [1e-5, 1e-3, 0.02, 0.2, 2.0]]
for k in range(2, 15):
    max_im = -1e9
    worst = None
    for tv in test_ts:
        yc = [math.comb(k, i) * eval_c(ht[(k, i)], tv) * ((-tv)**i) for i in range(k, -1, -1)]
        mx = max(abs(c) for c in yc)
        if mx == 0: continue
        yc = [c / mx for c in yc]
        mi = max(r.imag for r in np.roots(yc))
        if mi > max_im:
            max_im = mi; worst = tv
    print(f"  k={k:2d}: B_k^ht(t, -t*y) max Im(y) = {max_im:+.4e} at t={worst}")

print("\n=== TEST 4: Also test PLAIN (unweighted!) H_k(x, -x*y) = sum A_{k,i}(x) (-x)^i y^i and H_k^ht(t, -t*y)! ===")
for k in range(2, 8):
    max_im_x = -1e9
    for xv in test_xs:
        yc = [eval_c(A[(k, i)], xv) * ((-xv)**i) for i in range(k, -1, -1)]
        mx = max(abs(c) for c in yc)
        if mx == 0: continue
        yc = [c / mx for c in yc]
        max_im_x = max(max_im_x, max(r.imag for r in np.roots(yc)))
    print(f"  k={k}: plain sum A_{{k,i}}(x) (-x)^i y^i max Im(y) = {max_im_x:+.4e}")
