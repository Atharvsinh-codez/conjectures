"""THE MOMENT OF TRUTH FOR A_k(t, y) = sum_{i=0}^k binom(k,i) (D_{k+1}^i htilde_{k,i}(t)) y^i:
1. Check Newton inequality g_{k,i}(t)^2 - g_{k,i-1}(t) g_{k,i+1}(t) > 0 on ALL t > 0 (and x-coeff positivity!)
   at k = 10, 20, 35, 40, 50 (where unaligned A_{35,17}^2 - A_{35,16}A_{35,18} FAILS at x = +7e-17!!).
2. Check B_k(t, t*y) stability vs A_{k, m}(t, y) for m = k+1 vs m = 10*k vs m = inf!
3. Check the EXACT small-t expansion (up to t^6) of g_{k,i}(t)^2 - g_{k,i-1}(t)g_{k,i+1}(t) up to k = 80!
4. Check full x-coefficient positivity of (S_L^i A_{k,i})^2 - (S_L^{i-1} A_{k,i-1})(S_L^{i+1} A_{k,i+1})?
   Wait: in x-space, what IS D_{k+1}^i htilde_{k,i}(t)?!
   Remember: htilde_{k,i}(t) = (1+t)^{2k-i-1} A_{k,i}(t/(1+t)), and D_m((1+t)^d f(t/(1+t))) = (1+t)^{d+1} S_{m+d-1}(f)(t/(1+t))!
   So D_{k+1}^i htilde_{k,i}(t) = (1+t)^{2k-1} [ S_{3k-1} S_{3k-2} ... S_{3k-i} A_{k,i} ](x)!!
   Let's check ALL x-coefficients of g_{k,i}_x(x)^2 - g_{k,i-1}_x(x)*g_{k,i+1}_x(x) for k = 2..36!!
"""
import math
import mpmath as mp

# === Exact x-space operator S_L(f) = (1 + L*x)*f + x*(1-x)*f' ===
# [x^r] S_L(f) = (r + 1) f_r + (L - (r - 1)) f_{r-1}
def apply_S_int(L, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for r, c in enumerate(poly):
        res[r] += c * (r + 1)
        res[r + 1] += c * (L - r)
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def poly_mul(a, b):
    res = [0] * (len(a) + len(b) - 1)
    for i, va in enumerate(a):
        if va == 0: continue
        for j, vb in enumerate(b):
            res[i + j] += va * vb
    return res

def poly_sub(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] -= v
    return res

# Build exact A_{k, i}(x) up to k = 36 (fast pure Python big ints!)
# Recurrence: A_{0,0} = [1], A_{k,i} = S_{3k-2-i}( S_{3k-3-i} A_{k-1,i} + A_{k-1,i-1} )
A = [[None]*(k+1) for k in range(37)]
A[0][0] = [1]
for k in range(1, 37):
    for i in range(k + 1):
        u = [0]
        if i <= k - 1:
            u = poly_add(u, apply_S_int(3*k - 3 - i, A[k-1][i]))
        if i >= 1:
            u = poly_add(u, A[k-1][i-1])
        A[k][i] = apply_S_int(3*k - 2 - i, u)

# Sanity check: compare g_{k,i}_x = S_{3k-1} ... S_{3k-i} A_{k,i} with D_{k+1}^i htilde_{k,i} at k = 5!
# Remember: deg A_{k,i} = 2k - i - 1; htilde_{k,i}(t) = (1+t)^{2k-i-1} A_{k,i}(t/(1+t)).
# Let's check: D_m((1+t)^d f(x)) where x = t/(1+t), 1+t = 1/(1-x), t = x/(1-x), t(1+t)d/dt = x d/dx:
#   t(1+t) d/dt [(1-x)^{-d} f(x)] = (1-x)^{-d} [ x f'(x) + d x/(1-x) f(x) ]
#   (1 + m*t) (1-x)^{-d} f(x)     = (1-x)^{-d} [ (1 + (m-1)x)/(1-x) f(x) ]
#   Sum = (1-x)^{-(d+1)} [ (1 + (m + d - 1)x) f(x) + x(1-x) f'(x) ] = (1+t)^{d+1} S_{m+d-1} f(x)!!
# For m = k+1, d = 2k - i - 1: m + d - 1 = 3k - i - 1!!
# Wait: CHECK THAT! For i=1: d = 2k-2, m+d-1 = (k+1) + (2k-2) - 1 = 3k - 2!!
# Next step: d = 2k-1, m+d-1 = 3k - 1... wait: for i=2, first step has d=2k-3 => m+d-1 = 3k-3, second step has d=2k-2 => 3k-2!
# In general, applying D_{k+1} r times to htilde_{k,i} applies S_{3k-i-1}, then S_{3k-i}, ..., up to S_{3k-i+r-2}!
# So for r = i: applies S_{3k-i-1}, S_{3k-i}, ..., S_{3k-2}!!
# Let's verify at k = 4, i = 2 against direct D_5^2 htilde_{4,2}!
def apply_D_int(m, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + 1)
        res[r + 1] += coef * (r + m)
    return res

def to_t_from_x(f_x, deg):
    # (1+t)^deg f_x(t/(1+t)) = sum f_r t^r (1+t)^{deg-r}
    res = [0] * (deg + 1)
    for r, c in enumerate(f_x):
        for j in range(deg - r + 1):
            res[r + j] += c * math.comb(deg - r, j)
    return res

ht42 = to_t_from_x(A[4][2], 2*4 - 2 - 1)
lhs = apply_D_int(5, apply_D_int(5, ht42))
gx42 = apply_S_int(3*4 - 2, apply_S_int(3*4 - 3, A[4][2]))
rhs = to_t_from_x(gx42, 2*4 - 1)
assert lhs == rhs, f"MISMATCH:\n{lhs}\n{rhs}"
print("[EXACT IDENTITY VERIFIED: D_{k+1}^i htilde_{k,i}(t) <-> S_{3k-2} ... S_{3k-i-1} A_{k,i}(x)]", flush=True)

print("\n=== 1. CHECK x-COEFF POSITIVITY OF g_{k,i}_x(x)^2 - g_{k,i-1}_x(x)*g_{k,i+1}_x(x) FOR k = 2..36 ===")
# Also compare against unaligned A_{k,i}(x)^2 - A_{k,i-1}(x)*A_{k,i+1}(x) (which we know fails at k=19, i=8!)
aligned_neg_coeffs = []
unaligned_neg_coeffs = []
for k in range(2, 37):
    # Compute g_{k, i}_x for all i = 0..k
    G = []
    for i in range(k + 1):
        p = A[k][i]
        for s in range(i):
            p = apply_S_int(3*k - i - 1 + s, p)
        G.append(p)
    for i in range(1, k):
        diff_g = poly_sub(poly_mul(G[i], G[i]), poly_mul(G[i-1], G[i+1]))
        neg_g = [(r, v) for r, v in enumerate(diff_g) if v <= 0]
        if neg_g:
            aligned_neg_coeffs.append((k, i, len(neg_g), neg_g[0]))
        if k in [19, 35] and i in [8, 17]:
            diff_a = poly_sub(poly_mul(A[k][i], A[k][i]), poly_mul(A[k][i-1], A[k][i+1]))
            neg_a = [(r, v) for r, v in enumerate(diff_a) if v <= 0]
            if neg_a:
                unaligned_neg_coeffs.append((k, i, len(neg_a), neg_a[0][0]))
    if k in [10, 20, 30, 35, 36]:
        print(f"  through k={k:2d}: aligned g_i^2 - g_{'{i-1}'}g_{'{i+1}'} nonpos x-coeffs = {len(aligned_neg_coeffs)}", flush=True)
print(f"  Sanity unaligned fails (k=19,35): {unaligned_neg_coeffs}")
print(f"  TOTAL aligned g_i^2 - g_{'{i-1}'}g_{'{i+1}'} nonpos x-coeffs (k=2..36): {len(aligned_neg_coeffs)} {aligned_neg_coeffs[:5]}", flush=True)

print("\n=== 2. CHECK HIGHER m-ALIGNED (m = k+1 vs 2k vs 10k vs inf) AT k = 35, i = 17, x = +7e-17 ===")
k = 35; i = 17
mp.mp.dps = 60
def eval_mp(p, xv):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * xv + c
    return acc

for m_val in [k+1, k+2, 2*k, 5*k, 20*k, 100*k, 10**6]:
    # Applying D_m r times to htilde_{k,j} (deg 2k-j-1) applies S_{m + 2k - j - 2}, ..., S_{m + 2k - j - 2 + r - 1}
    # Wait: check at j=i-1, i, i+1!
    G3 = []
    for j in [i-1, i, i+1]:
        p = A[k][j]
        for s in range(j):
            p = apply_S_int(m_val + (2*k - j - 1) - 1 + s, p)
        G3.append(p)
    diff_g = poly_sub(poly_mul(G3[1], G3[1]), poly_mul(G3[0], G3[2]))
    neg_cnt = sum(1 for v in diff_g if v <= 0)
    # Check min on x in [0, 1e-12]
    min_val = min(eval_mp(diff_g, mp.mpf(s)) for s in ["1e-18", "3e-17", "7e-17", "1e-16", "5e-16", "1e-15", "1e-14"])
    print(f"  k=35, i=17, m={m_val:7d}: neg_x_coeffs={neg_cnt:3d}, min near 7e-17 = {mp.nstr(min_val, 5)}", flush=True)

# Also check B_{35}(x, x*y) (i.e. m = inf) at x = +7e-17 directly:
diff_inf = poly_sub(poly_mul(A[35][17], A[35][17]), poly_mul(A[35][16], A[35][18]))
val_inf = eval_mp(diff_inf, mp.mpf("7e-17"))
print(f"  k=35, i=17, m=    inf: neg_x_coeffs={sum(1 for v in diff_inf if v<=0):3d}, val at 7e-17 = {mp.nstr(val_inf, 5)}", flush=True)
