import numpy as np
import sympy as sp

# ------------------------------------------------------------
# Part 1: Exact GLZ A_{k,i}(x) and ALL 4 cross-pairs + minimal c*
# ------------------------------------------------------------
def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0]*n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def poly_deriv(p):
    if len(p) <= 1: return [0]
    return [i * p[i] for i in range(1, len(p))]

def poly_mul_x(p):
    return [0] + p

def apply_S(p, L):
    # S_L f = (1 + L x) f + x(1-x) f'
    # coeff of x^n: p[n] + L p[n-1] + n p[n] - (n-1) p[n-1] = (n+1) p[n] + (L - n + 1) p[n-1]
    d = len(p) - 1
    res = [0]*(d + 2)
    for n in range(d + 2):
        v = 0
        if n <= d:
            v += (n + 1) * p[n]
        if 1 <= n <= d + 1:
            v += (L - n + 1) * p[n - 1]
        res[n] = v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_P(p, c):
    # P_c f = (c - x d/dx) f -> coeff of x^n is (c - n) p[n]
    return [(c - n) * p[n] for n in range(len(p))]

def roots_real(p):
    # p is ascending coeffs
    c = np.array(p[::-1], dtype=np.float64)
    r = np.roots(c)
    assert np.max(np.abs(r.imag)) < 1e-6, f"Complex root: {r}"
    return np.sort(r.real)

def min_c_interlace(g, f):
    """
    Find minimal c > deg(f) such that g << P_c f, assuming g << f (deg g = deg f - 1).
    At each root r of g, we need P_c f(r) to have the SAME sign as f(r)!
    Wait: W(g, P_c f) >= 0 everywhere on R iff P_c f(r_m) has same sign as f(r_m)
    at every root r_m of g (since g << f, f(r_m) alternates sign, and deg P_c f = deg f
    with same leading sign when c > deg f).
    Let's check: P_c f(r) = c f(r) - r f'(r).
    Same sign as f(r) iff (c f(r) - r f'(r)) / f(r) > 0 iff c > r f'(r) / f(r)!!
    WAIT A SECOND! Is that NECESSARY AND SUFFICIENT for g << P_c f (for c > deg f)?!
    YES! Because by Lagrange interpolation / Hermite-Biehler: if g has d-1 distinct roots r_1 < ... < r_{d-1},
    and P_c f has degree d with same leading sign as f and same sign as f at all r_m,
    then P_c f has at least d-1 sign changes on (r_1, r_{d-1}) and 1 sign change on (-inf, r_1) and 1 on (r_{d-1}, +inf) -> wait:
    that's d+1 sign changes? No, d-2 intervals (r_m, r_{m+1}) + 2 outer intervals = d intervals = d roots, strictly interlaced by r_1..r_{d-1}!
    """
    rg = roots_real(g)
    # evaluate f(r) and r*f'(r)
    fp = poly_deriv(f)
    c_req = len(f) - 1  # need c > deg(f)
    ratios = []
    for r in rg:
        fv = sum(f[n] * (r**n) for n in range(len(f)))
        fpv = sum(fp[n] * (r**n) for n in range(len(fp)))
        ratios.append(r * fpv / fv)
    return max(c_req, max(ratios)), max(ratios), ratios

# Compute GLZ A_{k,i} up to k=14
KMAX = 14
A = {(0, 0): [1]}
for k in range(1, KMAX + 1):
    for i in range(k + 1):
        u1 = apply_S(A[(k-1, i)], 3*k - 3 - i) if i <= k - 1 else [0]
        u2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S(poly_add(u1, u2), 3*k - 2 - i)

print("=== 1. EXACT MINIMAL c* FOR g << P_c f (where c_GLZ = 3k - 1 - i) ===")
print("Recall: g << P_c f (for c > deg f, given g << f) <==> c > max_{r: g(r)=0} r*f'(r)/f(r)!")
print("Checking (a) Main: A_{k,i} << P_c A_{k,i-1}")
print("Checking (b) Pair(2,1): A_{k-1,i-1} << P_c S_{3k-2-i} A_{k-1,i-1}")
print("Checking (c) Pair(1,1): S_{3k-3-i} A_{k-1,i} << P_c S_{3k-2-i} A_{k-1,i-1}")
print("Checking (d) Pair(1,2): S_{3k-3-i} A_{k-1,i} << P_c A_{k-1,i-2}")
print("Checking (e) Pair(2,2): A_{k-1,i-1} << P_c A_{k-1,i-2}")

for k in range(3, 13):
    max_ratio_main = (-1e9, None)
    max_ratio_21 = (-1e9, None)
    max_ratio_11 = (-1e9, None)
    max_ratio_12 = (-1e9, None)
    max_ratio_22 = (-1e9, None)
    for i in range(1, k):
        c_glz = 3*k - 1 - i
        # (a) Main
        if len(A[(k, i)]) > 1:
            cstar, rmax, _ = min_c_interlace(A[(k, i)], A[(k, i-1)])
            margin = c_glz - cstar
            if rmax - c_glz > max_ratio_main[0]:
                max_ratio_main = (rmax - c_glz, (i, c_glz, cstar, rmax, len(A[(k, i-1)])-1))
        # (b) Pair(2,1): A_{k-1,i-1} << P_c S_{3k-2-i} A_{k-1,i-1}
        if len(A[(k-1, i-1)]) > 1:
            U2 = A[(k-1, i-1)]
            f21 = apply_S(U2, 3*k - 2 - i)
            cstar, rmax, _ = min_c_interlace(U2, f21)
            if rmax - c_glz > max_ratio_21[0]:
                max_ratio_21 = (rmax - c_glz, (i, c_glz, cstar, rmax, len(f21)-1))
        # (c) Pair(1,1): S_{3k-3-i} A_{k-1,i} << P_c S_{3k-2-i} A_{k-1,i-1}
        if i <= k - 1:
            U1 = apply_S(A[(k-1, i)], 3*k - 3 - i)
            f11 = apply_S(A[(k-1, i-1)], 3*k - 2 - i)
            cstar, rmax, _ = min_c_interlace(U1, f11)
            if rmax - c_glz > max_ratio_11[0]:
                max_ratio_11 = (rmax - c_glz, (i, c_glz, cstar, rmax, len(f11)-1))
        # (d) Pair(1,2): S_{3k-3-i} A_{k-1,i} << P_c A_{k-1,i-2}
        if 2 <= i <= k - 1:
            U1 = apply_S(A[(k-1, i)], 3*k - 3 - i)
            f12 = A[(k-1, i-2)]
            cstar, rmax, _ = min_c_interlace(U1, f12)
            if rmax - c_glz > max_ratio_12[0]:
                max_ratio_12 = (rmax - c_glz, (i, c_glz, cstar, rmax, len(f12)-1))
        # (e) Pair(2,2): A_{k-1,i-1} << P_c A_{k-1,i-2}
        if 2 <= i <= k - 1 and len(A[(k-1, i-1)]) > 1:
            U2 = A[(k-1, i-1)]
            f22 = A[(k-1, i-2)]
            cstar, rmax, _ = min_c_interlace(U2, f22)
            if rmax - c_glz > max_ratio_22[0]:
                max_ratio_22 = (rmax - c_glz, (i, c_glz, cstar, rmax, len(f22)-1))
    print(f"k={k:2d}: "
          f"Main(r-c)={max_ratio_main[0]:+7.3f}@{max_ratio_main[1][0]} "
          f"| P21={max_ratio_21[0]:+7.3f}@{max_ratio_21[1][0]} "
          f"| P11={max_ratio_11[0]:+7.3f}@{max_ratio_11[1][0]} "
          f"| P12={max_ratio_12[0]:+7.3f}@{max_ratio_12[1][0] if max_ratio_12[1] else '-'} "
          f"| P22={max_ratio_22[0]:+7.3f}@{max_ratio_22[1][0] if max_ratio_22[1] else '-'}")

print("\n=== 2. WAIT A SECOND! COMPUTE r * (S_{L-1} f)'(r) / (S_{L-1} f)(r) AT A ROOT r OF f!! ===")
x, L, c = sp.symbols('x L c')
# At a root r of f:
# S_{L-1} f(x) = (1 + (L-1)x) f(x) + x(1-x) f'(x)
# Evaluate at x = r where f(r) = 0:
# (S_{L-1} f)(r) = r(1-r) f'(r)
# (S_{L-1} f)'(x) = (L-1)f(x) + (1 + (L-1)x)f'(x) + (1-2x)f'(x) + x(1-x)f''(x)
# At x = r (f(r)=0):
# (S_{L-1} f)'(r) = (2 + (L-3)r) f'(r) + r(1-r) f''(r)
# THEREFORE:
# r * (S_{L-1} f)'(r) / (S_{L-1} f)(r) = r * [ (2 + (L-3)r) f'(r) + r(1-r) f''(r) ] / [ r(1-r) f'(r) ]
#                                      = (2 + (L-3)r) / (1 - r) + r * f''(r) / f'(r) !!
# CHECK c - r * (S_{L-1} f)'(r) / (S_{L-1} f)(r) when c = L:
expr = sp.simplify(L - (2 + (L-3)*x)/(1 - x))
print("L - (2 + (L-3)r)/(1-r) =", sp.factor(expr))
print("So at a root r < 0 of f:")
print("  P_L(S_{L-1} f)(r) / (S_{L-1} f)(r) = [ L - 2 - (2L - 3)r ] / (1 - r) - r * f''(r)/f'(r)")
print("  WAIT A SECOND! For r < 0, (1-r) > 1, so [L - 2 + (2L-3)|r|]/(1+|r|) IS ALWAYS POSITIVE (~ 2L-3 at -inf, L-2 at 0)!")
print("  AND what is - r * f''(r) / f'(r)?!")
print("  For f(x) = prod(x - r_j), f''(r_m)/f'(r_m) = 2 * sum_{j != m} 1 / (r_m - r_j)!")
print("  So -r_m * f''(r_m)/f'(r_m) = 2 * sum_{j != m} (-r_m) / (r_m - r_j) = 2 * sum_{j != m} |r_m| / (|r_j| - |r_m|)!")
print("  Wait: if r_j -> r_m from the LEFT (|r_j| > |r_m|), this is +infinity (helps)! From the RIGHT (|r_j| < |r_m|), it is -infinity (hurts)!")

print("\n=== 3. CHECKING P_c(S_{L-1} f)(r) / (S_{L-1} f)(r) FORMULA ON GLZ A_{k,i} ===")
for k in [5, 8, 11]:
    for i in [1, k//2, k-2]:
        f = A[(k-1, i-1)]
        if len(f) <= 1: continue
        Lval = 3*k - 1 - i
        sf = apply_S(f, Lval - 1)
        rf = roots_real(f)
        fp = poly_deriv(f)
        fpp = poly_deriv(fp)
        vals = []
        for r in rf:
            term1 = (Lval - 2 - (2*Lval - 3)*r) / (1 - r)
            fv_p = sum(fp[n]*r**n for n in range(len(fp)))
            fv_pp = sum(fpp[n]*r**n for n in range(len(fpp)))
            term2 = -r * fv_pp / fv_p
            vals.append((r, term1, term2, term1 + term2))
        worst = min(vals, key=lambda z: z[3])
        print(f"  k={k:2d}, i={i}: L={Lval:2d}, deg={len(f)-1:2d}, worst root r={worst[0]:8.4f}: term1={worst[1]:+7.3f}, -r*f''/f'={worst[2]:+7.3f}, total={worst[3]:+7.3f}")
