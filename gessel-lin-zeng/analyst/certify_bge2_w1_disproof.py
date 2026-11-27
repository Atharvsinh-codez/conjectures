import mpmath as mp
import sympy as sp
from fractions import Fraction
import math

# 1. Certify COMPLEX ROOTS (by exact SymPy real-root count / Sturm on exact Rationals!)
#    for SINGLE-operator prod_{j=1}^k (D_{b_j} + y)(1) with a_j = 1, w_j = 1, ALL b_j >= 2!
# 2. Test the EXACT SAME skewed-i attack (K=1000, i0=960 -> K-i0=40, w_eff=23.4)
#    on DOUBLE-operator prod_{j=1}^k D_{b_j}(D_{b_j} + y)(1) with a_j = 1, w_j = 1, ALL b_j >= 2!

def apply_D_int(b, p):
    # b is exact int, p is list of ints
    res = [0]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+b)
    return res

def add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

print("=== 1. EXACT INTEGER STURM CERTIFICATE FOR SINGLE-OPERATOR prod(D_{b_j}+y)(1), ALL b_j >= 2 (INTEGERS!) ===")
K = 1000; i0 = 960
# Precompute D_2^m(1) in EXACT integers for m = 0..44 (deg <= 44, takes 0.01s!)
D2 = [[1]]
for m in range(1, 45):
    D2.append(apply_D_int(2, D2[-1]))

# To keep integers smaller, divide C(K, r) by gcd of {C(K, r) : r in i0-3..i0+2}!
cs = {r: math.comb(K, r) for r in range(i0 - 3, i0 + 3)}
g = math.gcd(*cs.values())
P = {r: [ (cs[r] // g) * c for c in D2[K - r] ] for r in range(i0 - 3, i0 + 3)}

# Use our exact integer b1, b2 from the 120-dps hit!
b1 = 196200000000000  # 1.962e14 (EXACT INTEGER >= 2!)
b2 = 181330000000000  # 1.8133e14 (EXACT INTEGER >= 2!)

P1 = {r: add_int(apply_D_int(b1, P[r]), P[r-1]) for r in range(i0 - 2, i0 + 3)}
P2 = {r: add_int(apply_D_int(b2, P1[r]), P1[r-1]) for r in [i0 - 1, i0, i0 + 1]}

# Check exact interlacing break of P2[i0+1] and P2[i0], AND exact complex roots of P3[i0+1] = D_{b3} P2[i0+1] + P2[i0]!
t = sp.Symbol('t')
def count_real_roots_exact(coeffs):
    poly = sp.Poly(sum(c * t**r for r, c in enumerate(coeffs)), t, domain='ZZ')
    deg = poly.degree()
    # Make squarefree first if needed, or sp.count_roots handles multiplicities!
    n_real = poly.count_roots()
    sqf = (sp.gcd(poly, poly.diff()).degree() == 0)
    return deg, n_real, sqf

# First verify P2[i0] and P2[i0+1] are still real-rooted (or already complex!)
print("  Checking P2[i0] and P2[i0+1] exact real-root counts:")
for r in [i0, i0+1]:
    d, nr, sqf = count_real_roots_exact(P2[r])
    print(f"    P2[{r}]: deg={d}, real_roots={nr}, sqf={sqf}")

# Check exact sign of W(P2[i0+1], P2[i0]) at t_neg = -0.2 / s2 ~ -1.818e-13 = -1818 / 10^16
from fractions import Fraction
t_test = Fraction(-1818, 10**16)
def eval_frac(p, x):
    v = Fraction(0)
    for c in reversed(p): v = v * x + c
    return v
def deriv_int(p): return [r*p[r] for r in range(1, len(p))]
w_exact = eval_frac(P2[i0+1], t_test)*eval_frac(deriv_int(P2[i0]), t_test) - eval_frac(deriv_int(P2[i0+1]), t_test)*eval_frac(P2[i0], t_test)
print(f"  Exact W(P2[i0+1], P2[i0])(-1818/10^16) = {float(w_exact):.4e} (< 0? {w_exact < 0})")

# Now sweep exact integer b3 for P3[i0+1] = D_{b3} P2[i0+1] + P2[i0] to get COMPLEX ROOTS!
print("  Sweeping exact integer b3 for COMPLEX ROOTS of P3[i0+1] (deg 41):")
for b3 in [2, 10**10, 10**13, 10**14, 181330000000000, 5*10**14, 2*10**15, 10**16, 10**18]:
    P3_i = add_int(apply_D_int(b3, P2[i0+1]), P2[i0])
    d, nr, sqf = count_real_roots_exact(P3_i)
    print(f"    b3={b3:18d}: deg={d}, real_roots={nr}, sqf={sqf}" + ("  *** COMPLEX ROOTS CERTIFIED BY EXACT STURM!! ***" if nr < d else ""))
    if nr < d:
        break

print("\n=== 2. NOW TEST DOUBLE-OPERATOR prod_{j=1}^k D_{b_j}(D_{b_j} + y)(1), ALL b_j >= 2, a_j=1, w_j=1! ===")
# Wait: for DOUBLE-operator, D_2(D_2 + y) commutes with itself!
#   prod_{j=1}^K D_2(D_2 + y)(1) = sum_{i=0}^K C(K, i) D_2^{2K - i}(1) y^i !!
# Look at the slope ratio between i and i+1: STILL EXACTLY 2!!
# Wait: check what D_b(D_b + y) does at t = -eps where s_2*eps = alpha << 1 and 1 - b*eps = Phi!!
# Let's compute in 120-dps mpmath at K = 200, i0 = 185 (or K = 1000, i0 = 960? Wait: for DOUBLE, deg is 2K - i0!
# At K = 150, i0 = 135: K - i0 = 15, w_eff = 135/16 = 8.44 > 5.8, deg = 165 -> super fast!)
# Let's check K = 200, i0 = 180 (w_eff = 9.0, deg = 220) and K = 350, i0 = 325 (w_eff = 12.5, deg = 375)!
mp.mp.dps = 160
def apply_D_mp(b, p):
    res = [mp.mpf(0)]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+b)
    return res
def add_mp(a, b):
    res = [mp.mpf(0)]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res
def eval_mp(p, tv):
    v = mp.mpf(0)
    for c in reversed(p): v = v * tv + c
    return v
def deriv_mp(p): return [r*p[r] for r in range(1, len(p))]
def eval_W_mp(f, g, tv):
    return eval_mp(f, tv)*eval_mp(deriv_mp(g), tv) - eval_mp(deriv_mp(f), tv)*eval_mp(g, tv)

# Wait: for DOUBLE-operator at K = 400, i0 = 375 (K-i0 = 25, w_eff = 14.4, max deg = 428):
# Can we first derive the EXACT asymptotic M_i^{(w)} for DOUBLE-operator D_b(D_b + y) at t = -eps?!
# Let p_m(-eps) = 1 - 2^{2-m} alpha, 1 - b*eps = Phi!
# First u_m = D_b p_m + w p_{m-1} = (Phi + w) - (Phi + 1 + 2w) 2^{2-m} alpha!
# Second q_m = D_b u_m!
# At t = -eps: q_m(-eps) = Phi * u_m(-eps) - eps * u_m'(-eps)!
# We ALREADY computed eps * u_m'(-eps) = (1 - Phi) + (2 Phi + 2 + w) 2^{2-m} alpha!!
# SO q_m(-eps) = [Phi(Phi+w) - (1-Phi)] - [Phi(Phi+1+2w) + (2Phi+2+w)] 2^{2-m} alpha
#              = (Phi^2 + (w+1)Phi - 1) - (Phi^2 + (2w+3)Phi + (w+2)) 2^{2-m} alpha !!
# Look at q_1^2 - q_0 q_2 = A_2 * B_2 * alpha where A_2 = Phi^2 + (w+1)Phi - 1, B_2 = Phi^2 + (2w+3)Phi + (w+2)!
# Test numerically at K = 300, i0 = 275 (K-i0 = 25, w_eff = 10.6, deg = 328)!
K = 300; i0 = 275
D2m = [[mp.mpf(1)]]
for m in range(1, 2*K - i0 + 5):
    D2m.append(apply_D_mp(mp.mpf(2), D2m[-1]))
norm = mp.mpf(math.comb(K, i0))
P = {r: [(mp.mpf(math.comb(K, r)) / norm) * c for c in D2m[2*K - r]] for r in range(i0 - 3, i0 + 3)}
s2 = eval_mp(deriv_mp(P[i0+1]), mp.mpf(0)) / eval_mp(P[i0+1], mp.mpf(0))
w_eff = i0 / (K - i0 + 1)
print(f"  DOUBLE K={K}, i0={i0}, w_eff={w_eff:.2f}, s2={float(s2):.3e}")

best_M_dbl = (mp.mpf('1e100'), None)
for alpha in [0.05, 0.2, 0.5]:
    eps = mp.mpf(str(alpha)) / s2
    for Phi_f in [-30.0, -20.0, -15.0, -10.0, -8.0, -6.0, -4.0, -2.0, -1.0, -0.5]:
        Phi = mp.mpf(str(Phi_f))
        b1 = (1 - Phi) / eps
        P1 = {r: apply_D_mp(b1, add_mp(apply_D_mp(b1, P[r]), P[r-1])) for r in range(i0 - 2, i0 + 3)}
        tv = -eps
        v0, v1, v2 = eval_mp(P1[i0-1], tv), eval_mp(P1[i0], tv), eval_mp(P1[i0+1], tv)
        w21 = eval_W_mp(P1[i0+1], P1[i0], tv)
        M_val = (-tv)*(1+tv)*w21 + v1*v1 - v0*v2
        rel_M = M_val / (abs((-tv)*(1+tv)*w21) + v1*v1 + abs(v0*v2))
        if rel_M < best_M_dbl[0]:
            best_M_dbl = (rel_M, (alpha, Phi_f, b1, tv, P1))
print(f"  DOUBLE worst rel_M = {float(best_M_dbl[0]):+.6f} at alpha={best_M_dbl[1][0]}, Phi={best_M_dbl[1][1]}")
if best_M_dbl[0] < 0:
    print("  *** BOOM! M_i < 0 ALSO ACHIEVED FOR DOUBLE-OPERATOR prod D_{b_j}(D_{b_j}+y)(1) WITH ALL b_j>=2, a_j=1, w_j=1!! ***")
    _, (alpha, Phi_f, b1, t_neg, P1) = best_M_dbl
    v0, v1 = eval_mp(P1[i0-1], t_neg), eval_mp(P1[i0], t_neg)
    d1 = eval_mp(deriv_mp(P1[i0]), t_neg)
    b2_star = -(v1 + t_neg*(1+t_neg)*d1 + v0) / (t_neg * v1)
    print(f"  DOUBLE targeted b2* = {float(b2_star):.4e} (>=2? {b2_star >= 2})")
    if b2_star >= 2:
        P2 = {r: apply_D_mp(b2_star, add_mp(apply_D_mp(b2_star, P1[r]), P1[r-1])) for r in [i0-1, i0, i0+1]}
        min_w = mp.mpf('1e100')
        for mult in [0.9, 0.95, 0.99, 1.0, 1.01, 1.05, 1.1]:
            ww = eval_W_mp(P2[i0+1], P2[i0], t_neg * mp.mpf(str(mult)))
            if ww < min_w: min_w = ww
        print(f"  DOUBLE min W(P2[i0+1], P2[i0]) near t_neg = {float(min_w):.4e}")
