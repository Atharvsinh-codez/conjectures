"""Dual strategic diagnostic:
1. Verify at k = 35, i = 17 (where verifier found global-R failure at x = +7e-17 > 0) whether
   A_{35,17}(x)^2 > A_{35,16}(x) A_{35,18}(x) STILL HOLDS on ALL x < 0!!
2. Inspect Rossi-style alternating series dominance for psi_{k,i}(s) on s < 0 AND s > 0,
   and pole/principal-part dominance near s = 0, at k = 10, 20, 30!!
   Specifically:
   - How many of the d = 2k - i - 1 roots of psi_{k,i}(s) are in s < 0 vs s > 0?
   - On s < 0: at the natural sample points s_n where term n dominates sum_{m>=1} (-1)^{m-1} P_{k,i}(m) e^{m s},
     for which n = 1, 2, ..., n_max does ONE term dominate ALL other terms combined?!
   - On s > 0: at the natural sample points s_m where term m >= k+1 dominates sum_{m>=k+1} (-1)^m |P_{k,i}(-m)| e^{-m s},
     for which m >= k+1 does ONE term dominate ALL other terms combined?!
   - How many roots remain in the "bulk" window around s = 0?!
"""
import math
import numpy as np
import mpmath as mp

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_S(L, p):
    d = len(p) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(p):
        res[m] += c * (m + 1)
        res[m + 1] += c * (L - m)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def build_A(kmax):
    A = [[None] * (k + 1) for k in range(kmax + 1)]
    A[0][0] = [1]
    for k in range(1, kmax + 1):
        A[k][0] = apply_S(3*k - 2, apply_S(3*k - 3, A[k-1][0]))
        for i in range(1, k):
            L = 3*k - 2 - i
            A[k][i] = apply_S(L, poly_add(apply_S(L - 1, A[k-1][i]), A[k-1][i-1]))
        A[k][k] = apply_S(2*k - 2, A[k-1][k-1])
    return A

A = build_A(35)

# 1. Check A_{35,17}^2 - A_{35,16}*A_{35,18} on x < 0!
print("=== 1. CHECK A_{35,17}(x)^2 - A_{35,16}(x)*A_{35,18}(x) ON x < 0 vs x > 0 ===")
mp.mp.dps = 60
def eval_mp(poly, x):
    acc = mp.mpf(0)
    for c in reversed(poly): acc = acc * x + c
    return acc

min_neg = (1e300, None)
for s in np.linspace(-50.0, 20.0, 350):
    x = -mp.exp(mp.mpf(float(s)))
    v0 = eval_mp(A[35][16], x); v1 = eval_mp(A[35][17], x); v2 = eval_mp(A[35][18], x)
    rel = float((v1*v1 - v0*v2) / (v1*v1 + abs(v0*v2) + mp.mpf('1e-100')))
    if rel < min_neg[0]: min_neg = (rel, float(s))
print(f"  x < 0 (s in [-50, 20]): min rel (A17^2 - A16*A18) = {min_neg[0]:+.5e} at s={min_neg[1]:.2f}")
# Check verifier's positive x = 7e-17
x_pos = mp.mpf('7e-17')
v0 = eval_mp(A[35][16], x_pos); v1 = eval_mp(A[35][17], x_pos); v2 = eval_mp(A[35][18], x_pos)
print(f"  x = +7e-17 (verifier witness): A17^2 - A16*A18 = {mp.nstr(v1*v1 - v0*v2, 6)}")

# 2. Check Rossi-style alternating series dominance on s < 0 AND s > 0!
# Recall: sum_{n>=1} P_{k,i}(n) x^n = x A_{k,i}(x) / (1-x)^{3k-i+1},
# and     sum_{r>=0} |P_{k,i}(-(k+1+r))| u^r = A_{k,i}^{rev}(u) / (1-u)^{3k-i+1}!
def get_P_pos(A_poly, D_pow, nmax):
    """Exact integers P(1), P(2), ..., P(nmax) from x*A(x)/(1-x)^{D_pow}."""
    # P(n) = sum_{m=0}^{deg} a_m * binom(n - 1 - m + D_pow - 1, D_pow - 1)
    # Faster: compute via D_pow repeated prefix sums of A_poly padded with zeros!
    arr = [0] * nmax
    for m, c in enumerate(A_poly):
        if m < nmax: arr[m] = c
    for _ in range(D_pow):
        for j in range(1, nmax):
            arr[j] += arr[j - 1]
    return arr  # arr[n-1] = P(n)

def check_one_term_dominance(seq, label, n_check):
    """For alternating series sum_{j>=0} (-1)^j seq[j] e^{j * u} (u < 0),
    at optimal u_j maximizing term j (approx u_j = -0.5*(log(seq[j+1]/seq[j]) + log(seq[j]/seq[j-1]))),
    check if seq[j] e^{j u_j} > sum_{r != j} seq[r] e^{r u_j}!
    """
    dom_count = 0
    max_j_dom = -1
    ratios = []
    for j in range(0, min(n_check, len(seq) - 2)):
        if j == 0:
            # pick u very negative -> always dominates!
            dom_count += 1
            max_j_dom = 0
            continue
        r_prev = math.log(seq[j]) - math.log(seq[j - 1])
        r_next = math.log(seq[j + 1]) - math.log(seq[j])
        if r_prev <= r_next:
            # Not log-concave here!
            ratios.append((j, 999.0))
            continue
        u = -0.5 * (r_prev + r_next)
        if u >= 0:
            continue
        # Evaluate sum_{r != j} seq[r] e^{r u} / (seq[j] e^{j u})
        log_tj = math.log(seq[j]) + j * u
        err_sum = 0.0
        for r in range(max(0, j - 15), min(len(seq), j + 16)):
            if r == j: continue
            err_sum += math.exp(math.log(seq[r]) + r * u - log_tj)
        if err_sum < 1.0:
            dom_count += 1
            max_j_dom = j
        ratios.append((j, err_sum))
    return dom_count, max_j_dom, ratios

print("\n=== 2. ONE-TERM ALTERNATING SERIES DOMINANCE ON s < 0 AND s > 0 ===")
for k, i in [(10, 5), (20, 10), (30, 15), (30, 5), (30, 25)]:
    d = 2*k - i - 1
    D_pow = 3*k - i + 1
    # s < 0 series: P(n), n = 1, 2, ...
    seq_neg = get_P_pos(A[k][i], D_pow, 80)
    c_neg, mj_neg, r_neg = check_one_term_dominance(seq_neg, "s<0", 40)
    # s > 0 series: reciprocal A^{rev}(u) / (1-u)^{D_pow}
    seq_pos = get_P_pos(list(reversed(A[k][i])), D_pow, 80)
    c_pos, mj_pos, r_pos = check_one_term_dominance(seq_pos, "s>0", 40)
    # Also count exact roots in s < 0 (i.e. x in (-1, 0)) vs s > 0 (i.e. x in (-inf, -1))!
    # Wait: use sign changes / rough count via log-concavity
    print(f"  k={k:2d}, i={i:2d} (deg={d:2d}): s<0 1-term dom={c_neg:2d} (up to j={mj_neg:2d}) | "
          f"s>0 1-term dom={c_pos:2d} (up to j={mj_pos:2d}) | total_1term={c_neg+c_pos:2d}/{d+1:2d} "
          f"(bulk gap = {d + 1 - c_neg - c_pos:2d})")
