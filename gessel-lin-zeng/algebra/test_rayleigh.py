import math
import numpy as np
from itertools import combinations

def poly_trim(c):
    while len(c) > 1 and c[-1] == 0: c.pop()
    return c

def apply_T(a, b, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + a) * cm
        out[m + 1] += (m + b) * cm
    return poly_trim(out)

def mul_poly(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
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

def eval_grid(p, grid):
    acc = np.zeros_like(grid, dtype=float)
    for c in reversed(p):
        acc = acc * grid + float(c)
    return acc

t_int = np.linspace(-1.0, 0.0, 1001)
t_R = np.concatenate([np.linspace(-5.0, -1.0, 100), t_int, np.linspace(0.0, 5.0, 100)])

def compute_all_pS(k, ops, glz_double=False):
    """ops[j-1] = (a_j, b_j) for j=1..k.
    If glz_double=False: p_S = (prod_{j not in S} T_j)(1).
    If glz_double=True:  p_S = (prod_{j=1}^k T_j^{2 if j not in S else 1})(1)."""
    pS = {frozenset(): [1]}
    for j in range(1, k + 1):
        aj, bj = ops[j - 1]
        nxt = {}
        for S, poly in pS.items():
            # j not in S: apply T_j (or T_j^2)
            p_not = apply_T(aj, bj, poly)
            if glz_double:
                p_not = apply_T(aj, bj, p_not)
            nxt[S] = p_not
            # j in S: keep poly (or apply T_j once)
            p_in = apply_T(aj, bj, poly) if glz_double else poly
            nxt[S | {j}] = p_in
        pS = nxt
    return pS

# === TEST 1: Check Rayleigh property p_{I U {a}} p_{I U {b}} - p_I p_{I U {a,b}} on (-1, 0) ===
print("=== TEST 1: Check pairwise Rayleigh difference R = p_{I U a} p_{I U b} - p_I p_{I U {a,b}} on (-1, 0) ===")
for name, k, ops, dbl in [
    ("Inner GLZ k=5", 5, [(1, j+1) for j in range(1, 6)], False),
    ("Full GLZ k=5",  5, [(1, j+1) for j in range(1, 6)], True),
    ("Wild inner k=5", 5, [(1,2),(1,100),(1,2),(1,100),(1,2)], False),
    ("Wild full k=4", 4, [(1,2),(1,100),(1,2),(1,100)], True),
    ("Gen (a,b) k=4", 4, [(1,5),(10,12),(2,50),(3,4)], False),
]:
    pS = compute_all_pS(k, ops, dbl)
    # Test BOTH R = p_{I U a} p_{I U b} - p_I p_{I U ab} AND C_{I, {a,b}} = 2 p_{I U a} p_{I U b} - p_I p_{I U ab}
    fail_R_int = 0
    fail_C2_int = 0
    tot = 0
    worst_R = 999.0
    worst_C2 = 999.0
    for a, b in combinations(range(1, k + 1), 2):
        rem = [x for x in range(1, k + 1) if x not in (a, b)]
        for r in range(len(rem) + 1):
            for I_tup in combinations(rem, r):
                I = frozenset(I_tup)
                tot += 1
                prod1 = mul_poly(pS[I | {a}], pS[I | {b}])
                prod2 = mul_poly(pS[I], pS[I | {a, b}])
                R = sub_poly(prod1, prod2)
                C2 = sub_poly(add_poly(prod1, prod1), prod2)
                vR = eval_grid(R, t_int)
                vC2 = eval_grid(C2, t_int)
                if np.min(vR) < -1e-9:
                    fail_R_int += 1
                if np.min(vC2) < -1e-9:
                    fail_C2_int += 1
                worst_R = min(worst_R, float(np.min(vR)))
                worst_C2 = min(worst_C2, float(np.min(vC2)))
    print(f"  {name:18s} ({tot:3d} pairs): fail(R>=0 on (-1,0))={fail_R_int:3d} (min={worst_R:+.3e}) | fail(C2>=0 on (-1,0))={fail_C2_int:3d} (min={worst_C2:+.3e})")

# === TEST 2: Check ALL weight-monomial coefficients C_{I, J}(t) in q_i^2 - q_{i-1} q_{i+1} for |J| = 2, 4, 6! ===
print("\n=== TEST 2: Check ALL monomial coefficients C_{I, J}(t) of w^m in q_{k,i}^2 - q_{k,i-1}q_{k,i+1} on (-1, 0) ===")
for name, k, ops, dbl in [
    ("Inner GLZ k=6", 6, [(1, j+1) for j in range(1, 7)], False),
    ("Full GLZ k=6",  6, [(1, j+1) for j in range(1, 7)], True),
    ("Wild inner k=6", 6, [(1,2),(1,100),(1,2),(1,100),(1,2),(1,100)], False),
    ("Full GLZ k=7",  7, [(1, j+1) for j in range(1, 8)], True),
]:
    pS = compute_all_pS(k, ops, dbl)
    # For every even |J| = 2m >= 2 and disjoint I:
    for m in [1, 2, 3]:
        if 2*m > k: continue
        tot = 0
        fails_int = 0
        fails_R = 0
        for J_tup in combinations(range(1, k + 1), 2*m):
            J = set(J_tup)
            rem = [x for x in range(1, k + 1) if x not in J]
            for r in range(len(rem) + 1):
                for I_tup in combinations(rem, r):
                    I = frozenset(I_tup)
                    tot += 1
                    # sum_{U cup V = J, |U|=m} p_{I U U} p_{I U V} - sum_{|U|=m-1} p_{I U U} p_{I U V}
                    pos_sum = [0]
                    for U_tup in combinations(J_tup, m):
                        U = frozenset(U_tup)
                        V = frozenset(J - set(U_tup))
                        pos_sum = add_poly(pos_sum, mul_poly(pS[I | U], pS[I | V]))
                    neg_sum = [0]
                    for U_tup in combinations(J_tup, m - 1):
                        U = frozenset(U_tup)
                        V = frozenset(J - set(U_tup))
                        neg_sum = add_poly(neg_sum, mul_poly(pS[I | U], pS[I | V]))
                    diff = sub_poly(pos_sum, neg_sum)
                    if np.min(eval_grid(diff, t_int)) < -1e-7:
                        fails_int += 1
                    if np.min(eval_grid(diff, t_R)) < -1e-7:
                        fails_R += 1
        print(f"  {name:16s} |J|={2*m}: fails on (-1,0) = {fails_int:4d}/{tot:4d} | fails on R = {fails_R:4d}/{tot:4d}")
