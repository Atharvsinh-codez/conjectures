"""
Deep stress test (k up to 50 for x-coefficients, plus exact Sturm on R at k=19):
1. Check whether M_x^2 - X_x Y_x or W(A_{k,i-1}, A_{k,i}) EVER get negative x-coefficients up to k=50!
2. Check whether A_{19, 8}^2 - A_{19, 7} A_{19, 9} (whose [x^1] coeff is negative!) STILL has
   strictly positive values on all of R (exact Sturm), or if Turan on R also breaks!
3. Also check x-coefficients of M_x^2 - X_x Y_x and W(A_{k,i-1}, A_{k,i}) AND
   check [x^1] ratio formula for A_{k,i}.
4. Check Sturm on R for Turan A_{k,i}^2 - A_{k,i-1} A_{k,i+1} at k=14..20 (especially i = k//2).
"""

import sys
import time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import testbed as tb


def main():
    MAX_K = 50
    t0 = time.perf_counter()
    A = tb.compute_A_table(MAX_K)
    print(f"Computed exact A_{{k,i}} up to k={MAX_K} in {time.perf_counter()-t0:.2f}s.")

    # 1. Check M_x^2 - X_x Y_x x-coefficients up to k=50
    t1 = time.perf_counter()
    m2_neg = []
    for k in range(3, MAX_K + 1):
        for i in range(2, k):
            L = 3 * k - i - 2
            X_x = A[(k - 1, i - 2)]
            M_x = tb.op_S(L, A[(k - 1, i - 1)])
            Y_x = tb.op_S(L, tb.op_S(L - 1, A[(k - 1, i)]))
            P = tb.poly_sub(tb.poly_mul(M_x, M_x), tb.poly_mul(X_x, Y_x))
            neg_idx = [r for r, c in enumerate(P) if c < 0]
            if neg_idx:
                m2_neg.append((k, i, neg_idx[:5], P[neg_idx[0]]))
                if len(m2_neg) <= 5:
                    print(f"  [!] M_x^2 - X_x Y_x NEGATIVE x-coeff at k={k}, i={i}: idx={neg_idx[:5]}, val={P[neg_idx[0]]}")
    print(f"1. M_x^2 - X_x Y_x x-coeff scan (k=3..{MAX_K}, {time.perf_counter()-t1:.2f}s): total fails = {len(m2_neg)}")

    # 2. Check W(A_{k,i-1}, A_{k,i}) x-coefficients up to k=50
    t2 = time.perf_counter()
    w_neg = []
    for k in range(1, MAX_K + 1):
        for i in range(1, k + 1):
            P_big = A[(k, i - 1)]
            P_sml = A[(k, i)]
            W = tb.poly_sub(
                tb.poly_mul(tb.poly_deriv(P_big), P_sml),
                tb.poly_mul(P_big, tb.poly_deriv(P_sml)),
            )
            neg_idx = [r for r, c in enumerate(W) if c < 0]
            if neg_idx:
                w_neg.append((k, i, neg_idx[:5], W[neg_idx[0]]))
                if len(w_neg) <= 5:
                    print(f"  [!] W(A_{{k,i-1}}, A_{{k,i}}) NEGATIVE x-coeff at k={k}, i={i}: idx={neg_idx[:5]}, val={W[neg_idx[0]]}")
    print(f"2. W(A_{{k,i-1}}, A_{{k,i}}) x-coeff scan (k=1..{MAX_K}, {time.perf_counter()-t2:.2f}s): total fails = {len(w_neg)}")

    # 3. Check ALL 2x2 minors A_{k,i}A_{k,j} - A_{k,i-1}A_{k,j+1} first failure details
    first_offdiag = None
    first_diag = None
    for k in range(2, 30):
        for i in range(1, k):
            for j in range(i, k):
                P = tb.poly_sub(
                    tb.poly_mul(A[(k, i)], A[(k, j)]),
                    tb.poly_mul(A[(k, i - 1)], A[(k, j + 1)]),
                )
                neg_idx = [r for r, c in enumerate(P) if c < 0]
                if neg_idx:
                    if i == j and first_diag is None:
                        first_diag = (k, i, neg_idx, P[neg_idx[0]])
                    if i < j and first_offdiag is None:
                        first_offdiag = (k, i, j, neg_idx, P[neg_idx[0]])
        if first_diag and first_offdiag:
            break
    print(f"3. First negative x-coeff in 2x2 minors: diag={first_diag}, offdiag={first_offdiag}")

    # 4. Exact Sturm check on R for A_{k,i}^2 - A_{k,i-1} A_{k,i+1} at k=19, i=8, 9 (and k=18, i=9)
    # Wait: check sign at x = 0, x = 1/1000, x = 1/100, etc. FIRST (fast) before full Sturm!
    print("4. Checking Turan A_{k,i}^2 - A_{k,i-1} A_{k,i+1} on R at k=19, i=8,9 and k=25, i=12:")
    for k, i in [(19, 8), (19, 9), (20, 10), (25, 12), (30, 15)]:
        P = tb.poly_sub(tb.poly_mul(A[(k, i)], A[(k, i)]), tb.poly_mul(A[(k, i - 1)], A[(k, i + 1)]))
        neg_idx = [r for r, c in enumerate(P) if c < 0]
        # Evaluate at x = 1/q for small positive rational x to see if it EVER goes negative on (0, inf)!
        min_pos = None
        for q in [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000]:
            # P(1/q) * q^deg
            d = len(P) - 1
            val = sum(c * (q ** (d - r)) for r, c in enumerate(P))
            if val <= 0:
                min_pos = (f"1/{q}", val)
                break
        print(f"   ({k},{i}): deg={len(P)-1}, P[0]={P[0]}, P[1]={P[1]}, neg_idx={neg_idx[:5]}, neg_on_tested_pos_x={min_pos}")

    # Exact Sturm on (19, 8) and (19, 9)
    for k, i in [(19, 8), (19, 9)]:
        t3 = time.perf_counter()
        P = tb.poly_sub(tb.poly_mul(A[(k, i)], A[(k, i)]), tb.poly_mul(A[(k, i - 1)], A[(k, i + 1)]))
        ch = tb.sturm_chain(P, tb.poly_deriv(P))
        sqfree = (len(ch[-1]) == 1 and ch[-1][0] != 0)
        v_neg = tb.sign_variations(tb.chain_signs_at(ch, "-inf"))
        v_zero = tb.sign_variations(tb.chain_signs_at(ch, 0))
        v_pos = tb.sign_variations(tb.chain_signs_at(ch, "+inf"))
        print(f"   Exact Sturm for Turan ({k},{i}) ({time.perf_counter()-t3:.2f}s): sqfree={sqfree}, neg_roots={v_neg-v_zero}, pos_roots={v_zero-v_pos}")


if __name__ == "__main__":
    main()
