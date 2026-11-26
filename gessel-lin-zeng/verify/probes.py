"""
Adversarial exact-arithmetic probes for Gessel-Lin-Zeng Conjecture 15:
1. Check 4th cross-pair D_{k+1} htilde_{k-1,i} << htilde_{k-1,i-2} through k=20.
2. Stress-test the surviving bivariate stability candidate B_k(x, y) = sum binom(k,i) A_{k,i}(x) y^i
   and also sum binom(k,i) h_{k,i}(t) y^i along real specializations AND positive-slope lines (x0 + a*s, y0 + b*s).
"""

import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import testbed as tb


def real_root_count_squarefree_check(poly):
    poly = tb.poly_trim(poly)
    d = len(poly) - 1
    if d <= 1:
        return True, d, d
    ch = tb.sturm_chain(poly, tb.poly_deriv(poly))
    if len(ch[-1]) > 1:
        # Divide out gcd to check squarefree part? Actually just check roots of squarefree part
        return None, d, -1
    v_neg = tb.sign_variations(tb.chain_signs_at(ch, "-inf"))
    v_pos = tb.sign_variations(tb.chain_signs_at(ch, "+inf"))
    n_real = v_neg - v_pos
    return (n_real == d), d, n_real


def main():
    MAX_K = 20
    h = tb.compute_h_table(MAX_K)
    ht = {(k, i): tb.poly_div_exact_linear(h[(k, i)], -1)
          for k in range(1, MAX_K + 1) for i in range(k + 1)}

    p4_pass, p4_fail = 0, []
    for k in range(3, MAX_K + 1):
        for i in range(2, k):
            F = tb.poly_div_exact_linear(tb.op_D(k, h[(k - 1, i)]), -1)
            R = ht[(k - 1, i - 2)]
            ok, msg = tb.verify_strict_interlacing(F, R, -1, 0)
            if ok:
                p4_pass += 1
            else:
                p4_fail.append((k, i, msg))
    print(f"4th cross-pair D_{{k+1}} htilde_{{k-1,i}} << htilde_{{k-1,i-2}} (k=3..{MAX_K}): "
          f"PASS={p4_pass}, FAIL={len(p4_fail)} {p4_fail[:5]}")

    A = tb.compute_A_table(10)
    print("\n--- Testing B_k(x, y) = sum binom(k, i) A_{k,i}(x) y^i for real stability ---")
    for k in range(2, 11):
        bad = None
        # Test univariate specializations in y at integer x0
        for x0 in [1, 2, 3, 5, 10, 0, -1, -2, -3, -5, -10]:
            poly_y = [math.comb(k, i) * tb.poly_eval_int(A[(k, i)], x0) for i in range(k + 1)]
            ok, d, nr = real_root_count_squarefree_check(poly_y)
            if ok is False:
                bad = (f"in y at x={x0}", f"deg={d}, real={nr}")
                break
        if not bad:
            # Test univariate specializations in x at integer y0
            for y0 in [-1, -2, -3, -5, -10, -20, -100, 1, 2, 5, 10]:
                poly_x = [0]
                for i in range(k + 1):
                    coeff = math.comb(k, i) * (y0 ** i)
                    poly_x = tb.poly_add(poly_x, [coeff * c for c in A[(k, i)]])
                ok, d, nr = real_root_count_squarefree_check(poly_x)
                if ok is False:
                    bad = (f"in x at y={y0}", f"deg={d}, real={nr}")
                    break
        if not bad:
            # Test diagonal positive-slope lines x = s, y = a*s + b (a > 0)
            for a, b in [(1, 0), (1, 1), (1, -1), (2, -3), (1, -5), (3, -2), (1, -10), (5, -1)]:
                # Evaluate B_k(s, a*s + b) as polynomial in s
                poly_s = [0]
                lin_pow = [1]  # (b + a*s)^i
                lin = [b, a]
                for i in range(k + 1):
                    term = tb.poly_mul([math.comb(k, i)], tb.poly_mul(A[(k, i)], lin_pow))
                    poly_s = tb.poly_add(poly_s, term)
                    lin_pow = tb.poly_mul(lin_pow, lin)
                ok, d, nr = real_root_count_squarefree_check(poly_s)
                if ok is False:
                    bad = (f"line y={a}s+{b}", f"deg={d}, real={nr}")
                    break
        status = f"FAIL {bad}" if bad else "survived all tested lines"
        print(f"  B_{k}(x, y): {status}")

    print("\n--- Testing H_k^binom(t, y) = sum binom(k, i) h_{k,i}(t) y^i for real stability ---")
    for k in range(2, 9):
        bad = None
        for t0 in [1, 2, 5, -2, -3, -5]:
            poly_y = [math.comb(k, i) * tb.poly_eval_int(h[(k, i)], t0) for i in range(k + 1)]
            ok, d, nr = real_root_count_squarefree_check(poly_y)
            if ok is False:
                bad = (f"in y at t={t0}", f"deg={d}, real={nr}")
                break
        if not bad:
            for a, b in [(1, 0), (1, -1), (1, -5), (2, -3)]:
                poly_s = [0]
                lin_pow = [1]
                lin = [b, a]
                for i in range(k + 1):
                    term = tb.poly_mul([math.comb(k, i)], tb.poly_mul(h[(k, i)], lin_pow))
                    poly_s = tb.poly_add(poly_s, term)
                    lin_pow = tb.poly_mul(lin_pow, lin)
                ok, d, nr = real_root_count_squarefree_check(poly_s)
                if ok is False:
                    bad = (f"line y={a}s+{b}", f"deg={d}, real={nr}")
                    break
        status = f"FAIL {bad}" if bad else "survived all tested lines"
        print(f"  H_{k}^binom(t, y): {status}")

    print("\n--- Testing rational specializations in (-2, 1) for B_k(x,y) and H_k^binom(t,y) ---")
    def check_poly_y_rational(polys, weights, p, q):
        maxd = max(len(P) - 1 for P in polys)
        poly_y = []
        for i, P in enumerate(polys):
            val = sum(c * (p ** r) * (q ** (maxd - r)) for r, c in enumerate(P))
            poly_y.append(weights[i] * val)
        return real_root_count_squarefree_check(poly_y)

    for name, tbl in [("B_k(x,y) [binom]", A), ("H_k^binom(t,y) [binom]", h)]:
        print(f"  [{name}]")
        for k in range(2, 9):
            polys = [tbl[(k, i)] for i in range(k + 1)]
            w = [math.comb(k, i) for i in range(k + 1)]
            bad = None
            for q in [2, 3, 4, 5, 7, 10, 13, 20]:
                for p in range(-2 * q, q + 1):
                    if p == 0 or math.gcd(abs(p), q) != 1:
                        continue
                    ok, d, nr = check_poly_y_rational(polys, w, p, q)
                    if ok is False:
                        bad = (f"{p}/{q}", d, nr)
                        break
                if bad:
                    break
            print(f"    k={k}: {'FAIL at ' + str(bad) if bad else 'all rational specializations real-rooted in y'}")


if __name__ == "__main__":
    main()
