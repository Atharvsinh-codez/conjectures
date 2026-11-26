"""
Independent exact-arithmetic verification of #13 & #14 for Gessel-Lin-Zeng Conjecture 15:
1. Re-verify #13's exact disproof of B_{10}(x, y) = sum binom(10, i) A_{10, i}(x) y^i real stability:
   (a) Exact Fraction evaluation of Rayleigh determinant B_x B_y - B B_{xy} at (x, y) = (0, -3/4),
       plus exact check for k = 1..10.
   (b) Exact integer Sturm count of B_{10}(s, -2531/1000 + (130151/100) s) (cleared denominators).
   (c) Plain sum A_{k,i}(x) y^i at k=3, x = -1/100.
2. Stress-test x-coefficientwise nonnegativity AND strict positivity on R (via exact Sturm) up to high k:
   (a) M_x^2 - X_x Y_x where X_x = A_{k-1, i-2}, M_x = S_{3k-i-2} A_{k-1, i-1},
       Y_x = S_{3k-i-2} S_{3k-i-3} A_{k-1, i} (exact x-transform of M^2 - X Y).
   (b) All 2x2 minors A_{k, i} A_{k, j} - A_{k, i-1} A_{k, j+1} for 1 <= i <= j <= k-1.
   (c) Wronskians W(A_{k, i-1}, A_{k, i}) = A'_{k, i-1} A_{k, i} - A_{k, i-1} A'_{k, i}.
   (d) Also check 3x3 minors / Turan A_{k,i}^2 - A_{k,i-1} A_{k,i+1} on R via Sturm!
"""

from __future__ import annotations
from fractions import Fraction
import math
import sys
import time

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import testbed as tb


def poly_eval_frac(p, x: Fraction) -> Fraction:
    acc = Fraction(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc


def verify_B10_disproof(A):
    print("=== 1. Independent re-verification of B_{10}(x, y) real-stability disproof ===")
    # 1(a): Rayleigh determinant B_x * B_y - B * B_xy at (x0, y0) = (0, -3/4)
    k = 10
    x0 = Fraction(0)
    y0 = Fraction(-3, 4)

    B = Fraction(0)
    Bx = Fraction(0)
    By = Fraction(0)
    Bxy = Fraction(0)
    for i in range(k + 1):
        w = math.comb(k, i)
        Ai = A[(k, i)]
        dAi = tb.poly_deriv(Ai)
        val_A = poly_eval_frac(Ai, x0)
        val_dA = poly_eval_frac(dAi, x0)
        yi = y0 ** i
        B += w * val_A * yi
        Bx += w * val_dA * yi
        if i >= 1:
            yim1 = y0 ** (i - 1)
            By += w * i * val_A * yim1
            Bxy += w * i * val_dA * yim1

    rayleigh = Bx * By - B * Bxy
    expected = Fraction(-2649516238079970755, 4294967296)
    print(f"  Rayleigh(B_10) at (0, -3/4) = {rayleigh}")
    print(f"  Matches #13 value exactly?    {rayleigh == expected} (sign < 0: {rayleigh < 0})")
    assert rayleigh == expected and rayleigh < 0

    # Also check Rayleigh at (0, -3/4) for k=1..9 just for context
    for kk in range(1, 10):
        b0 = bx0 = by0 = bxy0 = Fraction(0)
        for i in range(kk + 1):
            w = math.comb(kk, i)
            Ai = A[(kk, i)]
            dAi = tb.poly_deriv(Ai)
            va = poly_eval_frac(Ai, x0)
            vda = poly_eval_frac(dAi, x0)
            b0 += w * va * (y0 ** i)
            bx0 += w * vda * (y0 ** i)
            if i >= 1:
                by0 += w * i * va * (y0 ** (i - 1))
                bxy0 += w * i * vda * (y0 ** (i - 1))
        r_kk = bx0 * by0 - b0 * bxy0
        print(f"    k={kk}: Rayleigh(B_{kk}) at (0, -3/4) = {float(r_kk):.6e} ({'POS' if r_kk > 0 else 'NONPOS'})")

    # 1(b): Exact integer Sturm count on positive-slope line
    # y(s) = -2531/1000 + (130151/100) * s = (-2531 + 1301510 * s) / 1000
    # Multiply B_{10}(s, y(s)) by 1000^10 so every term is an exact integer polynomial in s!
    lin_num = [-2531, 1301510]  # 1000 * y(s)
    poly_s = [0]
    lin_pow = [1]  # (1000 * y(s))^i
    for i in range(k + 1):
        scale = math.comb(k, i) * (1000 ** (k - i))
        term = tb.poly_mul([scale], tb.poly_mul(A[(k, i)], lin_pow))
        poly_s = tb.poly_add(poly_s, term)
        lin_pow = tb.poly_mul(lin_pow, lin_num)
    poly_s = tb.poly_primitive_pos(poly_s)
    deg_s = len(poly_s) - 1
    ch = tb.sturm_chain(poly_s, tb.poly_deriv(poly_s))
    sqfree = (len(ch[-1]) == 1 and ch[-1][0] != 0)
    n_real = tb.sign_variations(tb.chain_signs_at(ch, "-inf")) - tb.sign_variations(tb.chain_signs_at(ch, "+inf"))
    print(f"  B_10(s, -2531/1000 + (130151/100)s): deg={deg_s}, squarefree={sqfree}, real_roots={n_real}")
    assert deg_s == 19 and sqfree and n_real == 17

    # 1(c): Plain sum A_{3,i}(x) y^i at x = -1/100
    k3 = 3
    # Multiply by 100^5 (max deg of A_{3,i} is 5)
    p3_y = []
    for i in range(k3 + 1):
        val = sum(c * ((-1) ** r) * (100 ** (5 - r)) for r, c in enumerate(A[(k3, i)]))
        p3_y.append(val)
    ch3 = tb.sturm_chain(p3_y, tb.poly_deriv(p3_y))
    n_real3 = tb.sign_variations(tb.chain_signs_at(ch3, "-inf")) - tb.sign_variations(tb.chain_signs_at(ch3, "+inf"))
    print(f"  Plain sum A_{{3,i}}(-1/100) y^i: deg={len(p3_y)-1}, sqfree={len(ch3[-1])==1}, real_roots={n_real3}")
    assert len(p3_y) - 1 == 3 and len(ch3[-1]) == 1 and n_real3 == 1


def verify_exact_t_to_x_M2_XY(A, h, max_k_check=10):
    """Verify algebraically & exactly that M^2 - X Y in t equals (1+t)^{2(2k-i-1)} (M_x^2 - X_x Y_x)(t/(1+t))."""
    for k in range(3, max_k_check + 1):
        for i in range(2, k):
            X_t = tb.poly_div_exact_linear(h[(k - 1, i - 2)], -1)
            M_t = tb.poly_div_exact_linear(tb.op_D(k, h[(k - 1, i - 1)]), -1)
            Y_t = tb.poly_div_exact_linear(tb.op_D(k, tb.op_D(k, h[(k - 1, i)])), -1)
            diff_t = tb.poly_sub(tb.poly_mul(M_t, M_t), tb.poly_mul(X_t, Y_t))

            L = 3 * k - i - 2
            X_x = A[(k - 1, i - 2)]
            M_x = tb.op_S(L, A[(k - 1, i - 1)])
            Y_x = tb.op_S(L, tb.op_S(L - 1, A[(k - 1, i)]))
            diff_x = tb.poly_sub(tb.poly_mul(M_x, M_x), tb.poly_mul(X_x, Y_x))

            m = 2 * (2 * k - i - 1)
            diff_t_from_x = tb.transform_A_to_h(diff_x, m)
            assert diff_t == diff_t_from_x, f"t-to-x mismatch for M^2-XY at ({k},{i})"
    print(f"  Verified exact identity M(t)^2 - X(t)Y(t) = (1+t)^{{2(2k-i-1)}} [M_x^2 - X_x Y_x](t/(1+t)) for k=3..{max_k_check}.")


def stress_test_positivity(A, max_k_coeff=24, max_k_sturm=13):
    print(f"\n=== 2. Stress-testing exact x-coefficient nonnegativity (k <= {max_k_coeff}) & R-positivity (k <= {max_k_sturm}) ===")

    # (a) M_x^2 - X_x Y_x
    t0 = time.perf_counter()
    m2xy_checked = 0
    m2xy_neg_coeff = []
    m2xy_zero_coeff = []
    m2xy_real_root = []
    m2xy_sturm_checked = 0
    for k in range(3, max_k_coeff + 1):
        for i in range(2, k):
            L = 3 * k - i - 2
            X_x = A[(k - 1, i - 2)]
            M_x = tb.op_S(L, A[(k - 1, i - 1)])
            Y_x = tb.op_S(L, tb.op_S(L - 1, A[(k - 1, i)]))
            P = tb.poly_sub(tb.poly_mul(M_x, M_x), tb.poly_mul(X_x, Y_x))
            m2xy_checked += 1
            neg_idx = [r for r, c in enumerate(P) if c < 0]
            zero_idx = [r for r, c in enumerate(P) if c == 0]
            if neg_idx:
                m2xy_neg_coeff.append((k, i, neg_idx[:5], P[neg_idx[0]]))
            if zero_idx:
                m2xy_zero_coeff.append((k, i, zero_idx[:5]))
            if k <= max_k_sturm:
                m2xy_sturm_checked += 1
                ch = tb.sturm_chain(P, tb.poly_deriv(P))
                sqfree = (len(ch[-1]) == 1 and ch[-1][0] != 0)
                n_real = tb.sign_variations(tb.chain_signs_at(ch, "-inf")) - tb.sign_variations(tb.chain_signs_at(ch, "+inf"))
                if (not sqfree) or n_real > 0 or P[0] <= 0 or P[-1] <= 0:
                    m2xy_real_root.append((k, i, sqfree, n_real))
    print(f"  (a) M_x^2 - X_x Y_x (k=3..{max_k_coeff}, {m2xy_checked} pairs, {time.perf_counter()-t0:.2f}s):")
    print(f"      Negative x-coeffs: {len(m2xy_neg_coeff)} {m2xy_neg_coeff[:5]}")
    print(f"      Zero x-coeffs:     {len(m2xy_zero_coeff)} {m2xy_zero_coeff[:5]}")
    print(f"      Exact Sturm > 0 on ALL of R (k=3..{max_k_sturm}, {m2xy_sturm_checked} pairs): FAIL={len(m2xy_real_root)} {m2xy_real_root[:5]}")

    # (b) All 2x2 minors A_{k,i} A_{k,j} - A_{k,i-1} A_{k,j+1} (1 <= i <= j <= k-1)
    t1 = time.perf_counter()
    minors_checked = 0
    minors_neg_coeff = []
    minors_zero_coeff = []
    turan_sturm_checked = 0
    turan_real_root = []
    minor_sturm_checked = 0
    minor_real_root = []
    for k in range(2, max_k_coeff + 1):
        for i in range(1, k):
            for j in range(i, k):
                P = tb.poly_sub(
                    tb.poly_mul(A[(k, i)], A[(k, j)]),
                    tb.poly_mul(A[(k, i - 1)], A[(k, j + 1)]),
                )
                minors_checked += 1
                neg_idx = [r for r, c in enumerate(P) if c < 0]
                zero_idx = [r for r, c in enumerate(P) if c == 0]
                if neg_idx:
                    minors_neg_coeff.append((k, i, j, neg_idx[:5], P[neg_idx[0]]))
                if zero_idx:
                    minors_zero_coeff.append((k, i, j, zero_idx[:5]))
                if i == j and k <= max_k_sturm:
                    turan_sturm_checked += 1
                    ch = tb.sturm_chain(P, tb.poly_deriv(P))
                    sqfree = (len(ch[-1]) == 1 and ch[-1][0] != 0)
                    n_real = tb.sign_variations(tb.chain_signs_at(ch, "-inf")) - tb.sign_variations(tb.chain_signs_at(ch, "+inf"))
                    if (not sqfree) or n_real > 0 or P[0] <= 0 or P[-1] <= 0:
                        turan_real_root.append((k, i, sqfree, n_real))
                if i < j and k <= 10:
                    minor_sturm_checked += 1
                    ch = tb.sturm_chain(P, tb.poly_deriv(P))
                    sqfree = (len(ch[-1]) == 1 and ch[-1][0] != 0)
                    n_real = tb.sign_variations(tb.chain_signs_at(ch, "-inf")) - tb.sign_variations(tb.chain_signs_at(ch, "+inf"))
                    if (not sqfree) or n_real > 0 or P[0] <= 0 or P[-1] <= 0:
                        minor_real_root.append((k, i, j, sqfree, n_real))
    print(f"  (b) All 2x2 minors A_{{k,i}} A_{{k,j}} - A_{{k,i-1}} A_{{k,j+1}} (k=2..{max_k_coeff}, {minors_checked} triples, {time.perf_counter()-t1:.2f}s):")
    print(f"      Negative x-coeffs: {len(minors_neg_coeff)} {minors_neg_coeff[:5]}")
    print(f"      Zero x-coeffs:     {len(minors_zero_coeff)} {minors_zero_coeff[:5]}")
    print(f"      Turan (i=j) > 0 on ALL of R (k=2..{max_k_sturm}, {turan_sturm_checked} pairs): FAIL={len(turan_real_root)} {turan_real_root[:5]}")
    print(f"      Off-diag (i<j) > 0 on ALL of R (k=2..10, {minor_sturm_checked} triples): FAIL={len(minor_real_root)} {minor_real_root[:5]}")

    # (c) Wronskians W(A_{k, i-1}, A_{k, i}) = A'_{k, i-1} A_{k, i} - A_{k, i-1} A'_{k, i}
    # Check sign convention: deg(A_{k, i-1}) = d, deg(A_{k, i}) = d-1.
    # Leading term of A'_{k, i-1} A_{k, i} - A_{k, i-1} A'_{k, i} is (d - (d-1)) lc(A_{k, i-1}) lc(A_{k, i}) x^{2d-2} > 0!
    t2 = time.perf_counter()
    wron_checked = 0
    wron_neg_coeff = []
    wron_zero_coeff = []
    wron_sturm_checked = 0
    wron_real_root = []
    for k in range(1, max_k_coeff + 1):
        for i in range(1, k + 1):
            P_big = A[(k, i - 1)]
            P_sml = A[(k, i)]
            W = tb.poly_sub(
                tb.poly_mul(tb.poly_deriv(P_big), P_sml),
                tb.poly_mul(P_big, tb.poly_deriv(P_sml)),
            )
            wron_checked += 1
            neg_idx = [r for r, c in enumerate(W) if c < 0]
            zero_idx = [r for r, c in enumerate(W) if c == 0]
            if neg_idx:
                wron_neg_coeff.append((k, i, neg_idx[:5], W[neg_idx[0]]))
            if zero_idx:
                wron_zero_coeff.append((k, i, zero_idx[:5]))
            if k <= max_k_sturm:
                wron_sturm_checked += 1
                ch = tb.sturm_chain(W, tb.poly_deriv(W))
                sqfree = (len(ch[-1]) == 1 and ch[-1][0] != 0)
                n_real = tb.sign_variations(tb.chain_signs_at(ch, "-inf")) - tb.sign_variations(tb.chain_signs_at(ch, "+inf"))
                if (not sqfree) or n_real > 0 or W[0] <= 0 or W[-1] <= 0:
                    wron_real_root.append((k, i, sqfree, n_real))
    print(f"  (c) Wronskian W(A_{{k,i-1}}, A_{{k,i}}) = A'_{{k,i-1}} A_{{k,i}} - A_{{k,i-1}} A'_{{k,i}} (k=1..{max_k_coeff}, {wron_checked} pairs, {time.perf_counter()-t2:.2f}s):")
    print(f"      Negative x-coeffs: {len(wron_neg_coeff)} {wron_neg_coeff[:5]}")
    print(f"      Zero x-coeffs:     {len(wron_zero_coeff)} {wron_zero_coeff[:5]}")
    print(f"      Exact Sturm > 0 on ALL of R (k=1..{max_k_sturm}, {wron_sturm_checked} pairs): FAIL={len(wron_real_root)} {wron_real_root[:5]}")

    # (d) BONUS adversarial check: What about t-coefficients of M_t^2 - X_t Y_t,
    #     htilde_{k,i}^2 - htilde_{k,i-1} htilde_{k,i+1}, and Wronskian W(htilde_{k,i-1}, htilde_{k,i})?
    print(f"\n=== 3. Bonus check: t-coefficients & higher-distance Wronskians ===")
    wron_general_neg = []
    for k in range(2, 16):
        for i in range(k):
            for j in range(i + 1, k + 1):
                # W(A_{k,i}, A_{k,j}) = A'_{k,i} A_{k,j} - A_{k,i} A'_{k,j}, deg A_{k,i} > deg A_{k,j}
                W = tb.poly_sub(
                    tb.poly_mul(tb.poly_deriv(A[(k, i)]), A[(k, j)]),
                    tb.poly_mul(A[(k, i)], tb.poly_deriv(A[(k, j)])),
                )
                if any(c < 0 for c in W):
                    wron_general_neg.append((k, i, j))
    print(f"  General Wronskian W(A_{{k,i}}, A_{{k,j}}) (0 <= i < j <= k <= 15) negative x-coeffs: {len(wron_general_neg)} {wron_general_neg[:5]}")


def main():
    MAX_K = 24
    t0 = time.perf_counter()
    A = tb.compute_A_table(MAX_K)
    h = tb.compute_h_table(10)
    print(f"Computed exact A_{{k,i}} table up to k={MAX_K} in {time.perf_counter()-t0:.2f}s.\n")
    verify_B10_disproof(A)
    print()
    verify_exact_t_to_x_M2_XY(A, h, max_k_check=10)
    stress_test_positivity(A, max_k_coeff=MAX_K, max_k_sturm=13)


if __name__ == "__main__":
    main()
