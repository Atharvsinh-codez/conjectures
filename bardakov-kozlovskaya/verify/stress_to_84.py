"""
Push exact integer verification of Bardakov-Kozlovskaya Conjecture 3.3(2) & sign law
(-1)^{b(c+1)+c} through even n = 62..84 (plus exact odd-n check n = 1..31 and
small-c closed-form extraction / divisibility / sign-margin diagnostics).
"""

import math
import sys
import time
sys.path.insert(0, r"d:\solving-conjectures\work\bardakov-kozlovskaya\verify")
import engine as bk


def main():
    t0 = time.perf_counter()

    # 1. Exact check of odd n = 1..31 (verify c=0,b>0 vanish and odd sign law (-1)^{c(b+1)} for c>=1!)
    print("=== 1. Odd n = 1..31 exact check (c=0,b>0 zero; c>=1 sign law (-1)^(c(b+1))) ===")
    odd_zero_cge1 = []
    odd_sign_bad = []
    for n in range(1, 32, 2):
        coeffs = bk.compute_elementary_coeffs(n)
        for (a, b, c), val in coeffs.items():
            if c == 0:
                if b == 0:
                    assert val == 1
                else:
                    assert val == 0, f"Unexpected nonzero at odd n={n}, c=0, b={b}"
            else:
                exp_sgn = -1 if ((c * (b + 1)) % 2 == 1) else 1
                if val == 0:
                    odd_zero_cge1.append((n, a, b, c))
                elif (1 if val > 0 else -1) != exp_sgn:
                    odd_sign_bad.append((n, a, b, c, val, exp_sgn))
    print(f"  Odd n=1..31: c=0,b>0 ALL zero; c>=1 zeros={len(odd_zero_cge1)} {odd_zero_cge1[:5]}, "
          f"sign-law (-1)^(c(b+1)) violations={len(odd_sign_bad)} {odd_sign_bad[:5]}")

    # 2. Extract exact c=1 slice formula for even n=4..16!
    print("\n=== 2. Exact c=1 slice & (n-3, 0, 1) inspection (even n=4..20) ===")
    for n in range(4, 22, 2):
        coeffs = bk.compute_elementary_coeffs(n)
        val0 = coeffs[(n - 3, 0, 1)]
        # Check relation to 2^n, 3^n, binom, etc.
        # Note: at y=z=0, p_n(x,0,0)=x^n; what is [e_1^{n-3} e_3] p_n(x,y,z)?
        # Remember: at c=1, e_1^{n-3} e_3 is the ONLY term containing x^{n-2} y z that isn't in c=0!
        # Wait: c=0 is (e_1^2 - 4 e_2)^{n/2}! Its coeff of x^{n-2} y z is easy to compute!
        # Let's check exact ratios and prime factors!
        print(f"  n={n:2d}: c_{{n-3,0,1}} = {val0:15d},  -val/(n^2) = {-val0/(n*n):.4f},  -val/(4^n) = {-val0/(4**n):.6f}")

    # 3. Push even n = 62..84 (all even n!) and check nonvanishing + exact sign law (-1)^{b(c+1)+c}!
    print("\n=== 3. Pushing exact even n = 62..84 (stress-testing nonvanishing & sign law) ===")
    total_coeffs = 0
    zero_hits = []
    sign_bad = []
    for n in range(62, 85, 2):
        tn = time.perf_counter()
        coeffs = bk.compute_elementary_coeffs(n)
        n_zeros = 0
        n_sbad = 0
        for (a, b, c), val in coeffs.items():
            total_coeffs += 1
            exp_sgn = -1 if ((b * (c + 1) + c) % 2 == 1) else 1
            if val == 0:
                n_zeros += 1
                zero_hits.append((n, a, b, c))
            elif (1 if val > 0 else -1) != exp_sgn:
                n_sbad += 1
                sign_bad.append((n, a, b, c, val, exp_sgn))
        dt = time.perf_counter() - tn
        print(f"  n={n:2d}: {len(coeffs):4d} coeffs, zeros={n_zeros}, sign_fails={n_sbad} ({dt:.2f}s)")

    print(f"\nSUMMARY even n=62..84 ({total_coeffs} coeffs, total elapsed {time.perf_counter()-t0:.2f}s):")
    print(f"  Zero hits: {len(zero_hits)} {zero_hits[:5]}")
    print(f"  Sign-law violations: {len(sign_bad)} {sign_bad[:5]}")

    # 4. Also check 1 exact large spot check: n = 90 and n = 100!
    for n in [90, 100]:
        tn = time.perf_counter()
        coeffs = bk.compute_elementary_coeffs(n)
        n_zeros = sum(1 for v in coeffs.values() if v == 0)
        n_sbad = sum(1 for (a, b, c), v in coeffs.items() if (1 if v > 0 else -1) != (-1 if ((b * (c + 1) + c) % 2 == 1) else 1))
        print(f"  SPOT CHECK n={n:3d}: {len(coeffs):4d} coeffs, zeros={n_zeros}, sign_fails={n_sbad} ({time.perf_counter()-tn:.2f}s)")


if __name__ == "__main__":
    main()
