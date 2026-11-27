"""
Independent exact-arithmetic & high-precision audit of #50 for Bardakov-Kozlovskaya Conjecture 3.3(2):
1. Verify (-1)^j * Q_c^{(2m)}(1, u_j) > 0 for ALL even n = 4..60, ALL 1 <= c <= floor((2m-2)/3),
   ALL j = c..m - ceil(c/2), using exact engine.py integer polynomials + 100-dps mpmath.
2. Verify the exact 1D circle integral formula (including exact (-1)^{j+c} sign!) against exact
   (-1)^j * (1+t_j)^{2m-3c} * Q_c^{(2m)}(1, u_j) at n = 10, 12, 20 across even/odd c and corners.
3. Quantify the exact perturbation size R_z and cubic-root shift on the INTERIOR sample points
   (dropping j = c [replaced by u=0 via proved b=0] and even-c j = m-c/2 [replaced by u=+inf via proved a=0]).
"""

from __future__ import annotations
import math
import sys
import time
import mpmath as mp

sys.path.insert(0, r"d:\solving-conjectures\work\bardakov-kozlovskaya\verify")
import engine as bk


def main():
    mp.mp.dps = 100
    t0 = time.perf_counter()

    print("=== 1. Verifying (-1)^j Q_c^{(2m)}(1, u_j) > 0 for ALL even n = 4..60 ===")
    total_slices = 0
    total_points = 0
    fails = []
    for n in range(4, 61, 2):
        m = n // 2
        coeffs = bk.compute_elementary_coeffs(n)
        max_c = (n - 2) // 3
        for c in range(1, max_c + 1):
            d = (n - 3 * c) // 2
            if d < 1:
                continue
            total_slices += 1
            Q = [mp.mpf(coeffs[(n - 3 * c - 2 * b, b, c)]) for b in range(d + 1)]
            j_max = m - (c + 1) // 2
            assert j_max - c == d
            prev_u = None
            for j in range(c, j_max + 1):
                total_points += 1
                s = mp.sin(mp.pi * j / n) / mp.sin(mp.pi * (j + c) / n)
                tj = (s ** n) if (c % 2 == 0) else -(s ** n)
                uj = tj / ((1 + tj) ** 2)
                # Verify strict monotonicity
                if prev_u is not None:
                    if c % 2 == 0:
                        assert uj > prev_u
                    else:
                        assert uj < prev_u
                prev_u = uj
                # Evaluate Q(1, uj)
                val = mp.mpf(0)
                for coeff in reversed(Q):
                    val = val * uj + coeff
                signed = val if (j % 2 == 0) else -val
                if not (signed > 0):
                    fails.append((n, c, j, float(signed)))
    print(f"  Checked {total_slices} slices, {total_points} sample points across ALL even n=4..60 ({time.perf_counter()-t0:.2f}s):")
    print(f"  Sign (-1)^j Q_c(1, u_j) > 0 FAILS: {len(fails)} {fails[:5]}")
    assert len(fails) == 0

    print("\n=== 2. Exact numerical verification of 1D circle integral identity (n = 10, 12, 16) ===")
    mp.mp.dps = 60
    N_pts = 48
    for n, test_cases in [
        (10, [(2, 2), (2, 3), (2, 4), (1, 1), (1, 3)]),
        (12, [(2, 2), (2, 3), (2, 5), (3, 3), (3, 4)]),
        (16, [(2, 2), (2, 7), (3, 3), (3, 5), (4, 4), (4, 5), (4, 6)]),
    ]:
        m = n // 2
        coeffs = bk.compute_elementary_coeffs(n)
        for c, j in test_cases:
            d = (n - 3 * c) // 2
            Q = [mp.mpf(coeffs[(n - 3 * c - 2 * b, b, c)]) for b in range(d + 1)]
            s = mp.sin(mp.pi * j / n) / mp.sin(mp.pi * (j + c) / n)
            tj = (s ** n) if (c % 2 == 0) else -(s ** n)
            uj = tj / ((1 + tj) ** 2)
            Rz = (mp.sin(mp.pi * c / n) / mp.sin(mp.pi * (j + c) / n)) ** n

            # Exact LHS = (-1)^j * (1 + tj)^{n - 3c} * Q_c(1, uj)
            qval = mp.mpf(0)
            for coeff in reversed(Q):
                qval = qval * uj + coeff
            lhs = ((-1) ** j) * ((1 + tj) ** (n - 3 * c)) * qval

            # Compute 1D circle integral directly using P_n(1+tj, tj, tj*Rz*e^{i*phi}) via exact elementary coeffs!
            e1 = 1 + tj
            e2 = tj
            # Precompute Q_r(e1, e2) for all r = 0..n//3
            Qr_vals = []
            for r in range(n // 3 + 1):
                dr = (n - 3 * r) // 2
                qr = mp.mpf(0)
                for b in range(dr + 1):
                    a = n - 3 * r - 2 * b
                    qr += mp.mpf(coeffs[(a, b, r)]) * (e1 ** a) * (e2 ** b)
                Qr_vals.append(qr)

            acc = mp.mpc(0)
            acc_abs = mp.mpf(0)
            min_re = mp.mpf('1e100')
            for idx in range(N_pts):
                phi = -mp.pi + 2 * mp.pi * idx / N_pts
                eip = mp.exp(1j * phi)
                e3 = tj * Rz * eip
                P_val = mp.mpc(0)
                e3_pow = mp.mpc(1)
                for qr in Qr_vals:
                    P_val += qr * e3_pow
                    e3_pow *= e3
                F = ((-1) ** (j + c)) * P_val / ((abs(tj) * Rz * eip) ** c)
                acc += F
                acc_abs += abs(F)
                if F.real < min_re:
                    min_re = F.real
            rhs = (acc / N_pts).real
            ratio = rhs / lhs
            pos_ratio = (acc / N_pts).real / (acc_abs / N_pts)
            assert abs(ratio - 1) < mp.mpf('1e-35'), f"1D integral mismatch at n={n}, c={c}, j={j}: ratio={ratio}"
            print(f"  n={n:2d}, c={c}, j={j}: exact ratio = {mp.nstr(ratio, 12)}, Re(avg)/avg|F| = {float(pos_ratio):.5f}, minRe/avg = {float(min_re/rhs):+.5f}, Rz = {float(Rz):.3e}")

    print("\n=== 3. Checking max R_z on STRICTLY INTERIOR points (c+1 <= j <= m - ceil(c/2) - (1 if c even else 0)) ===")
    for n in [20, 40, 80, 160]:
        m = n // 2
        max_rz_all = (0.0, None)
        max_rz_int = (0.0, None)
        for c in range(2, (n - 2) // 3 + 1):
            j_lo = c + 1
            j_hi = m - (c + 1) // 2 - (1 if c % 2 == 0 else 0)
            for j in range(c, m - (c + 1) // 2 + 1):
                rz = (math.sin(math.pi * c / n) / math.sin(math.pi * (j + c) / n)) ** n
                if rz > max_rz_all[0]:
                    max_rz_all = (rz, (c, j))
                if j_lo <= j <= j_hi and rz > max_rz_int[0]:
                    max_rz_int = (rz, (c, j))
        print(f"  n={n:3d}: max R_z (all) = {max_rz_all[0]:.5f} at {max_rz_all[1]}; "
              f"max R_z (interior, dropping proved u=0 & even u=inf) = {max_rz_int[0]:.5f} at {max_rz_int[1]}")


if __name__ == "__main__":
    main()
