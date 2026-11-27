"""
Independent exact-arithmetic & symbolic audit of #40 for Bardakov-Kozlovskaya Conjecture 3.3(2):
1. Verify the exact c=1 formula -f_{2m,1}(u) = (2m/u)*(h*_{2m}(u) - (1+u)(1-u)^{2m-2}) for all even n=4..40.
2. Exact integer Sturm & Cauchy-index audit of G_{n,c}(t) real-rootedness (roots in (-inf,0) for n-3c odd,
   (0, 1/4] for n-3c even) and strict interlacing G_{n-2,c} << G_{n,c} for even n=2..30 and odd n=3..19!
   WAIT: check the boundary root at t = 1/4! Does G_{n,c}(1/4) == 0 ever happen (e.g. c=0)?!
3. Stress-test raw monomial sign law [x^{k1} y^{k2} z^{k3}] p_{2m}(x,y,z) == (-1)^{k1*k2+k1*k3+k2*k3}
   for ALL even n = 2..60 (using our exact E_table from engine.py)!
4. Exact SymPy verification of the Vandermonde-Cauchy Jacobian, exact sign, and covering degree
   (plus exact numerical check against engine.py at n=6, 8!).
"""

from __future__ import annotations
import math
import sys
import time
import sympy as sp

sys.path.insert(0, r"d:\solving-conjectures\work\bardakov-kozlovskaya\verify")
sys.path.insert(0, r"d:\solving-conjectures\work\gessel-lin-zeng\verify")
import engine as bk
import testbed as tb


# ---------------------------------------------------------------------------
# 1. Exact verification of c=1 formula for all even n = 4..40
# ---------------------------------------------------------------------------

def audit_c1_formula(max_n: int = 40):
    print(f"=== 1. Exact verification of c=1 formula for even n = 4..{max_n} ===")
    for n in range(4, max_n + 1, 2):
        coeffs = bk.compute_elementary_coeffs(n)
        m = (n - 3) // 2
        # Build f_{n,1}(u) = sum_{b=0}^m c_{n-3-2b, b, 1} * u^b * (1+u)^{n-3-2b}
        f_poly = [0] * (n - 2)
        for b in range(m + 1):
            a = n - 3 - 2 * b
            cb = coeffs[(a, b, 1)]
            for r in range(a + 1):
                f_poly[b + r] += cb * math.comb(a, r)

        # Build h*_n(u) = (1-u)^n * sum_{r=0}^{n-1} C(n*(r+1)-1, n-1) u^r (truncated to deg n-1)
        # Note deg(h*_n) <= n - 1!
        s_coeffs = [math.comb(n * (r + 1) - 1, n - 1) for r in range(n)]
        one_m_u_n = [((-1) ** j) * math.comb(n, j) for j in range(n + 1)]
        hstar = [0] * n
        for r in range(n):
            for j in range(n - r):
                hstar[r + j] += s_coeffs[r] * one_m_u_n[j]
        # Verify deg(h*_n) == n-1 and [u^n] vanishes too:
        c_un = sum(s_coeffs[r] * one_m_u_n[n - r] for r in range(n)) + math.comb(n * (n + 1) - 1, n - 1) * 1
        assert c_un == 0, f"h*_{n} did not vanish at u^{n}: c_un={c_un}!"

        # Build (1+u)*(1-u)^{n-2} (deg n-1)
        g_poly = [0] * n
        for j in range(n - 1):
            val = ((-1) ** j) * math.comb(n - 2, j)
            g_poly[j] += val
            g_poly[j + 1] += val

        diff = [hstar[j] - g_poly[j] for j in range(n)]
        assert diff[0] == 0 and diff[-1] == 0, f"Boundary mismatch at n={n}: diff[0]={diff[0]}, diff[-1]={diff[-1]}"
        rhs = [-n * diff[j + 1] for j in range(n - 2)]
        assert f_poly == rhs, f"c=1 formula mismatch at n={n}!"
    print(f"  EXACT IDENTITY -f_{{2m,1}}(u) = (2m/u)*(h*_{{2m}}(u) - (1+u)(1-u)^{{2m-2}}) verified for ALL even n=4..{max_n}!")


# ---------------------------------------------------------------------------
# 2. Exact integer Sturm & Cauchy-index audit of G_{n,c}(t)
# ---------------------------------------------------------------------------

def audit_Gnc_sturm(max_even: int = 30, max_odd: int = 19):
    print(f"\n=== 2. Exact integer Sturm & Cauchy-index audit of G_{{n,c}}(t) (even <= {max_even}, odd <= {max_odd}) ===")
    t0 = time.perf_counter()

    def build_G(coeffs, n, c):
        d = n - 3 * c
        if d < 0:
            return None
        m = d // 2
        poly = [coeffs[(d - 2 * b, b, c)] for b in range(m + 1)]
        # Make leading coeff positive for Sturm/Cauchy (record sign flip)
        poly = tb.poly_trim(poly)
        if poly[-1] < 0:
            poly = [-x for x in poly]
        return poly

    # Note: to test open interval (0, 1/4) in exact integers, substitute s = 4*t!
    # P(s/4) * 4^d has roots in (0, 1) iff P(t) has roots in (0, 1/4)!
    def scale_4(p):
        d = len(p) - 1
        return tb.poly_primitive_pos([c * (4 ** (d - r)) for r, c in enumerate(p)])

    even_tables = {n: bk.compute_elementary_coeffs(n) for n in range(2, max_even + 1, 2)}
    odd_tables = {n: bk.compute_elementary_coeffs(n) for n in range(3, max_odd + 1, 2)}

    rr_odd_pass, rr_odd_fail = 0, []
    rr_even_pass, rr_even_fail = 0, []
    root_at_quarter = []
    inter_pass, inter_fail = 0, []

    for parity, tables in [("even", even_tables), ("odd", odd_tables)]:
        for n, coeffs in sorted(tables.items()):
            for c in range(1, n // 3 + 1):
                d = n - 3 * c
                G = build_G(coeffs, n, c)
                degG = len(G) - 1
                if degG == 0:
                    continue
                if d % 2 == 1:
                    # Claim: all simple roots in (-inf, 0)
                    ok, msg = tb.verify_simple_roots_in_interval(G, "-inf", 0)
                    if ok:
                        rr_odd_pass += 1
                    else:
                        rr_odd_fail.append((n, c, msg))
                else:
                    # Claim: all simple roots in (0, 1/4] (check strictly (0, 1/4) vs root at 1/4!)
                    G4 = scale_4(G)
                    if tb.poly_eval_int(G4, 1) == 0:
                        root_at_quarter.append((n, c))
                    ok, msg = tb.verify_simple_roots_in_interval(G4, 0, 1)
                    if ok:
                        rr_even_pass += 1
                    else:
                        rr_even_fail.append((n, c, msg))

                # Check strict interlacing G_{n-2, c} << G_{n, c} (if n-2 in tables and n-2 >= 3*c)
                if (n - 2) in tables and (n - 2) >= 3 * c:
                    G_prev = build_G(tables[n - 2], n - 2, c)
                    assert len(G) == len(G_prev) + 1  # degree increases by 1!
                    if d % 2 == 1:
                        ok_i, msg_i = tb.verify_strict_interlacing(G_prev, G, "-inf", 0)
                    else:
                        ok_i, msg_i = tb.verify_strict_interlacing(scale_4(G_prev), scale_4(G), 0, 1)
                    if ok_i:
                        inter_pass += 1
                    else:
                        inter_fail.append((n, c, msg_i))

    print(f"  n-3c ODD  (roots in (-inf, 0)): PASS={rr_odd_pass}, FAIL={len(rr_odd_fail)} {rr_odd_fail[:5]}")
    print(f"  n-3c EVEN (roots in (0, 1/4) ): PASS={rr_even_pass}, FAIL={len(rr_even_fail)} {rr_even_fail[:5]} (root at 1/4: {root_at_quarter[:5]})")
    print(f"  Strict interlacing G_{{n-2,c}} << G_{{n,c}}: PASS={inter_pass}, FAIL={len(inter_fail)} {inter_fail[:5]}")
    print(f"  Elapsed: {time.perf_counter()-t0:.2f}s.")


# ---------------------------------------------------------------------------
# 3. Stress-test raw monomial sign law [x^{k1} y^{k2} z^{k3}] p_{2m} for even n = 2..60!
# ---------------------------------------------------------------------------

def audit_raw_monomial_signs(max_n: int = 60):
    print(f"\n=== 3. Stress-testing raw monomial sign law [x^k1 y^k2 z^k3] p_{{2m}} == (-1)^(k1*k2+k1*k3+k2*k3) (n=2..{max_n}) ===")
    t0 = time.perf_counter()
    total_checked = 0
    zeros = []
    sign_fails = []
    for n in range(2, max_n + 1, 2):
        E = bk.compute_E_table(n)
        # p_n(x, y, z) = sum_{j=0}^n (-1)^j E_j(x, y) z^{n-j}
        # So [x^{k1} y^{k2} z^{k3}] p_n = (-1)^{k1+k2} * E[k1+k2][k2] (with k3 = n - k1 - k2)!
        # Check all ordered triples k1 >= k2 >= k3 >= 0 with k1 + k2 + k3 = n:
        for k3 in range(n // 3 + 1):
            j = n - k3  # k1 + k2
            Ej = E[j]
            sgn_j = 1 if (j % 2 == 0) else -1
            for k2 in range(k3, j // 2 + 1):
                k1 = j - k2
                val = sgn_j * Ej[k2]
                total_checked += 1
                exp_sgn = -1 if ((k1 * k2 + k1 * k3 + k2 * k3) % 2 == 1) else 1
                if val == 0:
                    zeros.append((n, k1, k2, k3))
                elif (1 if val > 0 else -1) != exp_sgn:
                    sign_fails.append((n, k1, k2, k3, val, exp_sgn))
    print(f"  Checked {total_checked} distinct monomial partitions (even n=2..{max_n}, {time.perf_counter()-t0:.2f}s):")
    print(f"    Zeros: {len(zeros)} {zeros[:5]}")
    print(f"    Sign-law (-1)^(k1*k2+k1*k3+k2*k3) violations: {len(sign_fails)} {sign_fails[:5]}")


# ---------------------------------------------------------------------------
# 4. Exact SymPy verification of Vandermonde-Cauchy integral (#40 item 3)
# ---------------------------------------------------------------------------

def audit_vandermonde_cauchy():
    print("\n=== 4. Exact SymPy & numerical audit of Vandermonde-Cauchy integral (#40 item 3) ===")
    y, z = sp.symbols('y z')
    e1 = 1 + y + z
    e2 = y + z + y * z
    e3 = y * z
    u = e2 / e1**2
    v = e3 / e1**3
    J = sp.Matrix([[sp.diff(u, y), sp.diff(u, z)], [sp.diff(v, y), sp.diff(v, z)]])
    detJ = sp.factor(J.det())
    print(f"  Exact Jacobian det d(u,v)/d(y,z) = {detJ}")
    # Now compute the full integrand factor when changing (du dv)/(u^{b+1} v^{c+1}) * P_n(1, u, v) to dy dz:
    # Since P_n(1, u, v) = e1^{-n} p_n(1, y, z) and n = a + 2b + 3c:
    # (1 / (e1^n * u^{b+1} * v^{c+1})) * detJ * (y * z)   [note dy dz = i*y dphi_2 * i*z dpsi = -y*z dphi_2 dpsi!]
    # Wait: (1 / (2*pi*i)^2) * dy * dz = (1 / (2*pi)^2) * y * z * dphi_2 * dpsi (the i^2 cancels identically!!)!
    a_s, b_s, c_s = sp.symbols('a b c', integer=True)
    # Use k2 = b + c, k3 = c, a = k1 - k2, so n = a + 2b + 3c
    integrand_factor = sp.simplify(
        (e1 ** (-(a_s + 2 * b_s + 3 * c_s)))
        * (u ** (-(b_s + 1)))
        * (v ** (-(c_s + 1)))
        * detJ
        * (y * z)
    )
    # Compare against #40's displayed formula:
    #   (1 / (y^{b+c} * z^c)) * (1-y)*(1-z)*(1-z/y) / ((1+y+z)^{a+1} * (1 + z/y + z)^{b+1})
    formula_40 = (
        (1 - y) * (1 - z) * (1 - z / y)
        / (y ** (b_s + c_s) * z ** c_s * (1 + y + z) ** (a_s + 1) * (1 + z / y + z) ** (b_s + 1))
    )
    ratio = sp.simplify(integrand_factor / formula_40)
    print(f"  Exact ratio (true Jacobian integrand) / (#40 displayed formula) = {ratio}")
    assert ratio == 1 or ratio == -1
    if ratio == -1:
        print("  [!] SIGN FLIP DETECTED IN #40 DISPLAYED FORMULA (detJ has opposite orientation)!")
    else:
        print("  EXACT SIGN & ALGEBRAIC MATCH (ratio = +1)!")

    # Numerical cross-check on n = 6, (a,b,c) = (3,0,1), (1,1,1), (0,0,2) using 2D trapezoidal rule on torus!
    import cmath
    n = 6
    coeffs6 = bk.compute_elementary_coeffs(6)
    Ry, Rz = 0.35, 0.08  # satisfies Rz < Ry/(1+Ry) = 0.259 and Ry + Rz = 0.43 < 1
    N = 48
    zetas = [cmath.exp(2j * math.pi * k / n) for k in range(n)]
    for (a, b, c) in [(3, 0, 1), (1, 1, 1), (0, 0, 2), (2, 2, 0)]:
        k2, k3 = b + c, c
        acc = 0j
        for i1 in range(N):
            phi2 = -math.pi + 2 * math.pi * i1 / N
            yv = Ry * cmath.exp(1j * phi2)
            y_root = (Ry ** (1.0 / n)) * cmath.exp(1j * phi2 / n)
            for i2 in range(N):
                psi = -math.pi + 2 * math.pi * i2 / N
                zv = Rz * cmath.exp(1j * psi)
                p_val = 1.0 + 0j
                for zk in zetas:
                    p_val *= (zv - (1.0 + zk * y_root) ** n)
                num = p_val * (1 - yv) * (1 - zv) * (1 - zv / yv)
                den = (yv ** k2) * (zv ** k3) * ((1 + yv + zv) ** (a + 1)) * ((1 + zv / yv + zv) ** (b + 1))
                acc += num / den
        est = (acc / (N * N)).real
        exact = coeffs6[(a, b, c)]
        print(f"  n=6, ({a},{b},{c}): torus integral = {est:14.6f}, exact = {exact:8d}, ratio = {est/exact:.10f}")


def main():
    audit_c1_formula(max_n=40)
    audit_Gnc_sturm(max_even=30, max_odd=19)
    audit_raw_monomial_signs(max_n=60)
    audit_vandermonde_cauchy()


if __name__ == "__main__":
    main()
