"""
Independent exact-arithmetic audit of #31, #32, and #33 for Gessel-Lin-Zeng Conjecture 15:
1. Re-verify #31 item 1:
   (a) Exact integer Sturm & Cauchy-index audit of single-operator witness
       (K=1000 copies of b=2, b_1001=196200000000000, b_1002=181330000000000, i=960).
   (b) Exact SymPy derivation of single- and double-operator scaling limits M_sgl^{(w)} and M_dbl^{(w)}.
   (c) CRITICAL ADVERSARIAL CHECK: For the double-operator prod D_{b_j}(D_{b_j}+y)(1),
       does M_dbl < 0 actually break interlacing AFTER the outer D_{b_2} (i.e. W(D_{b_2} R_{i+1}, D_{b_2} R_i) < 0),
       or can the outer D_{b_2} heal the inner interlacing break? Test on exact scaling limit & finite K!
2. Audit #33: Exact signs of A, B, C, and J_raw in
     (D_{k+2} htilde_{k,i})^2 - (D_{k+2} htilde_{k,i-1})(D_{k+2} htilde_{k,i+1}) = A*Phi^2 + B*Phi + C
   with Phi = 1 + (k+2)*t on BOTH regimes (-1, -1/(k+2)] and [-1/(k+2), 0) for GLZ k = 2..30!
3. Independently verify #32 item 2: Ma-Wang Theorem 2 sharpness Laguerre certificates for
   i=4 (k=31 vs k=30), i=5 (k=16 vs k=15), and i=k-3 (k=28 vs k=27).
"""

from __future__ import annotations
from fractions import Fraction
import math
import sys
import time
import sympy as sp

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import testbed as tb


def poly_eval_frac(p: list[int], x: Fraction) -> Fraction:
    acc = Fraction(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc


# ---------------------------------------------------------------------------
# 1. Audit #31 item 1 (single-operator witness + double-operator scaling & outer D_{b_2})
# ---------------------------------------------------------------------------

def audit_31_item1():
    print("=== 1(a). Exact integer Sturm & Cauchy-index audit of #31 single-operator witness ===")
    K = 1000
    i0 = 960
    # Compute D_2^m(1) for m = 0..43 using exact Python integers
    D2 = [[1]]
    for _ in range(44):
        D2.append(tb.op_D(2, D2[-1]))

    # Note: [y^r] (D_2 + y)^K (1) = binom(K, r) * D_2^{K-r}(1).
    # Divide out common gcd of binomial coeffs for r in range(i0-3, i0+3) to keep exact ints compact
    # (does not change any roots or Wronskian signs!).
    cs = {r: math.comb(K, r) for r in range(i0 - 3, i0 + 3)}
    g = math.gcd(*cs.values())
    P0 = {r: [(cs[r] // g) * c for c in D2[K - r]] for r in range(i0 - 3, i0 + 3)}

    b1 = 196200000000000
    b2 = 181330000000000
    assert b1 >= 2 and b2 >= 2

    P1 = {r: tb.poly_add(tb.op_D(b1, P0[r]), P0[r - 1]) for r in range(i0 - 2, i0 + 3)}
    P2 = {r: tb.poly_add(tb.op_D(b2, P1[r]), P1[r - 1]) for r in (i0 - 1, i0, i0 + 1)}

    p60 = P2[i0]
    p61 = P2[i0 + 1]
    print(f"  deg(P2[960]) = {len(p60)-1}, deg(P2[961]) = {len(p61)-1}")
    ok60, msg60 = tb.verify_simple_roots_in_interval(p60, -1, 0)
    ok61, msg61 = tb.verify_simple_roots_in_interval(p61, -1, 0)
    print(f"  Exact Sturm P2[960] in (-1,0): {ok60} ({msg60})")
    print(f"  Exact Sturm P2[961] in (-1,0): {ok61} ({msg61})")
    assert ok60 and ok61

    ok_inter, msg_inter = tb.verify_strict_interlacing(p61, p60, -1, 0)
    print(f"  Exact Cauchy-index interlacing P2[961] << P2[960]: {ok_inter} ({msg_inter})")
    assert not ok_inter

    # Check exact W(P2[961], P2[960]) at t = -1818 / 10^16
    t_pt = Fraction(-1818, 10**16)
    v61 = poly_eval_frac(p61, t_pt)
    dv61 = poly_eval_frac(tb.poly_deriv(p61), t_pt)
    v60 = poly_eval_frac(p60, t_pt)
    dv60 = poly_eval_frac(tb.poly_deriv(p60), t_pt)
    W_val = v61 * dv60 - dv61 * v60
    print(f"  Exact W(P2[961], P2[960])(-1818/10^16) = {float(W_val):.6e} (EXACT < 0? {W_val < 0})")
    assert W_val < 0

    print("\n=== 1(b). Exact SymPy derivation of single- & double-operator scaling limits ===")
    Phi, w, a, u = sp.symbols('Phi w a u')
    B_val = (1 - Phi) / a
    def D_lim(f):
        return sp.expand((1 - B_val * u) * f + u * sp.diff(f, u))

    # p_m(u) = 1 - 2^{2-m} u  (m = 0,1,2,3 -> slopes 4, 2, 1, 1/2)
    pm = {m: 1 - sp.Rational(4, 2**m) * u for m in (0, 1, 2, 3)}
    # Single operator: q_m = D_lim(p_m) + w * p_{m-1}
    qs = {m: sp.expand(D_lim(pm[m]) + w * pm[m - 1]) for m in (1, 2, 3)}
    M_sgl = sp.expand(-u * (qs[3] * sp.diff(qs[2], u) - sp.diff(qs[3], u) * qs[2]) + w * (qs[2]**2 - qs[1] * qs[3])).subs(u, a)
    M_sgl_c1 = sp.collect(sp.expand(sp.series(M_sgl, a, 0, 2).removeO() / a), Phi)
    expected_sgl = (1 + w) * Phi**2 + (2 + 2*w + 3*w**2) * Phi + (1 + 4*w + 2*w**2 + 2*w**3)
    print(f"  Single-op M_sgl/alpha = {M_sgl_c1}")
    print(f"  Ratio M_sgl_c1 / expected_sgl = {sp.simplify(M_sgl_c1 / expected_sgl)}")
    disc_sgl = sp.expand((2 + 2*w + 3*w**2)**2 - 4*(1 + w)*(1 + 4*w + 2*w**2 + 2*w**3))
    print(f"  Single-op disc(Phi) = {sp.factor(disc_sgl)}; at w=5: {disc_sgl.subs(w,5)}, at w=6: {disc_sgl.subs(w,6)}")

    # Double operator: q_m = D_lim(D_lim(p_m) + w * p_{m-1})
    qd = {m: sp.expand(D_lim(D_lim(pm[m]) + w * pm[m - 1])) for m in (1, 2, 3)}
    M_dbl = sp.expand(-u * (qd[3] * sp.diff(qd[2], u) - sp.diff(qd[3], u) * qd[2]) + w * (qd[2]**2 - qd[1] * qd[3])).subs(u, a)
    M_dbl_c1 = sp.collect(sp.expand(sp.series(M_dbl, a, 0, 2).removeO() / a), Phi)
    expected_dbl = sp.Rational(1, 2)*(w + 1)*Phi**4 + (sp.Rational(3, 2)*w**2 + 3*w + 1)*Phi**3 + (w**3 + sp.Rational(9, 2)*w**2 + 3*w + 1)*Phi**2 + (w**3 - sp.Rational(3, 2)*w - 2)*Phi + w + sp.Rational(3, 2)
    print(f"  Double-op M_dbl/alpha = {M_dbl_c1}")
    print(f"  Ratio M_dbl_c1 / expected_dbl = {sp.simplify(M_dbl_c1 / expected_dbl)}")
    # Evaluate exact rational witness for M_dbl/alpha < 0 at w=5!
    print(f"  Double-op at w=5, Phi=-2: {M_dbl_c1.subs({w: 5, Phi: -2})} (<0? {M_dbl_c1.subs({w: 5, Phi: -2}) < 0})")

    print("\n=== 1(c). ADVERSARIAL CHECK: Does outer D_{b_2} preserve or heal the double-op interlacing break? ===")
    # After step 1 (D_{b_1}(D_{b_1} + y)) we have Q_r = qd[r] (times w^{2-r}).
    # At step 2, we apply D_{b_2}(D_{b_2} + y):
    #   R_m = D_{b_2}(Q_m) + w Q_{m-1} (inner sum, whose interlacing breaks when M_dbl < 0 at a root of R_2!)
    #   S_m = D_{b_2}(R_m)             (AFTER outer D_{b_2}!).
    # Wait: let's compute W_u(S_3, S_2) in the EXACT scaling limit at u = a when b_2 is chosen so R_2(a) = 0!!
    # At u = a, let D_2 f = (1 - B2*u)*f + u*f' with B2 = (1 - Phi2)/a.
    # For R_2(a) = 0: D_2(qd[2])(a) + w*qd[1](a) = 0 => Phi2*qd[2](a) + a*qd[2]'(a) + w*qd[1](a) = 0!
    # Wait: since qd[2](a) = Phi*(Phi + w) + O(a) and qd[1](a) = Phi*(Phi + w) + O(a),
    # R_2(a) = 0 gives Phi2 = -w + O(a)!!
    # Let's check exact W(S_3, S_2) at u = a (AND nearby u = a*(1 + eta)) for w=10, Phi=-5!
    # Even better: test an ACTUAL FINITE INTEGER double-operator witness!
    # Wait: for double-operator (D_2(D_2 + y))^K (1) = D_2^K (D_2 + y)^K (1)!
    # Notice D_2 commutes with itself, so [y^r] (D_2(D_2+y))^K (1) = binom(K, r) * D_2^{2K - r}(1)!!
    # Let's test K = 120, i0 = 110 (so 2K - i0 = 130, w_eff = 110/11 = 10, deg ~ 130, SUPER FAST in exact integers!)
    # Wait: first check in scaling limit whether W(S_3, S_2) goes NEGATIVE near u = a!
    Phi2 = sp.Symbol('Phi2')
    B2_val = (1 - Phi2) / a
    def D2_lim(f):
        return sp.expand((1 - B2_val * u) * f + u * sp.diff(f, u))

    # Substitute w = 10, Phi = -5 (where M_dbl/a is deeply negative!)
    qd10 = {m: sp.expand(qd[m].subs({w: 10, Phi: -5})) for m in (1, 2, 3)}
    R2 = sp.expand(D2_lim(qd10[2]) + 10 * qd10[1])
    R3 = sp.expand(D2_lim(qd10[3]) + 10 * qd10[2])
    # Solve R2(u=a) = 0 for Phi2 as series in a!
    R2_at_a = sp.expand(R2.subs(u, a))
    Phi2_sol = sp.solve(R2_at_a, Phi2)[0]
    # Now compute S2 = D2_lim(R2), S3 = D2_lim(R3), and check W_t(S3, S2) = -(S3 S2_u - S3_u S2) at u = a (with Phi2 = Phi2_sol)!
    # Wait: at u = a, R2(a) = 0, so S2(a) = a * R2'(a)!
    # Better: check W(R3, R2) vs W(S3, S2) at u = a * (1 + c*a) or u = a * lam!
    S2 = sp.expand(D2_lim(R2))
    S3 = sp.expand(D2_lim(R3))
    # Leading coeff in a of W_t(R3, R2)(u=a) and W_t(S3, S2)(u=a) at Phi2 = Phi2_sol!
    W_R = sp.expand(-(R3 * sp.diff(R2, u) - sp.diff(R3, u) * R2)).subs(u, a).subs(Phi2, Phi2_sol)
    W_S = sp.expand(-(S3 * sp.diff(S2, u) - sp.diff(S3, u) * S2)).subs(u, a).subs(Phi2, Phi2_sol)
    print(f"  Scaling limit (w=10, Phi=-5, R2(a)=0): leading W(R3,R2)(a) = {sp.series(W_R, a, 0, 2)}")
    print(f"  Scaling limit (w=10, Phi=-5, R2(a)=0): leading W(S3,S2)(a) = {sp.series(W_S, a, 0, 2)}")
    # Also scan Phi2 near -10 (since S2(a) = 0 happens at a SLIGHTLY different Phi2 than R2(a) = 0!)
    # Check W_S at u = a as a function of Phi2 (leading O(1/a) or O(1) term)!
    W_S_general = sp.series(sp.expand(-(S3 * sp.diff(S2, u) - sp.diff(S3, u) * S2)).subs(u, a), a, 0, 1).removeO()
    print(f"  Leading W(S3, S2)(a) as function of Phi2 (w=10, Phi=-5):")
    for p2v in [-12, -10.5, -10, -9.5, -8, -5, -2, -1, 0]:
        val = float(W_S_general.subs({Phi2: p2v, a: 1e-4}))
        print(f"    Phi2={p2v:5.1f}: W(S3,S2)(a) ~ {val:+.4e}")


# ---------------------------------------------------------------------------
# 2. Audit #33: Exact signs of A, B, C in GLZ D_{k+2}-Turan on both regimes (k=2..30)
# ---------------------------------------------------------------------------

def audit_33_ABC_signs(max_k: int = 30):
    print(f"\n=== 2. Audit #33: Signs of A, B, C in (D_{{k+2}} htilde_{{k,i}})^2 - (D_{{k+2}} htilde_{{k,i-1}})(D_{{k+2}} htilde_{{k,i+1}}) = A*Phi^2 + B*Phi + C (k=2..{max_k}) ===")
    t0 = time.perf_counter()
    h = tb.compute_h_table(max_k)
    ht = {(k, i): tb.poly_div_exact_linear(h[(k, i)], -1) for k in range(1, max_k + 1) for i in range(k + 1)}

    # Note:
    # A(t) = p1^2 - p0*p2
    # B(t) = t(1+t)*(2*p1*p1' - p0'*p2 - p0*p2') = t(1+t) * A'(t) !!
    # C(t) = t^2(1+t)^2 * (p1'^2 - p0'*p2')
    # WAIT A SECOND!! B(t) = t(1+t) * A'(t)!!
    # Let's check BOTH regimes:
    #   Left regime:  t in (-1, -1/(k+2)]   (where Phi = 1 + (k+2)*t <= 0)
    #   Right regime: t in [-1/(k+2), 0)    (where Phi = 1 + (k+2)*t >= 0)
    # Use exact Sturm on (-1, -1/(k+2)) and (-1/(k+2), 0) for k <= 12 (substituting t = s/(k+2) so -1/(k+2) is integer -1!),
    # AND dense exact Fraction grid for ALL k = 2..30 (especially near t = 0^-, t = -1/(k+2), and interior roots)!

    def eval_scaled(poly, num, den):
        # Evaluate poly(num/den) * den^deg as exact int
        d = len(poly) - 1
        acc = 0
        den_pow = 1
        # Horner with common denominator
        # Actually Fraction is plenty fast for k <= 30, or integer Horner:
        # sum c_r * num^r * den^{d-r}
        res = 0
        for c in reversed(poly):
            res = res * num + c * den_pow
            den_pow *= den
        # Wait: let's verify! For [c0, c1, c2]:
        # r=2: res = c2, den_pow = den
        # r=1: res = c2*num + c1*den, den_pow = den^2
        # r=0: res = c2*num^2 + c1*num*den + c0*den^2! EXACT!
        return res

    stats = {
        "A_left_neg": [], "A_right_neg": [],
        "B_left_neg": [], "B_left_pos": [],
        "B_right_neg": [], "B_right_pos": [],
        "C_left_neg": [], "C_right_neg": [],
        "J_left_lt1": [], "J_right_lt1": [],
        "full_left_neg": [], "full_right_neg": [],
    }

    for k in range(2, max_k + 1):
        m = k + 2  # Phi = 1 + m*t, split at t_split = -1/m
        # Sample exact rational points in (-1, -1/m] and [-1/m, 0)
        # Left points: -1 + 1/(10m), -9/10, -3/4, -1/2, -2/m, -11/(10m), -1/m
        left_pts = [(-10*m + 1, 10*m), (-9, 10), (-3, 4), (-1, 2), (-1, 3), (-2, m), (-11, 10*m), (-1, m)]
        left_pts = [(p, q) for p, q in left_pts if -q < p and p * m <= -q]
        # Right points: -1/m, -9/(10m), -1/(2m), -1/(5m), -1/(10m), -1/(100m), -1/(10^6), -1/(10^12)
        right_pts = [(-1, m), (-9, 10*m), (-1, 2*m), (-1, 5*m), (-1, 10*m), (-1, 100*m), (-1, 10**6), (-1, 10**12)]
        right_pts = [(p, q) for p, q in right_pts if p * m >= -q and p < 0]

        for i in range(1, k):
            p0 = ht[(k, i - 1)]
            p1 = ht[(k, i)]
            p2 = ht[(k, i + 1)]
            dp0, dp1, dp2 = tb.poly_deriv(p0), tb.poly_deriv(p1), tb.poly_deriv(p2)

            A_poly = tb.poly_sub(tb.poly_mul(p1, p1), tb.poly_mul(p0, p2))
            dA_poly = tb.poly_deriv(A_poly)
            # B(t) = t(1+t) * A'(t) (let's assert exact identity!)
            B_raw = tb.poly_sub(
                tb.poly_mul([2], tb.poly_mul(p1, dp1)),
                tb.poly_add(tb.poly_mul(dp0, p2), tb.poly_mul(p0, dp2)),
            )
            assert B_raw == dA_poly
            C_raw = tb.poly_sub(tb.poly_mul(dp1, dp1), tb.poly_mul(dp0, dp2))  # C(t) = t^2(1+t)^2 * C_raw(t)
            # Wronskians for J_raw = 4*W01*W12 - W02^2
            W01 = tb.poly_sub(tb.poly_mul(p0, dp1), tb.poly_mul(dp0, p1))
            W12 = tb.poly_sub(tb.poly_mul(p1, dp2), tb.poly_mul(dp1, p2))
            W02 = tb.poly_sub(tb.poly_mul(p0, dp2), tb.poly_mul(dp0, p2))
            J_diff = tb.poly_sub(tb.poly_mul([4], tb.poly_mul(W01, W12)), tb.poly_mul(W02, W02))

            # Full D_{k+2}-Turan
            Dp0, Dp1, Dp2 = tb.op_D(m, p0), tb.op_D(m, p1), tb.op_D(m, p2)
            Full_poly = tb.poly_sub(tb.poly_mul(Dp1, Dp1), tb.poly_mul(Dp0, Dp2))

            for regime, pts in [("left", left_pts), ("right", right_pts)]:
                for num, den in pts:
                    # Note t = num/den in (-1, 0), so t(1+t) < 0, t^2(1+t)^2 > 0!
                    va = eval_scaled(A_poly, num, den)
                    # B(t) = t(1+t) * dA(t), so sgn(B) = -sgn(dA) on (-1, 0)!
                    vdA = eval_scaled(dA_poly, num, den)
                    sgn_B = -tb.sign_int(vdA)
                    vc = eval_scaled(C_raw, num, den)
                    vj = eval_scaled(J_diff, num, den)
                    vf = eval_scaled(Full_poly, num, den)
                    if va <= 0 and len(stats[f"A_{regime}_neg"]) < 5:
                        stats[f"A_{regime}_neg"].append((k, i, f"{num}/{den}", va))
                    if sgn_B < 0 and len(stats[f"B_{regime}_neg"]) < 5:
                        stats[f"B_{regime}_neg"].append((k, i, f"{num}/{den}"))
                    if sgn_B > 0 and len(stats[f"B_{regime}_pos"]) < 5:
                        stats[f"B_{regime}_pos"].append((k, i, f"{num}/{den}"))
                    if vc <= 0 and len(stats[f"C_{regime}_neg"]) < 5:
                        stats[f"C_{regime}_neg"].append((k, i, f"{num}/{den}"))
                    if vj <= 0 and len(stats[f"J_{regime}_lt1"]) < 5:
                        stats[f"J_{regime}_lt1"].append((k, i, f"{num}/{den}"))
                    if vf <= 0 and len(stats[f"full_{regime}_neg"]) < 5:
                        stats[f"full_{regime}_neg"].append((k, i, f"{num}/{den}"))

    print(f"  Rational-grid audit across ALL k=2..{max_k}, 1<=i<=k-1 ({time.perf_counter()-t0:.2f}s):")
    for key, val in stats.items():
        print(f"    {key:16s}: count(capped@5)={len(val)}  samples={val[:3]}")

    # Exact Sturm on (-1, -1/(k+2)) and (-1/(k+2), 0) for k=2..11 to confirm NO hidden sign flips between grid points!
    print("  Exact Sturm sign verification on BOTH intervals (-1, -1/(k+2)) and (-1/(k+2), 0) for k=2..11:")
    def poly_rescale_m(p, m):
        # q(s) = m^deg * p(s / m) => s in (-m, -1) is t in (-1, -1/m); s in (-1, 0) is t in (-1/m, 0)
        d = len(p) - 1
        return tb.poly_primitive_pos([c * (m ** (d - r)) for r, c in enumerate(p)])

    def roots_in(p, lo, hi):
        ch = tb.sturm_chain(p, tb.poly_deriv(p))
        # Handle boundary roots safely by just reporting Sturm count on open (lo, hi) when non-vanishing
        vlo = tb.sign_variations(tb.chain_signs_at(ch, lo))
        vhi = tb.sign_variations(tb.chain_signs_at(ch, hi))
        return vlo - vhi, (len(ch[-1]) == 1)

    sturm_summary = {name: [0, 0] for name in ["A_L", "A_R", "B_L", "B_R", "C_L", "C_R", "J_L", "J_R"]}
    first_sturm = {}
    for k in range(2, 12):
        m = k + 2
        for i in range(1, k):
            p0, p1, p2 = ht[(k, i - 1)], ht[(k, i)], ht[(k, i + 1)]
            dp0, dp1, dp2 = tb.poly_deriv(p0), tb.poly_deriv(p1), tb.poly_deriv(p2)
            A_poly = tb.poly_sub(tb.poly_mul(p1, p1), tb.poly_mul(p0, p2))
            dA_poly = tb.poly_deriv(A_poly)
            C_raw = tb.poly_sub(tb.poly_mul(dp1, dp1), tb.poly_mul(dp0, dp2))
            W01 = tb.poly_sub(tb.poly_mul(p0, dp1), tb.poly_mul(dp0, p1))
            W12 = tb.poly_sub(tb.poly_mul(p1, dp2), tb.poly_mul(dp1, p2))
            W02 = tb.poly_sub(tb.poly_mul(p0, dp2), tb.poly_mul(dp0, p2))
            J_diff = tb.poly_sub(tb.poly_mul([4], tb.poly_mul(W01, W12)), tb.poly_mul(W02, W02))
            for tag, poly in [("A", A_poly), ("B", dA_poly), ("C", C_raw), ("J", J_diff)]:
                qs = poly_rescale_m(poly, m)
                rL, sqL = roots_in(qs, -m, -1)
                rR, sqR = roots_in(qs, -1, 0)
                if rL == 0: sturm_summary[f"{tag}_L"][0] += 1
                else:
                    sturm_summary[f"{tag}_L"][1] += 1
                    first_sturm.setdefault(f"{tag}_L", (k, i, rL))
                if rR == 0: sturm_summary[f"{tag}_R"][0] += 1
                else:
                    sturm_summary[f"{tag}_R"][1] += 1
                    first_sturm.setdefault(f"{tag}_R", (k, i, rR))
    for k_name, (p_cnt, f_cnt) in sturm_summary.items():
        print(f"    {k_name}: zero-root pairs={p_cnt}, sign-change pairs={f_cnt} (first={first_sturm.get(k_name)})")


# ---------------------------------------------------------------------------
# 3. Independently verify #32 item 2: Ma-Wang sharpness exact Laguerre certificates
# ---------------------------------------------------------------------------

def audit_32_ma_wang_sharpness():
    print("\n=== 3. Independent exact-Fraction audit of #32 Ma-Wang sharpness certificates ===")
    cases = [
        ("i=4 (k=31)", 31, tuple(range(28, 32)), tuple(range(1, 5)), 934, Fraction(-17936903, 778480)),
        ("i=5 (k=16)", 16, tuple(range(12, 17)), tuple(range(1, 6)), 244, Fraction(-1538333, 145228)),
        ("i=k-3 (k=28)", 28, tuple(range(4, 29)), tuple(range(1, 26)), 271, Fraction(-1885999, 258846)),
    ]
    for label, k, S_late, S_early, coeff, x_pt in cases:
        f_late = tb.A_k_S_via_operator(k, S_late)
        f_early = tb.A_k_S_via_operator(k, S_early)
        R = tb.poly_add([coeff * c for c in f_late], f_early)
        dR = tb.poly_deriv(R)
        ddR = tb.poly_deriv(dR)
        v0 = poly_eval_frac(R, x_pt)
        v1 = poly_eval_frac(dR, x_pt)
        v2 = poly_eval_frac(ddR, x_pt)
        lag = v1 * v1 - v0 * v2
        print(f"  {label}: deg(R)={len(R)-1}, R'(x)^2 - R(x)R''(x) = {float(lag):.6e} (EXACT < 0? {lag < 0})")
        assert lag < 0

    # Also verify exact Cauchy index that k-1 STILL passes A_{k-1, late} << A_{k-1, early} for each!
    prev_cases = [
        ("i=4 (k=30)", 30, tuple(range(27, 31)), tuple(range(1, 5))),
        ("i=5 (k=15)", 15, tuple(range(11, 16)), tuple(range(1, 6))),
        ("i=k-3 (k=27)", 27, tuple(range(4, 28)), tuple(range(1, 25))),
    ]
    for label, k, S_late, S_early in prev_cases:
        f_late = tb.A_k_S_via_operator(k, S_late)
        f_early = tb.A_k_S_via_operator(k, S_early)
        ok, msg = tb.verify_strict_interlacing(f_late, f_early, "-inf", 0)
        print(f"  Previous k check {label}: A_late << A_early = {ok} ({msg})")
        assert ok


def main():
    audit_31_item1()
    audit_33_ABC_signs(max_k=30)
    audit_32_ma_wang_sharpness()


if __name__ == "__main__":
    main()
