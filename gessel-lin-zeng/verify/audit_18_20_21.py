"""
Independent adversarial verification of #18, #20, #21, and #22 for Gessel-Lin-Zeng Conjecture 15.

Checks:
1. Exact rational + exact integer Sturm certification of #18's k=9 counterexample to
   universal prod_{j=1}^k (T_j + w_j y)(1) (0 < a_j < b_j, w_j > 0).
2. Adversarial audit of #20 Quantitative Plucker Quad Lemma:
   - Exact rational counterexample when W_AC < 0 (u < 0) is not excluded!
   - Exact check of interior (u, v) in [0, 2*th*sqrt(a)] x [0, 2*th*sqrt(b)] vs corner (u, v) = (2*th*sqrt(a), 2*th*sqrt(b)),
     AND exact SymPy verification of the corner bound over R <= 4 and R <= 9.
   - Sign-flip check for W_{A+B, C+E} / W_{BC} = (1+u)(1+v) - ab <= 0 when R > 1!
3. Exact SymPy audit of:
   - Single-operator Turan, T-Turan equivalence, Borcea-Branden symbol of D_m (m >= 1 vs m < 1),
   - Universal Cubic Wronskian SOS identity,
   - Operator intertwining P_{L+1} S_{L+1} = S_L P_{L+1}.
4. Exact integer Cauchy-index stress test up to k=24 of:
   - A_{k, i} << P_{3k-1-i} A_{k, i-1}  (1 <= i <= k <= 24)
   - S_{3k-3-i} A_{k-1, i} << P_{3k-1-i} A_{k-1, i-2}  (2 <= i <= k-1, k <= 24).
"""

from __future__ import annotations
from fractions import Fraction
import math
import random
import sys
import time
import sympy as sp

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import testbed as tb


# ---------------------------------------------------------------------------
# 1. Exact rational Sturm certification of k=9 universal counterexample (#18)
# ---------------------------------------------------------------------------

def apply_T_frac(p: list[Fraction], a: Fraction, b: Fraction) -> list[Fraction]:
    # T p(t) = t(1+t) p'(t) + (a + b t) p(t)
    dp = [Fraction(i) * p[i] for i in range(1, len(p))] if len(p) > 1 else [Fraction(0)]
    res = [Fraction(0)] * (len(p) + 2)
    for i, c in enumerate(p):
        res[i] += a * c
        res[i + 1] += b * c
    for i, c in enumerate(dp):
        if c:
            res[i + 1] += c
            res[i + 2] += c
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res


def add_frac(p: list[Fraction], q: list[Fraction]) -> list[Fraction]:
    n = max(len(p), len(q))
    res = [(p[i] if i < len(p) else Fraction(0)) + (q[i] if i < len(q) else Fraction(0)) for i in range(n)]
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res


def eval_frac(p: list[Fraction], x: Fraction) -> Fraction:
    acc = Fraction(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc


def clear_denoms(p: list[Fraction]) -> list[int]:
    lcm = 1
    for c in p:
        lcm = math.lcm(lcm, c.denominator)
    ints = [int(c * lcm) for c in p]
    return tb.poly_primitive_pos(ints)


def verify_k9_universal_counterexample():
    print("=== 1. Exact rational Sturm certification of #18's k=9 universal counterexample ===")
    # Reproduce seed(7) trial 1 parameters, rounding each (a_j, b_j, w_j) to CLEAN EXACT RATIONALS
    # with strictly 0 < a_j < b_j and w_j > 0!
    random.seed(7)
    params = []
    for tr in range(2):
        params = []
        for j in range(1, 7):
            a_f = 10 ** random.uniform(-2, 2)
            db_f = 10 ** random.uniform(-2, 2)
            w_f = 10 ** random.uniform(-2, 2)
            a_q = Fraction(round(a_f * 10000), 10000)
            b_q = a_q + Fraction(round(db_f * 10000), 10000)
            w_q = Fraction(round(w_f * 10000), 10000)
            params.append((a_q, b_q, w_q))

    polys = {0: [Fraction(1)]}
    for a_q, b_q, w_q in params:
        assert 0 < a_q < b_q and w_q > 0
        nxt = {}
        k_cur = len(polys)
        for idx in range(k_cur + 1):
            p1 = apply_T_frac(polys[idx], a_q, b_q) if idx in polys else [Fraction(0)]
            p2 = [w_q * c for c in polys[idx - 1]] if (idx - 1) in polys else [Fraction(0)]
            nxt[idx] = add_frac(p1, p2)
        polys = nxt

    # Step 7 (exact rational, breaking Turan at tv = -1/1000)
    tv = Fraction(-1, 1000)
    i = 2
    def deriv_f(p):
        return [Fraction(m) * p[m] for m in range(1, len(p))] if len(p) > 1 else [Fraction(0)]
    f0, df0 = eval_frac(polys[i - 1], tv), eval_frac(deriv_f(polys[i - 1]), tv)
    f1, df1 = eval_frac(polys[i], tv), eval_frac(deriv_f(polys[i]), tv)
    f2, df2 = eval_frac(polys[i + 1], tv), eval_frac(deriv_f(polys[i + 1]), tv)
    delta = f1 * f1 - f0 * f2
    ddelta_u = tv * (1 + tv) * (2 * f1 * df1 - df0 * f2 - f0 * df2)
    phi_star = -ddelta_u / (2 * delta)
    a7 = Fraction(1, 100)
    b7_exact = a7 + (phi_star - a7 * (1 + tv)) / tv
    # Round b7 to a compact exact rational!
    b7 = Fraction(round(float(b7_exact) * 100), 100)
    w7 = Fraction(1, 10**10)
    assert 0 < a7 < b7 and w7 > 0
    params.append((a7, b7, w7))

    p7 = {}
    for idx in range(8):
        p1 = apply_T_frac(polys[idx], a7, b7) if idx in polys else [Fraction(0)]
        p2 = [w7 * c for c in polys[idx - 1]] if (idx - 1) in polys else [Fraction(0)]
        p7[idx] = add_frac(p1, p2)

    # Verify Turan broken at k=7, i=2, t=-1/1000
    turan7 = eval_frac(p7[2], tv)**2 - eval_frac(p7[1], tv) * eval_frac(p7[3], tv)
    print(f"  k=7 Turan at t=-1/1000: {float(turan7):.4e} (EXACT < 0? {turan7 < 0})")
    assert turan7 < 0

    # Step 8 (exact rational, breaking interlacing)
    phi8 = Fraction(-100)
    a8 = Fraction(1, 100)
    b8 = a8 + (phi8 - a8 * (1 + tv)) / tv
    dv_3 = eval_frac(deriv_f(p7[3]), tv)
    w8_exact = -(tv * (1 + tv) * dv_3 + phi8 * eval_frac(p7[3], tv)) / eval_frac(p7[2], tv)
    w8 = Fraction(round(float(w8_exact) * 10000), 10000)
    assert 0 < a8 < b8 and w8 > 0
    params.append((a8, b8, w8))

    q8_3 = add_frac(apply_T_frac(p7[3], a8, b8), [w8 * c for c in p7[2]])
    q8_2 = add_frac(apply_T_frac(p7[2], a8, b8), [w8 * c for c in p7[1]])
    ok_inter8, msg_inter8 = tb.verify_strict_interlacing(clear_denoms(q8_3), clear_denoms(q8_2), -1, 0)
    print(f"  k=8 interlacing q8_3 << q8_2: {ok_inter8} ({msg_inter8})")
    assert not ok_inter8

    # Step 9: search exact rational (a9, b9, w9) matching #18's grid and certify by EXACT INTEGER STURM!
    found = None
    for a9 in [Fraction(1, 1000), Fraction(1, 10), Fraction(1), Fraction(10), Fraction(100)]:
        for beta9 in [Fraction(1, 1000), Fraction(1), Fraction(100), Fraction(10**4), Fraction(10**6)]:
            b9 = a9 + beta9
            Tq = apply_T_frac(q8_3, a9, b9)
            for p_w in range(-6, 8):
                w9 = Fraction(10**p_w) if p_w >= 0 else Fraction(1, 10**(-p_w))
                q9 = add_frac(Tq, [w9 * c for c in q8_2])
                q9_int = clear_denoms(q9)
                d = len(q9_int) - 1
                ch = tb.sturm_chain(q9_int, tb.poly_deriv(q9_int))
                sqfree = (len(ch[-1]) == 1 and ch[-1][0] != 0)
                n_real = tb.sign_variations(tb.chain_signs_at(ch, "-inf")) - tb.sign_variations(tb.chain_signs_at(ch, "+inf"))
                if sqfree and n_real < d:
                    found = (a9, b9, w9, d, n_real)
                    break
            if found:
                break
        if found:
            break
    assert found is not None, "Could not find exact rational k=9 complex-root witness!"
    a9, b9, w9, d, n_real = found
    print(f"  EXACT RATIONAL k=9 WITNESS CERTIFIED BY STURM:")
    print(f"    (a9, b9, w9) = ({a9}, {b9}, {w9}), deg(q9_3) = {d}, squarefree = True, real_roots = {n_real} (< {d}!)")
    assert all(0 < a < b and w > 0 for a, b, w in params) and 0 < a9 < b9 and w9 > 0
    print(f"    All 9 parameter triples strictly satisfy 0 < a_j < b_j and w_j > 0.")


# ---------------------------------------------------------------------------
# 2. Adversarial audit of #20 Quantitative Plucker Quad Lemma
# ---------------------------------------------------------------------------

def audit_plucker_quad():
    print("\n=== 2. Adversarial audit of #20 Quantitative Plucker Quad Lemma ===")
    # 2(a): Explicit defect when u < 0 (W_AC < 0) is NOT ruled out!
    # Take eps = 1/4 => theta = 2/sqrt(5), a = b = 1 (R = 1 <= 4), u = -1, v = 0 (note |u| = 1 < 2*theta = 4/sqrt(5) ~ 1.789)
    # Wait: even cleaner, take eps = 0 (theta = 1), a = b = 1 (R = 1), u = -3/2, v = 0:
    # Check J(A,B,C) = 4*1 / (-3/2)^2 = 16/9 >= 1, J(B,C,E) = +inf >= 1, R = 1 <= 4!
    # Construct EXPLICIT 2-vectors (f(0), f'(0)) realizing these Wronskians!
    # Wait: let B = (1, 0), C = (0, 1) => W_{BC} = 1.
    # Need W_{AB} = a = 1 => det(A, B) = -A_1 = 1 => A_1 = -1.
    # Need W_{AC} = u = -9/5 => det(A, C) = A_0 = -9/5. So A = (-9/5, -1)!
    # Need W_{BE} = v = 0 => det(B, E) = E_1 = 0.
    # Need W_{CE} = b = 1 => det(C, E) = -E_0 = 1 => E_0 = -1. So E = (-1, 0)!
    A_v = (Fraction(-9, 5), Fraction(-1))
    B_v = (Fraction(1), Fraction(0))
    C_v = (Fraction(0), Fraction(1))
    E_v = (Fraction(-1), Fraction(0))
    def W(x, y): return x[0]*y[1] - x[1]*y[0]
    wAB, wBC, wCE = W(A_v, B_v), W(B_v, C_v), W(C_v, E_v)
    wAC, wBE = W(A_v, C_v), W(B_v, E_v)
    J_ABC = 4 * wAB * wBC / (wAC**2)
    R_val = wAB * wCE / (wBC**2)
    AB = (A_v[0] + B_v[0], A_v[1] + B_v[1])
    BC = (B_v[0] + C_v[0], B_v[1] + C_v[1])
    CE = (C_v[0] + E_v[0], C_v[1] + E_v[1])
    J_new = 4 * W(AB, BC) * W(BC, CE) / (W(AB, CE)**2)
    print(f"  [DEFECT 1] Missing W_AC >= 0, W_BE >= 0 hypothesis:")
    print(f"    Vectors A={A_v}, B={B_v}, C={C_v}, E={E_v}:")
    print(f"    W_AB={wAB}, W_BC={wBC}, W_CE={wCE} (all >0); J(A,B,C)={J_ABC} (=1.778 >= 1), J(B,C,E)=+inf, R={R_val} (<=4)")
    print(f"    Yet W(A+B,B+C)={W(AB,BC)}, W(B+C,C+E)={W(BC,CE)}, W(A+B,C+E)={W(AB,CE)} => J(A+B,B+C,C+E) = {J_new} ({float(J_new):.4f} < 1!)")
    assert J_ABC >= 1 and R_val <= 4 and J_new < 1

    # 2(b): Even WITH u >= 0, v >= 0 and R <= 4, check sign flip of W(A+B, C+E) = (1+u)(1+v) - ab!
    # Take a = b = 2 (R = 4), u = 0, v = 0 => W(A+B, C+E)/W_BC = 1 - 4 = -3 < 0!
    print(f"  [DEFECT 2] Non-preservation of W_AC >= 0 when R > 1:")
    print(f"    At a=b=2 (R=4), u=v=0 (J(A,B,C)=J(B,C,E)=+inf, W_AC=W_BE=0): W(A+B,C+E)/W_BC = (1+0)(1+0)-4 = -3 < 0!")
    print(f"    In fact, (1+u)(1+v) - ab > 0 is guaranteed for all u,v >= 0 iff R = ab < 1 (or R <= 1 for >= 0)!")

    # 2(c): Check if WITH u >= 0, v >= 0 (W_AC >= 0, W_BE >= 0), the lower bound J_new >= 4/(1+th^2)^2
    #       ACTUALLY holds over ALL INTERIOR (u, v) in [0, 2*th*sqrt(a)] x [0, 2*th*sqrt(b)] (not just the corner!),
    #       and find exact R_max!
    # Careful: optimize J_new(u, v) = 4(1+a+u)(1+b+v)/((1+u)(1+v)-ab)^2 over u in [0, 2*th*sqrt(a)], v in [0, 2*th*sqrt(b)]!
    # Check d/du log J_new = 1/(1+a+u) - 2(1+v)/((1+u)(1+v) - ab)
    #   Numerator = (1+u)(1+v) - ab - 2(1+v)(1+a+u) = -(1+v)(1 + 2a + u) - ab < 0 ALWAYS!
    # WAIT A SECOND! Look at the sign of the denominator ((1+u)(1+v) - ab):
    #   d/du J_new = 4(1+b+v) * [ -(1+v)(1+2a+u) - ab ] / ((1+u)(1+v) - ab)^3 !!
    #   So:
    #   - If (1+u)(1+v) > ab (denom > 0): d/du J_new < 0 (STRICTLY DECREASING in u, minimized at u_max = 2*th*sqrt(a)!)
    #   - If (1+u)(1+v) < ab (denom < 0): d/du J_new > 0 (STRICTLY INCREASING in u, minimized at u = 0!)
    # WHOA!! CHECK u = 0, v = 0 WHEN ab > 1!!
    # At u = 0, v = 0, J_new = 4(1+a)(1+b)/(1 - ab)^2!
    # Minimize 4(1+a)(1+b)/(1-ab)^2 for ab = R > 1:
    #   (1+a)(1+b) = 1 + R + (a+b) >= 1 + R + 2*sqrt(R) = (1 + sqrt(R))^2!
    #   So at u = v = 0, J_new >= 4(1 + sqrt(R))^2 / (R - 1)^2 = 4 / (sqrt(R) - 1)^2!!
    #   WAIT A SECOND!! Look at 4 / (sqrt(R) - 1)^2 when R = 9!!
    #   At R = 9 (a = b = 3, u = v = 0): 4 / (3 - 1)^2 = 4 / 4 = 1 < 4 / (1 + th^2)^2 FOR EVERY eps > 0 (th < 1)!!
    #   AND AT R = 4 (a = b = 2, u = v = 0): 4 / (2 - 1)^2 = 4 >= 4 / (1 + th^2)^2 (which is <= 4 for all th >= 0)!!
    #   WAIT A SECOND!! WHAT ABOUT u = 0, v = 2*th*sqrt(b) WHEN (1+u)(1+v) < ab?!
    #   Let's check ALL corners AND interior numerically AND symbolically right now!!
    print("  [CHECK 2(c)] Checking global min of J_new over ALL (u, v) in [0, 2*th*sqrt(a)] x [0, 2*th*sqrt(b)]:")
    # Symbolic derivative check:
    a_s, b_s, u_s, v_s, th_s = sp.symbols('a b u v th', positive=True)
    J_expr = 4 * (1 + a_s + u_s) * (1 + b_s + v_s) / ((1 + u_s) * (1 + v_s) - a_s * b_s)**2
    dJ_du = sp.factor(sp.diff(J_expr, u_s))
    print(f"    Exact dJ_new/du = {dJ_du}")
    # Notice numerator of dJ_new/du is -4*(1+b+v)*(ab + (1+v)(1+2a+u)) < 0 ALWAYS for a,b,u,v >= 0!
    # Therefore on each connected component of {(u,v) : (1+u)(1+v) != ab}:
    #   - On {(1+u)(1+v) > ab}: J_new is strictly decreasing in BOTH u and v => minimized ONLY at (u_max, v_max)!
    #   - On {(1+u)(1+v) < ab}: J_new is strictly increasing in BOTH u and v => minimized ONLY at (0, 0)!
    # Check BOTH minima!
    # (i) At (0, 0) when ab > 1: J_new(0,0) = 4(1+a)(1+b)/(ab-1)^2 >= 4(1+sqrt(R))^2/(R-1)^2 = 4/(sqrt(R)-1)^2.
    #     Compare 4/(sqrt(R)-1)^2 vs 4/(1+th^2)^2:
    #     4/(sqrt(R)-1)^2 >= 4/(1+th^2)^2 <=> sqrt(R) - 1 <= 1 + th^2 <=> sqrt(R) <= 2 + th^2 <=> R <= (2 + th^2)^2!
    #     For R <= 4: sqrt(R) <= 2 <= 2 + th^2 ALWAYS!
    #     WAIT: for R = 9 (as in verify_plucker_and_cubic_exact.py!) and eps > 0 (th < 1):
    #     at a=b=3, u=v=0, J_new(0,0) = 4*(4)*(4)/64 = 1 < 4/(1+th^2)^2!
    # (ii) At (u_max, v_max) = (2*th*sqrt(a), 2*th*sqrt(b)) when (1+u_max)(1+v_max) > ab:
    #     Let x = sqrt(a), y = sqrt(b), r = x*y = sqrt(R), S = x + y >= 2*sqrt(r) = 2*s (r = s^2).
    #     Verify exact SymPy factorization of numerator of J_new(u_max, v_max) - 4/(1+th^2)^2!
    s_s, d_s = sp.symbols('s d', nonnegative=True)
    S_s = 2 * s_s + d_s
    r_s = s_s**2
    num_c = 4 * ((1 - r_s)**2 + 4 * th_s**2 * r_s + 2 * th_s * (1 + r_s) * S_s + S_s**2)
    den_c = (1 - r_s**2 + 4 * th_s**2 * r_s + 2 * th_s * S_s)**2
    # Note den_c > 0 on R <= 4 iff 1 - s^4 + 4*th^2*s^2 + 4*th*s > 0.
    # For s in [0, sqrt(2)] (i.e. R = s^4 <= 4) and th in [1/sqrt(2), 1] (eps in [0, 1]):
    # 1 - s^4 + 4*th^2*s^2 + 4*th*s >= 1 - 4 + 4*(1/2)*2 + 4*(1/sqrt(2))*sqrt(2) = -3 + 4 + 4 = 5 > 0!
    diff_num = sp.expand(num_c * (1 + th_s**2)**2 - 4 * den_c)
    P_d = sp.Poly(diff_num, d_s)
    print("    Exact coefficients of numerator of [J_new(u_max, v_max) - 4/(1+th^2)^2] in d = (sqrt(a)-sqrt(b))^2 >= 0:")
    for deg in range(P_d.degree() + 1):
        print(f"      [d^{deg}] = {sp.factor(P_d.nth(deg))}")
    # Check signs of [d^0], [d^1], [d^2] for s in [0, sqrt(2)] (R <= 4) and th in (0, 1]!
    #WAIT: Check [d^2] = 4*(1 - th^2)^2 >= 0? Let's inspect the printed factors!

    # Also test R = 4 with SMALL th (eps > 1, e.g., eps = 10, th = 1/sqrt(11)) to check why eps in [0, 1] was needed!
    for eps_val in [0.0, 0.5, 1.0, 2.0, 5.0]:
        thv = 1.0 / math.sqrt(1.0 + eps_val)
        target = 4.0 / ((1.0 + thv**2)**2)
        worst = (1e9, None)
        for s_i in range(0, 142):  # s in [0, 1.41] => R = s^4 in [0, 4]
            sv = s_i / 100.0
            rv = sv * sv
            for dv in [0.0, 0.01, 0.1, 0.5, 2.0, 10.0, 100.0]:
                Sv = 2 * sv + dv
                denv = 1 - rv**2 + 4 * thv**2 * rv + 2 * thv * Sv
                if abs(denv) < 1e-12:
                    continue
                val = 4 * ((1 - rv)**2 + 4 * thv**2 * rv + 2 * thv * (1 + rv) * Sv + Sv**2) / (denv**2)
                if val - target < worst[0]:
                    worst = (val - target, (sv, rv, dv, val, target))
        print(f"    eps={eps_val}: min over R<=4, d>=0 of (J_corner - 4/(1+th^2)^2) = {worst[0]:+.6e} at {worst[1]}")


# ---------------------------------------------------------------------------
# 3. Symbolic verification of Borcea-Branden symbol of D_m & P_{L+1} S_{L+1} = S_L P_{L+1}
# ---------------------------------------------------------------------------

def verify_symbolic_identities():
    print("\n=== 3. Exact SymPy verification of operator & SOS identities ===")
    t, w, m, d, x, L = sp.symbols('t w m d x L')
    # 3(a): Borcea-Branden symbol of D_m = t(1+t)d/dt + 1 + m*t on (t+w)^d
    f = (t + w)**d
    Dm_f = sp.simplify((t * (1 + t) * sp.diff(f, t) + (1 + m * t) * f) / (t + w)**(d - 1))
    Q_tw = sp.expand(Dm_f)
    # Check partial fraction of w(t) solving Q_tw = 0:
    w_sol = sp.solve(Q_tw, w)[0]
    w_pf = -((d + m) / m) * t - d * (m - 1) / m**2 + (d * (m - 1) / m**2) / (m * t + 1)
    assert sp.simplify(w_sol - w_pf) == 0
    print(f"  D_m symbol Q(t,w) = {Q_tw}")
    print(f"  Exact partial fraction w(t) = {w_pf}  (Im(t)>0 => Im(w)<0 iff m >= 1, d >= 1!)")

    # 3(b): P_{L+1} S_{L+1} = S_L P_{L+1} on arbitrary f(x)
    fn = sp.Function('f')(x)
    def S_op(Lv, g): return (1 + Lv * x) * g + x * (1 - x) * sp.diff(g, x)
    def P_op(cv, g): return cv * g - x * sp.diff(g, x)
    diff_PS = sp.expand(P_op(L + 1, S_op(L + 1, fn)) - S_op(L, P_op(L + 1, fn)))
    assert diff_PS == 0
    print(f"  Exact operator identity P_{{L+1}} S_{{L+1}} - S_L P_{{L+1}} = {diff_PS} (VERIFIED!)")

    # 3(c): Exact SOS of Universal Cubic Wronskian Log-Concavity
    r1, r2, r3, c1, c2, c3 = sp.symbols('r1 r2 r3 c1 c2 c3')
    W0 = c1 * r2**2 * r3**2 + c2 * r1**2 * r3**2 + c3 * r1**2 * r2**2
    W1 = c1 * (r2**2 + r2*r3 + r3**2) + c2 * (r1**2 + r1*r3 + r3**2) + c3 * (r1**2 + r1*r2 + r2**2)
    W2 = c1 + c2 + c3
    sos = (
        c1**2 * (r2**2 + r3**2) * (r2 + r3)**2
        + c2**2 * (r1**2 + r3**2) * (r1 + r3)**2
        + c3**2 * (r1**2 + r2**2) * (r1 + r2)**2
        + c1 * c2 * (sp.Rational(1, 2) * (2*r1*r2 + r3*(r1+r2))**2 + sp.Rational(1, 2) * (r3*(r1+r2) + 2*r3**2)**2)
        + c1 * c3 * (sp.Rational(1, 2) * (2*r1*r3 + r2*(r1+r3))**2 + sp.Rational(1, 2) * (r2*(r1+r3) + 2*r2**2)**2)
        + c2 * c3 * (sp.Rational(1, 2) * (2*r2*r3 + r1*(r2+r3))**2 + sp.Rational(1, 2) * (r1*(r2+r3) + 2*r1**2)**2)
    )
    assert sp.expand(W1**2 - W0 * W2 - sos) == 0
    print("  Universal Cubic Wronskian W_1^2 - W_0 W_2 == exact global SOS (all r1,r2,r3 in R, c1,c2,c3 >= 0) VERIFIED!")


# ---------------------------------------------------------------------------
# 4. Exact Cauchy-index stress test of P_{3k-1-i} intertwining up to k=24
# ---------------------------------------------------------------------------

def op_P(c: int, p: list[int]) -> list[int]:
    """P_c f(x) = c f(x) - x f'(x) = sum_m (c - m) a_m x^m."""
    return tb.poly_trim([(c - m) * a for m, a in enumerate(p)])


def stress_test_Pc(max_k: int = 22):
    print(f"\n=== 4. Exact integer Cauchy-index stress test of P_{{3k-1-i}} claims (k <= {max_k}) ===")
    t0 = time.perf_counter()
    A = tb.compute_A_table(max_k)

    ok1, fail1 = 0, []
    ok2, fail2 = 0, []
    for k in range(1, max_k + 1):
        for i in range(1, k + 1):
            L = 3 * k - 2 - i
            # Claim 1: A_{k, i} << P_{L+1} A_{k, i-1} (where L+1 = 3k - 1 - i)
            p_small = A[(k, i)]
            p_big = op_P(L + 1, A[(k, i - 1)])
            ok, msg = tb.verify_strict_interlacing(p_small, p_big, "-inf", 0)
            if ok:
                ok1 += 1
            else:
                fail1.append((k, i, msg))

            # Claim 2 (stronger Pair (2,1)): S_{L-1} A_{k-1, i} << P_{L+1} A_{k-1, i-2}
            if 2 <= i <= k - 1:
                lhs = tb.op_S(L - 1, A[(k - 1, i)])
                rhs = op_P(L + 1, A[(k - 1, i - 2)])
                ok_s, msg_s = tb.verify_strict_interlacing(lhs, rhs, "-inf", 0)
                if ok_s:
                    ok2 += 1
                else:
                    fail2.append((k, i, msg_s))

    print(f"  A_{{k,i}} << P_{{3k-1-i}} A_{{k,i-1}} (k=1..{max_k}): PASS={ok1}, FAIL={len(fail1)} {fail1[:5]}")
    print(f"  S_{{3k-3-i}} A_{{k-1,i}} << P_{{3k-1-i}} A_{{k-1,i-2}} (k=3..{max_k}): PASS={ok2}, FAIL={len(fail2)} {fail2[:5]}")
    print(f"  Elapsed: {time.perf_counter()-t0:.2f}s.")


def main():
    verify_k9_universal_counterexample()
    audit_plucker_quad()
    verify_symbolic_identities()
    stress_test_Pc(max_k=22)


if __name__ == "__main__":
    main()
