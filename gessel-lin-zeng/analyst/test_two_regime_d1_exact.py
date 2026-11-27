"""
Decisive test of the Two-Regime Split closure (Coordinator #33 item 1):
1. Exact d=1, d=2 analytic + 100-dps test of whether D_m preserves:
   (H1) J_raw(g,f,h) >= 1 on (-1, -1/m]
   (H2) A > 0, C0 > 0 on (-1, 0], and B >= 0 wherever J_raw < 1 on [-1/m, 0]
   Check BOTH output on (-1, -1/m] AND output on (-1, -1/(m+1)]!
2. Test the FULL GLZ step (first D_m, addition F_i = D_m h_i + h_{i-1}, second D_m)
   on 4-tuples satisfying (H1)+(H2): does it preserve (H1)+(H2) AND Plain Turan?
"""
import sympy as sp
import mpmath as mp
import numpy as np

mp.mp.dps = 100

def apply_Dm_mp(m, c):
    out = [mp.mpf(0)] * (len(c) + 1)
    for j, cj in enumerate(c):
        out[j] += (j + 1) * cj
        out[j + 1] += (j + m) * cj
    return out

def eval_mp(c, t):
    v = d = mp.mpf(0)
    for coeff in reversed(c):
        d = d*t + v
        v = v*t + coeff
    return v, d

def poly_from_roots_mp(roots):
    c = [mp.mpf(1)]
    for r in roots:
        inv_r = -mp.mpf(1) / r
        new_c = [mp.mpf(0)] * (len(c) + 1)
        for j, cj in enumerate(c):
            new_c[j] += cj
            new_c[j+1] += cj * inv_r
        c = new_c
    return c

def check_exact_d1():
    # d=1: g(t) = (1+a*t)(1+b*t), f(t) = c*(1+r*t), h(t) = 1, with 1 < a < r < b
    # Note: roots are -1/a, -1/b in (-1, 0) since a, b > 1.
    # J_raw(t) = 4*c^2*r*((a+b-r) + 2*a*b*t + a*b*r*t^2) / (a+b+2*a*b*t)^2.
    # For J_raw(t) to be < 1 near t=0 and >= 1 on (-1, -1/m]:
    # Notice J_raw(t) is STRICTLY DECREASING on (-1/(2a)...? Wait: W_gh = -(a+b+2abt) vanishes at t0 = -(a+b)/(2ab)!
    # Let's test explicit rational parameters!
    t = sp.Symbol('t')
    def D(m, p): return sp.expand(t*(1+t)*sp.diff(p, t) + (1 + m*t)*p)
    print("=== PART 1: Exact SymPy d=1 test of D_m on Two-Regime condition ===")
    cases = [
        # m, a, r, b, c (exact Rationals)
        (10, sp.Rational(2), sp.Rational(12), sp.Rational(100), sp.Rational(11, 10)),
        (10, sp.Rational(3), sp.Rational(15), sp.Rational(200), sp.Rational(11, 10)),
        (10, sp.Rational(5), sp.Rational(12), sp.Rational(500), sp.Rational(11, 10)),
        (10, sp.Rational(8), sp.Rational(12), sp.Rational(1000), sp.Rational(11, 10)),
        (20, sp.Rational(15), sp.Rational(25), sp.Rational(2000), sp.Rational(105, 100)),
        (5,  sp.Rational(4), sp.Rational(6),  sp.Rational(500), sp.Rational(105, 100)),
        (5,  sp.Rational(3,2), sp.Rational(5), sp.Rational(100), sp.Rational(105, 100)),
    ]
    for m, a, r, b, c in cases:
        g = sp.expand((1 + a*t)*(1 + b*t))
        f = sp.expand(c*(1 + r*t))
        h = sp.Integer(1)
        def min_J_and_check(G, F, H, t_hi):
            Wgf = sp.expand(G*sp.diff(F,t) - sp.diff(G,t)*F)
            Wfh = sp.expand(F*sp.diff(H,t) - sp.diff(F,t)*H)
            Wgh = sp.expand(G*sp.diff(H,t) - sp.diff(G,t)*H)
            # Check exact Sturm roots of num = 4*Wgf*Wfh - Wgh^2 on (-1, t_hi)!
            num = sp.Poly(sp.expand(4*Wgf*Wfh - Wgh**2), t)
            A = sp.Poly(sp.expand(F**2 - G*H), t)
            n_neg_J = sp.count_roots(num, -1, t_hi)
            # Evaluate at endpoints and midpoint
            val_hi = sp.N(num.eval(t_hi))
            val_m1 = sp.N(num.eval(-1))
            val_0  = sp.N(num.eval(0))
            n_A = sp.count_roots(A, -1, 0)
            return n_neg_J, float(val_m1), float(val_hi), float(val_0), n_A
        # First scale c so num(-1/m) > 0 and num(0) < 0
        # Compute c_crit at t = -1/m
        Wgf0 = sp.expand(g*sp.diff(1+r*t,t) - sp.diff(g,t)*(1+r*t))
        Wfh0 = sp.expand((1+r*t)*0 - r*1)
        Wgh0 = sp.expand(-sp.diff(g,t))
        tm = sp.Rational(-1, m)
        ratio_tm = sp.simplify(Wgh0.subs(t, tm)**2 / (4 * Wgf0.subs(t, tm) * Wfh0))
        if ratio_tm <= 0:
            continue
        # Set c^2 = ratio_tm * (1 + 1/100) so J_raw(-1/m) = 1.01
        c_val = sp.Rational(int(mp.ceil(mp.sqrt(float(ratio_tm)*1.01) * 1000)), 1000)
        f = sp.expand(c_val * (1 + r*t))
        in_info = min_J_and_check(g, f, h, tm)
        # Check B >= 0 where J < 1 on [tm, 0]: at t=0, A'(0) = 2*c^2*r - (a+b)
        Ap0 = float(2*c_val**2*r - (a+b))
        # Apply D_m
        gm, fm, hm = D(m, g), D(m, f), D(m, h)
        out_same = min_J_and_check(gm, fm, hm, tm)
        out_next = min_J_and_check(gm, fm, hm, sp.Rational(-1, m+1))
        # Also check exact min of J_raw(gm,fm,hm) on (-1, tm] and (-1, -1/(m+1)]
        print(f"m={m}, a={a}, r={r}, b={b}, c={float(c_val):.3f}, A'(0)={Ap0:.1f}:")
        print(f"  IN  (-1,-1/m]: roots(4WgfWfh-Wgh^2)={in_info[0]}, val(-1)={in_info[1]:.2e}, val(-1/m)={in_info[2]:.2e}, val(0)={in_info[3]:.2e}, roots(A)={in_info[4]}")
        print(f"  OUT (-1,-1/m]: roots={out_same[0]}, val(-1)={out_same[1]:.2e}, val(-1/m)={out_same[2]:.2e}, roots(A)={out_same[4]}")
        print(f"  OUT (-1,-1/(m+1)]: roots={out_next[0]}, val(-1/(m+1))={out_next[2]:.2e}")

check_exact_d1()

print("\n=== PART 2: Exact high-precision scan of d=1,2,3 for D_m dropping J_raw on (-1,-1/m] or (-1,-1/(m+1)] ===")
# Systematic sweep over (m, a, r, b) for d=1 where J_raw(-1/m) = 1.001, A > 0, A'(0) < 0!
t = sp.Symbol('t')
def D_sp(m, p): return sp.expand(t*(1+t)*sp.diff(p, t) + (1 + m*t)*p)

found_fail_same = None
found_fail_next = None
for m in [4, 10, 20, 50]:
    tm = -1.0 / m
    tm1 = -1.0 / (m + 1)
    for a in [1.1, 2.0, m*0.5, m*0.9, m*1.1, m*2.0]:
        for r_fac in [1.1, 1.5, 2.5, 5.0]:
            r = a * r_fac
            for b_fac in [1.5, 3.0, 10.0, 50.0]:
                b = r * b_fac
                # Exact mpmath computation
                ma, mr, mb = mp.mpf(a), mp.mpf(r), mp.mpf(b)
                g = [mp.mpf(1), ma + mb, ma*mb]
                f0 = [mp.mpf(1), mr]
                h = [mp.mpf(1)]
                # Check min J_raw(g, f0, h) on (-1, -1/m] (where c=1)
                grid_I = [-1 + mp.mpf(j)*(mp.mpf(-1)/m + 1)/60 for j in range(1, 61)]
                min_J0 = mp.inf
                for tv in grid_I:
                    gv, gd = eval_mp(g, tv); fv, fd = eval_mp(f0, tv); hv, hd = eval_mp(h, tv)
                    Wgf = gv*fd - gd*fv; Wfh = fv*hd - fd*hv; Wgh = gv*hd - gd*hv
                    J0 = 4*Wgf*Wfh / (Wgh*Wgh) if Wgh != 0 else mp.inf
                    if J0 < min_J0: min_J0 = J0
                if min_J0 <= 0: continue
                c = mp.sqrt(mp.mpf('1.0001') / min_J0)
                f = [c * cj for cj in f0]
                # Check A > 0 on (-1, 0] and A'(0) <= 0 if J(0) < 1
                A_min = min(eval_mp(f, tv)[0]**2 - eval_mp(g, tv)[0] for tv in grid_I + [mp.mpf(0), mp.mpf(-1)/(2*m)])
                if A_min <= 0: continue
                # Apply D_m
                gm, fm, hm = apply_Dm_mp(m, g), apply_Dm_mp(m, f), apply_Dm_mp(m, h)
                # Evaluate at -1/m and -1/(m+1) and grid
                min_same = mp.inf
                for tv in grid_I:
                    gv, gd = eval_mp(gm, tv); fv, fd = eval_mp(fm, tv); hv, hd = eval_mp(hm, tv)
                    Wgf = gv*fd - gd*fv; Wfh = fv*hd - fd*hv; Wgh = gv*hd - gd*hv
                    J = 4*Wgf*Wfh / (Wgh*Wgh) if Wgh != 0 else mp.inf
                    if J < min_same: min_same = J
                grid_next = [mp.mpf(-1)/m + mp.mpf(j)*(mp.mpf(-1)/(m+1) - mp.mpf(-1)/m)/20 for j in range(1, 21)]
                min_next = min_same
                for tv in grid_next:
                    gv, gd = eval_mp(gm, tv); fv, fd = eval_mp(fm, tv); hv, hd = eval_mp(hm, tv)
                    Wgf = gv*fd - gd*fv; Wfh = fv*hd - fd*hv; Wgh = gv*hd - gd*hv
                    J = 4*Wgf*Wfh / (Wgh*Wgh) if Wgh != 0 else mp.inf
                    if J < min_next: min_next = J
                # Check if A'(0) <= 0 (so Two-Regime condition holds!)
                Ap0 = float(2*c*c*mr - (ma+mb))
                # Check J_raw(0)
                J0_val = float(4*c*c*mr*(ma+mb-mr)/(ma+mb)**2)
                if min_same < 1 and (found_fail_same is None or min_same < found_fail_same[0]):
                    found_fail_same = (float(min_same), m, a, r, b, float(c), J0_val, Ap0)
                if min_next < 1 and (found_fail_next is None or min_next < found_fail_next[0]):
                    found_fail_next = (float(min_next), m, a, r, b, float(c), J0_val, Ap0)

print("Worst d=1 fail on (-1, -1/m]:", found_fail_same)
print("Worst d=1 fail on (-1, -1/(m+1)]:", found_fail_next)
