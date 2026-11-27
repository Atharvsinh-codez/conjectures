"""
1. Verify EVERY step of our 100% self-contained 12-line proof of c = 1 for m = 2..14:
   - Exact Rolle recurrence F_r(u) = n u(1-u) F'_{r-1}(u) + r(1 + (n-1)u) F_{r-1}(u) gives F_{n-1}(u)/(n-1)! == A_{n-1}^{(n)}(u)!
   - Exact formula gamma_{m-1} == (1/(2m)) sum_{k=0}^{m-1} (-1)^k cos(theta_k) / sin(theta_k)^{2m+1} (theta_k = (2k+1)pi/(4m))!
   - Verify alternating series lower bound > 4^{m-1} for ALL m >= 2!
2. Derive & verify exact closed form of c = 2 slice Q_2^{(2m)}(x+y, xy) from log p_{2m}(x, y, z) = 2m log(x-y) - sum z^r/r S_r(x, y)!
   Check whether Q_2 also admits a clean algebraic positivity proof (gamma-alternating / real-rooted in (0, 1/4))!
3. Check the exact formula for the e2 = 1 Vandermonde-Cauchy chart to share with @analyst!
"""
import math
import sympy as sp
import mpmath as mp

mp.mp.dps = 60

print("=== 1. VERIFY SELF-CONTAINED c = 1 PROOF (m = 2..12) ===")
for m in range(2, 13):
    n = 2 * m
    # Build F_r via exact recurrence F_0 = [1]
    F = [1]
    for r in range(1, n):
        # F_r(u) = n u(1-u) F'(u) + r(1 + (n-1)u) F(u)
        Fn = [0] * (len(F) + 1)
        for k, c in enumerate(F):
            # from n u(1-u) k c u^{k-1} = n k c (u^k - u^{k+1})
            Fn[k] += n * k * c
            Fn[k + 1] -= n * k * c
            # from r(1 + (n-1)u) c u^k = r c u^k + r(n-1) c u^{k+1}
            Fn[k] += r * c
            Fn[k + 1] += r * (n - 1) * c
        F = Fn
    fac = math.factorial(n - 1)
    assert all(v % fac == 0 for v in F)
    A_rec = [v // fac for v in F]
    A_dir = [sum(((-1)**r) * math.comb(n, r) * math.comb(n*(j - r + 1) - 1, n - 1) for r in range(j + 1)) for j in range(n)]
    assert A_rec == A_dir, f"Mismatch at m={m}!"

    # Compute exact gamma_{m-1} from B(u) = A(u)/(1+u) at u = -1: (-1)^{m-1} B(-1)
    # Equivalently (-1)^{m-1} A'(-1) = (-1)^{m-1} sum_{j=1}^{n-1} j (-1)^{j-1} A[j]!
    gam_exact = ((-1)**(m - 1)) * sum(j * ((-1)**(j - 1)) * A_dir[j] for j in range(1, n))

    # Compare with exact trig sum S_trig = (1/(2m)) sum_{k=0}^{m-1} (-1)^k cos(th_k) / sin(th_k)^{2m+1}
    S_trig = mp.mpf(0)
    t0_val = mp.mpf(0); t1_val = mp.mpf(0)
    for k in range(m):
        th = (2*k + 1) * mp.pi / (4 * m)
        term = mp.cos(th) / (mp.sin(th)**(2*m + 1))
        if k == 0: t0_val = term / (2*m)
        if k == 1: t1_val = term / (2*m)
        S_trig += ((-1)**k) * term
    S_trig /= (2 * m)
    rel_err = abs(S_trig - gam_exact) / abs(gam_exact)
    ratio_4 = float((t0_val - t1_val) / (4**(m - 1)))
    print(f"  m={m:2d}: F_{n-1}/(n-1)! == A OK, gam_{{m-1}}={gam_exact}, trig_rel_err={float(rel_err):.1e}, (t0-t1)/4^{{m-1}}={ratio_4:.2f}")

print("\n=== 2. EXACT FORMULA FOR c = 2 SLICE Q_2(x+y, xy) VIA VERONESE H_1, H_2 ===")
# Expand p_{2m}(x, y, z) = (x-y)^{2m} exp(-2m z H_1/(x-y)^{2m} - m z^2 H_2/(x-y)^{4m} + O(z^3))
# Coeff of z^2 in p_{2m}(x, y, z):
#   [z^2] p_{2m} = (1/(x-y)^{2m}) * [ 2 m^2 H_1(x,y)^2 - m H_2(x,y) ]!
# Wait: check at z^2 in P_{2m}(e1, e2, e3) = (e1^2 - 4 e2)^m + e3 Q_1(e1, e2) + e3^2 Q_2(e1, e2) + O(e3^3)!
# At e1 = s + z, e2 = p + s*z, e3 = p*z:
#   (e1^2 - 4 e2)^m = ((s^2 - 4p) - 2 s z + z^2)^m
#     z^2 coeff = m(s^2 - 4p)^{m-1} + 2 m(m-1) s^2 (s^2 - 4p)^{m-2}
#   e3 Q_1(e1, e2) = p z [ Q_1(s, p) + z (d/ds Q_1 + s d/dp Q_1) + O(z^2) ]
#     z^2 coeff = p (Q_{1,s}(s, p) + s Q_{1,p}(s, p))
#   e3^2 Q_2(e1, e2) = p^2 z^2 Q_2(s, p) + O(z^3)!
# Let's verify this exact c=2 formula at m = 3, 4, 5!
x, y = sp.symbols('x y')
for m in [3, 4, 5]:
    n = 2 * m
    # H_r(x, y) = sum_{j=0}^{2mr - r} [t^{2mj}] ((1-t^{2m})/(1-t))^{2mr} x^{2mr-r-j} y^j
    # Wait: (1-u)^{2mr} sum_{s>=0} C(2m(s+r)-1, 2mr-1) u^s:
    # Coeff of u^j is sum_{l=0}^j (-1)^l C(2mr, l) C(2m(j-l+r)-1, 2mr-1)
    def H_poly(r):
        deg = 2*m*r - r
        res = 0
        for j in range(deg + 1):
            c = sum(((-1)**l) * math.comb(2*m*r, l) * math.comb(2*m*(j - l + r) - 1, 2*m*r - 1) for l in range(j + 1))
            res += c * x**(deg - j) * y**j
        return sp.expand(res)
    H1 = H_poly(1)
    H2 = H_poly(2)
    z2_lhs = sp.cancel((2*m*m * H1**2 - m * H2) / (x - y)**(2*m))
    # Verify (x-y)^{2m} divides 2m^2 H1^2 - m H2!
    z2_lhs = sp.Poly(z2_lhs, x, y).as_expr()
    # Subtract c=0 z^2 term and c=1 z^2 term!
    s_xy, p_xy = x + y, x * y
    D_xy = (x - y)**2
    z2_c0 = sp.expand(m * D_xy**(m - 1) + 2 * m * (m - 1) * s_xy**2 * D_xy**(m - 2))
    # Compute Q1(s, p)
    Q1_xy = sp.cancel(-2*m * (H1 - s_xy * D_xy**(m - 1)) / p_xy)
    # Express Q1 in (s, p) to differentiate wrt s and p!
    # Wait: d/ds Q1(s, p) + s d/dp Q1(s, p) = ?
    # By chain rule: d/dx = d/ds + y d/dp, d/dy = d/ds + x d/dp!
    # So x d/dx + y d/dy = s d/ds + 2 p d/dp (homogeneity = deg * Q1),
    # and d/dx + d/dy = 2 d/ds + s d/dp!
    # Wait: what about d/ds + s d/dp?!
    # Notice: (x d/dy + y d/dx) = (x+y) d/ds + (x^2+y^2) d/dp = s d/ds + (s^2-2p) d/dp...
    # Simpler: in (x, y), keeping z->0: e1 = x+y+z, e2 = xy+(x+y)z, e3 = xyz.
    # Derivative of Q1(x+y+z, xy+(x+y)z) wrt z at z=0 is Q1_s + (x+y) Q1_p = ((x-y)Q1_x ...):
    # Solve: Q1_x = Q1_s + y Q1_p, Q1_y = Q1_s + x Q1_p => Q1_p = (Q1_y - Q1_x)/(x - y), Q1_s = (x Q1_x - y Q1_y)/(x - y)!
    # Hence Q1_s + (x+y) Q1_p = (x Q1_x - y Q1_y + (x+y)(Q1_y - Q1_x)) / (x - y) = (x Q1_y - y Q1_x) / (x - y)!!
    dQ1_x = sp.diff(Q1_xy, x); dQ1_y = sp.diff(Q1_xy, y)
    z2_c1 = sp.cancel(p_xy * (x * dQ1_y - y * dQ1_x) / (x - y))
    Q2_xy = sp.cancel((z2_lhs - z2_c0 - z2_c1) / (p_xy**2))
    # Evaluate Q2 at x=1, y=0 (s=1, p=0 -> coeff c_{n-6, 0, 2}!)
    c_lead = sp.simplify(Q2_xy.subs({x: 1, y: 0}))
    print(f"  m={m} (n={n}): Q2(1, 0) = c_{{{n-6},0,2}} = {c_lead}")
    # Also check Q2(x, y) factored or at (x=1, y=u)
    Q2_u = sp.Poly(sp.expand(Q2_xy.subs(x, 1)), y)
    print(f"    deg Q2(1, y) = {Q2_u.degree()}, Q2(1, 1) = {Q2_u.eval(1)}, Q2(1, -1) = {Q2_u.eval(-1)}")
