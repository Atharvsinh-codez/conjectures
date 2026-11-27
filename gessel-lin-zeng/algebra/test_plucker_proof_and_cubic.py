"""
1. Symbolic SymPy proof of the EXACT Stage 1 Plucker inequality:
   (1 + u^2)^2 (1 + 2 u p + p^2)(1 + 2 u q + q^2) - [(1 + 2 u p)(1 + 2 u q) - p^2 q^2]^2 >= 0!
   Find its EXACT factorization / SOS / discriminant in s = p + q, and exact threshold on x = p*q!
2. Check what real stability of the CUBIC C(t, y) = A - 3 B y + 3 C y^2 - E y^3 implies about
   (a) 4 W_AB W_BC / W_AC^2  and  (b) W_AB W_CE / W_BC^2!
3. Check if C_{k,i}(t, y) = A - 3 B y + 3 C y^2 - E y^3 (for A=h_{k-1,i-3}, B=D h_{k-1,i-2}, C=D^2 h_{k-1,i-1}, E=D^3 h_{k-1,i})
   is real-stable for GLZ up to k=20 (UHP test & Rayleigh test)!
4. Check if A - 3*c*B*y + 3*c*C*y^2 - E*y^3 (or without the 3's!) is real-stable!
"""
import math
import sympy as sp
import mpmath as mp
mp.mp.dps = 50

print("=== 1. SYMBOLIC PROOF OF STAGE 1 PLUCKER INEQUALITY ===", flush=True)
u, p, q, s, x, d = sp.symbols('u p q s x d', positive=True)
# Expression E(u, s, x) = (1+u^2)^2 * [(1+x)^2 + 2*u*(1+x)*s + s^2 - 4*(1-u^2)*x] - [1 + 2*u*s + 4*u^2*x - x^2]^2
W1W2 = (1 + x)**2 + 2*u*(1 + x)*s + s**2 - 4*(1 - u**2)*x
W3 = 1 + 2*u*s + 4*u**2*x - x**2
E_expr = sp.expand((1 + u**2)**2 * W1W2 - W3**2)
poly_s = sp.Poly(E_expr, s)
c2 = sp.factor(poly_s.nth(2))
c1 = sp.factor(poly_s.nth(1))
c0 = sp.factor(poly_s.nth(0))
print(f"  Quadratic in s = p + q:")
print(f"    c2 = {c2}")
print(f"    c1 = {c1}")
print(f"    c0 = {c0}")
disc_s = sp.factor(sp.expand(c1**2 - 4*c2*c0))
print(f"    disc_s = c1^2 - 4*c2*c0 = {disc_s}")

# Also substitute s = 2*r + d where r = sqrt(x) (so x = r^2, d >= 0)!
r = sp.Symbol('r', positive=True)
E_rd = sp.expand(E_expr.subs({x: r**2, s: 2*r + d}))
poly_d = sp.Poly(E_rd, d)
print(f"\n  Quadratic in d = (p + q) - 2*sqrt(p*q) >= 0 (with r = sqrt(p*q)):")
print(f"    d^2 coeff = {sp.factor(poly_d.nth(2))}")
print(f"    d^1 coeff = {sp.factor(poly_d.nth(1))}")
print(f"    d^0 coeff = {sp.factor(poly_d.nth(0))}")

# Check exact sign of d^1 coeff and disc_s!
# Wait: check numerically where E(u, 2*r + d, r^2) can EVER be negative for u in (0, 1], r > 0, d >= 0!
print("\n  Check exact min of E_expr over (u in (0,1], r in (0, 3], d >= 0):")
import random
random.seed(1)
for rmax in [1.0, 1.5, 2.0, 2.5, 3.0]:
    mn = 1e100
    worst = None
    for _ in range(100000):
        uv = random.uniform(0.01, 1.0)
        rv = random.uniform(0.01, rmax)
        # Check vertex in d if d_vert > 0!
        xv = rv*rv
        a_c = (1 - uv*uv)**2
        b_c = 2*uv*((1+uv*uv)**2*(1+xv) - 2*(1 + 4*uv*uv*xv - xv*xv))
        # min on s >= 2*rv:
        s_min = 2*rv
        if a_c > 1e-12:
            sv = -b_c / (2*a_c)
            if sv > s_min: s_min = sv
        val = (1+uv*uv)**2*((1+xv)**2 + 2*uv*(1+xv)*s_min + s_min**2 - 4*(1-uv*uv)*xv) - (1 + 2*uv*s_min + 4*uv*uv*xv - xv*xv)**2
        if val < mn:
            mn = val
            worst = (uv, rv, s_min)
    print(f"    r = sqrt(pq) <= {rmax}: min E = {mn:+.6e} at u={worst[0]:.4f}, r={worst[1]:.4f}, s={worst[2]:.4f}", flush=True)

print("\n=== 2. CHECK CUBIC C_{k,i}(t, y) = A - 3*B*y + 3*C*y^2 - E*y^3 FOR GLZ k=4..20 ===", flush=True)
def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p
def poly_add(a, b):
    n = max(len(a), len(b))
    return poly_trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def poly_deriv(a):
    return [i * a[i] for i in range(1, len(a))] if len(a) > 1 else [0]
def apply_Dj(p, j):
    dp = poly_deriv(p)
    res = [0] * (len(p) + 2)
    for i, c in enumerate(p):
        res[i] += c
        res[i+1] += j * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c
        res[i+2] += c
    return poly_trim(res)
def poly_eval(p, x):
    acc = mp.mpc(0) if isinstance(x, mp.mpc) else mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 21):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

# Test UHP stability of C_{k,i}(t, y) = A - 3*B*y + 3*C*y^2 - E*y^3 (and also without 3's, and with binomial coeff!)
t_uhp = [mp.mpc(-10, '1e-6'), mp.mpc(-1, '1e-6'), mp.mpc(-0.99, '1e-6'), mp.mpc(-0.5, '1e-6'),
         mp.mpc(-0.01, '1e-6'), mp.mpc(0, '1e-6'), mp.mpc(0.01, '1e-6'), mp.mpc(0.1, '1e-6'), mp.mpc(1, '1e-6'), mp.mpc(0, 1)]
for k in [5, 8, 12, 16, 20]:
    max_im = mp.mpf('-1e100')
    worst = None
    # Check middle & boundary i
    ilist = sorted(set([3, k//2, k//2+1, k-1]))
    for i in ilist:
        if i < 3 or i >= k: continue
        A = ht[k-1][i-3]
        B = apply_Dj(ht[k-1][i-2], k+1)
        C = apply_Dj(apply_Dj(ht[k-1][i-1], k+1), k+1)
        E = apply_Dj(apply_Dj(apply_Dj(ht[k-1][i], k+1), k+1), k+1)
        for tv in t_uhp:
            # C(t, y) = A - 3*B*y + 3*C*y^2 - E*y^3 => roots in y for t in UHP must have Im(y) <= 0!
            # Wait: check sign! For A << B << C << E, A - 3By + 3Cy^2 - Ey^3:
            # At k=4, i=3: A=(1+...), B=D(...): y-roots are POSITIVE!
            # As t increases, roots of A-Ey^3 move... wait: for Q = X - 2My + Yy^2,
            # Delta_Q = 2 W(X,M) - 2 W(X,Y) y + 2 W(M,Y) y^2 > 0 => Im(t) > 0 => Im(y) <= 0! Let's check BOTH min and max Im(y)!
            cy = [poly_eval(A, tv), -3*poly_eval(B, tv), 3*poly_eval(C, tv), -poly_eval(E, tv)]
            rts = mp.polyroots(list(reversed(cy)), maxsteps=200, extraprec=30)
            for rt in rts:
                if rt.imag > max_im:
                    max_im = rt.imag
                    worst = (i, complex(tv), complex(rt))
    print(f"  GLZ Cubic A - 3By + 3Cy^2 - Ey^3 at k={k:2d}: max Im(y) = {mp.nstr(max_im, 6)} at {worst[0], worst[1] if worst else None}", flush=True)

print("\n=== 3. CHECK MAX PLUCKER RATIO x^2 = W_AB W_CE / W_BC^2 FOR ANY REAL-STABLE CUBIC A - 3By + 3Cy^2 - Ey^3 ===", flush=True)
# Compute Rayleigh quadratic in y of C(t, y) = A - 3By + 3Cy^2 - Ey^3!
# C_t * C_y - C * C_ty = 3 * [W_AB - 2 W_AC y + (3 W_BC + W_AE) y^2 - 2 W_BE y^3 + W_CE y^4]!
# WAIT A SECOND!! CHECK THAT QUARTIC IN y!!
A_s, B_s, C_s, E_s, dA_s, dB_s, dC_s, dE_s, y_s = sp.symbols('A B C E dA dB dC dE y')
C_poly = A_s - 3*B_s*y_s + 3*C_s*y_s**2 - E_s*y_s**3
Ct_poly = dA_s - 3*dB_s*y_s + 3*dC_s*y_s**2 - dE_s*y_s**3
Cy_poly = sp.diff(C_poly, y_s)
Cty_poly = sp.diff(Ct_poly, y_s)
ray = sp.expand(Ct_poly * Cy_poly - C_poly * Cty_poly)
ray_p = sp.Poly(ray, y_s)
for deg in range(5):
    print(f"  [y^{deg}] Rayleigh(Cubic) = {sp.factor(ray_p.nth(deg))}")
