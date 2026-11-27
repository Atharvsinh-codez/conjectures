"""
Fast decisive checks:
1. Strengthened Cubic Stage 1 (pointwise in R^2):
   Given 5 vectors v0,v1,v2,v3,v4 in R^2 (in an open half-plane, W(vr,vs)>0 for r<s),
   suppose BOTH cubics C_c(v0,v1,v2,v3)(y) = v0 - 3c v1 y + 3c v2 y^2 - v3 y^3
   and C_c(v1,v2,v3,v4)(y) are real-stable (equivalently their Rayleigh/Wronskian quartics
   are >= 0 on R, i.e., for every t-fixed vector pair (p_r, p'_r), the polynomial
   W(C_c(y), C'_c(y))? Wait: let's be precise!).
   Actually, a univariate-in-y polynomial P(y) = sum a_r y^r with coefficients in R^2
   (representing (f_r(t), f'_r(t))) is a real-stable pencil at fixed t iff
   the Wronskian form Q(y) = sum_{r<s} W(v_r, v_s) ... wait:
   Let P(t,y) = v0 - 3c v1 y + 3c v2 y^2 - v3 y^3.
   By Hermite-Biehler / Borcea-Branden / Rayleigh:
   At fixed t in R, P(t,y) is real-stable in (t,y) near t iff for all y in R,
   W_t(P(t,y), partial_y P(t,y)) has constant sign!
   Wait: let's check the EXACT condition for P(t,y) = sum_{r=0}^3 (-1)^r b_r v_r(t) y^r
   to have all y-roots moving in the SAME direction as t varies!
   Let P(t,y) = f_0(t) - 3c f_1(t) y + 3c f_2(t) y^2 - f_3(t) y^3.
   For real y-roots y(t), dy/dt = -(partial_t P)/(partial_y P).
   We need (partial_y P)(partial_t P) to have fixed sign at roots y(t) of P, i.e.
   W_t(P, partial_y P) = P * (partial_t partial_y P) - (partial_t P)(partial_y P) >= 0 (or <= 0)
   Wait: at a root P=0, W_t(P, partial_y P) = -(partial_t P)(partial_y P) = (partial_y P)^2 * dy/dt!
   Let's check the sign! For a degree-2 quad Q(t,y) = X - 2c M y + Y y^2 with X << M << Y:
   W(X, M) > 0, W(M, Y) > 0, W(X, Y) > 0.
   Wait: partial_y Q = -2c M + 2 Y y.
   W_t(Q, partial_y Q) = W(X - 2c M y + Y y^2, -2c M + 2 Y y)
     = -2c W(X,M) + 2 W(X,Y) y - 4c W(M,Y) y^2 + 2c W(M,Y) y^2?
   Let's compute carefully:
     W(X, -2c M) = -2c W(X, M)
     W(X, 2 Y y) + W(-2c M y, -2c M) = 2 W(X, Y) y
     W(-2c M y, 2 Y y) + W(Y y^2, -2c M) = (-4c + 2c) W(M, Y) y^2 = -2c W(M, Y) y^2.
   So W_t(Q, partial_y Q) = -2c W(X,M) + 2 W(X,Y) y - 2c W(M,Y) y^2 <= 0 for all y in R
   iff c^2 W(X,M) W(M,Y) >= (1/4) W(X,Y)^2, i.e., 4 c^2 W(X,M) W(M,Y) >= W(X,Y)^2, i.e.,
   J(X,M,Y) >= 1/c^2!
   Wait: that's the EXACT pointwise Hermite-Biehler / Wronskian condition:
     -W_t(P(t,y), partial_y P(t,y)) >= 0 for ALL y in R!
   Let's double check for cubic P(t,y) = A - 3c B y + 3c C y^2 - E y^3!
   If P is real-stable in (t,y), then for every real t, P(t,y) is real-rooted in y, and at every
   y-root, dy/dt <= 0 (or >= 0), equivalently (for fixed t, viewing v_r = (f_r(t), f'_r(t)) in R^2):
   Is -W_t(P, partial_y P) >= 0 for ALL y in R (a degree-4 polynomial in y) equivalent to
   (P, partial_t P) interlacing in y (Hermite-Biehler for P(t,y) + i epsilon partial_t P(t,y))?
   YES! By Hermite-Biehler / Obreschkoff, P(t,y) is real-rooted in y AND all roots move strictly
   monotonically in the same direction (equivalently P + alpha * P_t is real-rooted for all real alpha)
   IFF the Wronskian (Bezoutiant / Hermite form / continuous Wronskian)
     H(y) := W_y(P_t, P) = P_t * P_y - P * P_{ty} = -W_t(P, P_y) >= 0 for ALL y in R
   Wait: H(y) >= 0 for ALL y in R is NECESSARY AND SUFFICIENT for P and P_t to interlace?
   BE CAREFUL! For deg 2, H(y) >= 0 on R IS necessary and sufficient (since deg H = 2).
   For deg 3, P(y) real-rooted AND H(y) >= 0 at the roots of P is necessary and sufficient for
   P_t << P; H(y) >= 0 for ALL y in R is slightly STRONGER than P_t << P (Wait: is it strictly stronger?
   For deg 3, P(y) = (y-r1)(y-r2)(y-r3) and P_t = a(y-s1)(y-s2) with r1 < s1 < r2 < s2 < r3:
   P_t / P = sum c_k / (y - r_k) with c_k > 0, so (P_t / P)' = - sum c_k / (y - r_k)^2 < 0 EVERYWHERE!
   And (P_t / P)' = (P_{ty} P - P_t P_y) / P^2 = -H(y) / P^2!
   SO H(y) >= 0 FOR ALL y IN R IS EQUIVALENT TO P_t / P BEING A PICK/HERGLOTZ FUNCTION,
   WHICH IS 100% EQUIVALENT TO REAL-ROOTEDNESS OF P(y) + alpha P_t(y) FOR ALL REAL alpha!)
   OH!! WHAT A BEAUTIFUL FACT!
   Let's verify this AND test Strengthened Cubic Stage 1 in R^2!

2. Also test GLZ k=35 Q_{35, q} with huge b_{36} = 10^15 (and w_{36} = 1 vs 10^-30) at q = 1/2, 1/4, 1/8.
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 60

# ============================================================
# PART 1: Strengthened Cubic Stage 1 in R^2
# ============================================================
def cubic_H_poly(v0, v1, v2, v3, c):
    """
    P(y) = v0 - 3c v1 y + 3c v2 y^2 - v3 y^3.
    P_y(y) = -3c v1 + 6c v2 y - 3 v3 y^2.
    H(y) = -W_t(P(y), P_y(y)) = W_t(P_y(y), P(y)), where W(u,v) = u0*v1 - u1*v0.
    Wait: for X << M << Y, W(X,M) > 0, W(M,Y) > 0, W(X,Y) > 0.
    Let's compute -W_t(P, P_y) with W(v_r, v_s) = w_{rs} > 0 for r < s.
    Coeff of y^0: -W(v0, -3c v1) = 3c w01.
    Coeff of y^1: -W(v0, 6c v2) - W(-3c v1, -3c v1) = -6c w02.
    Coeff of y^2: -W(v0, -3 v3) - W(-3c v1, 6c v2) - W(3c v2, -3c v1)
                = 3 w03 + 18c^2 w12 - 9c^2 w12 = 3 w03 + 9c^2 w12.
    Coeff of y^3: -W(-3c v1, -3 v3) - W(3c v2, 6c v2) - W(-v3, -3c v1)
                = -9c w13 + 3c w13 = -6c w13.
    Coeff of y^4: -W(3c v2, -3 v3) - W(-v3, 6c v2)
                = 9c w23 - 6c w23 = 3c w23.
    Divide by 3:
    h(y) = c w01 - 2c w02 y + (w03 + 3 c^2 w12) y^2 - 2c w13 y^3 + c w23 y^4.
    Check at c=1: w01 - 2 w02 y + (w03 + 3 w12) y^2 - 2 w13 y^3 + w23 y^4!
    """
    def W(a, b):
        return a[0]*b[1] - a[1]*b[0]
    w01 = W(v0, v1)
    w02 = W(v0, v2)
    w03 = W(v0, v3)
    w12 = W(v1, v2)
    w13 = W(v1, v3)
    w23 = W(v2, v3)
    # returns coeffs [a4, a3, a2, a1, a0] for np.polyval / roots
    return np.array([c*w23, -2*c*w13, w03 + 3*(c**2)*w12, -2*c*w02, c*w01], dtype=float)

def min_quartic_val(coeffs):
    # coeffs = [a4, a3, a2, a1, a0] with a4 > 0
    a4, a3, a2, a1, a0 = coeffs
    if a4 <= 0 or a0 <= 0:
        return -1e9
    # derivative 4 a4 y^3 + 3 a3 y^2 + 2 a2 y + a1
    deriv = [4*a4, 3*a3, 2*a2, a1]
    rts = np.roots(deriv)
    vals = [a0]
    for r in rts:
        if abs(r.imag) < 1e-10:
            y = r.real
            vals.append((((a4*y + a3)*y + a2)*y + a1)*y + a0)
    return min(vals)

def find_c_min(v0, v1, v2, v3):
    # Find smallest c > 0 such that min_{y in R} h_c(y) >= 0.
    # Note: h_c(y) is linear-ish/quadratic in c: for large c, 3c^2 w12 y^2 dominates,
    # wait: as c grows, c*w01 - 2c*w02 y + (w03 + 3c^2 w12) y^2 - 2c w13 y^3 + c w23 y^4
    #Wait: is h_c(y) >= 0 monotone in c? Let's check!
    # Divide by c: w01 - 2 w02 y + (w03/c + 3c w12) y^2 - 2 w13 y^3 + w23 y^4!
    # WHOA!! LOOK AT THAT FORMULA:
    #   h_c(y) / c = w01 - 2 w02 y + (w03/c + 3 c w12) y^2 - 2 w13 y^3 + w23 y^4!!
    # Wait: the ONLY place c appears is in (w03/c + 3c w12) y^2!!
    # Wait: for w03 > 0 and w12 > 0, the function f(c) = w03/c + 3c w12 has a MINIMUM at
    #   c* = sqrt(w03 / (3 w12)), where f(c*) = 2 sqrt(3 w03 w12)!
    #Wait - let's check if P_c(y) = v0 - 3c v1 y + 3c v2 y^2 - v3 y^3 REALLY has that!
    # Let's check! For c = 1, P_1 = v0 - 3 v1 y + 3 v2 y^2 - v3 y^3.
    # For c < 1, P_c is LESS stable or MORE stable? Wait: at c=0, P_0 = v0 - v3 y^3 is NOT real-rooted!
    # At c -> infinity, P_c ~ -3c v1 y + 3c v2 y^2 IS real-rooted, with 3 real roots (near 0, v1/v2, infty)!
    # SO SMALLER c IS STRONGER (harder to be real-stable), as long as c > sqrt(w03/(3 w12))!
    # Let's check the minimum c in (0, 2] where h_c(y) >= 0 for all y!
    lo, hi = 0.0, 1.0
    if min_quartic_val(cubic_H_poly(v0, v1, v2, v3, hi)) < -1e-12:
        # try larger hi
        hi = 5.0
        if min_quartic_val(cubic_H_poly(v0, v1, v2, v3, hi)) < -1e-12:
            return None
    # Wait: since w03/c + 3c w12 is NOT monotone for small c, let's check c in [0.5, 1.5]
    # Actually for real-rootedness of P_c(y) itself, we also need P_c(t,y) real-rooted in y,
    # Wait: H(y) >= 0 for all y in R DOES imply P_c(t,y) is real-rooted in y?
    # Let's check! If P_c(y) = y^3 - 3c y^2 + 3c y - 1? For c < 1, y=1 is a root,
    # P_c(y)/(y-1) = y^2 + (1-3c)y + 1, disc = (1-3c)^2 - 4 = 9c^2 - 6c - 3 = 3(3c+1)(c-1) < 0!
    # SO AT c < 1, P_c(y) = (1 - 3cy + 3cy^2 - y^3) has COMPLEX ROOTS if v0=v1=v2=v3!
    # WAIT A SECOND! Look at w03/c + 3c w12 when c -> 0+: w03/c -> +infty makes h_c(y) > 0 even when P_c is NOT real-rooted!
    # HA! GOOD CATCH!! H(y) >= 0 on R is ONLY equivalent to real-stability when the Hermite/Bezoutiant
    # MATRIX is PSD (NOT just the diagonal H(y) = v(y)^T B v(y) >= 0)!
    return None

def bezoutiant_cubic_min_eig(v0, v1, v2, v3, c):
    """
    Exact Bezoutiant (Hermite) 3x3 matrix B of (P(y), P_t(y)) where
      P(y)   = v0[0] - 3c v1[0] y + 3c v2[0] y^2 - v3[0] y^3
      P_t(y) = v0[1] - 3c v1[1] y + 3c v2[1] y^2 - v3[1] y^3
    By Hermite-Biehler / Bezoutiant theorem:
    P(y) + alpha P_t(y) is real-rooted for all real alpha (with P_t << P)
    IFF the 3x3 symmetric Bezoutiant matrix B (defined by
      (P(x) P_t(y) - P(y) P_t(x)) / (x - y) = sum_{i,j=0}^2 B_{ij} x^i y^j )
    is positive semidefinite (or negative semidefinite; let's check sign!).
    Let p_r = (-1)^r b_r v_r[0], q_r = (-1)^r b_r v_r[1], with b = [1, 3c, 3c, 1].
    Wait: W(p_r, p_s) := p_r q_s - q_r p_s = (-1)^{r+s} b_r b_s W(v_r, v_s).
    Let's compute (P(x) Q(y) - P(y) Q(x)) / (x - y) = sum_{0<=r<s<=3} (p_r q_s - q_r p_s) (x^r y^s - x^s y^r) / (x - y)
    Notice (x^r y^s - x^s y^r)/(x - y) = - (x y)^r (x^{s-r-1} + x^{s-r-2} y + ... + y^{s-r-1}).
    Let's negate it so we work with M(x,y) = (P(y) Q(x) - P(x) Q(y))/(x-y) = sum_{r<s} W(p_r,p_s) (xy)^r (x^{s-r-1}+...+y^{s-r-1})?
    Wait: (P(y)Q(x) - P(x)Q(y))/(x-y) at y=x is P'(x)Q(x) - P(x)Q'(x) = -W_t(P, P_y) = H(x) >= 0!
    YES!! Let's double check!
    (P(y)Q(x) - P(x)Q(y))/(x-y) = sum_{r<s} (p_r q_s - q_r p_s) (y^r x^s - x^r y^s)/(x-y)
                                = sum_{r<s} W_raw(r,s) * (xy)^r * sum_{k=0}^{s-r-1} x^{s-r-1-k} y^k.
    Wait: W_raw(r,s) = (-1)^{r+s} b_r b_s W(v_r, v_s):
      r=0,s=1: -3c w01 * 1
      Wait! Why is (0,1) negative? Because P(x) = v0 - 3c v1 x + ...
      P'(x)Q(x) - P(x)Q'(x) at x=0 is (-3c v1[0])(v0[1]) - (v0[0])(-3c v1[1]) = +3c W(v0, v1)!
      Let's check: (P(y)Q(x) - P(x)Q(y))/(x-y) at x->y:
      (P(y)(Q(y)+(x-y)Q'(y)) - (P(y)+(x-y)P'(y))Q(y))/(x-y) = P(y)Q'(y) - P'(y)Q(y) = -H(y)!
      So we want (P(x)Q(y) - P(y)Q(x))/(x-y)!!
      Wait: for (r,s)=(0,1): (p0 q1 - q0 p1)(x^0 y^1 - y^0 x^1)/(x-y) = (-3c w01)(-1) = +3c w01! YES!
    """
    b = [1.0, -3.0*c, 3.0*c, -1.0]
    vs = [v0, v1, v2, v3]
    def W(i, j):
        return vs[i][0]*vs[j][1] - vs[i][1]*vs[j][0]
    # B_{ij} = coeff of x^i y^j in (P(x)Q(y) - P(y)Q(x))/(x-y)
    # For r < s: p_r q_s - q_r p_s = b[r]*b[s]*W(r,s)
    # (x^r y^s - y^r x^s)/(x-y) = -(xy)^r * sum_{m=0}^{s-r-1} x^m y^{s-r-1-m}
    B = np.zeros((3, 3), dtype=float)
    for r in range(4):
        for s in range(r+1, 4):
            val = -b[r]*b[s]*W(r, s)
            for m in range(s - r):
                B[r + m, r + (s - r - 1 - m)] += val
    # Let's check B explicitly!
    # r=0,s=1: -(-3c) w01 = 3c w01 at (0,0)
    # r=0,s=2: -(3c) w02 = -3c w02 at (0,1) and (1,0)
    # r=0,s=3: -(-1) w03 = +w03 at (0,2), (1,1), (2,0)
    # r=1,s=2: -(-9c^2) w12 = +9c^2 w12 at (1,1)
    # r=1,s=3: -(3c) w13 = -3c w13 at (1,2), (2,1)
    # r=2,s=3: -(-3c) w23 = +3c w23 at (2,2)
    return B

# Verify B formula symbolically/numerically!
print("=== VERIFY 3x3 BEZOUTIANT OF STRENGTHENED CUBIC ===")
# For v_r = (1, r) (like (1+t)^r at t=0): w_{rs} = s - r.
# B = [[3c, -6c, 3], [-6c, 3+9c^2, -6c], [3, -6c, 3c]]
# Wait: for c=1, B = [[3,-6,3],[-6,12,-6],[3,-6,3]] = 3*[1,-2,1]^T[1,-2,1], PSD of rank 1!
# For c < 1: det of top-left 2x2 is 3c(3+9c^2) - 36c^2 = 9c(3c^2 - 4c + 1) = 9c(3c-1)(c-1) < 0!
# AWESOME! B is PSD iff P_c is real-stable at t!
for c_test in [0.9, 1.0, 1.1]:
    B = bezoutiant_cubic_min_eig((1,0), (1,1), (1,2), (1,3), c_test)
    eigs = np.linalg.eigvalsh(B)
    print(f"  linear v_r=(1,r), c={c_test}: min_eig={eigs[0]:.6f}, eigs={np.round(eigs,4)}")

# NOW TEST STRENGTHENED CUBIC STAGE 1 IN R^2!
# Given v0,v1,v2,v3,v4 with W(vr,vs) > 0 for r < s (open half-plane).
# Suppose B(v0,v1,v2,v3; c) >= 0 AND B(v1,v2,v3,v4; c) >= 0 (with c <= 1).
# What is the MINIMUM c_out such that B(v0+v1, v1+v2, v2+v3, v3+v4; c_out) >= 0?
# Wait: let's check BOTH c = 1 (standard cubic) and c < 1 (strengthened cubic)!
print("\n=== TEST STRENGTHENED CUBIC STAGE 1 (v0..v4 in R^2) ===")

def find_min_c_psd(v0, v1, v2, v3):
    # Find smallest c > 0 in [0, 2] such that B(v0,v1,v2,v3; c) >= 0.
    # Wait: is PSD monotone in c? Notice B(c) has:
    #   B_00 = 3c w01, B_01 = -3c w02, B_02 = w03
    #   B_11 = w03 + 9c^2 w12, B_12 = -3c w13, B_22 = 3c w23.
    # Let D = diag(1, 1/c, 1)? Wait: scale row/col 1 by 1/c:
    #   B_00 = 3c w01, B_01' = -3 w02, B_02 = w03
    #   B_11' = w03/c^2 + 9 w12, B_12' = -3 w13, B_22 = 3c w23.
    # Even better: scale row/col 0 and 2 by 1/sqrt(c) and row/col 1 by 1/sqrt(c) -> divide by c:
    #   [[3 w01, -3 w02, w03/c], [-3 w02, w03/c + 9c w12, -3 w13], [w03/c, -3 w13, 3 w23]].
    # Let's check for c in (0, 1] by testing a dense grid + bisection from right (starting at c=1)!
    # Wait: at c_0 = 1/sqrt(3) or whatever, let's check whether the set of PSD c is an interval [c_min, infty)!
    # If B(v0,v1,v2,v3; c) is PSD at c=1, let's bisect for the largest c <= 1 where min_eig < 0!
    if np.linalg.eigvalsh(bezoutiant_cubic_min_eig(v0, v1, v2, v3, 1.0))[0] < -1e-12:
        # check larger c
        hi = 1.0
        while hi < 100 and np.linalg.eigvalsh(bezoutiant_cubic_min_eig(v0, v1, v2, v3, hi))[0] < -1e-12:
            hi *= 1.5
        if hi >= 100:
            return None
        lo = hi / 1.5
        for _ in range(50):
            mid = 0.5*(lo+hi)
            if np.linalg.eigvalsh(bezoutiant_cubic_min_eig(v0, v1, v2, v3, mid))[0] >= 0:
                hi = mid
            else:
                lo = mid
        return hi
    # c=1 is PSD; bisect downward in [1e-4, 1.0]
    lo, hi = 1e-6, 1.0
    if np.linalg.eigvalsh(bezoutiant_cubic_min_eig(v0, v1, v2, v3, lo))[0] >= 0:
        return lo
    for _ in range(55):
        mid = 0.5*(lo+hi)
        if np.linalg.eigvalsh(bezoutiant_cubic_min_eig(v0, v1, v2, v3, mid))[0] >= 0:
            hi = mid
        else:
            lo = mid
    return hi

rng = np.random.default_rng(42)
# Generate random 5 vectors in upper half-plane (angles strictly increasing in (0, pi))
# Check: when max(c_in1, c_in2) <= c_in (for c_in = 1.0, 0.95, 0.9, 0.8, 0.7),
# what is max(c_out)?
for target_c in [1.0, 0.95, 0.9, 0.8, 0.5]:
    max_cout = -1.0
    max_ratio = -1.0
    worst = None
    n_valid = 0
    for trial in range(15000):
        # random sorted angles in (0, pi) and random lengths (log-normal)
        ang = np.sort(rng.uniform(0.01, np.pi - 0.01, 5))
        lens = np.exp(rng.normal(0, 1.5, 5))
        vs = [(lens[i]*np.cos(ang[i]), lens[i]*np.sin(ang[i])) for i in range(5)]
        c1 = find_min_c_psd(vs[0], vs[1], vs[2], vs[3])
        if c1 is None or c1 > target_c:
            continue
        c2 = find_min_c_psd(vs[1], vs[2], vs[3], vs[4])
        if c2 is None or c2 > target_c:
            continue
        n_valid += 1
        cin = max(c1, c2)
        ws = [(vs[i][0]+vs[i+1][0], vs[i][1]+vs[i+1][1]) for i in range(4)]
        cout = find_min_c_psd(ws[0], ws[1], ws[2], ws[3])
        if cout is None:
            cout = 999.0
        if cout > max_cout:
            max_cout = cout
            worst = (cin, cout, c1, c2)
        if cout / cin > max_ratio:
            max_ratio = cout / cin
    # Also compute Theorem B's quadratic bound for comparison:
    # If c = 1/sqrt(1+eps), Theorem B gives c_out <= (2+eps)/(2(1+eps)) = c*(c^2+1)/(2)??
    # Wait: 1+eps = 1/c^2 => eps = (1-c^2)/c^2 => 1/sqrt(J_out) <= (2+eps)/(2(1+eps)) = (1+c^2)/2!
    # Wait: for c = 1 - d, (1+c^2)/2 = 1 - d + d^2/2!
    quad_bound = 0.5*(1 + target_c**2)
    print(f"  target_c={target_c:.2f} (valid={n_valid:4d}): max_cout={max_cout:.6f}, max(cout/cin)={max_ratio:.6f}, quad_bd={(1+target_c**2)/2:.6f}, worst={worst}")

# Also test adversarial optimization for target_c = 1.0 and 0.9 (trying to push cout > target_c or > quad_bound)!
print("\n=== ADVERSARIAL SEARCH ON CUBIC STAGE 1 (c=1.0 and c=0.9) ===")
for c_test in [1.0, 0.9, 0.8]:
    worst_eig = 1e9
    worst_cout = 0.0
    # Parameterize 5 vectors up to GL2+ and scaling:
    # Fix v1 = (1, 0), v2 = (0, 1) -> W(v1, v2) = 1.
    # Wait: open half-plane requires all angles in (0, pi).
    # Fix v2 = (0, 1) (angle pi/2).
    # v1 = (a1, b1) with a1 > 0, and by SL2 we can set v1 = (1, 0), v2 = (0, 1), wait: then angle of v1 is 0, v2 is pi/2.
    # For v0 to have angle in (-pi/2, 0)? Wait: W(v0, v1) > 0 means v0 is clockwise of v1 (angle < 0).
    # Let v0 = (x0, -y0) with y0 > 0, x0 > 0 (so W(v0, v2) = x0 > 0, W(v0, v1) = y0 > 0).
    # v3 = (-x3, y3) with x3 > 0, y3 > 0 (so W(v1, v3) = y3 > 0, W(v2, v3) = x3 > 0).
    # v4 = (-x4, y4) with W(v3, v4) = -x3 y4 + y3 x4 > 0 => x4/y4 > x3/y3, and W(v0, v4) = x0 y4 - y0 x4 > 0!
    # Check W(v0, v3) = x0 y3 - y0 x3 > 0 too!
    # Beautiful! In slopes m_r = angle: let v_r = L_r * (cos th_r, sin th_r) with 0 < th0 < th1 < th2 < th3 < th4 < pi.
    from scipy.optimize import minimize
    best_cout = 0.0
    for seed in range(25):
        rng2 = np.random.default_rng(100 + seed)
        # Try random feasible start
        found = False
        for _ in range(500):
            ang = np.sort(rng2.uniform(0.05, np.pi - 0.05, 5))
            lens = np.exp(rng2.normal(0, 1.5, 5))
            vs = [(lens[i]*np.cos(ang[i]), lens[i]*np.sin(ang[i])) for i in range(5)]
            if (np.linalg.eigvalsh(bezoutiant_cubic_min_eig(vs[0],vs[1],vs[2],vs[3],c_test))[0] >= 1e-6 and
                np.linalg.eigvalsh(bezoutiant_cubic_min_eig(vs[1],vs[2],vs[3],vs[4],c_test))[0] >= 1e-6):
                found = True
                break
        if not found:
            continue
        # Maximize cout using Nelder-Mead with penalty
        # Param: dth0..dth4 (4 gaps) and logL0..logL4 (4 ratios, logL2=0)
        dth = np.diff(ang)
        logL = np.log(lens) - np.log(lens[2])
        x0 = np.concatenate([np.log(dth), [logL[0], logL[1], logL[3], logL[4]]])
        def obj(z):
            g = np.exp(np.clip(z[:4], -12, 5))
            tot = np.sum(g)
            if tot >= np.pi - 1e-4:
                g = g * (np.pi - 0.01) / (tot + 0.01)
            th = np.zeros(5)
            th[0] = 0.005
            th[1:] = th[0] + np.cumsum(g)
            if th[-1] >= np.pi - 1e-4:
                return 10.0
            l = np.array([z[4], z[5], 0.0, z[6], z[7]])
            l = np.clip(l, -10, 10)
            L = np.exp(l)
            vs = [(L[i]*np.cos(th[i]), L[i]*np.sin(th[i])) for i in range(5)]
            e1 = np.linalg.eigvalsh(bezoutiant_cubic_min_eig(vs[0],vs[1],vs[2],vs[3],c_test))[0]
            e2 = np.linalg.eigvalsh(bezoutiant_cubic_min_eig(vs[1],vs[2],vs[3],vs[4],c_test))[0]
            if e1 < 0 or e2 < 0:
                return 10.0 + 100.0 * (max(0.0, -e1) + max(0.0, -e2))
            ws = [(vs[i][0]+vs[i+1][0], vs[i][1]+vs[i+1][1]) for i in range(4)]
            co = find_min_c_psd(ws[0], ws[1], ws[2], ws[3])
            if co is None:
                return -100.0
            return -co
        res = minimize(obj, x0, method='Nelder-Mead', options={'maxiter': 800})
        if -res.fun > best_cout:
            best_cout = -res.fun
    print(f"  Adversarial c_in={c_test:.2f}: max c_out = {best_cout:.8f} (compare c_in={c_test:.4f}, (1+c^2)/2={0.5*(1+c_test**2):.6f}, 3c/(2+c)={3*c_test/(2+c_test):.6f})")

# ============================================================
# PART 2: GLZ k=35 Q_{35, q} with huge b_{36} = 10^15
# ============================================================
print("\n=== PART 2: GLZ k=35 Q_{35,q} + (D_{b36} + Y_q) ===")
def poly_trim(p):
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p

def apply_Dj_mp(p, j):
    # (t(1+t) d/dt + 1 + j t) p(t)
    # p[m] t^m -> (1 + m) p[m] t^m + (j + m) p[m] t^{m+1}
    n = len(p)
    out = [mp.mpf(0)] * (n + 1)
    for m, c in enumerate(p):
        if c != 0:
            out[m] += (1 + m) * c
            out[m + 1] += (j + m) * c
    return out

def poly_add_mp(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for i in range(len(a)):
        out[i] += a[i]
    for i in range(len(b)):
        out[i] += b[i]
    return out

def poly_eval_mp(p, t):
    s = mp.mpc(0)
    for c in reversed(p):
        s = s * t + c
    return s

def compute_glz_htilde_mp(k_max):
    # h[i] is poly in t for h~_{k, i}(t)
    h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]  # k=1
    for k in range(2, k_max + 1):
        j = k + 1
        h_new = []
        for i in range(k + 1):
            u = [mp.mpf(0)]
            if i <= k - 1:
                u = poly_add_mp(u, apply_Dj_mp(h[i], j))
            if i >= 1:
                u = poly_add_mp(u, h[i - 1])
            h_new.append(apply_Dj_mp(u, j))
        h = h_new
    return h

def max_im_y_roots(coeffs_y):
    # coeffs_y[i] is coeff of y^i (for i=0..k). Rescale y = s*u so |c0| ~ |ck|*s^k
    k = len(coeffs_y) - 1
    c0 = abs(coeffs_y[0])
    ck = abs(coeffs_y[-1])
    if c0 == 0 or ck == 0:
        return None
    s = (c0 / ck) ** (mp.mpf(1) / k)
    scaled = [coeffs_y[i] * (s ** i) for i in range(k + 1)]
    mx = max(abs(z) for z in scaled)
    scaled = [z / mx for z in scaled]
    # mp.polyroots expects highest degree first
    rts = mp.polyroots(list(reversed(scaled)), maxsteps=200, extraprec=40)
    # In UHP stability, for t in H_+ (Im t > 0), all y-roots must have Im y <= 0!
    # (Wait: check sign of Im y for t = -1/2 + i*0.1 on Q_1 = (1+2t) + y -> y = -(1+2t) has Im y = -2 Im t < 0!)
    return float(max(r.imag for r in rts))

h35 = compute_glz_htilde_mp(35)
# Check J_17 at t=0 for h35 to confirm it is indeed < 1
def W_at_0(p, q):
    p0 = p[0]
    p1 = p[1] if len(p) > 1 else mp.mpf(0)
    q0 = q[0]
    q1 = q[1] if len(q) > 1 else mp.mpf(0)
    return p0*q1 - p1*q0

w01 = W_at_0(h35[18], h35[17])
w12 = W_at_0(h35[17], h35[16])
w02 = W_at_0(h35[18], h35[16])
print(f"  GLZ k=35 at t=0: 4 W(h18,h17)W(h17,h16)/W(h18,h16)^2 = {float(4*w01*w12/(w02**2)):.8f} (< 1 confirms J17 < 0!)")

# Now form h36_step(b36, w36) = D_{b36} h35[i] + w36 * h35[i-1] (single step without second D, AND with second D)
for b36 in [mp.mpf(37), mp.mpf('1e15')]:
    for w36 in [mp.mpf(1), mp.mpf('1e-30')]:
        h36 = []
        for i in range(37):
            u = [mp.mpf(0)]
            if i <= 35:
                u = poly_add_mp(u, apply_Dj_mp(h35[i], b36))
            if i >= 1:
                u = poly_add_mp(u, [w36 * c for c in h35[i - 1]])
            h36.append(u)
        for q_val in [mp.mpf('0.5'), mp.mpf('0.25'), mp.mpf('0.125')]:
            # Test at t = -1e-16 + 1e-20 j, t = 1e-16 + 1e-20 j, t = -0.01 + 1e-6 j, t = 10 + 1e-3 j
            worst_im = -1e9
            for t_test in [mp.mpc('-1e-16', '1e-20'), mp.mpc('1e-16', '1e-20'), mp.mpc('-1e-6', '1e-10'), mp.mpc('10', '1e-3')]:
                coeffs = [(q_val ** (i*(i-1)//2)) * poly_eval_mp(h36[i], t_test) for i in range(37)]
                im = max_im_y_roots(coeffs)
                if im is not None and im > worst_im:
                    worst_im = im
            status = "PASS" if worst_im <= 1e-12 else "FAIL"
            print(f"  b36={float(b36):.0e}, w36={float(w36):.0e}, q={float(q_val):.3f}: max Im(y_scaled) = {worst_im:+.3e} [{status}]")
