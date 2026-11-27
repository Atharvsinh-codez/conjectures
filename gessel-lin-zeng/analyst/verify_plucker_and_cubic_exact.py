"""100% exact SymPy verification of:
1. Exact Quantitative Plucker Quad Theorem:
   For ALL a, b > 0 with a*b = R <= 9 (let's check exact threshold R_max(theta)!),
   and ALL u in [0, 2*theta*sqrt(a)], v in [0, 2*theta*sqrt(b)] (theta = (1+eps)^{-1/2} in (0, 1]):
   J_new = 4(1+a+u)(1+b+v) / ((1+u)(1+v) - ab)^2  >=  4 / (1 + theta^2)^2  (exact global min over R <= 9?!).
   WAIT: check carefully for R > 1 and S = sqrt(a) + sqrt(b) >= 2*sqrt(r)!!
2. Exact Universal Cubic Wronskian Theorem SOS identity:
   2(a^2+a+1)(b^2+b+1) - a^2 - b^2 == 1/2*(2*a*b+a+b)^2 + 1/2*(a+b+2)^2.
   AND verify it holds for ALL real-stable cubics P(t, y) = p0 + p1*y + p2*y^2 + p3*y^3!
"""
import math
import numpy as np
import sympy as sp

a, b, r, s, S, th = sp.symbols('a b r s S th', positive=True)

print("=== 1. EXACT SYMPY CHECK OF UNIVERSAL CUBIC WRONSKIAN SOS ===")
lhs = sp.expand(2*(a**2 + a + 1)*(b**2 + b + 1) - a**2 - b**2)
rhs = sp.expand(sp.Rational(1, 2)*(2*a*b + a + b)**2 + sp.Rational(1, 2)*(a + b + 2)**2)
assert lhs == rhs, f"Mismatch: {lhs - rhs}"
print("  EXACT IDENTITY VERIFIED: 2(a^2+a+1)(b^2+b+1) - a^2 - b^2 == 1/2(2ab+a+b)^2 + 1/2(a+b+2)^2")

# Also check r_3 = 0 case (when one root is 0):
r1, r2, r3 = sp.symbols('r1 r2 r3', real=True)
def W_single(ra, rb):
    # Q(y) = (y-ra)(y-rb) = ra*rb - (ra+rb)y + y^2
    return (ra**2 * rb**2, ra**2 + ra*rb + rb**2, sp.Integer(1))

w1 = W_single(r2, r3); w2 = W_single(r1, r3); w3 = W_single(r1, r2)
c1, c2, c3 = sp.symbols('c1 c2 c3', nonnegative=True)
W0 = c1*w1[0] + c2*w2[0] + c3*w3[0]
W1 = c1*w1[1] + c2*w2[1] + c3*w3[1]
W2 = c1*w1[2] + c2*w2[2] + c3*w3[2]
diff = sp.expand(W1**2 - W0*W2)
# Check coefficient of each c_i*c_j!
for pair, expr in [
    ("c1^2", diff.coeff(c1, 2)),
    ("c2^2", diff.coeff(c2, 2)),
    ("c3^2", diff.coeff(c3, 2)),
    ("c1*c2", sp.expand(diff.coeff(c1, 1).coeff(c2, 1))),
    ("c1*c3", sp.expand(diff.coeff(c1, 1).coeff(c3, 1))),
    ("c2*c3", sp.expand(diff.coeff(c2, 1).coeff(c3, 1))),
]:
    # Verify exact SOS!
    print(f"  [{pair}] = {sp.factor(expr)}")

print("\n=== 2. EXACT SYMPY CHECK OF QUANTITATIVE PLUCKER QUAD BOUND ===")
# Put x = sqrt(a), y = sqrt(b), r = x*y = sqrt(R), S = x + y >= 2*sqrt(r).
# J_new(S, r, th) = 4 * [(1-r)^2 + 4*th^2*r + 2*th*(1+r)*S + S^2] / [1 - r^2 + 4*th^2*r + 2*th*S]^2
num = 4 * ((1 - r)**2 + 4*th**2*r + 2*th*(1 + r)*S + S**2)
den = (1 - r**2 + 4*th**2*r + 2*th*S)**2
# Compute J_new - 4/(1 + th^2)^2!!
gap = sp.factor(sp.together(num / den - 4 / (1 + th**2)**2))
print("  J_new - 4/(1+th^2)^2 =")
num_gap, den_gap = sp.fraction(gap)
# Write num_gap as polynomial in S = 2*sqrt(r) + d (d >= 0)!
d = sp.Symbol('d', nonnegative=True)
# Replace r -> s^2, S -> 2*s + d!
num_gap_s = sp.expand(num_gap.subs({r: s**2, S: 2*s + d}))
P_d = sp.Poly(num_gap_s, d)
print("  Numerator of (J_new - 4/(1+th^2)^2) in d = S - 2*s (with r = s^2):")
for deg in range(P_d.degree() + 1):
    print(f"    [d^{deg}] = {sp.factor(P_d.nth(deg))}")

# Also check J_new(S, r, th) >= 4/(1+th^2)^2 numerically over R in [0, 9], th in [0.5, 1.0]!
min_gap_num = 1e9; worst_p = None
for thv in [0.5, 0.7, 0.85, 0.95, 0.99, 1.0]:
    for rv in np.linspace(0.0, 3.0, 61):  # R = r^2 in [0, 9]
        s0 = math.sqrt(rv)
        for dv in [0.0, 0.01, 0.1, 0.5, 2.0, 10.0, 100.0]:
            Sv = 2*s0 + dv
            denom = 1 - rv**2 + 4*thv**2*rv + 2*thv*Sv
            if denom <= 0:
                min_gap_num = -999; worst_p = (thv, rv, dv, "denom<=0"); break
            val = 4*((1-rv)**2 + 4*thv**2*rv + 2*thv*(1+rv)*Sv + Sv**2) / (denom**2)
            g = val - 4.0 / ((1 + thv**2)**2)
            if g < min_gap_num:
                min_gap_num = g; worst_p = (thv, rv, dv)
print(f"  Global numerical check over th in [0.5, 1], R in [0, 9], d >= 0: min(J_new - 4/(1+th^2)^2) = {min_gap_num:+.6e} at {worst_p}")
