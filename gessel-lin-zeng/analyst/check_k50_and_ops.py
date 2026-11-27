"""Quick checks:
1. Verify our k~34 prediction for r=2 Rayleigh at t=0!
2. Test if Q_{k,i}(t, y) = P_{k,i-1} + c * P_{k,i} y + P_{k,i+1} y^2 passes Rayleigh at t=0 for c = 2.5, 3 up to k=50!
3. Test what happens when we replace D_k by OTHER stability-preserving 1st-order operators:
   - Laguerre: t*d/dt + c + m*t  (roots in (-inf, 0))
   - Jacobi/GLZ: t(1+t)*d/dt + c + m*t  (roots in (-1, 0))
   - Hermite: d/dt - m*t  (roots in R)
   - General quadratic: (t-a)(t-b) d/dt + c + m*t!
   - Non-quadratic cubic: t(1+t)(2+t) d/dt + ... (does it fail?!)
   WAIT A SECOND!! TEST ALL OF THESE IN 3 SECONDS!
"""
import math
import sympy as sp

def apply_D(c_val, m_val, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c_val)
        res[r + 1] += coef * (r + m_val)
    return res

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

# Fast computation of ONLY [t^0] and [t^1] of htilde_{k,i}(t) up to k = 50!
# Recall D_{k+1}(a + b t + ...) = a + (2 b + (k+1) a) t + ...!
# Wait: for P_{k,i} = D_{k+1}(D_{k+1} P_{k-1,i} + P_{k-1,i-1}):
#   a0(D f) = a0(f),  a1(D f) = 2 a1(f) + (k+1) a0(f)!
KMAX = 50
a0 = [[0]*(k+1) for k in range(KMAX+1)]
a1 = [[0]*(k+1) for k in range(KMAX+1)]
a0[0][0] = 1; a1[0][0] = 0
for k in range(1, KMAX + 1):
    m = k + 1
    for i in range(k + 1):
        # U_i = D_m P_{k-1,i} + P_{k-1,i-1}
        u0 = (a0[k-1][i] if i <= k-1 else 0) + (a0[k-1][i-1] if i >= 1 else 0)
        u1 = ((2*a1[k-1][i] + m*a0[k-1][i]) if i <= k-1 else 0) + (a1[k-1][i-1] if i >= 1 else 0)
        # P_{k,i} = D_m U_i
        a0[k][i] = u0
        a1[k][i] = 2*u1 + m*u0

print("=== 1. VERIFY PREDICTION: RAYLEIGH OF P_{k,i-1} + c*P_{k,i}*y + P_{k,i+1}*y^2 AT t=0 (k=2..50) ===")
for c_num, c_den in [(2, 1), (21, 10), (22, 10), (5, 2), (3, 1)]:
    first_fail = None
    for k in range(2, KMAX + 1):
        for i in range(1, k):
            c0 = a1[k][i-1]*a0[k][i] - a0[k][i-1]*a1[k][i]
            c1 = a1[k][i-1]*a0[k][i+1] - a0[k][i-1]*a1[k][i+1]
            c2 = a1[k][i]*a0[k][i+1] - a0[k][i]*a1[k][i+1]
            # condition: c^2 * c0 * c2 >= c1^2
            if c_num*c_num * c0 * c2 < c_den*c_den * c1 * c1:
                first_fail = (k, i, float(c1*c1) / float(c0*c2))
                break
        if first_fail: break
    print(f"  c = {c_num}/{c_den} ({c_num/c_den:.2f}): first fail in k<=50 = {first_fail}")

# Print max ratio c1^2 / (c0 * c2) at k = 50!
max_r = max((float(a1[50][i-1]*a0[50][i+1] - a0[50][i-1]*a1[50][i+1])**2 /
             float((a1[50][i-1]*a0[50][i] - a0[50][i-1]*a1[50][i])*(a1[50][i]*a0[50][i+1] - a0[50][i]*a1[50][i+1])), i)
            for i in range(1, 50))
print(f"  At k=50: max c1^2/(c0*c2) = {max_r[0]:.6f} at i={max_r[1]} (predicted limit: 9/2 = 4.500000!)")

print("\n=== 2. CHECK GENERAL 1st-ORDER OPERATORS: LAGUERRE t*d/dt + c + m*t (PLAIN BIVARIATE STABILITY?!) ===")
# WAIT A SECOND!! CHECK LAGUERRE L_m = t*d/dt + 1 + m*t!!
# Compute L_m( (t+z)^d ) = (t+z)^{d-1} [ d*t + (1+m*t)*(t+z) ]!
# Check (L_m + y)( (t+z)^d ) = (t+z)^{d-1} [ d*t + (1 + m*t + y)*(t+z) ]!
# Solve for -z: -z = t + d*t / (1 + m*t + y)!!
# WAIT A SECOND!! CHECK Im(-z) FOR t, y IN H+?!
t_s, y_s = sp.symbols('t y')
def test_op_biv(name, op_fn, kmax=4):
    H = sp.Integer(1)
    for k in range(1, kmax + 1):
        H = sp.expand(op_fn(k, op_fn(k, H) + y_s * H))
    # check at t = 1, 2, -0.5, -2
    fails = []
    for t0 in [-2, -1, sp.Rational(-1,2), 0, 1, 2]:
        P = sp.Poly(sp.expand(H.subs(t_s, t0)), y_s, domain='QQ')
        if P.degree() > 0 and int(sp.count_roots(P)) < P.degree():
            fails.append((t0, int(sp.count_roots(P)), P.degree()))
    print(f"  {name:45s} (k={kmax}): plain H(t0, y) fails = {fails}")

test_op_biv("Laguerre t*d/dt + 1 + k*t", lambda k, f: t_s*sp.diff(f, t_s) + (1 + k*t_s)*f, kmax=4)
test_op_biv("Laguerre t*d/dt + k*t", lambda k, f: t_s*sp.diff(f, t_s) + k*t_s*f, kmax=4)
test_op_biv("Pure deriv + affine: d/dt + k*t", lambda k, f: sp.diff(f, t_s) + k*t_s*f, kmax=4)
test_op_biv("GLZ t(1+t)d/dt + 1 + k*t", lambda k, f: t_s*(1+t_s)*sp.diff(f, t_s) + (1 + k*t_s)*f, kmax=4)
