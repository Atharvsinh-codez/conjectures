import numpy as np
import sympy as sp

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0]*n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def poly_deriv(p):
    return [i * p[i] for i in range(1, len(p))] if len(p) > 1 else [0]

def apply_S(p, L):
    d = len(p) - 1
    res = [0]*(d + 2)
    for n in range(d + 2):
        v = 0
        if n <= d: v += (n + 1) * p[n]
        if 1 <= n <= d + 1: v += (L - n + 1) * p[n - 1]
        res[n] = v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_P(p, c):
    res = [(c - n) * p[n] for n in range(len(p))]
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def roots_real(p):
    c = np.array(p[::-1], dtype=np.float64)
    r = np.roots(c)
    assert np.max(np.abs(r.imag)) < 1e-6, f"Complex root: {r}"
    return np.sort(r.real)

def check_order_equal_deg(g, h):
    """Check if r_m(g) < r_m(h) < r_{m+1}(g) for all m (i.e. g << h equal degree, roots of g to the LEFT of h)."""
    rg = roots_real(g)
    rh = roots_real(h)
    assert len(rg) == len(rh)
    ok_LR = all(rg[m] < rh[m] for m in range(len(rg))) and all(rh[m] < rg[m+1] for m in range(len(rg)-1))
    ok_RL = all(rh[m] < rg[m] for m in range(len(rg))) and all(rg[m] < rh[m+1] for m in range(len(rg)-1))
    return ok_LR, ok_RL

KMAX = 14
A = {(0, 0): [1]}
for k in range(1, KMAX + 1):
    for i in range(k + 1):
        u1 = apply_S(A[(k-1, i)], 3*k - 3 - i) if i <= k - 1 else [0]
        u2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S(poly_add(u1, u2), 3*k - 2 - i)

print("=== 1. EQUAL-DEGREE INTERLACING: A_{k,i} vs P_d(A_{k,i-1}) AND vs A'_{k,i-1} ===")
fail_Pd = 0
fail_der_LR = 0
fail_der_RL = 0
total = 0
for k in range(2, 14):
    for i in range(1, k):
        g = A[(k, i)]
        f = A[(k, i-1)]
        d = len(f) - 1
        Pf = apply_P(f, d)
        df = poly_deriv(f)
        ok1_LR, ok1_RL = check_order_equal_deg(g, Pf)
        ok2_LR, ok2_RL = check_order_equal_deg(g, df)
        total += 1
        if not ok1_LR: fail_Pd += 1
        if not ok2_LR: fail_der_LR += 1
        if not ok2_RL: fail_der_RL += 1
print(f"Total pairs (k=2..13): {total}")
print(f"  A_{{k,i}} << P_d(A_{{k,i-1}}) (g left of P_d f): fails = {fail_Pd}")
print(f"  A_{{k,i}} << A'_{{k,i-1}}     (g left of f'):     fails = {fail_der_LR}")
print(f"  A'_{{k,i-1}} << A_{{k,i}}     (f' left of g):     fails = {fail_der_RL}")

print("\n=== 2. WAIT: WHAT ABOUT A_{k,i} vs D_m-DERIVATIVE OR (c - x d/dx) FOR SMALLER c < d? ===")
# Wait: for c < d, P_c f has a POSITIVE root (since leading coeff (c-d)a_d < 0, constant term c a_0 > 0)!
# Wait: what about A_{k,i} vs (f' + alpha * P_d f)? That's just (alpha*d + (1 - alpha*x)?) no, x f' + P_d f = d f!
# What about linear combination of f' and P_d f:
# Wait: f'(x) = sum n a_n x^{n-1} has deg d-1.
# P_d f(x) = sum (d-n) a_n x^n has deg d-1.
# Look at their roots: r(f') < r(P_d f)!
# Where does r(A_{k,i}) sit between r(f') and r(P_d f)?!
for k in [3, 5, 8, 11]:
    for i in [1, k//2, k-1]:
        g = A[(k, i)]
        f = A[(k, i-1)]
        d = len(f) - 1
        rg = roots_real(g)
        rdf = roots_real(poly_deriv(f))
        rpf = roots_real(apply_P(f, d))
        # compute relative position (rg - rdf)/(rpf - rdf)
        pos = (rg - rdf) / (rpf - rdf)
        print(f"  k={k:2d}, i={i:2d}: pos of r(A_{{k,i}}) in [r(f'), r(P_d f)]: min={pos.min():+.3f}, max={pos.max():+.3f}")

print("\n=== 3. CHECK GENERALIZED OPERATOR IDENTITY FOR (1-x)^a x^b d/dx WITH S_L! ===")
# Let's ask SymPy: what first-order differential operators Q = (p0 + p1*x) + (q0 + q1*x + q2*x^2)*d/dx
# satisfy Q * S_L = S_{L'} * Q' (or Q * S_L = S_{L-1} * Q)!!
x, L, Lp = sp.symbols('x L Lp')
a0, a1, b0, b1, b2 = sp.symbols('a0 a1 b0 b1 b2')
c0, c1, d0, d1, d2 = sp.symbols('c0 c1 d0 d1 d2')
f = sp.Function('f')(x)

def S_op(expr, Lv):
    return (1 + Lv*x)*expr + x*(1 - x)*sp.diff(expr, x)

def Q_op(expr, p0, p1, q0, q1, q2):
    return (p0 + p1*x)*expr + (q0 + q1*x + q2*x**2)*sp.diff(expr, x)

# Solve Q * S_L(f) == S_{Lp}(Q'(f)) for ALL f (collect coeffs of f, f', f'')
lhs = sp.expand(Q_op(S_op(f, L), a0, a1, b0, b1, b2))
rhs = sp.expand(S_op(Q_op(f, c0, c1, d0, d1, d2), Lp))
diff = sp.expand(lhs - rhs)
eqs = []
for deriv in [sp.diff(f, x, 2), sp.diff(f, x), f]:
    poly = sp.Poly(diff.coeff(deriv), x)
    eqs.extend(poly.coeffs())
# Fix Lp = L - 1 (or keep Lp free!)
sol_shift1 = sp.solve([e.subs(Lp, L - 1) for e in eqs], [a0, a1, b0, b1, b2, c0, c1, d0, d1, d2], dict=True)
print("Solutions to Q * S_L = S_{L-1} * Q' (1st order):")
for s in sol_shift1:
    print(" ", s)

sol_shift0 = sp.solve([e.subs(Lp, L) for e in eqs], [a0, a1, b0, b1, b2, c0, c1, d0, d1, d2], dict=True)
print("Solutions to Q * S_L = S_L * Q' (1st order):")
for s in sol_shift0:
    print(" ", s)
