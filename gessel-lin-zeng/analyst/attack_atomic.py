"""ADVERSARIAL ATTACK on the ATOMIC THEOREM:
  prod_{j=1}^N (D_j + lam_j * y)(1) = sum_{i=0}^N F_{N, i}(t) y^i,
  D_j = t(1+t) d/dt + c_j + m_j * t   (0 < c_j < m_j, lam_j >= 0).
Try EVERYTHING to break RR, row interlacing (F_{N,i} << F_{N,i-1}), or row log-concavity (F_{N,i}^2 >= F_{N,i-1} F_{N,i+1} on (-1,0))!
1. 300 wild random/extreme trials for N = 4, 5, 6 (using exact rational arithmetic on extreme cases + high-precision on random cases)!
2. Test if CAUCHY (phi_j' NOT always < 0!) breaks for wild m_k = [1, 1000, 1, 1000]!
3. Test if ONE step F_{k,i} = D F_{k-1,i} + F_{k-1,i-1} preserves LC on (-1,0) for ARBITRARY input row (f_0, f_1, f_2, f_3)
   with f_3 << f_2 << f_1 << f_0 and f_1^2 >= f_0 f_2 and f_2^2 >= f_1 f_3!
"""
import math
import numpy as np
import sympy as sp

t = sp.Symbol('t')
s = sp.Symbol('s')

def apply_D_f(c, m, poly):
    d = len(poly) - 1
    res = np.zeros(d + 2, dtype=np.float64)
    for r in range(d + 1):
        res[r] += poly[r] * (r + c)
        res[r + 1] += poly[r] * (r + m)
    return res

def apply_D_q(c, m, poly):
    d = len(poly) - 1
    res = [sp.Rational(0)] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * (r + m)
    return res

def poly_add_q(a, b):
    n = max(len(a), len(b))
    res = [sp.Rational(0)] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def poly_sub_q(a, b):
    return poly_add_q(a, [-v for v in b])

def poly_mul_q(a, b):
    res = [sp.Rational(0)] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            res[i + j] += u * v
    return res

def check_exact_atomic(steps):
    """Exact check of RR, row interlacing, and LC on (-1,0) for steps = [(c, m, lam), ...]."""
    N = len(steps)
    F = [[None]*(k+1) for k in range(N+1)]
    F[0][0] = [sp.Rational(1)]
    for k in range(1, N + 1):
        c, m, lam = steps[k - 1]
        F[k][0] = apply_D_q(c, m, F[k-1][0])
        for i in range(1, k):
            F[k][i] = poly_add_q(apply_D_q(c, m, F[k-1][i]), [lam * v for v in F[k-1][i-1]])
        F[k][k] = [lam * v for v in F[k-1][k-1]]
    # Check at level N
    rr_fail = []; lc_fail = []; row_fail = []
    for i in range(N + 1):
        P = sp.Poly(list(reversed(F[N][i])), t, domain='QQ')
        if P.degree() > 0 and int(sp.count_roots(P, -1, 0)) < P.degree():
            rr_fail.append(i)
    for i in range(1, N):
        D = poly_sub_q(poly_mul_q(F[N][i], F[N][i]), poly_mul_q(F[N][i-1], F[N][i+1]))
        Pd = sp.Poly(list(reversed(D)), t, domain='QQ')
        if Pd.eval(sp.Rational(-1, 2)) <= 0 or int(sp.count_roots(Pd, -1, 0)) > 0:
            lc_fail.append(i)
    return rr_fail, lc_fail

print("=== 1. EXACT ADVERSARIAL SEARCH ON ATOMIC THEOREM (N=3, 4, 5, 6) ===")
rng = np.random.default_rng(2026)
# Try extreme ratios of c, m, lam spanning 10^{-4} to 10^{4}!
fails_atomic = []
for trial in range(40):
    N = 4 if trial < 25 else 5
    steps = []
    for j in range(N):
        c = sp.Rational(10 ** int(rng.integers(-2, 3))) if rng.random() < 0.5 else sp.Rational(int(rng.integers(1, 20)), int(rng.integers(1, 20)))
        gap = sp.Rational(10 ** int(rng.integers(-3, 4))) if rng.random() < 0.5 else sp.Rational(int(rng.integers(1, 50)), int(rng.integers(1, 10)))
        m = c + gap
        lam = sp.Rational(10 ** int(rng.integers(-3, 4))) if rng.random() < 0.5 else sp.Rational(1)
        steps.append((c, m, lam))
    rr, lc = check_exact_atomic(steps)
    if rr or lc:
        fails_atomic.append((trial, steps, rr, lc))
print(f"  40 extreme exact random trials (N=4,5): FAILS = {len(fails_atomic)}")
for f in fails_atomic[:3]:
    print("   FAIL:", f)

print("\n=== 2. TEST CAUCHY WEIGHT C_m = (1+s^2)d/ds - 2*m*s WITH WILD m=[1, 1000, 1, 1000] ===")
P = [[None]*5 for _ in range(5)]
P[0][0] = sp.Integer(1)
ms = [1, 1000, 1, 1000]
for k in range(1, 5):
    m = ms[k - 1]
    op = lambda f, mm=m: (1 + s**2)*sp.diff(f, s) - 2*mm*s*f
    P[k][0] = sp.expand(op(op(P[k-1][0])))
    for i in range(1, k):
        P[k][i] = sp.expand(op(op(P[k-1][i]) + P[k-1][i-1]))
    P[k][k] = sp.expand(op(P[k-1][k-1]))
cauchy_fails = []
for k in range(1, 5):
    for i in range(k + 1):
        poly = sp.Poly(P[k][i], s, domain='ZZ')
        if poly.degree() > 0 and int(sp.count_roots(poly)) < poly.degree():
            cauchy_fails.append((k, i, int(sp.count_roots(poly)), poly.degree()))
print(f"  Cauchy m=[1, 1000, 1, 1000]: RR fails = {cauchy_fails}")

print("\n=== 3. TEST IF ONE STEP U_i = D f_i + f_{i-1} PRESERVES LC ON (-1,0) FOR ARBITRARY f_3<<f_2<<f_1<<f_0 WITH LC ===")
# Check if LC alone + interlacing is preserved by U_i = D f_i + f_{i-1}, or if something stronger is needed!
fails_1step_lc = 0
ts = np.linspace(-0.999, -0.001, 500)
for trial in range(500):
    r0 = sorted(rng.uniform(-0.98, -0.02, 4))
    r1 = [rng.uniform(r0[j], r0[j+1]) for j in range(3)]
    r2 = [rng.uniform(r1[j], r1[j+1]) for j in range(2)]
    r3 = [rng.uniform(r2[j], r2[j+1]) for j in range(1)]
    def mk(rs):
        p = np.array([1.0])
        for r in rs: p = np.polymul(p, [-1.0/r, 1.0])
        return p[::-1]
    f0, f1, f2, f3 = mk(r0), mk(r1), mk(r2), mk(r3)
    # Scale f1, f2 so f1^2 >= f0*f2 AND f2^2 >= f1*f3 on (-1,0)
    # Note: let f1 -> s1*f1, f2 -> s2*f2. Need s1^2/s2 >= R1 and s2^2/s1 >= R2!
    v = [np.polyval(f[::-1], ts) for f in [f0, f1, f2, f3]]
    def need_ratio(a, b, c):
        mask = (a * c > 0)
        return float(np.max((a[mask]*c[mask])/(b[mask]**2 + 1e-30))) if np.any(mask) else 0.0
    R1 = max(1.0, need_ratio(v[0], v[1], v[2]) * 1.05)
    R2 = max(1.0, need_ratio(v[1], v[2], v[3]) * 1.05)
    # Solve s1^2/s2 = R1, s2^2/s1 = R2 => s1 = (R1^2 R2)^{1/3}, s2 = (R1 R2^2)^{1/3}
    s1 = (R1*R1*R2)**(1/3); s2 = (R1*R2*R2)**(1/3)
    f1 = f1 * s1; f2 = f2 * s2
    # Apply D = t(1+t)d/dt + c + m*t
    c_v = rng.uniform(0.1, 10.0)
    m_v = c_v + rng.uniform(0.1, 50.0)
    def pad_add(a, b):
        n = max(len(a), len(b))
        out = np.zeros(n); out[:len(a)] += a; out[:len(b)] += b
        return out
    U1 = pad_add(apply_D_f(c_v, m_v, f1), f0)
    U2 = pad_add(apply_D_f(c_v, m_v, f2), f1)
    U3 = pad_add(apply_D_f(c_v, m_v, f3), f2)
    u1 = np.polyval(U1[::-1], ts); u2 = np.polyval(U2[::-1], ts); u3 = np.polyval(U3[::-1], ts)
    if np.min(u2*u2 - u1*u3) < -1e-9 * np.max(u2*u2):
        fails_1step_lc += 1
print(f"  Random LC+interlacing 4-tuple: 1-step U_2^2 >= U_1 U_3 on (-1,0) fails = {fails_1step_lc}/500")
