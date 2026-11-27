"""TEST: Does (r + t) F_{a,i}(t) + alpha F_{a,i-1}(t) REALLY have complex roots for some r in (0,1), alpha > 0?
(If YES, then finite M = 10^8 gives an EXACT counterexample to arbitrary-(c_j, m_j, lam_j) Atomic Theorem!
 AND we immediately check whether c_j = 1 (so r = 0, i.e., t * F + alpha * G, which ALWAYS preserves real roots!)
 STILL SURVIVES!)
"""
import math
import numpy as np
import sympy as sp

t = sp.Symbol('t')

def apply_GLZ(c, m, poly):
    d = len(poly) - 1
    res = [sp.Rational(0)] * (d + 2)
    for r, coef in enumerate(poly):
        res[r] += coef * (r + c)
        res[r + 1] += coef * (r + m)
    return res

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [sp.Rational(0)] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def poly_mul(a, b):
    res = [sp.Rational(0)] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            res[i + j] += u * v
    return res

D2 = [[sp.Rational(1)]]
for _ in range(6): D2.append(apply_GLZ(1, 2, D2[-1]))

# Take a = 3, i = 1: Fi = 3 D2^2(1) (deg 2), Fim1 = D2^3(1) (deg 3)
Fi = [3 * v for v in D2[2]]
Fim1 = D2[3]
r_i = sorted(np.roots([float(v) for v in reversed(Fi)]).real)
r_im1 = sorted(np.roots([float(v) for v in reversed(Fim1)]).real)
print("Roots of D2^2(1):", r_i)
print("Roots of D2^3(1):", r_im1)

# Pick -r strictly between r_im1[1] and r_i[1] (inside the SAME interval as r_i[1], creating a double-interval collision!)
# Wait: at t0 = (r_im1[1] + r_i[1])/2, both (t+r) and Fi have roots in (r_im1[1], r_im1[2])!
# So in (r_im1[1], r_im1[2]), (t+r)*Fi has TWO roots (same sign at both endpoints r_im1[1], r_im1[2]!),
# while alpha*Fim1 has the SAME sign as (t+r)*Fi on (r_im1[1], r_im1[2])? Let's check BOTH signs of (r_im1[1], r_i[1]) and (r_i[0], r_im1[1])!
found = None
for r_f in np.linspace(0.05, 0.95, 19):
    r_q = sp.Rational(int(round(r_f * 100)), 100)
    for exp_a in np.linspace(-4, 4, 33):
        a_f = 10 ** exp_a
        comb_f = np.polyadd(np.polymul([1.0, float(r_q)], [float(v) for v in reversed(Fi)]),
                            [a_f * float(v) for v in reversed(Fim1)])
        roots = np.roots(comb_f)
        if np.max(np.abs(roots.imag)) > 1e-4:
            found = (r_q, a_f, float(np.max(np.abs(roots.imag))))
            break
    if found: break

print("Limit counterexample (r+t)*Fi + alpha*Fim1:", found)
if found:
    r_q, a_f, _ = found
    a_q = sp.Rational(a_f).limit_denominator(1000)
    # Verify EXACTLY in sympy for limit AND for finite M = 10^6!!
    comb_lim = poly_add(poly_mul([r_q, sp.Rational(1)], Fi), [a_q * v for v in Fim1])
    Plim = sp.Poly(list(reversed(comb_lim)), t, domain='QQ')
    print(f"  EXACT limit (r={r_q}, alpha={a_q}): deg={Plim.degree()}, real_roots={int(sp.count_roots(Plim))}")
    M = 10**6
    c_big = r_q * M
    lam_big = a_q * M
    comb_fin = poly_add(apply_GLZ(c_big, M, Fi), [lam_big * v for v in Fim1])
    Pfin = sp.Poly(list(reversed(comb_fin)), t, domain='QQ')
    print(f"  EXACT finite M=10^6 (c={c_big}, m={M}, lam={lam_big}): deg={Pfin.degree()}, real_roots={int(sp.count_roots(Pfin))}")
