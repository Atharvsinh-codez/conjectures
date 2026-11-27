"""1. Fixed attack on Htilde_3(t,y) * (r + t + alpha*y)^b (fixing the i=0 multiple-root bug!)
   AND also attack directly on bare Atomic Theorem N=5,6,7 with 500 fast trials!
2. Inspect exact algebraic structure of F_{3,1}^2 - F_{3,0} F_{3,2} for general D_1, D_2, D_3!
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

P = [[None]*(k+1) for k in range(6)]
P[0][0] = [sp.Rational(1)]
for k in range(1, 6):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_GLZ(1, m, P[k-1][0])
    for i in range(1, k):
        U[i] = poly_add(apply_GLZ(1, m, P[k-1][i]), P[k-1][i-1])
    U[k] = P[k-1][k-1]
    for i in range(k+1):
        P[k][i] = apply_GLZ(1, m, U[i])

print("=== 1. FIXED ATTACK ON Htilde_k(t, y) * prod_{s=1}^b (r_s + t + alpha_s * y) ===")
# Use DISTINCT r_s so there are NO multiple roots!
rng = np.random.default_rng(2026)
rr_suspects = []
lc_suspects = []
for k_base in [3, 4]:
    for b in [1, 2, 3, 4, 6]:
        for trial in range(120):
            # Pick distinct r_s in (0.001, 0.999) and random alpha_s
            rs = rng.uniform(0.001, 0.999, b)
            alphas = 10 ** rng.uniform(-3, 3, b)
            F = [np.array([float(v) for v in P[k_base][i]], dtype=np.float64) for i in range(k_base + 1)]
            for s_idx in range(b):
                r_v, al_v = rs[s_idx], alphas[s_idx]
                newF = [None] * (len(F) + 1)
                newF[0] = np.polymul(F[0][::-1], [1.0, r_v])[::-1]
                for i in range(1, len(F)):
                    a1 = np.polymul(F[i][::-1], [1.0, r_v])
                    a2 = al_v * F[i-1][::-1]
                    newF[i] = np.polyadd(a1, a2)[::-1]
                newF[-1] = al_v * F[-1]
                F = newF
            # Check RR for i = 1..len(F)-2
            for i in range(1, len(F) - 1):
                # scale polynomial for conditioning
                p_hi = F[i][::-1]
                p_hi = p_hi / np.max(np.abs(p_hi))
                rts = np.roots(p_hi)
                if np.max(np.abs(rts.imag)) > 1e-3:
                    rr_suspects.append((k_base, b, i, float(np.max(np.abs(rts.imag))), rs[:3], alphas[:3]))
            # Check LC on (-1, 0)
            ts = np.linspace(-0.999, -0.001, 300)
            for i in range(1, len(F) - 1):
                v0 = np.polyval(F[i-1][::-1], ts)
                v1 = np.polyval(F[i][::-1], ts)
                v2 = np.polyval(F[i+1][::-1], ts)
                diff = v1*v1 - v0*v2
                if np.min(diff) < -1e-6 * np.max(v1*v1):
                    lc_suspects.append((k_base, b, i, float(np.min(diff)/np.max(v1*v1)), rs, alphas))
                    break
print(f"  rr_suspects (im > 1e-3): {len(rr_suspects)}, sample: {rr_suspects[:3]}")
print(f"  lc_suspects (min < -1e-6): {len(lc_suspects)}, sample: {[(x[0],x[1],x[2],x[3]) for x in lc_suspects[:3]]}")

# If any lc_suspect or rr_suspect, verify FIRST one in EXACT rational arithmetic!!
if lc_suspects:
    k_base, b, i_bad, rel_val, rs, alphas = lc_suspects[0]
    rs_q = [sp.Rational(int(round(v*1000)), 1000) for v in rs]
    al_q = [sp.Rational(v).limit_denominator(1000) for v in alphas]
    F = [[sp.Rational(v) for v in P[k_base][j]] for j in range(k_base + 1)]
    for s_idx in range(b):
        r_v, al_v = rs_q[s_idx], al_q[s_idx]
        newF = [None] * (len(F) + 1)
        newF[0] = [sp.Rational(0)] * (len(F[0]) + 1)
        for m, c in enumerate(F[0]):
            newF[0][m] += c * r_v; newF[0][m+1] += c
        for j in range(1, len(F)):
            a1 = [sp.Rational(0)] * (len(F[j]) + 1)
            for m, c in enumerate(F[j]):
                a1[m] += c * r_v; a1[m+1] += c
            newF[j] = poly_add(a1, [al_v * c for c in F[j-1]])
        newF[-1] = [al_v * c for c in F[-1]]
        F = newF
    # Exact check of LC and RR at i_bad!
    P_i = sp.Poly(list(reversed(F[i_bad])), t, domain='QQ')
    print(f"  EXACT CHECK of suspect k={k_base}, b={b}, i={i_bad}: deg={P_i.degree()}, real_roots(-1,0)={int(sp.count_roots(P_i, -1, 0))}")
    def pmul(a, b):
        res = [sp.Rational(0)]*(len(a)+len(b)-1)
        for ii, u in enumerate(a):
            for jj, v in enumerate(b): res[ii+jj] += u*v
        return res
    D = poly_add(pmul(F[i_bad], F[i_bad]), [-v for v in pmul(F[i_bad-1], F[i_bad+1])])
    Pd = sp.Poly(list(reversed(D)), t, domain='QQ')
    print(f"  EXACT LC roots in (-1,0): {int(sp.count_roots(Pd, -1, 0))}, eval(-1/2)={float(Pd.eval(sp.Rational(-1,2))):.3e}")
    # Also check ALL i for RR!
    for j in range(1, len(F) - 1):
        Pj = sp.Poly(list(reversed(F[j])), t, domain='QQ')
        rc = int(sp.count_roots(Pj, -1, 0))
        if rc < Pj.degree():
            print(f"  EXACT RR FAIL AT j={j}: rc={rc} < deg={Pj.degree()}!!")
