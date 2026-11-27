import mpmath as mp
import sympy as sp
from fractions import Fraction
import math

mp.mp.dps = 150

def apply_D_int(b, p):
    res = [0]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+b)
    return res

def add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def apply_D_mp(b, p):
    res = [mp.mpf(0)]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+b)
    return res

def add_mp(a, b):
    res = [mp.mpf(0)]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def eval_mp(p, tv):
    v = mp.mpf(0)
    for c in reversed(p): v = v * tv + c
    return v

def deriv_mp(p): return [r*p[r] for r in range(1, len(p))]
def eval_W_mp(f, g, tv):
    return eval_mp(f, tv)*eval_mp(deriv_mp(g), tv) - eval_mp(deriv_mp(f), tv)*eval_mp(g, tv)

print("=== 1. FAST SEARCH FOR COMPLEX ROOTS IN SINGLE-OPERATOR (a=1, w=1, b_j >= 2) AT K=1000, i0=960 ===")
K = 1000; i0 = 960
D2 = [[1]]
for m in range(1, 45): D2.append(apply_D_int(2, D2[-1]))
cs = {r: math.comb(K, r) for r in range(i0 - 3, i0 + 3)}
g = math.gcd(*cs.values())
P = {r: [ (cs[r] // g) * c for c in D2[K - r] ] for r in range(i0 - 3, i0 + 3)}

b1 = 196200000000000
b2 = 181330000000000
P1 = {r: add_int(apply_D_int(b1, P[r]), P[r-1]) for r in range(i0 - 2, i0 + 3)}
P2 = {r: add_int(apply_D_int(b2, P1[r]), P1[r-1]) for r in [i0 - 1, i0, i0 + 1]}

# Inspect the EXACT roots of P2[960] and P2[961] near t = -1.818e-13!
# Recall P2[960] has 2 roots close together where W(P2[961], P2[960]) < 0!
# Evaluate P2[960] and P2[961] on 40 points around t = -1.818e-13!
t_c = mp.mpf('-1.818e-13')
print("  Scanning P2[960] and P2[961] around t = -1.818e-13:")
p2_60 = [mp.mpf(c) for c in P2[i0]]
p2_61 = [mp.mpf(c) for c in P2[i0+1]]
# Find the exact sign changes of p2_60 and p2_61 in [-5e-13, -1e-14]!
grid = [-mp.mpf(10)**mp.mpf(str(e)) for e in [-14.0 + 0.05*m for m in range(45)]]
r60 = []; r61 = []
for idx in range(len(grid)-1):
    a, b = grid[idx], grid[idx+1]
    if eval_mp(p2_60, a) * eval_mp(p2_60, b) < 0:
        r60.append(mp.findroot(lambda x: eval_mp(p2_60, x), (a, b), solver='anderson'))
    if eval_mp(p2_61, a) * eval_mp(p2_61, b) < 0:
        r61.append(mp.findroot(lambda x: eval_mp(p2_61, x), (a, b), solver='anderson'))
print(f"  Rightmost roots of P2[960]: {[float(r) for r in r60]}")
print(f"  Rightmost roots of P2[961]: {[float(r) for r in r61]}")

# Notice if P2[960] has 2 roots r_a < r_b with NO root of P2[961] between them (e.g. r60[0] > r60[1] > r61[0] or similar!),
# then at the midpoint t_mid between those 2 roots of P2[960], P2[960](t_mid) has some sign,
# and to lift that lobe ABOVE 0 (killing BOTH real roots -> COMPLEX PAIR!!), we pick b3 so that
# D_{b3} P2[961](t_mid) has the OPPOSITE sign to P2[960](t_mid) and LARGER magnitude!
# Check sign of P2[960] and D_{b3} P2[961] at t_mid = (r60[0] + r60[1]) / 2!
if len(r60) >= 2:
    t_mid = (r60[0] + r60[1]) / 2
    v60 = eval_mp(p2_60, t_mid)
    v61 = eval_mp(p2_61, t_mid)
    dv61 = eval_mp(deriv_mp(p2_61), t_mid)
    print(f"  At t_mid={float(t_mid):.6e}: P2[960]={float(v60):.4e}, P2[961]={float(v61):.4e}")
    # D_{b3} P2[961](t_mid) = (1 + b3*t_mid)*v61 + t_mid*(1+t_mid)*dv61.
    # For P3[961](t_mid) = D_{b3} P2[961](t_mid) + v60 to flip sign (same sign as P2[960] OUTSIDE (r60[1], r60[0])):
    # Wait: v60 is the extremal value INSIDE the lobe; we want D_{b3} P2[961](t_mid) to have OPPOSITE sign to v60 and |D| > |v60|!
    # Solve (1 + b3*t_mid)*v61 + t_mid*(1+t_mid)*dv61 = -2 * v60 -> b3!
    b3_kill = ( (-2*v60 - t_mid*(1+t_mid)*dv61)/v61 - 1 ) / t_mid
    print(f"  Targeted b3_kill = {float(b3_kill):.4e} (>= 2? {b3_kill >= 2})")
    # Check sign of P3[961] on [r60[1]*1.5, r60[0]*0.5] for b3 near b3_kill!
    for mult in [0.5, 1.0, 2.0, 5.0, 10.0, 50.0]:
        b3_try = int(max(2, float(abs(b3_kill)) * mult))
        P3_try = [mp.mpf(c) for c in add_int(apply_D_int(b3_try, P2[i0+1]), P2[i0])]
        # Check if sign changes disappeared near t_mid!
        sc_cnt = 0
        subg = [r60[1] * mp.mpf(str(1.2 - 0.02*m)) for m in range(55)]
        for idx in range(len(subg)-1):
            if eval_mp(P3_try, subg[idx]) * eval_mp(P3_try, subg[idx+1]) < 0:
                sc_cnt += 1
        print(f"    b3={b3_try:20d}: sign_changes in lobe window = {sc_cnt}")
        if sc_cnt == 0:
            # Verify with EXACT SymPy count_roots!
            t_sym = sp.Symbol('t')
            poly = sp.Poly(sum(c * t_sym**r for r, c in enumerate(add_int(apply_D_int(b3_try, P2[i0+1]), P2[i0]))), t_sym, domain='ZZ')
            nr = poly.count_roots()
            print(f"    *** EXACT SYMPY STURM AT b3={b3_try}: deg={poly.degree()}, real_roots={nr} (< deg? {nr < poly.degree()})!! ***")
            if nr < poly.degree():
                break

print("\n=== 2. FAST 150-DPS TEST OF DOUBLE-OPERATOR prod D_{b_j}(D_{b_j} + y)(1) (a=1, w=1, b_j >= 2) ===")
# Check K = 300, i0 = 275 (w_eff = 10.6) AND K = 500, i0 = 465 (w_eff = 12.9)!
for K, i0 in [(300, 275), (500, 465)]:
    D2m = [[mp.mpf(1)]]
    for m in range(1, 2*K - i0 + 5):
        D2m.append(apply_D_mp(mp.mpf(2), D2m[-1]))
    norm = mp.mpf(math.comb(K, i0))
    Pd = {r: [(mp.mpf(math.comb(K, r)) / norm) * c for c in D2m[2*K - r]] for r in range(i0 - 3, i0 + 3)}
    s2 = eval_mp(deriv_mp(Pd[i0+1]), mp.mpf(0)) / eval_mp(Pd[i0+1], mp.mpf(0))
    w_eff = i0 / (K - i0 + 1)
    # Check asymptotic q_m(-eps) = (Phi^2 + (w+1)Phi - 1) - (Phi^2 + (2w+3)Phi + (w+2)) 2^{2-m} alpha!
    # Notice A_2(Phi) has roots near -w-1 and 0, B_2(Phi) has roots near -2w-2 and -1/2!
    # So A_2 * B_2 < 0 on Phi in (-2w-2, -w-1)!! Pick Phi around -1.5 * w_eff!
    best_M = (mp.mpf('1e100'), None)
    for alpha in [0.05, 0.2, 0.5]:
        eps = mp.mpf(str(alpha)) / s2
        for mult in [0.5, 1.0, 1.3, 1.5, 1.7, 2.0]:
            Phi = mp.mpf(str(-mult * w_eff))
            b1 = (1 - Phi) / eps
            P1d = {r: apply_D_mp(b1, add_mp(apply_D_mp(b1, Pd[r]), Pd[r-1])) for r in range(i0 - 2, i0 + 3)}
            tv = -eps
            v0, v1, v2 = eval_mp(P1d[i0-1], tv), eval_mp(P1d[i0], tv), eval_mp(P1d[i0+1], tv)
            w21 = eval_W_mp(P1d[i0+1], P1d[i0], tv)
            M_val = (-tv)*(1+tv)*w21 + v1*v1 - v0*v2
            rel_M = M_val / (abs((-tv)*(1+tv)*w21) + v1*v1 + abs(v0*v2))
            if rel_M < best_M[0]: best_M = (rel_M, (alpha, float(Phi), b1, tv, P1d))
    print(f"  DOUBLE K={K}, i0={i0} (w={w_eff:.1f}): worst rel_M = {float(best_M[0]):+.6f} at alpha={best_M[1][0]}, Phi={best_M[1][1]:.1f}")
    if best_M[0] < 0:
        _, (alpha, Phi_f, b1, t_neg, P1d) = best_M
        v0, v1 = eval_mp(P1d[i0-1], t_neg), eval_mp(P1d[i0], t_neg)
        d1 = eval_mp(deriv_mp(P1d[i0]), t_neg)
        b2_star = -(v1 + t_neg*(1+t_neg)*d1 + v0) / (t_neg * v1)
        print(f"    *** DOUBLE M_i < 0! b2* = {float(b2_star):.4e} (>=2? {b2_star >= 2}) ***")
        if b2_star >= 2:
            P2d = {r: apply_D_mp(b2_star, add_mp(apply_D_mp(b2_star, P1d[r]), P1d[r-1])) for r in [i0-1, i0, i0+1]}
            # Check W of BOTH inner u2 = D_{b2*} P1d + P1d_prev AND outer D_{b2*} u2!
            U2d = {r: add_mp(apply_D_mp(b2_star, P1d[r]), P1d[r-1]) for r in [i0, i0+1]}
            wu = eval_W_mp(U2d[i0+1], U2d[i0], t_neg)
            wp = eval_W_mp(P2d[i0+1], P2d[i0], t_neg)
            print(f"    W_inner(t_neg) = {float(wu):.4e}, W_outer(t_neg) = {float(wp):.4e}")
