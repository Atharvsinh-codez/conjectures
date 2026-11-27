import mpmath as mp
import sympy as sp
import math

# Quick 5-second check:
# 1. Check DOUBLE-operator at K = 800, i0 = 755 (w_eff = 16.4, deg = 845? Fast: only compute 6 polys!)
#    Wait: D_2^m(1) up to m = 850 takes ~2 seconds in mpmath! Let's do K = 700, i0 = 660 (w = 16.1, deg = 740, ~1s at dps=120)
# 2. Check SINGLE-operator K = 1000, i0 = 960: find exact integer b3 giving COMPLEX ROOTS of P3[961]
#    by scanning the sign-change count of P3[961] on ALL (-1, 0) (or checking where P3[961](t) has local min > 0)!

mp.mp.dps = 120

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

# Part 1: Find b3 making P3[961] = D_{b3} P2[961] + P2[960] have COMPLEX ROOTS!
# Notice P2[960] has roots r0 = -5.5216e-14, r1 = -1.8189e-13, r2 = -7.4132e-13.
# P2[961] has roots s0 = -7.2463e-14, s1 = -1.8089e-13, s2 = -1.1928e-12.
# Wait: on (r1, s1) = (-1.8189e-13, -1.8089e-13), P2[960] < 0 and P2[961] > 0 (check!).
# At t = -1.814e-13 (inside that tiny gap!): P2[960] < 0, P2[961] > 0, P2'[961] < 0!
# For D_{b3} P2[961](t) + P2[960](t) to be < 0 on BOTH sides of s1 AND at s1 (killing the 2 roots near -1.81e-13!):
# At s1 = -1.8089e-13, D_{b3} P2[961](s1) = s1*(1+s1)*P2'[961](s1) (INDEPENDENT of b3!)
# Check the sign of P2[960](s1) vs s1*(1+s1)*P2'[961](s1)!!
K = 1000; i0 = 960
D2 = [[1]]
for m in range(1, 45): D2.append(apply_D_int(2, D2[-1]))
cs = {r: math.comb(K, r) for r in range(i0 - 3, i0 + 3)}
g = math.gcd(*cs.values())
P = {r: [ (cs[r] // g) * c for c in D2[K - r] ] for r in range(i0 - 3, i0 + 3)}
b1 = 196200000000000; b2 = 181330000000000
P1 = {r: add_int(apply_D_int(b1, P[r]), P[r-1]) for r in range(i0 - 2, i0 + 3)}
P2 = {r: add_int(apply_D_int(b2, P1[r]), P1[r-1]) for r in [i0 - 1, i0, i0 + 1]}
p2_60 = [mp.mpf(c) for c in P2[i0]]
p2_61 = [mp.mpf(c) for c in P2[i0+1]]
# Check P3[961] at b3 = 2, 1e12, 1e14, 1e15, 1e16 around t = -1.81e-13!
s1 = mp.findroot(lambda x: eval_mp(p2_61, x), (mp.mpf('-1.815e-13'), mp.mpf('-1.80e-13')), solver='anderson')
r1 = mp.findroot(lambda x: eval_mp(p2_60, x), (mp.mpf('-1.825e-13'), mp.mpf('-1.815e-13')), solver='anderson')
print(f"s1 (root of P2[961]) = {float(s1):.8e}, r1 (root of P2[960]) = {float(r1):.8e}")
print(f"At s1: P2[960](s1) = {float(eval_mp(p2_60, s1)):.4e}, s1*P2'[961](s1) = {float(s1*eval_mp(deriv_mp(p2_61), s1)):.4e}")
# Notice BOTH have the SAME sign (negative!)! So P3[961](s1) is NEGATIVE for ALL b3!
# And on (-1.82e-13, -1.80e-13), P2[960] > 0 to the left of r1, < 0 to the right of r1;
# D_{b3} P2[961] is NEGATIVE at s1 and has a root to the left of s1; if b3 is chosen so D_{b3} P2[961] stays NEGATIVE at r1 too,
# let's check P3[961] on [-2.0e-13, -1.6e-13] for various b3!
for b3 in [2, 10**12, 5*10**12, 10**13, 3*10**13, 10**14, 10**15, 10**16]:
    p3 = [mp.mpf(c) for c in add_int(apply_D_int(b3, P2[i0+1]), P2[i0])]
    # Count sign changes on [-3e-13, -1e-13] (100 points)
    xs = [mp.mpf('-3e-13') + mp.mpf(m)*mp.mpf('2e-15') for m in range(101)]
    vals = [eval_mp(p3, x) for x in xs]
    sc = sum(1 for m in range(100) if vals[m]*vals[m+1] < 0)
    print(f"  b3={b3:18d}: sign_changes in [-3e-13, -1e-13] = {sc}, max_val={float(max(vals)):+.3e}, min_val={float(min(vals)):+.3e}")

# Part 2: Check DOUBLE-operator at K = 700, i0 = 660 (w_eff = 16.1)!
# Fast: compute ONLY D_2^m(1) using mp.mp.dps = 120 (or compute exact asymptotic formula for DOUBLE-operator!)
print("\n=== PART 2: EXACT ASYMPTOTIC FORMULA FOR DOUBLE-OPERATOR M_i^{(w)} ===")
# Compute exact polynomial in Phi, w for DOUBLE-operator M_i^{(w)} using SymPy!!
Phi, w, a, t = sp.symbols('Phi w a t')
# Let p_m(t) = 1 + 2^{2-m} s t (for m = 0, 1, 2, 3), with s*t = -a (so t = -a/s, b*t = Phi - 1, b = (1-Phi)*s/a).
# Compute q_m(t) = D_b(D_b p_m + w p_{m-1}) to O(a) as s -> inf (a fixed small, b*t = Phi - 1 fixed)!
# In variable u = s*t (so t = u/s, d/dt = s d/du, t d/dt = u d/du, 1+t -> 1, 1 + b t = 1 + (1-Phi)/a * (-u)):
# Wait: simpler! Set u = -s*t (so u = a at t = -eps, and 1 + b t = 1 - (1-Phi)*u/a)!
# Then D_b f(u) = (1 - (1-Phi)*u/a) f(u) + u f'(u) !! (Exact scaling limit as s -> inf with u = -s*t fixed!!)
u = sp.Symbol('u')
B_val = (1 - Phi) / a
def D_lim(f):
    return sp.expand((1 - B_val * u) * f + u * sp.diff(f, u))

# Start from p_m(u) = 1 - 2^{2-m} * u (for m = 0, 1, 2, 3 -> coeffs 4, 2, 1, 1/2):
pm = {m: 1 - sp.Rational(4, 2**m) * u for m in [0, 1, 2, 3]}
# Wait: for q_m = D_lim(D_lim(pm[m]) + w * pm[m-1]) (m = 1, 2, 3):
qm = {m: sp.expand(D_lim(D_lim(pm[m]) + w * pm[m-1])) for m in [1, 2, 3]}
# Compute M_lim(u=a) = u * (q_3 q_2' - q_3' q_2) + w * (q_2^2 - q_1 q_3) evaluated at u = a, to O(a)!
# Wait: check sign of W in u = -s*t! d/dt = -s d/du, (-t) = u/s, so (-t) W_t(q3, q2) = (u/s)(-s)(q3 q2_u - q3_u q2) = -u(q3 q2_u - q3_u q2)!
W_term = sp.expand(-u * (qm[3]*sp.diff(qm[2], u) - sp.diff(qm[3], u)*qm[2])).subs(u, a)
T_term = sp.expand(w * (qm[2]**2 - qm[1]*qm[3])).subs(u, a)
M_dbl = sp.series(W_term + T_term, a, 0, 2).removeO()
M_dbl_c1 = sp.factor(sp.expand(M_dbl / a))
print("  DOUBLE-operator M_i^{(w)} / alpha (exact leading O(alpha) coeff in Phi, w) =")
print("   ", sp.collect(sp.expand(M_dbl_c1), Phi))
# Evaluate minimum over Phi for w = 1, 5, 10, 15, 20, 30, 50!
for wv in [1, 5, 10, 15, 20, 30, 50, 100]:
    poly_Phi = sp.Poly(M_dbl_c1.subs(w, wv), Phi)
    # Find min over real Phi via roots of derivative
    dpoly = poly_Phi.diff()
    crit = [complex(r).real for r in sp.nroots(dpoly) if abs(complex(r).imag) < 1e-8]
    mn = min(float(poly_Phi.eval(c)) for c in crit)
    print(f"    w={wv:3d}: min_Phi (M_dbl/alpha) = {mn:+14.4f}  (at Phi={[round(c,2) for c in crit]})")
