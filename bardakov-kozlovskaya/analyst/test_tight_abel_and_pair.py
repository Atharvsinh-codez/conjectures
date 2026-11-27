"""
TEST TIGHTENED ELEMENTARY BOUNDS FOR MINOR ARCS (b_r^{pair} AND 2-STEP PAIRED ABEL)
AND CHECK THE EXACT SURVIVING MAJOR-ARC SET IN G* = {c*kap^2 < 42, c*q*kap > 0.85, kap < 2.5}!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd
from scipy.integrate import quad

# =============================================================================
# 1. EXACT UNIVERSAL BOUND ON |B_r - A_r| / c FOR tau_r = r*c*kap >= 0.85!
# =============================================================================
# Recall: for ANY r >= 1 with cos(2*pi*r*nu*m) = 1 (worst case, q | r):
#   B_r = 2 * sum_{m>=1} exp(-tau_r * (m/c) * (1 + m/c))
#   A_r =     sum_{m=1}^{c-1} exp(-tau_r * (m/c) * (1 - m/c))
# Notice:
#   - Continuum value: I(tau) := 2*int_0^inf e^{-tau u(1+u)} du - int_0^1 e^{-tau u(1-u)} du!
#   - Wait: EVEN BETTER! For general nu (when cos(2*pi*r*nu*m) varies), we bound |B_r - A_r|
#     by pairing m = 1..floor((c-1)/2) (using |e^{-tau u(1-u)} - e^{-tau u(1+u)}|!) PLUS interior tail m >= c/2!
#   Wait: CHECK THAT! If we pair WITHOUT cos, then:
#     sum_{m=1}^{floor((c-1)/2)} 2 * (e^{-tau u(1-u)} - e^{-tau u(1+u)}) + sum_{m >= c/2} 2 * e^{-tau u(1+u)}!
#     Wait: at tau = 1, the first sum is +2*int_0^{0.5}(e^{-u(1-u)} - e^{-u(1+u)})du = +0.28,
#     and the second sum is +2*int_{0.5}^inf e^{-u(1+u)} du = +0.52 -> TOTAL 0.80 (still 3x better than 2.2*sqrt(1/tau)!)
#     AND when q | r (D = 0, cos = +1), the two pieces have OPPOSITE SIGNS (+0.52 - 0.28 = +0.24!!)!
print("=== 1. EXACT TRIANGLE-PAIR BOUND vs EXACT (cos=1) VALUE OF |B_r - A_r| / c ===")
for tau in [0.5, 0.85, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0]:
    # Pair triangle bound (valid for ALL nu!):
    I_pair, _ = quad(lambda u: 2.0 * (exp(-tau*u*(1-u)) - exp(-tau*u*(1+u))), 0.0, 0.5)
    I_tail, _ = quad(lambda u: 2.0 * exp(-tau*u*(1+u)), 0.5, 0.5 + sqrt(35.0/tau))
    I_signed = I_tail - I_pair
    I_tri = I_pair + I_tail
    print(f"  tau={tau:5.2f}: signed={I_signed:+.5f} | tri={I_tri:.5f} | tri*tau^2={I_tri*tau*tau:.4f} | tri*sqrt(tau)={I_tri*sqrt(tau):.4f}")

# =============================================================================
# 2. CHECK 2-STEP PAIRED COSINE ABEL BOUND (4th-order difference -> 1 / sin^4(pi*r*nu)!)
# =============================================================================
# Recall: for ANY sequence f(m) (m >= 0) with f(0) = 0:
#   Step 1: 4*sin^2(th/2) * sum_{m>=1} f(m) cos(m*th) = f(1) - sum_{m>=1} d2f(m) cos(m*th),
#           where d2f(m) = f(m+1) - 2*f(m) + f(m-1).
#   Step 2: Apply the SAME identity to h(m) = d2f(m):
#           4*sin^2(th/2) * sum_{m>=1} h(m) cos(m*th) = -h(1)*(1 - cos(th)) + (h(1) - h(2))*cos(th)
#                                                       - sum_{m>=2} d2h(m) cos(m*th)!!
#   Divide by 16*sin^4(th/2):
#     sum_{m>=1} f(m) cos(m*th) = [ f(1) + h(1)/2 ] / (4*sin^2(th/2))
#                               + [ (h(2) - h(1))*cos(th) + sum_{m>=2} d4f(m) cos(m*th) ] / (16*sin^4(th/2))!!
#   WAIT A SECOND!! CHECK f(1) + h(1)/2 FOR f(m) = U_r(m) - V_r(m)!!
#     f(1) = f(1),  h(1) = f(2) - 2*f(1),  SO f(1) + h(1)/2 = f(2)/2!!
#   WAIT A MINUTE! CHECK THE SIGN OF h(1)*(1 - cos th)!
#   Let's verify the EXACT algebraic identity for 16*sin^4(th/2) * sum_{m>=1} f(m) cos(m*th) right now!
print("\n=== 2. VERIFY EXACT 2-STEP PAIRED COSINE ABEL IDENTITY & CONSTANT ===")
def check_2step_abel(c, jm, kap, r):
    nu = jm / c
    th = 2.0 * pi * r * nu
    s2 = 4.0 * (sin(th / 2.0) ** 2)
    if s2 < 1e-10: return None
    alpha = kap / c
    M = int(np.ceil(sqrt(40.0 / max(r * alpha, 1e-10)))) + 10
    # Construct EXACT paired sequence f(m) for m = 0..M such that B_r - A_r = 2*sum_{m>=1} f(m) cos(m*th) + mid!
    # Wait:
    #   A_r = sum_{m=1}^{c-1} exp(-r*alpha*m*(c-m)) cos(m*th)
    #   B_r = 2 * sum_{m=1}^{M} exp(-r*alpha*m*(c+m)) cos(m*th)
    # Wait: by symmetry m <-> c - m in A_r:
    #   cos((c - m)*th) = cos(2*pi*r*j - m*th) = cos(m*th) (EXACT for rational nu = j/c!!)
    #   So A_r = 2 * sum_{m=1}^{floor((c-1)/2)} exp(-r*alpha*m*(c-m)) cos(m*th) + [mid if c even]!
    # Even simpler (NO midpoint kink at m = c/2!!):
    #   KEEP A_r AS THE FULL SMOOTH SUM m = 1..c-1, AND NOTICE:
    #   sum_{m=1}^{c-1} exp(-r*alpha*m*(c-m)) cos(m*th) = (1/2) sum_{m=1}^{c-1} [exp(-r*alpha*m*(c-m)) + exp(-r*alpha*(c-m)*m)] cos(m*th)!
    #   On m = 1..c-1, f_e(m) = exp(-r*alpha*m*(c-m)) is SYMMETRIC and SMOOTH with NO kink at c/2,
    #   and f_e(0) = f_e(c) = 1!!
    #   Wait: f_e(m) has f_e(0) = 1, f_i(m) = exp(-r*alpha*m*(c+m)) has f_i(0) = 1!
    #   So setting f_e(m) = 0 for m >= c creates a kink at m = c? NO, at m = c: f_e(c) = 1 matches f_i(-c) = 1!
    #   OH!! MY!! GOD!!
    #   IN B_r = sum_{m=1}^M f_i(m) cos(m*th) + sum_{m=1}^M f_i(m) cos((c+m)*th):
    #     THE SECOND HALF OF B_r AT m' = c + m MEETS THE RIGHT END OF A_r AT m = c!!
    #   SO THE *SINGLE* DIFFERENCE F(m) = f_i(m) - f_e(m) (for m = 0..c) AND f_i(m) (for m > c)
    #   NO: JUST DEFINE f(m) = f_i(m) - f_e(m) for 0 <= m <= c/2, and wait: at m = c/2, f_e'(c/2) = 0!
    #   Check exact identity directly on F_e(m) = exp(-r*alpha*m*(c-m)) for m=0..c and F_i(m) = exp(-r*alpha*m*(c+m)) for m=0..M:
    ms_e = np.arange(0, c + 1, dtype=np.float64)
    Fe = np.exp(-r * alpha * ms_e * (c - ms_e))
    ms_i = np.arange(0, M + 1, dtype=np.float64)
    Fi = np.exp(-r * alpha * ms_i * (c + ms_i))
    # Exact B_r - A_r = 2 * sum_{m=1}^{M} Fi[m] cos(m*th) - sum_{m=1}^{c-1} Fe[m] cos(m*th)
    exact = (2.0 * float(np.sum(Fi[1:] * np.cos(ms_i[1:] * th)))
             - float(np.sum(Fe[1:c] * np.cos(ms_e[1:c] * th))))

    # Apply 2-step Abel to 2*Fi (m=1..M) AND -Fe (m=1..c-1) SEPARATELY, THEN CANCEL BOUNDARY AT m=0 AND m=c!!
    # Check exact 1-step (order-2) identity for finite sum S = sum_{m=1}^{N-1} g(m) cos(m*th) with g(0)=g(N)=1:
    #   (2 - 2*cos(th)) * S = -g(0) - g(N)*cos(N*th) + g(1)*... wait:
    #   (2 - 2*cos(th)) * sum_{m=1}^{N-1} g(m) cos(m*th)
    #     = g(1) + g(N-1)*cos(N*th) - g(0)*cos(th)... let's check directly:
    #   Expand (2 - 2*cos th) * sum_{m=1}^{N-1} g(m) cos(m th)
    #     = sum_{m=1}^{N-1} [2 g(m) - g(m-1) - g(m+1)] cos(m th)  (with g(0), g(N) boundary terms!)
    #     + g(0)*cos(th) - g(1) + g(N)*cos((N-1)*th) - g(N-1)*cos(N*th)!!
    #   Check! Since cos(c*th) = 1 and cos((c-1)*th) = cos(th) and Fe[0]=Fe[c]=1, Fe[1]=Fe[c-1]:
    #     For Fe: bdy_e = 2*cos(th) - 2*Fe[1] = 2*(cos(th) - Fe[1])!
    #     For 2*Fi: bdy_i = 2*(cos(th) - Fi[1])!
    #   SUBTRACT: bdy_i - bdy_e = 2*(Fe[1] - Fi[1])!! (THE 2*cos(th) CANCELS IDENTICALLY!!)
    d2Fe = Fe[2:] - 2.0 * Fe[1:-1] + Fe[:-2]  # length c - 1, indices m = 1..c-1
    d2Fi = Fi[2:] - 2.0 * Fi[1:-1] + Fi[:-2]  # length M - 1, indices m = 1..M-1
    rhs1 = (2.0 * (Fe[1] - Fi[1]) - 2.0 * float(np.sum(d2Fi * np.cos(ms_i[1:-1] * th)))
            + float(np.sum(d2Fe * np.cos(ms_e[1:-1] * th)))) / s2
    err1 = abs(exact - rhs1)

    # NOW APPLY A SECOND ORDER-2 ABEL STEP TO:
    #   T := -2 * sum_{m=1}^{M-1} d2Fi[m] cos(m*th) + sum_{m=1}^{c-1} d2Fe[m] cos(m*th)!!
    # Note: d2Fe is SYMMETRIC (d2Fe[m] = d2Fe[c-m]) and d2Fe[0] := 2*Fe[1] - 2*Fe[0]?
    # Define natural smooth extensions at m = 0:
    #   h_i(m) = Fi(m+1) - 2 Fi(m) + Fi(m-1) for m >= 1, and h_i(0) := 2*(Fi(1) - 1) [even extension]?
    #   NO! Just extend Fi(-1) := exp(-r*alpha*(-1)*(c-1)) = Fe(1)!!
    #   WAIT A SECOND!! STOP EVERYTHING!!
    #   Fi(m) = exp(-r*alpha*m*(c+m)) AT m = -1 IS *LITERALLY* Fe(1) = exp(-r*alpha*1*(c-1))!!
    #   SO d2Fi(0) = Fi(1) - 2*Fi(0) + Fi(-1) = Fi(1) - 2 + Fe(1) = d2Fe(0) IDENTICALLY!!
    #   WAIT A SECOND!! CHECK THAT!!
    #   Extend Fe(-1) := Fi(1) (the EXACT SAME analytic function exp(-r*alpha*m*(c-m)) at m = -1)!!
    #   THEN:
    #     h_i(0) = Fi(1) - 2 + Fe(1)
    #     h_e(0) = Fe(1) - 2 + Fi(1) = h_i(0)  IDENTICALLY!!
    #   SO THE cos(th) BOUNDARY TERM (2*h_i(0) - 2*h_e(0))*cos(th) = 0 VANISHES *IDENTICALLY*!!
    h0 = Fi[1] - 2.0 + Fe[1]
    he_ext = np.concatenate([[h0], d2Fe, [h0]]) # indices m = 0..c
    hi_ext = np.concatenate([[h0], d2Fi, [0.0]]) # indices m = 0..M
    d4Fe = he_ext[2:] - 2.0 * he_ext[1:-1] + he_ext[:-2] # m = 1..c-1
    d4Fi = hi_ext[2:] - 2.0 * hi_ext[1:-1] + hi_ext[:-2] # m = 1..M-1
    # Check identity: (2 - 2*cos(th)) * T = 2*(d2Fi[0] - d2Fe[0]) + 2*sum d4Fi cos - sum d4Fe cos!
    # Wait: check signs!
    #   (2 - 2*cos th) * sum_{m=1}^{c-1} d2Fe[m-1] cos(m th) = 2*(h0*cos th - d2Fe[0]) - sum d4Fe cos
    #   (2 - 2*cos th) * 2*sum_{m=1}^{M-1} d2Fi[m-1] cos(m th) = 2*(h0*cos th - d2Fi[0]) - 2*sum d4Fi cos
    # So (2 - 2*cos th) * T = 2*(d2Fi[0] - d2Fe[0]) + 2*sum d4Fi cos - sum d4Fe cos!!
    bdy1 = 2.0 * (Fe[1] - Fi[1])
    bdy2 = 2.0 * (d2Fi[0] - d2Fe[0])
    d4_sum = 2.0 * float(np.sum(d4Fi * np.cos(ms_i[1:-1] * th))) - float(np.sum(d4Fe * np.cos(ms_e[1:-1] * th)))
    rhs2 = bdy1 / s2 + (bdy2 + d4_sum) / (s2 * s2)
    err2 = abs(exact - rhs2)
    # Rigorous triangle upper bound on |B_r - A_r| / r:
    ub2 = (abs(bdy1) / s2 + (abs(bdy2) + 2.0*float(np.sum(np.abs(d4Fi))) + float(np.sum(np.abs(d4Fe)))) / (s2 * s2)) / r
    # Compare with 1-step bound:
    ub1 = (abs(bdy1) + 2.0*float(np.sum(np.abs(d2Fi))) + float(np.sum(np.abs(d2Fe)))) / (s2 * r)
    return abs(exact) / r, ub1, ub2, err1, err2

for c, jm, kap, r in [
    (100, 9, 0.002, 11), (100, 9, 0.002, 1), (100, 3, 0.005, 1), (100, 3, 0.02, 1),
    (50, 1, 0.02, 1), (50, 5, 0.02, 1), (30, 1, 0.05, 1), (30, 4, 0.05, 1),
    (20, 1, 0.10, 1), (10, 1, 0.15, 1), (6, 1, 0.25, 1), (200, 323%200, 0.005, 1),
]:
    ex, u1, u2, e1, e2 = check_2step_abel(c, jm, kap, r)
    print(f"  c={c:3d} jm={jm:2d} kap={kap:.3f} r={r:2d} (tau={r*c*kap:4.2f}): exact={ex:.4f} | ub1={u1:.4f} | ub2={u2:.4f} | id_err=({e1:.1e},{e2:.1e})")
