"""
CRITICAL PROOF-ENGINEERING TEST:
Test the EXACT Vandermonde-Cauchy integral in (y, z) (Chart I) AND the dual e2=1 chart (Chart II)
across ALL admissible (a, b, c) for n = 24 and n = 36!

1. Chart I (Vandermonde-Cauchy in (y, z) with x = 1):
   c_{a,b,c}^{(n)} = (1/(2pi)^2) int int Q(y, z) * W(y, z) d(phi2) d(psi)
   where:
     Q(y, z) = prod_{k=0}^{n-1} (z - (1 + zeta^k y^{1/n})^n) / (y^{k2} z^{k3})
     W(y, z) = (1 - y)(1 - z)(1 - z/y) / ((1 + y + z)^{a+1} (1 + z/y + z)^{b+1})
   Valid WHENEVER |y| + |z| < 1 and |z/y| + |z| < 1!
   Wait: at saddle radii Ry = r2^n, Rz = r3^n (r2 = sin(pi*k2/n)/sin(pi*k1/n), r3 = sin(pi*k3/n)/sin(pi*k1/n)):
   - First test the PURE MONOMIAL integrand Q(y, z) (which gives M_{k1,k2,k3} = [x^{k1} y^{k2} z^{k3}] p_n!)
     Wait: Q(y, z) has NO poles anywhere on y != 0, z != 0, so it works for EVERY (k1, k2, k3) INCLUDING k1=k2 (a=0) and k2=k3 (b=0)!!
     TEST IF (-1)^{k1*k2+k1*k3+k2*k3} * Re(Q(r2^n e^{i phi2}, r3^n e^{i psi})) >= -eps * avg(Q)
     FOR EVERY SINGLE (k1, k2, k3) AT n = 24 AND n = 36 (INCLUDING ALL CORNERS AND BOUNDARIES)!!
   - Second test the FULL integrand Q(y, z) * W(y, z) wherever Ry < 0.5 and Rz/Ry < 0.5!
   - Third test Chart II (small a: e2 = 1, e1 = R1*e^{i*phi1}, e3 = R3*e^{i*phi3}) on ALL small-a points
     especially (1, 10, 1) and (1, 7, 3) at n = 24!
"""
import mpmath as mp
import sys
sys.path.insert(0, r'd:\solving-conjectures\work\bardakov-kozlovskaya\verify')
import engine

mp.mp.dps = 50

def eval_Q_monomial(n, k2, k3, y, z):
    beta = y**(mp.mpf(1)/n)
    prod = mp.mpc(1)
    for k in range(n):
        prod *= (z - (1 + mp.exp(2*mp.pi*1j*k/n)*beta)**n)
    return prod / (y**k2 * z**k3)

print("=== PART 1: Test PURE MONOMIAL integrand Q(y, z) at Maillot saddle radii for ALL (a,b,c) at n=24 ===")
n = 24
E = engine.compute_E_table(n)
C = engine.compute_elementary_coeffs(n)

M_grid = 10
worst_mono_ratio = 10.0
worst_mono_minre = 10.0
worst_mono_info = None
for (a, b, c), A_val in sorted(C.items()):
    if c == 0: continue
    k1, k2, k3 = a+b+c, b+c, c
    M_val = ((-1)**(n - k3)) * E[n - k3][k2]
    sgn = 1 if M_val > 0 else -1
    assert sgn == (-1)**(k1*k2 + k1*k3 + k2*k3)
    r2 = mp.sin(mp.pi * k2 / n) / mp.sin(mp.pi * k1 / n)
    r3 = mp.sin(mp.pi * k3 / n) / mp.sin(mp.pi * k1 / n)
    Ry, Rz = r2**n, r3**n
    vals = []
    for j1 in range(M_grid):
        ph1 = 2*mp.pi*j1/M_grid
        y = Ry * mp.exp(1j*ph1)
        for j2 in range(M_grid):
            ph2 = 2*mp.pi*j2/M_grid
            z = Rz * mp.exp(1j*ph2)
            v = sgn * eval_Q_monomial(n, k2, k3, y, z) / abs(M_val)
            vals.append(v)
    avg = sum(vals)/len(vals)
    abs_avg = sum(abs(v) for v in vals)/len(vals)
    rat = float(abs(avg)/abs_avg)
    min_re = min(float(v.real) for v in vals)
    if rat < worst_mono_ratio:
        worst_mono_ratio = rat
        worst_mono_info = (a, b, c, rat, min_re, float(avg.real))
    if min_re < worst_mono_minre:
        worst_mono_minre = min_re
print(f"n=24 ALL {sum(1 for (a,b,c) in C if c>0)} (c>=1) monomials: worst |avg|/avg|Q| = {worst_mono_ratio:.4f} at {worst_mono_info}, global min Re(Q/avg) = {worst_mono_minre:+.4f}")

print("\n=== PART 2: Test Chart II (e2 = 1) on small-a points at n=24 (esp. (1,10,1) and (1,7,3)) ===")
# In Chart II: fix e2 = 1, vary e1 = R1 * e^{i*phi1}, e3 = R3 * e^{i*phi3}!
# Wait: what are the Maillot saddle radii for (e1, e2, e3) when e2 = 1?!
# For k1 >= k2 >= k3: x ~ r1^n, y ~ r2^n, z ~ r3^n with r1 = sin(t1), r2 = sin(t2), r3 = sin(t3).
# Normalize by e2 = x*y + ... -> scale so |x*y| = 1, i.e., divide all radii by sqrt(r1*r2)!
# Then:
#   |x| = (r1/r2)^{n/2}, |y| = (r2/r1)^{n/2}, |z| = (r3/sqrt(r1*r2))^n!
#   Wait: if r1 ~ r2 (small a!), then |x| ~ 1, |y| ~ 1, and e1 = x + y + z can be SMALL if x ~ -y (i.e. e2 = xy ~ -1)!
#   Wait: if e2 = -1, x ~ 1, y ~ -1, then e1 = x+y+z is small!
#   More simply, in (e1, e3) with e2 = 1:
#   R1 = |e1| ~ exp(d/da log|c|)? Check empirically: R1 ~ |c(a,b,c) / c(a+2,b-1,c)|^{-1/2} or R1 ~ (r1/r2)^{n/2}?!
#   Wait: at a=1, b=10, c=1 (n=24): c_{1,10,1} has e1^1 e2^{10} e3^1!
#   Notice that as a polynomial in e1, e1^1 is the LOWEST power of e1 (only a=1, 3, 5, ... appear for odd c=1)!
#   So taking R1 -> 0 gives EXACTLY the coefficient of e1^1 (with ZERO interference from a=3,5,...)!!
#   OH!! MY!! GOD!!
#   FOR a = 0: P(0, 1, e3) gives a=0 EXACTLY (1D integral in e3 only!)!
#   FOR a = 1: d/de1 P(0, 1, e3) gives a=1 EXACTLY (1D integral in e3 only!)!
#   CHECK THAT RIGHT NOW!!
def eval_Pn_from_e1e3(n, e1, e3):
    roots = mp.polyroots([mp.mpc(1), mp.mpc(-e1), mp.mpc(1), mp.mpc(-e3)], maxsteps=200, extraprec=60)
    x, y, z = roots
    alpha, beta = x**(mp.mpf(1)/n), y**(mp.mpf(1)/n)
    prod = mp.mpc(1)
    for k in range(n):
        prod *= (z - (alpha + mp.exp(2*mp.pi*1j*k/n)*beta)**n)
    return prod

for a, b, c in [(1, 10, 1), (1, 7, 3), (1, 4, 5), (1, 1, 7), (2, 8, 2), (2, 5, 4), (3, 9, 1)]:
    exact_c = C[(a, b, c)]
    sgn = 1 if exact_c > 0 else -1
    k1, k2, k3 = a+b+c, b+c, c
    r1, r2, r3 = mp.sin(mp.pi*k1/n), mp.sin(mp.pi*k2/n), mp.sin(mp.pi*k3/n)
    # Rescale so e2 = 1: R3 = (r3 / sqrt(r1*r2))^n (Wait: e3/e2^{3/2} = (r1*r2*r3 / (r1*r2)^{3/2})^n = (r3/sqrt(r1*r2))^n!)
    # And e1/e2^{1/2} = ((r1^n - r2^n) / (r1*r2)^{n/2})? Try R1 = (r1/r2)^{n/2} - (r2/r1)^{n/2} (or small R1 for a=1)!
    R3 = (r3 / mp.sqrt(r1 * r2))**n
    # For R1, try exact saddle ratio sqrt(|c(a,b,c)/c(a+2,b-1,c)|) vs formula (r1/r2)^{n/2} - (r2/r1)^{n/2}
    R1_pred = (r1/r2)**(n/2) - (r2/r1)**(n/2)
    if R1_pred <= 0: R1_pred = mp.mpf('0.1')
    # Also try small R1 = 0.01 for a=1!
    for tag, R1 in [("sinh_R1", R1_pred), ("small_R1", mp.mpf('0.05'))]:
        if tag == "small_R1" and a > 2: continue
        scale = (R1**a) * (R3**c) * abs(exact_c)
        M1 = 10
        vals = []
        for j1 in range(M1):
            ph1 = 2*mp.pi*j1/M1
            e1 = R1 * mp.exp(1j*ph1)
            for j3 in range(M1):
                ph3 = 2*mp.pi*j3/M1
                e3 = R3 * mp.exp(1j*ph3)
                v = sgn * eval_Pn_from_e1e3(n, e1, e3) * mp.exp(-1j*(a*ph1 + c*ph3)) / scale
                vals.append(v)
        avg = sum(vals)/len(vals)
        abs_avg = sum(abs(v) for v in vals)/len(vals)
        re_vals = [float(v.real) for v in vals]
        print(f"ChartII({tag:8s}) n={n} ({a},{b:2d},{c}) R1={float(R1):.4f} R3={float(R3):.2e}: avg={float(avg.real):+.5f} | Re in [{min(re_vals):+.3f}, {max(re_vals):+.3f}] | |avg|/avg|F|={float(abs(avg)/abs_avg):.4f}")
