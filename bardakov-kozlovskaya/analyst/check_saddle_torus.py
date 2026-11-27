"""
1. Fix syntax in Clausen Mahler-measure check & verify p_n(1, omega, omega^2) scaling.
2. CRITICAL ANALYTIC CHECK:
   Look at P_{2m}(1, e2, e3) on the torus (e2, e3) = (R2 * exp(i*phi2), R3 * exp(i*phi3))!
   Compute P_{2m}(1, e2, e3) DIRECTLY from the roots (x, y, z) of T^3 - T^2 + e2*T - e3 = 0
   using the symmetric product prod_{k=0}^{2m-1} (z - (alpha + zeta^k beta)^{2m})!
   Check:
   (a) What is the EXACT asymptotic rate function M(e2, e3) = lim_{n->inf} (1/n^2) log |P_n(1, e2, e3)|?!
       WAIT A SECOND!!
       When e2, e3 are FIXED (O(1)), x, y, z are O(1), so |x|^{1/n}, |y|^{1/n}, |z|^{1/n} -> 1!!
       WAIT: THAT MEANS M(e2, e3) IS CONSTANT = Cl_2(pi/3)/pi FOR ALL O(1) (e2, e3)!!
       TO SEE THE NONTRIVIAL RATE FUNCTION IN (e2, e3), HOW DO (e2, e3) SCALE WITH n?!
       CHECK THE SADDLE POINT SCALE:
       Since log|c_{a,b,c}| ~ n^2 * f(a/n, b/n, c/n),
       d/db log|c_{a,b,c}| ~ n * (df/dbeta)  ==>  R2 = |e2| ~ exp(-n * lambda_2)!!
       d/dc log|c_{a,b,c}| ~ n * (df/dgamma) ==>  R3 = |e3| ~ exp(-n * lambda_3)!!
   OH!! MY!! GOD!!
   IF e2 ~ exp(-n * lam2) AND e3 ~ exp(-n * lam3) WITH e1 = 1:
   THEN THE THREE ROOTS OF T^3 - T^2 + e2*T - e3 = 0 ARE:
     x ~ 1          ==>  |x|^{1/n} ~ 1 !!
     y ~ e2         ==>  |y|^{1/n} ~ exp(-lam2) =: r_2 !!
     z ~ e3 / e2    ==>  |z|^{1/n} ~ exp(-(lam3 - lam2)) =: r_3 !!
   HOLD EVERYTHING!! STOP!! FREEZE!!
   THAT MEANS AT THE SADDLE POINT (e1 = 1, e2 = +- r_2^n, e3 = +- (r_2 r_3)^n):
   THE THREE ROOTS (x, y, z) HAVE EXACT n-th ROOT MODULI (1, r_2, r_3) WITH 1 > r_2 > r_3 > 0!!
   CHECK THIS RIGHT NOW IN PYTHON!!
"""
import mpmath as mp
import numpy as np

mp.mp.dps = 50
cl2 = mp.clsin(2, mp.pi/3)
I_exact = mp.quad(lambda th: mp.log(2*mp.cos(th/2)), [-2*mp.pi/3, 2*mp.pi/3]) / (2*mp.pi)
print(f"Cl_2(pi/3)/pi = {float(cl2/mp.pi):.10f}, exact quad = {float(I_exact):.10f}")

# Check p_n(1, omega, omega^2)
omega = mp.exp(2*mp.pi*1j/3)
for n in [6, 12, 18, 24, 36, 60]:
    beta = mp.exp(2*mp.pi*1j / (3*n))
    prod = mp.mpc(1)
    for k in range(n):
        prod *= (omega**2 - (1 + mp.exp(2*mp.pi*1j*k/n) * beta)**n)
    val = prod.real
    print(f"n={n:3d}: sign={1 if val>0 else -1:+d}, log|p_n(1,w,w^2)|/n^2 = {float(mp.log(abs(val))/(n*n)):.10f}")

print("\n=== TEST THE EXPONENTIAL SADDLE-POINT SCALING: e1=1, e2 = s2*r2^n, e3 = s3*(r2*r3)^n ===")
# Compute exact P_{2m}(1, e2, e3) via cubic roots and product for n = 20, 40, 60!
def eval_Pn_from_e(n, e2, e3):
    # Roots of T^3 - T^2 + e2*T - e3 = 0 using mpmath polyroots
    roots = mp.polyroots([mp.mpc(1), mp.mpc(-1), mp.mpc(e2), mp.mpc(-e3)], maxsteps=200, extraprec=40)
    x, y, z = roots
    # alpha = x^{1/n}, beta = y^{1/n}
    alpha = x**(mp.mpf(1)/n)
    beta = y**(mp.mpf(1)/n)
    prod = mp.mpc(1)
    for k in range(n):
        zeta_k = mp.exp(2*mp.pi*1j*k / n)
        prod *= (z - (alpha + zeta_k * beta)**n)
    return prod

# Test Mahler measure of (1 + r2*U + r3*V) vs (1/n^2) log |P_n(1, e2, e3)| for r2 = 0.7, r3 = 0.4!
r2 = mp.mpf('0.7')
r3 = mp.mpf('0.4')
for n in [12, 20, 30, 40]:
    e2_pos = r2**n
    e3_pos = (r2*r3)**n
    # Evaluate at all 4 sign combos (e2, e3) = (+-e2_pos, +-e3_pos)!
    v_pp = eval_Pn_from_e(n,  e2_pos,  e3_pos).real
    v_pm = eval_Pn_from_e(n,  e2_pos, -e3_pos).real
    v_mp = eval_Pn_from_e(n, -e2_pos,  e3_pos).real
    v_mm = eval_Pn_from_e(n, -e2_pos, -e3_pos).real
    even_m = (v_mp + v_mm) / 2  # P_even(1, -R2, R3) -> should be > 0!
    odd_p  = -(v_pp - v_pm) / 2 # -P_odd(1, +R2, R3) -> should be > 0!
    print(f"n={n:2d} (r2=0.7, r3=0.4):")
    print(f"  P(1,+R2,+R3)={mp.nstr(v_pp,6)}, P(1,+R2,-R3)={mp.nstr(v_pm,6)} -> -P_odd(+R2,+R3)={mp.nstr(odd_p,6)} (log/n^2={float(mp.log(abs(odd_p))/(n*n)):.6f})")
    print(f"  P(1,-R2,+R3)={mp.nstr(v_mp,6)}, P(1,-R2,-R3)={mp.nstr(v_mm,6)} ->  P_even(-R2,+R3)={mp.nstr(even_m,6)} (log/n^2={float(mp.log(abs(even_m))/(n*n)):.6f})")

# NOW CHECK ON THE TORUS (e2 = R2 * e^{i*phi2}, e3 = R3 * e^{i*phi3}) at n = 20!
print("\n=== TORUS MODULUS TEST at n=20, r2=0.7, r3=0.4 ===")
n = 20
R2, R3 = r2**n, (r2*r3)**n
v_ev_ref = abs((eval_Pn_from_e(n, -R2, R3) + eval_Pn_from_e(n, -R2, -R3))/2)
v_od_ref = abs((eval_Pn_from_e(n,  R2, R3) - eval_Pn_from_e(n,  R2, -R3))/2)
max_ev_ratio = 0.0
max_od_ratio = 0.0
for j2 in range(12):
    phi2 = 2*mp.pi*j2/12
    for j3 in range(6):
        phi3 = mp.pi*j3/6
        e2c = R2 * mp.exp(1j*phi2)
        e3c = R3 * mp.exp(1j*phi3)
        vp = eval_Pn_from_e(n, e2c,  e3c)
        vm = eval_Pn_from_e(n, e2c, -e3c)
        ev_r = float(abs((vp + vm)/2) / v_ev_ref)
        od_r = float(abs((vp - vm)/2) / v_od_ref)
        if not (j2 == 6 and j3 == 0):
            if ev_r > max_ev_ratio: max_ev_ratio = ev_r
        if not (j2 == 0 and j3 == 0):
            if od_r > max_od_ratio: max_od_ratio = od_r
print(f"Max torus ratio away from peak: |P_even|/P_even(-R2,R3) = {max_ev_ratio:.6f}, |P_odd|/(-P_odd(R2,R3)) = {max_od_ratio:.6f}")
