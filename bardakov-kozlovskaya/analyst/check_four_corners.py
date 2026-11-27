"""
Check the EXACT phase/sign structure of P_{2m}(1, e2, e3) when (e2, e3) = (+-r2^n, +-r2^n*r3^n)!
WAIT A SECOND!! STOP EVERYTHING AND DO THE ALGEBRA OF THE ROOTS OF
  T^3 - T^2 + e2*T - e3 = 0
AT THE FOUR SIGN CORNERS:
  (1) e2 = +R2, e3 = +R3  (R2 = r2^n, R3 = r2^n r3^n, 1 > r2 > r3 > 0):
      roots are x ~ +1, y ~ +r2^n, z ~ +r3^n  (ALL THREE POSITIVE REAL!!)
  (2) e2 = +R2, e3 = -R3:
      roots are x ~ +1, y ~ +r2^n, z ~ -r3^n  (TWO POSITIVE, ONE NEGATIVE!!)
  (3) e2 = -R2, e3 = +R3:
      roots are x ~ +1, y ~ -r2^n, z ~ -r3^n  (ONE POSITIVE, TWO NEGATIVE!!)
  (4) e2 = -R2, e3 = -R3:
      roots are x ~ +1, y ~ -r2^n, z ~ +r3^n  (TWO POSITIVE, ONE NEGATIVE!!)

OH!! MY!! GOD!! STOP EVERYTHING!!
LOOK AT p_{2m}(x, y, z) = prod_{k=0}^{2m-1} (z - (alpha + zeta^k beta)^{2m}) AT ALL 4 CORNERS!!
Wait: when x > 0, y > 0, z > 0 (Corner 1):
  alpha = x^{1/2m} > 0, beta = y^{1/2m} > 0!
  Look at the 2m factors (z - (alpha + zeta^k beta)^{2m}):
  - k = 0: alpha + beta > 0 -> (alpha + beta)^{2m} > z -> NEGATIVE factor!
  - k = m: alpha - beta > 0 (since r2 < 1!) -> compare (alpha - beta)^{2m} with z = r3^{2m}:
    IF r3 < 1 - r2 (OUTSIDE triangle inequality: 1 > r2 + r3!!), then (alpha - beta)^{2m} > z -> NEGATIVE factor!
    IF r3 > 1 - r2 (INSIDE triangle inequality: 1 < r2 + r3!!), then (alpha - beta)^{2m} < z -> POSITIVE factor!!
  - ALL OTHER 2m - 2 factors come in COMPLEX CONJUGATE PAIRS (k and 2m-k), whose product is STRICTLY POSITIVE!!
HOLD EVERYTHING!! STOP!! FREEZE!!
WHAT HAPPENS AT CORNER 2 (x > 0, y > 0, z < 0)?!
  - k = 0: z - (alpha + beta)^{2m} < 0 (NEGATIVE!)
  - k = m: z - (alpha - beta)^{2m} < 0 (NEGATIVE! ALWAYS NEGATIVE because z < 0 and (alpha-beta)^{2m} > 0!!)
  - ALL OTHER 2m-2 factors come in conjugate pairs -> STRICTLY POSITIVE!!
  SO AT CORNER 2 (x>0, y>0, z<0):
    p_{2m}(x, y, z) = (-)*(-)*(+) = STRICTLY POSITIVE ALWAYS!!
  NOW CHECK CORNER 1 vs CORNER 2 WHEN 1 < r2 + r3 (INSIDE THE TRIANGLE!!):
    - Corner 1 (e2=+R2, e3=+R3): (-)*(+)*(+) = STRICTLY NEGATIVE!!
    - Corner 2 (e2=+R2, e3=-R3): (-)*(-)*(+) = STRICTLY POSITIVE!!
    SO P_{2m}(1, +R2, +R3) - P_{2m}(1, +R2, -R3) = (NEG) - (POS) = STRICTLY NEGATIVE (NO CANCELLATION AT ALL!!)!!
    WAIT A SECOND!! WHAT ABOUT WHEN 1 > r2 + r3 (OUTSIDE THE TRIANGLE!!)?!
    Let's check BOTH regimes (inside and outside triangle 1 vs r2+r3) AND ALSO Corners 3 & 4 RIGHT NOW!
"""
import mpmath as mp

mp.mp.dps = 50
def eval_Pn(n, x, y, z):
    alpha = mp.mpc(x)**(mp.mpf(1)/n)
    beta  = mp.mpc(y)**(mp.mpf(1)/n)
    prod = mp.mpc(1)
    for k in range(n):
        prod *= (z - (alpha + mp.exp(2*mp.pi*1j*k/n) * beta)**n)
    return prod

print("=== Check ALL 4 sign corners of (x, y, z) for even n=20 ===")
n = 20
for (r2, r3, desc) in [(0.7, 0.5, "INSIDE triangle (r2+r3=1.2 > 1)"),
                       (0.7, 0.4, "INSIDE triangle (r2+r3=1.1 > 1)"),
                       (0.4, 0.2, "OUTSIDE triangle (r2+r3=0.6 < 1)"),
                       (0.8, 0.1, "OUTSIDE triangle (r2+r3=0.9 < 1)")]:
    r2m, r3m = mp.mpf(r2), mp.mpf(r3)
    R2, R3 = r2m**n, (r2m*r3m)**n
    def ev_e(e2, e3):
        roots = mp.polyroots([mp.mpc(1), mp.mpc(-1), mp.mpc(e2), mp.mpc(-e3)], maxsteps=200, extraprec=40)
        return eval_Pn(n, roots[0], roots[1], roots[2]).real
    v_pp = ev_e( R2,  R3)
    v_pm = ev_e( R2, -R3)
    v_mp = ev_e(-R2,  R3)
    v_mm = ev_e(-R2, -R3)
    print(f"\n{desc}:")
    print(f"  (+,+): {mp.nstr(v_pp, 8)}   (+,-): {mp.nstr(v_pm, 8)}  --> -P_odd(+R2,+R3) = {mp.nstr(-(v_pp-v_pm)/2, 8)}")
    print(f"  (-,+): {mp.nstr(v_mp, 8)}   (-,-): {mp.nstr(v_mm, 8)}  -->  P_even(-R2,+R3)= {mp.nstr((v_mp+v_mm)/2, 8)}")
    print(f"  Ratio |P(+,+)|/|P(+,-)| = {mp.nstr(abs(v_pp/v_pm), 6)}, |P(-,+)|/|P(-,-)| = {mp.nstr(abs(v_mp/v_mm), 6)}")
