"""
1. Verify 3-piece split for Sub-regime 2A in 50-digit mpmath.iv -> get C1 down to ~22!
2. Check exact finite-n binomial c = 3 proof (6*e_3 = p_1^3 - 3*p_1*p_2 + 2*p_3) for ALL kap < 2.5!
3. Test the EXACT Major-Arc ($D = 0$ and $D != 0$) closed-form Bernoulli-B2 + Polya-Szego lower bound
   across EVERY surviving corner of G = {c >= 3, c*kap^2 < 22, c*q*kap > 0.14, X >= 1.45}:
   - D = 0: q=1 (c=4, 10, 50, 500), q=2 (c=20, 200), q=3 (c=30, 300), q=4 (c=200)
   - D != 0: D=1 (c=50, 200, 5000; tau=0.4..10), D=2 (c=200, 500), D=3 (c=200, 500), r0=3 D=1 (c=200, j=267),
             r0=2 D=-2 (c=200, j=299), r0=13 D=-1 (c=200, j=323, kap=0.001)!
"""
import numpy as np
import mpmath as mp
from mpmath import iv
from math import pi, sin, cos, log, exp, sqrt, lgamma

iv.dps = 50; mp.mp.dps = 50
PI = iv.pi; N0 = iv.mpf(132); KAP0 = iv.mpf('2.5')

# =============================================================================
# 1. 3-PIECE SPLIT FOR SUB-REGIME 2A IN 50-DIGIT mpmath.iv
# =============================================================================
alpha_max = iv.mpf('0.409284')
gam_ext = iv.sin(alpha_max / 4) / (alpha_max / 4)
gam_z = (alpha_max / 2) / iv.sin(alpha_max / 2)
zeta3_iv = iv.mpf(mp.zeta(3)) + iv.mpf([-1e-40, 1e-40])
t0_iv = iv.mpf(['1.59362426', '1.59362427'])
M_g_iv = (t0_iv ** 2) / (iv.exp(t0_iv) - 1)

def check_2A_3piece(C1_val):
    C1 = iv.mpf(C1_val)
    L_min = C1 / KAP0
    # Split [0, 0.5] into [0, f1] (pairing), [f1, f2] (pairing with lower exp), [f2, 0.5] (individual log)!
    # Even cleaner: [0, f1] (pairing, a1 = gam_ext*(1-f1)), [f1, 0.5] (individual log, OR 2-piece individual log [f1, f2] + [f2, 0.5]!)
    # Wait: why was S2 = 0.44 at C1=42? Because on [f1, 0.5], concavity slope b = gam_ext*(0.5 - f1) is smaller than at f1!
    # Split [f1, 0.5] at f2 = 0.35:
    #   On [f1, f2]: t(1-t) >= f1*(1-f1) + (1 - f1 - f2)*(t - f1)  (MUCH steeper slope 1 - f1 - f2!)
    #   On [f2, 0.5]: t(1-t) >= f2*(1-f2) + (0.5 - f2)*(t - f2)    (MUCH larger base f2*(1-f2)!)
    best = 1e9; best_p = None
    for f1_f in [0.14, 0.16, 0.18, 0.20, 0.22]:
        for f2_f in [0.30, 0.34, 0.38]:
            if f2_f <= f1_f + 0.05: continue
            f1 = iv.mpf(str(f1_f)); f2 = iv.mpf(str(f2_f))
            # (a) [0, f1]: pairing bound
            a1 = gam_ext * (1 - f1)
            S1 = 4 * gam_z * (2 * zeta3_iv / (C1 * (a1**3)) + M_g_iv / (L_min * (a1**2)))
            # (b) [f1, f2]: individual bound!
            # Note: for m in (f1*c, f2*c], 2*(log(1/(1-|a_{-m}^{-1}|)) + log(1/(1-|a_m|)))
            # Since |a_m| <= exp(-gam_ext * kap * m) <= exp(-gam_ext * f1 * L) AND |a_{-m}^{-1}| <= exp(-gam_ext * f1*(1-f1)*L):
            # Better: keep |a_{-m}^{-1}| and |a_m| SEPARATE! (|a_m| is MUCH smaller than |a_{-m}^{-1}|!)
            cA1 = gam_ext * f1 * (1 - f1)
            b1 = gam_ext * (1 - f1 - f2)
            A1_min = cA1 * L_min
            L1_pk = L_min if float((cA1 * L_min).a) >= 1.0 else (1 / cA1)
            # Exterior [f1, f2] (factor 2 for conjugate symmetry):
            S2a_ext = (2 / (1 - iv.exp(-A1_min))) * (iv.exp(-A1_min) + (L1_pk * iv.exp(-cA1 * L1_pk)) / (C1 * b1))
            # (c) [f2, 0.5]: exterior individual bound!
            cA2 = gam_ext * f2 * (1 - f2)
            b2 = gam_ext * (iv.mpf('0.5') - f2)
            A2_min = cA2 * L_min
            L2_pk = L_min if float((cA2 * L_min).a) >= 1.0 else (1 / cA2)
            S2b_ext = (2 / (1 - iv.exp(-A2_min))) * (iv.exp(-A2_min) + (L2_pk * iv.exp(-cA2 * L2_pk)) / (C1 * b2))
            # (d) ALL unpaired interior m > f1*c (including m in (f1*c, c/2] AND m > c/2!!):
            # Using |a_m| <= (1 + (4/(pi*N0)) * kap * m)^{-N0/2}:
            # Wait: for m <= c/2, m*pi/n < alpha_max/2 => sin(m*pi/n)/(m*pi/n) >= sinc(alpha_max/2) > 0.993!
            # Even using universal c_int = 2/pi:
            c_int = 2 / PI
            x_f1 = (2 * c_int * f1 * L_min) / N0
            a_f1 = (1 + x_f1) ** (-N0 / 2)
            L_ipk = N0 / (2 * c_int * f1 * (N0 / 2 - 2))
            L_use = L_min if float(L_min.a) >= float(L_ipk.b) else L_ipk
            x_use = (2 * c_int * f1 * L_use) / N0
            int_tail = (L_use / C1) * (N0 / (2 * c_int * (N0 / 2 - 1))) * ((1 + x_use) ** (-(N0 / 2 - 1)))
            S_int = (2 / (1 - a_f1)) * (a_f1 + int_tail)
            # (e) Midpoint (c even):
            S_mid = -iv.log(1 - iv.exp(-gam_ext * L_min / 4))
            tot = S1 + S2a_ext + S2b_ext + S_int + S_mid
            if float(tot.b) < best:
                best = float(tot.b); best_p = (f1_f, f2_f, float(S1.b), float(S2a_ext.b), float(S2b_ext.b), float(S_int.b), float(S_mid.b))
    return best, best_p

print("=== 1. 3-PIECE SUB-REGIME 2A CERTIFICATION IN 50-DIGIT mpmath.iv ===")
for C1_val in [24, 26, 28, 30, 32, 36]:
    b, p = check_2A_3piece(C1_val)
    print(f"  C1={C1_val:2d} (L_min={C1_val/2.5:5.2f}): B_2A <= {b:.5f} {'PASS!' if b < pi/2 else 'fail'} {p}")

# =============================================================================
# 2. CHECK EXACT FINITE-n BINOMIAL PROOF FOR c = 3 (ALL kap < 2.5)!
# =============================================================================
print("\n=== 2. EXACT FINITE-n c = 3 CHECK: (-1)^j * (p_1^3 - 3*p_1*p_2 + 2*p_3) FOR ALL kap < 2.5 ===")
for n in [132, 320, 1000, 10000]:
    c = 3
    min_rat_div3 = 1e9; min_rat_nondiv3 = 1e9
    for j in range(4, (n - 3) // 2):
        kap = 2*pi*sin(j*pi/n)*sin((j+3)*pi/n)/sin(3*pi/n)
        if kap >= 2.5: break
        rj = sin(j*pi/n)/sin((j+3)*pi/n)
        sj = sin(3*pi/n)/sin((j+3)*pi/n)
        # Note c=3 is ODD => p_r = n * sum_{t>=0} (-1)^t comb(n(r+t)-1, t*n) r_j^{tn}
        # Compute tilde_p_r = (-1)^{jr} sj^{nr} p_r!
        tp = []
        for r in [1, 2, 3]:
            t_c = int(round(r * j / 3.0))
            lo_t = max(0, t_c - 15); hi_t = t_c + 15
            terms = []
            for t in range(lo_t, hi_t + 1):
                lv = lgamma(n*(r+t)) - lgamma(t*n+1) - lgamma(n*r) + t*n*log(rj) + r*n*log(sj)
                terms.append((lv, (-1)**(t + j*r)))
            ml = max(x[0] for x in terms)
            s = sum(sg * exp(lv - ml) for lv, sg in terms)
            tp.append((ml + log(n), s))
        # Compare terms in 6*e_3(a) = tp_1^3 - 3*tp_1*tp_2 + 2*tp_3 (remember (-1)^{j(c+1)} = +1 for c=3!):
        # Wait: e_3(a_m) has sign (-1)^{j(c+1)} = +1 for c=3!
        # Check: a_m = (-1)^j sj^n lam_m => e_3(a) = (-1)^{3j} sj^{3n} e_3(lam) = (-1)^j sj^{3n} e_3(lam) > 0!
        # Let's normalize by max log of the 3 terms!
        l1 = 3*tp[0][0]; v1 = tp[0][1]**3
        l2 = tp[0][0] + tp[1][0]; v2 = -3.0 * tp[0][1] * tp[1][1]
        l3 = tp[2][0]; v3 = 2.0 * tp[2][1]
        ml = max(l1, l2, l3)
        t1 = v1 * exp(l1 - ml); t2 = v2 * exp(l2 - ml); t3 = v3 * exp(l3 - ml)
        tot = t1 + t2 + t3
        dom = tot / (abs(t1) + abs(t2) + abs(t3))
        if j % 3 == 0:
            min_rat_div3 = min(min_rat_div3, dom)
        else:
            min_rat_nondiv3 = min(min_rat_nondiv3, dom)
    print(f"  n={n:5d}: 3|j min(e3/sum|terms|)={min_rat_div3:+.6f} | 3!|j min(e3/sum|terms|)={min_rat_nondiv3:+.6f}")
