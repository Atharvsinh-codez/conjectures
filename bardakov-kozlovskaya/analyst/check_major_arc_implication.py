"""
1. Verify 3-line EXACT FINITE-n BINOMIAL PROOF of c = 2 for ALL kappa < 2.5 (n >= 132)!
   Recall for c = 2 (even!):
     p_r = n * sum_{t>=0} comb(n(r+t)-1, t*n) R_y^t > 0  (EVERY TERM STRICTLY POSITIVE!)
     2 * e_2(lambda_0, ..., lambda_{n-1}) = p_1^2 - p_2!
   Check exact ratio p_1^2 / p_2 across ALL (n, j) with c=2, kap < 2.5 (n = 132, 320, 1000, 10000)!
2. Check the EXACT small-c / D=0 / D!=0 boundary when we use 2-step paired Abel vs 1-step paired Abel
   for r with sin(pi*r*nu) != 0!
   Specifically: check the EXACT closed-form 2-step paired Abel constant!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, lgamma

# 1. Exact c = 2 check:
print("=== 1. EXACT FINITE-n c = 2 CHECK: p_1^2 / p_2 FOR ALL kap < 2.5 ===")
for n in [132, 320, 1000, 10000]:
    c = 2
    max_odd_rat = 0.0; worst_odd = None
    min_even_rat = 1e9; worst_even = None
    for j in range(3, n // 2):
        kap = 2*pi*sin(j*pi/n)*sin((j+2)*pi/n)/sin(2*pi/n)
        if kap >= 2.5: break
        rj = sin(j*pi/n)/sin((j+2)*pi/n)
        # Compute p_1 and p_2 in log-shifted space around peak t_* = r*j/2!
        vals = []
        for r in [1, 2]:
            t_c = int(round(r * j / 2.0))
            lo_t = max(0, t_c - 15)
            hi_t = t_c + 15
            lterms = [lgamma(n*(r+t)) - lgamma(t*n+1) - lgamma(n*r) + t*n*log(rj) for t in range(lo_t, hi_t+1)]
            m_lt = max(lterms)
            vals.append(log(n) + m_lt + log(sum(exp(x - m_lt) for x in lterms)))
        # log(p_1^2 / p_2) = 2*vals[0] - vals[1]
        rat = exp(2*vals[0] - vals[1])
        if j % 2 == 1:
            if rat > max_odd_rat: max_odd_rat = rat; worst_odd = (j, kap, rat)
        else:
            if rat < min_even_rat: min_even_rat = rat; worst_even = (j, kap, rat)
    print(f"  n={n:5d}: ODD j max(p1^2/p2)={max_odd_rat:.6f} at {worst_odd} | EVEN j min(p1^2/p2)={min_even_rat:.6f} at {worst_even}")

# 2. Check exact finite-n D=0 (q=1) positive-root sector:
# For ANY (n, c, j) with c | j (q = 1) and c*kap^2 < 45, kap < 2.5, n >= 132:
# What is max_m |arg(a_m)| over ALL active roots |a_m| > 1e-6?!
print("\n=== 2. CHECK MAX PHASE DEVIATION |arg(a_m)| FOR q = 1 (c | j) WITH c*kap^2 < 45 ===")
max_ang = 0.0; worst_ang = None
for n in [132, 320, 1000, 10000]:
    for c in range(2, int(0.08 * n) + 1):
        for mult in range(2, 20):
            j = mult * c
            if j > (n - c) // 2: break
            kap = 2*pi*sin(j*pi/n)*sin((j+c)*pi/n)/sin(c*pi/n)
            if kap >= 2.5 or c * kap * kap >= 45.0: continue
            rj = sin(j*pi/n)/sin((j+c)*pi/n)
            sj = sin(c*pi/n)/sin((j+c)*pi/n)
            th_star = pi * (n - c) / n
            # Check all roots with |log|a_m|| < 12:
            # Exterior: m = 1..c-1; Interior: m = 1..min(n-c-1, 400)
            ms_e = np.arange(1, c, dtype=np.float64)
            th_e = th_star + 2.0 * ms_e * pi / n
            z_e = ((1.0 + rj * np.exp(1j * th_e)) / (sj * np.exp(1j * j * pi / n))) ** n
            ang_e = np.abs(np.angle(z_e))
            M_i = min(n - c - 1, int(sqrt(30.0 * c / max(kap, 1e-9))) + 20)
            ms_i = np.arange(1, M_i + 1, dtype=np.float64)
            th_i = th_star - 2.0 * ms_i * pi / n
            lz_i = n * np.log((sj * np.exp(1j * j * pi / n)) / (1.0 + rj * np.exp(1j * th_i)))
            keep = lz_i.real > -12.0
            ang_i = np.abs(np.angle(np.exp(lz_i[keep]))) if np.any(keep) else np.array([0.0])
            ma = max(float(np.max(ang_e)), float(np.max(ang_i)))
            if ma > max_ang:
                max_ang = ma; worst_ang = (n, c, j, kap, c*kap*kap, ma)
print(f"  GLOBAL MAX |arg(a_m)| across ALL q=1 with c*kap^2 < 45: {max_ang:.5f} rad at {worst_ang}!!")
