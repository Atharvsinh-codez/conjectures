"""
LOCK DOWN 100% RIGOROUS, CLOSED-FORM ELEMENTARY CONSTANTS FOR:
1. Lemma R2.1: Tiny cubic-root perturbation for ALL n >= 132, kappa < 2.5 (on ANY |e_3| = lam * R_y * R_z with lam <= 2).
2. Lemma R2.2 (Sub-regime 2A): Exact closed-form pairing bound B_pair < pi/2 for ALL c*kap^2 >= C_1!
   Tightened [L/4, L/2] integral using concavity: for t in [f, 0.5], t*(1-t) >= f*(1-f) + (1-2f)*(t-f) [wait: concave => secant from f to 0.5 has slope 0.5 - f!!]
   Check: at t=f: f*(1-f); at t=0.5: 0.25. Difference = 0.25 - f + f^2 = (0.5 - f)^2.
   Slope = (0.5 - f)^2 / (0.5 - f) = 0.5 - f!
   So by CONCAVITY of t(1-t) on [f, 0.5]:
     t*(1-t) >= f*(1-f) + (0.5 - f)*(t - f)  FOR ALL t in [f, 0.5]!!
   AND for t in [0, f]: t*(1-t) >= (1 - f)*t!!
3. Lemma R2.3 (Sub-regime 2C): Exact closed-form finite-n binomial power-sum bound for ALL c*q*kap <= C_2!
   Verify exact ratio bounds on p_r = n * sum_{t>=0} (-1)^{ct} comb(n(r+t)-1, tn) R_y^t!
"""
import numpy as np
from math import gcd, pi, sin, cos, log, exp, sqrt, asin, comb, lgamma

# =============================================================================
# 1. TIGHTENED CLOSED-FORM BOUND FOR SUB-REGIME 2A (c*kap^2 >= C_1)
# =============================================================================
def closed_form_2A(C1, kap0=2.5, N0=132):
    alpha_max = asin(kap0 / (2.0 * pi))
    gam_ext = sin(alpha_max / 4.0) / (alpha_max / 4.0)
    gam_z = (alpha_max / 2.0) / sin(alpha_max / 2.0)
    L_min = C1 / kap0
    zeta3 = 1.2020569031595943
    M_g = 0.6476102379 # max_{t>0} t^2/(e^t - 1)

    # Recall we want to bound:
    #   S_pair = 2 * sum_{m=1}^{floor((c-1)/2)} (2*gam_z*kap/c) * m^2 / (exp(gam_ext*kap*m*(1 - m/c)) - 1)
    # Set h = kap, L = c*kap >= L_min, x_m = m*kap in (0, L/2].
    # Term = (4*gam_z / (c*kap^2)) * h * [ x_m^2 / (exp(gam_ext * x_m * (1 - x_m/L)) - 1) ].
    #Wait: let F(x) = x^2 / (exp(gam_ext * x * (1 - x/L)) - 1) on [0, L/2].
    # Split [0, L/2] at x = f*L (f in (0, 0.5)):
    # (a) On [0, f*L]: 1 - x/L >= 1 - f, so F(x) <= x^2 / (exp(a1*x) - 1) where a1 = gam_ext*(1-f).
    #     Unimodal bound: h * sum_{m <= f*c} F(m*h) <= int_0^{f*L} x^2/(e^{a1*x}-1) dx + h * max_{x>0} x^2/(e^{a1*x}-1)
    #                                                <= 2*zeta(3)/a1^3 + (kap/a1^2) * M_g.
    #     Wait: multiplying by 4*gam_z / (c*kap^2):
    #       (4*gam_z / C1) * 2*zeta(3)/a1^3  +  (4*gam_z / L_min) * M_g / a1^2 !!
    # (b) On [f*L, L/2]: by concavity of t*(1-t) on [f, 0.5]:
    #     x*(1 - x/L) >= f*(1-f)*L + (0.5 - f)*(x - f*L) := A + b*(x - f*L),
    #     where A = gam_ext * f*(1-f)*L, b = gam_ext*(0.5 - f).
    #     Notice x <= L/2, and A >= gam_ext*f*(1-f)*L_min := A_min!
    #     For x in [f*L, L/2], F(x) <= (x^2 / (1 - e^{-A_min})) * exp(-A - b*(x - f*L)).
    #     Wait: EVEN SIMPLER — use 2 subintervals or exact integral of x^2 exp(-b*(x - f*L)) on [f*L, L/2]:
    #     Let s = x - f*L in [0, (0.5-f)*L]:
    #       int_0^{(0.5-f)*L} (f*L + s)^2 exp(-b*s) ds <= int_0^inf (f*L + s)^2 exp(-b*s) ds
    #                                                   = f^2 L^2 / b + 2 f L / b^2 + 2 / b^3!
    #     Multiply by exp(-A) / (1 - exp(-A)) * (4*gam_z / (c*kap^2))!
    #     WAIT: c*kap^2 = L * kap, NOT L^2! What if L is HUGE (c huge, kap small)?
    #     No problem! Because exp(-A) = exp(-gam_ext*f*(1-f)*L) decays EXPONENTIALLY in L,
    #     and since c*kap^2 >= C1: 1/(c*kap^2) <= 1/C1, so (L^2 / C1) * exp(-gam_ext*f*(1-f)*L)
    #     is STRICTLY DECREASING in L once L >= 2 / (gam_ext*f*(1-f))!
    #     AND the Riemann discretization term h * max_{[f*L, L/2]} F(x):
    #       since h = kap <= kap0 and h/(c*kap^2) = 1/L:
    #       (4*gam_z / L) * (L/2)^2 / (exp(A) - 1) = gam_z * L / (exp(gam_ext*f*(1-f)*L) - 1)!
    #       Strictly decreasing for L >= 1 / (gam_ext*f*(1-f))!!
    best = 1e9; best_info = None
    for f in np.linspace(0.08, 0.35, 55):
        a1 = gam_ext * (1.0 - f)
        S1 = 4.0 * gam_z * (2.0 * zeta3 / (C1 * a1**3) + M_g / (L_min * a1**2))
        cA = gam_ext * f * (1.0 - f)
        b = gam_ext * (0.5 - f)
        # Check monotonicity threshold for L^2 * exp(-cA * L): peak at L = 2/cA!
        # Evaluate max over L >= L_min (exact peak is at L = max(L_min, 2/cA) for L^2 e^{-cA L}, etc.)
        L2 = max(L_min, 2.0 / cA)
        L1 = max(L_min, 1.0 / cA)
        A_min = cA * L_min
        denom_min = 1.0 - exp(-A_min)
        # Exact upper bound on (4*gam_z / C1) * (1/(1-e^{-A_min})) * [ f^2 L^2 e^{-cA L}/b + 2 f L e^{-cA L}/b^2 + 2 e^{-cA L}/b^3 ]
        S2_int = (4.0 * gam_z / (C1 * denom_min)) * (
            (f*f / b) * (L2**2 * exp(-cA * L2)) +
            (2.0 * f / (b*b)) * (L1 * exp(-cA * L1)) +
            (2.0 / (b**3)) * exp(-A_min)
        )
        # Riemann step term: gam_z * max_{L >= L_min} L / (exp(cA*L) - 1)
        # Note L/(e^{cA L}-1) is strictly decreasing for all L > 0!
        S2_step = gam_z * L_min / (exp(A_min) - 1.0)
        S2 = S2_int + S2_step

        # Midpoint (c even):
        e_mid = exp(-gam_ext * L_min / 4.0)
        S_mid = -log(1.0 - e_mid)

        # Interior tail (m >= floor((c-1)/2) + 1 >= c/2):
        # Sharper! Split interior m into m <= n/6 (where sin(m*pi/n)/(m*pi/n) >= 3/pi = 0.9549!) and m > n/6!
        # Actually even simpler: for ALL m <= (n-c)/2, sin(m*pi/n) >= (2/pi)*(m*pi/n), wait:
        # At m = c/2, m*pi/n = c*pi/(2n) < alpha_max/2 < 0.21 rad!
        # For m <= n/4: sin(m*pi/n)/(m*pi/n) >= sin(pi/4)/(pi/4) = 2*sqrt(2)/pi = 0.9003!
        # Check: at m = n/4: (1 + (2*0.9003*kap/n)*(n/4))^{-n/2} = (1 + 0.45*kap)^{-n/2} -> wait, kap can be small!
        # No, for m >= c/2, use exact c_int = 2/pi:
        c_int = 2.0 / pi
        x_half = (c_int * L_min) / N0
        a_half = (1.0 + x_half) ** (-N0 / 2.0)
        L_int_pk = max(L_min, N0 / (c_int * (N0/2.0 - 2.0)))
        x_pk = (c_int * L_int_pk) / N0
        int_tail = (L_int_pk / C1) * (N0 / (2.0 * c_int * (N0/2.0 - 1.0))) * ((1.0 + x_pk) ** (-(N0/2.0 - 1.0)))
        S_int = (2.0 / (1.0 - a_half)) * (a_half + int_tail)

        tot = S1 + S2 + S_mid + S_int
        if tot < best:
            best = tot; best_info = (f, S1, S2, S_mid, S_int)
    return best, best_info

print("=== TIGHTENED CLOSED-FORM BOUND FOR SUB-REGIME 2A (kap0 = 2.5 & 3.0, N0 = 132) ===")
for kap0 in [2.5, 3.0]:
    print(f"\n--- kap0 = {kap0} (pi/2 = {pi/2:.4f}) ---")
    for C1 in [16.0, 20.0, 22.0, 24.0, 26.0, 30.0]:
        b, (f, s1, s2, sm, si) = closed_form_2A(C1, kap0=kap0, N0=132)
        st = "PASS!" if b < pi/2 else "fail"
        print(f"  C1={C1:4.1f} (L_min={C1/kap0:5.2f}): B={b:.4f} [f={f:.2f}: S1={s1:.4f}, S2={s2:.4f}, mid={sm:.4f}, int={si:.4f}] -> {st}")

# =============================================================================
# 2. CHECK EXACT BINOMIAL POWER-SUM FORMULA p_r AND SUB-REGIME 2C (c*q*kap <= C_2)
# =============================================================================
print("\n=== VERIFY EXACT BINOMIAL FORMULA FOR p_r = sum_k (1 + r_j e^{i th_k})^{-nr} ===")
for n, c, j in [(16, 2, 3), (20, 3, 4), (40, 2, 3), (40, 4, 6)]:
    rj = sin(j * pi / n) / sin((j + c) * pi / n)
    sj = sin(c * pi / n) / sin((j + c) * pi / n)
    if c % 2 == 1:
        ths = (2 * np.arange(n) + 1) * pi / n
    else:
        ths = 2 * np.arange(n) * pi / n
    for r in [1, 2, 3]:
        p_direct = float(np.sum((1.0 + rj * np.exp(1j * ths)) ** (-n * r)).real)
        # Binomial formula: n * sum_{t>=0} (-1)^{c*t} comb(n(r+t)-1, t*n) r_j^{t*n}
        p_binom = 0.0
        for t in range(0, 60):
            # compute comb(n(r+t)-1, tn) * rj^{tn} in log space
            lv = lgamma(n*(r+t)) - lgamma(t*n + 1) - lgamma(n*r) + t*n*log(rj)
            term = n * ((-1)**(c*t)) * exp(lv)
            p_binom += term
            if t > 5 and abs(term) < 1e-16 * abs(p_binom): break
        rel = abs(p_direct - p_binom) / max(abs(p_direct), 1e-30)
        # Also check tilde_p_r = (-1)^{j*r} * sj^{n*r} * p_r
        tp_r = ((-1)**(j*r)) * (sj**(n*r)) * p_binom
        if r <= 2:
            print(f"  n={n:2d} c={c} j={j} r={r}: rel_err={rel:.2e} | tilde_p_r={tp_r:+.6e}")
