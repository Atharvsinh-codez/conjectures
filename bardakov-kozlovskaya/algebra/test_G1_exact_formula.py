"""
Discover and test the EXACT 100%-CLOSED-FORM Continued-Fraction Theorem for Piece G1!
Also check the true constant of F(gam, th) * gam * sin(th/2) on gam <= 1, and check
EXACT sum_{r>=1} |b_r^G|/r vs max_phi |Phi(phi)| on G1!
"""
import math

# 1. Check sharp constants of F(gam, th) = 4*gam*|Re z(1+z)/(1-z)^3|
print("=== 1. EXACT SHARP CONSTANTS FOR F(gam, th) ===")
for gam_max in [0.5, 1.0, 1.5, 2.0, 5.0]:
    kb_sin = 0.0; kb_th = 0.0; ka_sin = 0.0
    for ig in range(1, 200):
        gam = gam_max * ig / 200.0
        for it in range(1, 150):
            th = math.pi * it / 150.0
            z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
            F = 4.0 * gam * abs((z * (1.0 + z) / ((1.0 - z) ** 3)).real)
            s = math.sin(0.5 * th)
            kb_sin = max(kb_sin, F * gam * s)
            kb_th = max(kb_th, F * gam * (0.5 * th))
            ka_sin = max(ka_sin, F * s * s)
    print(f"  gam <= {gam_max:3.1f}: max F*sin^2 = {ka_sin:.5f}, max F*gam*sin = {kb_sin:.5f} (/pi={kb_sin/math.pi:.4f}), max F*gam*(th/2) = {kb_th:.5f} (/pi={kb_th/math.pi:.4f})")

# Wait! Check F(gam, th) vs EXACT 2-variable envelope:
# Look at the EXACT algebraic formula for F(gam, th) = 4*gam * |Re z(1+z)/(1-z)^3|!
# Let a = 2*sinh(gam/2) >= gam, b = 2*sin(th/2) in (0, 2]!
# Notice |1 - z|^2 = e^{-gam} * (4*sinh^2(gam/2) + 4*sin^2(th/2)) = e^{-gam} * (a^2 + b^2)!!
# WAIT A SECOND!! CHECK THE EXACT IDENTITY FOR Re( z(1+z)/(1-z)^3 ) IN TERMS OF (a, b)!!
# Write z = e^{-gam + i th}, multiply num & den by z^{-3/2} = e^{3gam/2 - 3i th/2}:
#   z(1+z) / (1-z)^3 = (z^{-1/2} + z^{1/2}) / (z^{-1/2} - z^{1/2})^3
#   = 2*cosh((gam - i th)/2) / (2*sinh((gam - i th)/2))^3 !!
#   = (1/4) * cosh(w) / sinh^3(w)  where  w := (gam - i th)/2 !!
# OH!! MY!! GOD!!
# Let's check THAT exact identity and its exact upper bound!
err_id = 0.0
max_ratio_exact = 0.0
for gam in [0.05, 0.2, 0.5, 1.0, 2.0, 4.0]:
    for th in [0.05, 0.2, 0.5, 1.0, 2.0, 3.0]:
        z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
        lhs = (z * (1.0 + z) / ((1.0 - z) ** 3)).real
        w = complex(0.5 * gam, -0.5 * th)
        rhs = 0.25 * (complex(math.cosh(w.real)*math.cos(w.imag), math.sinh(w.real)*math.sin(w.imag)) /
                      (complex(math.sinh(w.real)*math.cos(w.imag), math.cosh(w.real)*math.sin(w.imag)) ** 3)).real
        err_id = max(err_id, abs(lhs - rhs))
print(f"\n  EXACT identity Re[z(1+z)/(1-z)^3] = (1/4) Re[cosh(w)/sinh^3(w)] max err = {err_id:.2e}")

# NOW COMPUTE Re[cosh(u - iv) / sinh^3(u - iv)] (u = gam/2, v = th/2) EXPLICITLY!!
# Notice: (1/4) * cosh(w) / sinh^3(w) = -(1/8) * d/dw (1 / sinh^2(w))!
# Better: let S = sinh(u) cos(v), C = cosh(u) sin(v), so sinh(u - iv) = S - i C,
# and cosh(u - iv) = cosh(u) cos(v) - i sinh(u) sin(v)!
# Check the EXACT upper bound of F(gam, th) = gam * |Re[cosh(w)/sinh^3(w)]| vs simple rational function of (u, s = sin(v))!
# Notice |sinh(w)|^2 = sinh^2(u) + sin^2(v)!!
# And |cosh(w)| = sqrt(sinh^2(u) + cos^2(v)) <= cosh(u)!
# WAIT A SECOND! If you just take |cosh(w)| / |sinh(w)|^3, at v -> 0, u -> 0 it's 1/u^3 -> F <= 2 u / u^3 = 2/u^2 = 8/gam^2 (EXACT!)
# At u = 0.628 v (peak): |Re| is 3.032/(4*gam*v^2) while |cosh|/|sinh|^3 gives 1/(v^3*(1+0.394)^{3/2}) = 0.607/v^3 -> F <= 2*0.628*0.607/v^2 = 0.762/v^2 = 3.05/th^2 (ONLY 0.6% ABOVE 3.032!!)
# HOLD EVERYTHING!! CHECK THAT!!
# IS F(gam, th) <= gam * cosh(gam/2) / (sinh^2(gam/2) + sin^2(th/2))^{3/2}  EXACT FOR ALL gam > 0, th in (0, pi]?!
# YES!! BY TRIVIAL |Re(Z)| <= |Z| AND |cosh(u - iv)| = sqrt(cosh^2(u) - sin^2(v)) <= cosh(u)!!
max_triv_over_exact = 0.0
for gam in [0.01, 0.1, 0.3, 0.5, 1.0, 2.0]:
    u = 0.5 * gam
    for th in [0.01, 0.1, 0.3, 0.6, 1.2, 2.5, math.pi]:
        v = 0.5 * th
        z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
        F_ex = 4.0 * gam * abs((z * (1.0 + z) / ((1.0 - z) ** 3)).real)
        F_triv = gam * math.cosh(u) / ((math.sinh(u) ** 2 + math.sin(v) ** 2) ** 1.5)
        assert F_ex <= F_triv + 1e-12
# CHECK THE PEAK CONSTANT OF F_triv(gam, th) = gam * cosh(gam/2) / (sinh^2(gam/2) + sin^2(th/2))^{3/2}!!
# Notice: u * cosh(u) / (sinh^2(u) + s^2)^{3/2}:
#   Since u <= sinh(u) and cosh(u) = sqrt(1 + sinh^2(u)):
#   EVEN SHARPER: F_triv_2 = 2*u * sqrt(cosh^2(u) - s^2) / (sinh^2(u) + s^2)^{3/2}
#   Wait: max_u (2*u*cosh(u) / (sinh^2(u) + s^2)^{3/2}) * s^2:
#   At s -> 0: max_t (2*t / (t^2 + 1)^{3/2}) = 4/(3*sqrt(3)) = 0.76980 (vs exact 0.75806 — ONLY 1.5% LOOSER!!)
#   And for small u << s: 2*u / s^3 = gam / sin^3(th/2) (only 3rd power vs 4th power? Wait: at u << s, Re[cosh(w)/sinh^3(w)] has EXTRA factor u because cosh(-iv)/sinh^3(-iv) = i cos(v)/sin^3(v) IS PURELY IMAGINARY!!)
print("  VERIFIED: F(gam, th) <= gam * cosh(gam/2) / (sinh^2(gam/2) + sin^2(th/2))^{3/2} ALWAYS (1-line |Re Z| <= |Z|!),")
print("  AND note at u=0, cosh(-iv)/sinh^3(-iv) = i cos(v)/sin^3(v) is PURELY IMAGINARY (giving extra u^2 factor at small u)!")

# Check the EXACT closed-form formula for Re[cosh(u - iv) / sinh^3(u - iv)]!!
# Wait: expand Re[ cosh(u - iv) * sinh^3(u + iv) ] / |sinh(u - iv)|^6!
# Even simpler: check the 1-line upper bound
#   |Re[cosh(w)/sinh^3(w)]| <= u * cosh(u) * (sinh^2(u) + 3*sin^2(v)) / (sinh^2(u) + sin^2(v))^3 ??
# Let's test the EXACT formula for Re[cosh(u-iv) sinh^3(u+iv)]!
# Using cosh(u-iv) sinh(u+iv) = (1/2)(sinh(2u) + i sin(2v))
# and sinh^2(u+iv) = sinh^2(u) cos^2(v) - cosh^2(u) sin^2(v) + i * (1/2) sinh(2u) sin(2v)
#                  = (sinh^2(u) - sin^2(v)) + (i/2) sinh(2u) sin(2v)  (EXACT!!)!
# Multiply them! The REAL PART is:
#   (1/2) sinh(2u) * (sinh^2(u) - sin^2(v)) - (1/4) sinh(2u) sin^2(2v)
#   = (1/2) sinh(2u) * [ sinh^2(u) - sin^2(v) - 2 sin^2(v) cos^2(v) ]
#   = (1/2) sinh(2u) * [ sinh^2(u) - sin^2(v) * (1 + 2 cos^2(v)) ]  !!
# OH!! MY!! GOD!! CHECK THAT EXACT IDENTITY!!
err_re = 0.0
for gam in [0.02, 0.1, 0.5, 1.2, 3.0]:
    u = 0.5 * gam
    for th in [0.02, 0.2, 0.8, 1.7, 2.8, math.pi]:
        v = 0.5 * th
        z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
        F_ex = 4.0 * gam * (z * (1.0 + z) / ((1.0 - z) ** 3)).real
        S = math.sinh(u) ** 2
        T = math.sin(v) ** 2
        F_exact_formula = gam * 0.5 * math.sinh(gam) * (S - T * (1.0 + 2.0 * math.cos(v) ** 2)) / ((S + T) ** 3)
        err_re = max(err_re, abs(F_ex - F_exact_formula))
print(f"  EXACT CLOSED-FORM REAL PART IDENTITY MAX ERR = {err_re:.2e} !!")
print("  F(gam, th) = (gam/2) * sinh(gam) * |sinh^2(gam/2) - sin^2(th/2)*(1 + 2*cos^2(th/2))| / (sinh^2(gam/2) + sin^2(th/2))^3  EXACTLY!!")
