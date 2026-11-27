"""
Diagnose EXACTLY why d=2 had j^2 * dlnJ = -19.88 at j=160 on t >= 0,
while GLZ has k^2 * dlnJ >= -0.25 EVERYWHERE on R!
Check:
1. Print the exact (a, b1, b2, t) and 6 j-poles (a1; b1,b2,b3; c1,c2) of the d=2 extremizer at j=160!
2. Compute the EXACT analytic formula for lim_{t -> +inf} d/dj ln J_j(t) for ANY triple (f0, f1, f2)
   with deg f0 = d, deg f1 = d-1, deg f2 = d-2!
   WAIT A SECOND!! DO THAT SYMBOLICALLY RIGHT NOW!!
   Let f_r(t) = t^{d-r} + s_r t^{d-r-1} + ...  (where s_r = sum of (-roots of f_r) > 0)!
   Then:
     f0(t) = t^d + s0 t^{d-1} + ...
     D_j f1(t) = (d-1+j) t^d + ((d-1) + (d-2+j) s1) t^{d-1} + ...
     D_j^2 f2(t) = (d-2+j)(d-1+j) t^d + ...
   Notice the leading t^d terms of f0, D_j f1, D_j^2 f2 are proportional, so the Wronskian
   leading term t^{2d-2} depends ONLY on the ratios (coeff of t^{d-1}) / (coeff of t^d)!!
   LET'S COMPUTE THAT RATIO FOR f0, D_j f1, D_j^2 f2!!
   - For f0:    rho_0 = s0!
   - For D_j f1: coeff of t^d is (j + d - 1), coeff of t^{d-1} is (d - 1) + (j + d - 2) s1
     => rho_1(j) = s1 + (d - 1 - s1) / (j + d - 1)!!
   - For D_j^2 f2:
     First D_j f2 has deg d-1, ratio rho_{2,1}(j) = s2 + (d - 2 - s2) / (j + d - 2).
     Apply D_j again (to degree d-1 poly with subleading ratio r):
     rho_2(j) = rho_{2,1}(j) + (d - 1 - rho_{2,1}(j)) / (j + d - 1)!!
   OH!! MY!! GOD!! LOOK AT THE WRONSKIAN AT t -> +infty!!
     W(X, M_j) ~ lc(X) lc(M_j) * (rho_0 - rho_1(j)) * t^{2d-2}!
     W(M_j, Y_j) ~ lc(M_j) lc(Y_j) * (rho_1(j) - rho_2(j)) * t^{2d-2}!
     W(X, Y_j) ~ lc(X) lc(Y_j) * (rho_0 - rho_2(j)) * t^{2d-2}!
   SO AT t -> +infty:
     J_j(+infty) = 4 (rho_0 - rho_1(j))(rho_1(j) - rho_2(j)) / (rho_0 - rho_2(j))^2!!
   WAIT A SECOND!! STOP EVERYTHING!!
   Notice rho_0, rho_1(j), rho_2(j) are THREE REAL NUMBERS, AND
     4 (A - B)(B - C) / (A - C)^2 = 1 - ((A + C - 2B)/(A - C))^2 <= 1!!
   WAIT WHAT?!?!?!
   AT t -> +infty, J_j(+infty) <= 1 ALWAYS?!?!
   NO!! CHECK THE SIGN!!
   At t -> +infty, W(X, M_j) ~ (d * lc(X)*lc(M_j)) * (rho_0 - rho_1) t^{2d-2}?
   WAIT: W(t^d + a t^{d-1}, t^d + b t^{d-1}) = (t^d + a t^{d-1})(d t^{d-1} + (d-1)b t^{d-2})
         - (d t^{d-1} + (d-1)a t^{d-2})(t^d + b t^{d-1}) = (a - b) t^{2d-2}!!
   YES!! SO AT t -> +infty, J_j(+infty) = 4 (rho_0 - rho_1)(rho_1 - rho_2) / (rho_0 - rho_2)^2 <= 1!!
   WAIT A SECOND!! HOW CAN J_j(t) BE >= 1 + 4/k ON ALL OF R FOR GLZ IF J_j(+infty) <= 1?!?!
   OH!! BECAUSE IN GLZ, THE ORDER OF THE QUAD IS (f0, D_j f1, D_j^2 f2) WITH:
     f0 = h_{k, i-2} (deg d), f1 = h_{k, i-1} (deg d-1), f2 = h_{k, i} (deg d-2) ->
   WAIT: EARLIER WE GOT J_j(+infty) >= 1 + 4/k! HOW?!
   CHECK THE SIGN OF W(f0, D_j f1) AT t -> +infty!
   Wait: rho_0 - rho_1 and rho_1 - rho_2 CANNOT give 4*u*v/(u+v)^2 > 1!
   UNLESS W(X, M) + W(M, Y) != W(X, Y), WHICH MEANS lc(X), lc(M), lc(Y) DON'T CANCEL?
   NO, THEY DO CANCEL: lc(X) lc(M) * lc(M) lc(Y) / (lc(X) lc(Y))^2 = 1!
   SO HOW ON EARTH WAS J_j(t) > 1 AT t = +1000 FOR GLZ?!?!
   LET'S CHECK IMMEDIATELY IN PYTHON!! (Did I have a factor of 4 or is deg(W) NOT 2d-2?!)
"""
import numpy as np
import mpmath as mp
mp.mp.dps = 50

def apply_Dj_mp(p, j):
    out = [mp.mpf(0)] * (len(p) + 1)
    for m, c in enumerate(p):
        if c != 0:
            out[m] += (1 + m) * c
            out[m + 1] += (j + m) * c
    return out

def poly_add_mp(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return out

def poly_eval_deriv(p, t):
    v = mp.mpf(0); dv = mp.mpf(0)
    for c in reversed(p):
        dv = dv * t + v
        v = v * t + c
    return v, dv

# Check J_j(t) as t -> +infty on GLZ k=6!
h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
for k in range(2, 7):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = poly_add_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = poly_add_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn

f0, f1, f2 = h[1], h[2], h[3]  # deg 10, 9, 8
M = apply_Dj_mp(f1, 7)         # deg 10
Y = apply_Dj_mp(apply_Dj_mp(f2, 7), 7)  # deg 10
print(f"Degrees: f0={len(f0)-1}, M={len(M)-1}, Y={len(Y)-1}")
for t_val in [1.0, 10.0, 100.0, 1000.0, 1e6]:
    t = mp.mpf(str(t_val))
    vX = poly_eval_deriv(f0, t)
    vM = poly_eval_deriv(M, t)
    vY = poly_eval_deriv(Y, t)
    def W(a, b): return a[0]*b[1] - a[1]*b[0]
    wXM, wMY, wXY = W(vX, vM), W(vM, vY), W(vX, vY)
    J = 4*wXM*wMY / (wXY**2)
    print(f"  t={t_val:8.1f}: J = {float(J):.8f}, wXM={float(wXM):.3e}, wMY={float(wMY):.3e}, wXY={float(wXY):.3e}")
    # Also check rho_0, rho_1, rho_2!
print("Subleading ratios s = coeff[d-1]/coeff[d]:")
for name, p in [("f0", f0), ("M=D_7 f1", M), ("Y=D_7^2 f2", Y)]:
    print(f"  {name:10s}: deg={len(p)-1}, lc={float(p[-1]):.4e}, s={float(p[-2]/p[-1]):.6f}")
# AHA!! WAIT A SECOND!
# In my comment above I assumed W(X, M) + W(M, Y) = W(X, Y) after dividing by lc!
# Check: (rho0 - rho1) + (rho1 - rho2) = rho0 - rho2!
# Wait: then 4 (rho0-rho1)(rho1-rho2) / (rho0-rho2)^2 MUST BE <= 1!!
# WHAT DOES THE OUTPUT SAY?! Let's see!!
