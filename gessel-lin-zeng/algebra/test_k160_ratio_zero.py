"""
3-second check: At k=80 where (D_2+y)^80 + (D_{b*}+y) (FIXED w = 1!) BROKE plain Turan
(U_Turan/v1^2 = -0.293 at t* = -0.3 * 2^{-40}, b* = 1.023e13):
1. Does Q_{81, q}(t, y) at q = 0.5, 0.25, 0.125 PASS or FAIL?!
2. What is the ratio U_40(t*)^2 / (U_39(t*) U_41(t*))? (Remember q = 0.25 needs >= 1, q = 0.5 needs >= 2?
   NO: Q_{k,q} = sum q^{C(i,2)} U_i y^i has Turan ratio (q^{C(i,2)} U_i)^2 / (q^{C(i-1,2)} U_{i-1} q^{C(i+1,2)} U_{i+1})
   = (1/q) * U_i^2 / (U_{i-1} U_{i+1})! For real-rootedness in y, Newton requires >= 1, i.e., U_i^2/(U_{i-1} U_{i+1}) >= q!)
3. Can we push U_i^2 / (U_{i-1} U_{i+1}) BELOW 0.5 or BELOW 0.25 (or even NEGATIVE?!) at k = 80, 120, 160?!
"""
import mpmath as mp
mp.mp.dps = 120

p = [[mp.mpf(1)]]
for m in range(1, 90):
    prev = p[-1]
    nxt = [mp.mpf(0)] * (len(prev) + 1)
    for idx, c in enumerate(prev):
        nxt[idx] += (1 + idx) * c
        nxt[idx + 1] += (2 + idx) * c
    p.append(nxt)

def eval_pd(poly, tv):
    v = mp.mpf(0); dv = mp.mpf(0)
    for c in reversed(poly):
        dv = dv * tv + v
        v = v * tv + c
    return v, dv

for k in [80, 120, 160]:
    i = k // 2
    c_m1, c_0, c_1, c_2 = mp.binomial(k, i-2), mp.binomial(k, i-1), mp.binomial(k, i), mp.binomial(k, i+1)
    min_ratio = mp.mpf('1e100')
    best = None
    for alpha in [0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7]:
        t = -mp.mpf(str(alpha)) * (mp.mpf(2) ** (-(k - i)))
        vm1, dm1 = eval_pd(p[k - (i - 2)], t); vm1 *= c_m1; dm1 *= c_m1
        v0, d0 = eval_pd(p[k - (i - 1)], t); v0 *= c_0; d0 *= c_0
        v1, d1 = eval_pd(p[k - i], t); v1 *= c_1; d1 *= c_1
        v2, d2 = eval_pd(p[k - (i + 1)], t); v2 *= c_2; d2 *= c_2
        u0, u1, u2 = t*(1+t)*d0, t*(1+t)*d1, t*(1+t)*d2
        q0, q1, q2 = u0 + vm1, u1 + v0, u2 + v1
        # Choose Phi <= 1 + 2t so that U1 = q1 + Phi*v1 = 0!!
        # WAIT A SECOND!! IF U1(t) = 0, WHAT IS THE SIGN OF U0(t) * U2(t)??!
        # If U0(t) and U2(t) have the SAME SIGN at a root of U1(t), THEN U1^2 / (U0 U2) = 0 (< q FOR EVERY q > 0!!),
        # WHICH MEANS INTERLACING U_i << U_{i-1} IS ALREADY BROKEN!!
        if v1 != 0:
            Phi_zero = -q1 / v1
            if Phi_zero <= 1 + 2*t:
                U0z = q0 + Phi_zero * v0
                U2z = q2 + Phi_zero * v2
                if U0z * U2z > 0:
                    print(f"  BOOM!! k={k}, alpha={alpha}: U1=0 at Phi={float(Phi_zero):.4f} (b={(1-float(Phi_zero))/float(-t):.3e}), U0*U2 = {float(U0z*U2z):+.3e} > 0 => RATIO = 0 < q FOR ALL q>0!!")
        # Also minimize ratio U1^2 / (U0 U2) over a grid of Phi <= 1 + 2t where U0*U2 > 0!
        for phi_f in np.linspace(-10.0, 0.9, 200) if False else [-5.0, -3.0, -2.5, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.5, 0.0]:
            Phi = mp.mpf(str(phi_f))
            if Phi > 1 + 2*t: continue
            U0 = q0 + Phi * v0; U1 = q1 + Phi * v1; U2 = q2 + Phi * v2
            if U0 * U2 > 0:
                rat = (U1**2) / (U0 * U2)
                if rat < min_ratio:
                    min_ratio = rat
                    best = (alpha, phi_f, float((1 - Phi)/(-t)))
    print(f"  k={k:3d}: min U_i^2 / (U_{{i-1}} U_{{i+1}}) (with U_{{i-1}} U_{{i+1}} > 0) = {float(min_ratio):.6f} at {best}")
