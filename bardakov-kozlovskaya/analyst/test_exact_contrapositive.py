"""
DIAGNOSE THE EXACT CONTRAPOSITIVE OF MINOR ARCS + FIXED D=0 MULTISECTION + D!=0 MAJOR ARCS!
1. Fix the constant i*pi*(q-1)/q offset in test_unified_breakthrough.py and measure the TRUE
   variation of Im R_{non-q} and Re R_{non-q} (especially on the peak |theta| <= pi/2 and globally!)
   for D = 0 (q >= 2) when c*q*kap > 0.85 AND |Phi@1| >= 1.45!
2. Compute an explicit closed-form Minor-Arc upper bound:
     B_minor(c, j, kap) = sum_{r=1}^{R} min(b_r^{pair}, b_r^{Abel}) + tail_{R+1}
   and find EVERY surviving point in G* = {c*kap^2 < 42, c*q*kap > 0.85, kap < 2.5}
   where B_minor >= 1.45!
3. Check the EXACT parameter region of surviving D = 0 and D != 0 points!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd

# =============================================================================
# 1. FIXED D = 0 (q >= 2) MULTISECTION WHEN c*q*kap > 0.85 AND |Phi@1| >= 1.45!
# =============================================================================
def test_D0_fixed(c, j, kap):
    g = gcd(c, j); q = c // g
    alpha = kap / c
    M_cut = int(np.ceil(sqrt(45.0 / max(alpha, 1e-12)))) + 10
    ms = np.arange(-(c + M_cut), M_cut + 1, dtype=np.float64)
    re_a = -alpha * ms * (c + ms)
    im_a = -2.0 * pi * ms * j / c
    keep = re_a > -40.0
    re_a = re_a[keep]; im_a = im_a[keep]

    # Solve exact saddle for positive-root main term:
    lo, hi = -80.0, 40.0
    for _ in range(60):
        vm = 0.5 * (lo + hi)
        z = np.clip(q * re_a + vm, -70.0, 70.0)
        if float(np.sum(1.0 / (1.0 + np.exp(-z)))) > c: hi = vm
        else: lo = vm
    V = 0.5 * (lo + hi)
    lam = exp(V / q)
    # Check if lam is very close to any root |a_k|^{-1}; if so, shift slightly to midpoint!
    rv = re_a + log(lam)
    # Actually, let's check BOTH unshifted lam_* AND midpoint-shifted lam!
    if float(np.min(np.abs(rv))) < 0.005 * kap:
        # Shift log(lam) to sit midway between the two nearest root log-moduli!
        s_ra = np.sort(-re_a)
        idx = int(np.searchsorted(s_ra, log(lam)))
        if 0 < idx < len(s_ra):
            lam = exp(0.5 * (s_ra[idx - 1] + s_ra[idx]))
        else:
            lam *= exp(0.01 * kap)
        V = q * log(lam)
        rv = re_a + log(lam)

    z = np.clip(q * re_a + V, -70.0, 70.0)
    p = 1.0 / (1.0 + np.exp(-z))
    H_q = float(np.sum(p * (1.0 - p))) / q

    # Compute exact R_{non-q}(psi) = sum_{r>=1, q!|r} -(A_r e^{-ir psi} + B_r e^{ir psi})/r
    # Notice: (1/q)*sum_{t=0}^{q-1} log(1 - u e^{i(psi+2pi t/q)}) = (1/q)*log(1 - u^q e^{iq psi})
    # ONLY has ZERO constant offset! Wait: why did test_unified_breakthrough.py have 1.5708?
    # Let's check! Because ub = exp(-rv - i im_a) had (ub)^q = exp(-q rv - i q im_a) = +exp(-q rv) > 0,
    # Wait: for u = e^{-rv - i im_a} with |u| < 1, log(1 - u e^{-i psi}) has ZERO branch cuts,
    # SO WHERE DID 1.5708 COME FROM?!
    # AHA!! In test_unified_breakthrough.py, ub had len(ub) = ext (an ODD number!),
    # wait: why would sum log(1 - ub e^{-i psi}) - (1/q) sum log(1 - ub^q e^{-iq psi}) equal 1.5708*i?!
    # CHECK AT psi = pi/q: 1 - ub^q e^{-i pi} = 1 + |ub|^q > 1 (REAL!), so its imag part is 0!
    # What about sum_k log(1 - ub_k e^{-i pi/q})?!
    # Wait: the roots ub_k are symmetric in CONJUGATE PAIRS EXCEPT m = -c/2 (when c is even)!
    # For m = -c/2 (when c is even, e.g. c=30, j=45, q=2):
    #   im_a[-c/2] = -2*pi*(-c/2)*j/c = +pi*j! (for j=45 odd: im_a = pi, so ub = -|ub| < 0!)
    #   At psi = pi/2: 1 - ub * e^{-i pi/2} = 1 - (-|ub|)(-i) = 1 - i|ub| -> angle = -atan(|ub|) ~ -pi/4!
    #   No, wait: why was it EXACTLY pi/2 = 1.5708?!
    # Let's print R_vals at psi = pi/q and psi = 0 directly to see!
    psis = np.linspace(-pi, pi, 512, endpoint=False)
    big = rv > 0.0; sm = ~big
    ub = np.exp(-rv[big] - 1j * im_a[big]) if np.any(big) else np.array([], dtype=np.complex128)
    us = np.exp(rv[sm] + 1j * im_a[sm]) if np.any(sm) else np.array([], dtype=np.complex128)
    R_vals = np.zeros(512, dtype=np.complex128)
    eip = np.exp(1j * psis); epm = np.conj(eip)
    eiqp = np.exp(1j * q * psis); eiqm = np.conj(eiqp)
    if len(ub) > 0:
        for k0 in range(0, len(ub), 200):
            ch = ub[k0:k0+200]
            R_vals += np.sum(np.log(1.0 - ch[None, :] * epm[:, None]), axis=1)
            R_vals -= (1.0 / q) * np.sum(np.log(1.0 - (ch**q)[None, :] * eiqm[:, None]), axis=1)
    if len(us) > 0:
        for k0 in range(0, len(us), 200):
            ch = us[k0:k0+200]
            R_vals += np.sum(np.log(1.0 - ch[None, :] * eip[:, None]), axis=1)
            R_vals -= (1.0 / q) * np.sum(np.log(1.0 - (ch**q)[None, :] * eiqp[:, None]), axis=1)

    # Check mean of R_vals (must be 0!) and value at psi = pi/q:
    idx_pk = int(np.argmin(np.abs(psis - pi / q)))
    th_pk = np.abs(np.angle(np.exp(1j * (q * psis - pi))))
    pk = th_pk <= pi / 2
    # Also check R_vals - R_vals[idx_pk] (pure oscillation around peak!):
    dIm_pk = float(np.max(np.abs(R_vals[pk].imag - R_vals[idx_pk].imag)))
    dIm_all = float(np.max(np.abs(R_vals.imag - R_vals[idx_pk].imag)))
    return q, g, lam, H_q, int(np.sum(big)), R_vals[idx_pk], float(np.max(np.abs(R_vals.real))), dIm_pk, dIm_all

print("=== 1. INSPECT R_{non-q}(psi = pi/q) AND VARIATION FOR D = 0 (q >= 2) ===")
for c, j, cqk in [
    (12, 18, 1.5), (12, 18, 3.0), (15, 20, 1.5), (20, 30, 1.5), (20, 30, 3.0),
    (20, 25, 1.5), (20, 25, 3.0), (30, 45, 1.5), (30, 45, 3.0), (30, 45, 6.0),
    (40, 50, 1.5), (40, 50, 3.0), (50, 75, 2.0), (50, 75, 5.0),
]:
    q = c // gcd(c, j); kap = cqk / (c * q)
    _, g, lam, Hq, next_, Rpk, mre, dim_pk, dim_all = test_D0_fixed(c, j, kap)
    print(f"  c={c:2d} j={j:2d} (q={q},g={g:2d}) cqk={cqk:.1f}: lam={lam:.4f} ext={next_:2d} Hq={Hq:5.2f} | R(pi/q)={Rpk.real:+.4f}{Rpk.imag:+.4f}i | max|ReR|={mre:.4f} | dIm_pk={dim_pk:.4f} dIm_all={dim_all:.4f}")

# =============================================================================
# 2. CHECK EXACT CONTRAPOSITIVE OF MINOR-ARC BOUND B_minor = sum min(b_r^pair, b_r^Abel) >= 1.45
# =============================================================================
print("\n=== 2. MAP SURVIVING POINTS OF CLOSED-FORM MINOR-ARC BOUND IN G* (c*q*kap > 0.85, c*kap^2 < 42) ===")
# Exact simple model bounds for each harmonic r >= 1 on |w| = 1:
#   b_r^{pair}: Recall from Lemma R2.2:
#     For r*c*kap >= 4: |A_r - B_r| / r <= (4 / (r^3 * c * kap^2)) * 1.15  [check exact constant!]
#     For all r >= 1:  |A_r - B_r| / r <= (1 / r) * min(4.6 / (r^2 * c * kap^2), 2.5 * sqrt(c / (r * kap)))
#   b_r^{Abel}:
#     For sin(pi*r*nu) != 0: |A_r - B_r| / r <= (1.35 / sin^2(pi*r*nu)) * min(1/(r*c), sqrt(kap/(r*c)))
# Wait: FIRST verify b_r^{pair} <= (1/r)*min(4.6/(r^2*c*kap^2), 2.5*sqrt(c/(r*kap))) on exact model!
def exact_br_model(c, jm, kap, r):
    alpha = kap / c
    M_cut = int(np.ceil(sqrt(35.0 / max(r * alpha, 1e-8)))) + 5
    ms_e = np.arange(-(c - 1), 0, dtype=np.float64)
    A_r = float(np.sum(np.exp(r * alpha * ms_e * (c + ms_e)) * np.cos(2*pi*r*ms_e*jm/c)))
    ms_i = np.concatenate([np.arange(1, M_cut + 1, dtype=np.float64),
                           np.arange(-(c + M_cut), -c, dtype=np.float64)])
    lw = -r * alpha * ms_i * (c + ms_i)
    keep = lw > -35.0
    B_r = float(np.sum(np.exp(lw[keep]) * np.cos(2*pi*r*ms_i[keep]*jm/c)))
    return abs(B_r - A_r) / r

def B_minor_bound(c, jm, kap, K_abel=1.35, C_p1=4.8, C_p2=2.2):
    nu = jm / c
    tot = 0.0
    # Sum r = 1..60 + tail r >= 61 using C_p1 / (r^3 * c * kap^2) -> tail <= C_p1 / (2 * 60^2 * c * kap^2)
    for r in range(1, 61):
        bp = min(C_p1 / ((r**3) * c * kap * kap), C_p2 * sqrt(c / (r * kap)) / r)
        s2 = sin(pi * r * nu) ** 2
        if s2 > 1e-12:
            ba = (K_abel / s2) * min(1.0 / (r * c), sqrt(kap / (r * c)))
            tot += min(bp, ba)
        else:
            tot += bp
    tot += C_p1 / (2.0 * 3600.0 * c * kap * kap)
    return tot

# Test across c in [2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50, 100, 200]
kaps = np.geomspace(0.002, 2.49, 60)
for c in [2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50, 100]:
    surv_D0 = []
    surv_Dne0 = []
    for jm in range(0, c // 2 + 1):
        g = gcd(c, jm) if jm > 0 else c
        q = c // g
        for kap in kaps:
            if c * q * kap <= 0.85 or c * kap * kap >= 42.0: continue
            bm = B_minor_bound(c, jm, kap)
            if bm >= 1.45:
                # Find best small rational p0/r0 minimizing |r0*jm - c*p0| with r0 <= sqrt(c):
                best_r0 = q; best_D = 0; best_score = 1e9
                for r0 in range(1, min(q, int(sqrt(c)) + 2) + 1):
                    p0 = round(r0 * jm / c)
                    D = abs(r0 * jm - c * p0)
                    if D == 0:
                        best_r0 = r0; best_D = 0; break
                    # Score: D / (c / r0)^{1/2}
                    sc = D * sqrt(r0) / sqrt(c)
                    if sc < best_score:
                        best_score = sc; best_r0 = r0; best_D = D
                if best_D == 0:
                    surv_D0.append((jm, q, g, kap, bm))
                else:
                    surv_Dne0.append((jm, q, best_r0, best_D, c / best_r0, kap, c*best_r0*kap, bm))
    d0_qs = sorted(set((x[1], x[2]) for x in surv_D0))
    dne0_min_ct = min((x[4] / (x[3]**2) for x in surv_Dne0), default=999.0)
    dne0_sample = sorted(set((x[2], x[3], round(x[4], 1)) for x in surv_Dne0))[:6]
    print(f"  c={c:3d}: D=0 surv={len(surv_D0):3d} (q,g)={d0_qs[:6]} | D!=0 surv={len(surv_Dne0):3d} min(ct/D^2)={dne0_min_ct:6.1f} (r0,D,ct)={dne0_sample}")
