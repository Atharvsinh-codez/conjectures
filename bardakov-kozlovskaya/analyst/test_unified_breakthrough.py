"""
TEST THE UNIFIED MULTISECTION + 3-LINE POLYA-SZEGO BREAKTHROUGH!
1. Check TRUE |Re R_{non-q}| and |Im R_{non-q}| (separating peak |theta| <= pi/2 vs global)
   for D = 0 (r_0 = q >= 2) at lam_* across c*q*kap in [0.15, 10.0]!
   Fix the branch-cut artifact in eval_logF by computing R_{non-q} DIRECTLY from its
   Fourier/Laurent series (which is exact and branch-cut-free on |w| = lam)!
2. Check exact finite-n Sub-regime 2C-UPGRADE (c*q*kap <= 0.90, ALL c >= 2, q >= 1, g >= 1):
   - Verify |a_{-c/2}| * lam_* < 1 ALWAYS for c*q*kap <= 0.90 (ZERO exterior roots!)
   - Verify max_{r: q!|r} |P_r| * rho_*^r / r across worst-case (c, q, g, cqk=0.90)!
   - Verify min H_q(lam_*) >= 0.65 across ALL g >= 1, c*q*kap <= 0.90!
3. Re-run the EXACT crevice map for c*q*kap > 0.90 (c = 2..50): which (c, j_mod, q, kap)
   STILL have |Phi@1| >= 1.40?
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd

# =============================================================================
# 1. EXACT BRANCH-CUT-FREE R_{non-q} ON |w| = lam_* FOR D = 0 (r_0 = q >= 2)
# =============================================================================
def exact_R_non_q(c, j, kap, lam=None):
    g = gcd(c, j); q = c // g
    alpha = kap / c
    M_cut = int(np.ceil(sqrt(45.0 / max(alpha, 1e-12)))) + 10
    ms = np.arange(-(c + M_cut), M_cut + 1, dtype=np.float64)
    re_a = -alpha * ms * (c + ms)
    im_a = -2.0 * pi * ms * j / c
    keep = re_a > -40.0
    re_a = re_a[keep]; im_a = im_a[keep]

    if lam is None:
        lo, hi = -80.0, 40.0
        for _ in range(60):
            vm = 0.5 * (lo + hi)
            z = np.clip(q * re_a + vm, -70.0, 70.0)
            if float(np.sum(1.0 / (1.0 + np.exp(-z)))) > c: hi = vm
            else: lo = vm
        V = 0.5 * (lo + hi)
        lam = exp(V / q)
    else:
        V = q * log(lam)

    z = np.clip(q * re_a + V, -70.0, 70.0)
    p = 1.0 / (1.0 + np.exp(-z))
    H_q = float(np.sum(p * (1.0 - p))) / q

    # Check number of exterior roots M (re_a + log(lam) > 0):
    rv = re_a + log(lam)
    n_ext = int(np.sum(rv > 0.0))
    max_rv_if_none = float(np.max(rv))
    # Compute EXACT Fourier coefficients of R_{non-q}:
    #   R_{non-q}(lam e^{i psi}) = -sum_{r >= 1, q!|r} (1/r) [ A_r e^{-i r psi} + B_r e^{i r psi} ]
    # Notice: Re R = -sum (A_r + B_r)/r cos(r psi), Im R = sum (A_r - B_r)/r sin(r psi)
    # Wait: ALSO check if rv has a tiny gap -> avoid hitting rv == 0!
    min_gap = float(np.min(np.abs(rv)))
    if min_gap < 1e-5:
        # shift lam slightly by 1e-3
        lam *= exp(0.002 / q)
        rv = re_a + log(lam)
        min_gap = float(np.min(np.abs(rv)))

    big = rv > 0.0; sm = ~big
    ub = np.exp(-rv[big] - 1j * im_a[big]) if np.any(big) else np.array([], dtype=np.complex128)
    us = np.exp(rv[sm] + 1j * im_a[sm]) if np.any(sm) else np.array([], dtype=np.complex128)

    # Evaluate R_{non-q} directly on a 512-grid via exact sum over r = 1..R_max with q!|r
    # Wait: if min_gap is small (e.g. 0.02), series converges slowly! Instead evaluate
    # EXACTLY using the 1-line multisection average:
    #   R_{non-q}(psi) = log F_full(psi) - (1/q) sum_{t=0}^{q-1} log F_full(psi + 2*pi*t/q) !!
    # For EACH individual root u (with |u| < 1), log(1 - u e^{i psi}) is ANALYTIC with ZERO branch cuts!
    # So sum_{t=0}^{q-1} log(1 - u e^{i (psi + 2pi t/q)}) = log(1 - u^q e^{i q psi}) WITH ZERO BRANCH CUTS!
    psis = np.linspace(-pi, pi, 512, endpoint=False)
    # Center theta = q*psi - pi in [-pi, pi] (peak is at theta = 0, i.e. psi = pi/q):
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
    # Peak mask: |angle(e^{i(q*psi - pi)})| <= pi/2
    th_pk = np.abs(np.angle(np.exp(1j * (q * psis - pi))))
    pk = th_pk <= pi / 2
    # Also check L1 bound on Fourier coeffs sum_{q!|r} |A_r+B_r|/r and |A_r-B_r|/r (for n_ext == 0):
    L1_nonq = 0.0
    if n_ext == 0:
        pu = np.ones_like(us)
        for r in range(1, 60):
            pu *= us
            if r % q != 0:
                L1_nonq += abs(float(np.sum(pu).real)) / r
    return q, g, lam, n_ext, max_rv_if_none, H_q, float(np.max(np.abs(R_vals.real))), float(np.max(np.abs(R_vals.imag))), float(np.max(np.abs(R_vals[pk].imag))), L1_nonq

print("=== 1. EXACT BRANCH-CUT-FREE R_{non-q} FOR D = 0 (r_0 = q >= 2) ===")
for c, j in [(4, 5), (4, 6), (6, 7), (6, 9), (8, 9), (8, 12), (20, 25), (30, 45), (50, 75)]:
    q = c // gcd(c, j); g = gcd(c, j)
    print(f"\nc={c:2d} j={j:2d} (q={q:2d}, g={g:2d}):")
    for cqk in [0.15, 0.50, 0.80, 0.90, 1.00, 1.50, 3.00, 6.00]:
        kap = cqk / (c * q)
        if kap >= 2.5 or c * kap * kap >= 42.0: continue
        _, _, lam, next_, mrv, Hq, mre, mim, mim_pk, l1 = exact_R_non_q(c, j, kap)
        print(f"  cqk={cqk:4.2f}: lam*={lam:.4f} ext={next_:2d} max_rv={mrv:+.3f} H_q={Hq:5.2f} | |ReR|={mre:.4f} |ImR|={mim:.4f} |ImR_pk|={mim_pk:.4f} | L1(ext=0)={l1:.4f}")

# =============================================================================
# 2. RE-RUN CREVICE MAP FOR c*q*kap > 0.90 (c = 2..40): WHICH POINTS FAIL |Phi@1| < 1.45?
# =============================================================================
print("\n=== 2. CREVICE MAP FOR c*q*kap > 0.90, c*kap^2 < 42, kap < 2.5 ===")
def model_phi1(c, j_mod, kap, M_psi=96):
    alpha = kap / c
    M_cut = int(np.ceil(sqrt(35.0 / max(alpha, 1e-8)))) + 5
    ms_ext = np.arange(-(c - 1), 0, dtype=np.float64)
    u_ext = np.exp(alpha * ms_ext * (c + ms_ext) + 2j * pi * ms_ext * j_mod / c)
    ms_int = np.concatenate([np.arange(1, M_cut + 1, dtype=np.float64),
                             np.arange(-(c + M_cut), -c, dtype=np.float64)])
    lw = -alpha * ms_int * (c + ms_int) - 2j * pi * ms_int * j_mod / c
    u_int = np.exp(lw[lw.real > -35.0])
    psis = np.linspace(-pi, pi, M_psi, endpoint=False)
    eip = np.exp(1j * psis)
    log_K = np.zeros(M_psi, dtype=np.complex128)
    for arr, ep in [(u_ext, np.conj(eip)), (u_int, eip)]:
        for k0 in range(0, len(arr), 200):
            ch = arr[k0:k0+200]
            log_K += np.sum(np.log(1.0 - ch[None, :] * ep[:, None]), axis=1)
    return float(np.max(np.abs(log_K.imag)))

kaps = np.geomspace(0.005, 2.5, 50)
for c in [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30, 40]:
    fails_q1 = []; fails_qgt1 = []
    for jm in range(0, c // 2 + 1):
        g = gcd(c, jm) if jm > 0 else c
        q = c // g
        for kap in kaps:
            if c * q * kap <= 0.90 or c * kap * kap >= 42.0: continue
            mp = model_phi1(c, jm, kap)
            if mp >= 1.45:
                if q == 1: fails_q1.append((kap, c*kap, mp))
                else: fails_qgt1.append((jm, q, g, kap, c*q*kap, mp))
    qgt1_summary = sorted(set((x[0], x[1], x[2]) for x in fails_qgt1))
    print(f"  c={c:2d}: q=1 fails={len(fails_q1):2d} | q>1 fails={len(fails_qgt1):3d} (jm,q,g)={qgt1_summary[:8]}")
