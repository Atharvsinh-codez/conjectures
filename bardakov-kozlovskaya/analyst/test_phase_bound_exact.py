"""
TEST THE EXACT LOG-KERNEL PHASE BOUND AT lambda = 1:
On |w| = 1, the EXACT finite-n integrand is:
  K_exact(psi) = (2 - 2*cos(psi)) * exp( -sum_{r>=1} (A_r + B_r)/r * cos(r*psi) ) * exp( i * Phi_exact(psi) )
where:
  A_r = sum_{m = -(c-1)}^{-1} a_m^{-r}   (EXACTLY REAL by conjugate symmetry!)
  B_r = sum_{m not in [-c, 0]} a_m^{+r}  (EXACTLY REAL by conjugate symmetry!)
  Phi_exact(psi) = -sum_{r>=1} (B_r - A_r)/r * sin(r*psi).

If |Phi_exact(psi)| < pi/2 for ALL psi in [-pi, pi], then Re(K_exact(psi)) >= 0 POINTWISE EVERYWHERE,
which IMMEDIATELY PROVES (-1)^j * s_j^{n*c} * e_c(lam_0,...,lam_{n-1}) > 0!!

Let's check:
1. For what (n, c, j) does L1_phase := sum_{r>=1} |B_r - A_r|/r < pi/2 hold?
2. For what (n, c, j) does max_psi |Phi_exact(psi)| < pi/2 hold?
3. What about ALL kappa (not just kappa < 3.5! Could this ALSO prove Regime 1?!)?!
4. When L1_phase >= pi/2, what is c * q * kappa and kappa?
"""
import numpy as np
from math import gcd

def kappa_val(n, c, j):
    return 2.0 * np.pi * np.sin(j * np.pi / n) * np.sin((j + c) * np.pi / n) / np.sin(c * np.pi / n)

def exact_phase_stats(n, c, j, R_max=60, M_psi=256):
    if c % 2 == 1:
        ths = ((2*np.arange(n) + 1) * np.pi / n + np.pi) % (2*np.pi) - np.pi
    else:
        ths = (2*np.arange(n) * np.pi / n + np.pi) % (2*np.pi) - np.pi
    th_star = np.pi * (n - c) / n
    ext_mask = np.abs(ths) > th_star + 1e-12
    int_mask = np.abs(ths) < th_star - 1e-12
    rj = np.sin(j * np.pi / n) / np.sin((j + c) * np.pi / n)
    sj = np.sin(c * np.pi / n) / np.sin((j + c) * np.pi / n)
    sgn_j = 1.0 if (j % 2 == 0) else -1.0

    # u_ext = a_m^{-1} for ext, u_int = a_m for int (all have |u| < 1)
    w_ext = (1.0 + rj * np.exp(1j * ths[ext_mask])) / sj
    nz = np.abs(w_ext) > 0
    lw_ext = n * np.log(w_ext[nz]) + (1j*np.pi if sgn_j < 0 else 0.0)
    u_ext = np.exp(lw_ext) # |u_ext| < 1

    w_int = sj / (1.0 + rj * np.exp(1j * ths[int_mask]))
    lw_int = n * np.log(w_int) + (1j*np.pi if sgn_j < 0 else 0.0)
    keep = lw_int.real > -40.0
    u_int = np.exp(lw_int[keep]) # |u_int| < 1

    # Compute A_r, B_r for r = 1..R_max
    rs = np.arange(1, R_max + 1)
    # Vectorized powers
    A_r = np.zeros(R_max)
    B_r = np.zeros(R_max)
    p_ext = np.ones_like(u_ext)
    p_int = np.ones_like(u_int)
    for idx in range(R_max):
        p_ext = p_ext * u_ext
        p_int = p_int * u_int
        A_r[idx] = float(np.sum(p_ext).real)
        B_r[idx] = float(np.sum(p_int).real)
    L1_phase = float(np.sum(np.abs(B_r - A_r) / rs))
    L1_sum_only = float(np.sum((np.abs(A_r) + np.abs(B_r)) / rs))

    # Exact Phi_exact(psi) via exact log(1 - u e^{+-i psi})
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    eip = np.exp(1j * psis)
    log_K = np.zeros(M_psi, dtype=np.complex128)
    if len(u_ext) > 0:
        log_K += np.sum(np.log(1.0 - u_ext[None, :] * np.conj(eip)[:, None]), axis=1)
    if len(u_int) > 0:
        log_K += np.sum(np.log(1.0 - u_int[None, :] * eip[:, None]), axis=1)
    max_phi = float(np.max(np.abs(log_K.imag)))
    # Also exact ratio on |w|=1
    mid = 2.0 - 2.0 * np.cos(psis)
    sh = float(np.max(log_K.real))
    integ = mid * np.exp(log_K - sh)
    rat = float(np.mean(integ).real) / float(np.mean(np.abs(integ)))
    return max_phi, L1_phase, L1_sum_only, rat, A_r[:4], B_r[:4]

print("=== CHECK EXACT PHASE |Phi_exact(psi)| AT lambda = 1 ACROSS ALL INTERIOR POINTS ===")
for n in [40, 80, 160, 320]:
    tot = 0
    ok_phi = 0
    ok_L1 = 0
    ok_naive = 0
    max_phi_all = 0.0; worst_phi = None
    max_L1_kap2 = 0.0; worst_L1_kap2 = None
    fails_phi = []
    for c in range(2, (2*n - 3)//3 + 1):
        j_max = n//2 - (c + 1)//2 - (1 if c % 2 == 0 else 0)
        for j in range(c + 1, j_max + 1):
            tot += 1
            kap = kappa_val(n, c, j)
            g = gcd(c, j); q = c // g
            mp, l1, l1n, rat, Ar, Br = exact_phase_stats(n, c, j, R_max=40, M_psi=128)
            if mp < np.pi / 2: ok_phi += 1
            else: fails_phi.append((c, j, q, kap, c*q*kap, mp, l1, rat))
            if l1 < np.pi / 2: ok_L1 += 1
            if l1n < np.pi / 2: ok_naive += 1
            if mp > max_phi_all: max_phi_all = mp; worst_phi = (c, j, q, kap, c*q*kap, mp, l1, rat)
            if kap >= 2.0 and l1 > max_L1_kap2:
                max_L1_kap2 = l1; worst_L1_kap2 = (c, j, kap, mp, l1, l1n)
    print(f"n={n:3d} ({tot:5d} pts): |Phi|<pi/2: {ok_phi}/{tot} | L1<pi/2: {ok_L1}/{tot} | naive<pi/2: {ok_naive}/{tot}")
    print(f"   worst |Phi|: {worst_phi}")
    print(f"   max L1 for kap>=2.0: {max_L1_kap2:.4f} at {worst_L1_kap2}")
    if fails_phi:
        fails_phi.sort(key=lambda x: -x[4]) # sort by c*q*kap descending
        print(f"   |Phi|>=pi/2 count={len(fails_phi)}, MAX c*q*kap among fails={fails_phi[0][4]:.4f} at (c={fails_phi[0][0]}, j={fails_phi[0][1]}, q={fails_phi[0][2]}, kap={fails_phi[0][3]:.3f}, |Phi|={fails_phi[0][5]:.3f})")
        fails_phi.sort(key=lambda x: -x[3]) # sort by kap descending
        print(f"   MAX kap among |Phi|>=pi/2 fails={fails_phi[0][3]:.4f} at (c={fails_phi[0][0]}, j={fails_phi[0][1]}, q={fails_phi[0][2]}, c*q*kap={fails_phi[0][4]:.3f})")
