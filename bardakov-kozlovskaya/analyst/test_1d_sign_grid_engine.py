"""
CHECK THE 1D SIGN-GRID PROOF ENGINE FOR BARDAKOV-KOZLOVSKAYA CONJECTURE 3.3(2)!!
1. Check whether u = 0 has the EXACT SAME SIGN as u_c (so we ONLY need j = c+1 .. floor((n-c)/2),
   strictly inside r_j > s_j and r_j <= 1!) for ALL even n = 4..50 and ALL c >= 1!
   Wait: check the sign at u_c vs u=0 carefully!
2. Check the 1D integral kernel K(psi) = prod_{|th_k|>th*}(1 - R_k/z) * prod_{|th_k|<th*}(1 - z/R_k)
   at z = (-1)^j s_j^n e^{i psi} across ALL (n, c, j) for n = 8, 20, 40, 80, 160!
   Specifically compute:
   - I_{n,c,j} = (1/2pi) int 4 sin^2(psi/2) K(psi) d(psi) (MUST be > 0!)
   - min_psi Re(K(psi)), max_psi |K(psi) - 1|, and max_k |ratio_k|!
3. Check the EXACT 1D contour integral for f_{n,c}(u_j) on |z| = s_j^n (or |e3| = r_j^n s_j^n):
   compare the full integrand for f_{n,c}(u_j) with the pure E_{n-c}(1, u_j)/u_j^c integrand!
"""
import mpmath as mp
import numpy as np
import sys
sys.path.insert(0, r'd:\solving-conjectures\work\bardakov-kozlovskaya\verify')
import engine

mp.mp.dps = 60

def eval_poly_mp(cs, u):
    acc = mp.mpf(0)
    for c in reversed(cs): acc = acc * u + c
    return acc

def eval_fnc(C, n, c, u):
    d = n - 3 * c
    if abs(1 + u) < mp.mpf('1e-30'):
        return mp.mpf(0) if d % 2 == 1 else C[(0, d//2, c)] * (u**(d//2))
    t = u / ((1 + u)**2)
    g = [C[(d - 2*b, b, c)] for b in range(d // 2 + 1)]
    return eval_poly_mp(g, t) * ((1 + u)**d)

print("=== PART 1: Check sign(f_{n,c}(0)) vs sign(f_{n,c}(u_c)) and sign(f_{n,c}(u_{c+1})) ===")
same_as_uc = 0; opp_of_uc1 = 0; tot = 0
for n in [6, 8, 10, 12, 16, 20, 30, 40]:
    C = engine.compute_elementary_coeffs(n)
    for c in range(1, n // 3 + 1):
        d = n - 3 * c
        if d < 1: continue
        tot += 1
        s0 = 1 if C[(d, 0, c)] > 0 else -1
        rc = mp.sin(c * mp.pi / n) / mp.sin(2 * c * mp.pi / n)
        uc = ((-1)**c) * (rc**n)
        sc = 1 if eval_fnc(C, n, c, uc) > 0 else -1
        if s0 == sc: same_as_uc += 1
        if d >= 2:
            rc1 = mp.sin((c + 1) * mp.pi / n) / mp.sin((2 * c + 1) * mp.pi / n)
            uc1 = ((-1)**c) * (rc1**n)
            sc1 = 1 if eval_fnc(C, n, c, uc1) > 0 else -1
            if s0 == -sc1: opp_of_uc1 += 1
        else:
            opp_of_uc1 += 1
print(f"sign(f(0)) == sign(f(u_c)): {same_as_uc}/{tot} | sign(f(0)) == -sign(f(u_{{c+1}})): {opp_of_uc1}/{tot} !!")

print("\n=== PART 2: Inspect the 1D integral I_{n,c,j} = (1/2pi) int 4 sin^2(psi/2) K(psi) dpsi ===")
# Test across ALL (c, j) with c >= 1 and c+1 <= j <= (n-c)//2 for n = 8, 20, 40, 80, 160!
for n in [8, 20, 40, 80, 160]:
    min_I = 1e9
    max_I = -1e9
    min_ReK = 1e9
    max_mod_ratio = 0.0
    worst_I_info = None
    worst_K_info = None
    M_psi = 64
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    sin2 = 4.0 * np.sin(psis / 2.0)**2 # (M,)
    eip = np.exp(1j * psis) # (M,)
    for c in range(1, n // 3 + 1):
        d = n - 3 * c
        if d <= 0: continue
        # Check j from c+1 to (n-c)//2 (and also j = (n-c)//2 where r_j = 1 for even c!)
        j_start = c + 1
        j_end = (n - c) // 2
        if j_start > j_end: continue
        # Angles theta_k in (-pi, pi]:
        # For odd c: theta_k = (2k+1)*pi/n (k=0..n-1), wrapped to (-pi, pi]
        # For even c: theta_k = 2k*pi/n (k=0..n-1), wrapped to (-pi, pi]
        if c % 2 == 1:
            ths = ((2*np.arange(n) + 1) * np.pi / n + np.pi) % (2*np.pi) - np.pi
        else:
            ths = (2*np.arange(n) * np.pi / n + np.pi) % (2*np.pi) - np.pi
        th_star = np.pi * (n - c) / n
        # Strict masks (exclude the 2 boundary points |th| == th_star within 1e-12)
        ext_mask = np.abs(ths) > th_star + 1e-12  # size c - 1
        int_mask = np.abs(ths) < th_star - 1e-12  # size n - c - 1
        assert np.sum(ext_mask) == c - 1 and np.sum(int_mask) == n - c - 1
        for j in range(j_start, j_end + 1):
            rj = np.sin(j * np.pi / n) / np.sin((j + c) * np.pi / n)
            sj = np.sin(c * np.pi / n) / np.sin((j + c) * np.pi / n)
            # Compute log(R_k / ((-1)^j sj^n)) safely!
            # R_k = (1 + rj * exp(i * th_k))^n
            w_ext = (1.0 + rj * np.exp(1j * ths[ext_mask])) / sj
            w_int = sj / (1.0 + rj * np.exp(1j * ths[int_mask]))
            # Note |w_ext| < 1 and |w_int| < 1!
            # Special case: if rj == 1 and th == pi (even c, j=(n-c)/2), 1 + rj*e^{i*pi} = 0 -> w_ext = 0!
            rho_max = max(np.max(np.abs(w_ext)) if len(w_ext) else 0.0,
                          np.max(np.abs(w_int)) if len(w_int) else 0.0)
            if rho_max**n > max_mod_ratio: max_mod_ratio = float(rho_max**n)
            # Compute exact R_k / ((-1)^j sj^n) = ((-1)^j) * w_ext^n
            sgn_j = 1.0 if (j % 2 == 0) else -1.0
            # (M, c-1) and (M, n-c-1)
            if len(w_ext):
                fac_ext = 1.0 - (sgn_j * (w_ext**n))[None, :] / eip[:, None]
                prod_ext = np.prod(fac_ext, axis=1)
            else:
                prod_ext = np.ones(M_psi, dtype=np.complex128)
            fac_int = 1.0 - (sgn_j * (w_int**n))[None, :] * eip[:, None]
            prod_int = np.prod(fac_int, axis=1)
            K = prod_ext * prod_int
            I_val = float(np.mean(sin2 * K).real)
            rk_min = float(np.min(K.real))
            if I_val < min_I:
                min_I = I_val; worst_I_info = (c, j, I_val, float(rho_max**n))
            if I_val > max_I: max_I = I_val
            if rk_min < min_ReK:
                min_ReK = rk_min; worst_K_info = (c, j, rk_min, float(rho_max**n))
    print(f"n={n:3d}: I in [{min_I:+.5f}, {max_I:+.5f}] (worst={worst_I_info}) | min Re(K)={min_ReK:+.5f} (worst={worst_K_info}) | max |ratio|^n={max_mod_ratio:.5f}")
