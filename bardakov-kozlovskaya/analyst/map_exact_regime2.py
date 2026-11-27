"""
MAP THE EXACT FINITE-n REGIME 2 (kappa < 3.5) LANDSCAPE AT n = 1000 AND n = 5000!
Check:
1. For ALL interior points with kappa < 3.5 at n = 1000 (1850 pts):
   - How many satisfy max_psi |Phi_exact(psi)| < pi/2 at lambda = 1?
   - How many satisfy L1_phase := sum |B_r - A_r|/r < pi/2 at lambda = 1?
   - What is the EXACT boundary of failures in terms of (c, q, kap, c*kap, c*q*kap)?
2. For extreme finite-n points at n = 10000 (including c = 200, j = 201 [c/n = 0.02, kap = 0.79],
   and c = 2..20 with tiny kap, and c = 100, j = 200 [q=1, kap=0.39]):
   - What is max|Phi_exact| at lambda = 1?
   - What is rat_exact at lambda = 1 vs optimal lambda?
"""
import numpy as np
from math import gcd

def kappa_val(n, c, j):
    return 2.0 * np.pi * np.sin(j * np.pi / n) * np.sin((j + c) * np.pi / n) / np.sin(c * np.pi / n)

def exact_roots_log(n, c, j):
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

    w_ext = (1.0 + rj * np.exp(1j * ths[ext_mask])) / sj
    nz = np.abs(w_ext) > 0
    lw_ext = n * np.log(w_ext[nz]) + (1j*np.pi if sgn_j < 0 else 0.0)
    u_ext = np.exp(lw_ext)

    w_int = sj / (1.0 + rj * np.exp(1j * ths[int_mask]))
    lw_int = n * np.log(w_int) + (1j*np.pi if sgn_j < 0 else 0.0)
    keep = lw_int.real > -40.0
    u_int = np.exp(lw_int[keep])
    return u_ext, u_int

def eval_exact_at_1(n, c, j, M_psi=128):
    u_ext, u_int = exact_roots_log(n, c, j)
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    eip = np.exp(1j * psis)
    log_K = np.zeros(M_psi, dtype=np.complex128)
    for arr, ep in [(u_ext, np.conj(eip)), (u_int, eip)]:
        for k0 in range(0, len(arr), 200):
            ch = arr[k0:k0+200]
            log_K += np.sum(np.log(1.0 - ch[None, :] * ep[:, None]), axis=1)
    max_phi = float(np.max(np.abs(log_K.imag)))
    mid = 2.0 - 2.0 * np.cos(psis)
    sh = float(np.max(log_K.real))
    integ = mid * np.exp(log_K - sh)
    rat = float(np.mean(integ).real) / float(np.mean(np.abs(integ)))
    return max_phi, rat

print("=== SCAN ALL REGIME-2 POINTS (kap < 3.5) AT n = 1000 ===")
n = 1000
tot = 0; ok_phi = 0; fails = []
worst_phi = 0.0; worst_info = None
for c in range(2, int(0.11 * n) + 2):
    for j in range(c + 1, int(0.20 * n) + 2):
        if j > (n - c) // 2: continue
        kap = kappa_val(n, c, j)
        if kap >= 3.5: continue
        tot += 1
        mp, rat = eval_exact_at_1(n, c, j, M_psi=128)
        g = gcd(c, j); q = c // g
        if mp < np.pi / 2: ok_phi += 1
        else: fails.append((c, j, q, kap, c*kap, c*q*kap, mp, rat))
        if mp > worst_phi: worst_phi = mp; worst_info = (c, j, q, kap, c*kap, c*q*kap, mp, rat)
print(f"n=1000: {ok_phi}/{tot} have |Phi_exact|<pi/2 at lam=1! (fails={len(fails)})")
print(f"  worst |Phi|: c={worst_info[0]}, j={worst_info[1]}, q={worst_info[2]}, kap={worst_info[3]:.4f}, c*k={worst_info[4]:.3f}, c*q*k={worst_info[5]:.3f}, |Phi|={worst_info[6]:.4f}, rat@1={worst_info[7]:.4f}")
if fails:
    fails.sort(key=lambda x: -x[3])
    print("  Top 8 fails by kap:")
    for f in fails[:8]:
        print(f"    c={f[0]:3d} j={f[1]:3d} q={f[2]:3d} kap={f[3]:.4f} c*k={f[4]:.3f} c*q*k={f[5]:.3f} |Phi|={f[6]:.3f} rat@1={f[7]:.4f}")

print("\n=== EXTREME FINITE-n CORNERS AT n = 10000 ===")
n = 10000
test_pts = [
    (2, 3), (2, 4), (2, 5), (2, 10), (2, 20), (2, 40),
    (3, 4), (3, 6), (3, 9), (3, 15),
    (4, 5), (4, 8), (4, 12), (4, 20),
    (10, 11), (10, 20), (10, 30), (10, 50),
    (20, 21), (20, 40), (20, 60),
    (50, 51), (50, 100), (50, 150),
    (100, 101), (100, 150), (100, 200),
    (200, 201), (200, 299), (200, 323), (200, 400),
    (500, 501), (500, 750), (500, 809), (500, 1000),
]
for c, j in test_pts:
    kap = kappa_val(n, c, j)
    g = gcd(c, j); q = c // g
    mp, rat = eval_exact_at_1(n, c, j, M_psi=256)
    print(f"  c={c:3d} j={j:4d} q={q:3d} kap={kap:7.4f} c*k={c*kap:8.2f} | max|Phi@1|={mp:8.4f} | rat@1={rat:8.4f}")
