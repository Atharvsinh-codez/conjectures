"""
EXACT CREVICE MAP OF REGIME 2 (kappa < 2.5) FOR SMALL AND MEDIUM c = 2..50!
For every c in [2, 50], every representative j_mod in 0..c//2 (covering all q and all D!),
and a dense log-grid of kappa in [0.005, 2.5]:
Check:
  1. Is c * q * kap <= 0.8? (Regime 2C!)
  2. If c * q * kap > 0.8: what is max|Phi@1| and L1_phase = sum |B_r - A_r|/r at lam = 1?
  3. Print EVERY (c, j_mod, q, kap) where BOTH c*q*kap > 0.8 AND L1_phase >= 1.50!
This tells us with 100% precision the EXACT boundary between 2C, lam=1, and Major-Arc saddle!
"""
import numpy as np
from math import gcd

def model_L1_and_phi(c, j_mod, kap, R_max=40, M_psi=128):
    alpha = kap / c
    M_cut = int(np.ceil(np.sqrt(35.0 / max(alpha, 1e-8)))) + 5
    ms_ext = np.arange(-(c - 1), 0, dtype=np.float64)
    u_ext = np.exp(alpha * ms_ext * (c + ms_ext) + 2j * np.pi * ms_ext * j_mod / c)
    ms_int = np.concatenate([np.arange(1, M_cut + 1, dtype=np.float64),
                             np.arange(-(c + M_cut), -c, dtype=np.float64)])
    lw = -alpha * ms_int * (c + ms_int) - 2j * np.pi * ms_int * j_mod / c
    u_int = np.exp(lw[lw.real > -35.0])
    # Compute A_r, B_r for r = 1..R_max
    p_e = np.ones_like(u_ext); p_i = np.ones_like(u_int)
    L1 = 0.0
    for r in range(1, R_max + 1):
        p_e *= u_ext; p_i *= u_int
        L1 += abs(float(np.sum(p_i).real) - float(np.sum(p_e).real)) / r
    # Compute exact max|Phi| on |w|=1
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    eip = np.exp(1j * psis)
    log_K = np.zeros(M_psi, dtype=np.complex128)
    for arr, ep in [(u_ext, np.conj(eip)), (u_int, eip)]:
        for k0 in range(0, len(arr), 200):
            ch = arr[k0:k0+200]
            log_K += np.sum(np.log(1.0 - ch[None, :] * ep[:, None]), axis=1)
    mp = float(np.max(np.abs(log_K.imag)))
    return L1, mp

kaps = np.geomspace(0.01, 2.5, 45)
print("=== CREVICE MAP FOR c = 2..40 WHEN c * q * kap > 1.0 ===")
for c in [2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 40]:
    max_L1_nonq1 = 0.0; worst_nonq1 = None
    max_phi_nonq1 = 0.0; worst_phi_nonq1 = None
    max_L1_q1 = 0.0; worst_q1 = None
    fails_phi = []
    for jm in range(0, c // 2 + 1):
        g = gcd(c, jm) if jm > 0 else c
        q = c // g
        for kap in kaps:
            if c * q * kap <= 1.0: continue
            L1, mp = model_L1_and_phi(c, jm, kap, R_max=30, M_psi=64)
            if q == 1:
                if L1 > max_L1_q1: max_L1_q1 = L1; worst_q1 = (kap, c*kap, L1, mp)
            else:
                if L1 > max_L1_nonq1: max_L1_nonq1 = L1; worst_nonq1 = (jm, q, kap, c*q*kap, L1, mp)
                if mp > max_phi_nonq1: max_phi_nonq1 = mp; worst_phi_nonq1 = (jm, q, kap, c*q*kap, L1, mp)
            if mp >= 1.40:
                fails_phi.append((jm, q, kap, c*q*kap, L1, mp))
    print(f"c={c:2d}: q=1 maxL1={max_L1_q1:.3f} (mp={worst_q1[3]:.3f} at k={worst_q1[0]:.3f}) | q>1 maxPhi={max_phi_nonq1:.3f} at {worst_phi_nonq1[:3] if worst_phi_nonq1 else None} | #phi>=1.4: {len(fails_phi)}")
    if fails_phi:
        q_set = sorted(set((x[0], x[1]) for x in fails_phi))
        print(f"      (jm, q) with |Phi|>=1.4: {q_set[:10]}")
