"""
TEST THE EXACT BOUNDARY OF THE TWO ANALYTIC REGIMES IN THE (n, c, j) PLANE:
1. Check T(alpha; j, c) = [w^c] prod_{m in Z} (1 + e^{-2pi i m j/c} e^{-alpha m(c+m)} w)
   for ALL c = 2..12, j = 1..c, and alpha in [0.001, 2.0]!
   Verify if (-1)^j * T(alpha; j, c) > 0 ALWAYS (and compute exact ratios!).
2. For ALL (c, j) at n = 1000, 5000 with kappa < 4.0:
   Compare:
   - Regime A (unscaled |z| = s_j^n): when is unscaled rat ALWAYS >= 0.80?!
   - Regime B (Poisson power-sum dominance: q^2 * alpha <= C): when does the leading
     q-power term (1/g!) * ((-1)^{q-1} P_q / q)^g dominate ALL other partitions of c?!
   Check if Regime A and Regime B OVERLAP WITH HUGE MARGIN across the ENTIRE (n, c, j) plane!!
"""
import numpy as np
from math import gcd, factorial

print("=== PART 1: Verify (-1)^j * T(alpha; j, c) > 0 AND check exact Poisson power-sum dominance ===")
def eval_T_and_Poisson(alpha, j, c):
    # Compute exact power sums P_r = sum_{m=-M}^M e^{-2pi i r m j / c} e^{-alpha r m(c+m)} for r=1..c
    M = int(np.ceil(c/2 + np.sqrt(50.0 / max(alpha, 1e-6)))) + 10
    ms = np.arange(-M, M + 1, dtype=np.float64)
    P = np.zeros(c + 1, dtype=np.complex128)
    for r in range(1, c + 1):
        # shift exponent by -alpha*r*c^2/4 so all exponents <= 0 (NO overflow!)
        ex = -alpha * r * (ms + c/2.0)**2
        ph = -2.0 * np.pi * r * ms * j / c
        P[r] = np.sum(np.exp(ex + 1j * ph))
    # Compute e_c via Newton identities (all shifted by exp(+alpha*c^3/4)):
    # e_0 = 1, k*e_k = sum_{r=1}^k (-1)^{r-1} P_r e_{k-r}
    e = np.zeros(c + 1, dtype=np.complex128)
    e[0] = 1.0
    for k in range(1, c + 1):
        s = 0.0 + 0.0j
        for r in range(1, k + 1):
            s += ((-1)**(r - 1)) * P[r] * e[k - r]
        e[k] = s / k
    g = gcd(c, j); q = c // g
    # Leading q-only term: (1/g!) * ((-1)^{q-1} P_q / q)^g
    lead = ((((-1)**(q - 1)) * P[q] / q) ** g) / factorial(g)
    sgn_j = 1.0 if (j % 2 == 0) else -1.0
    return float((sgn_j * e[c]).real), float((sgn_j * lead).real), float(abs(e[c] - lead) / abs(lead))

# Sweep alpha from 0.001 to 1.0 for various (c, j)
for c, j in [(2, 1), (2, 2), (3, 1), (3, 3), (4, 1), (4, 2), (4, 4), (5, 1), (6, 1), (6, 2), (6, 3), (6, 6), (10, 1), (10, 2), (10, 5), (10, 10)]:
    g = gcd(c, j); q = c // g
    min_T = 1e9
    # Test small alpha where alpha * c * q <= 1.0 (Poisson regime) and larger alpha!
    pois_errs = []
    for alpha in [0.001, 0.005, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]:
        Tv, Lv, err = eval_T_and_Poisson(alpha, j, c)
        if Tv < min_T: min_T = Tv
        if alpha * c * q <= 0.5:
            pois_errs.append((alpha, alpha*c*q, err))
    max_perr = max((x[2] for x in pois_errs), default=-1.0)
    print(f"c={c:2d}, j={j:2d} (g={g:2d}, q={q:2d}): min ((-1)^j T) = {min_T:+.4e} | max Poisson rel_err (for alpha*c*q<=0.5) = {max_perr:.4f}")

print("\n=== PART 2: Check OVERLAP between unscaled |z|=s_j^n AND Poisson regime alpha*c*q <= C ===")
# Wait: check what unscaled |z|=s_j^n ratio is when alpha*c*q >= 0.5 (i.e. q*kappa >= 0.5)!!
def unscaled_rat(n, c, j, M_psi=128):
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    sin2 = 4.0 * np.sin(psis / 2.0)**2
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
    log_K = np.zeros(M_psi, dtype=np.complex128)
    if np.any(ext_mask):
        w_ext = (1.0 + rj * np.exp(1j * ths[ext_mask])) / sj
        nz = np.abs(w_ext) > 0
        if np.any(nz):
            lw = n * np.log(w_ext[nz]) + (1j*np.pi if sgn_j < 0 else 0.0)
            log_K += np.sum(np.log(1.0 - np.exp(lw[None, :] - 1j*psis[:, None])), axis=1)
    w_int = sj / (1.0 + rj * np.exp(1j * ths[int_mask]))
    lw = n * np.log(w_int) + (1j*np.pi if sgn_j < 0 else 0.0)
    keep = lw.real > -40.0
    if np.any(keep):
        log_K += np.sum(np.log(1.0 - np.exp(lw[keep][None, :] + 1j*psis[:, None])), axis=1)
    sh = np.max(log_K.real)
    integ = sin2 * np.exp(log_K - sh)
    val = np.mean(integ)
    abs_val = np.mean(np.abs(integ))
    return float(val.real / abs_val), float(val.real * np.exp(min(sh, 50.0)))

# Scan ALL (c, j) at n = 500 with kappa < 4.0 and check both q*kappa and unscaled_rat and Poisson err!
n = 500
min_rat_above = {0.5: 1e9, 1.0: 1e9, 1.5: 1e9, 2.0: 1e9}
worst_above = {}
max_perr_below = {0.5: 0.0, 1.0: 0.0, 1.5: 0.0, 2.0: 0.0}
worst_below = {}
for c in range(2, 60):
    for j in range(c + 1, 120):
        kap = 2.0 * np.pi * np.sin(j*np.pi/n) * np.sin((j+c)*np.pi/n) / np.sin(c*np.pi/n)
        if kap >= 4.0: continue
        g = gcd(c, j); q = c // g; alpha = kap / c
        qkap = q * kap # = alpha * c * q!
        rat, Iv = unscaled_rat(n, c, j, M_psi=64)
        _, _, perr = eval_T_and_Poisson(alpha, j, c)
        for thr in [0.5, 1.0, 1.5, 2.0]:
            if qkap >= thr:
                if rat < min_rat_above[thr]:
                    min_rat_above[thr] = rat
                    worst_above[thr] = (c, j, q, kap, qkap, rat)
            else:
                if perr > max_perr_below[thr]:
                    max_perr_below[thr] = perr
                    worst_below[thr] = (c, j, q, kap, qkap, perr)

print(f"At n={n}, scanning ALL (c,j) with kappa < 4.0:")
for thr in [0.5, 1.0, 1.5, 2.0]:
    print(f"  threshold q*kappa={thr:.1f}:")
    print(f"    ABOVE (unscaled |z|=s_j^n): min rat = {min_rat_above[thr]:.4f} at {worst_above.get(thr)}")
    print(f"    BELOW (Poisson leading q):  max rel_err = {max_perr_below[thr]:.4f} at {worst_below.get(thr)}")
