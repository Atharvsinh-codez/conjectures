"""
Map the EXACT boundary of |Phi_exact(psi)| < pi/2 and L1_phase < pi/2 across ALL (c, q, kappa)
at VERY LARGE n (up to n = 10^6), including the extreme scaling corner:
  q = 1, c = 100..1000, kappa = 0.01..3.5!
Also check what happens on the EXACT contour |w| = lam_* (and whether the integrand has
positive real part everywhere or only in a saddle peak, and what scaling parameter controls it!).
"""
import numpy as np
from math import gcd

def kappa_val(n, c, j):
    return 2.0 * np.pi * np.sin(j * np.pi / n) * np.sin((j + c) * np.pi / n) / np.sin(c * np.pi / n)

def test_model_phase(alpha, c, j, M_psi=512):
    # Pure Gaussian-phase model at lam = 1:
    # u_ext(m) = exp(alpha*m*(c+m) + 2pi i m j/c) for m = -(c-1)..-1
    # u_int(m) = exp(-alpha*m*(c+m) - 2pi i m j/c) for m >= 1 and m <= -(c+1)
    M_cut = int(np.ceil(np.sqrt(45.0 / max(alpha, 1e-12)))) + 5
    ms_ext = np.arange(-(c - 1), 0, dtype=np.float64)
    u_ext = np.exp(alpha * ms_ext * (c + ms_ext) + 2j * np.pi * ms_ext * j / c)
    ms_pos = np.arange(1, M_cut + 1, dtype=np.float64)
    ms_neg = np.arange(-(c + M_cut), -c, dtype=np.float64)
    ms_int = np.concatenate([ms_pos, ms_neg])
    lw_int = -alpha * ms_int * (c + ms_int) - 2j * np.pi * ms_int * j / c
    u_int = np.exp(lw_int[lw_int.real > -45.0])

    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    eip = np.exp(1j * psis)
    # Chunked sum to avoid huge arrays
    log_K = np.zeros(M_psi, dtype=np.complex128)
    for arr, sgn in [(u_ext, -1), (u_int, +1)]:
        ep = np.conj(eip) if sgn < 0 else eip
        for k0 in range(0, len(arr), 200):
            chunk = arr[k0:k0+200]
            log_K += np.sum(np.log(1.0 - chunk[None, :] * ep[:, None]), axis=1)
    max_phi = float(np.max(np.abs(log_K.imag)))
    # Also compute exact rat at lam = 1
    mid = 2.0 - 2.0 * np.cos(psis)
    sh = float(np.max(log_K.real))
    integ = mid * np.exp(log_K - sh)
    rat = float(np.mean(integ).real) / float(np.mean(np.abs(integ)))
    # First few A_r, B_r
    A1 = float(np.sum(u_ext).real)
    B1 = float(np.sum(u_int).real)
    g = gcd(c, j); q = c // g
    Aq = float(np.sum(u_ext**q).real)
    Bq = float(np.sum(u_int**q).real)
    return max_phi, rat, A1, B1, Aq, Bq

print("=== EXTREME SCALING CORNERS OF THE (c, q, kappa) SPACE AT lam = 1 ===")
print(f"{'c':>5} {'j':>6} {'q':>4} {'kap':>7} {'c*q*k':>8} {'c*q^2*k^2':>10} {'max|Phi|':>9} {'rat@1':>8} {'Aq':>9} {'Bq':>9} {'Bq-Aq':>9}")
cases = [
    # Small c, small kap (Poisson corner)
    (2, 3, 0.5), (2, 3, 0.2), (2, 3, 0.05),
    (2, 4, 1.5), (2, 4, 0.5), (2, 4, 0.1),
    (3, 4, 0.5), (3, 4, 0.1), (3, 6, 1.0), (3, 6, 0.3),
    (4, 8, 1.5), (4, 8, 0.5), (4, 8, 0.1),
    # Large c, q=1, moderate & small kap!
    (20, 40, 1.0), (20, 40, 0.3), (20, 40, 0.1), (20, 40, 0.05),
    (100, 200, 1.0), (100, 200, 0.3), (100, 200, 0.1), (100, 200, 0.03), (100, 200, 0.01),
    (500, 1000, 0.1), (500, 1000, 0.03), (500, 1000, 0.01), (500, 1000, 0.003),
    # Large c, q=2, small kap
    (100, 150, 0.3), (100, 150, 0.1), (100, 150, 0.03), (100, 150, 0.01),
    # Large c, q=c (coprime!), small kap
    (50, 51, 1.0), (50, 51, 0.3), (50, 51, 0.1), (50, 51, 0.03), (50, 51, 0.01),
    (200, 201, 0.3), (200, 201, 0.05), (200, 201, 0.01),
]
for c, j, kap in cases:
    g = gcd(c, j); q = c // g
    alpha = kap / c
    mp, rat, A1, B1, Aq, Bq = test_model_phase(alpha, c, j)
    print(f"{c:5d} {j:6d} {q:4d} {kap:7.4f} {c*q*kap:8.3f} {c*q*q*kap*kap:10.4f} {mp:9.4f} {rat:8.4f} {Aq:9.2f} {Bq:9.2f} {Bq-Aq:9.3f}")
