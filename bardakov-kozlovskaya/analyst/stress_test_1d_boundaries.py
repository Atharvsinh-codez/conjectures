"""
STRESS-TEST THE TWO CRITICAL ASYMPTOTIC BOUNDARIES OF THE 1D SIGN-GRID PROOF:
1. Edge Case 1 (d = 2, u_j = +1): n = 3c + 2, j = c + 1 = (n-c)/2, for c = 10, 26, 50, 100, 200 (n up to 602)!
   Also d = 3 (u_j close to -1): n = 3c + 3, j = c + 1 = (n-c-1)/2!
2. Edge Case 2 (Fixed small c, j as n -> inf):
   Compute I_{n,c,j} = (1/2pi) int 4 sin^2(psi/2) K(psi) dpsi for (c, j) = (2,3), (2,4), (2,5), (3,4), (3,5), (4,5), (5,6)
   at HUGE n = 160, 320, 640, 1280, 2560, 5120!!
   (Note: for fixed (c, j) and large n, u_j = (-1)^c (j/(j+c))^n and z = (c/(j+c))^n are BOTH exponentially tiny,
    so x(z) = 1 + O(e^{-an}), y(z) = u_j(1 + O(e^{-bn})), meaning I_full == I_{n,c,j} to exponential precision!)
"""
import numpy as np

def test_I_kernel(n, c, j, M_psi=128):
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    sin2 = 4.0 * np.sin(psis / 2.0)**2
    eip = np.exp(1j * psis)
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
    # Use log-sum to avoid any float64 underflow/overflow at n=5120!
    log_K = np.zeros(M_psi, dtype=np.complex128)
    if np.any(ext_mask):
        w_ext = (1.0 + rj * np.exp(1j * ths[ext_mask])) / sj
        # skip exact 0 (r_j=1, th=pi)
        nz = np.abs(w_ext) > 0
        if np.any(nz):
            lw = n * np.log(w_ext[nz]) + (1j*np.pi if sgn_j < 0 else 0.0) # shape (k,)
            # factor: 1 - exp(lw - i*psi)
            fac = 1.0 - np.exp(lw[None, :] - 1j*psis[:, None])
            log_K += np.sum(np.log(fac), axis=1)
    w_int = sj / (1.0 + rj * np.exp(1j * ths[int_mask]))
    lw = n * np.log(w_int) + (1j*np.pi if sgn_j < 0 else 0.0)
    # Only include terms with Re(lw) > -40 (others contribute < 1e-17)
    keep = lw.real > -40.0
    if np.any(keep):
        fac = 1.0 - np.exp(lw[keep][None, :] + 1j*psis[:, None])
        log_K += np.sum(np.log(fac), axis=1)
    K = np.exp(log_K)
    val = np.mean(sin2 * K)
    abs_val = np.mean(sin2 * np.abs(K))
    return float(val.real), float(val.real / abs_val), float(np.max(np.abs(K)))

print("=== EDGE CASE 2: Fixed small (c, j) at HUGE n = 160..5120 ===")
for (c, j) in [(1, 2), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5), (5, 6)]:
    parts = []
    for n in [160, 320, 640, 1280, 2560, 5120]:
        I_v, rat, mxK = test_I_kernel(n, c, j, M_psi=128)
        parts.append(f"n={n}:{I_v:+.3e}({rat:+.2f})")
    print(f"(c={c},j={j}): " + "  ".join(parts))

print("\n=== EDGE CASE 1: d = 2 (n = 3c+2, j = c+1, u_j = +1) AND d = 3 (n = 3c+3, j = c+1) up to n = 602 ===")
def test_edge1_full(n, c, j, M_psi=128):
    psis = np.linspace(-np.pi, np.pi, M_psi, endpoint=False)
    eip = np.exp(1j * psis)
    zetas = np.exp(2j * np.pi * np.arange(n) / n)
    if c % 2 == 1:
        ths = ((2*np.arange(n) + 1) * np.pi / n + np.pi) % (2*np.pi) - np.pi
    else:
        ths = (2*np.arange(n) * np.pi / n + np.pi) % (2*np.pi) - np.pi
    th_star = np.pi * (n - c) / n
    int_mask = np.abs(ths) < th_star - 1e-12
    rj = np.sin(j * np.pi / n) / np.sin((j + c) * np.pi / n)
    sj = np.sin(c * np.pi / n) / np.sin((j + c) * np.pi / n)
    uj = ((-1)**c) * (rj**n)
    sgn_j = 1.0 if (j % 2 == 0) else -1.0
    z_v = sgn_j * (sj**n) * eip
    b_q = 1.0 + uj - z_v
    c_q = (1.0 - z_v) * (uj - z_v)
    disc = np.sqrt((1.0 - uj)**2 + 2.0*(1.0 + uj)*z_v - 3.0*z_v*z_v)
    x_v = 0.5 * (b_q + disc)
    swap = np.abs(x_v) < np.abs(0.5 * (b_q - disc))
    x_v = np.where(swap, 0.5 * (b_q - disc), x_v)
    y_v = c_q / x_v
    ax = np.exp(np.log(x_v) / n)
    by = np.exp(np.log(y_v) / n)
    w = ax[:, None] + by[:, None] * zetas[None, :]
    log_w = n * np.log(w)
    log_z = np.log(z_v)[:, None]
    use_z = log_z.real > log_w.real
    log_fac = np.where(use_z,
                       log_z + np.log(1.0 - np.exp(np.where(use_z, log_w - log_z, -100.0))),
                       log_w + 1j*np.pi + np.log(1.0 - np.exp(np.where(~use_z, log_z - log_w, -100.0))))
    log_pn = np.sum(log_fac, axis=1)
    e3 = z_v * (1.0 - z_v) * (uj - z_v)
    e3p = (1.0 - z_v)*(uj - z_v) - z_v*(uj - z_v) - z_v*(1.0 - z_v)
    log_meas = np.log(z_v) + np.log(e3p) - (c + 1) * np.log(e3)
    log_ref = np.sum(n * np.log(np.abs(1.0 + rj * np.exp(1j * ths[int_mask])))) + n * np.log(sj) - c * n * np.log(rj)
    exp_sgn = 1.0 if (j % 2 == 0) else -1.0
    integ = np.exp(log_pn + log_meas - log_ref + (1j*np.pi if exp_sgn < 0 else 0.0))
    return float(np.mean(integ).real), float(np.mean(integ).real / np.mean(np.abs(integ))), float(np.min(integ.real))

for c in [10, 26, 50, 100, 200]:
    n = 3*c + 2 # even since c is even!
    I_f, rat, mn = test_edge1_full(n, c, c + 1)
    print(f"d=2 (n={n:3d}, c={c:3d}, j={c+1:3d}, u=+1): I_full={I_f:+.5f}, ratio={rat:+.4f}, minRe={mn:+.4f}")

for c in [9, 25, 49, 99, 199]:
    n = 3*c + 3 # even since c is odd!
    I_f, rat, mn = test_edge1_full(n, c, c + 1)
    print(f"d=3 (n={n:3d}, c={c:3d}, j={c+1:3d}, u<0 ): I_full={I_f:+.5f}, ratio={rat:+.4f}, minRe={mn:+.4f}")
