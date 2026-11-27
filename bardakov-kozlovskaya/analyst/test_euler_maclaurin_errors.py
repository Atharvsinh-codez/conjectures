"""
MEASURE EXACT EULER-MACLAURIN & JACOBI-THETA ERROR CONSTANTS!
In our exact decomposition of Q(y, z) at y = r2^n e^{i phi2}, z = r3^n e^{i psi} (r2=sin(t2)/sin(t1), r3=sin(t3)/sin(t1)):
  1. Interior phase error: Delta_phi(phi2) = sum_{k in K} n*arg(1 + r2 e^{i theta_k}) - k2*phi2 (mod 2pi).
     We proved BOTH B1 and B2 terms vanish identically -> O(1/n)!
     Measure max |Delta_phi(phi2)| across ALL (k1, k2, k3) for n = 12, 24, 48, 80!
  2. Look at the exact phase of Pi_exact(phi2, psi) = prod_{k in K}(1 - z/(1+r2 e^{i th_k})^n) * prod_{k notin K}(1 - (1+r2 e^{i th_k})^n/z)!
     WAIT A SECOND!!
     WHY DID Re(sgn * Q(y, z)) STAY >= 0 EVEN AT n = 4, 6, 8, 10, 12?!
     COMPUTE THE EXACT ARGUMENT arg(sgn * Q(y, z)) ACROSS THE ENTIRE TORUS AT n = 4, 6, 12, 24!
     If |arg(sgn * Q(y, z))| <= pi/2 EVERYWHERE on the torus, then Re(sgn * Q) >= 0 EVERYWHERE!
     Let's check max |arg(sgn * Q(y, z))| (in units of pi/2)!!
"""
import numpy as np

print("=== Check max |arg(sgn * Q(y, z))| / (pi/2) across ALL (k1, k2, k3) at n = 4, 6, 8, 12, 20, 30 ===")
for n in [4, 6, 8, 12, 20, 30]:
    zetas = np.exp(2j * np.pi * np.arange(n) / n)
    max_arg_all = 0.0
    max_arg_info = None
    max_phase_err = 0.0
    max_phase_info = None
    g = np.linspace(-np.pi + 0.0137, np.pi + 0.0137, 25, endpoint=False)
    P2, P3 = np.meshgrid(g, g, indexing='ij')
    p2f, p3f = P2.ravel(), P3.ravel()
    for k3 in range(1, n // 3 + 1):
        for k2 in range(k3, (n - k3) // 2 + 1):
            k1 = n - k2 - k3
            sgn = (-1)**(k1*k2 + k1*k3 + k2*k3)
            r2 = np.sin(np.pi * k2 / n) / np.sin(np.pi * k1 / n)
            r3 = np.sin(np.pi * k3 / n) / np.sin(np.pi * k1 / n)
            y_v = (r2**n) * np.exp(1j * p2f)
            z_v = (r3**n) * np.exp(1j * p3f)
            beta_v = r2 * np.exp(1j * p2f / n)
            factors = z_v[:, None] - (1.0 + beta_v[:, None] * zetas[None, :])**n
            # Use log sum to avoid overflow at n=30!
            # Scale each factor by max(|z|, |1+beta*zeta|^n)
            terms2 = (1.0 + beta_v[:, None] * zetas[None, :])**n
            # Actually at n=30, |1+r2|^30 <= 2^30 = 1e9, completely safe in float64!
            prods = np.prod(factors, axis=1)
            vals = sgn * prods / (y_v**k2 * z_v**k3)
            angs = np.abs(np.angle(vals))
            m_ang = float(np.max(angs))
            if m_ang > max_arg_all:
                max_arg_all = m_ang
                idx = int(np.argmax(angs))
                max_arg_info = (k1, k2, k3, float(p2f[idx]), float(p3f[idx]), float(np.abs(vals[idx])/np.mean(np.abs(vals))))
            # Also check interior phase error Delta_phi(phi2)
            th_star = np.pi * (n - k3) / n
            for ph2 in g:
                ths = ((2*np.pi*np.arange(n) + ph2) / n + np.pi) % (2*np.pi) - np.pi
                mask = np.abs(ths) < th_star
                if np.sum(mask) == n - k3:
                    S = n * np.sum(np.angle(1.0 + r2 * np.exp(1j * ths[mask])))
                    err = abs(((S - k2 * ph2) + np.pi) % (2*np.pi) - np.pi)
                    if err > max_phase_err:
                        max_phase_err = err
                        max_phase_info = (k1, k2, k3, float(ph2))
    print(f"n={n:2d}: max |arg(sgn*Q)|/(pi/2) = {max_arg_all/(np.pi/2):.6f} at k={max_arg_info[:3]} (|Q|/avg={max_arg_info[5]:.2e}) | max |Delta_phi| = {max_phase_err:.5f} rad at k={max_phase_info[:3]}")
