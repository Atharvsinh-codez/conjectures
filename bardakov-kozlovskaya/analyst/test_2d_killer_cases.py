"""
CRITICAL COMPARISON:
1. Check the EXACT 2D torus integral for the monomial coefficient M_{k1, k2, k3} (and e-coeff c_{a,b,c})
   at the EXACT killer Regime-2 cases where 1D had Diophantine resonances:
     - (k1, k2, k3) = (995, 3, 2)  [c=2, j=3, n=1000, kap=0.148]
     - (k1, k2, k3) = (994, 4, 2)  [c=2, j=4, n=1000, kap=0.237]
     - (k1, k2, k3) = (9995, 3, 2) [c=2, j=3, n=10000, kap=0.015]
     - (k1, k2, k3) = (9994, 4, 2) [c=2, j=4, n=10000, kap=0.024]
     - (k1, k2, k3) = (9970, 20, 10) [c=10, j=20, n=10000, kap=0.118]
     - (k1, k2, k3) = (9599, 201, 200) [c=200, j=201, n=10000, kap=0.793]
   Recall the EXACT 2D torus integral for M_{k1, k2, k3} (with x=1, y=r^n e^{i phi2}, z=s^n e^{i psi}):
     (-1)^{k1*k2 + k1*k3 + k2*k3} M_{k1, k2, k3} * r^{n k2} s^{n k3}
     = (1/(2pi)^2) int_{-pi}^{pi} int_{-pi}^{pi} e^{-i k2 phi2 - i k3 psi} p_n(1, r^n e^{i phi2}, s^n e^{i psi}) dphi2 dpsi
   where r = sin(k2*pi/n)/sin(k1*pi/n), s = sin(k3*pi/n)/sin(k1*pi/n)!
   Wait: check the EXACT sign and phase at phi2 = c*pi (or k3*pi), psi = k2*pi!
2. Also inspect the peak structure of 1D F(lam_L2 e^{i psi}) at c=200, j=267, kap=0.01!
"""
import numpy as np

def test_2d_monomial(n, k2, k3, M=64):
    k1 = n - k2 - k3
    r = np.sin(k2 * np.pi / n) / np.sin(k1 * np.pi / n)
    s = np.sin(k3 * np.pi / n) / np.sin(k1 * np.pi / n)
    # Wait: p_n(1, y, z) = prod_{m=0}^{n-1} (z - (1 + y^{1/n} e^{2pi i m/n})^n)
    # Set y = r^n e^{i phi2} => y^{1/n} = r e^{i phi2/n}.
    # Set z = s^n e^{i psi}.
    # Shift phi2 = k3*pi + dphi2, psi = k2*pi + dpsi?
    # Let's just evaluate on the FULL torus (phi2, psi) in [-pi, pi)^2!
    phi2s = np.linspace(-np.pi, np.pi, M, endpoint=False)
    psis = np.linspace(-np.pi, np.pi, M, endpoint=False)
    ms = np.arange(n, dtype=np.float64)
    # For each phi2, compute log|...| and phase of prod_m (s^n e^{i psi} - (1 + r e^{i(phi2 + 2pi m)/n})^n)
    # Wait: factor out (1 + r e^{i th})^n when |...| > s^n, and s^n e^{i psi} when |...| < s^n!
    log_vals = np.zeros((M, M), dtype=np.complex128)
    for i, phi2 in enumerate(phi2s):
        ths = (phi2 + 2.0 * np.pi * ms) / n
        # complex numbers w_m = (1 + r e^{i th}) / s
        w = (1.0 + r * np.exp(1j * ths)) / s
        lw = n * np.log(w) # log((1 + r e^{i th})^n / s^n)
        # For each m: z - R_m = s^n * (e^{i psi} - e^{lw})
        # If lw.real > 0: factor out -e^{lw} => log(-e^{lw}) + log(1 - e^{-lw + i psi})
        # If lw.real <= 0: factor out e^{i psi} => i*psi + log(1 - e^{lw - i psi})
        big = lw.real > 0.0
        # Constant in psi:
        # Wait: add n*log(s) later (drops out in ratio!)
        acc = np.zeros(M, dtype=np.complex128)
        if np.any(big):
            lb = lw[big]
            acc += np.sum(lb) + 1j * np.pi * len(lb)
            inv_e = np.exp(-lb)
            # keep only inv_e with |inv_e| > 1e-16
            kp = np.abs(inv_e) > 1e-15
            if np.any(kp):
                ep = np.exp(1j * psis)
                # chunk
                arr = inv_e[kp]
                for k0 in range(0, len(arr), 200):
                    ch = arr[k0:k0+200]
                    acc += np.sum(np.log(1.0 - ch[None, :] * ep[:, None]), axis=1)
        sm = ~big
        if np.any(sm):
            ls = lw[sm]
            acc += len(ls) * 1j * psis
            e_s = np.exp(ls)
            kp = np.abs(e_s) > 1e-15
            if np.any(kp):
                epm = np.exp(-1j * psis)
                arr = e_s[kp]
                for k0 in range(0, len(arr), 200):
                    ch = arr[k0:k0+200]
                    # handle near-1 safely
                    near1 = np.abs(ch) > 1.0 - 1e-9
                    if np.any(near1):
                        p1 = np.prod(1.0 - ch[near1][None, :] * epm[:, None], axis=1)
                        with np.errstate(divide='ignore'):
                            acc += np.where(np.abs(p1) > 0, np.log(p1), -1e9)
                    far = ~near1
                    if np.any(far):
                        acc += np.sum(np.log(1.0 - ch[far][None, :] * epm[:, None]), axis=1)
        # Multiply by e^{-i k2 phi2 - i k3 psi} and target sign (-1)^{k1*k2 + k1*k3 + k2*k3}:
        sgn_par = (k1*k2 + k1*k3 + k2*k3) % 2
        log_vals[i, :] = acc - 1j * k2 * phi2 - 1j * k3 * psis + (1j * np.pi if sgn_par else 0.0)

    sh = float(np.max(log_vals.real))
    F = np.exp(log_vals - sh)
    val = float(np.mean(F).real)
    im = float(np.mean(F).imag)
    l1 = float(np.mean(np.abs(F)))
    min_re = float(np.min(F.real))
    # Also check peak location and active phase!
    active = (log_vals.real - sh) > -4.6
    max_arg_act = float(np.max(np.abs(np.angle(F[active]))))
    return val / l1, min_re / l1, max_arg_act, float(np.mean(active)), abs(im) / l1

print("=== EXACT 2D TORUS INTEGRAL FOR MONOMIALS AT KILLER REGIME-2 CASES ===")
print(f"{'n':>6} {'k2':>4} {'k3':>4} {'kap':>7} | {'rat_2D':>9} {'minRe/L1':>10} {'argAct':>8} {'fracAct':>8} {'im_err':>8}")
for n, k2, k3 in [
    (60, 17, 1),
    (1000, 3, 2),
    (1000, 4, 2),
    (1000, 6, 3),
    (1000, 20, 10),
    (10000, 3, 2),
    (10000, 4, 2),
    (10000, 20, 10),
    (10000, 201, 200),
]:
    kap = 2*np.pi*np.sin(k2*np.pi/n)*np.sin((k2+k3)*np.pi/n)/np.sin(k3*np.pi/n)
    M = 48 if n <= 1000 else 32
    rat, minr, aa, fa, ime = test_2d_monomial(n, k2, k3, M=M)
    print(f"{n:6d} {k2:4d} {k3:4d} {kap:7.4f} | {rat:9.5f} {minr:10.5f} {aa:8.4f} {fa:8.4f} {ime:8.1e}")
