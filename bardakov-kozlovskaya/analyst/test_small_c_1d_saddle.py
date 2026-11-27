"""
CRITICAL TEST OF THE TWO REMAINING THIN STRIPS (small c AND small a) AT HUGE n = 80, 160, 320, 640!

PART 1 (Small c = 2, 3, 4, 5, 10 at obtuse k1 = 0.7n and very obtuse k1 = 0.9n):
  In our Vandermonde-Cauchy integral:
    c_{a,b,c}^{(n)} = (1/(2pi)^2) int int Q(y, z) W(y, z) d(phi2) d(psi).
  Notice that if we do the exact psi-integral FIRST (or use a BETTER saddle radius Ry* != r2^n),
  Wait: for c = 2, [z^2] p_n(1, y, z) = (n/2) * (n*H1(y)^2 - H2(y)) / (1-y)^n = e_2(1/(1+zeta^k y^{1/n})^n) * (1-y)^n!
  More generally, for ANY c >= 1 and ANY fixed y = Ry * exp(i*phi2), we can compute
    E_{n-c}(1, y) = (-1)^c * (1-y)^n * e_c(lambda_0, ..., lambda_{n-1})
  in O(n*c) time via DP (elementary symmetric polynomial of lambda_k = 1/(1 + zeta^k y^{1/n})^n)!!
  WAIT A SECOND!!
  CHECK E_{n-c}(1, Ry * exp(i*phi2)) * exp(-i*k2*phi2) AS A FUNCTION OF phi2 IN [-pi, pi]!!
  1. For ODD c: check if its modulus peaks SHARPLY at phi2 = 0 with sign -1 (all coeffs < 0? Wait: E_{n-c}(1, y) has sign (-1)^{k2*(c+1)+c} -> for odd c, ALL POSITIVE coeffs in y -> sharp Gaussian peak at phi2 = 0!!)!
  2. For EVEN c: E_{n-c}(1, y) has ALTERNATING coeffs (-1)^{k2} in y -> sharp peak at phi2 = pi (i.e. y = -Ry < 0)!!
  OH!! MY!! GOD!!
  CHECK RIGHT NOW WHETHER AT FIXED c = 2, 3, 4, 5, 10 AND n = 80, 160, 320, 640,
  E_{n-c}(1, Ry* e^{i phi2}) * W_c(Ry* e^{i phi2}) HAS A SHARP GAUSSIAN PEAK WITH |avg|/avg|F| -> O(1/sqrt(n))
  OR EVEN BETTER: CHECK THE EXACT WIDTH AND TAIL DECAY AROUND THE PEAK phi2 = 0 (odd c) / phi2 = pi (even c)!!

PART 2 (Small a = 1, 2, 3 at c = 2, 3, 5, n/4 for n = 40, 80, 160):
  Test Chart II (e2 = -1 for even c, e2 = +1 for odd c) at the right R1, R3!
"""
import numpy as np

def eval_Ec_from_y(n, c, y):
    # Compute (-1)^c * [z^c] p_n(1, y, z) = (1 - y)^n * e_c(1 / (1 + zeta^k y^{1/n})^n)
    # Wait: return ACTUAL [z^c] p_n(1, y, z) = (-1)^c * (1 - y)^n * e_c(lambda_k)!
    # Use log-rescaling to avoid overflow/underflow for n up to 640!
    r = abs(y)**(1.0 / n)
    ph = np.angle(y)
    zetas = np.exp(2j * np.pi * np.arange(n) / n)
    w = 1.0 + r * np.exp(1j * ph / n) * zetas # shape (n,)
    log_lam = -n * np.log(w)
    # Shift all log_lam by max Re(log_lam) so exp(log_lam - shift) has max modulus 1!
    sh = np.max(log_lam.real)
    lams = np.exp(log_lam - sh)
    # Elementary symmetric polynomial e_0..e_c via DP
    e = np.zeros(c + 1, dtype=np.complex128)
    e[0] = 1.0
    for lam in lams:
        # e[j] += lam * e[j-1] for j = c down to 1
        e[1:] = e[1:] + lam * e[:-1]
    # Exact log of [z^c] p_n(1, y, z):
    # (-1)^c * (1 - y)^n * e[c] * exp(c * sh)
    log_val = (1j * np.pi * (c % 2)) + n * np.log(1.0 - y) + np.log(e[c]) + c * sh
    return log_val

def find_saddle_Ry_1d(n, k2, c):
    # Find saddle point Ry > 0 maximizing/minimizing Re(log_Ec(s * Ry)) - k2 * log(Ry)
    # where s = +1 if c is odd (peak at y = +Ry), s = -1 if c is even (peak at y = -Ry)!
    s = 1.0 if (c % 2 == 1) else -1.0
    # Start near Maillot r2^n
    k1 = n - k2 - c
    r2 = np.sin(np.pi * k2 / n) / np.sin(np.pi * k1 / n)
    u = n * np.log(max(min(r2, 0.999), 1e-6)) # log(Ry)
    for _ in range(25):
        h = 1e-4
        f0 = eval_Ec_from_y(n, c, s * np.exp(u)).real - k2 * u
        fp = eval_Ec_from_y(n, c, s * np.exp(u + h)).real - k2 * (u + h)
        fm = eval_Ec_from_y(n, c, s * np.exp(u - h)).real - k2 * (u - h)
        g = (fp - fm) / (2 * h)
        H = (fp - 2*f0 + fm) / (h * h)
        if abs(g) < 1e-8 or H <= 1e-12:
            break
        u -= np.clip(g / H, -2.0, 2.0)
    return float(np.exp(u)), float(H)

print("=== PART 1: Exact 1D y-circle profile of [z^c] p_n(1, y, z) at optimal saddle Ry* ===")
for n in [80, 160, 320, 640]:
    print(f"\n--- n = {n} ---")
    for lbl, frac1, c in [("c=1 obt 0.7", 0.7, 1), ("c=2 obt 0.7", 0.7, 2), ("c=2 v.obt 0.9", 0.9, 2),
                          ("c=3 obt 0.7", 0.7, 3), ("c=5 obt 0.7", 0.7, 5), ("c=10 obt 0.7", 0.7, 10)]:
        k1 = int(frac1 * n)
        k2 = n - k1 - c
        if k2 < c: continue
        Ry, curv = find_saddle_Ry_1d(n, k2, c)
        sgn_peak = 1.0 if (c % 2 == 1) else -1.0
        phi0 = 0.0 if (c % 2 == 1) else np.pi
        # Expected sign of M_{k1,k2,c} = (-1)^{k1*k2 + k1*c + k2*c}
        sgn_M = (-1)**(k1*k2 + k1*c + k2*c)
        # Evaluate on 200-point grid around circle phi2 in [phi0 - pi, phi0 + pi)
        phis = phi0 + np.linspace(-np.pi, np.pi, 200, endpoint=False)
        log_vals = np.array([eval_Ec_from_y(n, c, Ry * np.exp(1j * ph)) - 1j * k2 * ph for ph in phis])
        if sgn_M < 0: log_vals += 1j * np.pi
        # Peak is at index 100 (phi = phi0)
        pk = log_vals[100].real
        vals = np.exp(log_vals - pk)
        avg = np.mean(vals)
        abs_avg = np.mean(np.abs(vals))
        # Check modulus at phi0 + pi/2 and phi0 + pi (opposite side of circle!)
        mod_90 = float(np.abs(vals[150]))
        mod_180 = float(np.abs(vals[0]))
        # Check exact phase at peak (should be 0!)
        ang_pk = float(np.angle(vals[100]))
        print(f"  {lbl:14s} (k2={k2:3d}, Ry*={Ry:.4e}, H={curv:6.1f}): avg={avg.real:+.4f}, |avg|/avg|F|={abs(avg)/abs_avg:.4f}, |F(90deg)|={mod_90:.3e}, |F(180deg)|={mod_180:.3e}, ang_pk={ang_pk:+.2e}")
