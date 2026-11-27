"""
Deep probe of Gap G corners (including r0 = 1, 2, 3, 5 with D = 0 and D != 0, and huge c = 5000)!
Inspects:
  1. Exact lam_L2 and ratio Re(int)/int|I| at lam_L2 vs lam = 1 vs explicit formula lam_exp.
  2. Number and location of peaks of log|I(phi)| on [-pi, pi].
  3. Exact phase |Phi(phi)| across the support of the peaks (where log|I(phi)| >= max - 15).
  4. Whether A_r(lam) - B_r(lam) (exact Fourier phase expansion!) or Poisson summation on log(1 - a_k w)
     gives a 100% unified elementary bound WITHOUT major/minor arcs!
"""
import math
import cmath

def gaussian_a_list(c, j, kap):
    alpha = math.pi * kap / c
    x_max = int(math.sqrt(0.25 * c * c + 45.0 / alpha)) + 5
    shift = 0.0 if (c % 2 == 0) else 0.5
    roots = []
    sgn = 1.0 if (j % 2 == 0) else -1.0
    for k in range(-x_max, x_max + 1):
        x = k + shift
        log_r = alpha * (0.25 * c * c - x * x)
        ang = 2.0 * math.pi * x * j / c + (0.0 if sgn > 0 else math.pi)
        roots.append((log_r, ang))
    return roots

def eval_F_log(roots, c, lam, phi):
    log_lam = math.log(lam)
    log_abs = -c * log_lam
    tot_ang = -c * phi
    for log_r, ang in roots:
        lm = log_r + log_lam
        th = ang + phi
        if lm < -40.0:
            continue
        elif lm > 40.0:
            log_abs += lm
            tot_ang += th + math.pi
        else:
            r = math.exp(lm)
            z = complex(1.0 - r * math.cos(th), -r * math.sin(th))
            log_abs += math.log(abs(z))
            tot_ang += math.atan2(z.imag, z.real)
    return log_abs, tot_ang

def inspect_case(c, j, kap, M=512):
    roots = gaussian_a_list(c, j, kap)
    sgn_pref = 1.0 if ((j * (c + 1) + c) % 2 == 0) else -1.0
    # Find lam_L2:
    def obj(log_lam):
        lam = math.exp(log_lam)
        vals = [eval_F_log(roots, c, lam, 2.0 * math.pi * (q + 0.5) / 128)[0] for q in range(128)]
        mx = max(vals)
        return 2.0 * mx + math.log(sum(math.exp(2.0 * (v - mx)) for v in vals) / 128)
    grid = [-6.0 + 12.0 * i / 24 for i in range(25)]
    bi = min(range(len(grid)), key=lambda i: obj(grid[i]))
    lo, hi = grid[max(0, bi - 1)], grid[min(len(grid) - 1, bi + 1)]
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    for _ in range(30):
        m1, m2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
        if obj(m1) < obj(m2): hi = m2
        else: lo = m1
    lam = math.exp(0.5 * (lo + hi))

    data = [eval_F_log(roots, c, lam, 2.0 * math.pi * (k + 0.5) / M) for k in range(M)]
    mx = max(v[0] for v in data)
    # Check max |Phi| where weight >= exp(-10) of peak!
    max_phi_dev = 0.0
    s_re = 0.0
    s_abs = 0.0
    peaks = []
    for k in range(M):
        la, ph = data[k]
        # Shift phase by pi if sgn_pref < 0:
        eff_ph = (ph + (0.0 if sgn_pref > 0 else math.pi) + math.pi) % (2.0 * math.pi) - math.pi
        w = math.exp(la - mx) if (la - mx > -700) else 0.0
        s_re += w * math.cos(eff_ph)
        s_abs += w
        if la >= mx - 8.0:
            max_phi_dev = max(max_phi_dev, abs(eff_ph))
        # Local peak?
        prev_la = data[(k - 1) % M][0]
        next_la = data[(k + 1) % M][0]
        if la >= prev_la and la >= next_la and la >= mx - 1.0:
            peaks.append(round(2.0 * (k + 0.5) / M, 3))

    # Also compute Fourier coefficients A_r, B_r for r = 1..6 at lam!
    log_lam = math.log(lam)
    n_in = sum(1 for lr, _ in roots if lr + log_lam > 0)
    AB = []
    for r in range(1, 7):
        Ar = sum(math.exp(r * (lr + log_lam)) * math.cos(r * ang) for lr, ang in roots if lr + log_lam < 0 and r * (lr + log_lam) > -40)
        Br = sum(math.exp(-r * (lr + log_lam)) * math.cos(r * ang) for lr, ang in roots if lr + log_lam > 0 and -r * (lr + log_lam) > -40)
        AB.append((r, round(Ar + Br, 2), round(Ar - Br, 2)))

    print(f"c={c:4d} j={j:4d} kap={kap:.4f} c*kap={c*kap:6.2f} | lam*={lam:.4f} (n_in={n_in:4d}, diff={n_in-c:+4d}) "
          f"| rat={s_re/s_abs:+.5f} max|Phi|(w>e-8)={max_phi_dev:.4f} peaks/pi={peaks[:6]}")
    print(f"    (r, A_r+B_r, A_r-B_r) = {AB}")

cases = [
    (500, 1000, 0.01),
    (200, 201, 0.01),
    (200, 300, 0.02),
    (300, 400, 0.01),   # r0=3, D=0 (q=3)
    (300, 401, 0.01),   # r0=3, D=1 (j/c = 4/3 + 1/300)
    (300, 401, 0.005),  # r0=3, D=1, tiny c*kap = 1.5
    (300, 401, 0.002),  # r0=3, D=1, c*kap = 0.6!
    (200, 201, 0.002),  # r0=1, D=1, c*kap = 0.4!
    (200, 201, 0.0005), # r0=1, D=1, c*kap = 0.1!
    (500, 701, 0.004),  # r0=5, D=1 (7/5 + 1/500), c*kap = 2.0
    (500, 701, 0.001),  # r0=5, D=1, c*kap = 0.5
    (100, 141, 0.02),   # j/c = 1.41 ~ sqrt(2) (irrational / minor arc!)
    (100, 162, 0.02),   # j/c = 1.62 ~ phi (golden ratio!)
    (100, 162, 0.005),  # golden ratio, c*kap = 0.5!
]
for c, j, kap in cases:
    inspect_case(c, j, kap)
