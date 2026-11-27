"""
Fast dissection of the Regime 2 Gap G = { c*kappa^2 < 22, c*q*kappa > 0.8, kappa < 2.5 }!
Tests the exact Gaussian-limit model AND exact finite-n model on extreme corners of G:
  - c=500, j=1000, kap=0.01 (D=0, q=1)
  - c=200, j=201, kap=0.01  (D=1, q=200)
  - c=5000, j=5001, kap=0.0001..0.05
  - c=50, j=53, kap=0.2     (D=3, q=50)
  - c=20, j=23, kap=0.5     (moderate c, moderate kap)
  - c=10, j=13, kap=1.0     (small c, moderate kap)
  - c=7, j=14, kap=1.5      (small c, D=0)
Computes:
  1. Exact saddle point lam_* of L^2 norm (and L^1 / real-axis / log|F| saddle),
  2. How many roots lie inside |w| = lam_* vs outside,
  3. Exact ratio Re(int I) / int |I| at lam_*,
  4. Exact phase Phi(phi) and log|I(phi)| profiles around the peak(s)!
"""
import math
import cmath

def gaussian_a_list(c, j, kap, K_cut=None):
    """Exact Gaussian-limit roots a_k for even/odd c:
       a_k = (-1)^j * exp((pi/4)*c*kap - (pi*kap/c)*x^2) * exp(2*pi*i * x * j / c)
       where x = k (even c) or k + 1/2 (odd c).
    """
    alpha = math.pi * kap / c
    # Cutoff where alpha*x^2 - (pi/4)*c*kap > 40 => x > sqrt(c^2/4 + 40/alpha)
    x_max = int(math.sqrt(0.25 * c * c + 45.0 / alpha)) + 5
    if K_cut is not None:
        x_max = min(x_max, K_cut)
    shift = 0.0 if (c % 2 == 0) else 0.5
    roots = []
    sgn = 1.0 if (j % 2 == 0) else -1.0
    for k in range(-x_max, x_max + 1):
        x = k + shift
        log_r = alpha * (0.25 * c * c - x * x)
        ang = 2.0 * math.pi * x * j / c
        if sgn < 0:
            ang += math.pi
        roots.append((log_r, ang))
    return roots

def eval_F_log(roots, c, lam, phi):
    """Return (log|w^{-c} prod(1 - a_k w)|, total_phase) at w = lam * e^{i phi}."""
    log_lam = math.log(lam)
    log_abs = -c * log_lam
    tot_ang = -c * phi
    for log_r, ang in roots:
        # a_k * w has log_mod = log_r + log_lam, angle = ang + phi
        lm = log_r + log_lam
        th = ang + phi
        if lm < -40.0:
            continue
        elif lm > 40.0:
            # 1 - a_k w ~ -a_k w
            log_abs += lm
            # phase of (-e^{i th}) = th + pi
            tot_ang += th + math.pi
        else:
            r = math.exp(lm)
            z = complex(1.0 - r * math.cos(th), -r * math.sin(th))
            log_abs += math.log(abs(z))
            tot_ang += math.atan2(z.imag, z.real)
    return log_abs, tot_ang

def find_lam_L2(roots, c, j, M=256):
    """Find lam minimizing log(avg |w^{-c} F(w)|^2) via golden-section search in log(lam)."""
    def obj(log_lam):
        lam = math.exp(log_lam)
        vals = [eval_F_log(roots, c, lam, 2.0 * math.pi * (q + 0.5) / M)[0] for q in range(M)]
        mx = max(vals)
        return 2.0 * mx + math.log(sum(math.exp(2.0 * (v - mx)) for v in vals) / M)
    # Bracket around log_lam = 0:
    lo, hi = -15.0, 15.0
    # Narrow bracket roughly first:
    grid = [lo + (hi - lo) * i / 30 for i in range(31)]
    best_i = min(range(len(grid)), key=lambda i: obj(grid[i]))
    lo = grid[max(0, best_i - 1)]
    hi = grid[min(len(grid) - 1, best_i + 1)]
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    for _ in range(35):
        m1 = hi - gr * (hi - lo)
        m2 = lo + gr * (hi - lo)
        if obj(m1) < obj(m2):
            hi = m2
        else:
            lo = m1
    return math.exp(0.5 * (lo + hi))

def analyze_case(c, j, kap, M=512):
    roots = gaussian_a_list(c, j, kap)
    d = math.gcd(c, j)
    q = c // d
    # Check at lam = 1 and at lam = lam_L2:
    lam_opt = find_lam_L2(roots, c, j, M=min(M, 256))
    res = {}
    for name, lam in [("1", 1.0), ("opt", lam_opt)]:
        data = [eval_F_log(roots, c, lam, 2.0 * math.pi * (k + 0.5) / M) for k in range(M)]
        mx = max(v[0] for v in data)
        sgn_pref = 1.0 if ((j * (c + 1) + c) % 2 == 0) else -1.0  # wait: check exact sign of [w^c] prod(1 - a_k w)!
        # WAIT: we know [w^c] prod(1 - a_k w) = (-1)^c e_c(a) = (-1)^{c + jc} ... let's print raw sign too!
        # Remember a_k = (-1)^j Rz / w_k, so e_c(a) = (-1)^{jc} Rz^c e_c(w^{-1}),
        # and we want (-1)^j e_c(w^{-1}) > 0 <=> (-1)^{j + jc + c} [w^c] prod(1 - a_k w) > 0!
        s_re = 0.0
        s_abs = 0.0
        n_inside = sum(1 for lr, _ in roots if lr + math.log(lam) > 0)
        for la, ph in data:
            w = math.exp(la - mx) if (la - mx > -700) else 0.0
            s_re += w * math.cos(ph)
            s_abs += w
        rat = (sgn_pref * s_re) / s_abs if s_abs > 0 else 0.0
        # Find where max |I| occurs:
        best_idx = max(range(M), key=lambda k: data[k][0])
        best_phi = 2.0 * math.pi * (best_idx + 0.5) / M
        res[name] = (lam, rat, n_inside, best_phi)
    # Also compute log|a_0|:
    log_a0 = 0.25 * math.pi * c * kap
    print(f"c={c:4d} j={j:4d} d={d:3d} q={q:4d} kap={kap:.4f} c*kap={c*kap:6.2f} c*k^2={c*kap*kap:5.2f} |a0|={math.exp(log_a0):8.2e} | "
          f"lam=1: rat={res['1'][1]:+.4f} (in={res['1'][2]:4d}) | "
          f"lam*={res['opt'][0]:8.4f} (lam*|a0|={res['opt'][0]*math.exp(log_a0):7.3f}, in={res['opt'][2]:4d}): "
          f"rat={res['opt'][1]:+.4f} peak_phi/pi={res['opt'][3]/math.pi:.3f}")

print("=== DISSECTION OF REGIME 2 GAP G ===")
test_cases = [
    (500, 1000, 0.01),
    (200, 201, 0.01),
    (200, 201, 0.05),
    (200, 201, 0.20),
    (200, 203, 0.02),
    (200, 300, 0.02),
    (50, 51, 0.05),
    (50, 51, 0.30),
    (50, 53, 0.20),
    (50, 75, 0.20),
    (50, 100, 0.20),
    (20, 21, 0.50),
    (20, 23, 0.50),
    (20, 27, 0.50),
    (10, 11, 1.00),
    (10, 13, 1.00),
    (10, 15, 1.00),
    (7, 14, 1.50),
    (7, 10, 1.50),
    (6, 7, 1.80),
    (4, 5, 2.20),
    (3, 4, 2.20),
]
for c, j, kap in test_cases:
    analyze_case(c, j, kap)
