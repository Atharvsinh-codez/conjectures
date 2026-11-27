"""
Systematic scan of the Regime 2 Gap G = { kappa < 2.5, c*kappa^2 < 22, c*q*kappa > 0.8 }:
1. Test the exact boundary of where |Phi(phi; lam=1)| < 1.45 ALREADY holds (requiring ZERO rescaling!).
2. For ALL cases where |Phi(phi; lam=1)| >= 1.45, check:
   - Does there ALWAYS exist a unique small r0 >= 1 with D = r0*j - c*p0 small?
   - What is max(r0*kappa) among ALL failures of |Phi@1| < 1.45?
   - What is min(c / (r0 * D^2)) among ALL D != 0 failures?
   - What is min(c / r0) among ALL D = 0 failures?
3. Test an EXPLICIT 1-line closed-form choice of lam (and explicit Hayman/Bernoulli lower bound!) on ALL failures!
"""
import math
import cmath

def gaussian_a_list(c, j, kap):
    alpha = math.pi * kap / c
    x_max = int(math.sqrt(0.25 * c * c + 40.0 / alpha)) + 5
    shift = 0.0 if (c % 2 == 0) else 0.5
    roots = []
    sgn = 1.0 if (j % 2 == 0) else -1.0
    for k in range(-x_max, x_max + 1):
        x = k + shift
        log_r = alpha * (0.25 * c * c - x * x)
        ang = 2.0 * math.pi * x * j / c + (0.0 if sgn > 0 else math.pi)
        roots.append((log_r, ang))
    return roots

def max_phi_at_1(c, j, kap, M=256):
    """Compute exact max |Phi(phi; lam=1)| over phi in (0, 2*pi) (excluding the boundary zero at phi=0)."""
    roots = gaussian_a_list(c, j, kap)
    # At lam = 1, factor out the 2 boundary roots at x = +-c/2 (which have a = +1, giving (1 - e^{i phi})^2 = -e^{i phi} * 4 sin^2(phi/2))!
    # All remaining roots have |a| != 1: small (|a| > 1, count c-1) and large (|a| < 1).
    # For each small root, write 1 - a e^{i phi} = (-a e^{i phi}) * (1 - a^{-1} e^{-i phi}).
    # Total phase of (-1)^{j(c+1)+c} * e^{-i c phi} * prod(1 - a e^{i phi}) / (4 sin^2(phi/2)) is:
    #   Phi(phi) = sum_{small} arg(1 - a^{-1} e^{-i phi}) + sum_{large} arg(1 - a e^{i phi})!
    small = []
    large = []
    for lr, ang in roots:
        if abs(lr) < 1e-12:
            continue
        elif lr > 0:
            if lr < 40.0:
                small.append((math.exp(-lr), -ang))  # a^{-1} = exp(-lr) * exp(-i ang)
        else:
            if lr > -40.0:
                large.append((math.exp(lr), ang))
    mx_phi = 0.0
    for q in range(M):
        phi = 2.0 * math.pi * (q + 0.5) / M
        ph = 0.0
        for r, a in small:
            th = a - phi
            ph += math.atan2(-r * math.sin(th), 1.0 - r * math.cos(th))
        for r, a in large:
            th = a + phi
            ph += math.atan2(-r * math.sin(th), 1.0 - r * math.cos(th))
        if abs(ph) > mx_phi:
            mx_phi = abs(ph)
    return mx_phi

def best_dirichlet(c, j, max_r0=20):
    """Find r0 in 1..max_r0 minimizing r0 * D^2 (or |D|)."""
    best = None
    for r0 in range(1, min(max_r0, c) + 1):
        p0 = round(r0 * j / c)
        if math.gcd(r0, p0) != 1:
            continue
        D = r0 * j - c * p0
        score = (0, r0) if D == 0 else (1, r0 * D * D)
        if best is None or score < best[0]:
            best = (score, r0, p0, D)
    return best[1], best[2], best[3]

print("=== SCANNING WHERE |Phi@1| >= 1.40 IN GAP G ===")
failures = []
for c in [3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 50, 80, 120, 200]:
    # Test representative j values: multiples, near-multiples, half-integers, thirds, golden ratio
    js = sorted(set([
        c + 1, c + 2, c + 3, 2 * c,
        int(1.5 * c), int(1.5 * c) + 1,
        int(4 * c // 3), int(4 * c // 3) + 1,
        int(1.6180339 * c),
    ]))
    for j in js:
        if j <= c: continue
        d = math.gcd(c, j)
        q = c // d
        for kap in [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 1.8]:
            if c * kap * kap >= 22.0: continue
            if c * q * kap <= 0.8: continue
            mp = max_phi_at_1(c, j, kap, M=192)
            if mp >= 1.40:
                r0, p0, D = best_dirichlet(c, j, max_r0=20)
                failures.append((mp, c, j, kap, r0, D, r0 * kap, c / (r0 * D * D) if D != 0 else c / r0))

failures.sort(reverse=True)
print(f"Total |Phi@1| >= 1.40 hits: {len(failures)}")
if failures:
    d0 = [f for f in failures if f[5] == 0]
    dn = [f for f in failures if f[5] != 0]
    print(f"  D == 0 hits: {len(d0)}, max(r0*kap) = {max(f[6] for f in d0) if d0 else 0:.4f}, min(c/r0) = {min(f[7] for f in d0) if d0 else 0:.2f}")
    print(f"  D != 0 hits: {len(dn)}, max(r0*kap) = {max(f[6] for f in dn) if dn else 0:.4f}, min(c/(r0*D^2)) = {min(f[7] for f in dn) if dn else 0:.2f}")
    print("  Top 15 worst |Phi@1| hits:")
    for f in failures[:15]:
        print(f"    |Phi@1|={f[0]:6.3f}  c={f[1]:3d} j={f[2]:3d} kap={f[3]:.3f}  r0={f[4]:2d} D={f[5]:+3d}  r0*kap={f[6]:.4f}  c/(r0*D^2 or r0)={f[7]:6.2f}")
    print("  Worst D != 0 with smallest c/(r0*D^2):")
    for f in sorted(dn, key=lambda x: x[7])[:10]:
        print(f"    |Phi@1|={f[0]:6.3f}  c={f[1]:3d} j={f[2]:3d} kap={f[3]:.3f}  r0={f[4]:2d} D={f[5]:+3d}  r0*kap={f[6]:.4f}  c/(r0*D^2)={f[7]:6.2f}")
    print("  Worst D == 0 with smallest c/r0:")
    for f in sorted(d0, key=lambda x: x[7])[:10]:
        print(f"    |Phi@1|={f[0]:6.3f}  c={f[1]:3d} j={f[2]:3d} kap={f[3]:.3f}  r0={f[4]:2d} D={f[5]:+3d}  r0*kap={f[6]:.4f}  c/r0={f[7]:6.2f}")
