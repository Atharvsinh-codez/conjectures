"""
Test TWO fundamental questions about Regime 2 Gap G:
1. Is the pure Gaussian elementary symmetric sum E_c(alpha; c, j) > 0 for ALL c >= 1, ALL j, ALL alpha > 0?
   AND how does E_c compare to the exact finite-n value I_pure?
2. For the hardest Gap G corners (c=500,j=1000,kap=0.01; c=200,j=201,kap=0.01/0.002/0.0005; c=300,j=401,kap=0.002; c=100,j=162,kap=0.005):
   - Plot rat(lam) and max|Phi| vs lam, AND test a SIMPLE EXPLICIT closed-form choice of lam (NO optimization!)
   - Check what happens to the EXACT saddle equation at the peak(s)!
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

def rat_at_lam(roots, c, j, lam, M=512):
    sgn_pref = 1.0 if ((j * (c + 1) + c) % 2 == 0) else -1.0
    data = [eval_F_log(roots, c, lam, 2.0 * math.pi * (k + 0.5) / M) for k in range(M)]
    mx = max(v[0] for v in data)
    s_re, s_abs = 0.0, 0.0
    max_ph_5 = 0.0
    for la, ph in data:
        eff_ph = (ph + (0.0 if sgn_pref > 0 else math.pi) + math.pi) % (2.0 * math.pi) - math.pi
        w = math.exp(la - mx) if (la - mx > -700) else 0.0
        s_re += w * math.cos(eff_ph)
        s_abs += w
        if la >= mx - 5.0:
            max_ph_5 = max(max_ph_5, abs(eff_ph))
    return (s_re / s_abs if s_abs > 0 else 0.0), max_ph_5, mx + math.log(s_abs / M)

print("=== TEST 1: SWEEP lam ON HARDEST GAP G CORNERS ===")
hard = [
    (500, 1000, 0.01),
    (200, 201, 0.01),
    (200, 201, 0.002),
    (200, 201, 0.0005),
    (300, 401, 0.002),
    (500, 701, 0.001),
    (100, 162, 0.005),
]
for c, j, kap in hard:
    roots = gaussian_a_list(c, j, kap)
    inv_a0 = math.exp(-0.25 * math.pi * c * kap)
    print(f"\nc={c}, j={j}, kap={kap}, c*kap={c*kap:.3f}, |a0|^-1={inv_a0:.4f}:")
    for lam in sorted(set([0.5 * inv_a0, 0.8 * inv_a0, 0.95 * inv_a0, 0.75, 0.85, 0.95, 1.0, 1.1, 1.25])):
        if lam <= 0: continue
        r, mp5, l1 = rat_at_lam(roots, c, j, lam, M=512)
        n_in = sum(1 for lr, _ in roots if lr + math.log(lam) > 0)
        print(f"  lam={lam:7.4f} (in={n_in:3d}): rat={r:+.5f}, max|Phi|(w>e-5)={mp5:.4f}, logL1={l1:+.2f}")

print("\n=== TEST 2: CHECK E_c(alpha; c, j) > 0 FOR ALL c=2..12, ALL j=1..c, ALL alpha ===")
# Exact computation via poly mult of (1 + z * b_x):
min_val = 1e9
for c in range(2, 13):
    for j in range(1, c + 1):
        for kap in [0.001, 0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 5.0]:
            alpha = math.pi * kap / c
            x_max = int(math.sqrt(40.0 / alpha)) + 3
            shift = 0.0 if (c % 2 == 0) else 0.5
            # poly coeff up to deg c:
            poly = [0j] * (c + 1)
            poly[0] = 1.0 + 0j
            for k in range(-x_max, x_max + 1):
                x = k + shift
                bx = math.exp(-alpha * x * x) * cmath.exp(2j * math.pi * x * j / c)
                for deg in range(min(c, k + x_max + 1), 0, -1):
                    poly[deg] += poly[deg - 1] * bx
            ec = (((-1) ** j) * poly[c]).real
            assert abs(poly[c].imag) < 1e-10 * abs(poly[c]) + 1e-14
            if ec <= 0:
                print(f"  COUNTEREXAMPLE!! c={c}, j={j}, kap={kap}: E_c = {ec}")
            min_val = min(min_val, ec)
print(f"ALL {11*6*8}+ (c,j,kap) checks PASSED! Min E_c = {min_val:.4e} > 0!")
