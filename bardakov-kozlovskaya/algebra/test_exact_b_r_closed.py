"""
Test subtracting the 1-line exact closed-form geometric sum:
    b_r^{main} = 4 * r * (kap / c) * Re( z * (1 + z) / (1 - z)^3 ),   z = exp(-r*kap + 2*pi*i * r * j / c)!
1. Check how close b_r^{main} is to exact finite-n b_r across G*!
2. Check the EXACT maximum of |b_r^{main}| / r as a function of (c, j, kap, r), AND check
   sum_{r>=1} |b_r^{main}| / r on ALL Minor Arcs outside:
     - G2 (D=0): d >= 5, d*q^4*kap^2 <= 25 (or whatever exact cutoff!),
     - G3 (D!=0): c/(r0*D^2) >= 15, r0*kap <= 0.25!
"""
import math
import cmath

def br_main(c, j, kap, r):
    z = math.exp(-r * kap) * cmath.exp(2j * math.pi * r * j / c)
    return 4.0 * r * (kap / c) * (z * (1.0 + z) / ((1.0 - z) ** 3)).real

def sum_br_main(c, j, kap, R_max=300):
    tot = 0.0
    max_term = (0.0, 1)
    for r in range(1, R_max + 1):
        v = abs(br_main(c, j, kap, r)) / r
        tot += v
        if v > max_term[0]:
            max_term = (v, r)
    # Tail r > R_max bounded by sum 8/(r^3 c kap^2):
    tot += 8.0 / (c * kap * kap * 2.0 * R_max * R_max)
    return tot, max_term[1], max_term[0]

# Check exact analytic max of |4*r*(kap/c)*Re(z(1+z)/(1-z)^3)| / r over ALL r*kap > 0 for fixed theta = 2*pi*r*j/c!
print("=== 1. EXACT ANALYTIC BOUND ON |b_r^{main}| / r AS A FUNCTION OF theta = 2*pi*r*j/c ===")
for deg in [1, 2, 5, 10, 20, 45, 90, 180]:
    th = deg * math.pi / 180.0
    # Maximize 4*gamma * |Re(z(1+z)/(1-z)^3)| * (2*sin(th/2))^2 over gamma > 0:
    # Note: for small th, max over gamma of 4*gamma*|Re(...)| * th^2 should be ~ 8 * 0.27 = 2.16 -> Wait: 0.055 * (2*pi)^2 = 2.17!
    best = 0.0
    for k in range(1, 500):
        gam = th * (k / 100.0)
        z = math.exp(-gam) * cmath.exp(1j * th)
        val = 4.0 * gam * abs((z * (1.0 + z) / ((1.0 - z) ** 3)).real) * (4.0 * (math.sin(th / 2.0) ** 2))
        if val > best: best = val
    print(f"  theta={deg:3d} deg: max_gam [4*gam*|Re z(1+z)/(1-z)^3| * 4 sin^2(th/2)] = {best:.5f}  (-> coeff of c/(r*D^2) is {best/(4*math.pi*math.pi):.5f}!)")

print("\n=== 2. SCAN sum |b_r^{main}| / r ACROSS ALL MINOR ARCS OUTSIDE G2 (d>=5) AND G3 (c/(r0*D^2)>=15, r0*kap<=0.25) ===")
worst_minor = (0.0, None)
hits_over_135 = []
for c in [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 50, 80, 120, 200]:
    for j in range(c + 1, min(2 * c + 1, c + 45)):
        d = math.gcd(c, j); q = c // d
        for kap in [0.005, 0.01, 0.02, 0.04, 0.08, 0.15, 0.30, 0.60, 1.20, 2.00]:
            if c * kap * kap >= 42.0 or c * q * kap <= 0.85: continue
            # Check if inside G2 (D=0): d >= 5 and d*q^4*kap^2 <= 25:
            if d >= 5 and d * ((q * q * kap) ** 2) <= 25.0:
                continue
            # Check if inside G3 (D!=0): exists r0 >= 1 coprime to p0 with D != 0, c/(r0*D^2) >= 15, r0*kap <= 0.25:
            in_G3 = False
            for r0 in range(1, min(c, int(0.25 / kap) + 1) + 1):
                if r0 * kap > 0.25 + 1e-12: break
                p0 = round(r0 * j / c)
                D = r0 * j - c * p0
                if D != 0 and math.gcd(r0, p0) == 1 and c / (r0 * D * D) >= 15.0:
                    in_G3 = True
                    break
            if in_G3:
                continue
            tot, r_dom, t_dom = sum_br_main(c, j, kap, R_max=200)
            if tot > worst_minor[0]:
                worst_minor = (tot, (c, j, d, q, kap, r_dom, round(t_dom, 3)))
            if tot > 1.35:
                hits_over_135.append((tot, c, j, d, q, kap, r_dom, round(t_dom, 3)))

print(f"  WORST Minor-Arc sum |b_r^main|/r = {worst_minor[0]:.4f} at {worst_minor[1]}")
print(f"  Total Minor-Arc hits > 1.35: {len(hits_over_135)}")
for h in sorted(hits_over_135, reverse=True)[:15]:
    print(f"    sum={h[0]:.4f}: c={h[1]:3d} j={h[2]:3d} d={h[3]:2d} q={h[4]:3d} kap={h[5]:.3f} r_dom={h[6]:3d} t_dom={h[7]:.3f}")
