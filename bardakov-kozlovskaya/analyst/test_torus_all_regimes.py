"""
TEST THE NONNEGATIVE-SQUARE TORUS MIRACLE ACROSS EVERY SINGLE REGIME OF THE SIMPLEX!
For n = 24 and n = 36, test ALL boundary & corner & interior types:
  - Deep interior: (4,4,4) at n=24; (6,6,6) at n=36
  - Near c=1: (n-3, 0, 1), (n-5, 1, 1), (1, (n-4)//2, 1)
  - Near c=2: (n-6, 0, 2), (n-8, 1, 2), (0, (n-6)//2, 2)
  - Near b=0: (n-3c, 0, c) for c = 1, 2, 3, n//6, n//3
  - Near b=1: (n-2-3c, 1, c)
  - Near a=0: (0, b, c) with 2b+3c=n (c even!)
  - Near a=1: (1, b, c) with 2b+3c=n-1 (c odd!)
  - Near equilateral corner: (0, 0, n//3), (1, 1, n//3 - 1), (0, 3, n//3 - 2), (2, 2, n//3 - 2)

For each (a, b, c), evaluate
  F(phi2, phi3) = sgn * P_n(1, R2 e^{i*phi2}, R3 e^{i*phi3}) * e^{-i*b*phi2 - i*c*phi3} / (R2^b R3^c |c_{a,b,c}|)
on a 12 x 12 grid (actually check BOTH the exact saddle radius R2, R3 from Maillot AND
optimized/adjusted R2, R3 when a=0 or b=0 where r2->1 or r2=r3!).
Report: avg(Re F) [should be +1.000], min Re(F), max Re(F), and ratio |avg|/avg|F|!
"""
import mpmath as mp
import sys
sys.path.insert(0, r'd:\solving-conjectures\work\bardakov-kozlovskaya\verify')
import engine

mp.mp.dps = 50

def eval_Pn_from_e(n, e2, e3):
    roots = mp.polyroots([mp.mpc(1), mp.mpc(-1), mp.mpc(e2), mp.mpc(-e3)], maxsteps=200, extraprec=50)
    x, y, z = roots
    alpha = x**(mp.mpf(1)/n)
    beta = y**(mp.mpf(1)/n)
    prod = mp.mpc(1)
    for k in range(n):
        prod *= (z - (alpha + mp.exp(2*mp.pi*1j*k/n)*beta)**n)
    return prod

def test_point(n, a, b, c, coeffs, M=12, label=""):
    exact_c = coeffs[(a, b, c)]
    sgn = 1 if exact_c > 0 else -1
    k1, k2, k3 = a + b + c, b + c, c
    # Saddle radii from Maillot:
    # Note: for a=0, k1=k2 => r2=1 => R2=1! Wait: at e1=1, e2=R2*e^{i*phi2}, if R2=1 then x, y are O(1)!
    # Better: use exact saddle point of log P_n (or for a=0 just test R2=1 or R2=2, and note for a=0 we can just take R2 -> inf!)
    # Let's first test raw Maillot radii r2 = sin(pi*k2/n)/sin(pi*k1/n), r3 = sin(pi*k3/n)/sin(pi*k1/n)!
    # If a=0, k1=k2, r2=1; if we nudge slightly or use r2=1, let's see what happens!
    r2 = mp.sin(mp.pi * k2 / n) / mp.sin(mp.pi * k1 / n)
    r3 = mp.sin(mp.pi * k3 / n) / mp.sin(mp.pi * k1 / n)
    # Handle a=0 (r2=1): at e1=1, the saddle in e2 is at infinity! So scale R2 = 10, R3 = 10*r3^n? No,
    # if a=0, e1^0 means we can set e1 = 0 directly (P_n(0, e2, e3)!)!
    if a == 0:
        # Test exact 1D circle for P_n(0, 1, R3 e^{i*phi3}) (since a=0 => 2b+3c=n, c=2k)!
        # Here e1=0, e2=1, e3=R3*e^{i*phi3} -> roots of T^3 + T - e3 = 0 -> x ~ i, y ~ -i, z ~ e3!
        # Wait: better test R2 = 100.0 (large e2/e1^2) OR direct saddle for (a=0)!
        # Let's test direct Maillot with e1=0: x,y,z roots of T^3 + e2*T - e3 = 0 with e2 = 1, e3 = r3^n?
        # Wait: for a=0, k1 = k2 = (n-c)/2, k3 = c, so r1 = r2 = 1, r3 = sin(pi*c/n)/sin(pi*(n-c)/(2n)) = 2*sin(pi*c/(2n))!
        R3 = r3**n
        vals = []
        for j3 in range(2*M):
            phi3 = -mp.pi + 2*mp.pi*j3/(2*M)
            e3 = R3 * mp.exp(1j * phi3)
            roots = mp.polyroots([mp.mpc(1), mp.mpc(0), mp.mpc(1), mp.mpc(-e3)], maxsteps=200, extraprec=50)
            x, y, z = roots
            alpha, beta = x**(mp.mpf(1)/n), y**(mp.mpf(1)/n)
            prod = mp.mpc(1)
            for k in range(n):
                prod *= (z - (alpha + mp.exp(2*mp.pi*1j*k/n)*beta)**n)
            v = sgn * prod * mp.exp(-1j * c * phi3) / (R3**c * abs(exact_c))
            vals.append(v)
        avg = sum(vals)/len(vals)
        abs_avg = sum(abs(v) for v in vals)/len(vals)
        re_vals = [float(v.real) for v in vals]
        print(f"[a=0 1D] {label:16s} n={n} ({a:2d},{b:2d},{c:2d}): avg={float(avg.real):+.5f} | Re in [{min(re_vals):+.3f}, {max(re_vals):+.3f}] | |avg|/avg|F|={float(abs(avg)/abs_avg):.4f}")
        return

    R2 = r2**n
    R3 = (r2 * r3)**n
    scale = (R2**b) * (R3**c) * abs(exact_c)
    vals = []
    for j2 in range(M):
        phi2 = -mp.pi + 2*mp.pi*j2/M
        for j3 in range(M):
            phi3 = -mp.pi + 2*mp.pi*j3/M
            e2 = R2 * mp.exp(1j * phi2)
            e3 = R3 * mp.exp(1j * phi3)
            v = sgn * eval_Pn_from_e(n, e2, e3) * mp.exp(-1j * (b*phi2 + c*phi3)) / scale
            vals.append(v)
    avg = sum(vals) / len(vals)
    abs_avg = sum(abs(v) for v in vals) / len(vals)
    re_vals = [float(v.real) for v in vals]
    print(f"         {label:16s} n={n} ({a:2d},{b:2d},{c:2d}): avg={float(avg.real):+.5f} | Re in [{min(re_vals):+.3f}, {max(re_vals):+.3f}] | |avg|/avg|F|={float(abs(avg)/abs_avg):.4f}")

n = 24
coeffs = engine.compute_elementary_coeffs(n)
tests = [
    ("deep interior",  4, 4, 4),
    ("c=1, b=0",      21, 0, 1),
    ("c=1, b=1",      19, 1, 1),
    ("c=1, b=mid",    11, 5, 1),
    ("c=1, a=1",       1, 10, 1),
    ("c=2, b=0",      18, 0, 2),
    ("c=2, b=1",      16, 1, 2),
    ("c=2, b=mid",    10, 4, 2),
    ("c=2, a=0",       0, 9, 2),
    ("c=3, a=1",       1, 7, 3),
    ("b=0, c=mid",    12, 0, 4),
    ("b=1, c=mid",    10, 1, 4),
    ("a=0, c=mid",     0, 6, 4),
    ("a=1, c=mid",     1, 4, 5),
    ("near corner 1",  1, 1, 7),
    ("near corner 2",  2, 2, 6),
    ("near corner 3",  3, 0, 7),
    ("corner a=0,b=0", 0, 0, 8),
]
for lbl, a, b, c in tests:
    assert a + 2*b + 3*c == n
    test_point(n, a, b, c, coeffs, M=10, label=lbl)
