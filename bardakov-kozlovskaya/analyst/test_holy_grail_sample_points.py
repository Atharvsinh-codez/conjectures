"""
THE HOLY GRAIL TEST!!
For even n = 2m and EVERY c = 1, 2, ..., floor(n/3) with d = n - 3c >= 1:
Define the d + 1 EXPLICIT sample points:
  r_j = sin(j*pi/n) / sin((j+c)*pi/n),    u_j = (-1)^c * r_j^n,    j = c, c+1, ..., n - 2c.
(Note: for even c, u_j > 0; for odd c, u_j < 0; and r_{n-c-j} = 1 / r_j!)

TEST:
1. Does f_{n, c}(u_j) = Q_{n, c}(1 + u_j, u_j) STRICTLY ALTERNATE SIGN across j = c, c+1, ..., n - 2c
   for ALL even n = 4, 6, ..., 50 and ALL c = 1..floor(n/3)?!
2. Also test the SAME sample points on E_{n-c}(1, u_j) = (-1)^c [z^c] p_n(1, u_j, z)
   for j = 1, 2, ..., n - c - 1 (all n - c + 1 sample points? Wait: j = 1..n-c-1 is n-c-1 points, plus 0 and inf gives n-c+1 -> ALL n-c roots of E_{n-c}(1, u)!!)!
3. Compare f_{n, c}(u_j) against (-1)^c E_{n-c}(1, u_j) / u_j^c at EVERY sample point u_j!!
   CHECK THE RATIO  f_{n, c}(u_j) / ( (-1)^c E_{n-c}(1, u_j) / u_j^c ) !!
"""
import mpmath as mp
import sys
sys.path.insert(0, r'd:\solving-conjectures\work\bardakov-kozlovskaya\verify')
import engine

mp.mp.dps = 80

def eval_poly_mp(coeffs_low_to_high, u):
    acc = mp.mpf(0)
    for c in reversed(coeffs_low_to_high):
        acc = acc * u + c
    return acc

def eval_fnc(C, n, c, u):
    # f_{n,c}(u) = sum_{b=0}^{(n-3c)//2} c_{n-3c-2b, b, c} u^b (1+u)^{n-3c-2b}
    d = n - 3 * c
    # Evaluate via G_{n,c}(u/(1+u)^2) * (1+u)^d (avoid u=-1 division!)
    if abs(1 + u) < mp.mpf('1e-30'):
        return mp.mpf(0) if d % 2 == 1 else C[(0, d//2, c)] * (u**(d//2))
    t = u / ((1 + u)**2)
    g = [C[(d - 2*b, b, c)] for b in range(d // 2 + 1)]
    return eval_poly_mp(g, t) * ((1 + u)**d)

print("=== TEST 1: Explicit sample points u_j = (-1)^c * (sin(j*pi/n)/sin((j+c)*pi/n))^n for f_{n,c}(u) ===")
total_polys = 0
pass_polys = 0
fails = []
for n in [4, 6, 8, 10, 12, 14, 16, 20, 24, 30, 40, 50]:
    C = engine.compute_elementary_coeffs(n)
    E = engine.compute_E_table(n)
    n_ok = 0; n_tot = 0
    min_ratio_f_E = mp.inf
    max_ratio_f_E = -mp.inf
    worst_fail_n = None
    for c in range(1, n // 3 + 1):
        d = n - 3 * c
        if d <= 0: continue
        n_tot += 1; total_polys += 1
        # Sample points j = c, c+1, ..., n - 2c (total d + 1 points!)
        # Wait: at j = (n-c)/2 (when n-c is even, i.e. c is even!), r_j = 1 -> u_j = +1 (for even c)!
        # What if c is odd and j = (n-c)/2? Not an integer! Good!
        vals = []
        ratios = []
        for j in range(c, n - 2*c + 1):
            rj = mp.sin(j * mp.pi / n) / mp.sin((j + c) * mp.pi / n)
            uj = ((-1)**c) * (rj**n)
            fv = eval_fnc(C, n, c, uj)
            # Also evaluate (-1)^c * E_{n-c}(1, uj) / uj^c = [z^c] p_n(1, uj, z) / uj^c
            # Wait: [z^c] p_n(1, u, z) = (-1)^{n-c} E_{n-c}(1, u) = (-1)^c E_{n-c}(1, u)
            Ev = ((-1)**c) * eval_poly_mp(E[n - c], uj) / (uj**c)
            vals.append(fv)
            if Ev != 0:
                rat = fv / Ev
                ratios.append(rat)
                if rat < min_ratio_f_E: min_ratio_f_E = rat
                if rat > max_ratio_f_E: max_ratio_f_E = rat
        # Check strict sign alternation!
        signs = [1 if v > 0 else (-1 if v < 0 else 0) for v in vals]
        alt = all(signs[idx] != 0 and signs[idx] == -signs[idx - 1] for idx in range(1, len(signs)))
        if alt:
            n_ok += 1; pass_polys += 1
        else:
            fails.append((n, c, d, signs))
            if worst_fail_n is None:
                worst_fail_n = (c, d, signs[:10])
    print(f"n={n:2d}: PASS = {n_ok:2d}/{n_tot:2d} | f(u_j) / ([z^c]p_n/u_j^c) in [{float(min_ratio_f_E):+.4f}, {float(max_ratio_f_E):+.4f}] | first_fail={worst_fail_n}")
print(f"\nOVERALL: {pass_polys}/{total_polys} polynomials strictly alternate at explicit Maillot sample points!")
if fails:
    print("First 8 fails:", fails[:8])
