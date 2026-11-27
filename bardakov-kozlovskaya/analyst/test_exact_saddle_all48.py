"""
1. Test RANDOM IRRATIONAL ANGLES on the pure monomial integrand Q(r2^n e^{i phi2}, r3^n e^{i psi})
   at n = 6, 12, 24 to see if Re(sgn * Q) is literally >= 0 everywhere or has O(1/n) negative dips!
2. Test the EXACT saddle point (R2*, R3*) (found by minimizing log(P_even(1,-R2,R3)/(R2^b R3^c))
   or log(-P_odd(1,R2,R3)/(R2^b R3^c))) for ALL 48 c>=1 coefficients of n=24, especially (1,10,1) & (2,8,2),
   AND compare (R2*, R3*) with the closed-form Maillot + Vandermonde radii!
"""
import mpmath as mp
import numpy as np
import sys
sys.path.insert(0, r'd:\solving-conjectures\work\bardakov-kozlovskaya\verify')
import engine

mp.mp.dps = 50

def eval_Q_monomial(n, k2, k3, y, z):
    beta = y**(mp.mpf(1)/n)
    prod = mp.mpc(1)
    for k in range(n):
        prod *= (z - (1 + mp.exp(2*mp.pi*1j*k/n)*beta)**n)
    return prod / (y**k2 * z**k3)

print("=== PART 1: Random irrational-angle stress test of pure monomial integrand Q(y, z) ===")
rng = np.random.default_rng(2026)
for n in [6, 12, 24]:
    E = engine.compute_E_table(n)
    C = engine.compute_elementary_coeffs(n)
    min_re_all = 1e9
    worst_pt = None
    for (a, b, c) in sorted(C.keys()):
        if c == 0: continue
        k1, k2, k3 = a+b+c, b+c, c
        M_val = ((-1)**(n - k3)) * E[n - k3][k2]
        sgn = 1 if M_val > 0 else -1
        r2 = mp.sin(mp.pi * k2 / n) / mp.sin(mp.pi * k1 / n)
        r3 = mp.sin(mp.pi * k3 / n) / mp.sin(mp.pi * k1 / n)
        Ry, Rz = r2**n, r3**n
        # Test 25 random angles
        for _ in range(25):
            ph1 = mp.mpf(float(rng.uniform(-np.pi, np.pi)))
            ph2 = mp.mpf(float(rng.uniform(-np.pi, np.pi)))
            v = sgn * eval_Q_monomial(n, k2, k3, Ry*mp.exp(1j*ph1), Rz*mp.exp(1j*ph2)) / abs(M_val)
            rv = float(v.real)
            if rv < min_re_all:
                min_re_all = rv
                worst_pt = (a, b, c, float(ph1), float(ph2), rv, float(v.imag))
    print(f"n={n:2d}: global min Re(sgn*Q / |M|) over random irrational angles = {min_re_all:+.6f} at {worst_pt}")

print("\n=== PART 2: Exact saddle point (R2*, R3*) for ALL tricky cases at n=24 ===")
# Use exact polynomial P_24(1, e2, e3) to find exact saddle point (R2*, R3*) via Newton's method
# on log-convex P_even(1, -R2, R3)/(R2^b R3^c) [even c] or -P_odd(1, R2, R3)/(R2^b R3^c) [odd c]!
n = 24
C = engine.compute_elementary_coeffs(n)

def eval_P_poly(C, e2, e3):
    # Evaluate P_n(1, e2, e3) directly from dict C (100x faster than polyroots!)
    acc = mp.mpc(0)
    # Precompute powers
    p2 = [mp.mpc(1)]
    for _ in range(n//2): p2.append(p2[-1] * e2)
    p3 = [mp.mpc(1)]
    for _ in range(n//3): p3.append(p3[-1] * e3)
    for (a, b, c), v in C.items():
        acc += v * p2[b] * p3[c]
    return acc

def find_exact_saddle(C, n, a, b, c):
    # Minimize f(u, v) = log( sum |c_{a',b',c'}| exp(b'*u + c'*v) ) - b*u - c*v over parity of c
    parity = c % 2
    terms = [(b1, c1, abs(v)) for (a1, b1, c1), v in C.items() if c1 % 2 == parity]
    # Start near Maillot guess
    k1, k2, k3 = a+b+c, b+c, c
    r2 = float(mp.sin(mp.pi*k2/n)/mp.sin(mp.pi*k1/n)) if a > 0 else 1.5
    r3 = float(mp.sin(mp.pi*k3/n)/mp.sin(mp.pi*k1/n)) if a > 0 else float(mp.sin(mp.pi*k3/n))
    u = n * np.log(max(r2, 1e-6))
    v = n * np.log(max(r2*r3, 1e-12))
    # Handle boundary cases b=0 or a=0 cleanly (saddle goes to -inf/+inf; clamp or use finite step)
    if b == 0: u = -60.0
    if a == 0: u = +20.0
    for it in range(40):
        # log-sum-exp for stability
        exps = [b1*u + c1*v + float(mp.log(val)) for b1, c1, val in terms]
        mx = max(exps)
        ws = [np.exp(e - mx) for e in exps]
        Z = sum(ws)
        eb = sum(w*b1 for w, (b1, c1, _) in zip(ws, terms)) / Z
        ec = sum(w*c1 for w, (b1, c1, _) in zip(ws, terms)) / Z
        vbb = sum(w*(b1-eb)**2 for w, (b1, c1, _) in zip(ws, terms)) / Z
        vcc = sum(w*(c1-ec)**2 for w, (b1, c1, _) in zip(ws, terms)) / Z
        vbc = sum(w*(b1-eb)*(c1-ec) for w, (b1, c1, _) in zip(ws, terms)) / Z
        gb, gc = eb - b, ec - c
        if abs(gb) < 1e-9 and abs(gc) < 1e-9:
            break
        det = vbb*vcc - vbc*vbc
        if det < 1e-14:
            # 1D step in v only
            if vcc > 1e-14: v -= np.clip(gc / vcc, -5, 5)
            break
        du = (vcc*gb - vbc*gc) / det
        dv = (-vbc*gb + vbb*gc) / det
        step = min(1.0, 5.0 / max(abs(du), abs(dv), 1e-12))
        u -= step * du
        v -= step * dv
        u = np.clip(u, -100, 100)
        v = np.clip(v, -200, 200)
    return mp.exp(mp.mpf(u)), mp.exp(mp.mpf(v))

# Check ALL 48 (c>=1) coefficients of n=24 at their exact saddle (R2*, R3*)!
worst_rat = 10.0
worst_minre = 10.0
worst_case = None
print("Testing ALL 48 (c>=1) coeffs of n=24 at exact saddle (R2*, R3*) with 16x16 torus grid (P_even or P_odd):")
for (a, b, c), exact_c in sorted(C.items()):
    if c == 0: continue
    sgn = 1 if exact_c > 0 else -1
    R2, R3 = find_exact_saddle(C, n, a, b, c)
    # Evaluate P_n^{parity}(1, R2 e^{i phi2}, R3 e^{i phi3}) * e^{-i(b phi2 + c phi3)}
    # Note: P_even/P_odd = (P(e2,e3) +- P(e2,-e3))/2!
    scale = (R2**b) * (R3**c) * abs(exact_c)
    M = 12
    vals = []
    for j2 in range(M):
        ph2 = 2*mp.pi*j2/M
        e2 = R2 * mp.exp(1j*ph2)
        for j3 in range(M):
            ph3 = 2*mp.pi*j3/M
            e3 = R3 * mp.exp(1j*ph3)
            vp = eval_P_poly(C, e2, e3)
            vm = eval_P_poly(C, e2, -e3)
            p_par = (vp + vm)/2 if c % 2 == 0 else (vp - vm)/2
            v = sgn * p_par * mp.exp(-1j*(b*ph2 + c*ph3)) / scale
            vals.append(v)
    avg = sum(vals)/len(vals)
    abs_avg = sum(abs(v) for v in vals)/len(vals)
    rat = float(abs(avg)/abs_avg)
    min_re = min(float(v.real) for v in vals)
    if rat < worst_rat:
        worst_rat = rat
        worst_case = (a, b, c, rat, min_re, float(mp.log(R2)/n), float(mp.log(R3)/n))
    if min_re < worst_minre:
        worst_minre = min_re
    if (a, b, c) in [(1,10,1), (2,8,2), (1,7,3), (4,4,4), (0,9,2), (21,0,1), (0,0,8), (1,1,7)]:
        print(f"  ({a:2d},{b:2d},{c}): logR2/n={float(mp.log(R2)/n):+.4f}, logR3/n={float(mp.log(R3)/n):+.4f} | Re in [{min_re:+.4f}, {max(float(v.real) for v in vals):+.4f}] | |avg|/avg|F|={rat:.4f}")
print(f"OVERALL n=24 (ALL 48 c>=1): worst |avg|/avg|F| = {worst_rat:.4f} at {worst_case}, global min Re = {worst_minre:+.4f}")
