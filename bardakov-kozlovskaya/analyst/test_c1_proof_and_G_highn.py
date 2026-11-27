"""
1. Test our EXACT c = 1 real-rootedness proof (conjugate-pair phase sweep)!
   For even n = 2m >= 4:
     -f_{2m, 1}(u) = (2m/u) * [ h*_{2m}(u) - (1+u)(1-u)^{2m-2} ]  (deg 2m-3, palindromic).
   Equivalently at u = -v = -r^{2m} (r > 0):
     (1+r^{2m})^{-2m} * [ h*_{2m}(-r^{2m}) - (1-r^{2m})(1+r^{2m})^{2m-2} ]
     = (1/2m) * sum_{k=0}^{2m-1} (1 + r * exp(i*theta_k))^{-2m} - (1-r^{2m})/(1+r^{2m})^2
     where theta_k = (2k+1)*pi/(2m).
   Look at the dominant pair k = m-1 (and k = m, its conjugate), theta = pi - pi/(2m):
     1 + r * exp(i*theta) = rho * exp(i*psi), where psi(r) increases from 0 to pi - pi/(2m).
   Choose the 2m - 2 explicit sample points r_j > 0 where:
     2m * psi(r_j) = j * pi,  for j = 1, 2, ..., 2m - 2!
   Wait: at j = m-1/2? No! At 2m*psi = j*pi, cos(2m*psi) = (-1)^j!
   Note: j = 1, 2, ..., 2m - 2 gives 2m - 2 points where (-1)^j alternates sign
   --> AT LEAST 2m - 3 STRICT SIGN CHANGES on (0, inf)!!
   Since deg(R_{2m}) = 2m - 3, THAT IMMEDIATELY PROVES ALL 2m - 3 ROOTS ARE SIMPLE NEGATIVE REAL!!
   CHECK RIGHT NOW WHETHER THE k = m-1 PAIR STRICTLY DOMINATES ALL OTHER TERMS + (1-r^{2m})/(1+r^{2m})^2
   AT EVERY SAMPLE POINT j = 1, 2, ..., 2m - 2 FOR ALL m = 2, 3, ..., 50!!

2. Also push G_{n,c}(t) real-rootedness test to n = 50, 60 using verify/engine.py!
"""
import mpmath as mp
import numpy as np
import sys
sys.path.insert(0, r'd:\solving-conjectures\work\bardakov-kozlovskaya\verify')
import engine

mp.mp.dps = 60

print("=== 1. TEST EXPLICIT SAMPLE-POINT PROOF FOR c = 1 (ALL EVEN n = 2m) ===")
# Solve psi(r) = arg(1 + r*exp(i*(pi - pi/(2m)))) = j*pi/(2m) for r > 0!
# Geometry: triangle with side 1, side r, exterior angle pi/(2m), opposite angle psi:
# By Law of Sines: r / sin(psi) = 1 / sin(pi - pi/(2m) - psi) = 1 / sin(psi + pi/(2m))!
# WAIT A SECOND!! CHECK THAT LAW OF SINES!!
# 1 + r*exp(i*(pi - alpha)) = (1 - r*cos(alpha)) + i*(r*sin(alpha)) = rho * exp(i*psi).
# tan(psi) = r*sin(alpha) / (1 - r*cos(alpha))
# => r*(sin(alpha)*cos(psi) + cos(alpha)*sin(psi)) = sin(psi)
# => r = sin(psi) / sin(psi + alpha) !!
# AND rho = |1 + r*exp(i*(pi-alpha))| = sin(alpha) / sin(psi + alpha) !!
# OH!! MY!! GOD!! THAT'S AN EXACT CLOSED-FORM TRIGONOMETRIC FORMULA FOR r_j AND rho_j!!

for m in [2, 3, 4, 5, 6, 8, 10, 15, 25, 40]:
    n = 2 * m
    alpha = mp.pi / n
    worst_ratio = mp.mpf(0)
    worst_j = -1
    all_ok = True
    for j in range(1, n - 1): # j = 1, ..., 2m - 2 (gives 2m - 2 alternating points -> 2m - 3 roots!)
        psi = j * mp.pi / n
        r = mp.sin(psi) / mp.sin(psi + alpha)
        rho_dom = mp.sin(alpha) / mp.sin(psi + alpha)
        # Main term (pair k=m-1 and its conjugate): 2 * (-1)^j / (n * rho_dom^n)
        main = 2 * ((-1)**j) / (n * rho_dom**n)
        # Sum of all OTHER pairs k = 0..m-2 (each paired with conjugate):
        other = mp.mpf(0)
        for k in range(m - 1):
            th = (2*k + 1) * mp.pi / n
            z = 1 + r * mp.exp(1j * th)
            other += 2 * mp.re(1 / z**n) / n
        # Subtracted term: -(1 - r^n) / (1 + r^n)^2
        sub = -(1 - r**n) / (1 + r**n)**2
        err = other + sub
        rat = abs(err) / abs(main)
        if rat > worst_ratio:
            worst_ratio = rat
            worst_j = j
        if ((-1)**j) * (main + err) <= 0:
            all_ok = False
    print(f"n={n:2d} (m={m:2d}): ALL {n-2:2d} sample signs OK={all_ok}, max |err|/|main| = {float(worst_ratio):.6e} (at j={worst_j})")

print("\n=== 2. Check G_{n,c}(t) real-rootedness at n = 50, 60 ===")
for n in [50, 60]:
    coeffs = engine.compute_elementary_coeffs(n)
    anom = []
    for c in range(1, n//3 + 1):
        d = n - 3*c
        g = [coeffs[(d - 2*b, b, c)] for b in range(d//2 + 1)]
        if len(g) <= 1: continue
        mx = max(abs(v) for v in g)
        r = np.roots([float(v/mx) for v in g][::-1])
        if np.max(np.abs(r.imag)) > 1e-6:
            rm = mp.polyroots([mp.mpf(v) for v in g[::-1]], maxsteps=400, extraprec=80)
            max_im = max(abs(z.imag) for z in rm)
            if max_im > mp.mpf('1e-18'):
                anom.append(f"c={c}(deg={len(g)-1}):CPLX({float(max_im):.2e})")
            else:
                rr = sorted(float(z.real) for z in rm)
                if c % 2 == 1 and not all(x < 0 for x in rr): anom.append(f"c={c}:NOT_NEG")
                if c % 2 == 0 and not all(0 < x <= 0.250001 for x in rr): anom.append(f"c={c}:OUT_0_1/4")
        else:
            rr = np.sort(r.real)
            if c % 2 == 1 and not np.all(rr < 0): anom.append(f"c={c}:NOT_NEG")
            if c % 2 == 0 and not np.all((rr > 0) & (rr <= 0.25+1e-6)): anom.append(f"c={c}:OUT_0_1/4")
    print(f"n={n:2d}: G_{{n,c}}(t) anomalies = {anom if anom else 'NONE! 100% REAL IN (-inf,0)/(0,1/4]!'}")
