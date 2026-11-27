"""
Verify our exact Maillot / hyperbolic-volume asymptotic formula for |A_{k1,k2,k3}|
and inspect the EXACT torus behavior of P_{2m}(1, R2 e^{i phi2}, R3 e^{i phi3}) at the saddle point!

Exact Maillot theorem for m(a + b*U + c*V) where (a, b, c) form a triangle with opposite angles (t1, t2, t3):
  m(a + b*U + c*V) = (1/pi) * [ t1*log(a) + t2*log(b) + t3*log(c) + (1/2)*(Cl2(2*t1) + Cl2(2*t2) + Cl2(2*t3)) ]?
  Wait: check with a = b = c = 1 (t1=t2=t3=pi/3):
  (1/(2*pi)) * 3 * Cl2(2*pi/3) = (3/(2*pi)) * (2/3)*Cl2(pi/3) = Cl2(pi/3)/pi = 0.3230659472! YES!
  Wait: check a=1, b=r2, c=r3 with Legendre transform:
  L(k1, k2, k3) = m(1 + r2*U + r3*V) - (k2/n)*log(r2) - (k3/n)*log(r3)
                = (1/(2*pi)) * [ Cl2(2*pi*k1/n) + Cl2(2*pi*k2/n) + Cl2(2*pi*k3/n) ] !!
  WAIT A SECOND!! Notice that for an OBTUSE triangle (k1 > n/2, so 2*pi*k1/n > pi),
  Cl2(2*pi*k1/n) is NEGATIVE, wait: is Maillot's formula Cl2(2*t1)+Cl2(2*t2)+Cl2(2*t3) still valid for obtuse triangles?!
  YES! Because 2*t1 = 2*pi - 2*(t2+t3), and Cl2(2*pi*k1/n) + Cl2(2*pi*k2/n) + Cl2(2*pi*k3/n) = 2 * Vol(ideal tetrahedron) > 0 ALWAYS for t1+t2+t3=pi!
  Wait: check why n=36 val_mid gave 0.250673 earlier! Let's check exact indices!
"""
import mpmath as mp
import numpy as np
import sys
sys.path.insert(0, r'd:\solving-conjectures\work\bardakov-kozlovskaya\verify')
import engine

mp.mp.dps = 50

def maillot_legendre(k1, k2, k3, n):
    t1, t2, t3 = mp.pi*k1/n, mp.pi*k2/n, mp.pi*k3/n
    # Note: in Bloch-Wigner dilogarithm, Vol = (1/2)*(Cl2(2t1) + Cl2(2t2) + Cl2(2t3))
    return (mp.clsin(2, 2*t1) + mp.clsin(2, 2*t2) + mp.clsin(2, 2*t3)) / (2*mp.pi)

print("=== 1. Verify Maillot hyperbolic-volume formula against exact engine coefficients (n=36, 60) ===")
for n in [36, 60]:
    coeffs = engine.compute_elementary_coeffs(n) # check return format first!
    # Wait: let's inspect coeffs type safely
    if isinstance(coeffs, dict):
        sample_k = list(coeffs.keys())[:3]
        print(f"n={n}: engine.compute_elementary_coeffs dict sample keys = {sample_k}")
        break

# Inspect engine functions and compare!
for n in [24, 36, 60]:
    coeffs = engine.compute_elementary_coeffs(n)
    print(f"\nn={n}:")
    # Pick 6 diverse rays (a, b, c) with a+2b+3c = n
    rays = [
        (0, 0, n//3),           # equilateral (1/3, 1/3, 1/3)
        (n//6, n//6, n//6),     # k1=n/2, k2=n/3, k3=n/6 (30-60-90 right triangle!)
        (n//3, n//6, n//9 if n%9==0 else n//6),
        (n//2, n//12, n//9 if n%9==0 else n//6),
        (0, n//4, n//6),        # a=0
        (n//2, 0, n//6),        # b=0
    ]
    for a, b, c in rays:
        if a + 2*b + 3*c != n or a < 0 or b < 0 or c < 0:
            # adjust a
            a = n - 2*b - 3*c
            if a < 0: continue
        val = coeffs[(a, b, c)]
        k1, k2, k3 = a + b + c, b + c, c
        emp = float(mp.log(abs(val)) / (n*n))
        pred = float(maillot_legendre(k1, k2, k3, n))
        print(f"  (a,b,c)=({a:2d},{b:2d},{c:2d}) [k=({k1:2d},{k2:2d},{k3:2d})]: log|c|/n^2 = {emp:.6f} vs Maillot = {pred:.6f} (diff={emp-pred:+.6f})")

print("\n=== 2. Torus amplitude & phase at EXACT saddle point for right triangle (k1,k2,k3)=(1/2, 1/3, 1/6) ===")
# k1/n = 1/2, k2/n = 1/3, k3/n = 1/6 -> r1=1, r2=sqrt(3)/2, r3=1/2!
# Wait: r2 = sin(pi/3)/sin(pi/2) = sqrt(3)/2 = 0.866025, r3 = sin(pi/6)/sin(pi/2) = 0.5!
r2 = mp.sqrt(3) / 2
r3 = mp.mpf('0.5')

def eval_Pn_from_e(n, e2, e3):
    roots = mp.polyroots([mp.mpc(1), mp.mpc(-1), mp.mpc(e2), mp.mpc(-e3)], maxsteps=200, extraprec=50)
    x, y, z = roots
    alpha = x**(mp.mpf(1)/n)
    beta = y**(mp.mpf(1)/n)
    prod = mp.mpc(1)
    for k in range(n):
        prod *= (z - (alpha + mp.exp(2*mp.pi*1j*k/n)*beta)**n)
    return prod

for n in [12, 24, 36]:
    R2 = r2**n
    R3 = (r2 * r3)**n
    # Scan phi2 in [0, 2*pi) with phi3 = 0, and phi3 in [0, 2*pi) with phi2 = pi (even) / 0 (odd)
    v_ev_0 = abs((eval_Pn_from_e(n, -R2, R3) + eval_Pn_from_e(n, -R2, -R3))/2)
    v_od_0 = abs((eval_Pn_from_e(n,  R2, R3) - eval_Pn_from_e(n,  R2, -R3))/2)
    print(f"\nn={n} (r2=sqrt(3)/2, r3=1/2): P_even(-R2,R3)={mp.nstr(v_ev_0,5)}, -P_odd(R2,R3)={mp.nstr(v_od_0,5)}")
    ev_prof = []
    od_prof = []
    for j in range(9):
        phi2 = mp.pi * j / 8
        vp = eval_Pn_from_e(n, R2*mp.exp(1j*phi2),  R3)
        vm = eval_Pn_from_e(n, R2*mp.exp(1j*phi2), -R3)
        ev_prof.append(f"{float(abs((vp+vm)/2)/v_ev_0):.3f}")
        od_prof.append(f"{float(abs((vp-vm)/2)/v_od_0):.3f}")
    print("  |P_even|/peak for phi2 = 0..pi (step pi/8):", " ".join(ev_prof))
    print("  |P_odd |/peak for phi2 = 0..pi (step pi/8):", " ".join(od_prof))
