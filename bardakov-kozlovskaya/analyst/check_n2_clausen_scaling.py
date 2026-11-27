"""
Test EXACT integral / analytic scaling of p_{2m}(x, y, z) when (x, y, z) = (r1^n, r2^n, r3^n)
AND check the exact constant lim_{n->inf} (1/n^2) log |c_{0,0,n/3}^{(n)}|!!
WAIT A SECOND!!
For 6 | n: c_{0,0,n/3}^{(n)} = p_n(1, omega, omega^2) = prod_{k=0}^{n-1} (omega^2 - (1 + zeta^k * omega^{1/n})^n)!
WAIT: |(1 + e^{i*theta})| = 2 |cos(theta/2)|!
So log |p_n(1, omega, omega^2)| ~ n * (n / (2*pi)) * integral_{-2*pi/3}^{2*pi/3} log(2 cos(theta/2)) d(theta)!!
HOLD ON!! COMPUTE THAT INTEGRAL RIGHT NOW:
  (1/(2*pi)) * integral_{-2*pi/3}^{2*pi/3} log(2 cos(theta/2)) d(theta)
  = (2/pi) * Cl_2(pi/3) / 2 ... = (1/pi) * Cl_2(pi/3) = 0.3230659472... !!
CHECK IF log|c_{0,0,n/3}| / n^2 CONVERGES TO Cl_2(pi/3) / pi = 0.3230659472... !!
"""
import mpmath as mp

mp.mp.dps = 50
# Clausen function Cl_2(pi/3) = Im(Li_2(e^{i*pi/3}))
cl2 = mp.clsin(2, mp.pi/3)
I_exact = mp.quad(lambda th: mp.log(2*mp.cos(th/2)), [-2*mp.pi/3, 2*mp.pi/3]) / (2*mp.pi)
print(f"Cl_2(pi/3)/pi = {float(cl2/mp.pi):.10f}")
print(f"Exact integral (1/2pi) int_{-2pi/3}^{2pi/3} log(2cos(th/2)) dth = {float(I_exact):.10f}")

# Check exact p_n(1, omega, omega^2) for n = 6, 12, 18, 24, 30, 36, 48, 60!
omega = mp.exp(2*mp.pi*1j/3)
for n in [6, 12, 18, 24, 30, 36, 48, 60, 120]:
    # alpha = 1, beta = exp(2*pi*i / (3n)), z = omega^2
    beta = mp.exp(2*mp.pi*1j / (3*n))
    prod = mp.mpc(1)
    for k in range(n):
        zeta_k = mp.exp(2*mp.pi*1j*k / n)
        prod *= (omega**2 - (1 + zeta_k * beta)**n)
    val = prod.real
    assert abs(prod.imag) < mp.mpf('1e-20') * abs(val)
    print(f"n={n:3d}: sign={1 if val>0 else -1:+d}, log|p_n(1,w,w^2)|/n^2 = {float(mp.log(abs(val))/(n*n)):.10f}")
