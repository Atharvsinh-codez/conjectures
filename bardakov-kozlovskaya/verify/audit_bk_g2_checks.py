"""
Independent 30-digit adversarial verification & counterexample checks for updated algebra/certify_G2_D0.py:
  1. Line 747-748 & 775: U_y(x, 0) = -dpsi'(x) is INCREASING (not decreasing!) on x >= c/2!
     Test U(x, y) - U(x, 0) vs C_psi*kap^2*y + C_u2*(kap/c)*y^2 on [c/2, 0.060*n] and [0.15*n, 0.25*n]!
  2. Line 1083-1085: u_dead_top at Q(x) = 1.00 and u_side at x = 0.25*n vs dU_top!
  3. Line 1016: U(0, H) - U(0, 0) vs 0.5*(c_y_lo - C_psi*kap_cap)*Y^2!
  4. Line 960-974: sign swap of Lam_tilde < 1 vs Lam_tilde > 1 in q_dpsi_x0!
  5. Lemma G2.7 (L1311-1352): verify sin_se_coeff, c_se_g2, F_pref, and multiplicative E_cubic!
"""
import math
import mpmath as mp

mp.mp.dps = 30

def check_g2_point(label, d, q, p0, n):
    c = d * q
    j = d * p0
    assert math.gcd(p0, q) == 1 and d >= 2 and q >= 2 and j > c and n % 2 == 0 and n >= 132
    assert c <= (n - 3) // 3
    th_c = c * mp.pi / n
    psi_j = j * mp.pi / n
    al1 = (j + c) * mp.pi / n
    r1, r2, r3 = mp.sin(al1), mp.sin(psi_j), mp.sin(th_c)
    s = r2 / r1
    kap = 2 * mp.pi * r1 * r2 / r3
    q2k = q * q * kap
    T = mp.mpf('0.5') * c * q * kap
    ck2 = c * kap * kap
    print(f"\n[{label}] d={d}, q={q}, p0={p0}, c={c}, j={j}, n={n}")
    print(f"  kap={float(kap):.6f}, q^2*kap={float(q2k):.4f}, T={float(T):.4f}, c*kap^2={float(ck2):.3f}, al1/pi={float(al1/mp.pi):.5f}")
    assert kap < 2.5 and ck2 < 42 and q2k <= 2.0 and c * q * kap >= 0.85

    def log_b(z):
        eta = 2 * mp.pi * z / n
        return n * mp.log((r3 / r1) / (1 - s * mp.exp(1j * eta))) - 2j * mp.pi * (mp.mpf(j) / c) * (z - mp.mpf(c) / 2)

    def U(x, y): return mp.re(log_b(x + 1j * y))
    def V(x, y): return mp.im(log_b(x + 1j * y))
    def U_y_at_0(x):
        eta = 2 * mp.pi * x / n
        D = 1 - 2 * s * mp.cos(eta) + s * s
        return 2 * mp.pi * (mp.mpf(j) / c - s * (mp.cos(eta) - s) / D)

    # Constants from certify_G2_D0.py for q >= 2
    C_psi_q2 = mp.mpf('0.080404')
    C_u2_q2 = mp.mpf('1.00645')
    # Check U_y(x, 0) at x = c/2, x = 0.05*n, x = 0.06*n, x = 0.15*n, x = 0.25*n
    cpsi_k2 = C_psi_q2 * kap * kap
    print(f"  C_psi*kap^2 = {float(cpsi_k2):.6e}")
    for name, x in [("0", mp.mpf(0)), ("c/2", mp.mpf(c)/2), ("0.02n", mp.mpf('0.02')*n),
                    ("0.06n", mp.mpf('0.06')*n), ("0.15n", mp.mpf('0.15')*n), ("0.25n", mp.mpf('0.25')*n)]:
        uy = U_y_at_0(x)
        Q = (1 - 2*s*mp.cos(2*mp.pi*x/n) + s*s) / ((r3/r1)**2) - 1
        print(f"    x={name:5s}: Q={float(Q):8.4f}, U(x,0)={float(U(x,0)):+12.3f}, U_y(x,0)={float(uy):+.6e} (vs C_psi*kap^2={float(cpsi_k2):.3e})")

    # Find actual Lam_* = 1.0 approx (or exact saddle for small n, or just lam=1.0 -> log lam = 0)
    # Find actual H_pi where U(x_*, H) = 0 and -V(x_*, H) = pi!
    # Start guess: x_* = c/2, H = pi / (kap)
    try:
        sol = mp.findroot(lambda x, y: (U(x, y), V(x, y) + mp.pi), (mp.mpf(c) / 2, mp.pi / kap))
        xs, H_pi = sol[0], sol[1]
        Y = mp.sqrt(2 * kap / c) * H_pi
        X = mp.sqrt(2 * kap / c) * xs
        dU_0 = U(0, H_pi) - U(0, 0)
        # Compare dU_0 with 0.5 * c_y_eff * Y^2!
        # Also compare U(x, H_pi) - U(x, 0) at x = c/2, x_Q02, x_Q1, 0.25n with C_psi*kap^2*H + C_u2*(kap/c)*H^2!
        quad_claim = C_psi_q2 * kap**2 * H_pi + C_u2_q2 * (kap / c) * H_pi**2
        # Find exact x where Q(x) = 0.20 and Q(x) = 1.00
        def find_x_Q(Q_tgt):
            # 1 - 2s cos(eta) + s^2 = (1 + Q_tgt)*(r3/r1)^2
            rhs = (1 + s*s - (1 + Q_tgt)*(r3/r1)**2) / (2 * s)
            return n * mp.acos(rhs) / (2 * mp.pi)
        x_Q02 = find_x_Q(mp.mpf('0.20'))
        x_Q1  = find_x_Q(mp.mpf('1.00'))
        print(f"  At H_pi={float(H_pi):.4f} (x_*={float(xs):.3f}, X={float(X):.3f}, Y={float(Y):.4f}, 2pi*H/q={float(2*mp.pi*H_pi/q):.3f}):")
        print(f"    dU(0, H_pi)      = {float(dU_0):+.6e} vs 0.5*Y^2 = {float(0.5*Y**2):+.6e} (ratio dU/(0.5*Y^2) = {float(dU_0/(0.5*Y**2)):.4f})")
        for nm, xv in [("c/2", mp.mpf(c)/2), ("x_Q02", x_Q02), ("x_Q1", x_Q1), ("0.25n", mp.mpf('0.25')*n)]:
            dU_x = U(xv, H_pi) - U(xv, 0)
            UxH = U(xv, H_pi)
            print(f"    x={nm:5s}: dU(x,H)={float(dU_x):+10.4f} vs claimed quad={float(quad_claim):8.4f} | U(x,H)={float(UxH):+12.3f}")
    except Exception as ex:
        print(f"  findroot failed: {ex}")

# Test 1: moderate Box G2 point (q=2, d=2, kap=0.50 -> T=2.0)
# n ~ 2*pi^2*j*(j+c)/(c*kap): d=2, q=2, c=4, p0=3, j=6, kap=0.49 => n ~ 2*pi^2*6*10/(4*0.49) = 604
check_g2_point("q=2, d=2, kap~0.49 (T~2.0)", d=2, q=2, p0=3, n=604)

# Test 2: small T=0.45, q=2, d=2, LARGE p0=101 (j/c = 50.5!)
# kap=0.1125 => n ~ 2*pi^2*202*206/(4*0.1125) = 1825300
check_g2_point("q=2, d=2, T~0.45, LARGE j/c=50.5", d=2, q=2, p0=101, n=1825300)

# Test 3: q=5, LARGE c*kap^2 ~ 41.0, small kap (breaks L1016 cy_eff!)
# Want c*kap^2 = 41, q=5, kap = 0.005 => c = 41/0.000025 = 1640000 => d = 328000, p0 = 6, j = 6*d = 1968000
# n ~ 2*pi^2*j*(j+c)/(c*kap) = 2*pi^2*1.968e6*3.608e6 / (1.64e6 * 0.005) = 1.709e10
check_g2_point("q=5, c*kap^2~41 (L1016 test)", d=328000, q=5, p0=6, n=17090000000)

# Test 4: q=2, kap=0.49, high j/n near alp1_q2 = 0.02537*pi, small T -> test Q(x)=1.00!
# Wait: for q=2, d=2, T = 4*kap. To make T small AND alp1 large: alp1 = (j+c)*pi/n ~ kap/(2*pi) * (r3/r2) <= kap/(2*pi),
# so alp1_max = 0.02537*pi requires kap ~ 0.50 (T >= 2.0)! Let's test q=2, d=2, p0=3, kap=0.495 (n=600) AND d=2, q=2, p0=3, T=0.45 (n=2640)!
check_g2_point("q=2, d=2, T~0.45, p0=3", d=2, q=2, p0=3, n=2640)
