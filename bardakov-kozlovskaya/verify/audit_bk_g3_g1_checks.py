"""
Fast 30-50 digit independent adversarial verification & counterexample checks for:
  1. G3-ABC Sub-Lemma G3.3: exact Poisson identity + BOTH y>0 and y<0 poles!
  2. G1 Lemma G1.1: Master HS identity exact convention check!
  3. G1 Lemma G1.2: S_q(d, X) non-monotonicity in d>=2, M_hi counterexample, B1/B2 d>=2 bug, B2b t2/t1 < 1.5*a1!
  4. G1 Lemma G1.3: constants 1.72 & 0.90 and c>=16 cout_l clamp!
"""
import math
import mpmath as mp

mp.mp.dps = 30

print("=== CHECK 1: EXACT POISSON IDENTITY FOR G3.3 (y > 0 vs y < 0 POLES!) ===")
# Test on a moderate-n VALID Piece G3-A point (or G3-C point) where we can compute exact log F_pure - ct*Psi_fin!
# Take a clean G3-A1/A2 or G3-C point with smallish n:
# G3-A2: D=2, r0=1, c=40, p0=1, j=42, kap=0.40 => n ~ 2*pi^2*42*82/(40*0.40) = 4248!
# Check G3 conditions: ct=40>=40, r0*kap=0.40<=0.45, r0^3*c*kap^2=6.4<=60, ct*tau^2/16 = 40*256/16 = 640 >= 40!
def test_g3_poisson(label, n, c, j, r0, p0, D, u0, theta):
    th_c = c * mp.pi / n
    psi_j = j * mp.pi / n
    al1 = (j + c) * mp.pi / n
    r1, r2, r3 = mp.sin(al1), mp.sin(psi_j), mp.sin(th_c)
    sj = r2 / r1
    kap = 2 * mp.pi * r1 * r2 / r3
    ct = mp.mpf(c) / r0
    tau = c * r0 * kap
    assert r0 * j == c * p0 + D and math.gcd(p0, r0) == 1
    assert c * kap**2 < 42 and ct >= 10 * D**2 and r0 * kap <= mp.mpf('0.45')
    assert r0**3 * c * kap**2 <= 60 and ct * tau**2 / (D**4) >= 40
    print(f"\n[{label}] n={n}, c={c}, j={j}, r0={r0}, D={D}, ct={float(ct):.1f}, kap={float(kap):.4f}, tau={float(tau):.3f}")

    def zeta(u):
        w = sj * mp.exp(2j * u * th_c)
        return r0 * n * mp.log((r3 / r1) / (1 - w)) + 1j * mp.pi * r0 * j - 2j * mp.pi * c * p0 * (u + mp.mpf('0.5'))

    re0 = mp.re(zeta(u0))
    Lam = mp.exp(-re0)
    # Compute exact discrete sum over ALL grid points u_k = (k + shift)/c in (-n/(2c), n/(2c)]
    # Here r0 = 1 so F_pure(lam * e^{i(psi_0 + theta)}) = (-w^{-1} (1 - (-1)^j w)^2) * prod_{inner} (1 - a_k w) * prod_{outer} (1 - a_k^{-1} w^{-1})
    # Note: for u_k = +/- 1/2, a_k = (-1)^j, so the boundary factor (1 - (-1)^j w)^2 corresponds to u_k = -1/2 and +1/2!
    # Wait: let's compute the EXACT discrete log-sum sum_{k} log(1 + Lam * e^{zeta(u_k) + i*theta}) (for r0=1)
    # minus the exact continuum integral ct * int_{-Ln}^{Ln} log(1 + Lam * e^{zeta(u) + i*theta}) du!
    assert r0 == 1
    shift = mp.mpf(0) if c % 2 == 0 else mp.mpf('0.5')
    S_disc = mp.mpc(0)
    for k in range(-n // 2, n // 2 + 1):
        x = k + shift
        if x <= -n / 2 + 1e-9 or x > n / 2 + 1e-9:
            continue
        u_k = x / c
        # For |x| == c/2, factor is 1 + Lam*e^{i*theta} (if we include both -1/2 and +1/2)
        z_k = Lam * mp.exp(zeta(u_k) + 1j * theta)
        S_disc += mp.log(1 + z_k)
    # Continuum integral: split at [-Ln, -2, -u0, u0, 2, Ln] for fast high precision
    Ln = mp.mpf(n) / (2 * c)
    f_int = lambda u: mp.log(1 + Lam * mp.exp(zeta(u) + 1j * theta))
    S_cont = ct * (mp.quad(f_int, [-Ln, -5, -1, -u0, 0, u0, 1, 5, Ln]))
    diff_exact = S_disc - S_cont
    print(f"  EXACT (S_disc - ct*int) = {float(mp.re(diff_exact)):+.12e} {float(mp.im(diff_exact)):+.12e}i")

    # Find ALL poles near +/- u0 in BOTH y > 0 AND y < 0!
    im_p = mp.im(zeta(u0))
    im_m = mp.im(zeta(-u0))
    def find_poles_near(u_sign, branch_y_sign):
        u_ref = u_sign * u0
        im_ref = mp.im(zeta(u_ref))
        # Find nearest 3 odd multiples of pi on the appropriate side!
        # For D > 0:
        #   at +u0: y>0 goes inward (x < u0 => Im decreases!), y<0 goes outward (x > u0 => Im increases!)
        #   at -u0: y>0 goes inward (x > -u0 => Im increases!), y<0 goes outward (x < -u0 => Im decreases!)
        # Wait: direction of Im change along Re=re0:
        #   sign(d Im / d|y|) = -u_sign * branch_y_sign * sign(D)
        dir_im = -u_sign * branch_y_sign * (1 if D > 0 else -1)
        base = int(mp.floor((im_ref + theta) / mp.pi))
        if dir_im > 0:
            first_odd = base + 1 if base % 2 == 0 else base + 2
            targets = [(first_odd + 2 * m) * mp.pi for m in range(3)]
        else:
            first_odd = base if base % 2 != 0 else base - 1
            targets = [(first_odd - 2 * m) * mp.pi for m in range(3)]
        found = []
        for idx, tgt in enumerate(targets):
            dx = (tgt - im_ref - theta) / (2 * mp.pi * D)
            x0 = u_ref + dx
            eta0 = abs(tau * (x0**2 - u0**2) / (2 * mp.pi * D))
            u_g = x0 + 1j * branch_y_sign * min(max(eta0, mp.mpf('0.002')), mp.mpf('0.4'))
            try:
                root = mp.findroot(lambda u: zeta(u) - re0 + 1j * theta - 1j * tgt, u_g)
                if abs(zeta(root) - re0 + 1j * theta - 1j * tgt) < mp.mpf('1e-20'):
                    if branch_y_sign * mp.im(root) > 0:
                        # Contribution: -log(1 - e^{2*pi*i*sigma*ct*(u_p + 1/2)}) or +log(...)? Check sign!
                        sig = branch_y_sign
                        # Note: for even c, grid is k/c so e^{2*pi*i*c*u_k} = 1; for odd c, u_k = (k+1/2)/c so e^{2*pi*i*c*(u_k+1/2)} = 1!
                        phase_shift = mp.mpf(0) if c % 2 == 0 else mp.mpf('0.5')
                        val = mp.log(1 - mp.exp(2j * mp.pi * sig * ct * (root + phase_shift)))
                        found.append((idx + 1, root, val))
            except Exception as ex:
                pass
        return found

    sum_ypos = mp.mpc(0)
    sum_yneg = mp.mpc(0)
    for us, lab_u in [(+1, "+u0"), (-1, "-u0")]:
        for ys, lab_y in [(+1, "y>0"), (-1, "y<0")]:
            poles = find_poles_near(us, ys)
            for idx, rt, val in poles:
                exp_v = 2 * mp.pi * ct * abs(mp.im(rt))
                if idx <= 2:
                    print(f"    [{lab_u}, {lab_y}, m={idx}] u_p={float(mp.re(rt)):+.5f}{float(mp.im(rt)):+.5f}i, 2pi*ct*|y|={float(exp_v):6.2f}, log={float(mp.re(val)):+.6e}{float(mp.im(val)):+.6e}i")
                if ys > 0: sum_ypos += val
                else:      sum_yneg += val
    print(f"  Sum(y>0 ONLY)        = {float(mp.re(sum_ypos)):+.12e} {float(mp.im(sum_ypos)):+.12e}i")
    print(f"  Sum(y<0 ONLY)        = {float(mp.re(sum_yneg)):+.12e} {float(mp.im(sum_yneg)):+.12e}i")
    print(f"  Sum(y>0) + Sum(y<0)  = {float(mp.re(sum_ypos+sum_yneg)):+.12e} {float(mp.im(sum_ypos+sum_yneg)):+.12e}i")
    print(f"  Diff vs -(y>0 + y<0) = {float(abs(diff_exact + (sum_ypos + sum_yneg))):.3e}")
    print(f"  Diff vs +(y>0 + y<0) = {float(abs(diff_exact - (sum_ypos + sum_yneg))):.3e}")
    print(f"  Diff vs +(y>0) ONLY  = {float(abs(diff_exact - sum_ypos)):.3e}")
    print(f"  Diff vs -(y>0) ONLY  = {float(abs(diff_exact + sum_ypos)):.3e}")

test_g3_poisson("G3-A (D=2, r0=1, c=40, kap~0.40, theta=0)", n=4248, c=40, j=42, r0=1, p0=1, D=2, u0=mp.mpf('0.45'), theta=mp.mpf('0.0'))
test_g3_poisson("G3-A (D=2, r0=1, c=40, kap~0.40, theta=1.0)", n=4248, c=40, j=42, r0=1, p0=1, D=2, u0=mp.mpf('0.45'), theta=mp.mpf('1.0'))

print("\n=== CHECK 2: G1 LEMMA G1.1 MASTER HS IDENTITY (EXACT CONVENTION) ===")
def HS_rhs(t, nu):
    eps = mp.sqrt(t / (mp.pi * nu))
    def F(u):
        D1 = ((t + u)**2 + 1) * ((t - u)**2 + 1)
        D2 = (t**2 + (1 + u)**2) * (t**2 + (1 - u)**2)
        return 8 * t**2 * ((t**2 + 1) * (t**2 - 3) - u**4) / (D1 * D2)
    return mp.quad(lambda z: (2 / mp.sqrt(2 * mp.pi)) * z**2 * mp.exp(-mp.mpf('0.5') * z**2) * F(eps * z), [0, mp.inf])

# Test exact D(T, nu) = int_{-1}^1 cos(pi*nu*u)*exp(-0.5*T*(1-u^2)) du - 2*int_0^inf cos(pi*nu*(1+y))*exp(-T*(y+y^2/2)) dy
# vs factor 2*pi*nu and T vs 2*T!
for t_val, nu_val in [(mp.mpf('0.5'), 1), (mp.mpf('1.2'), 2)]:
    T_val = t_val * mp.pi * nu_val
    rhs = HS_rhs(t_val, nu_val)
    # Convention A: cos(2*pi*nu*u) with exp(-T*(1-u^2))? Or cos(pi*nu*u) with exp(-0.5*T*(1-u^2))?
    for name, cfreq, Tfac_in, Tfac_out in [
        ("pi*nu, T/2", mp.pi, mp.mpf('0.5'), mp.mpf('1')),
        ("pi*nu, T",   mp.pi, mp.mpf('1'),   mp.mpf('2')),
        ("2pi*nu, T",  2*mp.pi, mp.mpf('1'), mp.mpf('2')),
        ("2pi*nu, 2T", 2*mp.pi, mp.mpf('2'), mp.mpf('4')),
    ]:
        I_in = mp.quad(lambda u: mp.cos(cfreq * nu_val * u) * mp.exp(-Tfac_in * T_val * (1 - u**2)), [-1, 1])
        I_out = 2 * mp.quad(lambda y: mp.cos(cfreq * nu_val * (1 + y)) * mp.exp(-Tfac_out * T_val * (y + mp.mpf('0.5') * y**2)), [0, mp.inf])
        lhs = 2 * mp.pi**2 * nu_val**2 * (I_in - I_out)
        if abs(lhs) > 1e-12:
            print(f"  (t={float(t_val)},nu={nu_val}) [{name:11s}]: LHS={float(lhs):+.10f}, RHS={float(rhs):+.10f}, ratio={float(rhs/lhs):+.8f}, diff={float(abs(lhs-rhs)):.2e}")

print("\n=== CHECK 3: FAST CLOSED-FORM D0(T) AND S_q(d, X) NON-MONOTONICITY IN d >= 2 ===")
def D0_fast(T):
    # int_{-1}^1 exp(-T/2*(1-u^2)) du = exp(-T/2) * sqrt(2*pi/T) * erfi(sqrt(T/2))
    # 2*int_0^inf exp(-T*(y+y^2/2)) dy = exp(T/2) * sqrt(2*pi/T) * erfc(sqrt(T/2))
    s = mp.sqrt(T / 2)
    I_in = mp.exp(-T / 2) * mp.sqrt(2 * mp.pi / T) * mp.erfi(s)
    I_out = mp.exp(T / 2) * mp.sqrt(2 * mp.pi / T) * mp.erfc(s)
    return I_in - I_out

def S_q_fast(d, X):
    s = mp.mpf(0)
    for l in range(1, 300):
        s += (mp.mpf(d) / (2 * l)) * abs(D0_fast(mp.mpf('0.5') * l * d * X))
    return s

print("  D0(T) scan:")
for T_s in ['0.5', '1.0', '1.5', '2.0', '2.2', '2.5', '3.0', '4.0', '6.0', '10.0', '15.0', '30.0']:
    T = mp.mpf(T_s)
    v = D0_fast(T)
    print(f"    T={T_s:4s}: D0(T) = {float(v):+.7f}, T*|D0(T)| = {float(T*abs(v)):.6f}, T^2*|D0(T)| = {float(T*T*abs(v)):.5f}")

print("  S_q(d, X) for d = 1, 2, 3, 4, 5, 8:")
for X_test in [mp.mpf('2.0'), mp.mpf('2.2'), mp.mpf('2.5'), mp.mpf('3.0'), mp.mpf('4.0'), mp.mpf('19.0')]:
    vals = [(d, float(S_q_fast(d, X_test))) for d in [1, 2, 3, 4, 5, 8]]
    print(f"    X={float(X_test):5.1f}: " + ", ".join(f"d={d}:{v:.5f}" for d, v in vals))

print("\n=== CHECK 4: G1 LEMMA G1.3 CONSTANTS 1.72 & 0.90 ===")
f1 = lambda u: 2 * (u**mp.mpf('1.5')) / mp.sinh(mp.mpf('0.99') * u)
f2 = lambda u: (u**2) / mp.sinh(mp.mpf('0.99') * u)
print(f"  max 2*u^(3/2)/sinh(0.99*u)  = {float(max(f1(mp.mpf(k)/500) for k in range(1, 3000))):.6f}  (script L1486 uses 1.72!)")
print(f"  max u^2/sinh(0.99*u)        = {float(max(f2(mp.mpf(k)/500) for k in range(1, 3000))):.6f}  (script L1513 uses 0.90!)")
f1s = lambda u: 2 * (u**mp.mpf('1.5')) / mp.sinh(mp.mpf('0.998') * u)
f2s = lambda u: (u**2) / mp.sinh(mp.mpf('0.998') * u)
print(f"  max 2*u^(3/2)/sinh(0.998*u) = {float(max(f1s(mp.mpf(k)/500) for k in range(1, 3000))):.6f}  (script L1686 uses 1.72!)")
print(f"  max u^2/sinh(0.998*u)       = {float(max(f2s(mp.mpf(k)/500) for k in range(1, 3000))):.6f}  (script L1693 uses 0.90!)")
