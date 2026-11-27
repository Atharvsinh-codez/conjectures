"""
Verify the EXACT 1-line Heat-Kernel Identity for b^G(c, D, gam):
  b^G(c, D, gam) = (gam / (2*c)) * int_0^1 E_Z [ M(gam/2 + sig*sqrt(lam)*Z, pi*D/c)
                                               + M(gam/2, pi*D/c + sig*sqrt(lam)*Z) ] dlam
where sig = sqrt(gam / (2*c)), Z ~ N(0, 1), and M(u, v) = Re[ cosh(u + i*v) / sinh^3(u + i*v) ]!
Check to 1e-12 across small and large (c, D, gam)!
"""
import math
import cmath

def b_r_gauss(c, D, gam):
    th = 2.0 * math.pi * D / c
    B = sum(math.exp(-gam * m * (1.0 - m / c)) * math.cos(m * th) for m in range(1, c) if gam * m * (1.0 - m / c) < 40)
    A = 0.0; m = 1
    while True:
        ex = gam * m * (1.0 + m / c)
        if ex > 40: break
        A += 2.0 * math.exp(-ex) * math.cos(m * th); m += 1
    return B - A

def M_func(u, v):
    # Re[ cosh(u + iv) / sinh^3(u + iv) ]
    w = complex(u, v)
    sh = cmath.sinh(w)
    ch = cmath.cosh(w)
    return (ch / (sh * sh * sh)).real

# Also verify the direct Hubbard-Stratonovich identity for D(T, nu) first!
def D_HS(T, nu):
    # (2 / sqrt(2*pi)) int_{-inf}^{inf} e^{-s^2/2} [ (T + s*sqrt(T))/((T+s*sqrt(T))^2 + (pi*nu)^2) - T/(T^2 + (pi*nu + s*sqrt(T))^2) ] ds
    N = 400; smax = 8.0; ds = 2.0 * smax / N
    b = math.pi * nu; rt = math.sqrt(T)
    tot = 0.0
    for i in range(N):
        s = -smax + (i + 0.5) * ds
        w = math.exp(-0.5 * s * s) / math.sqrt(2.0 * math.pi)
        h = s * rt
        t1 = (T + h) / ((T + h) ** 2 + b * b)
        t2 = T / (T * T + (b + h) ** 2)
        tot += 2.0 * w * (t1 - t2) * ds
    return tot

def b_heat_kernel(c, D, gam):
    sig = math.sqrt(gam / (2.0 * c))
    u0 = 0.5 * gam
    v0 = math.pi * D / c
    N_lam = 40; N_s = 200; smax = 7.5; ds = 2.0 * smax / N_s
    tot = 0.0
    for il in range(N_lam):
        lam = (il + 0.5) / N_lam
        st = sig * math.sqrt(lam)
        e_val = 0.0
        for j in range(N_s):
            s = -smax + (j + 0.5) * ds
            w = math.exp(-0.5 * s * s) / math.sqrt(2.0 * math.pi)
            e_val += w * (M_func(u0 + st * s, v0) + M_func(u0, v0 + st * s)) * ds
        tot += e_val / N_lam
    return (gam / (2.0 * c)) * tot

# Wait: also check b_HS directly (without integrating over lam!):
# sum_k (1/((a+h) + i(b + k*pi*c))) = (1/c) coth((a+h+ib)/c)!
# So b^G(c, D, gam) = (1/2) * E_Z [ Re coth(gam/2 + sig*Z + i*pi*D/c) - Re coth(gam/2 + i*(pi*D/c + sig*Z)) ] !!
# WAIT A SECOND!! CHECK THAT 0-INTEGRAL COTH IDENTITY TOO!!
def b_coth_HS(c, D, gam):
    sig = math.sqrt(gam / (2.0 * c))
    u0 = 0.5 * gam; v0 = math.pi * D / c
    N_s = 300; smax = 8.0; ds = 2.0 * smax / N_s
    tot = 0.0
    for j in range(N_s):
        s = -smax + (j + 0.5) * ds
        w = math.exp(-0.5 * s * s) / math.sqrt(2.0 * math.pi)
        c1 = (cmath.cosh(complex(u0 + sig * s, v0)) / cmath.sinh(complex(u0 + sig * s, v0))).real
        c2 = (cmath.cosh(complex(u0, v0 + sig * s)) / cmath.sinh(complex(u0, v0 + sig * s))).real
        tot += w * (c1 - c2) * ds
    return tot

print("=== VERIFY EXACT COTH-HS AND HEAT-KERNEL IDENTITIES FOR b^G(c, D, gam) ===")
for c, D, gam in [(6, 1, 0.5), (10, 3, 0.2), (15, 1, 0.35), (20, 1, 0.5), (29, 3, 0.15), (100, 4, 0.15)]:
    bg = b_r_gauss(c, D, gam)
    bc = b_coth_HS(c, D, gam)
    bh = b_heat_kernel(c, D, gam)
    # Also check Poisson sum with D_HS!
    T = 0.5 * c * gam
    bp = 0.5 * c * sum(D_HS(T, abs(D + k * c)) for k in range(-20, 21))
    print(f"  c={c:3d} D={D} gam={gam:.2f}: b^G={bg:+.8f} | coth-HS={bc:+.8f} (err={abs(bg-bc):.2e}) | heat={bh:+.8f} (err={abs(bg-bh):.2e}) | D_HS={bp:+.8f}")
