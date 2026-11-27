"""Quick adversarial check of Lemma G2.4 strip pole-crossing in Box G2."""
import math, cmath

def check_point(d, q, p0, q2k, psi):
    c = d * q; j = d * p0; kt = q2k / (q * q)
    T = 0.5 * d * q2k
    n = max(132, int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt))) // 2 * 2)
    th_c = c * math.pi / n; psi_j = j * math.pi / n; alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1; kap = 2.0 * math.pi * r1 * r2 / r3
    # Solve saddle Lam_* for G_q^+
    shift = 0.0 if (c % 2 == 0) else 0.5
    lrs = []
    log_Rz = n * math.log(r3 / r1)
    for k in range(-n // 2, n // 2 + 1):
        x = k + shift
        if x <= -n / 2 + 1e-9 or x > n / 2 + 1e-9: continue
        eta = 2.0 * math.pi * x / n
        lr = n * math.log(abs(complex(1.0 - s * math.cos(eta), -s * math.sin(eta)))) - log_Rz
        if lr < 45.0: lrs.append(-lr)
    lx = 0.0
    for _ in range(40):
        val = -float(d); der = 0.0
        for lr in lrs:
            lmq = q * (lr + lx)
            if lmq < -45.0: continue
            elif lmq > 45.0: val += 1.0 / q
            else:
                eu = math.exp(lmq); t = eu / (1.0 + eu); val += t / q; der += t / (1.0 + eu)
        if abs(val) < 1e-13: break
        lx -= val / der
    Lam = math.exp(q * lx)
    A = math.pi - abs(psi)
    th_t = A / q
    c0 = 0.78 if A < 1.8 else 0.72
    H_t = c0 * th_t / kap
    # Exact log(u(z)) = lx + n*log((r3/r1)/(1 - s*exp(i*2*pi*z/n))) - 2*pi*i*(j/c)*(z - c/2) + i*pi*j + i*th_t
    # Solve log(u(z)) = 0 by Newton starting near x = c/2 * 1.5, y = 0.5 * th_t / kap
    def log_u(z):
        eta = 2.0 * math.pi * z / n
        w = 1.0 - s * cmath.exp(1j * eta)
        # unwrap near real axis
        val = lx + n * (math.log(r3 / r1) - cmath.log(w)) - 2j * math.pi * (j / c) * (z - c / 2.0) + 1j * th_t
        if j % 2 == 1:
            val += 1j * math.pi
        # adjust 2*pi*i multiple to nearest 0
        im = (val.imag + math.pi) % (2.0 * math.pi) - math.pi
        return complex(val.real, im)
    z = complex(0.75 * c, 0.5 * th_t / kap)
    for _ in range(40):
        f = log_u(z)
        h = 1e-6 * max(1.0, abs(z))
        df = (log_u(z + h) - log_u(z - h)) / (2 * h)
        z -= f / df
    err = abs(log_u(z))
    print(f"d={d}, q={q}, q2k={q2k}, psi={psi:.3f} (n={n}, T={T:.2f}, Lam={Lam:.3f}): "
          f"POLE z_p = {z.real:.2f} + {z.imag:.2f}i (x_p/(c/2)={z.real/(c/2):.3f}), "
          f"y_p = {z.imag:.3f} vs H_t = {H_t:.3f} ==> y_p/H_t = {z.imag/H_t:.4f} {'< 1 POLE INSIDE STRIP!!' if z.imag < H_t else '>= 1 safe'} (err={err:.1e})")

for d, q, p0, q2k in [(2, 2, 3, 0.45), (2, 2, 3, 1.0), (2, 2, 3, 1.95), (2, 4, 5, 0.45), (2, 4, 5, 1.95), (3, 2, 3, 0.45), (4, 2, 3, 1.50)]:
    check_point(d, q, p0, q2k, 2.0 * math.pi / 3.0)
    check_point(d, q, p0, q2k, 1.25)
