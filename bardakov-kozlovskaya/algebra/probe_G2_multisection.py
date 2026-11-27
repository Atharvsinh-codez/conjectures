"""
Probe the EXACT q-multisection identity on EXACT FINITE-n D = 0 Major Arcs:
    F(w) = G_q(w^q)^{1/q} * exp(R_{non-q}(w)),   G_q(w^q) := prod_{k=0}^{n-1} (1 - a_k^q w^q).
Check for ALL d >= 5, q = 1..10, d*q^2*kap >= 0.85, d*(q^2*kap)^2 <= 60:
  1. Exact saddle Lam_* = lam_*^q of the POSITIVE-ROOT polynomial G_q(W) (W = -w^q)!
  2. Exact curvature H_q = (1/q) * sum_k Re(a_k^q Lam_* / (1 + a_k^q Lam_*^2...)) of G_q(W)^{1/q} in psi = q*phi - pi!
  3. Exact max |Re(R_{non-q})| and max |Im(R_{non-q})| across the ENTIRE circle |w| = lam_* (and on the peak support)!
  4. Wait: check if R_{non-q} blows up ONLY near the zeros a_k^{-1} (where G_q ALSO vanishes!!),
     and compare |F(w)| vs |G_q(w^q)|^{1/q} directly!!
"""
import math
import cmath

def exact_finite_n_roots(n, c, j):
    th_c = c * math.pi / n
    psi_j = j * math.pi / n
    alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1
    log_Rz = n * math.log(r3 / r1)
    roots = []
    shift = 0.0 if (c % 2 == 0) else 0.5
    sgn_ang = 0.0 if (j % 2 == 0) else math.pi
    k = 0
    while True:
        added = False
        klist = [0] if (k == 0 and c % 2 == 0) else ([0, -1] if (k == 0 and c % 2 == 1) else [k, -k - (1 if c % 2 == 1 else 0)])
        for kk in klist:
            eta = 2.0 * (kk + shift) * math.pi / n
            if abs(eta) > math.pi + 1e-12:
                continue
            z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
            log_r = log_Rz - n * math.log(abs(z))
            ang = sgn_ang - n * math.atan2(z.imag, z.real)
            if log_r > -55.0:
                roots.append((log_r, ang))
                added = True
        if not added and k > c:
            break
        k += 1
        if k > n // 2 + 2:
            break
    return roots, 2.0 * math.pi * r1 * r2 / r3

def eval_F_and_Gq(roots, c, q, lam, phi):
    """Return (log|w^{-c}F|, arg(w^{-c}F), log|W^{-d}Gq^{1/q}|, arg(W^{-d}Gq^{1/q})) at w = lam*e^{i phi}."""
    d = c // q
    log_lam = math.log(lam)
    la_F = -c * log_lam
    ph_F = -c * phi
    la_G = -c * log_lam
    ph_G = -c * phi
    # Wait: for G_q(w^q) = prod (1 - a_k^q w^q), track each factor continuously in phi?
    # Better: compute log(1 - a_k w) and (1/q)*log(1 - a_k^q w^q) directly for each k!
    for lr, ang in roots:
        lm = lr + log_lam
        th = ang + phi
        if lm < -45.0:
            continue
        elif lm > 45.0:
            la_F += lm
            ph_F += th + math.pi
            la_G += lm
            ph_G += th + math.pi / q  # wait: phase of (-a_k^q w^q)^{1/q} mod 2*pi/q!
        else:
            r = math.exp(lm)
            z = complex(1.0 - r * math.cos(th), -r * math.sin(th))
            la_F += math.log(abs(z))
            ph_F += math.atan2(z.imag, z.real)
            # q-th power:
            rq = r ** q
            thq = q * th
            if rq > 1e18:
                la_G += lm
                ph_G += th + math.pi / q
            elif rq < 1e-18:
                pass
            else:
                zq = complex(1.0 - rq * math.cos(thq), -rq * math.sin(thq))
                la_G += math.log(abs(zq)) / q
                ph_G += math.atan2(zq.imag, zq.real) / q
    return la_F, ph_F, la_G, ph_G

def solve_Gq_saddle(roots, c, q):
    """Solve (1/q) * sum_k Re(a_k^q * Lam / (1 + a_k^q * Lam)) = d for Lam = lam^q > 0."""
    d = c // q
    x = 0.0  # x = log(lam)
    for _ in range(40):
        val = -float(d)
        der = 0.0
        for lr, ang in roots:
            lmq = q * (lr + x)
            angq = q * ang
            if lmq < -45.0:
                continue
            elif lmq > 45.0:
                val += 1.0 / q
            else:
                u = math.exp(lmq) * complex(math.cos(angq), math.sin(angq))
                # W = -w^q = +Lam, so 1 - a_k^q w^q = 1 + u
                t = u / (1.0 + u)
                val += t.real / q
                der += (t / (1.0 + u)).real  # d/dx [t/q] = q * (u/(1+u)^2) / q = u/(1+u)^2!
        if abs(val) < 1e-12 or der <= 0:
            break
        x -= max(-0.2, min(0.2, val / der))
    # Compute H_psi = (1/q) * sum_k Re(u / (1+u)^2) (exact curvature in psi = q*phi - pi!):
    H_psi = 0.0
    for lr, ang in roots:
        lmq = q * (lr + x)
        angq = q * ang
        if abs(lmq) > 45.0:
            continue
        u = math.exp(lmq) * complex(math.cos(angq), math.sin(angq))
        H_psi += (u / ((1.0 + u) ** 2)).real / q
    return math.exp(x), H_psi

def mk_n(c, j, kt):
    n = int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt)))
    return max(132, n + (n % 2))

print("=== EXACT q-MULTISECTION PROBE ON D = 0 MAJOR ARCS (d >= 5) ===")
grid = [
    (1, 2, 5, 0.50), (1, 2, 5, 1.50), (1, 2, 6, 0.80), (1, 2, 10, 0.30), (1, 2, 50, 0.05),
    (2, 3, 5, 0.25), (2, 3, 5, 0.80), (2, 3, 6, 0.50), (2, 3, 10, 0.15), (2, 3, 50, 0.03),
    (3, 4, 5, 0.15), (3, 4, 5, 0.60), (3, 4, 10, 0.10), (3, 4, 30, 0.03),
    (5, 7, 5, 0.10), (5, 7, 5, 0.30), (5, 7, 15, 0.04),
    (7, 9, 5, 0.06), (10, 11, 5, 0.03), (10, 11, 5, 0.15),
]
for q, p0, d, kt in grid:
    c = d * q
    j = d * p0
    n = mk_n(c, j, kt)
    roots, kap = exact_finite_n_roots(n, c, j)
    lam, H_psi = solve_Gq_saddle(roots, c, q)
    sgn_pref = 1.0 if ((j * (c + 1) + c) % 2 == 0) else -1.0
    # Evaluate F vs G_q^{1/q} across M points on the circle:
    M = max(256, 64 * q)
    # Specifically measure max |log|F| - log|G_q^{1/q}|| and max |arg(F) - arg(G_q^{1/q})| where |G_q^{1/q}| >= exp(-6) * peak!
    # Also check peak value of |F/G_q^{1/q} - 1|!
    laG_list = []
    laF_list = []
    phF_list = []
    phG_list = []
    for i in range(M):
        phi = 2.0 * math.pi * (i + 0.5) / M
        laF, phF, laG, phG = eval_F_and_Gq(roots, c, q, lam, phi)
        laF_list.append(laF); phF_list.append(phF)
        laG_list.append(laG); phG_list.append(phG)
    mxG = max(laG_list)
    mxF = max(laF_list)
    # Check max |laF - laG| on |G| >= exp(-6)*peak AND globally (laF - laG upper bound!!):
    max_re_diff_pk = max(abs(laF_list[i] - laG_list[i]) for i in range(M) if laG_list[i] >= mxG - 6.0)
    max_re_diff_all = max(laF_list[i] - laG_list[i] for i in range(M))
    # Also check exact |F(phi_nu)/G_q(phi_nu)^{1/q} - 1| at the peak phi0 = pi/q:
    phi0 = math.pi if (q % 2 == 1) else (math.pi / q)
    laF0, phF0, laG0, phG0 = eval_F_and_Gq(roots, c, q, lam, phi0)
    # Exact phase difference at peak (mod 2*pi):
    epF0 = (phF0 + (0.0 if sgn_pref > 0 else math.pi) + math.pi) % (2.0 * math.pi) - math.pi
    # Check rat_F:
    s_re = sum(math.exp(laF_list[i] - mxF) * math.cos(phF_list[i] + (0.0 if sgn_pref > 0 else math.pi)) for i in range(M) if laF_list[i] - mxF > -700)
    s_abs = sum(math.exp(laF_list[i] - mxF) for i in range(M) if laF_list[i] - mxF > -700)
    print(f"q={q:2d} d={d:2d} kap={kap:.3f} dq^2k={d*q*q*kap:5.2f} dq^4k^2={d*(q*q*kap)**2:6.2f} | "
          f"lam={lam:.4f} H_psi={H_psi:6.2f} | |Re R|_pk={max_re_diff_pk:.4f} max(Re R)_all={max_re_diff_all:+.4f} "
          f"|epF0|={abs(epF0):.4f} | rat={s_re/s_abs:+.5f}")
