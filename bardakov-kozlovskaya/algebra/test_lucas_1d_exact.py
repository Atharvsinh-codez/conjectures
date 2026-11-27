"""
Test the EXACT Lucas-Chebyshev factorization of our 1D Circle Integrand:
  P_{2m}(1+t, t, p) = x^{2m} prod_{k=0}^{2m-1} (L_{2m}(zeta^k (p/x^3)^{1/(2m)}) - (1+t-x)/x)
vs the frozen-x (x = 1!) approximation:
  P_{Lucas, x=1}(t, p) = prod_{k=0}^{2m-1} (L_{2m}(zeta^k p^{1/(2m)}) - t)
at j = c (c=2, j=2; c=3, j=3; top_left c=12, j=12 at n=40; top_left c=26, j=26 at n=80)
and bulk (n=80, c=13, j=20)!
"""
import mpmath as mp

mp.mp.dps = 120

def L2m_mp(m, w):
    sq = mp.sqrt(1 + 4 * w)
    return ((1 + sq) / 2)**(2 * m) + ((1 - sq) / 2)**(2 * m)

def eval_x_root_mp(t, p):
    # Find the root x near 1 of T^3 - (1+t) T^2 + t T - p = 0 by Newton!
    r = mp.mpc(1)
    s1 = 1 + t
    for _ in range(10):
        fv = ((r - s1) * r + t) * r - p
        dfv = (3 * r - 2 * s1) * r + t
        r -= fv / dfv
    return r

def test_lucas_vs_exact(n, c, j, N=48):
    m = n // 2
    r1 = mp.sin((j + c) * mp.pi / n)
    r2 = mp.sin(j * mp.pi / n)
    r3 = mp.sin(c * mp.pi / n)
    s = r2 / r1
    Ry = s ** n
    Rz = (r3 / r1) ** n
    t = Ry if (c % 2 == 0) else -Ry
    log_denom_r = c * (mp.log(Ry) + mp.log(Rz))

    half = N // 2
    ex_lvs = []; fr_lvs = []
    max_x_dev = mp.mpf(0)
    for idx in range(half):
        phi = mp.pi * (idx + mp.mpf('0.5')) / half
        p = t * Rz * mp.exp(1j * phi)
        xv = eval_x_root_mp(t, p)
        if abs(xv - 1) > max_x_dev: max_x_dev = abs(xv - 1)

        # EXACT via Lucas-Chebyshev product!
        p_ex = p / (xv ** 3)
        s_ex = (1 + t - xv) / xv
        rp_ex = abs(p_ex) ** (mp.mpf(1) / n); ap_ex = mp.arg(p_ex) / n
        sl_ex = n * mp.log(xv)
        for k in range(n):
            wk = rp_ex * mp.exp(1j * (ap_ex + 2 * mp.pi * k / n))
            sl_ex += mp.log(L2m_mp(m, wk) - s_ex)
        lv_ex = sl_ex - log_denom_r - 1j * c * phi
        if c % 2 == 1: lv_ex += 1j * mp.pi
        if j % 2 == 1: lv_ex += 1j * mp.pi

        # FROZEN x = 1 Lucas-Chebyshev product!
        rp_fr = abs(p) ** (mp.mpf(1) / n); ap_fr = mp.arg(p) / n
        sl_fr = mp.mpc(0)
        for k in range(n):
            wk = rp_fr * mp.exp(1j * (ap_fr + 2 * mp.pi * k / n))
            sl_fr += mp.log(L2m_mp(m, wk) - t)
        lv_fr = sl_fr - log_denom_r - 1j * c * phi
        if c % 2 == 1: lv_fr += 1j * mp.pi
        if j % 2 == 1: lv_fr += 1j * mp.pi

        ex_lvs.append(lv_ex); fr_lvs.append(lv_fr)

    sc = max(lv.real for lv in ex_lvs)
    ev = [mp.exp(lv - sc) for lv in ex_lvs]
    fv = [mp.exp(lv - sc) for lv in fr_lvs]
    avg_ex = sum(v.real for v in ev) / half
    avg_fr = sum(v.real for v in fv) / half
    max_diff = max(abs(ev[i] - fv[i]) for i in range(half)) / avg_ex
    avg_diff = abs(avg_ex - avg_fr) / avg_ex
    min_re_fr = min(v.real for v in fv) / avg_fr
    return float(max_x_dev), float(max_diff), float(avg_diff), float(min_re_fr)

print("=== EXACT vs FROZEN-x=1 LUCAS-CHEBYSHEV 1D PRODUCT ===")
print(f"{'case':15s} {'n':>4s} {'c':>3s} {'j':>3s} | {'max|x-1|':>10s} {'max|E-L|/avg':>12s} {'|avgE-avgL|/avg':>15s} {'minRe(L)/avg':>12s}")
for tag, n, c, j in [
    ("c=2,j=2",     40,  2,  2), ("c=2,j=2",     80,  2,  2), ("c=2,j=2",    160,  2,  2),
    ("c=3,j=3",     40,  3,  3), ("c=3,j=3",     80,  3,  3),
    ("bulk",        80, 13, 20),
    ("top_left",    40, 12, 12), ("top_left",    80, 26, 26),
    ("top_mid",     80, 25, 26),
]:
    xd, md, ad, mrf = test_lucas_vs_exact(n, c, j, N=48)
    print(f"{tag:15s} {n:4d} {c:3d} {j:3d} | {xd:10.2e} {md:12.3e} {ad:15.3e} {mrf:+12.4f}")
