"""
Check the TRUE error budget of the 1D Circle Integral:
1. Measure max_phi |F_exact(phi) - F_frozen(phi)| / avg(F_exact) across ALL regimes (fixing the 0/0 bug!).
2. Specifically inspect TOP CORNER (c=26, j=26,27 at n=80) to see WHY F_exact has minRe/avg = +0.72 and +0.23!
3. Inspect BOTTOM-LEFT CORNER (c=2,3,4, j=c..c+3 at n=80, 160, 320) to see the EXACT mechanism keeping Re/absRe ~ 1!
"""
import mpmath as mp

mp.mp.dps = 120

def eval_cubic_roots_mp(t, e3_val):
    s1 = 1 + t
    if abs(e3_val) < mp.mpf('1e-8') * abs(t) * max(mp.mpf('1e-12'), abs(1 - t)**2):
        rts = [mp.mpc(1), mp.mpc(t), e3_val / t]
        for idx in range(3):
            r = rts[idx]
            for _ in range(6):
                fv = ((r - s1) * r + t) * r - e3_val
                dfv = (3 * r - 2 * s1) * r + t
                if dfv != 0: r -= fv / dfv
            rts[idx] = r
    else:
        p = t - s1 * s1 / 3
        q = -2 * s1**3 / 27 + s1 * t / 3 - e3_val
        disc = (q / 2)**2 + (p / 3)**3
        sd = mp.sqrt(disc)
        u3 = -q / 2 + sd; v3 = -q / 2 - sd
        if abs(v3) > abs(u3): u3, v3 = v3, u3
        u = u3 ** (mp.mpf(1) / 3)
        om = mp.exp(2j * mp.pi / 3)
        rts = []
        for k in range(3):
            uk = u * (om ** k)
            vk = -p / (3 * uk) if uk != 0 else mp.mpc(0)
            r = uk + vk + s1 / 3
            for _ in range(5):
                fv = ((r - s1) * r + t) * r - e3_val
                dfv = (3 * r - 2 * s1) * r + t
                if dfv != 0: r -= fv / dfv
            rts.append(r)
    rts.sort(key=lambda w: abs(w), reverse=True)
    return rts[0], rts[1], rts[2]

def compare_exact_vs_frozen(n, c, j, N=48):
    r1 = mp.sin((j + c) * mp.pi / n)
    r2 = mp.sin(j * mp.pi / n)
    r3 = mp.sin(c * mp.pi / n)
    s = r2 / r1
    Ry = s ** n
    Rz = (r3 / r1) ** n
    t = Ry if (c % 2 == 0) else -Ry
    log_denom_r = c * (mp.log(Ry) + mp.log(Rz))

    w_frozen = []
    for k in range(n):
        eta = (2*k * mp.pi / n) if (c % 2 == 0) else ((2*k + 1) * mp.pi / n)
        if eta > mp.pi: eta -= 2 * mp.pi
        w_frozen.append((1 - s * mp.exp(1j * eta)) ** n)

    half = N // 2
    # Evaluate both exact and frozen at phi = pi*(idx + 0.5)/half (midpoint rule avoids exact 0!)
    F_ex = []; F_fr = []
    for idx in range(half):
        phi = mp.pi * (idx + mp.mpf('0.5')) / half
        e3_val = t * Rz * mp.exp(1j * phi)
        xv, yv, zv = eval_cubic_roots_mp(t, e3_val)
        yr = yv / xv; zr = zv / xv
        ry = abs(yr) ** (mp.mpf(1) / n); ay = mp.arg(yr) / n
        s_log = mp.mpc(0)
        for k in range(n):
            wk = 1 + ry * mp.exp(1j * (ay + 2 * mp.pi * k / n))
            s_log += mp.log(zr - (wk ** n))
        lv = n * mp.log(xv) + s_log - log_denom_r - 1j * c * phi
        if c % 2 == 1: lv += 1j * mp.pi
        if j % 2 == 1: lv += 1j * mp.pi

        s_frz = mp.mpc(0)
        for wk in w_frozen:
            s_frz += mp.log(Rz * mp.exp(1j * phi) - wk)
        lv_frz = s_frz - log_denom_r - 1j * c * phi
        if c % 2 == 1: lv_frz += 1j * mp.pi
        if j % 2 == 1: lv_frz += 1j * mp.pi

        F_ex.append(lv); F_fr.append(lv_frz)

    scale = max(lv.real for lv in F_ex)
    ex_v = [mp.exp(lv - scale) for lv in F_ex]
    fr_v = [mp.exp(lv - scale) for lv in F_fr]
    avg_ex = sum(v.real for v in ex_v) / half
    avg_fr = sum(v.real for v in fr_v) / half
    abs_re_ex = sum(abs(v.real) for v in ex_v) / half
    max_diff = max(abs(ex_v[i] - fr_v[i]) for i in range(half)) / avg_ex
    avg_diff = abs(avg_ex - avg_fr) / avg_ex
    # Also check how much of the integral comes from the frozen leading 4*sin^2((phi+j*pi)/2)
    return float(avg_ex / abs_re_ex), float(max_diff), float(avg_diff), float(Rz), float(Ry)

print("=== EXACT vs FROZEN COMPARISON (NO 0/0 BUG!) ===")
print(f"{'case':15s} {'n':>4s} {'c':>3s} {'j':>3s} | {'Re/absRe':>9s} {'max|E-F|/avg':>12s} {'|avgE-avgF|/avg':>15s} | {'Rz':>9s} {'Ry':>9s}")
for tag, n, c, j in [
    ("c=2,j=2",     40,  2,  2), ("c=2,j=2",     80,  2,  2), ("c=2,j=2",    160,  2,  2),
    ("c=2,j=3",     80,  2,  3), ("c=2,j=4",     80,  2,  4), ("c=2,j=5",     80,  2,  5),
    ("c=3,j=3",     80,  3,  3), ("c=3,j=4",     80,  3,  4), ("c=3,j=5",     80,  3,  5),
    ("c=2,j=m/2",   80,  2, 20), ("c=2,u=1/4",   80,  2, 39),
    ("bulk",        80, 13, 20), ("bulk_rt",     80, 14, 33),
    ("top_left",    40, 12, 12), ("top_left",    80, 26, 26),
    ("top_mid",     80, 25, 26), ("top_rt",      80, 26, 27),
]:
    r_re, md, ad, Rz, Ry = compare_exact_vs_frozen(n, c, j, N=48)
    print(f"{tag:15s} {n:4d} {c:3d} {j:3d} | {r_re:+9.6f} {md:12.3e} {ad:15.3e} | {Rz:9.2e} {Ry:9.2e}")
