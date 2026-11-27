"""
3-second 250-dps mpmath check of the exact n=80, 120, 200, 300 behavior of the 1D circle integral
at the ONLY places 50-dps underflowed (c=2, j=2; c=3, j=3; bulk c=m//3, j=m//2)!
Also test the SIMPLER frozen-root approximation (x=1, y=t, z=Rz e^{i phi}) vs EXACT cubic roots!
"""
import mpmath as mp

mp.mp.dps = 250

def eval_cubic_roots_mp(t, e3_val):
    # For |e3_val| << |t|*|1-t|, Newton from (1, t, e3_val/t) converges to 250 dps in 4 iterations!
    # For general e3_val, use Cardano + 5 Newton steps!
    s1 = 1 + t
    if abs(e3_val) < mp.mpf('1e-10') * abs(t) * max(mp.mpf('1e-10'), abs(1 - t)**2):
        r0 = mp.mpc(1)
        r1 = mp.mpc(t)
        r2 = e3_val / t
        rts = [r0, r1, r2]
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
        u3 = -q / 2 + sd
        v3 = -q / 2 - sd
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

def eval_1d_250(n, c, j, N=48):
    r1 = mp.sin((j + c) * mp.pi / n)
    r2 = mp.sin(j * mp.pi / n)
    r3 = mp.sin(c * mp.pi / n)
    s = r2 / r1
    Ry = s ** n
    Rz = (r3 / r1) ** n
    t = Ry if (c % 2 == 0) else -Ry
    log_denom_r = c * (mp.log(Ry) + mp.log(Rz))
    half = N // 2
    lvs = []; ws = []
    for idx in range(half + 1):
        phi = mp.pi * idx / half
        w_tr = mp.mpf('0.5') if (idx == 0 or idx == half) else mp.mpf(1)
        e3_val = t * Rz * mp.exp(1j * phi)
        xv, yv, zv = eval_cubic_roots_mp(t, e3_val)
        yr = yv / xv; zr = zv / xv
        ry = abs(yr) ** (mp.mpf(1) / n); ay = mp.arg(yr) / n
        s_log = mp.mpc(0)
        for k in range(n):
            wk = 1 + ry * mp.exp(1j * (ay + 2 * mp.pi * k / n))
            fac = zr - (wk ** n)
            if abs(fac) == 0:
                s_log = mp.mpc(-1e9); break
            s_log += mp.log(fac)
        lv = n * mp.log(xv) + s_log - log_denom_r - 1j * c * phi
        if c % 2 == 1: lv += 1j * mp.pi
        if j % 2 == 1: lv += 1j * mp.pi
        lvs.append(lv); ws.append(w_tr)
    max_re = max(lv.real for lv in lvs)
    F = [mp.exp(lv - max_re) for lv in lvs]
    tot_w = sum(ws)
    avg_re = sum(ws[i] * F[i].real for i in range(len(F))) / tot_w
    abs_avg = sum(ws[i] * abs(F[i]) for i in range(len(F))) / tot_w
    min_re = min(f.real for f in F)
    return float(avg_re / abs_avg), float(min_re / avg_re)

for tag, cases in [
    ("c=2, j=2", [(80, 2, 2), (160, 2, 2), (240, 2, 2)]),
    ("c=3, j=3", [(80, 3, 3), (160, 3, 3), (240, 3, 3)]),
    ("c=2, j=3", [(80, 2, 3), (160, 2, 3), (240, 2, 3)]),
    ("bulk",     [(80, 13, 20), (160, 26, 40), (200, 33, 50)]),
]:
    print(f"--- {tag} ---")
    for (n, c, j) in cases:
        r, mr = eval_1d_250(n, c, j, N=48)
        print(f"  n={n:3d} (c={c:2d}, j={j:2d}): Re(avg)/avg|F| = {r:+.6f}, minRe/avg = {mr:+.4f}")
