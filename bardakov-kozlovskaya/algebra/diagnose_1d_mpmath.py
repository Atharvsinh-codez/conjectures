"""
Exact high-precision (50-dps mpmath) diagnostic of the 1D Circle Integral
  I_{m,c,j}(phi; Rv) = (-1)^{j+c} P_{2m}(1 + t_j, t_j, t_j * Rz * e^{i phi}) / (|t_j| * Rz * e^{i phi})^c
at ALL regimes for n = 20, 40, 80, 200, 500:
  - Fixed small c = 2, 3, 4 at j = c (left corner) and j = m/2 (mid) and j = j_end (right corner)!
  - Bulk (c ~ n/6, j ~ n/4)
  - Top corner (c = floor((n-2)/3), j = c and j = j_end)
Check BOTH the Maillot radius Rz_M = (r3/r1)^n AND scaled radii alpha * Rz_M (e.g. alpha = 0.25, 0.5, 1, 2, 4)!
"""
import math
import mpmath as mp

mp.mp.dps = 50

def eval_cubic_roots_mp(t, e3_val):
    # Solve T^3 - (1+t) T^2 + t T - e3_val = 0 cleanly in 50-dps mpmath via exact Cardano + 2 Newton polishes!
    # Depress: T = W + (1+t)/3 -> W^3 + p W + q = 0
    s1 = 1 + t
    p = t - s1 * s1 / 3
    q = -2 * s1**3 / 27 + s1 * t / 3 - e3_val
    disc = (q / 2)**2 + (p / 3)**3
    sd = mp.sqrt(disc)
    u3 = -q / 2 + sd
    v3 = -q / 2 - sd
    # Pick larger cube root to avoid cancellation
    if abs(v3) > abs(u3):
        u3, v3 = v3, u3
    if u3 == 0:
        r0 = s1 / 3
        rts = [r0, r0, r0]
    else:
        u = u3 ** (mp.mpf(1) / 3)
        om = mp.exp(2j * mp.pi / 3)
        rts = []
        for k in range(3):
            uk = u * (om ** k)
            vk = -p / (3 * uk) if uk != 0 else mp.mpc(0)
            r = uk + vk + s1 / 3
            # 2 Newton steps for full 50-dps accuracy even when e3 ~ 1e-200!
            for _ in range(3):
                fv = ((r - s1) * r + t) * r - e3_val
                dfv = (3 * r - 2 * s1) * r + t
                if dfv != 0:
                    r -= fv / dfv
            rts.append(r)
    rts.sort(key=lambda w: abs(w), reverse=True)
    return rts[0], rts[1], rts[2]

def eval_1d_mp(n, c, j, alpha=1.0, N=64):
    m = n // 2
    r1 = mp.sin((j + c) * mp.pi / n)
    r2 = mp.sin(j * mp.pi / n)
    r3 = mp.sin(c * mp.pi / n)
    s = r2 / r1
    Ry = s ** n
    Rz = mp.mpf(alpha) * ((r3 / r1) ** n)
    t = Ry if (c % 2 == 0) else -Ry
    log_denom_r = c * (mp.log(Ry) + mp.log(Rz))

    # Symmetry: integrand(-phi) = conj(integrand(phi)), so only compute phi in [0, pi] (N/2+1 pts) with trapezoidal weights!
    half = N // 2
    lvs = []
    weights = []
    for idx in range(half + 1):
        phi = mp.pi * idx / half
        w_tr = mp.mpf('0.5') if (idx == 0 or idx == half) else mp.mpf(1)
        e3_val = t * Rz * mp.exp(1j * phi)
        xv, yv, zv = eval_cubic_roots_mp(t, e3_val)
        yr = yv / xv
        zr = zv / xv
        ry = abs(yr) ** (mp.mpf(1) / n)
        ay = mp.arg(yr) / n
        # sum_{k=0}^{n-1} log(zr - (1 + ry e^{i(ay + 2pi k/n)})^n)
        s_log = mp.mpc(0)
        for k in range(n):
            wk = 1 + ry * mp.exp(1j * (ay + 2 * mp.pi * k / n))
            fac = zr - (wk ** n)
            if abs(fac) == 0:
                s_log = mp.mpc(-1e9)
                break
            s_log += mp.log(fac)
        lv = n * mp.log(xv) + s_log - log_denom_r - 1j * c * phi
        if c % 2 == 1: lv += 1j * mp.pi
        if j % 2 == 1: lv += 1j * mp.pi
        lvs.append(lv)
        weights.append(w_tr)

    max_re = max(lv.real for lv in lvs)
    F = [mp.exp(lv - max_re) for lv in lvs]
    tot_w = sum(weights)
    avg_re = sum(weights[i] * F[i].real for i in range(len(F))) / tot_w
    abs_avg = sum(weights[i] * abs(F[i]) for i in range(len(F))) / tot_w
    min_re = min(f.real for f in F)
    return float(avg_re / abs_avg), float(min_re / avg_re)

print("=== 50-DPS MPMATH 1D CIRCLE INTEGRAL STRESS TEST (n = 20, 40, 80, 200) ===")
for tag, fn in [
    ("c=2, j=2 (left bot)",       lambda n: (2, 2)),
    ("c=2, j=m//2 (mid bot)",     lambda n: (2, n//4)),
    ("c=2, j=m-1 (u=1/4 bot)",    lambda n: (2, n//2 - 1)),
    ("c=3, j=3 (left bot)",       lambda n: (3, 3)),
    ("c=3, j=m-2 (right bot)",    lambda n: (3, n//2 - 2)),
    ("c=m//3, j=m//2 (bulk)",     lambda n: (max(2, n//6), n//4)),
    ("c=max_c, j=c (top left)",   lambda n: ((n-2)//3, (n-2)//3)),
    ("c=max_c, j=j_end (top rt)", lambda n: ((n-2)//3, n//2 - ((n-2)//3 + 1)//2)),
]:
    print(f"\n--- {tag} ---")
    for n in [20, 40, 80, 200]:
        c, j = fn(n)
        r1, mr1 = eval_1d_mp(n, c, j, alpha=1.0, N=48)
        extra = ""
        if r1 < 0.90:
            # Also try alpha = 0.5, 2.0, 4.0 to see if saddle radius shift helps!
            best = (r1, 1.0, mr1)
            for al in [0.25, 0.5, 2.0, 4.0]:
                ra, mra = eval_1d_mp(n, c, j, alpha=al, N=48)
                if ra > best[0]: best = (ra, al, mra)
            extra = f" | best alpha={best[1]}: ratio={best[0]:+.5f}, minRe/avg={best[2]:+.4f}"
        print(f"  n={n:3d} (c={c:2d}, j={j:2d}): alpha=1 -> Re(avg)/avg|F|={r1:+.6f}, minRe/avg={mr1:+.4f}{extra}")
