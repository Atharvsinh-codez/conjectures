"""
Deep quantitative analysis of the 1D Circle Integral error budget across ALL regimes of (c, j):
  I_{m, c, j}(phi) = (-1)^{j+c} P_{2m}(1+t_j, t_j, t_j R_z e^{i phi}) / (|t_j| R_z e^{i phi})^c.

Check:
1. True positivity ratio: int_0^pi Re(F) / int_0^pi |Re(F)| (NOT avg(|F|), since Im(F) integrates to 0!)
   at c=2, j=2 and c=3, j=3 for n = 40, 80, 160, 240, 400!
2. EXACT decomposition of F(phi) = F_frozen(phi) * E_cubic(phi), where:
   - F_frozen(phi) = (-1)^{j+c} prod_{k=0}^{2m-1} (R_z e^{i phi} - (1 + xi_k s)^{2m}) / (R_z e^{i phi})^c
     (NOTE: its exact integral is (-1)^{j+c} [z^c] p_{2m}(1, t_j, z) = (-1)^j E_{2m-c}(1, t_j)!!)
   - E_cubic(phi) = F(phi) / F_frozen(phi) = cubic-root perturbation factor!
   WAIT A SECOND!! CHECK |E_cubic(phi) - 1| ACROSS ALL REGIMES!!
   And check the EXACT bound on the off-boundary factors of F_frozen(phi):
     F_frozen(phi) = C_{m,c,j} * 4 sin^2((phi + j pi)/2) * prod_{small} (1 - w_k / (R_z e^{i phi})) * prod_{large} (1 - R_z e^{i phi} / w_k)!
   CHECK sum_{small} |w_k / R_z| + sum_{large} |R_z / w_k| across ALL regimes!
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

def analyze_case(n, c, j, N=48):
    r1 = mp.sin((j + c) * mp.pi / n)
    r2 = mp.sin(j * mp.pi / n)
    r3 = mp.sin(c * mp.pi / n)
    s = r2 / r1
    Ry = s ** n
    Rz = (r3 / r1) ** n
    t = Ry if (c % 2 == 0) else -Ry
    log_denom_r = c * (mp.log(Ry) + mp.log(Rz))

    # Frozen factors w_k = (1 - s e^{i eta_k})^n
    # Even c: eta_k = 2*pi*k/n (k=0..n-1). Odd c: eta_k = (2k+1)*pi/n (k=0..n-1).
    # Compute sum_{small} |w_k|/Rz + sum_{large} Rz/|w_k|!
    sum_small = mp.mpf(0); sum_large = mp.mpf(0)
    w_frozen = []
    for k in range(n):
        eta = (2*k * mp.pi / n) if (c % 2 == 0) else ((2*k + 1) * mp.pi / n)
        # normalize eta to (-pi, pi]
        if eta > mp.pi: eta -= 2 * mp.pi
        wk = (1 - s * mp.exp(1j * eta)) ** n
        w_frozen.append((eta, wk))
        ae = abs(eta)
        th_c = c * mp.pi / n
        if abs(ae - th_c) < mp.mpf('1e-30'):
            continue  # boundary pair!
        elif ae < th_c:
            sum_small += abs(wk) / Rz
        else:
            sum_large += Rz / abs(wk)

    half = N // 2
    lvs = []; lvs_frz = []; ws = []
    max_cubic_dev = mp.mpf(0)
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

        # Frozen log value (exact same normalization):
        # P_frozen(phi) = prod_{k=0}^{n-1} (Rz e^{i phi} - w_k), divided by (Rz e^{i phi})^c!
        # Wait: in exact P_{2m}, e3 = t Rz e^{i phi}, so z ~ Rz e^{i phi}, and we divided by |t|^c Rz^c e^{ic phi} * (-1)^c;
        # wait: why was t^c in denominator? Because at frozen roots, P_{2m}(1, t, Rz e^{i phi}) is NOT multiplied by t^c!
        # Wait: in our exact formula, P_{2m}(1+t, t, t Rz e^{i phi}) has roots x~1, y~t, z~Rz e^{i phi},
        # so P_{2m} ~ prod (Rz e^{i phi} - (1+xi_k s)^n), and dividing by (t Rz e^{i phi})^c has an extra t^{-c} = Ry^{-c} (-1)^c!
        s_frz = mp.mpc(0)
        for eta, wk in w_frozen:
            fac = Rz * mp.exp(1j * phi) - wk
            if abs(fac) == 0:
                s_frz = mp.mpc(-1e9); break
            s_frz += mp.log(fac)
        lv_frz = s_frz - log_denom_r - 1j * c * phi
        if c % 2 == 1: lv_frz += 1j * mp.pi
        if j % 2 == 1: lv_frz += 1j * mp.pi
        lvs_frz.append(lv_frz)
        if s_log.real > -1e8 and s_frz.real > -1e8:
            dev = abs(mp.exp(lv - lv_frz) - 1)
            if dev > max_cubic_dev: max_cubic_dev = dev

    max_re = max(lv.real for lv in lvs)
    F = [mp.exp(lv - max_re) for lv in lvs]
    tot_w = sum(ws)
    avg_re = sum(ws[i] * F[i].real for i in range(len(F))) / tot_w
    abs_re_avg = sum(ws[i] * abs(F[i].real) for i in range(len(F))) / tot_w
    min_re = min(f.real for f in F)
    return float(avg_re / abs_re_avg), float(min_re / avg_re), float(sum_small), float(sum_large), float(max_cubic_dev)

print("=== ERROR BUDGET ANALYSIS OF 1D CIRCLE INTEGRAL ===")
print(f"{'case':16s} {'n':>4s} {'c':>3s} {'j':>3s} | {'Re/absRe':>9s} {'minRe/avg':>9s} | {'sum_sm':>9s} {'sum_lg':>9s} {'max|Ec-1|':>10s}")
for tag, n, c, j in [
    ("c=2,j=2", 40, 2, 2), ("c=2,j=2", 80, 2, 2), ("c=2,j=2", 160, 2, 2), ("c=2,j=2", 320, 2, 2),
    ("c=3,j=3", 40, 3, 3), ("c=3,j=3", 80, 3, 3), ("c=3,j=3", 160, 3, 3),
    ("c=2,j=m/2", 40, 2, 10), ("c=2,j=m/2", 80, 2, 20),
    ("c=2,u=1/4", 40, 2, 19), ("c=2,u=1/4", 80, 2, 39),
    ("bulk", 40, 6, 10), ("bulk", 80, 13, 20),
    ("top_left", 40, 12, 12), ("top_left", 80, 26, 26),
    ("top_rt", 40, 12, 14), ("top_rt", 80, 26, 27),
]:
    r_re, mr, s_sm, s_lg, ec = analyze_case(n, c, j, N=48)
    print(f"{tag:16s} {n:4d} {c:3d} {j:3d} | {r_re:+9.6f} {mr:+9.4f} | {s_sm:9.3e} {s_lg:9.3e} {ec:10.3e}")
