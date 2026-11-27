"""Deep check #2 for Piece G3 (Cases ABC & D): y>0 outer poles up to L_n, d_m/mlog fallbacks, straddle bin, d_lam_th_ub, 0.39*Del0."""
import math
import cmath
import sys
import numpy as np
import mpmath as mp
from mpmath import iv

mp.mp.dps = 50
iv.dps = 50
PI = iv.pi

sys.path.insert(0, r"d:\solving-conjectures\work\bardakov-kozlovskaya\analyst")
import certify_G3_D_full_iv as g3d

def make_point(c, j, r0, target_kap):
    n_est = int(round(2.0 * math.pi**2 * j * (j + c) / (c * target_kap)))
    n = max(132, 2 * ((n_est + 1) // 2))
    th_c = c * math.pi / n
    psi_j = j * math.pi / n
    al_1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(al_1), math.sin(psi_j), math.sin(th_c)
    sj = r2 / r1
    kap = 2.0 * math.pi * r1 * r2 / r3
    tau = c * r0 * kap
    p0 = round(r0 * j / c)
    D = r0 * j - c * p0
    ct = c / r0
    return dict(n=n, c=c, j=j, r0=r0, p0=p0, D=D, ct=ct, kap=kap, tau=tau,
                th_c=th_c, psi_j=psi_j, al_1=al_1, r1=r1, r2=r2, r3=r3, sj=sj)

def zeta_fn(pt, u):
    w = pt['sj'] * cmath.exp(2j * u * pt['th_c'])
    return pt['r0'] * pt['n'] * cmath.log((pt['r3'] / pt['r1']) / (1.0 - w)) + 1j * math.pi * pt['r0'] * pt['j'] - 2j * math.pi * pt['c'] * pt['p0'] * (u + 0.5)

def T_fn(pt, u):
    return 4.0 * pt['sj'] / ((1.0 - pt['sj'])**2) * (math.sin(u * pt['th_c'])**2)

print("=== CHECK A: y > 0 outer poles across FULL x in [u0, L_n] (both D=+1 and D=-1, up to y = 0.50*r0) ===")
for c, j, r0, kap in [
    (12, 13, 1, 0.44),     # D=+1, high kap, small ct
    (12, 23, 1, 0.44),     # D=-1, high kap, small ct
    (20, 21, 1, 0.44),     # D=+1, tau=8.8
    (30, 31, 1, 0.44),     # D=+1, tau=13.2
    (89, 90, 1, 0.449),    # D=+1, tau=40
    (40, 42, 1, 0.40),     # D=+2, tau=16
    (20, 21, 2, 0.22),     # r0=2, D=+2? wait 2*21 - 20*2 = +2
    (21, 32, 2, 0.22),     # r0=2, 2*32 - 21*3 = +1
]:
    pt = make_point(c, j, r0, kap)
    tau, ct, D = pt['tau'], pt['ct'], pt['D']
    u0 = 0.48
    re0 = zeta_fn(pt, u0).real
    L_n = pt['n'] / (2.0 * c)
    G = g3d.box_geom(iv.mpf(str(pt['kap'])))
    Del0 = G['C_D0'] * iv.mpf(str(r0 * c * pt['kap']**2))
    d_lo = abs(D) - Del0
    E_out_0 = float(((4 * PI**2 * d_lo / G['te_str']) * (iv.mpf(ct) / iv.mpf(str(tau)))).a)
    # Sign of y on the "same-side-as-inner" outer branch is sgn(D)
    sgn = 1 if D > 0 else -1
    # Check where Re zeta(x + i*sgn*y) = re0 for x in (u0, L_n], y in (0, 0.50*r0]
    min_E = 1e9; best = None
    # Sample logarithmically from u0 to L_n
    xs = np.concatenate([np.linspace(u0 + 0.01, min(L_n, 2.0), 150),
                         np.geomspace(min(L_n, 2.0), L_n, 250)])
    for x in xs:
        # At y=0, Re zeta(x) - re0 < 0. Check max over y in (0, 0.50*r0]
        ys = np.linspace(0.0, 0.50 * r0, 120)
        vals = [zeta_fn(pt, x + 1j * sgn * y).real - re0 for y in ys]
        if max(vals) >= 0:
            for idx in range(len(ys) - 1):
                if vals[idx] < 0 and vals[idx + 1] >= 0:
                    ya, yb = ys[idx], ys[idx + 1]
                    for _ in range(30):
                        ym = 0.5 * (ya + yb)
                        if zeta_fn(pt, x + 1j * sgn * ym).real - re0 < 0: ya = ym
                        else: yb = ym
                    ym = 0.5 * (ya + yb)
                    Ev = 2 * math.pi * ct * ym
                    if Ev < min_E:
                        min_E = Ev; best = (x, ym, T_fn(pt, x))
                    break
    print(f"  c={c:2d} j={j:2d} r0={r0} D={D:+d} tau={tau:5.2f} L_n={L_n:6.1f}: E_out_0(code)={E_out_0:6.2f} | min exact E(x>u0, sgn*y>0)={min_E:9.2f} at {best}")

print("\n=== CHECK B: Why did y>0 outer branch stay > E_out_0 even when T(x)>0? Check Re zeta(x+iy) - Re zeta(x)! ===")
pt = make_point(89, 90, 1, 0.449)
for x in [0.5, 1.0, 3.0, 7.0, 15.0]:
    # Check d/dy Re zeta(x+iy) at y=0 and Re zeta(x + 0.5i) - Re zeta(u0)
    re0 = zeta_fn(pt, 0.48).real
    rex = zeta_fn(pt, x).real
    d_re_0 = (zeta_fn(pt, x + 1e-6j).real - rex) / 1e-6
    re_top = zeta_fn(pt, x + 0.5j).real
    print(f"  x={x:4.1f} T(x)={T_fn(pt, x):.4f}: Re(x)-Re(u0)={rex-re0:+8.2f}, d/dy|_{{0}}={d_re_0:+7.2f}, Re(x+0.5i)-Re(u0)={re_top-re0:+8.2f}")

print("\n=== CHECK C: Does d_m <= 0.05 or x_f_lo <= 1e-6 ACTUALLY trigger in _eval_g1_fin_single_iv and eval_Psi_diff_W1_iv? ===")
# Instrument by inspecting bins for representative boxes across D1, D2, D3
for lbl, ta_f, tb_f, u0_a, u0_b, ct_lo_f, ct_hi_f, th_a, th_b in [
    ("D1-T1 peak",   1.05, 1.25, 0.42, 0.42, 25.6, 48.0, 0.0, 0.083),
    ("D1-T3 valley", 1.05, 1.25, 0.26, 0.26, 80.0, 115.0, 1.5, 1.8),
    ("D2-T5 peak",   1.65, 1.82, 0.185, 0.185, 185.0, 360.0, 0.0, 0.04),
    ("D2-B valley",  1.65, 1.82, 0.16, 0.185, 360.0, None, 1.5, 1.8),
    ("D3-A1 peak",   2.80, 2.98, 0.42, 0.42, 10.0, 14.0, 0.0, 0.1),
    ("D3-B valley",  9.50, 11.2, 0.495, 0.510, 55.0, None, 2.8, 3.14),
]:
    ta = iv.mpf(str(ta_f)); tb = iv.mpf(str(tb_f))
    u0_iv = iv.mpf([iv.mpf(str(u0_a)).a, iv.mpf(str(u0_b)).b])
    ct_lo = iv.mpf(str(ct_lo_f)); ct_hi = None if ct_hi_f is None else iv.mpf(str(ct_hi_f))
    G, c_min_iv, r0ck2_hi, Del0_hi = g3d.make_params(ta, tb, ct_lo, ct_hi)
    u0_lo = iv.mpf(u0_iv.a); u0_hi = iv.mpf(u0_iv.b)
    # Check _eval_g1_fin_single_iv bins for d_m <= 0.05 or x_f_lo <= 1e-6:
    L_val = float(iv.sqrt(u0_hi ** 2 + iv.mpf('11.0') / (G['c_re5'] * ta)).b)
    L_iv = iv.mpf(str(round(L_val, 2)))
    u_lo_f = float(u0_lo.a); u_hi_f = float(u0_hi.b); L_f = float(L_iv.b)
    u_act = min(L_f, math.sqrt(u_hi_f * u_hi_f + 4.5 / max(float(ta.a), 0.5)))
    n_bins = 56
    pts = [u_lo_f * k / float(n_bins) for k in range(n_bins + 1)]
    if u_hi_f > u_lo_f + 1e-10: pts.append(u_hi_f)
    for k in range(1, n_bins + 1): pts.append(u_hi_f + (u_act - u_hi_f) * ((k / float(n_bins)) ** 1.25))
    for k in range(1, 9): pts.append(u_act + (L_f - u_act) * (k / 8.0))
    hits_dm = []; hits_xlo = []
    for i in range(len(pts) - 1):
        ua = iv.mpf(str(pts[i])); ub = iv.mpf(str(pts[i + 1]))
        if (ub - ua).b <= mp.mpf('1e-14'): continue
        if ub.b <= u0_lo.a + mp.mpf('1e-12'):
            x_f_lo = G['c_rein'] * ta * (u0_lo ** 2 - ub ** 2)
            r_f_ub = iv.exp(-x_f_lo)
            dz_e = g3d.dzeta_err_bound(G, ua, ub, u0_lo, u0_hi, tb, Del0_hi, r0ck2_hi, True)
            ph_lim = g3d.iv_max_hi(abs(PI - 2 * PI * ua), abs(PI - 2 * PI * ub)) + dz_e
            d_ang = iv.mpf(1) if ph_lim.b <= (PI / 2).a else (iv.sin(PI - ph_lim) if ph_lim.b < PI.a - mp.mpf('0.03') else iv.mpf(0))
            d_m = iv.mpf(max((1 - r_f_ub).a, d_ang.a))
            if d_m.a <= mp.mpf('0.05'): hits_dm.append(('in', float(ua.a), float(ub.b), float(d_m.a), float(x_f_lo.a)))
            if x_f_lo.a <= mp.mpf('1e-6'): hits_xlo.append(('in', float(ua.a), float(ub.b), float(d_m.a)))
        elif ua.a >= u0_hi.b - mp.mpf('1e-12'):
            x_f_lo = G['c_re5'] * ta * (ua ** 2 - u0_hi ** 2)
            r_f_ub = iv.exp(-x_f_lo)
            dz_e = g3d.dzeta_err_bound(G, ua, ub, u0_lo, u0_hi, tb, Del0_hi, r0ck2_hi, False)
            d_ang = iv.mpf(0)
            if ub.b <= mp.mpf('0.90'):
                ph_lim = g3d.iv_max_hi(abs(2 * PI * ua - PI), abs(2 * PI * ub - PI)) + dz_e
                d_ang = iv.mpf(1) if ph_lim.b <= (PI / 2).a else (iv.sin(PI - ph_lim) if ph_lim.b < PI.a - mp.mpf('0.03') else iv.mpf(0))
            d_m = iv.mpf(max((1 - r_f_ub).a, d_ang.a))
            if d_m.a <= mp.mpf('0.05'): hits_dm.append(('out', float(ua.a), float(ub.b), float(d_m.a), float(x_f_lo.a)))
            if x_f_lo.a <= mp.mpf('1e-6'): hits_xlo.append(('out', float(ua.a), float(ub.b), float(d_m.a)))
    print(f"  {lbl:14s}: d_m<=0.05 hits={hits_dm[:3]} (total {len(hits_dm)}) | x_f_lo<=1e-6 hits={len(hits_xlo)} (used in W1 valley mlog_hi=14.0!)")

print("\n=== CHECK D: Exact max |h(u)| / Del0 on [0, 0.515] vs hardcoded 0.39 ===")
G0 = g3d.box_geom(iv.mpf('0.45'))
# Exact Delta_2(t) = A * sin^2(t*th_c) / (1 + T(t)), with T(t) <= gam * t^2 <= 0.0028 on [0, 0.515], th_c <= 0.036
# On [0, 0.5]: sin^2(t*th_c)/(t*th_c)^2 in [S_half^2, 1], 1/(1+T(t)) in [1/(1+T_in), 1]
# Let f(t) = Delta_2(t)/A_0 = r(t) * t^2 where r(t) in [r_lo, 1], r_lo = S_half^2 / (1 + T_in)
S_half = iv.sin(iv.mpf('0.5') * G0['thc_hi']) / (iv.mpf('0.5') * G0['thc_hi'])
r_lo = (S_half ** 2) / (1 + G0['T_in'])
print(f"  r_lo = {float(r_lo.a):.6f} (ratio of Delta_2(t)/t^2 on [0, 0.5] varies by <= {float((1-r_lo).b)*100:.3f}%)")
# Since int_0^{1/2} f(t) dt = (1/2) Del0, f(t) <= (1/r_lo) * 12 * Del0 * t^2 and >= r_lo * 12 * Del0 * t^2
# So for u in [0, 0.5]: h(u) = int_0^u (Del0 - f(t)) dt <= Del0 * (u - 4 * r_lo * u^3), max at u = 1/sqrt(12*r_lo)
u_pk = 1 / iv.sqrt(12 * r_lo)
h_max_over_Del0 = u_pk - 4 * r_lo * (u_pk ** 3)
# And for u in [0.5, 0.515]: -h(0.515) = int_{0.5}^{0.515} (f(t) - Del0) dt <= Del0 * (4/r_lo * (0.515^3 - 0.5^3) - 0.015)
h_neg_over_Del0 = (4 / r_lo) * (iv.mpf('0.515')**3 - iv.mpf('0.5')**3) - iv.mpf('0.015')
print(f"  PROVED max h(u)/Del0 on [0,0.5] <= {float(h_max_over_Del0.b):.6f}, max (-h(u))/Del0 on [0.5,0.515] <= {float(h_neg_over_Del0.b):.6f}, osc <= {float((h_max_over_Del0 + h_neg_over_Del0).b):.6f} << 0.39!")

print("\n=== CHECK E: Straddle bin skipped parabola integral -int_{u0_lo}^{u0_hi} 2 Re(NE_even * W2_par) du ===")
# Test at tau = 2.0, u0 in [0.20, 0.225] (max_du = 0.025 in D2!) and tau = 3.0, u0 in [0.32, 0.335] (D3!)
for tau_f, ua_f, ub_f, th_f in [(2.0, 0.20, 0.225, 2.5), (3.0, 0.32, 0.335, 3.0), (10.0, 0.495, 0.510, 3.1)]:
    u0_mid = 0.5 * (ua_f + ub_f)
    def integrand(u):
        x = tau_f * abs(u0_mid**2 - u**2)
        if u < u0_mid:
            ph0 = math.pi - 2 * math.pi * u
        else:
            ph0 = math.pi + 2 * math.pi * u
        r = math.exp(-x)
        # EVEN part of Li_3(-r e^{i ph})
        z0 = -r * cmath.exp(1j * ph0)
        zp = -r * cmath.exp(1j * (ph0 + th_f))
        zm = -r * cmath.exp(1j * (ph0 - th_f))
        ev = 0.5 * (complex(mp.polylog(3, zp)) + complex(mp.polylog(3, zm))) - complex(mp.polylog(3, z0))
        W2 = (0.75 * tau_f**2) / ((tau_f * u - 1j * math.pi)**4)
        return (-2.0 * ev * W2).real
    xs = np.linspace(ua_f, ub_f, 40)
    val = np.trapezoid([integrand(x) for x in xs], xs)
    # Compare with b_st in code:
    W1_abs = 1.0 / (2 * tau_f * ua_f**2 + 2 * math.pi**2 / tau_f)
    b_st = 4 * (math.pi**2 / 6) * W1_abs * (ub_f - ua_f)
    print(f"  tau={tau_f:4.1f} u=[{ua_f},{ub_f}] th={th_f}: SKIPPED par_ne_re = {val:+.6f} (code b_st={b_st:.6f}, note W1_par on straddle ~ {b_st/2:.6f})")

print("\n=== CHECK F: Exact d_lam_th_ub (ABC L898, D L1418/1490) vs true finite-n Re(d/dloglam Psi_fin(th)) ===")
for c, j, r0, kap, u0 in [(60, 61, 1, 0.015, 0.42), (40, 41, 1, 0.03, 0.25), (25, 26, 1, 0.08, 0.30)]:
    pt = make_point(c, j, r0, kap)
    tau = pt['tau']
    # Compute exact Re(-1 + int_{-L_n}^{L_n} z/(1+z) du) at th = pi (valley!) for u0
    def dlam_exact(th):
        re0 = zeta_fn(pt, u0).real
        L_u = min(pt['n'] / (2 * c), 5.0)
        from scipy.integrate import quad
        def f(u):
            z = cmath.exp(zeta_fn(pt, u) - re0 + 1j * th)
            return (z / (1.0 + z)).real
        val, _ = quad(f, -L_u, L_u, points=[-u0, u0], limit=200, epsabs=1e-9)
        return -1.0 + val
    # Compare with code's d_lam_th_ub:
    G, Del0, _, rel_B1, mult_B2, _, _, _ = g3d.curv_W3_small_u_iv(iv.mpf(str(tau)), iv.mpf(str(tau)), iv.mpf(str(u0)), iv.mpf(str(pt['ct'])))
    s_max = iv.mpf(str(tau)) * iv.mpf(str(u0)) / PI
    d_lam_code = float((-(1 - s_max**2 - 2*Del0) + (1 + rel_B1)*(iv.mpf(str(tau))/(6*PI)) + mult_B2*((iv.mpf(str(tau))**2)/(8*PI**2))*(2*iv.mpf(str(u0)) + iv.sqrt(PI/(G['c_re5']*iv.mpf(str(tau)))))).b)
    ex_pi = dlam_exact(math.pi)
    print(f"  c={c} j={j} tau={tau:.2f} u0={u0}: exact d_lam(pi)={ex_pi:+.5f} <= code d_lam_th_ub={d_lam_code:+.5f} {'OK' if ex_pi <= d_lam_code else 'VIOLATION!!'}")
