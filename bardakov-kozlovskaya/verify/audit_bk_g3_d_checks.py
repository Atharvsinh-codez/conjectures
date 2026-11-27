"""Adversarial mathematical audit checks for Piece G3 Case D (and shared G3-ABC helpers)."""
import math
import cmath
import numpy as np
import mpmath as mp
from mpmath import iv

mp.mp.dps = 50
iv.dps = 50
PI = iv.pi

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

def exact_funcs(pt):
    c, r0, p0, sj, th_c, n, r1, r3, D, tau = (
        pt['c'], pt['r0'], pt['p0'], pt['sj'], pt['th_c'],
        pt['n'], pt['r1'], pt['r3'], pt['D'], pt['tau']
    )
    B = 4.0 * math.pi * c * r0 * th_c
    def zeta(u):
        w = sj * cmath.exp(2j * u * th_c)
        return r0 * n * cmath.log((r3 / r1) / (1.0 - w)) + 1j * math.pi * r0 * pt['j'] - 2j * math.pi * c * p0 * (u + 0.5)
    def A_fn(u):
        w = sj * cmath.exp(2j * u * th_c)
        return 2j * math.pi * c * ((r0 + p0) * w - p0)
    def P_fn(u):
        return A_fn(u) / (1.0 - sj)
    def P_par(u):
        return -2.0 * tau * u + 2j * math.pi * D
    def dzeta(u):
        w = sj * cmath.exp(2j * u * th_c)
        return A_fn(u) / (1.0 - w)
    def W1_fin(u):
        w = sj * cmath.exp(2j * u * th_c)
        return B * w / (A_fn(u) ** 2)
    def W2_fin(u):
        w = sj * cmath.exp(2j * u * th_c)
        A = A_fn(u)
        return (3.0 * B**2 * w**2 + A * B * (2j * th_c) * w * (1.0 + w)) / (A ** 4)
    def W2_par(u):
        return 12.0 * tau**2 / (P_par(u) ** 4)
    def W3_fin(u):
        w = sj * cmath.exp(2j * u * th_c)
        A = A_fn(u)
        return (15.0 * B**3 * w**3 + 10.0 * A * B**2 * (2j * th_c) * w**2 * (1.0 + w)
                - A**2 * B * (4.0 * th_c**2) * w * (1.0 + 4.0 * w + w**2)) / (A ** 6)
    def T_fn(u):
        return 4.0 * sj / ((1.0 - sj) ** 2) * (math.sin(u * th_c) ** 2)
    return zeta, A_fn, P_fn, P_par, dzeta, W1_fin, W2_fin, W2_par, W3_fin, T_fn

import sys
sys.path.insert(0, r"d:\solving-conjectures\work\bardakov-kozlovskaya\analyst")
import certify_G3_D_full_iv as g3d

print("=== CHECK 1: rel_W2, CW3_free_D1, cert_B1_B2_pert_D1 (missing 1/(1-s_j)!) ===")
test_pts = [
    make_point(20, 21, 1, 0.25),   # D=+1, j=c+1, r0=1, tau=5.0, ct=20
    make_point(12, 13, 1, 0.40),   # D=+1, j=c+1, r0=1, tau=4.8, ct=12
    make_point(15, 29, 1, 0.40),   # D=-1, j=2c-1, r0=1, tau=6.0, ct=15
    make_point(40, 41, 1, 0.03),   # D=+1, j=c+1, tau=1.20, ct=40 (Band D1!)
    make_point(25, 26, 1, 0.08),   # D=+1, j=c+1, tau=2.00, ct=25 (Band D2!)
    make_point(60, 61, 1, 0.015),  # D=+1, j=c+1, tau=0.90, ct=60 (Case G3-C!)
    make_point(40, 42, 1, 0.30),   # D=+2, j=c+2, tau=12.0 (Case G3-A!)
    make_point(89, 90, 1, 0.449),  # D=+1, j=c+1, tau=40.0 (Case G3-B!)
]
for pt in test_pts:
    zeta, A_fn, P_fn, P_par, dzeta, W1_fin, W2_fin, W2_par, W3_fin, T_fn = exact_funcs(pt)
    tau = pt['tau']; ct = pt['ct']; kap = pt['kap']; c = pt['c']; r0 = pt['r0']; D = pt['D']
    # Build tightest possible box around this exact point:
    G = g3d.box_geom(iv.mpf(str(kap)))
    r0ck2 = iv.mpf(str(r0 * c * kap**2))
    Del0 = G['C_D0'] * r0ck2
    c_min = iv.mpf(c)
    # Test W2 relative error at u = 0.0 and u = 0.40:
    for u in [0.0, 0.40, 1.0]:
        w2_rel_exact = abs(W2_fin(u) / W2_par(u) - 1.0)
        ep = g3d.eps_P_bound(G, iv.mpf(str(tau)), iv.mpf(str(tau)), iv.mpf(str(u)), iv.mpf(str(u)), c_min, Del0)
        inv4 = 1 / ((1 - ep) ** 4)
        P_over_2tau = iv.sqrt(iv.mpf(str(u))**2 + (PI / iv.mpf(str(tau)))**2) * (1 + ep)
        rel_W2_code = ((G['te_hi']**2 - 1) + 4 * iv.mpf(str(u)) * G['thc_hi']) * inv4 + (inv4 - 1) + inv4 * (iv.mpf(4)/3) * G['te_hi'] * G['thc_hi'] * P_over_2tau
        if w2_rel_exact > float(rel_W2_code.b) and abs(D) == 1:
            print(f"  [VIOLATION rel_W2!] c={c} j={pt['j']} D={D:+d} tau={tau:.2f} kap={kap:.3f} u={u}: EXACT={w2_rel_exact:.6f} > CODE={float(rel_W2_code.b):.6f}")
        # Also test W3 ratio:
        w3_ratio_exact = abs(W3_fin(u)) * (abs(P_fn(u)) ** 6) / (120.0 * tau**3)
        T_loc = iv.mpf(str(max(T_fn(u), float(G['T_in'].b))))
        cw3_code = g3d.CW3_free_D1(G, c_min, Del0, T_loc)
        if w3_ratio_exact > float(cw3_code.b) and abs(D) == 1:
            print(f"  [VIOLATION CW3_free_D1!] c={c} j={pt['j']} D={D:+d} tau={tau:.2f} kap={kap:.3f} u={u}: EXACT={w3_ratio_exact:.6f} > CODE={float(cw3_code.b):.6f}")
    # Also test |W2_fin(u)/W2_par(u)| vs mult_W2 in cert_B1_B2_pert_D1 at u0 = 0.01, 0.25:
    for u0 in [0.01, 0.25]:
        rel_B1, mult_B2 = g3d.cert_B1_B2_pert_D1(G, c_min, iv.mpf(str(u0)), Del0, r0ck2)
        b2_ratio_exact = abs((W2_fin(u0) / dzeta(u0)) / (W2_par(u0) / P_par(u0)))
        if b2_ratio_exact > float(mult_B2.b) and abs(D) == 1:
            print(f"  [VIOLATION mult_B2!] c={c} j={pt['j']} D={D:+d} tau={tau:.2f} kap={kap:.3f} u0={u0}: EXACT={b2_ratio_exact:.6f} > CODE={float(mult_B2.b):.6f}")

print("\n=== CHECK 2: Exact |W2_fin| and |W3_fin| second-term ratio (true vs code) ===")
for pt in test_pts:
    zeta, A_fn, P_fn, P_par, dzeta, W1_fin, W2_fin, W2_par, W3_fin, T_fn = exact_funcs(pt)
    kap = pt['kap']; c = pt['c']; r0 = pt['r0']; th_c = pt['th_c']; sj = pt['sj']
    G = g3d.box_geom(iv.mpf(str(kap)))
    # Term 2 of W2_fin / W2_par at u=0:
    u = 0.0
    t2_w2_true = abs(A_fn(u) * (4*math.pi*c*r0*th_c) * (2j*th_c) * sj * (1 + sj) / (A_fn(u)**4)) / abs(W2_par(u))
    ep = g3d.eps_P_bound(G, iv.mpf(str(pt['tau'])), iv.mpf(str(pt['tau'])), iv.mpf(0), iv.mpf(0), iv.mpf(c), G['C_D0']*iv.mpf(str(r0*c*kap**2)))
    P_over_2tau = (PI / iv.mpf(str(pt['tau']))) * (1 + ep)
    t2_w2_code = float(((1/((1-ep)**4)) * (iv.mpf(4)/3) * G['te_hi'] * G['thc_hi'] * P_over_2tau).b)
    # P_th exact vs code at u=0 and u=1.0:
    Pth_exact_0 = abs(P_fn(0.0)) * th_c / (pt['tau'] * (1.0 - sj))
    a_p_code = float(((1 + G['C_D0']*iv.mpf(str(r0*c*kap**2))) / (2 * G['S']**2 * iv.mpf(c))).b)
    # T-term of P_th:
    u_t = 1.0
    T1 = T_fn(u_t)
    b_p_exact = math.sqrt(abs(P_fn(u_t))**2 - abs(P_fn(0.0))**2) * th_c / (pt['tau'] * (1.0 - sj) * math.sqrt(T1))
    f_p_hi = (1 + 1/iv.mpf(c)) * (1 + 1/(2*iv.mpf(c)))
    b_p_code = float(iv.sqrt(f_p_hi * G['te_hi'] / (2 * G['S']**4)).b)
    print(f"  c={c:2d} j={pt['j']:2d} D={pt['D']:+d} kap={kap:.3f}: T2_W2 true={t2_w2_true:.5f} vs code={t2_w2_code:.5f} | a_p true={Pth_exact_0:.5f} vs code={a_p_code:.5f} | b_p true={b_p_exact:.5f} vs code={b_p_code:.5f}")

print("\n=== CHECK 3: |P(u)|^2 vs 4*((tau*u)^2 + (d_lo*pi)^2) in exact_R3_W3_D1 ===")
for pt in test_pts:
    zeta, A_fn, P_fn, P_par, dzeta, W1_fin, W2_fin, W2_par, W3_fin, T_fn = exact_funcs(pt)
    tau = pt['tau']; kap = pt['kap']; c = pt['c']; r0 = pt['r0']; D = pt['D']
    G = g3d.box_geom(iv.mpf(str(kap)))
    Del0 = float((G['C_D0'] * iv.mpf(str(r0 * c * kap**2))).b)
    d_lo = abs(D) - Del0
    for u in [0.3, 1.0, 2.0, 4.0]:
        P2_exact = abs(P_fn(u)) ** 2
        P2_code_lb = 4.0 * ((tau * u)**2 + (d_lo * math.pi)**2)
        if P2_exact < P2_code_lb:
            print(f"  [VIOLATION |P(u)|^2 LB!] c={c} j={pt['j']} D={D:+d} tau={tau:.2f} kap={kap:.3f} u={u}: EXACT={P2_exact:.4f} < CODE_LB={P2_code_lb:.4f} (ratio={P2_exact/P2_code_lb:.4f}, 1/ratio^3={(P2_code_lb/P2_exact)**3:.4f})")

print("\n=== CHECK 4: y > 0 outer poles (|x_p| > u_0, y_p > 0) vs E_out_0 and c_out ===")
# Solve 1 + exp(zeta(u) - Re zeta(u0) + i*th) = 0 numerically in y > 0 for x > u0!
for pt in [make_point(20, 21, 1, 0.25), make_point(12, 13, 1, 0.40), make_point(30, 31, 1, 0.40), make_point(40, 41, 1, 0.03), make_point(89, 90, 1, 0.449)]:
    zeta, A_fn, P_fn, P_par, dzeta, W1_fin, W2_fin, W2_par, W3_fin, T_fn = exact_funcs(pt)
    tau = pt['tau']; kap = pt['kap']; c = pt['c']; r0 = pt['r0']; ct = pt['ct']
    u0 = 0.45 if tau > 3 else 0.25
    re0 = zeta(u0).real
    G = g3d.box_geom(iv.mpf(str(kap)))
    Del0 = G['C_D0'] * iv.mpf(str(r0 * c * kap**2))
    d_lo = 1 - Del0; d_hi = 1 + Del0
    C_x = iv.sqrt(G['te_str'] / G['c_re50'])
    C_quad = (3 + d_hi / d_lo) * G['te_str'] * C_x
    E_out_0 = float(((4 * PI**2 * d_lo / G['te_str']) * (iv.mpf(ct) / iv.mpf(str(tau)))).a)
    c_out = float((2 * PI * iv.sqrt((2 * PI * iv.mpf(ct) * (iv.mpf(ct) / iv.mpf(str(tau)))) / C_quad)).a)
    # Find where Re zeta(x + iy) = re0 in x > u0, y in (0, 0.50*r0]
    min_E_pos = 1e9; worst_xy = None
    for x in np.linspace(u0 + 0.02, min(pt['n']/(2*c) - 0.1, 15.0), 400):
        # Check if Re zeta(x + iy) - re0 ever crosses 0 for y in (0, 0.50*r0]
        ys = np.linspace(0.0, 0.50 * r0, 100)
        vals = [zeta(x + 1j * y).real - re0 for y in ys]
        for idx in range(len(ys) - 1):
            if vals[idx] * vals[idx + 1] <= 0:
                # bisect y
                ya, yb = ys[idx], ys[idx + 1]
                for _ in range(30):
                    ym = 0.5 * (ya + yb)
                    if (zeta(x + 1j * ym).real - re0) * vals[idx] <= 0:
                        yb = ym
                    else:
                        ya = ym
                ym = 0.5 * (ya + yb)
                E_val = 2 * math.pi * ct * ym
                if E_val < min_E_pos:
                    min_E_pos = E_val; worst_xy = (x, ym, T_fn(x))
                break
    print(f"  c={c:2d} j={pt['j']:2d} tau={tau:5.2f} kap={kap:.3f}: E_out_0(code)={E_out_0:.2f}, c_out={c_out:.2f} | min exact E(x>u0, y>0)={min_E_pos:.2f} at {worst_xy}")

print("\n=== CHECK 5: Re P(u) vs dP_re_rel (L489) for D=-1 ===")
for c, j, r0, kap in [(15, 29, 1, 0.40), (10, 19, 1, 0.44), (20, 39, 1, 0.30), (100, 199, 1, 0.20)]:
    pt = make_point(c, j, r0, kap)
    zeta, A_fn, P_fn, P_par, dzeta, *_ = exact_funcs(pt)
    G = g3d.box_geom(iv.mpf(str(pt['kap'])))
    Del0 = G['C_D0'] * iv.mpf(str(r0 * c * pt['kap']**2))
    for u in [0.05, 0.3, 0.5, 1.0]:
        re_err_exact = abs(P_fn(u).real - P_par(u).real) / abs(P_par(u))
        rat_re = (2 * iv.mpf(str(pt['tau'])) * iv.mpf(str(u))) / iv.sqrt((2 * iv.mpf(str(pt['tau'])) * iv.mpf(str(u)))**2 + (2*PI)**2)
        dP_re_code = float((((G['te_hi'] - 1) + (1 + Del0) / (2 * iv.mpf(c)) + ((2 * iv.mpf(str(u)) * G['thc_hi'])**2)/6) * rat_re).b)
        if re_err_exact > dP_re_code:
            print(f"  [VIOLATION dP_re_rel!] c={c} j={j} D={pt['D']} u={u}: EXACT={re_err_exact:.7f} > CODE={dP_re_code:.7f}")

print("\n=== CHECK 6: Does gam * L_iv^2 <= 0.05 and L_iv * thc_hi <= 0.080 hold across D1/D2/D3? ===")
max_gamL2 = 0.0; max_Lthc = 0.0; worst_box = None
for ta_f, tb_f, ct_lo_f in [(1.05, 1.25, 25.6), (1.50, 1.65, 14.7), (2.60, 2.80, 10.0), (2.80, 2.98, 10.0), (4.38, 4.60, 10.0), (4.60, 4.82, 10.0), (9.5, 11.2, 10.0), (32.0, 40.0, 17.07)]:
    G, c_min, r0ck2, Del0 = g3d.make_params(iv.mpf(str(ta_f)), iv.mpf(str(tb_f)), iv.mpf(str(ct_lo_f)))
    L_val = float(iv.sqrt(iv.mpf('0.5')**2 + iv.mpf('11.0') / (G['c_re5'] * iv.mpf(str(ta_f)))).b)
    gL2 = float(G['gam'].b) * L_val**2
    Lthc = L_val * float(G['thc_hi'].b)
    if gL2 > max_gamL2: max_gamL2 = gL2; worst_box = (ta_f, tb_f, ct_lo_f, L_val, gL2, Lthc)
    if Lthc > max_Lthc: max_Lthc = Lthc
    print(f"  tau=[{ta_f:4.2f},{tb_f:4.2f}] ct>={ct_lo_f:5.2f}: L={L_val:.3f}, gam*L^2={gL2:.5f} (<=0.05?), L*thc={Lthc:.5f} (<=0.080?)")

print("\n=== CHECK 7: Exact err_HB1 in curv_W3_small_u_iv vs true finite-n H_fin0 - Re B1_fin ===")
for c, j, r0, kap in [(40, 41, 1, 0.03), (25, 26, 1, 0.07), (15, 16, 1, 0.10), (15, 29, 1, 0.10), (40, 79, 1, 0.03)]:
    pt = make_point(c, j, r0, kap)
    zeta, A_fn, P_fn, P_par, dzeta, W1_fin, W2_fin, W2_par, W3_fin, T_fn = exact_funcs(pt)
    tau = pt['tau']; D = pt['D']
    G, c_min, r0ck2, Del0 = g3d.make_params(iv.mpf(str(tau)), iv.mpf(str(tau)), iv.mpf(str(pt['ct'])))
    u0 = 0.25
    rel_B1, mult_B2 = g3d.cert_B1_B2_pert_D1(G, c_min, iv.mpf(str(u0)), Del0, r0ck2)
    err_HB1_code = float(((1 - G['r_re0']) + (G['te_hi'] - 1) + rel_B1 + 2 * Del0).b)
    # Exact H_fin0 - Re B1_fin(0):
    def d2z(u):
        w = pt['sj'] * cmath.exp(2j * u * pt['th_c'])
        return -(4.0 * math.pi * c * r0 * pt['th_c'] * w) / ((1.0 - w)**2)
    def wrap(x): return (x + math.pi) % (2 * math.pi) - math.pi
    ap0 = wrap(zeta(u0).imag); am0 = wrap(zeta(-u0).imag)
    H_fin0_ex = -2.0 * (1.0 / dzeta(u0)).real
    B1_fin_ex = (1j * (ap0 * d2z(u0) / (dzeta(u0)**3) - am0 * d2z(-u0) / (dzeta(-u0)**3))).real
    # Parabola value:
    s = tau * u0 / math.pi
    par_val = (2*u0/(1+s**2) + (1-2*u0)*(1-3*s**2)/((1+s**2)**3))
    ex_val = (H_fin0_ex - B1_fin_ex) / (tau / (2 * math.pi**2))
    diff = abs(ex_val - par_val)
    print(f"  c={c} j={j} D={D:+d} tau={tau:.2f}: |ex - par|={diff:.6f} vs err_HB1_code={err_HB1_code:.6f} {'[VIOLATION!!]' if diff > err_HB1_code else 'OK'}")
