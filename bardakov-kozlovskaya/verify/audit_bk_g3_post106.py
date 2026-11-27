"""Post-#106 adversarial verification of all 5 fixed formulas in certify_G3_ABC_full_iv.py and certify_G3_D_full_iv.py."""
import math
import cmath
import sys
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

def exact_funcs(pt):
    c, r0, p0, sj, th_c, n, r1, r3, D, tau = (
        pt['c'], pt['r0'], pt['p0'], pt['sj'], pt['th_c'],
        pt['n'], pt['r1'], pt['r3'], pt['D'], pt['tau']
    )
    B = 4.0 * math.pi * c * r0 * th_c
    def A_fn(u):
        w = sj * cmath.exp(2j * u * th_c)
        return 2j * math.pi * c * ((r0 + p0) * w - p0)
    def P_fn(u): return A_fn(u) / (1.0 - sj)
    def P_par(u): return -2.0 * tau * u + 2j * math.pi * D
    def dzeta(u):
        w = sj * cmath.exp(2j * u * th_c)
        return A_fn(u) / (1.0 - w)
    def W2_fin(u):
        w = sj * cmath.exp(2j * u * th_c)
        A = A_fn(u)
        return (3.0 * B**2 * w**2 + A * B * (2j * th_c) * w * (1.0 + w)) / (A ** 4)
    def W2_par(u): return 12.0 * tau**2 / (P_par(u) ** 4)
    def W3_fin(u):
        w = sj * cmath.exp(2j * u * th_c)
        A = A_fn(u)
        return (15.0 * B**3 * w**3 + 10.0 * A * B**2 * (2j * th_c) * w**2 * (1.0 + w)
                - A**2 * B * (4.0 * th_c**2) * w * (1.0 + 4.0 * w + w**2)) / (A ** 6)
    def T_fn(u): return 4.0 * sj / ((1.0 - sj) ** 2) * (math.sin(u * th_c) ** 2)
    return A_fn, P_fn, P_par, dzeta, W2_fin, W2_par, W3_fin, T_fn

test_pts = [
    make_point(20, 21, 1, 0.25),
    make_point(12, 13, 1, 0.40),
    make_point(15, 29, 1, 0.40),
    make_point(40, 41, 1, 0.03),
    make_point(25, 26, 1, 0.08),
    make_point(60, 61, 1, 0.015),
    make_point(89, 90, 1, 0.449),
    make_point(10, 11, 1, 0.44),
    make_point(10, 19, 1, 0.44),
    make_point(21, 32, 2, 0.22),
]

violations = 0
for pt in test_pts:
    A_fn, P_fn, P_par, dzeta, W2_fin, W2_par, W3_fin, T_fn = exact_funcs(pt)
    tau, ct, kap, c, r0, D, th_c, sj = pt['tau'], pt['ct'], pt['kap'], pt['c'], pt['r0'], pt['D'], pt['th_c'], pt['sj']
    G = g3d.box_geom(iv.mpf(str(kap)))
    r0ck2 = iv.mpf(str(r0 * c * kap**2))
    Del0 = G['C_D0'] * r0ck2
    c_min = iv.mpf(c)
    f_p_lo = (1 - 1 / c_min) * (1 - 1 / (2 * c_min))
    mob_in  = f_p_lo * (G['S'] ** 4) * G['te_lo'] * G['sinc_in_sq']
    mob_out = f_p_lo * (G['S'] ** 4) * G['te_lo'] * G['sinc_c5_sq']
    d_lo = 1 - Del0

    # 1. thc_1ms >= th_c / (1 - sj)
    if th_c / (1.0 - sj) > float(G['thc_1ms'].b):
        print(f"  [FAIL thc_1ms] c={c} j={pt['j']}: {th_c/(1-sj)} > {float(G['thc_1ms'].b)}")
        violations += 1

    # 2. rel_W2, CW3_free_D1, |P(u)|^2 with mob_out, dP_re_rel
    for u in [0.0, 0.3, 0.45, 1.0, 2.0]:
        if T_fn(u) > 0.05: continue
        ep = g3d.eps_P_bound(G, iv.mpf(str(tau)), iv.mpf(str(tau)), iv.mpf(str(u)), iv.mpf(str(u)), c_min, Del0)
        inv4 = 1 / ((1 - ep) ** 4)
        P_over_2tau = iv.sqrt(iv.mpf(str(u))**2 + (PI / iv.mpf(str(tau)))**2) * (1 + ep)
        rel_W2 = ((G['te_hi']**2 - 1) + 4 * iv.mpf(str(u)) * G['thc_hi']) * inv4 + (inv4 - 1) + inv4 * (iv.mpf(4)/3) * G['te_hi'] * G['thc_1ms'] * P_over_2tau
        w2_rel_ex = abs(W2_fin(u) / W2_par(u) - 1.0)
        if w2_rel_ex > float(rel_W2.b):
            print(f"  [FAIL rel_W2] c={c} j={pt['j']} u={u}: {w2_rel_ex} > {float(rel_W2.b)}")
            violations += 1
        # CW3_free_D1
        T_loc = iv.mpf(str(max(T_fn(u), float(G['T_in'].b))))
        cw3 = g3d.CW3_free_D1(G, c_min, Del0, T_loc)
        w3_ex = abs(W3_fin(u)) * (abs(P_fn(u))**6) / (120.0 * tau**3)
        if w3_ex > float(cw3.b):
            print(f"  [FAIL CW3_free_D1] c={c} j={pt['j']} u={u}: {w3_ex} > {float(cw3.b)}")
            violations += 1
        # |P(u)|^2 >= 4*(mob * tau^2 * u^2 + (d_lo*pi)^2)
        mob = mob_in if u <= 0.52 else mob_out
        P2_lb = float((4 * (mob * (iv.mpf(str(tau)) * iv.mpf(str(u)))**2 + (d_lo * PI)**2)).a)
        if abs(P_fn(u))**2 < P2_lb:
            print(f"  [FAIL |P|^2] c={c} j={pt['j']} u={u}: {abs(P_fn(u))**2} < {P2_lb}")
            violations += 1

    # 3. a_p and b_p individually
    a_p_ex = abs(P_fn(0.0)) * th_c / (tau * (1.0 - sj))
    a_p_bd = float(((1 + Del0) / (G['S'] * iv.cos(G['al']) * c_min)).b)
    u_t = 0.5
    b_p_ex = math.sqrt(abs(P_fn(u_t))**2 - abs(P_fn(0.0))**2) * th_c / (tau * (1.0 - sj) * math.sqrt(T_fn(u_t)))
    b_p_bd = float((iv.sqrt((1 + 1/c_min) * (2 + 1/c_min)) / (G['S'] * iv.cos(G['al']))).b)
    if a_p_ex > a_p_bd or b_p_ex > b_p_bd:
        print(f"  [FAIL a_p/b_p] c={c} j={pt['j']}: a_p {a_p_ex} vs {a_p_bd}, b_p {b_p_ex} vs {b_p_bd}")
        violations += 1

    # 4. far_T5 C_W3_T5 check at T(u) >= 0.05 (sample u up to L_n)
    far_r5, C_W3_T5 = g3d.far_T5_helpers_D1(G, c_min, Del0)
    L_n = pt['n'] / (2.0 * c)
    for frac in [0.1, 0.25, 0.5, 0.8, 0.99]:
        u = frac * L_n
        if T_fn(u) <= 0.05: continue
        w3_norm = abs(W3_fin(u)) * (abs(P_fn(u))**4) / tau
        if w3_norm > float(C_W3_T5.b):
            print(f"  [FAIL C_W3_T5] c={c} j={pt['j']} u={u:.2f} T={T_fn(u):.3f}: {w3_norm} > {float(C_W3_T5.b)}")
            violations += 1

    # 5. int_mlog_ub check vs numerical quad near u0
    u0 = 0.35
    ua = 0.348; ub = 0.35
    alpha = G['c_rein'] * iv.mpf(str(tau)) * iv.mpf(str(u0 + ua))
    bd = float(g3d.int_mlog_ub(alpha, iv.mpf(0), iv.mpf('0.002'), iv.mpf(0), iv.mpf(1), iv.mpf('0.002')).b)
    from scipy.integrate import quad
    ex_int, _ = quad(lambda u: -math.log(1.0 - math.exp(-float(G['c_rein'].a) * tau * (u0**2 - u**2))), ua, ub, epsabs=1e-12)
    if ex_int > bd:
        print(f"  [FAIL int_mlog_ub] c={c}: exact={ex_int} > bd={bd}")
        violations += 1

print(f"\nPOST-#106 ADVERSARIAL VERIFICATION COMPLETE: {violations} violations across all {len(test_pts)} extreme Piece G3 points!")
assert violations == 0
