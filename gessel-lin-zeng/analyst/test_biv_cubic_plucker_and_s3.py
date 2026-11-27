"""DECISIVE TEST OF:
1. BIVARIATE CUBIC PLUCKER CONJECTURE:
   Given TWO overlapping real-stable cubics:
     P1(t, y) = A + 3 B y + 3 C y^2 + E y^3
     P2(t, y) = B + 3 C y + 3 E y^2 + F y^3
   Is P_new(t, y) = (A+B) + 3(B+C)y + 3(C+E)y^2 + (E+F)y^3 ALWAYS REAL-STABLE?!
   (Test at fixed t_0 in R: P(t_0, y) and P_t(t_0, y) interlacing!
    Wait: by Hermite-Biehler, P(t, y) is real-stable iff at EVERY t_0 in R,
    P(t_0, y) has all real roots in y AND P_t(t_0, y) interlaces P(t_0, y) in the SAME orientation!!
    THAT'S A PURELY POINTWISE LINEAR-IN-t CONDITION ON 10 REAL NUMBERS
    (a,b,c,e,f) AND (a',b',c',e',f') AT A SINGLE POINT t_0!!
    LET'S TEST 10,000 RANDOM + BOUNDARY REAL-STABLE OVERLAPPING CUBIC JETS!!)
2. FIXED Stage 3 pole-gap test at d = 2 (fix W sign!) for j = 20, 40, 80!
"""
import numpy as np

print("=== 1. POINTWISE RAYLEIGH / INTERLACING TEST OF BIVARIATE CUBIC PLUCKER (20,000 TRIALS!) ===")
# At a fixed t_0 in R, a cubic jet (P(y), Q(y)) = (p0 + 3 p1 y + 3 p2 y^2 + p3 y^3, q0 + 3 q1 y + 3 q2 y^2 + q3 y^3)
# satisfies the real-stability jet condition (Q << P in same orientation, i.e. Q(y_j)*P'(y_j) >= 0 at all 3 roots of P)
# iff P has 3 real roots r1 <= r2 <= r3 AND Q(y) = c0*P(y) + sum_{m=1}^3 c_m * P(y)/(y - r_m) with c_1, c_2, c_3 >= 0!
#Wait: for TWO OVERLAPPING cubics (a, b, c, e) and (b, c, e, f) with derivatives (a', b', c', e') and (b', c', e', f'):
# They share (b, c, e) and (b', c', e')!
# How to generate ALL overlapping real-stable cubic jets?
# Start from a real-stable quartic? No, MORE GENERAL:
# 1. Pick P1 = a + 3 b y + 3 c y^2 + e y^3 with 3 negative roots -r1, -r2, -r3 (r_i > 0), and Q1 = a' + 3b'y + 3c'y^2 + e'y^3
#    with c0 in R and c1, c2, c3 >= 0!
# 2. That fixes (b, c, e) AND (b', c', e')!
# 3. Now pick ANY f > 0 such that P2 = b + 3 c y + 3 e y^2 + f y^3 has 3 negative roots s1 < s2 < s3 < 0!
# 4. Now check the interval of f' in R such that Q2 = b' + 3 c' y + 3 e' y^2 + f' y^3 satisfies Q2(s_m) * P2'(s_m) >= 0
#    at all 3 roots s1 < s2 < s3 of P2!!
#    Notice Q2(s_m) = (b' + 3 c' s_m + 3 e' s_m^2) + f' * s_m^3 is LINEAR in f'!!
#    So the valid f' range is an EXACT interval [fp_lo, fp_hi]!!
# 5. Test BOTH endpoints fp_lo, fp_hi AND midpoint of f' for P_new = (a+b) + 3(b+c)y + 3(c+e)y^2 + (e+f)y^3
#    and Q_new = (a'+b') + 3(b'+c')y + 3(c'+e')y^2 + (e'+f')y^3!!

def cubic_roots_and_rayleigh(p0, p1, p2, p3, q0, q1, q2, q3):
    """Return (all_real, min_QPprime_sign_check) for P = p0 + 3 p1 y + 3 p2 y^2 + p3 y^3, Q = q0 + 3 q1 y + 3 q2 y^2 + q3 y^3.
    Orientation check: Q(y_j) * P'(y_j) >= 0 at all 3 roots y_j of P!"""
    rts = np.roots([p3, 3*p2, 3*p1, p0])
    if np.max(np.abs(rts.imag)) > 1e-7:
        return False, -1e9
    rts = np.sort(rts.real)
    vals = []
    for yj in rts:
        dp = 3*p3*yj*yj + 6*p2*yj + 3*p1
        qv = q3*yj**3 + 3*q2*yj*yj + 3*q1*yj + q0
        vals.append(dp * qv)
    scale = max(abs(v) for v in vals) + 1e-300
    return True, min(vals) / scale

rng = np.random.default_rng(2026)
valid_jets = 0; fails_jet = 0; min_rel_ray = 1e9; first_fail = None
for trial in range(25000):
    r1, r2, r3 = 10 ** rng.uniform(-1.5, 1.5, 3)
    # P1(y) = (y + r1)(y + r2)(y + r3) = r1 r2 r3 + (r1r2+r1r3+r2r3)y + (r1+r2+r3)y^2 + y^3
    a = r1 * r2 * r3
    b = (r1*r2 + r1*r3 + r2*r3) / 3.0
    c = (r1 + r2 + r3) / 3.0
    e = 1.0
    # Q1(y) = c0*P1(y) + c1*(y+r2)(y+r3) + c2*(y+r1)(y+r3) + c3*(y+r1)(y+r2) with c1,c2,c3 >= 0!
    c0 = rng.normal(0, 2.0)
    c1, c2, c3 = 10 ** rng.uniform(-2, 2, 3)
    if rng.random() < 0.3: c1 = 0.0  # boundary!
    if rng.random() < 0.3: c2 = 0.0
    q_c0 = c0*a + c1*r2*r3 + c2*r1*r3 + c3*r1*r2
    q_c1 = c0*3*b + c1*(r2+r3) + c2*(r1+r3) + c3*(r1+r2)
    q_c2 = c0*3*c + c1 + c2 + c3
    q_c3 = c0*e
    ap, bp, cp, ep = q_c0, q_c1/3.0, q_c2/3.0, q_c3

    # Valid range for f > 0 such that P2 = b + 3 c y + 3 e y^2 + f y^3 is real-rooted:
    A_q = -b*b; B_q = 6*b*c*e - 4*c**3; C_q = 3*c*c*e*e - 4*b*e**3
    dd = B_q*B_q - 4*A_q*C_q
    if dd <= 0: continue
    sd = np.sqrt(dd)
    f_lo = max(1e-9, (-B_q + sd) / (2*A_q))
    f_hi = (-B_q - sd) / (2*A_q)
    if f_hi <= f_lo: continue
    f = rng.uniform(f_lo, f_hi)
    rts2 = np.roots([f, 3*e, 3*c, b])
    if np.max(np.abs(rts2.imag)) > 1e-7: continue
    rts2 = np.sort(rts2.real)
    # Find valid interval [fp_lo, fp_hi] for f' so that Q2(s_m) * P2'(s_m) >= 0 for all m = 0,1,2!
    # Q2(s) = (bp + 3 cp s + 3 ep s^2) + fp * s^3.
    # Require dp * (const + fp * s^3) >= 0!
    fp_lo, fp_hi = -1e100, 1e100
    ok_fp = True
    for s in rts2:
        dp = 3*f*s*s + 6*e*s + 3*c
        cst = bp + 3*cp*s + 3*ep*s*s
        coef = dp * (s**3)
        rhs = -dp * cst
        if abs(coef) < 1e-14:
            if rhs > 1e-12: ok_fp = False; break
        elif coef > 0:
            fp_lo = max(fp_lo, rhs / coef)
        else:
            fp_hi = min(fp_hi, rhs / coef)
    if not ok_fp or fp_lo > fp_hi + 1e-12: continue
    if fp_lo < -1e50 or fp_hi > 1e50: continue
    valid_jets += 1
    for fp in [fp_lo, fp_hi, 0.5*(fp_lo + fp_hi)]:
        ok_r, mval = cubic_roots_and_rayleigh(a+b, b+c, c+e, e+f, ap+bp, bp+cp, cp+ep, ep+fp)
        if not ok_r or mval < -1e-6:
            fails_jet += 1
            if first_fail is None:
                # Also compute Plucker R of (a,b,c,e) and (b,c,e,f)!
                w_ab = ap*b - a*bp; w_bc = bp*c - b*cp; w_ce = cp*e - c*ep; w_ef = ep*f - e*fp
                first_fail = (mval, w_ab*w_ce/(w_bc**2 + 1e-30), w_bc*w_ef/(w_ce**2 + 1e-30))
            break
        if mval < min_rel_ray: min_rel_ray = mval

print(f"  Valid overlapping real-stable cubic jets = {valid_jets}: FAILS = {fails_jet}, min rel ray = {min_rel_ray:+.3e}, first_fail = {first_fail}")

print("\n=== 2. FIXED STAGE 3 POLE-GAP TEST AT d = 2 (CHECK W SIGN FIRST!) ===")
def apply_D_flt(m, p):
    res = np.zeros(len(p) + 1)
    for r, c in enumerate(p): res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def eval_W_flt(f, g, ts):
    fp = np.arange(1, len(f))*f[1:]; gp = np.arange(1, len(g))*g[1:]
    return np.polyval(fp[::-1], ts)*np.polyval(g[::-1], ts) - np.polyval(f[::-1], ts)*np.polyval(gp[::-1], ts)

ts_r = np.concatenate([-10**np.linspace(2, -4, 30), np.linspace(-0.99, -0.01, 60), [0.0], 10**np.linspace(-4, 2, 30)])

# Check sign of W(f0, D_r f1) on actual GLZ k=2: f2=[1], f1=[1, 2], f0=D_3(1+2t)=[1, 7, 8]?
# Wait: for f0 = (1+2t)(1+4t), f1 = (1+3t), f2 = 1: W(f0, D_r f1)(0) = f0'(0)*1 - 1*(D_r f1)'(0) = 6 - (6 + r) = -r < 0!!
# AHA!! W(f0, D_r f1) IS NEGATIVE (f' g - f g' < 0 when f < g in (-1, 0)!!)!
# So pole-gap condition is: max_t W(f0, D_r f1) < 0 AND max_t W(D_r f1, D_r^2 f2) < 0!
for j in [20, 40, 80]:
    max_drop = -1e9; valid_s3 = 0; best_cfg = None
    for _ in range(5000):
        u = np.sort(1.0 + 10 ** rng.uniform(-1, 3, 3))
        b2, a, b1 = u[0], u[1], u[2]
        f2 = np.array([1.0]); f1 = np.array([1.0, a]); f0 = np.array([1.0, b1+b2, b1*b2])
        r_cut = 0.75 * j
        Mr = apply_D_flt(r_cut, f1); Yr = apply_D_flt(r_cut, apply_D_flt(r_cut, f2))
        if np.max(eval_W_flt(f0, Mr, ts_r)) >= 0 or np.max(eval_W_flt(Mr, Yr, ts_r)) >= 0:
            continue
        Mj = apply_D_flt(j, f1); Yj = apply_D_flt(j, apply_D_flt(j, f2))
        w1 = eval_W_flt(f0, Mj, ts_r); w2 = eval_W_flt(Mj, Yj, ts_r); w12 = eval_W_flt(f0, Yj, ts_r)
        if np.max(w1) >= 0 or np.max(w2) >= 0: continue
        Jj = float(np.min(4*w1*w2 / (w12*w12 + 1e-300)))
        if Jj < 1.0: continue
        valid_s3 += 1
        Mj1 = apply_D_flt(j+1, f1); Yj1 = apply_D_flt(j+1, apply_D_flt(j+1, f2))
        w1p = eval_W_flt(f0, Mj1, ts_r); w2p = eval_W_flt(Mj1, Yj1, ts_r); w12p = eval_W_flt(f0, Yj1, ts_r)
        Jj1 = float(np.min(4*w1p*w2p / (w12p*w12p + 1e-300)))
        drop = Jj - Jj1
        if drop > max_drop:
            max_drop = drop; best_cfg = (round(b2,2), round(a,2), round(b1,2), round(Jj,4), round(Jj1,4))
    print(f"  j={j:2d} (r*<=0.75j, valid={valid_s3}): max(J(j)-J(j+1))={max_drop:+.6f} (j*drop={j*max_drop:+.4f}, j^2*drop={j*j*max_drop:+.3f}), cfg={best_cfg}")
