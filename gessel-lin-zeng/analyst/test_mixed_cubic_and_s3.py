"""1. Test Bivariate Cubic Plucker for MIXED-SIGN y-roots (t_0 INSIDE (-1, 0) where A,B,C,E,F can have ANY real signs!!).
2. Fix Stage 3 d=2 test (check why valid=0 — probably sign of W(f0, D_r f1) or sign of b1, b2!)."""
import numpy as np

print("=== 1. BIVARIATE CUBIC PLUCKER WITH ARBITRARY REAL (MIXED-SIGN!) y-ROOTS ===")
def cubic_roots_and_rayleigh(p0, p1, p2, p3, q0, q1, q2, q3):
    if abs(p3) < 1e-12: return None, 0.0
    rts = np.roots([p3, 3*p2, 3*p1, p0])
    if np.max(np.abs(rts.imag)) > 1e-6:
        return False, -1e9
    rts = np.sort(rts.real)
    vals = []
    for yj in rts:
        dp = 3*p3*yj*yj + 6*p2*yj + 3*p1
        qv = q3*yj**3 + 3*q2*yj*yj + 3*q1*yj + q0
        vals.append(dp * qv)
    scale = max(abs(v) for v in vals) + 1e-300
    # Check ALL same sign (here positive since P1, P2 had positive orientation!)
    return True, min(vals) / scale

rng = np.random.default_rng(2026)
valid_mixed = 0; fails_mixed = 0; first_fail = None
for _ in range(20000):
    # Pick 3 ARBITRARY REAL roots r1, r2, r3 (positive AND negative!) for P1(y) = (y - r1)(y - r2)(y - r3)
    r1, r2, r3 = rng.normal(0, 2.0, 3)
    # P1(y) = (-r1 r2 r3) + (r1 r2 + r1 r3 + r2 r3) y - (r1 + r2 + r3) y^2 + y^3
    a = -r1 * r2 * r3
    b = (r1*r2 + r1*r3 + r2*r3) / 3.0
    c = -(r1 + r2 + r3) / 3.0
    e = 1.0
    c0 = rng.normal(0, 2.0)
    c1, c2, c3 = 10 ** rng.uniform(-2, 2, 3)
    # Q1(y) = c0*P1(y) + sum c_m P1(y)/(y - r_m)
    # P1(y)/(y - r1) = r2 r3 - (r2 + r3) y + y^2
    q_c0 = c0*a + c1*r2*r3 + c2*r1*r3 + c3*r1*r2
    q_c1 = c0*3*b - c1*(r2+r3) - c2*(r1+r3) - c3*(r1+r2)
    q_c2 = c0*3*c + c1 + c2 + c3
    q_c3 = c0*e
    ap, bp, cp, ep = q_c0, q_c1/3.0, q_c2/3.0, q_c3

    # Valid range for f in R such that P2 = b + 3 c y + 3 e y^2 + f y^3 has 3 real roots:
    # Disc/27 = -b^2 f^2 + (6 b c e - 4 c^3) f + (3 c^2 e^2 - 4 b e^3) >= 0
    A_q = -b*b; B_q = 6*b*c*e - 4*c**3; C_q = 3*c*c*e*e - 4*b*e**3
    dd = B_q*B_q - 4*A_q*C_q
    if dd <= 0 or abs(b) < 1e-6: continue
    sd = np.sqrt(dd)
    f_lo = (-B_q + sd) / (2*A_q); f_hi = (-B_q - sd) / (2*A_q)
    if f_hi <= f_lo: continue
    f = rng.uniform(f_lo, f_hi)
    if abs(f) < 1e-4 or abs(e + f) < 1e-4: continue
    rts2 = np.roots([f, 3*e, 3*c, b])
    if np.max(np.abs(rts2.imag)) > 1e-6: continue
    rts2 = np.sort(rts2.real)
    # Check orientation: we need Q2(s_m) * P2'(s_m) >= 0 (same sign as Q1(r_m)*P1'(r_m) >= 0!)
    # Wait: check sign of P1'(r_m)*Q1(r_m)!
    # At r1: P1'(r1) = (r1-r2)(r1-r3), Q1(r1) = c1*(r1-r2)(r1-r3) => P1'(r1)*Q1(r1) = c1*(r1-r2)^2*(r1-r3)^2 >= 0!
    # Wait: for P2 with leading coeff f (could be NEGATIVE!), if f < 0, P2(y) = f*(y-s1)(y-s2)(y-s3) still has
    # P2'(s1)*Q2(s1) = f^2 * c1' * ... >= 0 ONLY IF we keep orientation W(P2, Q2) >= 0, i.e., P2' Q2 - P2 Q2' >= 0!
    # At a root s of P2, (P2_y Q2 - P2 Q2_y)(s) = P2'(s) Q2(s) >= 0 REGARDLESS OF SIGN OF f!!
    fp_lo, fp_hi = -1e100, 1e100
    ok_fp = True
    for s in rts2:
        dp = 3*f*s*s + 6*e*s + 3*c
        cst = bp + 3*cp*s + 3*ep*s*s
        coef = dp * (s**3)
        rhs = -dp * cst
        if abs(coef) < 1e-12:
            if rhs > 1e-10: ok_fp = False; break
        elif coef > 0: fp_lo = max(fp_lo, rhs / coef)
        else: fp_hi = min(fp_hi, rhs / coef)
    if not ok_fp or fp_lo > fp_hi + 1e-10 or fp_lo < -1e40 or fp_hi > 1e40: continue
    # Skip if all roots of P1 and P2 happen to be negative (we specifically want mixed-sign!)
    if max(r1, r2, r3) < 0 and max(rts2) < 0: continue
    valid_mixed += 1
    for fp in [fp_lo, fp_hi, 0.5*(fp_lo + fp_hi)]:
        ok_r, mval = cubic_roots_and_rayleigh(a+b, b+c, c+e, e+f, ap+bp, bp+cp, cp+ep, ep+fp)
        if ok_r is False or (ok_r is True and mval < -1e-5):
            fails_mixed += 1
            if first_fail is None:
                first_fail = (ok_r, mval, (round(r1,2), round(r2,2), round(r3,2)), round(f,3))
            break
print(f"  Mixed-sign overlapping real-stable cubic jets = {valid_mixed}: FAILS = {fails_mixed}, first_fail = {first_fail}")

print("\n=== 2. CHECK d=2 STAGE 3 SIGN AND POLE-GAP ===")
def apply_D_flt(m, p):
    res = np.zeros(len(p) + 1)
    for r, c in enumerate(p): res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def eval_W_flt(f, g, ts):
    fp = np.arange(1, len(f))*f[1:]; gp = np.arange(1, len(g))*g[1:]
    return np.polyval(fp[::-1], ts)*np.polyval(g[::-1], ts) - np.polyval(f[::-1], ts)*np.polyval(gp[::-1], ts)

ts_r = np.concatenate([-10**np.linspace(2, -4, 30), np.linspace(-0.99, -0.01, 60), [0.0], 10**np.linspace(-4, 2, 30)])
# Print one GLZ-like d=2 example to see exact signs of W!
# f2 = [1], f1 = [1, 3], f0 = [1, 7, 8] (roots -0.78, -0.16; f1 root -0.33 => f2 << f1 << f0!)
f2 = np.array([1.0]); f1 = np.array([1.0, 3.0]); f0 = np.array([1.0, 7.0, 8.0])
M5 = apply_D_flt(5, f1); Y5 = apply_D_flt(5, apply_D_flt(5, f2))
w1 = eval_W_flt(f0, M5, ts_r); w2 = eval_W_flt(M5, Y5, ts_r); w12 = eval_W_flt(f0, Y5, ts_r)
print(f"  Sanity d=2 (f2=1, f1=1+3t, f0=1+7t+8t^2, j=5): W_XM in [{w1.min():.2f}, {w1.max():.2f}], W_MY in [{w2.min():.2f}, {w2.max():.2f}]")
# AHA! At j=5, let's see why W_MY was positive or negative!
# Now run the d=2 pole-gap search properly with the actual sign of W!
sgn = 1.0 if w1.min() > 0 else -1.0
for j in [20, 40, 80]:
    max_drop = -1e9; valid_s3 = 0; best_cfg = None
    for _ in range(10000):
        # Need f1 root -1/a and f0 roots -1/b2 < -1/a < -1/b1 < 0, AND D_{0.75 j} f2 = 1 + 0.75 j t (root -1/(0.75 j))
        # For f1 < D_r f2, need -1/a < -1/r => a < r!!
        # THAT'S WHY valid=0 before: a was up to 1000 >> r!!
        r_cut = 0.75 * j
        a = rng.uniform(1.05, r_cut * 0.95)
        b2 = rng.uniform(1.01, a * 0.99)
        b1 = rng.uniform(a * 1.01, a * 20.0)
        f2 = np.array([1.0]); f1 = np.array([1.0, a]); f0 = np.array([1.0, b1+b2, b1*b2])
        Mr = apply_D_flt(r_cut, f1); Yr = apply_D_flt(r_cut, apply_D_flt(r_cut, f2))
        if np.min(sgn * eval_W_flt(f0, Mr, ts_r)) <= 0 or np.min(sgn * eval_W_flt(Mr, Yr, ts_r)) <= 0:
            continue
        Mj = apply_D_flt(j, f1); Yj = apply_D_flt(j, apply_D_flt(j, f2))
        w1 = eval_W_flt(f0, Mj, ts_r); w2 = eval_W_flt(Mj, Yj, ts_r); w12 = eval_W_flt(f0, Yj, ts_r)
        if np.min(sgn * w1) <= 0 or np.min(sgn * w2) <= 0: continue
        Jj = float(np.min(4*w1*w2 / (w12*w12 + 1e-300)))
        if Jj < 1.0: continue
        valid_s3 += 1
        Mj1 = apply_D_flt(j+1, f1); Yj1 = apply_D_flt(j+1, apply_D_flt(j+1, f2))
        w1p = eval_W_flt(f0, Mj1, ts_r); w2p = eval_W_flt(Mj1, Yj1, ts_r); w12p = eval_W_flt(f0, Yj1, ts_r)
        Jj1 = float(np.min(4*w1p*w2p / (w12p*w12p + 1e-300)))
        drop = Jj - Jj1
        if drop > max_drop:
            max_drop = drop; best_cfg = (round(b2,2), round(a,2), round(b1,2), round(Jj,5), round(Jj1,5))
    print(f"  j={j:2d} (r*<=0.75j, valid={valid_s3}): max(J(j)-J(j+1))={max_drop:+.6f} (j*drop={j*max_drop:+.4f}, j^2*drop={j*j*max_drop:+.3f}), cfg={best_cfg}")
