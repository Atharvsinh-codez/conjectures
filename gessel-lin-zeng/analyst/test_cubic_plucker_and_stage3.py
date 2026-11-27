"""Two high-leverage tests:
1. CUBIC PLUCKER CONJECTURE:
   Suppose P1(t,y) = A + 3 B y + 3 C y^2 + E y^3 and P2(t,y) = B + 3 C y + 3 E y^2 + F y^3 are BOTH real-stable!
   (a) Univariate (constants a,b,c,e,f > 0): does real-rootedness of a+3by+3cy^2+ey^3 and b+3cy+3ey^2+fy^3
       ALWAYS imply real-rootedness of (a+b) + 3(b+c)y + 3(c+e)y^2 + (e+f)y^3?!
   (b) Bivariate (deg d=2,3 in t): does real stability of P1 and P2 ALWAYS imply real stability of
       P_new(t,y) = (A+B) + 3(B+C)y + 3(C+E)y^2 + (E+F)y^3?!
   (c) Check GLZ aligned cubic window W^{(3)}_{k, i}(t, y) = g_{k,i-3} + 3 g_{k,i-2} y + 3 g_{k,i-1} y^2 + g_{k,i} y^3
       for FULL UHP stability for k = 3..9!
2. STAGE 3 SHIFT TEST WITH POLE-GAP CONDITION r^* <= 3/4 * j:
   For d=2,3 polynomials (f0, f1, f2) in RZ(-1, 0) with f2 << f1 << f0, J(j) >= 1 + 2/j, AND
   r_{XM}, r_{MY} <= 0.75 * j (i.e. f0 < D_{0.75 j} f1 and D_{0.75 j} f1 < D_{0.75 j}^2 f2):
   Can J(j) - J(j+1) STILL be Omega(1/j) (> 1/j^2), or is it forced to be O(1/j^2)?!
"""
import numpy as np
from scipy.optimize import minimize

print("=== 1(a). UNIVARIATE CUBIC PLUCKER TEST (50,000 RANDOM TRIALS + BOUNDARY) ===")
def cubic_disc(p0, p1, p2, p3):
    # Discriminant of p0 + 3 p1 y + 3 p2 y^2 + p3 y^3 (divided by 27):
    # For a + 3b y + 3c y^2 + d y^3: disc/27 = 3 b^2 c^2 - 4 a c^3 - 4 b^3 d + 6 a b c d - a^2 d^2 >= 0!
    return 3*p1*p1*p2*p2 - 4*p0*p2**3 - 4*p1**3*p3 + 6*p0*p1*p2*p3 - p0*p0*p3*p3

rng = np.random.default_rng(2026)
uni_fails = 0; min_disc_rel = 1e9; worst_uni = None
for _ in range(50000):
    # Parametrize real-rooted cubic P1 by 3 positive roots r1, r2, r3 (i.e. roots -r1, -r2, -r3):
    r1, r2, r3 = 10 ** rng.uniform(-2, 2, 3)
    a = r1 * r2 * r3
    b = (r1*r2 + r1*r3 + r2*r3) / 3.0
    c = (r1 + r2 + r3) / 3.0
    e = 1.0
    # Now pick f > 0 such that b + 3 c y + 3 e y^2 + f y^3 is real-rooted!
    # Discriminant in f is quadratic: -b^2 f^2 + (6bce - 4c^3) f + (3c^2 e^2 - 4b e^3) >= 0!
    # Wait: for b + 3cy + 3ey^2 + fy^3, p0=b, p1=c, p2=e, p3=f:
    # A*f^2 + B*f + C >= 0 with A = -b^2 < 0, B = 6 b c e - 4 c^3, C = 3 c^2 e^2 - 4 b e^3.
    A_q = -b*b; B_q = 6*b*c*e - 4*c**3; C_q = 3*c*c*e*e - 4*b*e**3
    dd = B_q*B_q - 4*A_q*C_q
    if dd < 0: continue
    sd = np.sqrt(dd)
    f_lo = max(0.0, (-B_q + sd) / (2*A_q))  # since A_q < 0, roots are (-B +/- sd)/(2A), lo is +sd!
    f_hi = (-B_q - sd) / (2*A_q)
    if f_hi <= f_lo: continue
    for f in [f_lo + 1e-12*(f_hi-f_lo), f_hi - 1e-12*(f_hi-f_lo), rng.uniform(f_lo, f_hi)]:
        d_new = cubic_disc(a+b, b+c, c+e, e+f)
        scale = ((b+c)*(c+e))**2 + 1e-300
        if d_new < -1e-10 * scale:
            uni_fails += 1
            if uni_fails <= 3:
                print(f"  UNI CUBIC FAIL: a={a:.3g},b={b:.3g},c={c:.3g},e={e:.3g},f={f:.3g}, disc={d_new/scale:.3e}")
        if d_new / scale < min_disc_rel:
            min_disc_rel = d_new / scale; worst_uni = (a,b,c,e,f)
print(f"  Univariate Cubic Plucker fails: {uni_fails}, min rel disc = {min_disc_rel:+.3e}")

print("\n=== 1(b). GLZ ALIGNED CUBIC WINDOW W^{(3)}_{k,i}(t,y) = g_{k,i-3}+3g_{k,i-2}y+3g_{k,i-1}y^2+g_{k,i}y^3 UHP (k=3..9) ===")
def apply_D_int(m, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly): res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

H = [[None]*(k+1) for k in range(10)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, 10):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D_int(m, H[k-1][0])
    for i in range(1, k): U[i] = poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1): H[k][i] = apply_D_int(m, U[i])

test_ys = [complex(rx, ix) for rx in [-20, -3, -1, -0.5, 0, 0.5, 3, 20] for ix in [1e-4, 1e-2, 0.2, 2.0]]
for k in range(3, 10):
    m = k + 1
    G = []
    for i in range(k+1):
        p = list(H[k][i])
        for _ in range(i): p = apply_D_int(m, p)
        G.append(p)
    d = len(G[0]) - 1
    min_im, max_im = 1e9, -1e9
    for i in range(3, k + 1):
        g0, g1, g2, g3 = G[i-3], G[i-2], G[i-1], G[i]
        for yv in test_ys:
            c = np.array([float(g0[r]) + 3*yv*float(g1[r]) + 3*(yv**2)*float(g2[r]) + (yv**3)*float(g3[r]) for r in range(d+1)], dtype=np.complex128)
            if abs(c[-1]) < 1e-30: continue
            rts = np.roots(c[::-1])
            min_im = min(min_im, float(np.min(rts.imag)))
            max_im = max(max_im, float(np.max(rts.imag)))
    print(f"  k={k}: Cubic window UHP Im(t) in [{min_im:+.3e}, {max_im:+.3e}] (stable={min_im >= -1e-9})", flush=True)

print("\n=== 2. STAGE 3 SHIFT LOSS AT d=2 WITH POLE-GAP CONDITION r_{XM}, r_{MY} <= 0.75 * j ===")
# f2 = 1, f1 = 1 + a*t (root -1/a in (-1,0) => a > 1), f0 = (1 + b1*t)(1 + b2*t) with b1 > a > b2 > 1!
# Wait: roots of f0 are -1/b2 < -1/a < -1/b1 < 0 => b2 < a < b1!
def apply_D_flt(m, p):
    res = np.zeros(len(p) + 1)
    for r, c in enumerate(p): res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def eval_W_flt(f, g, ts):
    fp = np.arange(1, len(f))*f[1:]; gp = np.arange(1, len(g))*g[1:]
    return np.polyval(fp[::-1], ts)*np.polyval(g[::-1], ts) - np.polyval(f[::-1], ts)*np.polyval(gp[::-1], ts)

ts_r = np.concatenate([-10**np.linspace(2, -4, 30), np.linspace(-0.99, -0.01, 60), [0.0], 10**np.linspace(-4, 2, 30)])

for j in [20, 40, 80]:
    max_drop = -1e9; best_cfg = None
    # Random + targeted search over 1 < b2 < a < b1 with r_XM, r_MY <= 0.75 * j!
    for _ in range(3000):
        # Pick 1 < b2 < a < b1
        u = np.sort(1.0 + 10 ** rng.uniform(-1, 3, 3))
        b2, a, b1 = u[0], u[1], u[2]
        f2 = np.array([1.0])
        f1 = np.array([1.0, a])
        f0 = np.array([1.0, b1 + b2, b1 * b2])
        # Check pole gap at r = 0.75 * j: need W(f0, D_r f1) > 0 AND W(D_r f1, D_r^2 f2) > 0 on ALL t in R!
        r_cut = 0.75 * j
        Mr = apply_D_flt(r_cut, f1); Yr = apply_D_flt(r_cut, apply_D_flt(r_cut, f2))
        if np.min(eval_W_flt(f0, Mr, ts_r)) <= 0 or np.min(eval_W_flt(Mr, Yr, ts_r)) <= 0:
            continue
        # Compute J(j) and J(j+1)
        Mj = apply_D_flt(j, f1); Yj = apply_D_flt(j, apply_D_flt(j, f2))
        w1 = eval_W_flt(f0, Mj, ts_r); w2 = eval_W_flt(Mj, Yj, ts_r); w12 = eval_W_flt(f0, Yj, ts_r)
        Jj = float(np.min(4*w1*w2 / (w12*w12 + 1e-300)))
        if Jj < 1.0: continue
        Mj1 = apply_D_flt(j+1, f1); Yj1 = apply_D_flt(j+1, apply_D_flt(j+1, f2))
        w1p = eval_W_flt(f0, Mj1, ts_r); w2p = eval_W_flt(Mj1, Yj1, ts_r); w12p = eval_W_flt(f0, Yj1, ts_r)
        Jj1 = float(np.min(4*w1p*w2p / (w12p*w12p + 1e-300)))
        drop = Jj - Jj1
        if drop > max_drop:
            max_drop = drop; best_cfg = (b2, a, b1, Jj, Jj1)
    print(f"  j={j:2d} (r*<=0.75j): max(J(j)-J(j+1)) = {max_drop:+.6f} (j*drop={j*max_drop:+.4f}, j^2*drop={j*j*max_drop:+.3f}), cfg={best_cfg}")
