"""
Three surgical tests (< 15s total):
1. Route I Link 2: Adversarial Nelder-Mead check of whether R_in <= R_cap (for R_cap = 2, 3, 4)
   with J_in >= 1 CLOSES under Stage 1 Pascal addition of 5 vectors (open half-plane AND without W_{04}>0)!
   Also check if Cubic Bezoutiant PSD (with J >= 1+eps) forces R <= 2 or R <= 4 under Pascal addition!
2. Route I Pole Shift Lemma:
   Given f1 << f0 in (-1, 0) and f0 << D_r f1 (0 < r < m),
   what is max (r_out - r) such that D_m f0 << D_{r_out} D_m f1?
   Wait: check if D_m D_r f1 << D_{r'} D_m f1 ALREADY holds universally for ALL f1 in RZ(-1,0)!
3. Route I Universal Stage 3 under (j*delta <= 3) AND (j^2*|delta_1 - delta_0| <= C):
   Adversarial optimization on d=2 and d=3 triples (f0, f1, f2) in RZ(-1,0) with:
     J_j(t) >= 1 + 2/j on R,
     j * delta_0(t), j * delta_1(t) <= 3.0 on R,
     j^2 * |delta_1(t) - delta_0(t)| <= 3.0 on R!
   What is the WORST possible j^2 * (ln J_{j+1} - ln J_j)?!
   If worst > -1.0 (beating Stage 1's +1/j^2 surplus!!), WE HAVE FOUND THE EXACT CLOSING INVARIANT!
"""
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(101)

# ============================================================
# TEST 1: ADVERSARIAL CHECK OF R <= R_cap CLOSURE UNDER PASCAL ADDITION
# ============================================================
print("=== TEST 1: ADVERSARIAL R_out FOR R_in <= R_cap (J_in >= 1) UNDER PASCAL ADDITION ===")
def W2(a, b): return a[0]*b[1] - a[1]*b[0]

def get_JR(vs):
    w01 = W2(vs[0], vs[1]); w12 = W2(vs[1], vs[2]); w23 = W2(vs[2], vs[3])
    w02 = W2(vs[0], vs[2]); w13 = W2(vs[1], vs[3])
    if w01 <= 0 or w12 <= 0 or w23 <= 0 or w02 <= 0 or w13 <= 0:
        return None
    J1 = 4*w01*w12 / (w02**2)
    J2 = 4*w12*w23 / (w13**2)
    R = (w01 * w23) / (w12**2)
    return min(J1, J2), R

for R_cap in [1.0, 2.0, 3.0, 4.0]:
    best_Rout = 0.0
    for seed in range(25):
        rng2 = np.random.default_rng(600 + seed)
        found = False
        for _ in range(1000):
            ang = np.sort(rng2.uniform(0.05, np.pi - 0.05, 5))
            lens = np.exp(rng2.normal(0, 1.5, 5))
            vs = [(lens[i]*np.cos(ang[i]), lens[i]*np.sin(ang[i])) for i in range(5)]
            r1 = get_JR(vs[:4]); r2 = get_JR(vs[1:])
            if r1 and r2 and r1[0] >= 1.0 and r2[0] >= 1.0 and r1[1] <= R_cap and r2[1] <= R_cap:
                found = True
                break
        if not found: continue
        x0 = np.concatenate([np.log(np.diff(ang)), np.log(lens[[0,1,3,4]] / lens[2])])
        def objR(z):
            g = np.exp(np.clip(z[:4], -12, 4))
            tot = np.sum(g)
            if tot >= np.pi - 1e-4: g = g * (np.pi - 0.01) / (tot + 0.01)
            th = np.zeros(5); th[0] = 0.005; th[1:] = th[0] + np.cumsum(g)
            L = np.exp(np.clip([z[4], z[5], 0.0, z[6], z[7]], -10, 10))
            vs = [(L[i]*np.cos(th[i]), L[i]*np.sin(th[i])) for i in range(5)]
            r1 = get_JR(vs[:4]); r2 = get_JR(vs[1:])
            if r1 is None or r2 is None: return 10.0
            pen = 0.0
            if r1[0] < 1.0: pen += 200*(1.0 - r1[0])
            if r2[0] < 1.0: pen += 200*(1.0 - r2[0])
            if r1[1] > R_cap: pen += 200*(r1[1] - R_cap)
            if r2[1] > R_cap: pen += 200*(r2[1] - R_cap)
            if pen > 0: return 10.0 + pen
            ws = [(vs[i][0]+vs[i+1][0], vs[i][1]+vs[i+1][1]) for i in range(4)]
            ro = get_JR(ws)
            return -ro[1] if ro else 10.0
        res = minimize(objR, x0, method='Nelder-Mead', options={'maxiter': 1500})
        if -res.fun > best_Rout:
            best_Rout = -res.fun
    print(f"  Adversarial R_in <= {R_cap:.1f} (J_in >= 1, open half-plane): max R_out = {best_Rout:.6f} ({'CLOSES!' if best_Rout <= R_cap else 'EXCEEDS'})")

# ============================================================
# TEST 2: UNIVERSAL D_m COMMUTATOR POLE SHIFT: D_m D_r f << D_{r'} D_m f
# ============================================================
print("\n=== TEST 2: Smallest shift (r' - r) so D_m D_r f << D_{r'} D_m f for random f in RZ(-1,0) ===")
def apply_Dj_f(p, j):
    out = np.zeros(len(p) + 1)
    for m, c in enumerate(p):
        out[m] += (1 + m) * c
        out[m + 1] += (j + m) * c
    return out

def roots_f(p):
    r = np.roots(p[::-1])
    if np.max(np.abs(r.imag)) > 1e-6: return None
    return np.sort(r.real)

def eq_interlaces(r1, r2):
    # equal degree, r1 << r2: r1[0] < r2[0] < r1[1] < ... < r1[d-1] < r2[d-1]
    if r1 is None or r2 is None or len(r1) != len(r2): return False
    return all(r1[q] < r2[q] for q in range(len(r1))) and all(r2[q] < r1[q+1] for q in range(len(r1)-1))

for d, m, r in [(2, 10, 6.0), (4, 15, 10.0), (6, 20, 13.0)]:
    max_shift = -1e9
    for _ in range(300):
        rts = -rng.uniform(0.01, 0.99, d)
        f = np.poly(rts)[::-1]
        lhs = roots_f(apply_Dj_f(apply_Dj_f(f, r), m))
        Dmf = apply_Dj_f(f, m)
        # bisect shift s in [-1, 5] such that D_m D_r f << D_{r+s} D_m f
        lo, hi = 0.0, 10.0
        if not eq_interlaces(lhs, roots_f(apply_Dj_f(Dmf, r + hi))):
            max_shift = 999.0
            break
        for _ in range(30):
            mid = 0.5*(lo + hi)
            if eq_interlaces(lhs, roots_f(apply_Dj_f(Dmf, r + mid))):
                hi = mid
            else:
                lo = mid
        if hi > max_shift:
            max_shift = hi
    print(f"  d={d}, m={m}, r={r}: max shift (r' - r) for D_m D_r f << D_{{r'}} D_m f = {max_shift:.5f} (compare (m-r)/m = {(m-r)/m:.4f})")

# ============================================================
# TEST 3: UNIVERSAL STAGE 3 ON d=2 AND d=3 WITH j*delta <= 3 AND j^2*|delta1-delta0| <= 3!
# ============================================================
print("\n=== TEST 3: ADVERSARIAL d=2 STAGE 3 UNDER j*delta <= 3 AND j^2*|delta_1 - delta_0| <= C ===")
# For d=2: f2 = 1, f1 = t + a, f0 = (t + b1)(t + b2) with 0 < b1 < a < b2 < 1.
# delta_0(t) = W(f0, t f1) / W(f0, D_j f1)
# delta_1(t) = W(f1, t f2) / W(f1, D_j f2)
def eval_d2_all(a, b1, b2, j, t):
    X  = b1*b2 + (b1+b2)*t + t*t; Xp = (b1+b2) + 2*t
    f1 = a + t; f1p = 1.0
    M  = a + (2 + a*j)*t + (j + 1)*t*t; Mp = (2 + a*j) + 2*(j + 1)*t
    Y  = 1 + 3*j*t + j*(j + 1)*t*t; Yp = 3*j + 2*j*(j + 1)*t
    tf1 = a*t + t*t; tf1p = a + 2*t
    tf2 = t; tf2p = 1.0
    Df2 = 1 + j*t; Df2p = j
    wXM = X*Mp - Xp*M; wMY = M*Yp - Mp*Y; wXY = X*Yp - Xp*Y
    J = 4*wXM*wMY / (wXY**2)
    d0 = (X*tf1p - Xp*tf1) / wXM
    d1 = (f1*tf2p - f1p*tf2) / (f1*Df2p - f1p*Df2)
    return J, wXM, wMY, wXY, d0, d1

ts_g = np.concatenate([np.linspace(-5.0, -1.05, 15), np.linspace(-0.99, -0.005, 40), np.linspace(0.0, 5.0, 20), [20.0, 100.0]])
for C_lip in [1.0, 3.0, 5.0, 1000.0]:
    for j in [40, 120]:
        worst_dlnJ = 1e9
        valid = 0
        for _ in range(15000):
            r = np.sort(np.exp(rng.uniform(-7, 0, 3)))
            if rng.random() < 0.5: r = 1.0 - r[::-1]
            b1, a, b2 = r[0], r[1], r[2]
            ok = True
            for tc in ts_g:
                J, w1, w2, w3, d0, d1 = eval_d2_all(a, b1, b2, j, tc)
                if w1 <= 0 or w2 <= 0 or w3 <= 0 or J < 1.0:
                    ok = False; break
                if j*d0 > 3.0 or j*d1 > 3.0 or (j**2)*abs(d1 - d0) > C_lip:
                    ok = False; break
            if not ok: continue
            valid += 1
            # compute min j^2 * ln(J_{j+1} / J_j) over ts_g
            md = min((j**2)*np.log(eval_d2_all(a,b1,b2,j+1,tc)[0] / eval_d2_all(a,b1,b2,j,tc)[0]) for tc in ts_g)
            if md < worst_dlnJ:
                worst_dlnJ = md
        print(f"  C_lip={C_lip:6.1f}, j={j:3d} (valid={valid:4d}): worst j^2 * ln(J_{{j+1}}/J_j) = {worst_dlnJ:+.4f}")
