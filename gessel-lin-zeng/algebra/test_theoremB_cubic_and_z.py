"""
Test 3 breakthrough algebraic questions:
1. Complete the symbolic proof of Theorem B: verify on u in (0,1], x in (0,1] that
   whenever d1 < 0, disc_s <= 0 (so E(u, s, x) >= 0 for ALL s >= 2*sqrt(x), 0 < u <= 1, 0 < x <= 1)!
   Also check exact upper bound on x = pq implied by Cubic real stability!
2. Test CUBIC Stage 1 propagation:
   Given 5-tuple (A,B,C,E,F) where BOTH Cubics (A,B,C,E) and (B,C,E,F) are real-stable,
   is (A+B, B+C, C+E, E+F) ALWAYS real-stable?!
   And check GLZ Cubic C_{k,i}(t,y) = A - 3By + 3Cy^2 - Ey^3 and C_alpha up to k=26!
3. Test Q_{k,q} (q=1/4, 1/8): is D_{k+1} Q_{k-1}(t, y) + z * y * Q_{k-1}(t, q*y) real-stable in (t,y)
   for ALL z > 0 (z in {1e-4, 1e-2, 1, 1e2, 1e4})?!
"""
import math, random
import sympy as sp
import numpy as np
import mpmath as mp
mp.mp.dps = 50

print("=== 1. VERIFY CASE 2 OF THEOREM B: d1 < 0 => disc_s <= 0 on (u, r) in (0,1]^2 ===", flush=True)
# Recall d1/2 = 2*r^4*u + r^2*u^5 - 6*r^2*u^3 + r^2*u + 2*r*u^4 - 4*r*u^2 + 2*r + u^5 + 2*u^3 - u
# and disc_s <= 0  <==>  (r^2 + 1 - 2*u^2)^2 <= (1 - u^2)^3
# Check on dense 200x200 grid of (u, r) in (0, 1]^2!
u_vals = np.linspace(0.001, 1.0, 400)
r_vals = np.linspace(0.001, 1.0, 400)
U, R = np.meshgrid(u_vals, r_vals)
X = R**2
d1_half = 2*R**4*U + R**2*U**5 - 6*R**2*U**3 + R**2*U + 2*R*U**4 - 4*R*U**2 + 2*R + U**5 + 2*U**3 - U
disc_factor = (X + 1 - 2*U**2)**2 - (1 - U**2)**3
mask = d1_half < 0
print(f"  Grid points with d1 < 0: {np.sum(mask)} / {mask.size}", flush=True)
if np.sum(mask) > 0:
    print(f"  On d1 < 0: max(disc_factor) = {np.max(disc_factor[mask]):+.6e} (NEED <= 0!)", flush=True)
# WAIT A SECOND! Check if c1_half < 0 already implies disc_factor <= 0!
# Actually check algebraic relation between d1_half and disc_factor!
# Let's check for which (u, r) in (0,1]^2 we have disc_factor > 0 AND d1_half < 0!
bad = (d1_half < 0) & (disc_factor > 0)
print(f"  Bad points (d1 < 0 AND disc_s > 0): {np.sum(bad)}", flush=True)

# Also check exact upper bound on x = r^2 from Cubic Rayleigh quartic at y = (p/q)^{1/2}:
# 3 - 4*(u+v)*r + (2 + 4*u*v)*r^2 - r^4 >= 0  for u,v in [u_min, 1]!
# Check max r^2 = p*q over ALL (u, v) in [u_min, 1]^2!
for umin in [0.0, 0.5, 0.8, 0.9, 0.95, 0.99]:
    # Max of r^2 happens at u = v = umin!
    # Solve 3 - 8*u*r + (2 + 4*u^2)*r^2 - r^4 = 0 for largest positive root r!
    rts = np.roots([-1, 0, 2 + 4*umin**2, -8*umin, 3])
    real_pos = sorted([z.real for z in rts if abs(z.imag) < 1e-8 and z.real > 0])
    # Check where poly >= 0! At r=0 val=3 > 0; at r=1 val=3-8u+2+4u^2-1 = 4(1-u)^2 >= 0!
    print(f"  Cubic Rayleigh bound (u,v >= {umin:.2f}): positive roots in r = {[round(z,4) for z in real_pos]} => max x=pq <= {real_pos[-1]**2:.4f}", flush=True)

print("\n=== 2. TEST CUBIC STAGE 1 PROPAGATION ON RANDOM 5-TUPLES (A,B,C,E,F) ===", flush=True)
# Test if (A,B,C,E) real-stable Cubic AND (B,C,E,F) real-stable Cubic => (A+B, B+C, C+E, E+F) real-stable Cubic!
# Even simpler: check BOTH (a) pointwise Rayleigh quartic (using Plucker!) AND (b) actual polynomials!
def poly_mul(a, b):
    res = [0.0]*(len(a)+len(b)-1)
    for ii, ai in enumerate(a):
        for jj, bj in enumerate(b): res[ii+jj] += ai*bj
    return res
def pfromr(rts):
    p = [1.0]
    for r in rts: p = poly_mul(p, [-r, 1.0])
    return p
def poly_deriv(a):
    return [i * a[i] for i in range(1, len(a))] if len(a) > 1 else [0]
def poly_add_f(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0.0) + (b[i] if i < len(b) else 0.0) for i in range(n)]

# Check cubic real stability via UHP t -> max Im(y) <= 0 (for A - 3By + 3Cy^2 - Ey^3)
uhp_pts = [complex(-5, 1e-4), complex(-1, 1e-4), complex(-0.5, 1e-4), complex(-0.1, 1e-4),
           complex(0, 1e-4), complex(0.1, 1e-4), complex(1, 1e-4), complex(0, 0.5), complex(-0.5, 0.5)]
def is_cubic_stable(A, B, C, E):
    for tv in uhp_pts:
        a = np.polyval(A[::-1], tv)
        b = np.polyval(B[::-1], tv)
        c = np.polyval(C[::-1], tv)
        e = np.polyval(E[::-1], tv)
        rts = np.roots([-e, 3*c, -3*b, a])
        if np.max(np.imag(rts)) > 1e-7:
            return False, float(np.max(np.imag(rts)))
    return True, 0.0

# Also test via Rayleigh quartic in Plucker variables (EXACT 5-tuple algebraic check!):
# For 5 polynomials (A,B,C,E,F) at a fixed t:
# Each polynomial is represented by a 2-vector (f(t), f'(t)) in R^2!
# WAIT A SECOND!! THAT'S SO SMART!!
# W(f, g) = det([f, g])! So ANY valid tuple of Wronskians at a fixed t is JUST 5 vectors v_0, v_1, v_2, v_3, v_4 in R^2
# in counterclockwise angular order!!
# NO NEED TO GUESS PLUCKER RELATIONS — JUST SAMPLE 5 VECTORS IN R^2!!
random.seed(123)
tested_cub = 0
fail_cub_add = 0
worst_cub = None
for _ in range(50000):
    if tested_cub >= 2000: break
    # Sample 5 angles 0 < th0 < th1 < th2 < th3 < th4 < pi and radii r0..r4 > 0
    ths = sorted([random.uniform(0.01, math.pi - 0.01) for _ in range(5)])
    rads = [10 ** random.uniform(-1, 1) for _ in range(5)]
    vecs = [(r*math.cos(th), r*math.sin(th)) for r, th in zip(rads, ths)]
    def W(u, v): return u[0]*v[1] - u[1]*v[0]
    def check_rayleigh_cubic(v0, v1, v2, v3):
        # Rayleigh quartic >= 0 for all y in R:
        # w01 - 2*w02*y + (3*w12 + w03)*y^2 - 2*w13*y^3 + w23*y^4 >= 0
        c0 = W(v0, v1)
        c1 = -2*W(v0, v2)
        c2 = 3*W(v1, v2) + W(v0, v3)
        c3 = -2*W(v1, v3)
        c4 = W(v2, v3)
        if c0 <= 0 or c4 <= 0: return False
        # Check no real roots of odd multiplicity (i.e. min over R >= 0):
        # Critical points from derivative 4*c4*y^3 + 3*c3*y^2 + 2*c2*y + c1 = 0
        crts = np.roots([4*c4, 3*c3, 2*c2, c1])
        for z in crts:
            if abs(z.imag) < 1e-7:
                y = z.real
                val = (((c4*y + c3)*y + c2)*y + c1)*y + c0
                if val < -1e-9 * max(c0, c4):
                    return False
        return True
    if not check_rayleigh_cubic(vecs[0], vecs[1], vecs[2], vecs[3]): continue
    if not check_rayleigh_cubic(vecs[1], vecs[2], vecs[3], vecs[4]): continue
    tested_cub += 1
    # Now add adjacent pairs: u0=v0+v1, u1=v1+v2, u2=v2+v3, u3=v3+v4!
    u = [(vecs[m][0]+vecs[m+1][0], vecs[m][1]+vecs[m+1][1]) for m in range(4)]
    if not check_rayleigh_cubic(u[0], u[1], u[2], u[3]):
        fail_cub_add += 1
        if worst_cub is None: worst_cub = (ths, rads)
print(f"  POINTWISE CUBIC RAYLEIGH ADDITION (2000 valid 5-tuples in R^2): fails = {fail_cub_add}/{tested_cub}!!", flush=True)

# WAIT A SECOND! ALSO TEST QUARTIC AND DEGREE-n REAL-STABLE ADDITION IN R^2!!
# What about Quad (4 vectors in R^2): if (v0,v1,v2) and (v1,v2,v3) satisfy Quad Rayleigh (with W(v0,v3) > 0, since th0<th3<pi!),
# does (v0+v1, v1+v2, v2+v3) ALWAYS satisfy Quad Rayleigh?!
tested_q = 0
fail_q_add = 0
for _ in range(10000):
    if tested_q >= 2000: break
    ths = sorted([random.uniform(0.01, math.pi - 0.01) for _ in range(4)])
    rads = [10 ** random.uniform(-1, 1) for _ in range(4)]
    vecs = [(r*math.cos(th), r*math.sin(th)) for r, th in zip(rads, ths)]
    def W(u, v): return u[0]*v[1] - u[1]*v[0]
    def check_quad(v0, v1, v2):
        return 4*W(v0,v1)*W(v1,v2) - W(v0,v2)**2 >= 0
    if not (check_quad(vecs[0], vecs[1], vecs[2]) and check_quad(vecs[1], vecs[2], vecs[3])): continue
    tested_q += 1
    u = [(vecs[m][0]+vecs[m+1][0], vecs[m][1]+vecs[m+1][1]) for m in range(3)]
    if not check_quad(u[0], u[1], u[2]):
        fail_q_add += 1
print(f"  POINTWISE QUAD RAYLEIGH ADDITION (with 0<th0<th3<pi, i.e. W_AE>0!): fails = {fail_q_add}/{tested_q}!!", flush=True)

print("\n=== 3. TEST Q_{k,q} (q=1/4, 1/8): IS D_{k+1} Q_{k-1}(t,y) + z * y * Q_{k-1}(t, q*y) STABLE FOR ALL z > 0? ===", flush=True)
def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p
def poly_add(a, b):
    n = max(len(a), len(b))
    return poly_trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def apply_Dj(p, j):
    dp = poly_deriv(p)
    res = [0] * (len(p) + 2)
    for i, c in enumerate(p):
        res[i] += c
        res[i+1] += j * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c
        res[i+2] += c
    return poly_trim(res)
def poly_eval_mp(p, x):
    acc = mp.mpc(0)
    for c in reversed(p): acc = acc * x + c
    return acc

ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 9):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

t_pts = [mp.mpc(-10, '1e-5'), mp.mpc(-1, '1e-5'), mp.mpc(-0.5, '1e-5'), mp.mpc(-0.01, '1e-5'),
         mp.mpc(0, '1e-5'), mp.mpc(0.01, '1e-5'), mp.mpc(1, '1e-5'), mp.mpc(0, 1)]
for q_str, q_val in [("1/2", mp.mpf('0.5')), ("1/4", mp.mpf('0.25')), ("1/8", mp.mpf('0.125'))]:
    for z_val in [mp.mpf('1e-4'), mp.mpf('0.01'), mp.mpf('1'), mp.mpf('100'), mp.mpf('1e4')]:
        max_im = mp.mpf('-1e100')
        for k in [2, 3, 5, 7]:
            D_ht = {i: apply_Dj(ht[k-1][i], k+1) for i in range(k)}
            for tv in t_pts:
                # Coeff of y^i in D_{k+1} Q_{k-1}(t,y) + z * y * Q_{k-1}(t, q*y):
                cy = []
                for i in range(k + 1):
                    v1 = poly_eval_mp(D_ht[i], tv) * (q_val**(i*(i-1)//2)) if i <= k-1 else mp.mpc(0)
                    v2 = z_val * poly_eval_mp(ht[k-1][i-1], tv) * (q_val**((i-1)*(i-2)//2)) * (q_val**(i-1)) if i >= 1 else mp.mpc(0)
                    cy.append(v1 + v2)
                rts = mp.polyroots(list(reversed(cy)), maxsteps=200, extraprec=30)
                for r in rts:
                    if r.imag > max_im: max_im = r.imag
        if max_im > mp.mpf('1e-12'):
            print(f"  FAIL q={q_str}, z={mp.nstr(z_val,3)}: max Im(y) = {mp.nstr(max_im, 6)}", flush=True)
        else:
            print(f"  PASS q={q_str}, z={mp.nstr(z_val,3)}: max Im(y) = {mp.nstr(max_im, 5)}", flush=True)
