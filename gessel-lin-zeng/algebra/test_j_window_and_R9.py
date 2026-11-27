"""
Fast 10-second checks:
1. Check exact max R = (pq)^2 for which Theorem B (Phi >= 4/(1+u^2)^2) holds for ALL u in [1/sqrt(2), 1] (i.e. eps in [0, 1]):
   Does it hold for R <= 4, R <= 9, or even R <= 16?!
2. Check the EXACT stability window [j_min(k), j_max(k)] of GLZ Quad (htilde_{k,i-2}, D_j htilde_{k,i-1}, D_j^2 htilde_{k,i})
   at k = 10, 20, 35 (including t = 1e-16..1e-1 where k=35 fails at j=inf)!
3. Check the EXACT Stage 3 drop on random triples that ALSO satisfy Theorem C's condition (r* <= alpha * j, e.g. h_1 << D_{0.7 j} h_2 and h_0 << D_{0.7 j} h_1)!
"""
import math, random
import numpy as np
import mpmath as mp
mp.mp.dps = 50

print("=== 1. EXACT R = (pq)^2 THRESHOLD FOR THEOREM B WHEN eps in [0, 1] (u in [1/sqrt(2), 1]) ===", flush=True)
# Recall E(u, s, x) = (1+u^2)^2 * [(1+x)^2 + 2*u*(1+x)*s + s^2 - 4*(1-u^2)*x] - [1 + 2*u*s + 4*u^2*x - x^2]^2
# where x = p*q = sqrt(R), s = p+q >= 2*sqrt(x).
for R_max in [1.0, 2.0, 4.0, 6.0, 9.0, 12.0, 16.0]:
    x_max = math.sqrt(R_max)
    mn = 1e100
    worst = None
    for uv in np.linspace(1.0/math.sqrt(2.0), 1.0, 200):
        for xv in np.linspace(1e-4, x_max, 200):
            # Minimize quadratic in s >= 2*sqrt(xv):
            a_c = (1 - uv*uv)**2
            b_c = 2*uv*((1+uv*uv)**2*(1+xv) - 2*(1 + 4*uv*uv*xv - xv*xv))
            s_min = 2*math.sqrt(xv)
            if a_c > 1e-12:
                sv = -b_c / (2*a_c)
                if sv > s_min: s_min = sv
            val = (1+uv*uv)**2*((1+xv)**2 + 2*uv*(1+xv)*s_min + s_min**2 - 4*(1-uv*uv)*xv) - (1 + 2*uv*s_min + 4*uv*uv*xv - xv*xv)**2
            if val < mn:
                mn = val; worst = (uv, xv, s_min)
    print(f"  R <= {R_max:4.1f} (x=pq <= {x_max:.3f}), u in [0.707, 1]: min E = {mn:+.6e} at u={worst[0]:.4f}, x={worst[1]:.4f}, s={worst[2]:.4f}", flush=True)

print("\n=== 2. GLZ QUAD STABILITY WINDOW [j_min(k), j_max(k)] AT k=10, 20, 35 ===", flush=True)
def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p
def poly_add(a, b):
    n = max(len(a), len(b))
    return poly_trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def poly_deriv(a):
    return [i * a[i] for i in range(1, len(a))] if len(a) > 1 else [0]
def apply_Dj(p, j):
    dp = poly_deriv(p)
    res = [0] * (len(p) + 2)
    for i, c in enumerate(p):
        res[i] += c; res[i+1] += j * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c; res[i+2] += c
    return poly_trim(res)
def poly_eval(p, x):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 36):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

tgrid = [mp.mpf(v) for v in [-10, -2, -1, -0.5, -0.1, 0, 1e-18, 1e-16, 1e-14, 1e-10, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.5, 2, 20]]
for k in [10, 20, 35]:
    i = k // 2 + 1 # middle index (note at k=35, i=16 gives h_{35,14}, h_{35,15}, h_{35,16}!)
    if k == 35: i = 16
    h0, h1, h2 = ht[k][i-2], ht[k][i-1], ht[k][i]
    print(f"  --- GLZ k={k}, i={i} (actual step uses j={k+2}) ---", flush=True)
    for j in [int(0.7*k), int(0.8*k), k, k+2, 2*k, 10*k, 1000*k, 10**10]:
        M = apply_Dj(h1, j)
        Y = apply_Dj(apply_Dj(h2, j), j)
        dX, dM, dY = poly_deriv(h0), poly_deriv(M), poly_deriv(Y)
        mn = mp.mpf('1e100'); wt = None; wneg = False
        for tv in tgrid:
            x, dx = poly_eval(h0, tv), poly_eval(dX, tv)
            mv, dm = poly_eval(M, tv), poly_eval(dM, tv)
            yv, dy = poly_eval(Y, tv), poly_eval(dY, tv)
            w1 = x*dm - dx*mv; w2 = mv*dy - dm*yv; w3 = x*dy - dx*yv
            if w1 <= 0 or w2 <= 0: wneg = True; break
            g = (4*w1*w2)/(w3*w3) if w3 != 0 else mp.mpf('1e100')
            if g < mn: mn = g; wt = float(tv)
        print(f"    j={j:10d}: min_gamma = {'W<=0' if wneg else mp.nstr(mn, 6)} at t={wt}", flush=True)

print("\n=== 3. RE-TEST d=2 ADVERSARIAL STAGE 3 WITH THEOREM C GAP CONDITION (r* <= 0.75 * j)! ===", flush=True)
# At d=2: h2=1, h1=t+a, h0=(t+b1)(t+b2).
# Require h1 << D_{r} h2 (i.e. a >= 1/(0.75*j) = 4/(3*j)) AND h0 << D_{0.75*j} h1!
def W_quad(A, B):
    return (A[0]*B[1] - A[1]*B[0], 2*(A[0]*B[2] - A[2]*B[0]), A[1]*B[2] - A[2]*B[1])
def min_gamma_d2(a, b1, b2, j):
    X = (b1*b2, b1+b2, 1.0)
    M = (a, 2.0 + a*j, j + 1.0)
    Y = (1.0, 3.0*j, j*(j + 1.0))
    w1, w2, w3 = W_quad(X, M), W_quad(M, Y), W_quad(X, Y)
    if w1[2] <= 0 or w1[1]**2 - 4*w1[0]*w1[2] >= 0: return -1.0
    if w2[2] <= 0 or w2[1]**2 - 4*w2[0]*w2[2] >= 0: return -1.0
    p1, p2, p3 = np.array([w1[2],w1[1],w1[0]]), np.array([w2[2],w2[1],w2[0]]), np.array([w3[2],w3[1],w3[0]])
    p12 = np.polymul(p1, p2)
    crit = np.polysub(np.polymul(np.polyder(p12), p3), 2.0*np.polymul(p12, np.polyder(p3)))
    cands = [z.real for z in np.roots(crit) if abs(z.imag) < 1e-7] + [-1e6, 0.0, 1e6]
    mn = (4.0*w1[2]*w2[2])/(w3[2]**2) if abs(w3[2]) > 1e-30 else 1e100
    for tv in cands:
        v1, v2, v3 = np.polyval(p1, tv), np.polyval(p2, tv), np.polyval(p3, tv)
        val = (4.0*v1*v2)/(v3*v3) if abs(v3) > 1e-30 else 1e100
        if val < mn: mn = val
    return mn

random.seed(9)
for j in [12, 24, 48, 96]:
    r_cut = 0.75 * j
    N = j + 1
    best_sc = 0.0; best_info = None
    for _ in range(30000):
        a = random.uniform(1.0/r_cut + 1e-4, 0.999)
        # Roots of D_{r_cut}(t+a) = (r_cut+1)t^2 + (2+a*r_cut)t + a:
        disc = (2 + a*r_cut)**2 - 4*(r_cut+1)*a
        if disc <= 0: continue
        rM1 = (-(2+a*r_cut) - math.sqrt(disc))/(2*(r_cut+1))
        rM2 = (-(2+a*r_cut) + math.sqrt(disc))/(2*(r_cut+1))
        b2 = random.uniform(-rM1 + 1e-5, 0.9999)
        b1 = random.uniform(-rM2 + 1e-5, a - 1e-5)
        if not (0 < b1 < b2 < 1): continue
        g1 = min_gamma_d2(a, b1, b2, j)
        if g1 < 1.0: continue
        g2 = min_gamma_d2(a, b1, b2, j + 1)
        sc = (N**2)*(1.0 - g2/g1)
        if sc > best_sc: best_sc = sc; best_info = (a, b1, b2, g1, g2)
    print(f"  d=2 with r*<=0.75*j, j={j:2d} (N={N:2d}): MAX N^2*(1 - g_{{j+1}}/g_j) = {best_sc:.4f} (g1={best_info[3]:.4f} -> g2={best_info[4]:.4f})", flush=True)
