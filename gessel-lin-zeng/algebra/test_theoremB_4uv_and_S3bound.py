"""
1. Test if Theorem B holds for ALL u in (0, 1] and ALL 0 < x <= 4*u^2 (i.e. W(A, E) >= 0)!
   Check numerically on dense grid AND symbolically!
2. Test d=2 adversarial Stage 3 drop when BOTH gamma_{0.75 j} >= 1 AND gamma_j >= 1!
   Does N^2 * (1 - gamma_{j+1}/gamma_j) stay BOUNDED as j = 12, 24, 48, 96, 192?!
   ALSO test when gamma_{j-1} >= 1 + 2/j and gamma_j >= 1 + 2/j (two consecutive integers)!
"""
import math, random
import numpy as np

print("=== 1. TEST THEOREM B FOR ALL u in (0, 1] AND ALL 0 < x <= 4*u^2 (W(A,E) >= 0) ===", flush=True)
# E(u, s, x) = (1+u^2)^2 * [(1+x)^2 + 2*u*(1+x)*s + s^2 - 4*(1-u^2)*x] - [1 + 2*u*s + 4*u^2*x - x^2]^2
mn = 1e100
worst = None
for uv in np.linspace(0.005, 1.0, 400):
    for frac in np.linspace(1e-4, 1.0, 400):
        xv = frac * 4.0 * uv * uv
        a_c = (1.0 - uv*uv)**2
        b_c = 2.0*uv*((1.0+uv*uv)**2*(1.0+xv) - 2.0*(1.0 + 4.0*uv*uv*xv - xv*xv))
        s_min = 2.0*math.sqrt(xv)
        if a_c > 1e-12:
            sv = -b_c / (2.0*a_c)
            if sv > s_min: s_min = sv
        val = (1.0+uv*uv)**2*((1.0+xv)**2 + 2.0*uv*(1.0+xv)*s_min + s_min**2 - 4.0*(1.0-uv*uv)*xv) - (1.0 + 2.0*uv*s_min + 4.0*uv*uv*xv - xv*xv)**2
        if val < mn:
            mn = val; worst = (uv, xv, s_min)
print(f"  ALL u in (0,1], ALL 0 < x <= 4*u^2 (160,000 grid pts): min E = {mn:+.6e} at u={worst[0]:.4f}, x={worst[1]:.4f}, s={worst[2]:.4f}", flush=True)

# Also check ASYMMETRIC u != v with x <= 4*u*v!
random.seed(42)
mn_asym = 1e100
worst_a = None
for _ in range(300000):
    u = random.uniform(0.01, 1.0)
    v = random.uniform(0.01, 1.0)
    u0 = max(u, v)
    x = random.uniform(1e-4, 4.0 * u * v)
    ratio = 10 ** random.uniform(-3, 3)
    p = math.sqrt(x * ratio); q = math.sqrt(x / ratio)
    W1 = 1 + 2*u*p + p*p
    W2 = 1 + 2*v*q + q*q
    W3 = (1 + 2*u*p)*(1 + 2*v*q) - x*x
    val = (1 + u0*u0)**2 * W1 * W2 - W3*W3
    if val < mn_asym:
        mn_asym = val; worst_a = (u, v, p, q)
print(f"  ASYMMETRIC u,v in (0,1], 0 < pq <= 4*u*v (300,000 pts): min = {mn_asym:+.6e} at {tuple(round(z,4) for z in worst_a)}", flush=True)

print("\n=== 2. TEST d=2 STAGE 3 DROP WHEN gamma_{0.75 j} >= 1 AND gamma_j >= 1 ===", flush=True)
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
for j in [12, 24, 48, 96, 192]:
    r_cut = 0.75 * j
    N = j + 1
    best_sc = 0.0; best_info = None
    tested = 0
    for tr in range(40000):
        mode = tr % 3
        if mode == 0:
            a = random.uniform(1.0/r_cut + 1e-5, 0.999)
        elif mode == 1:
            a = 1.0/r_cut + 10**random.uniform(-4, -0.5)*(1 - 1.0/r_cut)
        else:
            a = random.uniform(1.0/r_cut + 1e-5, min(0.999, 5.0/r_cut))
        disc = (2 + a*r_cut)**2 - 4*(r_cut+1)*a
        if disc <= 0: continue
        rM1 = (-(2+a*r_cut) - math.sqrt(disc))/(2*(r_cut+1))
        rM2 = (-(2+a*r_cut) + math.sqrt(disc))/(2*(r_cut+1))
        b2 = random.uniform(-rM1 + 1e-6, 0.9999)
        b1 = random.uniform(-rM2 + 1e-6, a - 1e-6)
        if not (0 < b1 < b2 < 1): continue
        g_r = min_gamma_d2(a, b1, b2, r_cut)
        if g_r < 1.0: continue
        g1 = min_gamma_d2(a, b1, b2, j)
        if g1 < 1.0: continue
        tested += 1
        g2 = min_gamma_d2(a, b1, b2, j + 1)
        sc = (N**2)*(1.0 - g2/g1)
        if sc > best_sc: best_sc = sc; best_info = (a, b1, b2, g_r, g1, g2)
    print(f"  j={j:3d} (N={N:3d}, {tested:5d} valid): MAX N^2*(1-g_{{j+1}}/g_j) = {best_sc:8.4f} | g_r={best_info[3]:.4f}, g_j={best_info[4]:.4f} -> g_{{j+1}}={best_info[5]:.4f} (a={best_info[0]:.5f})", flush=True)
