"""
Test the SELF-CLOSING PLUCKER QUAD INDUCTION:
1. Exact Stage 1 inequality:
   Phi(u, v, p, q) = 4*(1 + 2*u*p + p^2)*(1 + 2*v*q + q^2) / ((1 + 2*u*p)*(1 + 2*v*q) - p^2*q^2)^2
   Check if Phi(u, v, p, q) >= 4 / (1 + max(u,v)^2)^2 for ALL u,v in (0,1] and p*q <= 1 (or p*q <= X_max)!
2. Stage 1 Plucker ratio x_out^2 = W(A+B,B+C)*W(C+E,E+F) / W(B+C,C+E)^2:
   Check max x_out^2 over all 5-tuples (A,B,C,E,F) satisfying u,v,w <= 1 and input Plucker ratios <= 1!
3. Does D_j preserve x^2 = W_AB W_CE / W_BC^2 <= 1 on random 4-tuples (or GLZ 4-tuples)?
   Track max_{t in R} x^2 across Stage 0 -> Stage 1 -> Stage 2 (D) -> Stage 3 (shift) for GLZ k=6..18!
4. Quick check of Pair (3,2) for k=4..14 and Family I/II with wild b_j (rescued from probe_4_directions.py).
"""
import math, random
import numpy as np
import mpmath as mp
mp.mp.dps = 50

print("=== 1. TEST STAGE 1 QUANTITATIVE INEQUALITY Phi(u,v,p,q) >= 4/(1 + max(u,v)^2)^2 ===", flush=True)
random.seed(42)
for pq_bound in [1.0, 1.5, 2.0, 2.5]:
    min_diff = 1e100
    worst = None
    for _ in range(200000):
        u = random.uniform(0.5, 0.999)
        v = u if random.random() < 0.3 else random.uniform(0.5, 0.999)
        umax = max(u, v)
        # sample p*q <= pq_bound
        pq = random.uniform(1e-4, pq_bound)
        ratio = 10 ** random.uniform(-3, 3)
        p = math.sqrt(pq * ratio)
        q = math.sqrt(pq / ratio)
        W1 = 1 + 2*u*p + p*p
        W2 = 1 + 2*v*q + q*q
        W3 = (1 + 2*u*p)*(1 + 2*v*q) - p*p*q*q
        if abs(W3) < 1e-12: continue
        phi = 4 * W1 * W2 / (W3 * W3)
        target = 4 / ((1 + umax*umax)**2)
        diff = phi - target
        if diff < min_diff:
            min_diff = diff
            worst = (u, v, p, q, phi, target)
    print(f"  p*q <= {pq_bound}: min(Phi - 4/(1+umax^2)^2) = {min_diff:+.6e}  (worst u={worst[0]:.4f}, v={worst[1]:.4f}, p={worst[2]:.4f}, q={worst[3]:.4f}, phi={worst[4]:.6f}, tgt={worst[5]:.6f})", flush=True)

print("\n=== 2. TEST STAGE 1 PLUCKER RATIO PROPAGATION ON 5-TUPLES (A,B,C,E,F) ===", flush=True)
# Normalize W_BC = 1, W_CE = c, W_AB = a, W_EF = f.
# Input Plucker ratios: x1^2 = W_AB*W_CE/W_BC^2 = a*c <= X1, x2^2 = W_BC*W_EF/W_CE^2 = f/c^2 <= X2!
# Wait: W_1 = W(A+B, B+C) = a + 2*u*sqrt(a) + 1
#       W_2 = W(B+C, C+E) = 1 + 2*v*sqrt(c) + c
#       W_3 = W(C+E, E+F) = c + 2*w*sqrt(c*f) + f
# So x_out^2 = W_1 * W_3 / W_2^2!
# Notice W_3 / c = 1 + 2*w*sqrt(f/c) + f/c!
# Let p = sqrt(a), q = sqrt(c), r = sqrt(f)/c (so x1 = p*q, x2 = r*q ... wait: x2 = sqrt(f)/c = r! No, f/c^2 = x2^2 => sqrt(f)/c = x2!)
# WAIT A SECOND! Look at W_1 * W_3 / W_2^2:
# W_1 = 1 + 2*u*p + p^2
# W_2 = 1 + 2*v*q + q^2
# W_3 = c*(1 + 2*w*sqrt(f/c) + f/c) = q^2 * (1 + 2*w*(x2*q) + (x2*q)^2) = q^2 + 2*w*x2*q^3 + x2^2*q^4!
# Substitute p = x1 / q!
# W_1 * W_3 / W_2^2 = (1 + 2*u*x1/q + x1^2/q^2) * q^2 * (1 + 2*w*x2*q + x2^2*q^2) / (1 + 2*v*q + q^2)^2
#                   = (q^2 + 2*u*x1*q + x1^2) * (1 + 2*w*x2*q + x2^2*q^2) / (1 + 2*v*q + q^2)^2 !!
# OH!! MY!! GOD!! CHECK THAT FORMULA WHEN x1 <= 1, x2 <= 1, u,v,w <= 1!!
# Wait: if u = v = w = 1 and x1 = x2 = 1:
# Numerator = (q+1)^2 * (1+q)^2 = (1+q)^4, Denominator = (1+q)^4 => RATIO = 1.000000 IDENTICALLY!!
for x_bound in [1.0, 1.5, 2.0]:
    max_xout = 0.0
    worst = None
    for _ in range(100000):
        u = random.uniform(0.5, 1.0)
        v = random.uniform(0.5, 1.0)
        w = random.uniform(0.5, 1.0)
        x1 = random.uniform(0.01, x_bound)
        x2 = random.uniform(0.01, x_bound)
        q = 10 ** random.uniform(-2, 2)
        num = (q*q + 2*u*x1*q + x1*x1) * (1 + 2*w*x2*q + x2*x2*q*q)
        den = (1 + 2*v*q + q*q)**2
        val = math.sqrt(num / den)
        if val > max_xout:
            max_xout = val
            worst = (u, v, w, x1, x2, q)
    print(f"  x1,x2 <= {x_bound} (u,v,w in [0.5,1]): max x_out = {max_xout:.6f} at {tuple(round(z,3) for z in worst)}", flush=True)

# ALSO check if u, v, w are ALL close to each other (e.g. u,w <= v + eps, or v >= 0.85, u,w <= 1):
for vmin in [0.9, 0.95, 0.99]:
    max_xout = 0.0
    for _ in range(50000):
        u = random.uniform(vmin, 1.0)
        v = random.uniform(vmin, 1.0)
        w = random.uniform(vmin, 1.0)
        x1 = random.uniform(0.01, 1.0)
        x2 = random.uniform(0.01, 1.0)
        q = 10 ** random.uniform(-2, 2)
        num = (q*q + 2*u*x1*q + x1*x1) * (1 + 2*w*x2*q + x2*x2*q*q)
        den = (1 + 2*v*q + q*q)**2
        val = math.sqrt(num / den)
        if val > max_xout: max_xout = val
    print(f"  x1,x2 <= 1.0, u,v,w in [{vmin}, 1]: max x_out = {max_xout:.6f}", flush=True)

print("\n=== 3. TRACK PLUCKER RATIO x^2 ACROSS STAGES 0->1->2->3 FOR GLZ k=8, 12, 16 ===", flush=True)
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
        res[i] += c
        res[i+1] += j * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c
        res[i+2] += c
    return poly_trim(res)
def poly_eval(p, x):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 18):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

grid_R = [mp.mpf(v) for v in [-100, -10, -2, -1.1, -1, -0.9, -0.7, -0.5, -0.3, -0.2, -0.1, -0.01, 0, 0.05, 0.2, 1, 10]]
def max_plucker(A, B, C, E):
    dA, dB, dC, dE = poly_deriv(A), poly_deriv(B), poly_deriv(C), poly_deriv(E)
    mx = mp.mpf(0)
    wt = None
    for tv in grid_R:
        a, da = poly_eval(A, tv), poly_eval(dA, tv)
        b, db = poly_eval(B, tv), poly_eval(dB, tv)
        c, dc = poly_eval(C, tv), poly_eval(dC, tv)
        e, de = poly_eval(E, tv), poly_eval(dE, tv)
        wAB = a*db - da*b
        wBC = b*dc - db*c
        wCE = c*de - dc*e
        r = (wAB*wCE)/(wBC*wBC)
        if r > mx:
            mx = r
            wt = float(tv)
    return mx, wt

for k in [8, 12, 16]:
    i = k // 2 + 1
    U = {r: poly_add(apply_Dj(ht[k-1][r], k+1) if r <= k-1 else [0], ht[k-1][r-1] if r >= 1 else [0]) for r in range(k+1)}
    def Dpow(p, j, m):
        for _ in range(m): p = apply_Dj(p, j)
        return p
    p0a = max_plucker(ht[k-1][i-3], Dpow(ht[k-1][i-2], k+1, 1), Dpow(ht[k-1][i-1], k+1, 2), Dpow(ht[k-1][i], k+1, 3))
    p1  = max_plucker(U[i-3], Dpow(U[i-2], k+1, 1), Dpow(U[i-1], k+1, 2), Dpow(U[i], k+1, 3))
    p2  = max_plucker(ht[k][i-3], Dpow(ht[k][i-2], k+1, 1), Dpow(ht[k][i-1], k+1, 2), Dpow(ht[k][i], k+1, 3))
    p3  = max_plucker(ht[k][i-3], Dpow(ht[k][i-2], k+2, 1), Dpow(ht[k][i-1], k+2, 2), Dpow(ht[k][i], k+2, 3))
    print(f"  k={k:2d}, i={i}: x^2 S0={mp.nstr(p0a[0],6)}@{p0a[1]} -> S1(add)={mp.nstr(p1[0],6)}@{p1[1]} -> S2(D)={mp.nstr(p2[0],6)}@{p2[1]} -> S3(shift)={mp.nstr(p3[0],6)}@{p3[1]}", flush=True)

print("\n=== 4. QUICK CHECK OF PAIR (3,2) FOR k=4..16 ===", flush=True)
for k in [4, 6, 8, 10, 12, 16]:
    fails = 0
    for idx in range(3, k):
        lhs = apply_Dj(apply_Dj(ht[k-1][idx], k+1), k+2)
        rhs = ht[k-1][idx-3]
        rL = sorted([r.real for r in mp.polyroots(list(reversed(lhs)), maxsteps=300, extraprec=30)])
        rR = sorted([r.real for r in mp.polyroots(list(reversed(rhs)), maxsteps=300, extraprec=30)])
        if not all(rR[m] < rL[m] < rR[m+1] for m in range(len(rL))):
            fails += 1
    print(f"  k={k:2d}: Pair (3,2) fails = {fails}/{k-3}", flush=True)
