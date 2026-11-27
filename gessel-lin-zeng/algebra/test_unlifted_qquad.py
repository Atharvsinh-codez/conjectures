"""
CRITICAL 5-SECOND TEST OF THE STAGE-3-FREE UNLIFTED q-QUAD INDUCTION!
For unlifted polynomials E << C << B << A in (-1, 0) (degrees d-1, d, d+1, d+2):
  Q1^(q)(t, y) = A + 2*B*y + q*C*y^2  (real-stable <=> 4*W(A,B)*W(B,C) >= q * W(A,C)^2 on R!)
  Q2^(q)(t, y) = B + 2*C*y + q*E*y^2  (real-stable <=> 4*W(B,C)*W(C,E) >= q * W(B,E)^2 on R!)
Test:
1. For GLZ k = 4..35, what is min_{i, t in R} [4 W(h_{i-1}, h_i) W(h_i, h_{i+1}) / W(h_{i-1}, h_{i+1})^2]?
2. For random valid 4-tuples (A, B, C, E) satisfying Q1^(q) and Q2^(q) real-stable at q = 1/2, 1/4, 1/8, 1/16:
   Does Q_new^(q) := Q1^(q) + D_b Q2^(q) = (A + D_b B) + 2*(B + D_b C)*y + q*(C + D_b E)*y^2
   STAY REAL-STABLE (for b >= 2)?!
   Also test WITH the outer D_b applied: D_b(Q1^(q) + D_b Q2^(q))!
"""
import math, random
import numpy as np
import mpmath as mp
mp.mp.dps = 50

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
def poly_eval_mp(p, x):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

print("=== 1. GLZ UNLIFTED min_{t in R} 4*W(h_{i-1},h_i)*W(h_i,h_{i+1}) / W(h_{i-1},h_{i+1})^2 ===", flush=True)
ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 36):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

tgrid = [mp.mpf(v) for v in [-100, -10, -2, -1.1, -1, -0.95, -0.5, -0.1, -0.01, -1e-6, 0, 1e-16, 1e-6, 0.01, 0.1, 1, 10, 100]]
for k in [5, 10, 15, 20, 35]:
    dh = {i: poly_deriv(ht[k][i]) for i in range(k+1)}
    mn = mp.mpf('1e100'); wt = None
    for i in range(1, k):
        for tv in tgrid:
            f0, df0 = poly_eval_mp(ht[k][i-1], tv), poly_eval_mp(dh[i-1], tv)
            f1, df1 = poly_eval_mp(ht[k][i], tv),   poly_eval_mp(dh[i], tv)
            f2, df2 = poly_eval_mp(ht[k][i+1], tv), poly_eval_mp(dh[i+1], tv)
            w01 = f0*df1 - df0*f1; w12 = f1*df2 - df1*f2; w02 = f0*df2 - df0*f2
            # Note h_{i+1} << h_i << h_{i-1} => W(h_{i-1}, h_i) < 0, W(h_i, h_{i+1}) < 0, product > 0!
            r = (4*w01*w12)/(w02*w02) if w02 != 0 else mp.mpf('1e100')
            if r < mn: mn = r; wt = (i, float(tv))
    print(f"  GLZ k={k:2d}: min_R 4*W01*W12/W02^2 = {mp.nstr(mn, 7)} at {wt}", flush=True)

print("\n=== 2. TEST 1-STEP UNLIFTED q-QUAD PROPAGATION ON RANDOM 4-TUPLES (A, B, C, E) ===", flush=True)
def pfromr_np(rts):
    p = np.array([1.0])
    for r in rts: p = np.polymul(p, np.array([1.0, -r]))
    return p[::-1]
def apply_D_np(j, p):
    d = len(p) - 1
    out = np.zeros(d + 2)
    for m, cm in enumerate(p):
        out[m] += (m + 1)*cm; out[m+1] += (m + j)*cm
    return out
def wronsk_np(f, g):
    fh, gh = f[::-1], g[::-1]
    return np.polysub(np.polymul(fh, np.polyder(gh)), np.polymul(np.polyder(fh), gh))
def poly_add_np(a, b):
    n = max(len(a), len(b))
    out = np.zeros(n)
    out[:len(a)] += a; out[:len(b)] += b
    return out

ts = np.concatenate([np.linspace(-20, -1.01, 100), np.linspace(-0.999, -0.001, 250), np.linspace(0.001, 20, 100)])
def check_qquad(X, M, Y, q):
    # Check Y << M << X (so W(X, M) < 0, W(M, Y) < 0!) and 4 W(X,M) W(M,Y) >= q W(X,Y)^2 on R!
    w1 = np.polyval(wronsk_np(X, M), ts)
    w2 = np.polyval(wronsk_np(M, Y), ts)
    w3 = np.polyval(wronsk_np(X, Y), ts)
    if np.max(w1) >= 0 or np.max(w2) >= 0: return -1.0
    return float(np.min((4*w1*w2) / (w3*w3 + 1e-300)))

random.seed(42)
for q in [0.5, 0.25, 0.125, 0.0625, 0.01]:
    tested = 0; fails_add = 0; fails_int = 0; worst = 1e100
    for _ in range(5000):
        if tested >= 200: break
        # Generate E << C << B << A in (-1, 0) of degrees 1, 2, 3, 4
        rA = sorted([random.uniform(-0.98, -0.02) for _ in range(4)])
        rB = [random.uniform(rA[m]+1e-3*(rA[m+1]-rA[m]), rA[m+1]-1e-3*(rA[m+1]-rA[m])) for m in range(3)]
        rC = [random.uniform(rB[m]+1e-3*(rB[m+1]-rB[m]), rB[m+1]-1e-3*(rB[m+1]-rB[m])) for m in range(2)]
        rE = [random.uniform(rC[0]+1e-3*(rC[1]-rC[0]), rC[1]-1e-3*(rC[1]-rC[0]))]
        A = pfromr_np(rA)
        B = pfromr_np(rB) * (10**random.uniform(-1, 1))
        C = pfromr_np(rC) * (10**random.uniform(-1, 1))
        E = pfromr_np(rE) * (10**random.uniform(-1, 1))
        # Check Q1^(q) and Q2^(q) are real-stable (ratio >= q)
        r1 = check_qquad(A, B, C, q)
        if r1 < q * 1.01: continue
        r2 = check_qquad(B, C, E, q)
        if r2 < q * 1.01: continue
        tested += 1
        b = random.choice([2.0, 5.0, 15.0, 50.0])
        X_new = poly_add_np(A, apply_D_np(b, B))
        M_new = poly_add_np(B, apply_D_np(b, C))
        Y_new = poly_add_np(C, apply_D_np(b, E))
        rn = check_qquad(X_new, M_new, Y_new, q)
        if rn < 0: fails_int += 1
        if rn < q: fails_add += 1
        if rn < worst: worst = rn
    print(f"  q={q:6.4f} ({tested:3d} tested): fails_int={fails_int}, fails_qquad={fails_add}, min_ratio_out={worst:.5f} (need >= {q})", flush=True)
