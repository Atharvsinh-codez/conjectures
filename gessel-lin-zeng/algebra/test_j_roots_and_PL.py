"""
Test two brand-new exact algebraic structures:
1. Theorem C & j-roots of W(h0, D_j h1), W(D_j h1, D_j^2 h2), W(h0, D_j^2 h2) for GLZ:
   - What is the smallest r = r(k) such that htilde_{k, i-1} << D_r htilde_{k, i} (so W(f0, t f1)/W(f0, D_{k+1} f1) < 1/(k+1-r))?
   - What are the exact j-roots of W(h0, D_j h1) (deg 1), W(h0, D_j^2 h2) (deg 2), W(D_j h1, D_j^2 h2) (deg 3) at fixed t in R?
   - Do the j-roots of W(h0, D_j h1) * W(D_j h1, D_j^2 h2) INTERLACE or LIE TO THE RIGHT of the j-roots of W(h0, D_j^2 h2)?!
2. P_{L+1} intertwining route:
   - Does A_{k, i} << P_{3k-1-i} A_{k, i-1} hold through k = 16?
   - Does S_{L-1} A_{k-1, i} << P_{L+1} A_{k-1, i-2} (stronger Pair (2,1)) hold through k = 16?
   - Does f << P_{L+1}(S_L f) hold for ALL f in RZ(-inf, 0) with deg f < L (random test)?
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
def apply_SL(p, L):
    dp = poly_deriv(p)
    res = [0] * (len(p) + 2)
    for i, c in enumerate(p):
        res[i] += c; res[i+1] += L * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c; res[i+2] -= c
    return poly_trim(res)
def poly_eval(p, x):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 17):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

print("=== 1. THEOREM C: MAX j-ROOT r*(k) OF W(htilde_{k,i-1}, D_j htilde_{k,i}) = 0 OVER ALL t in R ===", flush=True)
# Note W(f0, D_j f1) = W(f0, D_0 f1) + j * W(f0, t f1) = 0 <=> j = -W(f0, D_0 f1) / W(f0, t f1)!
# Since W(f0, t f1) > 0 everywhere on R, r*(k) = max_{i, t in R} (-W(f0, D_0 f1) / W(f0, t f1))!
tgrid = [mp.mpf(v) for v in [-100, -10, -2, -1.1, -1, -0.95, -0.8, -0.5, -0.2, -0.05, 0, 0.05, 0.2, 0.5, 1, 5, 100]]
for k in [4, 6, 8, 10, 12, 14, 16]:
    max_r = mp.mpf('-1e100')
    min_r = mp.mpf('1e100')
    worst = None
    for i in range(1, k + 1):
        f0 = ht[k][i-1]
        f1 = ht[k][i]
        D0f1 = apply_Dj(f1, 0) # t(1+t)f1' + f1
        tf1 = [0] + f1
        df0, dD0, dtf1 = poly_deriv(f0), poly_deriv(D0f1), poly_deriv(tf1)
        for tv in tgrid:
            v0, dv0 = poly_eval(f0, tv), poly_eval(df0, tv)
            vD, dvD = poly_eval(D0f1, tv), poly_eval(dD0, tv)
            vt, dvt = poly_eval(tf1, tv), poly_eval(dtf1, tv)
            wD0 = v0*dvD - dv0*vD
            wt = v0*dvt - dv0*vt
            rval = -wD0 / wt
            if rval > max_r: max_r = rval; worst = (i, float(tv))
            if rval < min_r: min_r = rval
    print(f"  k={k:2d} (k+1={k+1:2d}): r*(k) = max j-root of W(h_{{i-1}}, D_j h_i) = {mp.nstr(max_r, 6)} (ratio r*/k = {mp.nstr(max_r/k, 5)}) at {worst}; min={mp.nstr(min_r, 5)}", flush=True)

print("\n=== 2. INSPECT ALL j-ROOTS OF W_XM(j), W_MY(j), W_XY(j) AT k=10, i=6 ===", flush=True)
k = 10; i = 6
h0, h1, h2 = ht[k][i-2], ht[k][i-1], ht[k][i]
# Expand D_j h1 = D0(h1) + j * (t h1)
# Expand D_j^2 h2 = D0^2(h2) + j * (2 t D0(h2) + t(1+t) h2) + j^2 * (t^2 h2)
D0h1 = apply_Dj(h1, 0); th1 = [0] + h1
D0h2 = apply_Dj(h2, 0); D02h2 = apply_Dj(D0h2, 0)
Y1 = poly_add([0] + [2*c for c in D0h2], poly_add([0]+h2, [0,0]+h2))
Y2 = [0, 0] + h2
for tv in [mp.mpf(v) for v in [-2, -0.5, 0.2, 0.5, 10]]:
    def evalW(f, g):
        return poly_eval(f, tv)*poly_eval(poly_deriv(g), tv) - poly_eval(poly_deriv(f), tv)*poly_eval(g, tv)
    # W_XM(j) = a0 + a1 j
    a0, a1 = evalW(h0, D0h1), evalW(h0, th1)
    # W_XY(j) = c0 + c1 j + c2 j^2
    c0, c1, c2 = evalW(h0, D02h2), evalW(h0, Y1), evalW(h0, Y2)
    # W_MY(j) = b0 + b1 j + b2 j^2 + b3 j^3
    b0 = evalW(D0h1, D02h2)
    b1 = evalW(D0h1, Y1) + evalW(th1, D02h2)
    b2 = evalW(D0h1, Y2) + evalW(th1, Y1)
    b3 = evalW(th1, Y2)
    r_XM = [-a0/a1]
    r_XY = mp.polyroots([c2, c1, c0])
    r_MY = mp.polyroots([b3, b2, b1, b0])
    # Compare sum of numerator j-roots vs 2 * sum of denominator j-roots!
    sum_num = (-a0/a1) + (-b2/b3)
    sum_den2 = 2 * (-c1/c2)
    print(f"  t={float(tv):5.2f}: r_XM={[round(float(z.real),2) for z in r_XM]}, r_MY={[complex(round(float(z.real),2), round(float(z.imag),2)) for z in r_MY]}, r_XY={[complex(round(float(z.real),2), round(float(z.imag),2)) for z in r_XY]} | sum_num - 2*sum_den = {float(sum_num - sum_den2):+.4f}", flush=True)

print("\n=== 3. TEST f << P_{L+1}(S_L f) FOR RANDOM f in RZ(-inf, 0) WITH deg f < L ===", flush=True)
random.seed(7)
fails_fPS = 0
for tr in range(300):
    d = random.randint(1, 5)
    L = random.randint(d + 1, d + 10)
    rts = sorted([-10**random.uniform(-2, 2) for _ in range(d)])
    f = [1.0]
    for r in rts:
        nf = [0.0]*(len(f)+1)
        for idx, c in enumerate(f):
            nf[idx] += -r*c; nf[idx+1] += c
        f = nf
    Sf = apply_SL(f, L)
    PSf = [((L+1) - m)*c for m, c in enumerate(Sf)]
    rPS = np.roots(list(reversed(PSf)))
    if np.max(np.abs(np.imag(rPS))) > 1e-6:
        fails_fPS += 1; continue
    rPS = np.sort(np.real(rPS))
    if len(rPS) != d + 1 or not all(rPS[m] < rts[m] < rPS[m+1] for m in range(d)):
        fails_fPS += 1
print(f"  Random f in RZ(-inf,0) (300 trials): f << P_{{L+1}}(S_L f) fails = {fails_fPS}/300", flush=True)

print("\n=== 4. TEST STRONGER PAIR (2,1): S_{L-1} A_{k-1,i} << P_{L+1} A_{k-1,i-2} FOR k=4..14 ===", flush=True)
A = {0: {0: [1]}}
for k in range(1, 15):
    A[k] = {}
    for idx in range(k+1):
        p1 = apply_SL(A[k-1][idx], 3*k-3-idx) if idx <= k-1 else [0]
        p2 = A[k-1][idx-1] if idx >= 1 else [0]
        A[k][idx] = apply_SL(poly_add(p1, p2), 3*k-2-idx)

for k in [4, 6, 8, 10, 12, 14]:
    ok1 = 0; ok2 = 0; tot1 = 0; tot2 = 0
    for idx in range(1, k+1):
        L = 3*k - 2 - idx
        # Check A_{k,i} << P_{L+1} A_{k,i-1}
        f = A[k][idx-1]; Pf = [((L+1)-m)*c for m,c in enumerate(f)]; g = A[k][idx]
        rP = sorted([r.real for r in mp.polyroots(list(reversed(Pf)), maxsteps=300, extraprec=30)])
        rg = sorted([r.real for r in mp.polyroots(list(reversed(g)), maxsteps=300, extraprec=30)]) if len(g)>1 else []
        if all(rP[m] < rg[m] < rP[m+1] for m in range(len(rg))): ok1 += 1
        tot1 += 1
        # Check S_{L-1} A_{k-1,i} << P_{L+1} A_{k-1,i-2} (for 2 <= i <= k-1)
        if 2 <= idx <= k-1:
            lhs = apply_SL(A[k-1][idx], L - 1)
            rhs = [((L+1)-m)*c for m,c in enumerate(A[k-1][idx-2])]
            rL = sorted([r.real for r in mp.polyroots(list(reversed(lhs)), maxsteps=300, extraprec=30)])
            rR = sorted([r.real for r in mp.polyroots(list(reversed(rhs)), maxsteps=300, extraprec=30)])
            if all(rR[m] < rL[m] < rR[m+1] for m in range(len(rL))): ok2 += 1
            tot2 += 1
    print(f"  k={k:2d}: A_{{k,i}} << P_{{L+1}} A_{{k,i-1}}: {ok1}/{tot1};  S_{{L-1}} A_{{k-1,i}} << P_{{L+1}} A_{{k-1,i-2}}: {ok2}/{tot2}", flush=True)
