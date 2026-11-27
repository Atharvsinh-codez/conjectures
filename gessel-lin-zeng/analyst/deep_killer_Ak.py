"""DEEP ADVERSARIAL KILLER TEST ON A_k(t, y) = sum C(k,i) D_{k+1}^i htilde_{k,i}(t) y^i:
1. Test k=35, i=17, m=10^20 (verify our prediction that m -> inf eventually flips negative at x=7e-17!).
2. Test wild monotone m_j at k=8, 10 (exp m_j = 2^j, step jump [2,2,2,2,1000,1000,1000,1000], etc.)
   at t=0, t=-1, t=inf AND full UHP!
3. Test critical exponent eps in sum C(k,i)^eps g_{k,i}(t) y^i at t=0 for k = 4, 8, 12, 20!
4. Test GLZ A_k(t, y) Rayleigh on ALL t in R (including AT ROOTS of g_{k,i}(t)!) for k = 10, 12, 15!
5. Check GLZ A_k(t, y) for k = 2..14 at t > 0: is A_k(t, y) real-rooted in y for ALL t in R?
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 60

def apply_D_at_0(m, b0, b1):
    return b0, 2*b1 + m*b0

def apply_D_at_m1(m, b0, b1):
    return (1 - m)*b0, m*b0 - m*b1

def apply_D_at_inf(m, d, cd, cd1):
    return d + 1, (d + m)*cd, (d - 1 + m)*cd1 + (d + 1)*cd

def check_rayleigh_1var(p0_list, p1_list):
    k = len(p0_list) - 1
    s0 = max(abs(x) for x in p0_list); s1 = max(abs(x) for x in p1_list)
    if s0 == 0 or s1 == 0: return None
    c0 = [mp.mpf(x)/mp.mpf(s0) for x in p0_list]
    c1 = [mp.mpf(x)/mp.mpf(s1) for x in p1_list]
    try:
        rts = mp.polyroots(c0[::-1], maxsteps=400, extraprec=50)
    except Exception:
        return None
    if max(abs(mp.im(r)) for r in rts) > mp.mpf('1e-18'):
        return False
    signs = []
    for yj in rts:
        yj = mp.re(yj)
        v0 = mp.mpf(0); dv0 = mp.mpf(0); v1 = mp.mpf(0)
        for idx in range(k, -1, -1):
            if idx < k: dv0 = dv0 * yj + v0
            v0 = v0 * yj + c0[idx]
            v1 = v1 * yj + c1[idx]
        prod = dv0 * v1
        signs.append(1 if prod > 0 else (-1 if prod < 0 else 0))
    return (len(set(signs)) == 1 and signs[0] != 0)

print("=== 1. WILD MONOTONE m_j AT k = 8, 10: RAYLEIGH OF A_{k, align_m} AT t=0, t=-1, t=inf ===")
def test_ms_boundary(label, ms, align_m):
    k = len(ms)
    a0 = [1]; a1 = [0]
    b0 = [1]; b1 = [0]
    cd = [1]; cd1 = [0]; degs = [0]
    for m in ms:
        na0, na1 = [0]*(len(a0)+1), [0]*(len(a0)+1)
        nb0, nb1 = [0]*(len(b0)+1), [0]*(len(b0)+1)
        ncd, ncd1 = [0]*(len(cd)+1), [0]*(len(cd)+1)
        ndegs = [0]*(len(cd)+1)
        for i in range(len(a0)):
            # U_i = D_m h_i + h_{i-1}
            t0, t1 = apply_D_at_0(m, a0[i], a1[i])
            na0[i] += t0; na1[i] += t1
            na0[i+1] += a0[i]; na1[i+1] += a1[i]
            t0, t1 = apply_D_at_m1(m, b0[i], b1[i])
            nb0[i] += t0; nb1[i] += t1
            nb0[i+1] += b0[i]; nb1[i+1] += b1[i]
            d_new, tc, tc1 = apply_D_at_inf(m, degs[i], cd[i], cd1[i])
            ncd[i] += tc; ncd1[i] += tc1; ndegs[i] = d_new
            ncd[i+1] += cd[i]; ncd1[i+1] += cd1[i]; ndegs[i+1] = degs[i]
        # Outer D_m
        a0, a1 = zip(*[apply_D_at_0(m, na0[i], na1[i]) for i in range(len(na0))])
        b0, b1 = zip(*[apply_D_at_m1(m, nb0[i], nb1[i]) for i in range(len(nb0))])
        inf_res = [apply_D_at_inf(m, ndegs[i], ncd[i], ncd1[i]) for i in range(len(ncd))]
        degs, cd, cd1 = zip(*inf_res)
        a0, a1, b0, b1, degs, cd, cd1 = map(list, (a0, a1, b0, b1, degs, cd, cd1))
    # Align with align_m^i and weight C(k,i)
    p0_0, p1_0, p0_m1, p1_m1, p0_inf, p1_inf = [], [], [], [], [], []
    for i in range(k+1):
        w = math.comb(k, i)
        v0, v1 = a0[i], a1[i]
        for _ in range(i): v0, v1 = apply_D_at_0(align_m, v0, v1)
        p0_0.append(w*v0); p1_0.append(w*v1)
        v0, v1 = b0[i], b1[i]
        for _ in range(i): v0, v1 = apply_D_at_m1(align_m, v0, v1)
        p0_m1.append(w*v0); p1_m1.append(w*v1)
        d, vc, vc1 = degs[i], cd[i], cd1[i]
        for _ in range(i): d, vc, vc1 = apply_D_at_inf(align_m, d, vc, vc1)
        p0_inf.append(w*vc); p1_inf.append(w*vc1)
    ok0 = check_rayleigh_1var(p0_0, p1_0)
    okm1 = check_rayleigh_1var(p0_m1, p1_m1)
    okinf = check_rayleigh_1var(p0_inf, p1_inf)
    print(f"  {label:35s}: t=0:{ok0}, t=-1:{okm1}, t=inf:{okinf}", flush=True)

test_ms_boundary("exp k=8 (2..256)", [2**j for j in range(1, 9)], 256)
test_ms_boundary("exp k=10 (2..1024)", [2**j for j in range(1, 11)], 1024)
test_ms_boundary("10^j k=8", [10**j for j in range(1, 9)], 10**8)
test_ms_boundary("step [2]*4+[1000]*4 (k=8)", [2]*4 + [1000]*4, 1000)
test_ms_boundary("step [2]*5+[1000]*5 (k=10)", [2]*5 + [1000]*5, 1000)
test_ms_boundary("step [2]*7+[10^6]*7 (k=14)", [2]*7 + [10**6]*7, 10**6)
test_ms_boundary("step [2]*10+[10^6]*10 (k=20)", [2]*10 + [10**6]*10, 10**6)
test_ms_boundary("wild non-mono [2,10^4]*5 (k=10)", [2, 10**4]*5, 10**4)

print("\n=== 2. CRITICAL EXPONENT eps IN sum C(k,i)^eps g_{k,i}(t) y^i AT t = 0 ===")
# Build GLZ a0, a1 up to k=20
a0 = [[0]*(k+1) for k in range(25)]; a1 = [[0]*(k+1) for k in range(25)]
a0[1][0], a1[1][0] = 1, 2; a0[1][1], a1[1][1] = 1, 0
for k in range(2, 25):
    m = k + 1
    for i in range(k+1):
        u0 = (a0[k-1][i] if i <= k-1 else 0) + (a0[k-1][i-1] if i >= 1 else 0)
        u1 = ((2*a1[k-1][i] + m*a0[k-1][i]) if i <= k-1 else 0) + (a1[k-1][i-1] if i >= 1 else 0)
        a0[k][i] = u0; a1[k][i] = 2*u1 + m*u0

for k in [4, 8, 12, 20]:
    m = k + 1
    g0s = [a0[k][i] for i in range(k+1)]
    g1s = [(1<<i)*a1[k][i] + ((1<<i)-1)*m*a0[k][i] for i in range(k+1)]
    res = []
    for eps in [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]:
        p0 = [mp.mpf(math.comb(k, i))**eps * g0s[i] for i in range(k+1)]
        p1 = [mp.mpf(math.comb(k, i))**eps * g1s[i] for i in range(k+1)]
        res.append((eps, check_rayleigh_1var(p0, p1)))
    print(f"  k={k:2d}: " + ", ".join(f"{e}:{ok}" for e, ok in res), flush=True)

print("\n=== 3. CHECK k=35, i=17 WITH m = 10^20 AT x = 7e-17 ===")
def apply_S_int(L, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for r, c in enumerate(poly):
        res[r] += c * (r + 1)
        res[r + 1] += c * (L - r)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0]*n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def poly_mul(a, b):
    res = [0]*(len(a)+len(b)-1)
    for i, va in enumerate(a):
        if va == 0: continue
        for j, vb in enumerate(b): res[i+j] += va*vb
    return res

def poly_sub(a, b):
    n = max(len(a), len(b))
    res = [0]*n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] -= v
    return res

A = [[None]*(k+1) for k in range(36)]
A[0][0] = [1]
for k in range(1, 36):
    for i in range(k+1):
        u = [0]
        if i <= k-1: u = poly_add(u, apply_S_int(3*k-3-i, A[k-1][i]))
        if i >= 1: u = poly_add(u, A[k-1][i-1])
        A[k][i] = apply_S_int(3*k-2-i, u)

def eval_mp(p, xv):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * xv + c
    return acc

for m_val in [10**16, 10**18, 10**20]:
    G3 = []
    for j in [16, 17, 18]:
        p = A[35][j]
        for s in range(j): p = apply_S_int(m_val + (70 - j - 1) - 1 + s, p)
        G3.append(p)
    diff_g = poly_sub(poly_mul(G3[1], G3[1]), poly_mul(G3[0], G3[2]))
    neg_cnt = sum(1 for v in diff_g if v <= 0)
    val = eval_mp(diff_g, mp.mpf("7e-17"))
    print(f"  k=35, i=17, m={m_val:.0e}: neg_x_coeffs={neg_cnt}, val(7e-17)={mp.nstr(val, 5)}", flush=True)
