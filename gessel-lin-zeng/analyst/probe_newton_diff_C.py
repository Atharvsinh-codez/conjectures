"""Probe the Newton-difference basis C_{k, r}(t) = binom(k, r) Delta^r ghat_{k, 0}(t)
(equivalently sum_{i=0}^k g_{k,i}(t) y^i = sum_{r=0}^{k-1} C_{k, r}(t) y^r (1+y)^{k-r})!
Check for k = 2..10:
1. All r = 0..k-1: are ALL C_{k, r}(t) real-rooted in [-1, 0], with all positive coefficients?
2. What is the exact interlacing relation between C_{k, r}(t) and C_{k, r+1}(t)?
3. What is the exact interlacing relation between Delta^r ghat_{k, i} across i AND across r?
4. Check if C_{k, r}(t) / t (for r >= 1) is an INTERLACING SEQUENCE in r!
"""
import math
import mpmath as mp

mp.mp.dps = 80

def apply_D(m, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def poly_sub(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] -= v
    return res

KMAX = 10
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D(m, H[k-1][0])
    for i in range(1, k): U[i] = poly_add(apply_D(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1): H[k][i] = apply_D(m, U[i])

def roots_mp(p):
    lo = 0
    while lo < len(p) and p[lo] == 0: lo += 1
    q = p[lo:]
    while len(q) > 1 and q[-1] == 0: q.pop()
    d = len(q) - 1
    if d <= 0: return lo, [], 0.0
    g = (abs(mp.mpf(q[0])) / abs(mp.mpf(q[-1]))) ** (mp.mpf(1) / d)
    cs = [mp.mpf(q[r]) * (g ** r) for r in range(d + 1)]
    rts = mp.polyroots(cs[::-1], maxsteps=400, extraprec=50)
    max_im = float(max(abs(mp.im(r)) for r in rts))
    reals = sorted([mp.re(r) * g for r in rts])
    return lo, reals, max_im

def interlace_dir(r1, r2):
    if len(r1) == len(r2):
        d = len(r1)
        if all(r1[j] < r2[j] for j in range(d)) and all(r2[j] < r1[j+1] for j in range(d-1)): return "r1<r2"
        if all(r2[j] < r1[j] for j in range(d)) and all(r1[j] < r2[j+1] for j in range(d-1)): return "r2<r1"
        return "NO_eq"
    elif len(r1) == len(r2) + 1:
        d = len(r2)
        if all(r1[j] < r2[j] < r1[j+1] for j in range(d)): return "r2<<r1"
        return "NO_1"
    return f"deg({len(r1)},{len(r2)})"

for k in range(2, 10):
    m = k + 1; kfact = math.factorial(k)
    Ghat = []
    for i in range(k+1):
        p = list(H[k][i])
        for _ in range(i): p = apply_D(m, p)
        scale = kfact // math.comb(k, i)
        Ghat.append([c * scale for c in p])
    # Compute all forward differences Delta^r ghat_{k, 0} for r = 0..k
    cur = Ghat
    Cs = []
    for r in range(k + 1):
        Cs.append(cur[0])
        cur = [poly_sub(cur[j+1], cur[j]) for j in range(len(cur) - 1)]
    # Check each r = 0..k-1
    summary = []
    rts_list = []
    for r in range(k):
        p = Cs[r]
        all_pos = all(c >= 0 for c in p)
        lo, rts, mi = roots_mp(p)
        ok_rz = (mi < 1e-20 and (not rts or (rts[0] > -1 and rts[-1] < 0)))
        rts_list.append((lo, rts))
        summary.append(f"r={r}:(t^{lo},pos={all_pos},rz={ok_rz})")
    # Check interlacing between consecutive r and r+1 (for r >= 0)
    inters = []
    for r in range(k - 1):
        inters.append(f"{r}->{r+1}:{interlace_dir(rts_list[r][1], rts_list[r+1][1])}")
    print(f"k={k}: " + " ".join(summary))
    print(f"     interlace: " + ", ".join(inters), flush=True)
