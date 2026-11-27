"""
Quick 2-second check:
On GLZ D-lifted family E_{k, i} = D_{k+1}^i h~_{k-1, i} (i = 0..k-1, deg 2k-3):
Wait: check the direction of interlacing!
At k=2: E_{2,0} = h~_{1,0} = 1+2t (root -1/2), E_{2,1} = D_3 h~_{1,1} = 1+3t (root -1/3).
So E_{2,0} << E_{2,1}, i.e. W(E_{2,0}, E_{2,1}) > 0!
Check on GLZ for k = 4..16:
Does W(E_{k, i}, E_{k, i+r}) > 0 on ALL of R for r = 1, 2, 3, 4?!
Also check W(E_{k, i}, E_{k, i+r}) on (-1, 0) and on R!
"""
import numpy as np

def apply_Dj(p, j):
    out = [0.0]*(len(p)+1)
    for m, c in enumerate(p):
        out[m] += (1+m)*c
        out[m+1] += (j+m)*c
    return out

def poly_add(a, b):
    n = max(len(a), len(b))
    out = [0.0]*n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return out

def poly_deriv(p):
    return [m*p[m] for m in range(1, len(p))] or [0.0]

def poly_mul(a, b):
    out = [0.0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out

def poly_sub(a, b):
    n = max(len(a), len(b))
    out = [0.0]*n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] -= b[i]
    return out

def wronskian(f, g):
    # f g' - f' g
    return poly_sub(poly_mul(f, poly_deriv(g)), poly_mul(poly_deriv(f), g))

def poly_eval(p, t):
    s = 0.0
    for c in reversed(p):
        s = s*t + c
    return s

h = [[1.0, 2.0], [1.0]]
ts = np.concatenate([np.linspace(-2.0, -1.01, 40), np.linspace(-0.999, -0.001, 300), np.linspace(0.0, 5.0, 60)])

for k in range(2, 13):
    j = k + 1
    # Build E_{k, i} = D_j^i h_{k-1, i} for i = 0..k-1
    E = []
    for i in range(k):
        p = h[i]
        for _ in range(i):
            p = apply_Dj(p, j)
        E.append(p)
    if k >= 4:
        res = {}
        for r in [1, 2, 3, 4]:
            if r >= k:
                continue
            min_ratio = 1e99
            neg_count = 0
            for i in range(k - r):
                w = wronskian(E[i], E[i + r])
                w1 = wronskian(E[i], E[i + 1])
                vals = [poly_eval(w, t) / max(1e-300, poly_eval(w1, t)) for t in ts]
                mv = min(vals)
                if mv < min_ratio:
                    min_ratio = mv
                if mv < 0:
                    neg_count += 1
            res[r] = (neg_count, min_ratio)
        print(f"k={k:2d}: " + ", ".join(f"r={r}: neg={v[0]} (min W_r/W_1={v[1]:+.4e})" for r, v in res.items()))
    # step to h_k
    hn = []
    for i in range(k + 1):
        u = [0.0]
        if i <= k - 1: u = poly_add(u, apply_Dj(h[i], j))
        if i >= 1: u = poly_add(u, h[i - 1])
        hn.append(apply_Dj(u, j))
    h = hn
