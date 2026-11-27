"""
CRITICAL TEST: Are D_{k+1} Q_{k-1}(t, y) and Q_{k-1}(t, q*y) in PROPER POSITION in H_+^2?!
i.e., is Im( D_{k+1} Q_{k-1}(t, y) / Q_{k-1}(t, q*y) ) > 0 for ALL (t, y) in H_+^2?!
Equivalently: is D_{k+1} Q_{k-1}(t, y) + z * Q_{k-1}(t, q*y) REAL STABLE in ALL THREE variables (t, y, z)?!

If YES, setting z = y immediately gives real stability of D_{k+1} Q_{k-1}(t,y) + y Q_{k-1}(t,qy)
and applying D_{k+1} gives real stability of Q_k(t, y)!!

Also test:
- Exact range of q for k=2, 3, 4, 5, 6, 8, 10!
- Also test if D_b F(t, y) / F(t, q*y) has Im > 0 on H_+^2 for general products F(t, y) = prod (t + a_r + b_r y)
  or Hutchinson-separated real-stable F(t, y)!
"""
import math, random
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
        res[i] += c
        res[i+1] += j * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c
        res[i+2] += c
    return poly_trim(res)
def poly_eval(p, x):
    acc = mp.mpc(0)
    for c in reversed(p): acc = acc * x + c
    return acc

ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 11):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

random.seed(42)
# Build a VERY adversarial test set of (t, y) in H_+^2!
# Including tiny Im(t), tiny Im(y), large |t|, large |y|, near roots of Q_{k-1}(t, q*y)!
pts = []
for _ in range(400):
    tr = random.choice([-1000, -100, -10, -2, -1.01, -0.99, -0.5, -0.1, -0.01, 0, 0.01, 0.1, 1, 10, 100]) + random.uniform(-0.3, 0.3)
    ti = 10 ** random.uniform(-6, 2)
    yr = random.choice([-1000, -100, -10, -1, -0.1, -0.01, 0, 0.01, 0.1, 1, 10, 100, 1000]) + random.uniform(-0.3, 0.3)
    yi = 10 ** random.uniform(-6, 2)
    pts.append((mp.mpc(tr, ti), mp.mpc(yr, yi)))

print("=== 1. TEST Im( D_{k+1} Q_{k-1}(t,y) / Q_{k-1}(t, q*y) ) > 0 ON H_+^2 ===", flush=True)
for q_str, q_val in [("0.6", mp.mpf('0.6')), ("1/2", mp.mpf('0.5')), ("1/3", mp.mpf(1)/3), ("1/4", mp.mpf('0.25')), ("1/8", mp.mpf('0.125')), ("1/16", mp.mpf('0.0625'))]:
    print(f"  --- q = {q_str} ---", flush=True)
    for k in [2, 3, 4, 5, 6, 8, 10]:
        D_ht = {i: apply_Dj(ht[k-1][i], k+1) for i in range(k)}
        neg_cnt = 0
        min_im = mp.mpf('1e100')
        worst = None
        for tv, yv in pts:
            P = sum(poly_eval(D_ht[i], tv) * (q_val**(i*(i-1)//2)) * (yv**i) for i in range(k))
            S = sum(poly_eval(ht[k-1][i], tv) * (q_val**(i*(i-1)//2)) * ((q_val*yv)**i) for i in range(k))
            if S == 0: continue
            im = (P / S).imag
            if im < 0:
                neg_cnt += 1
            if im < min_im:
                min_im = im
                worst = (complex(tv), complex(yv), float(im))
        # ALSO test directly AT roots of Q_{k-1}(t, q*y) for real t = -10..10 (adversarial near-pole test!)
        pole_neg = 0
        for tr_val in [mp.mpf(v) for v in [-10, -2, -0.5, -0.05, 0, 0.05, 1, 10]]:
            # Find y-roots of S(tr_val, y) = 0, then nudge t -> tr_val + 1e-6j, y -> yr + 1e-6j!
            cy = [poly_eval(ht[k-1][i], tr_val) * (q_val**(i*(i-1)//2 + i)) for i in range(k)]
            if k <= 6:
                rts = mp.polyroots(list(reversed(cy)), maxsteps=200, extraprec=30)
                for rt in rts:
                    tv = tr_val + mp.mpc(0, '1e-6')
                    yv = rt.real + mp.mpc(0, '1e-6')
                    P = sum(poly_eval(D_ht[i], tv) * (q_val**(i*(i-1)//2)) * (yv**i) for i in range(k))
                    S = sum(poly_eval(ht[k-1][i], tv) * (q_val**(i*(i-1)//2)) * ((q_val*yv)**i) for i in range(k))
                    if S != 0 and (P / S).imag < 0:
                        pole_neg += 1
                        if (P/S).imag < min_im:
                            min_im = (P/S).imag
                            worst = (complex(tv), complex(yv), float((P/S).imag))
        print(f"    k={k:2d}: neg={neg_cnt}/400, pole_neg={pole_neg}, min_Im = {mp.nstr(min_im, 5)} {'' if neg_cnt+pole_neg==0 else 'WORST='+str((round(worst[0].real,3), round(worst[1].real,3)))}", flush=True)
