"""
Fast tests of:
1. Is D_{k+1} Q_{k-1}(t, y) + z * y * Q_{k-1}(t, q*y) real-stable in (t, y, z)
   (i.e. proper position: Im(D_{k+1} Q_{k-1}(t,y) / (y Q_{k-1}(t, q y))) has CONSTANT sign on H_+^2)?
   Test q = 1/2, 1/4, 1/8, 1/16 for k = 2, 3, 4, 5, 6!
2. Combinatorial refinements of JSP_{k,i} at k = 3, 4, 5:
   Test last_letter, last_base, first_base, end_dir, last_k_pos, first_k_pos for ALL-PAIRS interlacing!
"""
import math, random
from collections import defaultdict
import numpy as np
import mpmath as mp
mp.mp.dps = 40

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
for k in range(2, 9):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

print("=== 1. PROPER POSITION OF D_{k+1} Q_{k-1}(t,y) AND y*Q_{k-1}(t, q*y) ON H_+^2 ===", flush=True)
# Test 200 random (t, y) in H_+^2 (covering tiny/medium/large Im and Re on (-inf, inf))!
random.seed(5)
pts = []
for _ in range(300):
    tr = random.choice([-100, -10, -2, -1.01, -0.99, -0.5, -0.01, 0, 0.01, 1, 10]) + random.uniform(-0.5, 0.5)
    ti = 10 ** random.uniform(-4, 1)
    yr = random.choice([-100, -10, -1, -0.1, 0, 0.1, 1, 10, 100]) + random.uniform(-0.5, 0.5)
    yi = 10 ** random.uniform(-4, 1)
    pts.append((mp.mpc(tr, ti), mp.mpc(yr, yi)))

for q_str, q_val in [("1/2", mp.mpf('0.5')), ("1/4", mp.mpf('0.25')), ("1/8", mp.mpf('0.125')), ("1/16", mp.mpf('0.0625'))]:
    print(f"  --- q = {q_str} ---", flush=True)
    for k in [2, 3, 4, 5, 6]:
        # P(t,y) = D_{k+1} Q_{k-1}(t,y), R(t,y) = y * Q_{k-1}(t, q*y)
        D_ht = {i: apply_Dj(ht[k-1][i], k+1) for i in range(k)}
        pos_cnt = 0
        neg_cnt = 0
        max_pos = mp.mpf('-1e100')
        min_neg = mp.mpf('1e100')
        for tv, yv in pts:
            P = sum(poly_eval(D_ht[i], tv) * (q_val**(i*(i-1)//2)) * (yv**i) for i in range(k))
            R = yv * sum(poly_eval(ht[k-1][i], tv) * (q_val**(i*(i-1)//2)) * ((q_val*yv)**i) for i in range(k))
            if R == 0: continue
            im = (P / R).imag
            if im > 0:
                pos_cnt += 1
                if im > max_pos: max_pos = im
            else:
                neg_cnt += 1
                if im < min_neg: min_neg = im
        print(f"    k={k}: Im(D Q_{{k-1}} / (y Q_{{k-1}}(t,qy))) -> pos={pos_cnt}, neg={neg_cnt} (max_pos={mp.nstr(max_pos,4) if pos_cnt else 0})", flush=True)

print("\n=== 2. COMBINATORIAL REFINEMENTS OF JSP_{k,i} AT k=3, 4 (AND k=5, i=2) ===", flush=True)
# Generate JSP_{k,i} by iterative insertion of (bar_m, m, m) or (m, m)
# Represent letters as integers: bar_m = 2*m - 1, m = 2*m (so bar_1=1 < 1=2 < bar_2=3 < 2=4 < ...)
def gen_jsp(k):
    # returns dict i -> list of words (tuples)
    words = {0: [()]}
    for m in range(1, k + 1):
        bar_m = 2*m - 1
        unb_m = 2*m
        new_words = defaultdict(list)
        # Option 1 (m not in S): insert bar_m anywhere, then insert (unb_m, unb_m) together anywhere!
        for i, wlist in words.items():
            for w in wlist:
                L = len(w)
                # Insert bar_m at pos p1 in 0..L
                for p1 in range(L + 1):
                    w1 = w[:p1] + (bar_m,) + w[p1:]
                    # Insert (unb_m, unb_m) at pos p2 in 0..L+1
                    for p2 in range(L + 2):
                        w2 = w1[:p2] + (unb_m, unb_m) + w1[p2:]
                        new_words[i].append(w2)
                # Option 2 (m in S): insert (unb_m, unb_m) at pos p2 in 0..L
                for p2 in range(L + 1):
                    w2 = w[:p2] + (unb_m, unb_m) + w[p2:]
                    new_words[i + 1].append(w2)
        words = new_words
    return words

def des(w):
    return sum(1 for a, b in zip(w, w[1:]) if a > b)

def check_interlace_pair(p, q):
    # Check if p << q OR q << p (roots in (-inf, 0), deg diff <= 1)
    # Trim trailing zeros
    p = list(p); q = list(q)
    while len(p) > 1 and p[-1] == 0: p.pop()
    while len(q) > 1 and q[-1] == 0: q.pop()
    dp, dq = len(p)-1, len(q)-1
    if abs(dp - dq) > 1: return None
    rp = sorted(np.roots(p[::-1]).real) if dp >= 1 else []
    rq = sorted(np.roots(q[::-1]).real) if dq >= 1 else []
    if dp >= 1 and np.max(np.abs(np.roots(p[::-1]).imag)) > 1e-6: return "nonreal"
    if dq >= 1 and np.max(np.abs(np.roots(q[::-1]).imag)) > 1e-6: return "nonreal"
    # Check p << q (rightmost root of q >= rightmost root of p, or deg q = deg p + 1)
    def leq_interlace(a, b):
        # a << b: deg b == deg a => a[0]<=b[0]<=a[1]<=b[1]...; deg b == deg a + 1 => b[0]<=a[0]<=b[1]<=...<=b[-1]
        if len(b) == len(a):
            return all(a[m] <= b[m] + 1e-7 for m in range(len(a))) and all(b[m] <= a[m+1] + 1e-7 for m in range(len(a)-1))
        elif len(b) == len(a) + 1:
            return all(b[m] <= a[m] + 1e-7 and a[m] <= b[m+1] + 1e-7 for m in range(len(a)))
        return False
    p_le_q = leq_interlace(rp, rq)
    q_le_p = leq_interlace(rq, rp)
    if p_le_q and q_le_p: return "both"
    if p_le_q: return "<"
    if q_le_p: return ">"
    return "FAIL"

def test_stats_on(k, words):
    stats = {
        "last_letter": lambda w: w[-1],
        "last_base": lambda w: (w[-1] + 1) // 2,
        "first_base": lambda w: (w[0] + 1) // 2,
        "end_dir": lambda w: -1 if w[-2] > w[-1] else (0 if w[-2] == w[-1] else 1),
        "start_dir": lambda w: -1 if w[0] > w[1] else (0 if w[0] == w[1] else 1),
        "pos_last_k": lambda w: len(w) - 1 - w[::-1].index(2*k),
        "dist_last_k_end": lambda w: w[::-1].index(2*k),
        "pos_first_k": lambda w: w.index(2*k),
        "num_ltr_max": lambda w: sum(1 for idx in range(len(w)) if all(w[idx] > w[j] for j in range(idx))),
        "num_rtl_max": lambda w: sum(1 for idx in range(len(w)) if all(w[idx] > w[j] for j in range(idx+1, len(w)))),
        "num_plateaus": lambda w: sum(1 for a, b in zip(w, w[1:]) if a == b),
        "last_des_pos": lambda w: max([idx+1 for idx in range(len(w)-1) if w[idx] > w[idx+1]], default=0),
    }
    print(f"  --- k={k} ---", flush=True)
    for sname, sfn in stats.items():
        all_ok = True
        fail_info = None
        for i in range(1, k):
            buckets = defaultdict(lambda: defaultdict(int))
            for w in words[i]:
                buckets[sfn(w)][des(w)] += 1
            keys = sorted(buckets.keys())
            polys = []
            for key in keys:
                md = max(buckets[key].keys())
                p = [buckets[key][d] for d in range(md + 1)]
                polys.append((key, p))
            # Check ALL pairs!
            # Also check if there is a consistent total order (all pairs interlace!)
            for a in range(len(polys)):
                for b in range(a + 1, len(polys)):
                    rel = check_interlace_pair(polys[a][1], polys[b][1])
                    if rel in (None, "nonreal", "FAIL"):
                        all_ok = False
                        fail_info = (i, polys[a][0], polys[b][0], rel)
                        break
                if not all_ok: break
            if not all_ok: break
        print(f"    {sname:18s}: all-pairs interlace = {all_ok} {'' if all_ok else str(fail_info)}", flush=True)

w3 = gen_jsp(3)
test_stats_on(3, w3)
w4 = gen_jsp(4)
test_stats_on(4, w4)
