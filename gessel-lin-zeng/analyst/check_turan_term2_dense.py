import mpmath as mp

mp.mp.dps = 80

def apply_D_int(m, p):
    res = [0]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v

KMAX = 10
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    H[k][0] = apply_D_int(m, apply_D_int(m, H[k-1][0]))
    for i in range(1, k):
        H[k][i] = apply_D_int(m, poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1]))
    H[k][k] = apply_D_int(m, H[k-1][k-1])

print("=== DENSE 500-POINT CHECK ON (-1, 0) OF Term2 = 2 h_{i-1} D h_i - h_{i-2} D h_{i+1} - h_i D h_{i-1} ===")
ts = [mp.mpf(-m) / mp.mpf(500) for m in range(1, 500)]
for k in [4, 6, 8, 10]:
    m = k + 1
    min_t2 = (mp.mpf('1e100'), None)
    neg_cnt = 0
    for i in range(2, k - 1):
        hm2, hm1, h0, hp1 = H[k-1][i-2], H[k-1][i-1], H[k-1][i], H[k-1][i+1]
        Dm1, D0, Dp1 = apply_D_int(m, hm1), apply_D_int(m, h0), apply_D_int(m, hp1)
        for tv in ts:
            a, b, c = eval_mp(hm2, tv), eval_mp(hm1, tv), eval_mp(h0, tv)
            da, db, dc = eval_mp(Dm1, tv), eval_mp(D0, tv), eval_mp(Dp1, tv)
            t2 = 2*b*db - a*dc - c*da
            sc = abs(2*b*db) + abs(a*dc) + abs(c*da) + mp.mpf('1e-300')
            rel = t2 / sc
            if rel < min_t2[0]: min_t2 = (rel, (i, float(tv)))
            if t2 < 0: neg_cnt += 1
    print(f"  k={k:2d}: min_rel_Term2 = {float(min_t2[0]):+.6f} at (i,t)={min_t2[1]}, neg_hits={neg_cnt}")
