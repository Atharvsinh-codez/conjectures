import mpmath as mp

# 1. Verify k = 72 at dps = 300 (checking if k=72 t=-0.99 was just 80-dps precision loss!)
# 2. Check the 3 terms of plain Turan u_i^2 - u_{i-1} u_{i+1} on GLZ for k = 10, 20, 35, 50:
#      Term1 = (D h_i)^2 - (D h_{i-1})(D h_{i+1})
#      Term2 = 2 h_{i-1} D h_i - h_{i-2} D h_{i+1} - h_i D h_{i-1}
#      Term3 = h_{i-1}^2 - h_{i-2} h_i
# 3. Also check: is A_{k,i}(x)^2 - A_{k,i-1}(x) A_{k,i+1}(x) >= 0 on ALL of R (incl. x > 0)?!

mp.mp.dps = 300

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

def deriv_int(p): return [r*p[r] for r in range(1, len(p))] if len(p) > 1 else [0]
def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v

KMAX = 72
row = [[1, 2], [1]]

for k in range(2, KMAX + 1):
    m = k + 1
    prev_row = row
    new_row = [None]*(k + 1)
    new_row[0] = apply_D_int(m, apply_D_int(m, row[0]))
    for i in range(1, k):
        new_row[i] = apply_D_int(m, poly_add_int(apply_D_int(m, row[i]), row[i-1]))
    new_row[k] = apply_D_int(m, row[k-1])
    row = new_row

    if k in [10, 25, 45, 72]:
        # Check k=72 at t=-0.99, -0.9999, -0.5, -0.01, -1e-6, 0, AND x>0 (t > 0 or t < -1)!
        ts_in = [mp.mpf(s) for s in ["-0.9999", "-0.99", "-0.95", "-0.8", "-0.5", "-0.2", "-0.05", "-1e-3", "-1e-6", "0"]]
        ts_out = [mp.mpf(s) for s in ["-10", "-2", "-1.01", "0.01", "1", "10"]]
        is_chk = sorted(set([2, 3, k//2, k-2]))
        is_chk = [i for i in is_chk if 2 <= i <= k - 2]
        min_tur_in = mp.mpf('1e100'); min_tur_out = mp.mpf('1e100')
        min_t1 = mp.mpf('1e100'); min_t2 = mp.mpf('1e100'); min_J_left = mp.mpf('1e100')
        min_Dtur = mp.mpf('1e100')
        for i in is_chk:
            D_next = {r: apply_D_int(k+2, row[r]) for r in [i-1, i, i+1]}
            d_cur  = {r: deriv_int(row[r]) for r in [i-1, i, i+1]}
            # prev_row for Term 1, Term 2 of u_i = D_{k+1} prev[i] + prev[i-1]
            Dp = {r: apply_D_int(k+1, prev_row[r]) for r in [i-1, i, i+1]}
            for tv in ts_in:
                v0, v1, v2 = eval_mp(row[i-1], tv), eval_mp(row[i], tv), eval_mp(row[i+1], tv)
                sc = v1*v1 + abs(v0*v2) + mp.mpf('1e-300')
                r_tur = (v1*v1 - v0*v2) / sc
                if r_tur < min_tur_in: min_tur_in = r_tur

                dv0, dv1, dv2 = eval_mp(D_next[i-1], tv), eval_mp(D_next[i], tv), eval_mp(D_next[i+1], tv)
                r_Dtur = (dv1*dv1 - dv0*dv2) / (dv1*dv1 + abs(dv0*dv2) + mp.mpf('1e-300'))
                if r_Dtur < min_Dtur: min_Dtur = r_Dtur

                if tv <= -mp.mpf(1)/(k+2):
                    p0, p1, p2 = eval_mp(d_cur[i-1], tv), eval_mp(d_cur[i], tv), eval_mp(d_cur[i+1], tv)
                    w1 = v1*p0 - p1*v0; w2 = v2*p1 - p2*v1; w12 = v2*p0 - p2*v0
                    J = 4*w1*w2/(w12*w12) if w12 != 0 else mp.mpf('1e10')
                    if J < min_J_left: min_J_left = J

                # Check Term 1 and Term 2 of u_i^2 - u_{i-1} u_{i+1}
                hm2, hm1, h0 = eval_mp(prev_row[i-2], tv), eval_mp(prev_row[i-1], tv), eval_mp(prev_row[i], tv)
                dm1, d0, dp1 = eval_mp(Dp[i-1], tv), eval_mp(Dp[i], tv), eval_mp(Dp[i+1], tv)
                t1 = d0*d0 - dm1*dp1
                t2 = 2*hm1*d0 - hm2*dp1 - h0*dm1
                sc1 = d0*d0 + abs(dm1*dp1) + mp.mpf('1e-300')
                sc2 = abs(2*hm1*d0) + abs(hm2*dp1) + abs(h0*dm1) + mp.mpf('1e-300')
                if t1/sc1 < min_t1: min_t1 = t1/sc1
                if t2/sc2 < min_t2: min_t2 = t2/sc2

            for tv in ts_out:
                v0, v1, v2 = eval_mp(row[i-1], tv), eval_mp(row[i], tv), eval_mp(row[i+1], tv)
                sc = v1*v1 + abs(v0*v2) + mp.mpf('1e-300')
                r_tur = (v1*v1 - v0*v2) / sc
                if r_tur < min_tur_out: min_tur_out = r_tur

        print(f"k={k:2d} (300 dps): plain_Tur_in>={float(min_tur_in):+.5f} | plain_Tur_out>={float(min_tur_out):+.5f} | D_next_Tur>={float(min_Dtur):+.5f}")
        print(f"              J_raw(-1,-1/m]>={float(min_J_left):+.5f} | Term1>={float(min_t1):+.5f} | Term2(cross)>={float(min_t2):+.5f}")
