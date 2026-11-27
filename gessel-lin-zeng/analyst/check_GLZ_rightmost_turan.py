import mpmath as mp

# Check GLZ plain Turan, D_{k+2} Turan, Term1, Term2, and J_raw AT THE EXACT RIGHTMOST-ROOT SCALE
# t = -alpha / s_{k, i}, where s_{k, i} = h'_{k, i}(0) / h_{k, i}(0), for k = 35, 50, 72!
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

    if k in [35, 50, 72]:
        print(f"\n=== GLZ k = {k} AT RIGHTMOST-ROOT SCALE t = -alpha / s_{{k,i}} ===")
        for i in [2, k//4, k//2, 3*k//4, k-2]:
            h0, h1, h2 = row[i-1], row[i], row[i+1]
            s1 = mp.mpf(h1[1]) / mp.mpf(h1[0])
            # Also D_{k+2} on h0, h1, h2
            Dh0, Dh1, Dh2 = apply_D_int(k+2, h0), apply_D_int(k+2, h1), apply_D_int(k+2, h2)
            # And prev_row for Term1, Term2
            pm2, pm1, p0, p1 = prev_row[i-2], prev_row[i-1], prev_row[i], prev_row[i+1]
            Dpm1, Dp0, Dpp1 = apply_D_int(k+1, pm1), apply_D_int(k+1, p0), apply_D_int(k+1, p1)
            # Sweep alpha across [1e-4, 1e4] logarithmically (30 points) around 1/s1!
            min_tur = (mp.mpf('1e100'), None)
            min_Dtur = (mp.mpf('1e100'), None)
            min_t2 = (mp.mpf('1e100'), None)
            min_J = (mp.mpf('1e100'), None)
            dh0, dh1, dh2 = deriv_int(h0), deriv_int(h1), deriv_int(h2)
            for exp in [-4, -2, -1, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 1.0, 1.3, 1.7, 2.0, 2.5, 3.0, 4.0]:
                alpha = mp.power(10, mp.mpf(str(exp)))
                tv = -alpha / s1
                if tv <= -1: continue
                v0, v1, v2 = eval_mp(h0, tv), eval_mp(h1, tv), eval_mp(h2, tv)
                r_tur = (v1*v1 - v0*v2) / (v1*v1 + abs(v0*v2) + mp.mpf('1e-300'))
                if r_tur < min_tur[0]: min_tur = (r_tur, float(alpha))
                dv0, dv1, dv2 = eval_mp(Dh0, tv), eval_mp(Dh1, tv), eval_mp(Dh2, tv)
                r_Dtur = (dv1*dv1 - dv0*dv2) / (dv1*dv1 + abs(dv0*dv2) + mp.mpf('1e-300'))
                if r_Dtur < min_Dtur[0]: min_Dtur = (r_Dtur, float(alpha))
                # Term 2 of u_i = D_{k+1} p_i + p_{i-1}
                hm2, hm1, hp0 = eval_mp(pm2, tv), eval_mp(pm1, tv), eval_mp(p0, tv)
                dm1, d0, dp1 = eval_mp(Dpm1, tv), eval_mp(Dp0, tv), eval_mp(Dpp1, tv)
                t2 = 2*hm1*d0 - hm2*dp1 - hp0*dm1
                sc2 = abs(2*hm1*d0) + abs(hm2*dp1) + abs(hp0*dm1) + mp.mpf('1e-300')
                if t2/sc2 < min_t2[0]: min_t2 = (t2/sc2, float(alpha))
                # J_raw
                q0, q1, q2 = eval_mp(dh0, tv), eval_mp(dh1, tv), eval_mp(dh2, tv)
                w1 = v1*q0 - q1*v0; w2 = v2*q1 - q2*v1; w12 = v2*q0 - q2*v0
                J = 4*w1*w2/(w12*w12) if w12 != 0 else mp.mpf('1e10')
                if J < min_J[0]: min_J = (J, float(alpha))
            print(f"  i={i:2d} (1/s1={float(1/s1):.2e}): Tur>={float(min_tur[0]):+.5f}@{min_tur[1]:.1e} | D_Tur>={float(min_Dtur[0]):+.5f}@{min_Dtur[1]:.1e} | Term2>={float(min_t2[0]):+.5f}@{min_t2[1]:.1e} | J_raw>={float(min_J[0]):.5f}@{min_J[1]:.1e}")
