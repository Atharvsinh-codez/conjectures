"""
Fast 10-second check of two key questions:
1. Fixed stress-test of Q_{k, q}(t, y) for b_j >= 2 at q = 1/4, 1/8 (normalize cy by max|cy| + dps=180).
2. For GLZ k = 8, 12, 16: check range [min_t, max_t] of u_{k,i}(t) = W_XY / (2 sqrt(W_XM W_MY)) on R,
   check if u_{i-1}(t) + u_{i+1}(t) <= 2 u_i(t) holds pointwise,
   and check if D_j preserves R = W_AB W_CE / W_BC^2 <= 1 (or <= 4) on random 4-tuples!
"""
import math, random
import mpmath as mp
mp.mp.dps = 150

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
def poly_eval(p, x):
    acc = mp.mpc(0) if isinstance(x, mp.mpc) else mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

print("=== 1. Q_{k,q} ON EXTREME b_j >= 2 (k=8, rescaled polyroots) ===", flush=True)
t_uhp = [mp.mpc('-100', '1e-6'), mp.mpc('-1.001', '1e-7'), mp.mpc('-0.999', '1e-7'),
         mp.mpc('-0.5', '1e-6'), mp.mpc('-1e-6', '1e-10'), mp.mpc('0', '1e-7'),
         mp.mpc('1e-6', '1e-10'), mp.mpc('1', '1e-6'), mp.mpc('100', '1e-6'), mp.mpc('0', '1')]
cases = [
    ("const [2]*10", [2]*10),
    ("step-up [2]*4+[1e6]*4", [2]*4 + [10**6]*4),
    ("step-dn [1e6]*4+[2]*4", [10**6]*4 + [2]*4),
    ("spike [2]*3+[1e6]+[2]*4", [2]*3 + [10**6] + [2]*4),
    ("dip [1e6]*3+[2]+[1e6]*4", [10**6]*3 + [2] + [10**6]*4),
    ("geom 10^j (j=1..7)", [10**j for j in range(1, 8)]),
    ("rev-geom 10^(8-j)", [10**(8-j) for j in range(1, 8)]),
]

for fam_name, doubled in [("Family I (single)", False), ("Family II (doubled)", True)]:
    print(f"  --- {fam_name} ---", flush=True)
    for cname, bs in cases:
        polys = {0: [mp.mpf(1)]}
        for step, b in enumerate(bs, start=1):
            b = mp.mpf(b)
            new_p = {}
            for i in range(step + 1):
                p1 = apply_Dj(polys[i], b) if i <= step - 1 else [mp.mpf(0)]
                p2 = polys[i-1] if i >= 1 else [mp.mpf(0)]
                u = poly_add(p1, p2)
                new_p[i] = apply_Dj(u, b) if doubled else u
            polys = new_p
        k = len(bs)
        out = []
        for q_str, q_val in [("1/2", mp.mpf('0.5')), ("1/4", mp.mpf('0.25')), ("1/8", mp.mpf('0.125'))]:
            max_im = mp.mpf('-1e100')
            err_cnt = 0
            for tv in t_uhp:
                cy = [poly_eval(polys[i], tv) * (q_val**(i*(i-1)//2)) for i in range(k+1)]
                # Rescale y -> s*y so |cy[0]| ~ |cy[k]| (geometric scaling!)
                s = (abs(cy[0]) / abs(cy[-1])) ** (mp.mpf(1) / k) if abs(cy[-1]) > 0 else mp.mpf(1)
                cy_sc = [cy[i] * (s**i) for i in range(k+1)]
                mx = max(abs(c) for c in cy_sc)
                cy_sc = [c / mx for c in cy_sc]
                try:
                    rts = mp.polyroots(list(reversed(cy_sc)), maxsteps=400, extraprec=60)
                    for r in rts:
                        im = (r * s).imag
                        if im > max_im: max_im = im
                except Exception:
                    err_cnt += 1
            tag = f"FAIL({mp.nstr(max_im,3)})" if max_im > mp.mpf('1e-12') else (f"ERR({err_cnt})" if err_cnt else "OK")
            out.append(f"q={q_str}:{tag}")
        print(f"    {cname:25s}: " + "  ".join(out), flush=True)

print("\n=== 2. CHECK RANGE OF u_{k,i}(t) AND CONCAVITY u_{i-1}+u_{i+1} <= 2 u_i FOR GLZ ===", flush=True)
mp.mp.dps = 50
ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 17):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

tgrid = [mp.mpf(v) for v in [-100, -10, -2, -1.1, -1, -0.9, -0.7, -0.5, -0.3, -0.1, -0.01, 0, 0.05, 0.2, 1, 10, 100]]
for k in [8, 12, 16]:
    # Compute u_{k,i}(t) for i = 2..k-1
    D0 = ht[k-1]
    D1 = {idx: apply_Dj(D0[idx], k+1) for idx in range(k)}
    D2 = {idx: apply_Dj(D1[idx], k+1) for idx in range(k)}
    dD0 = {idx: poly_deriv(D0[idx]) for idx in range(k)}
    dD1 = {idx: poly_deriv(D1[idx]) for idx in range(k)}
    dD2 = {idx: poly_deriv(D2[idx]) for idx in range(k)}
    u_tab = {}
    for i in range(2, k):
        u_tab[i] = []
        for tv in tgrid:
            x, dx = poly_eval(D0[i-2], tv), poly_eval(dD0[i-2], tv)
            m, dm = poly_eval(D1[i-1], tv), poly_eval(dD1[i-1], tv)
            y, dy = poly_eval(D2[i], tv),   poly_eval(dD2[i], tv)
            w1 = x*dm - dx*m; w2 = m*dy - dm*y; w3 = x*dy - dx*y
            u_tab[i].append(w3 / (2*mp.sqrt(w1*w2)))
    min_u = min(min(v) for v in u_tab.values())
    max_u = max(max(v) for v in u_tab.values())
    # Check concavity u_{i-1} + u_{i+1} <= 2 u_i pointwise!
    conc_fails = 0
    max_viol = mp.mpf('-1e100')
    for i in range(3, k - 1):
        for idx_t in range(len(tgrid)):
            viol = u_tab[i-1][idx_t] + u_tab[i+1][idx_t] - 2*u_tab[i][idx_t]
            if viol > 0: conc_fails += 1
            if viol > max_viol: max_viol = viol
    print(f"  GLZ k={k:2d}: u in [{mp.nstr(min_u, 5)}, {mp.nstr(max_u, 5)}]; concavity u_{{i-1}}+u_{{i+1}}<=2u_i fails={conc_fails} (max viol={mp.nstr(max_viol, 5)})", flush=True)
