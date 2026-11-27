"""
Fast 15-second algebraic tests:
1. Does Q_{k, q}(t, y) for Family I prod(D_{b_j} + Y_q)(1) and Family II prod D_{b_j}(D_{b_j} + Y_q)(1)
   stay real-stable for EXTREME b_j ([2]*5+[1e4]*5, [1e4]*5+[2]*5, [2,1e4]*5, random b_j) at q=1/2, 1/4, 1/8?!
2. Does Family I or Family II EVER break interlacing or Turan(-1,0) on extreme b_j (k=10)?
3. Check P_{L+1} A_{k, i-1} vs A_{k, i} (from exact intertwining P_{L+1} S_{L+1} = S_L P_{L+1})!
4. Check if Q_{k, q}(t, y) at q=1/4 is STILL real-stable for GLZ at k=35 (near t=0, t=1e-16, t=-0.5, t=10)!
"""
import math, random
import mpmath as mp
mp.mp.dps = 100

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
    acc = mp.mpc(0) if isinstance(x, mp.mpc) else mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

print("=== 1. TEST Q_{k,q} ON EXTREME b_j (a_j=1, w_j=1 FIXED!) ===", flush=True)
t_uhp = [mp.mpc(-10, '1e-6'), mp.mpc(-1, '1e-6'), mp.mpc(-0.99, '1e-6'), mp.mpc(-0.5, '1e-6'),
         mp.mpc(-0.01, '1e-6'), mp.mpc(-1e-5, '1e-8'), mp.mpc(0, '1e-6'), mp.mpc(1e-5, '1e-8'),
         mp.mpc(0.01, '1e-6'), mp.mpc(1, '1e-6'), mp.mpc(0, 1)]

cases = [
    ("step-up [2]*4+[10000]*4", [2]*4 + [10000]*4),
    ("step-dn [10000]*4+[2]*4", [10000]*4 + [2]*4),
    ("alt [2,10000]*4", [2, 10000]*4),
    ("alt [10000,2]*4", [10000, 2]*4),
    ("tiny-hi [1.01]*4+[10000]*4", [mp.mpf('1.01')]*4 + [10000]*4),
    ("hi-tiny [10000]*4+[1.01]*4", [10000]*4 + [mp.mpf('1.01')]*4),
]

for fam_name, doubled in [("Family I (single D_b+Y_q)", False), ("Family II (doubled D_b(D_b+Y_q))", True)]:
    print(f"  --- {fam_name} ---", flush=True)
    for cname, bs in cases:
        polys = {0: [mp.mpf(1)]}
        for step, b in enumerate(bs, start=1):
            new_p = {}
            for i in range(step + 1):
                p1 = apply_Dj(polys[i], b) if i <= step - 1 else [mp.mpf(0)]
                p2 = polys[i-1] if i >= 1 else [mp.mpf(0)]
                u = poly_add(p1, p2)
                new_p[i] = apply_Dj(u, b) if doubled else u
            polys = new_p
        k = len(bs)
        res_str = []
        for q_str, q_val in [("1/2", mp.mpf('0.5')), ("1/4", mp.mpf('0.25')), ("1/8", mp.mpf('0.125'))]:
            max_im = mp.mpf('-1e100')
            for tv in t_uhp:
                cy = [poly_eval(polys[i], tv) * (q_val**(i*(i-1)//2)) for i in range(k+1)]
                rts = mp.polyroots(list(reversed(cy)), maxsteps=300, extraprec=50)
                for r in rts:
                    if r.imag > max_im: max_im = r.imag
            res_str.append(f"q={q_str}:{'FAIL('+mp.nstr(max_im,3)+')' if max_im > mp.mpf('1e-12') else 'OK'}")
        print(f"    {cname:28s}: " + "  ".join(res_str), flush=True)

print("\n=== 2. CHECK P_{L+1} A_{k, i-1} vs A_{k, i} (WHERE L = 3k-2-i) ===", flush=True)
A = {0: {0: [1]}}
for k in range(1, 11):
    A[k] = {}
    for i in range(k+1):
        p1 = apply_SL(A[k-1][i], 3*k-3-i) if i <= k-1 else [0]
        p2 = A[k-1][i-1] if i >= 1 else [0]
        A[k][i] = apply_SL(poly_add(p1, p2), 3*k-2-i)

for k in [4, 6, 8, 10]:
    ok_all = 0; tot = 0
    for i in range(1, k+1):
        L = 3*k - 2 - i
        # P_{L+1} f = (L+1) f - x f'
        f = A[k][i-1]
        Pf = [((L+1) - m) * c for m, c in enumerate(f)]
        g = A[k][i]
        rP = sorted([r.real for r in mp.polyroots(list(reversed(Pf)), maxsteps=300, extraprec=30)])
        rg = sorted([r.real for r in mp.polyroots(list(reversed(g)), maxsteps=300, extraprec=30)]) if len(g)>1 else []
        # Check g << P_{L+1} f (deg P = deg g + 1: rP[0] < rg[0] < rP[1] < ... < rP[-1])
        ok = all(rP[m] < rg[m] < rP[m+1] for m in range(len(rg)))
        if ok: ok_all += 1
        tot += 1
    print(f"  k={k:2d}: A_{{k,i}} << P_{{L+1}} A_{{k,i-1}}: {ok_all}/{tot}", flush=True)

print("\n=== 3. CHECK Q_{k, q} FOR GLZ AT k=24, 30, 35 WITH q=1/2 AND q=1/4 ===", flush=True)
mp.mp.dps = 180
ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 36):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

# Test at t = 1e-16 + 1e-20j (where A_i^2 - A_{i-1}A_{i+1} broke!), t = -0.5 + 1e-6j, t = 0 + 1e-6j, t = 1 + 1e-6j!
t_check = [mp.mpc('1e-16', '1e-20'), mp.mpc('1e-6', '1e-10'), mp.mpc('-0.5', '1e-6'), mp.mpc('-0.99', '1e-6'), mp.mpc('1', '1e-6')]
for q_str, q_val in [("1/2", mp.mpf('0.5')), ("1/4", mp.mpf('0.25'))]:
    for k in [24, 30, 35]:
        max_im = mp.mpf('-1e100')
        for tv in t_check:
            cy = [poly_eval(ht[k][i], tv) * (q_val**(i*(i-1)//2)) for i in range(k+1)]
            try:
                rts = mp.polyroots(list(reversed(cy)), maxsteps=500, extraprec=80)
                for r in rts:
                    if r.imag > max_im: max_im = r.imag
            except Exception as e:
                print(f"    Q_{k} q={q_str} at t={complex(tv)}: err {e}", flush=True)
        print(f"  GLZ Q_{k}(t,y) q={q_str}: max Im(y) = {mp.nstr(max_im, 5)}", flush=True)
