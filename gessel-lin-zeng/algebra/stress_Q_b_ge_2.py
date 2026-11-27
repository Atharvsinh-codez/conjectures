"""
Ruthlessly stress-test Q_{k, q}(t, y) = prod_{j=1}^k (D_{b_j} + Y_q)(1) and doubled prod D_{b_j}(D_{b_j} + Y_q)(1)
for ALL b_j >= 2 at q = 1/2, 1/4, 1/8, 1/16!
Can we BREAK Q_{k, q} for b_j >= 2 (e.g. b_j in [2, 1e8], k=2..12)?
Also test k=3 with extreme b_1, b_2, b_3 >= 2, and test t = -1 + eps, t = 0 + eps, t = 1e-10 + 1e-14j!
"""
import math, random
import mpmath as mp
mp.mp.dps = 120

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

t_uhp = [
    mp.mpc('-100', '1e-6'), mp.mpc('-2', '1e-6'), mp.mpc('-1.0001', '1e-8'),
    mp.mpc('-1', '1e-8'), mp.mpc('-0.9999', '1e-8'), mp.mpc('-0.5', '1e-6'),
    mp.mpc('-0.01', '1e-6'), mp.mpc('-1e-6', '1e-10'), mp.mpc('-1e-10', '1e-14'),
    mp.mpc('0', '1e-8'), mp.mpc('1e-10', '1e-14'), mp.mpc('1e-6', '1e-10'),
    mp.mpc('0.01', '1e-6'), mp.mpc('1', '1e-6'), mp.mpc('100', '1e-6'), mp.mpc('0', '1')
]

cases = [
    ("const [2]*12", [2]*12),
    ("const [1000]*12", [1000]*12),
    ("step-up [2]*6+[1e8]*6", [2]*6 + [10**8]*6),
    ("step-dn [1e8]*6+[2]*6", [10**8]*6 + [2]*6),
    ("spike [2]*4+[1e8]+[2]*4", [2]*4 + [10**8] + [2]*4),
    ("dip [1e8]*4+[2]+[1e8]*4", [10**8]*4 + [2] + [10**8]*4),
    ("geom 2^j (j=1..10)", [2**j for j in range(1, 11)]),
    ("rev-geom 2^(11-j)", [2**(11-j) for j in range(1, 11)]),
    ("super-geom 10^j (j=1..8)", [10**j for j in range(1, 9)]),
    ("rev-super 10^(9-j)", [10**(9-j) for j in range(1, 9)]),
]

for fam_name, doubled in [("Family I: prod(D_{b_j}+Y_q)(1)", False), ("Family II: prod D_{b_j}(D_{b_j}+Y_q)(1)", True)]:
    print(f"=== {fam_name} ===", flush=True)
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
            worst_t = None
            for tv in t_uhp:
                cy = [poly_eval(polys[i], tv) * (q_val**(i*(i-1)//2)) for i in range(k+1)]
                rts = mp.polyroots(list(reversed(cy)), maxsteps=400, extraprec=60)
                for r in rts:
                    if r.imag > max_im:
                        max_im = r.imag; worst_t = complex(tv)
            if max_im > mp.mpf('1e-12'):
                out.append(f"q={q_str}:FAIL({mp.nstr(max_im,3)}@{worst_t.real:.1e})")
            else:
                out.append(f"q={q_str}:OK")
        print(f"  {cname:26s}: " + "  ".join(out), flush=True)
