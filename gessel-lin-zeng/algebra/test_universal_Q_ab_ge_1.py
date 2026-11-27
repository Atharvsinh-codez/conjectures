"""
CRITICAL TEST OF THE alpha_j >= 1, beta_j >= 1 UNIVERSAL Q-STABILITY HYPOTHESIS!
1. Check our exact endpoint prediction: for [1.01]^4 + [10^4]^4 (beta = 0.01 => beta/(beta+1) = 1/101),
   does q = 1/250 PASS?!
2. Ruthlessly try to BREAK Q_{k, q}(t, y) = prod_{j=1}^k (T_j + w_j Y_q)(1) at q = 1/4, 1/8, 1/16
   with alpha_j >= 1, beta_j >= 1 (i.e. a_j >= 1, b_j >= a_j + 1), arbitrary w_j > 0!
   Test BOTH (a) random wild (a_j, b_j, w_j) spanning [1, 1e6], AND (b) adversarial k=6..10!
3. Also test if a SINGLE-STEP operator F(t,y) -> T F(t,y) + w y F(t, q y) preserves a simple
   quantitative stability class (e.g. real-stable + Hutchinson ratio on R, or sector-stable) for alpha>=1, beta>=1!
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
def apply_T(p, a, b):
    dp = poly_deriv(p)
    res = [0] * (len(p) + 2)
    for i, c in enumerate(p):
        res[i] += a * c; res[i+1] += b * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c; res[i+2] += c
    return poly_trim(res)
def poly_eval(p, x):
    acc = mp.mpc(0) if isinstance(x, mp.mpc) else mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

def check_uhp_Q(polys, k, q_val, t_pts):
    max_im = mp.mpf('-1e100')
    worst_t = None
    for tv in t_pts:
        cy = [poly_eval(polys[i], tv) * (q_val**(i*(i-1)//2)) for i in range(k+1)]
        s = (abs(cy[0]) / abs(cy[-1])) ** (mp.mpf(1) / k) if abs(cy[-1]) > 0 else mp.mpf(1)
        cy_sc = [cy[i] * (s**i) for i in range(k+1)]
        mx = max(abs(c) for c in cy_sc)
        if mx == 0: continue
        cy_sc = [c / mx for c in cy_sc]
        rts = mp.polyroots(list(reversed(cy_sc)), maxsteps=400, extraprec=60)
        for r in rts:
            im = (r * s).imag
            if im > max_im: max_im = im; worst_t = complex(tv)
    return max_im, worst_t

t_uhp = [mp.mpc('-100', '1e-6'), mp.mpc('-1.001', '1e-7'), mp.mpc('-0.999', '1e-7'),
         mp.mpc('-0.5', '1e-6'), mp.mpc('-0.01', '1e-6'), mp.mpc('-1e-6', '1e-10'),
         mp.mpc('0', '1e-7'), mp.mpc('1e-6', '1e-10'), mp.mpc('0.01', '1e-6'),
         mp.mpc('1', '1e-6'), mp.mpc('100', '1e-6'), mp.mpc('0', '1')]

print("=== 1. VERIFY ENDPOINT PREDICTION: [1.01]^4 + [10^4]^4 AT q = 1/250 ===", flush=True)
polys = {0: [mp.mpf(1)]}
bs = [mp.mpf('1.01')]*4 + [mp.mpf(10000)]*4
for step, b in enumerate(bs, start=1):
    polys = {i: poly_add(apply_T(polys[i], mp.mpf(1), b) if i <= step-1 else [mp.mpf(0)],
                         polys[i-1] if i >= 1 else [mp.mpf(0)]) for i in range(step+1)}
for q_str, q_val in [("1/50", mp.mpf(1)/50), ("1/150", mp.mpf(1)/150), ("1/250", mp.mpf(1)/250), ("1/500", mp.mpf(1)/500)]:
    mi, wt = check_uhp_Q(polys, 8, q_val, t_uhp)
    print(f"  [1.01]^4+[1e4]^4 at q={q_str:5s}: max Im(y) = {mp.nstr(mi, 4)} {'FAIL' if mi > mp.mpf('1e-12') else 'OK'} @{wt.real:.2g}", flush=True)

print("\n=== 2. STRESS-TEST UNIVERSAL T_j WITH alpha_j >= 1, beta_j >= 1, WILD w_j > 0 (30 trials, k=6) ===", flush=True)
random.seed(2026)
for q_str, q_val in [("1/2", mp.mpf('0.5')), ("1/4", mp.mpf('0.25')), ("1/8", mp.mpf('0.125')), ("1/16", mp.mpf('0.0625'))]:
    fails = 0
    worst = None
    for tr in range(25):
        k = 6
        polys = {0: [mp.mpf(1)]}
        params = []
        for j in range(1, k+1):
            alpha = mp.mpf(1) + mp.mpf(10 ** random.uniform(-2, 4)) if random.random() < 0.7 else mp.mpf(1)
            beta  = mp.mpf(1) + mp.mpf(10 ** random.uniform(-2, 4)) if random.random() < 0.7 else mp.mpf(1)
            w_j   = mp.mpf(10 ** random.uniform(-3, 3))
            params.append((float(alpha), float(beta), float(w_j)))
            polys = {i: poly_add(apply_T(polys[i], alpha, alpha + beta) if i <= j-1 else [mp.mpf(0)],
                                 [w_j*c for c in polys[i-1]] if i >= 1 else [mp.mpf(0)]) for i in range(j+1)}
        mi, wt = check_uhp_Q(polys, k, q_val, t_uhp)
        if mi > mp.mpf('1e-12'):
            fails += 1
            if worst is None: worst = (tr, float(mi), wt.real, params)
    print(f"  Universal (alpha>=1, beta>=1, wild w_j) at q={q_str:4s}: fails = {fails}/25 {'' if not worst else 'ex: Im='+str(round(worst[1],4))+'@'+str(round(worst[2],3))}", flush=True)

print("\n=== 3. STRESS-TEST FIXED a_j=1, w_j=1, RANDOM b_j >= 2 (30 trials, k=8) ===", flush=True)
random.seed(777)
for q_str, q_val in [("1/2", mp.mpf('0.5')), ("1/4", mp.mpf('0.25')), ("1/8", mp.mpf('0.125'))]:
    fails = 0
    for tr in range(25):
        k = 8
        polys = {0: [mp.mpf(1)]}
        for j in range(1, k+1):
            b = mp.mpf(2) + mp.mpf(10 ** random.uniform(-2, 5)) if random.random() < 0.8 else mp.mpf(2)
            polys = {i: poly_add(apply_T(polys[i], mp.mpf(1), b) if i <= j-1 else [mp.mpf(0)],
                                 polys[i-1] if i >= 1 else [mp.mpf(0)]) for i in range(j+1)}
        mi, wt = check_uhp_Q(polys, k, q_val, t_uhp)
        if mi > mp.mpf('1e-12'): fails += 1
    print(f"  Fixed a=1, w=1, random b_j>=2 (k=8) at q={q_str:4s}: fails = {fails}/25", flush=True)
