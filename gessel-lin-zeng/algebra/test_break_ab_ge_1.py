"""
5-second check:
1. Can we break universal Q-stability (and even real-rootedness!) with alpha_j >= 1, beta_j >= 1, w_j > 0?
   Test if J_i(t) < 0 on (-1,0) already occurs at k=6 with alpha_j >= 1, beta_j >= 1, w_j > 0,
   and if so, push to k=7 Turan break -> k=8 interlacing break -> k=9 complex roots!
2. Also check at k=34 GLZ + step 35 with a_35=1, w_35=1, large b_35: does Turan break at t=-1e-16?
"""
import math, random
import mpmath as mp
mp.mp.dps = 60

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
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

print("=== 1. CAN J_i(t) < 0 ON (-1,0) FOR alpha_j >= 1, beta_j >= 1, w_j > 0 AT k=6? ===", flush=True)
random.seed(123)
tgrid = [-mp.mpf(v) for v in ['0.9999', '0.99', '0.9', '0.7', '0.5', '0.3', '0.1', '0.01', '0.0001']]
found = None
for tr in range(80):
    k = 6
    polys = {0: [mp.mpf(1)]}
    for j in range(1, k+1):
        alpha = mp.mpf(1) + (mp.mpf(10 ** random.uniform(0, 4)) if random.random() < 0.5 else mp.mpf(0))
        beta  = mp.mpf(1) + (mp.mpf(10 ** random.uniform(0, 4)) if random.random() < 0.5 else mp.mpf(0))
        w_j   = mp.mpf(10 ** random.uniform(-4, 4))
        polys = {i: poly_add(apply_T(polys[i], alpha, alpha+beta) if i <= j-1 else [mp.mpf(0)],
                             [w_j*c for c in polys[i-1]] if i >= 1 else [mp.mpf(0)]) for i in range(j+1)}
    dp = {i: poly_deriv(polys[i]) for i in range(k+1)}
    for i in range(1, k):
        for tv in tgrid:
            f0, df0 = poly_eval(polys[i-1], tv), poly_eval(dp[i-1], tv)
            f1, df1 = poly_eval(polys[i], tv),   poly_eval(dp[i], tv)
            f2, df2 = poly_eval(polys[i+1], tv), poly_eval(dp[i+1], tv)
            w01 = f0*df1 - df0*f1; w12 = f1*df2 - df1*f2; w02 = f0*df2 - df0*f2
            if 4*w01*w12 - w02*w02 < 0:
                found = (tr, i, tv, float(4*w01*w12/(w02*w02)), polys)
                break
        if found: break
    if found: break

if found:
    tr, i, tv, ratio, polys = found
    print(f"  FOUND J_i < 0 at trial {tr}, i={i}, t={mp.nstr(tv,6)}: 4W01W12/W02^2 = {ratio:.5f} < 1!", flush=True)
    # Step 7 with alpha7 >= 1, beta7 >= 1!
    dp = {idx: poly_deriv(polys[idx]) for idx in range(7)}
    f0, df0 = poly_eval(polys[i-1], tv), poly_eval(dp[i-1], tv)
    f1, df1 = poly_eval(polys[i], tv),   poly_eval(dp[i], tv)
    f2, df2 = poly_eval(polys[i+1], tv), poly_eval(dp[i+1], tv)
    delta = f1*f1 - f0*f2
    ddelta_u = tv*(1+tv)*(2*f1*df1 - df0*f2 - f0*df2)
    phi_star = -ddelta_u / (2*delta)
    # Solve alpha*(1+tv) + beta*tv = phi_star with alpha >= 1, beta >= 1:
    # (alpha - 1)*(1+tv) + (beta - 1)*tv = phi_star - (1 + 2*tv)
    rhs = phi_star - (1 + 2*tv)
    if rhs >= 0:
        beta7 = mp.mpf(1); alpha7 = mp.mpf(1) + rhs / (1 + tv)
    else:
        alpha7 = mp.mpf(1); beta7 = mp.mpf(1) + rhs / tv
    w7 = mp.mpf('1e-12')
    p7 = {idx: poly_add(apply_T(polys[idx], alpha7, alpha7+beta7) if idx <= 6 else [mp.mpf(0)],
                        [w7*c for c in polys[idx-1]] if idx >= 1 else [mp.mpf(0)]) for idx in range(8)}
    v7 = poly_eval(p7[i], tv)**2 - poly_eval(p7[i-1], tv)*poly_eval(p7[i+1], tv)
    print(f"  Step 7 (alpha7={mp.nstr(alpha7,5)}>=1, beta7={mp.nstr(beta7,5)}>=1): Turan = {mp.nstr(v7,6)} < 0!", flush=True)
    # Step 8: break interlacing with alpha8 >= 1, beta8 >= 1!
    dv_ip1 = poly_eval(poly_deriv(p7[i+1]), tv)
    v_ip1 = poly_eval(p7[i+1], tv)
    v_i = poly_eval(p7[i], tv)
    for sign_phi in [1, -1]:
        broke = False
        for mag in [10, 100, 1000, 1e5, 1e7]:
            phi8 = mp.mpf(sign_phi * mag)
            rhs8 = phi8 - (1 + 2*tv)
            if rhs8 >= 0:
                beta8 = mp.mpf(1); alpha8 = mp.mpf(1) + rhs8/(1+tv)
            else:
                alpha8 = mp.mpf(1); beta8 = mp.mpf(1) + rhs8/tv
            w8 = -(tv*(1+tv)*dv_ip1 + phi8*v_ip1) / v_i
            if w8 <= 0: continue
            q8_ip1 = poly_add(apply_T(p7[i+1], alpha8, alpha8+beta8), [w8*c for c in p7[i]])
            q8_i   = poly_add(apply_T(p7[i], alpha8, alpha8+beta8),   [w8*c for c in p7[i-1]])
            r1 = mp.polyroots(list(reversed(q8_ip1)), maxsteps=200, extraprec=30)
            r0 = mp.polyroots(list(reversed(q8_i)),   maxsteps=200, extraprec=30)
            if max(abs(z.imag) for z in r1) > mp.mpf('1e-10') or max(abs(z.imag) for z in r0) > mp.mpf('1e-10'):
                print(f"  COMPLEX ROOTS ALREADY AT k=8!", flush=True); broke = True; break
            r1 = sorted([z.real for z in r1]); r0 = sorted([z.real for z in r0])
            if not all(r0[m] < r1[m] < r0[m+1] for m in range(len(r1))):
                print(f"  INTERLACING BROKEN AT k=8 (alpha8={mp.nstr(alpha8,4)}>=1, beta8={mp.nstr(beta8,4)}>=1, w8={mp.nstr(w8,4)})!", flush=True)
                # Push to k=9 complex roots!
                for b9 in [mp.mpf(2), mp.mpf(10), mp.mpf(1000), mp.mpf(1e6)]:
                    for a9 in [mp.mpf(1), mp.mpf(100), mp.mpf(1e5)]:
                        Tq = apply_T(q8_ip1, a9, a9+b9)
                        for pw in range(-8, 10):
                            w9 = mp.mpf(10)**pw
                            q9 = poly_add(Tq, [w9*c for c in q8_i])
                            rts9 = mp.polyroots(list(reversed(q9)), maxsteps=200, extraprec=30)
                            mx = max(abs(z.imag) for z in rts9)
                            if mx > mp.mpf('1e-10'):
                                print(f"  BOOM!! COMPLEX ROOTS AT k=9 WITH ALL alpha>=1, beta>=1! (a9={a9}, b9={a9+b9}, w9={w9}, max|Im|={mp.nstr(mx,5)})", flush=True)
                                broke = True; break
                        if broke: break
                    if broke: break
                break
        if broke: break
else:
    print("  No J_i < 0 found in 80 trials at k=6!", flush=True)
