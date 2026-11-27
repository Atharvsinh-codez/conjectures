"""
CRITICAL ADVERSARIAL TESTS:
(A) Can we break Universal Turan AND Universal Interlacing / Real-Rootedness of prod_{j=1}^k (T_j + w_j y)(1)
    by picking adversarial (a_j, b_j, w_j) where J_i(t) < 0?
(B) Check GLZ Turan htilde_{k,i}(t)^2 >= htilde_{k,i-1}(t) htilde_{k,i+1}(t) on (-1, 0) —
    ESPECIALLY near t = -1/(k+1) (where Phi = 1 + (k+1)t = 0!) — for k = 35, 40, 45, 50!
(C) Check GLZ J_i(t) on (-1, 0) at t = -1/(k+1) for k = 35..50, AND check E_{k,i} = D_{k+1}^i htilde_{k-1,i}
    J_i^E(t) = 4 W(E_{i-1},E_i) W(E_i,E_{i+1}) - W(E_{i-1},E_{i+1})^2 on ALL of R for k=10..35!
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
        res[i] += a * c
        res[i+1] += b * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c
        res[i+2] += c
    return poly_trim(res)
def poly_eval(p, x):
    acc = mp.mpf(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc

print("=== (A) ADVERSARIAL ATTACK ON UNIVERSAL prod_{j=1}^k (T_j + w_j y)(1) ===", flush=True)
# Replay trial 41 from test_wronskian_turan_invariant.py (or find any k=6 with J_i(t) < 0)
random.seed(7)
found_bad = None
for tr in range(50):
    k = 6
    polys = {0: [mp.mpf(1)]}
    params = []
    for j in range(1, k+1):
        a_j = mp.mpf(10 ** random.uniform(-2, 2))
        b_j = a_j + mp.mpf(10 ** random.uniform(-2, 2))
        w_j = mp.mpf(10 ** random.uniform(-2, 2))
        params.append((a_j, b_j, w_j))
        new_p = {}
        for i in range(j+1):
            p1 = apply_T(polys[i], a_j, b_j) if i in polys else [mp.mpf(0)]
            p2 = [w_j * c for c in polys[i-1]] if (i-1) in polys else [mp.mpf(0)]
            new_p[i] = poly_add(p1, p2)
        polys = new_p
    dp = {i: poly_deriv(polys[i]) for i in range(k+1)}
    for i in range(1, k):
        for tv in [-mp.mpf('0.9999'), -mp.mpf('0.9'), -mp.mpf('0.5'), -mp.mpf('0.1'), -mp.mpf('0.001')]:
            f0, df0 = poly_eval(polys[i-1], tv), poly_eval(dp[i-1], tv)
            f1, df1 = poly_eval(polys[i], tv),   poly_eval(dp[i], tv)
            f2, df2 = poly_eval(polys[i+1], tv), poly_eval(dp[i+1], tv)
            w01 = f0*df1 - df0*f1
            w12 = f1*df2 - df1*f2
            w02 = f0*df2 - df0*f2
            if 4*w01*w12 - w02*w02 < 0:
                # Compute optimal Phi* = -Delta'(u) / (2 Delta(u))
                # Delta = f1^2 - f0*f2, Delta'_u = t(1+t)*(2*f1*df1 - df0*f2 - f0*df2)
                delta = f1*f1 - f0*f2
                ddelta_u = tv*(1+tv)*(2*f1*df1 - df0*f2 - f0*df2)
                phi_star = -ddelta_u / (2*delta)
                # Solve a + b*tv = phi_star with 0 < a < b:
                # write b = a + beta (a > 0, beta > 0): a*(1+tv) + beta*tv = phi_star
                if phi_star >= 0:
                    beta = mp.mpf('0.01')
                    a7 = (phi_star - beta*tv) / (1 + tv)
                    b7 = a7 + beta
                else:
                    a7 = mp.mpf('0.01')
                    beta = (phi_star - a7*(1+tv)) / tv
                    b7 = a7 + beta
                # Check (T7 p_i)^2 - (T7 p_{i-1})(T7 p_{i+1}) at tv!
                q0 = apply_T(polys[i-1], a7, b7)
                q1 = apply_T(polys[i], a7, b7)
                q2 = apply_T(polys[i+1], a7, b7)
                v_turan = poly_eval(q1, tv)**2 - poly_eval(q0, tv)*poly_eval(q2, tv)
                if v_turan < 0:
                    found_bad = (tr, i, tv, a7, b7, v_turan, polys)
                    break
        if found_bad: break
    if found_bad: break

if found_bad:
    tr, i, tv, a7, b7, v_turan, polys = found_bad
    print(f"  BOOM! UNIVERSAL TURAN BROKEN AT k=7! (from trial {tr}, i={i}, t={mp.nstr(tv,8)})", flush=True)
    print(f"  a7={mp.nstr(a7,6)}, b7={mp.nstr(b7,6)} (0<a7<b7: {0<a7<b7})", flush=True)
    print(f"  At w7 -> 0+: q1(t)^2 - q0(t)*q2(t) = {mp.nstr(v_turan, 8)} < 0!", flush=True)
    # Apply step 7 with w7 = 1e-10 and verify Turan fails strictly at k=7!
    w7 = mp.mpf('1e-10')
    p7 = {}
    for idx in range(8):
        p_a = apply_T(polys[idx], a7, b7) if idx in polys else [mp.mpf(0)]
        p_b = [w7 * c for c in polys[idx-1]] if (idx-1) in polys else [mp.mpf(0)]
        p7[idx] = poly_add(p_a, p_b)
    v7 = poly_eval(p7[i], tv)**2 - poly_eval(p7[i-1], tv)*poly_eval(p7[i+1], tv)
    print(f"  With w7=1e-10: q_{{7,{i}}}(t)^2 - q_{{7,{i-1}}}(t)*q_{{7,{i+1}}}(t) = {mp.nstr(v7, 8)}", flush=True)
    # NOW TRY TO BREAK INTERLACING AT STEP 8!
    # Look for t* where q_{7,i}(t)^2 < q_{7,i-1}(t)*q_{7,i+1}(t) AND pick a8, b8, w8 so q_{8,i} does NOT interlace q_{8,i-1}!
    # Remember from Theorem A proof: at a root r of q_{8,i} = T8 p7[i] + w8 p7[i-1],
    # we want p7[i-1](r) * (T8 p7[i-1](r) + w8 p7[i-2](r)) to have WRONG sign!
    # Wait: pick t* = tv where p7[i]^2 < p7[i-1]*p7[i+1] (so p7[i-1], p7[i+1] have SAME sign!)
    # Let's scan w8 and (a8, b8) to make T8 p7[i+1](tv) + w8 p7[i](tv) = 0 (so tv is a root of q_{8, i+1}!),
    # then evaluate q_{8, i}(tv) = T8 p7[i](tv) + w8 p7[i-1](tv)!!
    # Specifically, if we set w8 = - T8 p7[i+1](tv) / p7[i](tv) (need w8 > 0!),
    # then at r = tv, p7[i](r) * q_{8, i}(r) = p7[i](r)*T8 p7[i](r) - p7[i-1](r)*T8 p7[i+1](r)
    # = t(1+t)*(p7[i]*p7[i]' - p7[i-1]*p7[i+1]') + Phi8*(p7[i]^2 - p7[i-1]*p7[i+1])!!
    # SINCE p7[i]^2 - p7[i-1]*p7[i+1] < 0, BY MAKING Phi8 VERY LARGE POSITIVE (or NEGATIVE), WE CAN FLIP THE SIGN OF p7[i](r)*q_{8,i}(r)!!
    v_im1 = poly_eval(p7[i-1], tv)
    v_i   = poly_eval(p7[i], tv)
    v_ip1 = poly_eval(p7[i+1], tv)
    dv_im1 = poly_eval(poly_deriv(p7[i-1]), tv)
    dv_i   = poly_eval(poly_deriv(p7[i]), tv)
    dv_ip1 = poly_eval(poly_deriv(p7[i+1]), tv)
    print(f"  At tv={mp.nstr(tv,6)}: p7[i-1]={mp.nstr(v_im1,4)}, p7[i]={mp.nstr(v_i,4)}, p7[i+1]={mp.nstr(v_ip1,4)}", flush=True)
    # Try sweeping Phi8 in [-1e6, 1e6] so that w8 = -(T8 p7[i+1](tv))/p7[i](tv) > 0, and check roots of q_{8, i+1} vs q_{8, i}!
    broke_int = False
    broke_rr = False
    for sign_phi in [1, -1]:
        for mag in [1, 10, 100, 1000, 10000, 1e6]:
            phi8 = mp.mpf(sign_phi * mag)
            if phi8 >= 0:
                beta8 = mp.mpf('0.01')
                a8 = (phi8 - beta8*tv) / (1 + tv)
                b8 = a8 + beta8
            else:
                a8 = mp.mpf('0.01')
                beta8 = (phi8 - a8*(1+tv)) / tv
                b8 = a8 + beta8
            T8_ip1_val = tv*(1+tv)*dv_ip1 + phi8*v_ip1
            w8 = -T8_ip1_val / v_i
            if w8 <= 0: continue
            # Build q_{8, i+1} and q_{8, i} and q_{8, i-1}!
            q8_ip1 = poly_add(apply_T(p7[i+1], a8, b8), [w8*c for c in p7[i]])
            q8_i   = poly_add(apply_T(p7[i], a8, b8),   [w8*c for c in p7[i-1]])
            rts_ip1 = mp.polyroots(list(reversed(q8_ip1)), maxsteps=200, extraprec=30)
            rts_i   = mp.polyroots(list(reversed(q8_i)),   maxsteps=200, extraprec=30)
            im_ip1 = max(abs(r.imag) for r in rts_ip1)
            im_i   = max(abs(r.imag) for r in rts_i)
            if im_ip1 > mp.mpf('1e-10') or im_i > mp.mpf('1e-10'):
                print(f"  WHOA! NONREAL ROOT AT k=8! phi8={phi8}, im_ip1={im_ip1}, im_i={im_i}", flush=True)
                broke_rr = True
            else:
                r1 = sorted([r.real for r in rts_ip1])
                r0 = sorted([r.real for r in rts_i])
                # Check r0[0] < r1[0] < r0[1] < ... < r1[-1] < r0[-1]
                ok = True
                for m in range(len(r1)):
                    if not (r0[m] < r1[m] < r0[m+1]):
                        ok = False; break
                if not ok:
                    print(f"  BOOM! UNIVERSAL INTERLACING BROKEN AT k=8! phi8={float(phi8)}, a8={mp.nstr(a8,5)}, b8={mp.nstr(b8,5)}, w8={mp.nstr(w8,5)}", flush=True)
                    print(f"    r(q8_{i+1}) = {[float(x) for x in r1]}", flush=True)
                    print(f"    r(q8_{i})   = {[float(x) for x in r0]}", flush=True)
                    broke_int = True
                    # NOW CHECK STEP 9 REAL-ROOTEDNESS! Since q8_{i+1} does NOT interlace q8_i,
                    # T9 q8_{i+1} + w9 q8_i CAN HAVE COMPLEX ROOTS?! Wait: q8_{i+1} << T9 q8_{i+1}, if q8_{i+1} and q8_i don't interlace, try varying w9!
                    for w9 in [mp.mpf(10)**p for p in range(-6, 7)]:
                        q9 = poly_add(apply_T(q8_ip1, mp.mpf(1), mp.mpf(2)), [w9*c for c in q8_i])
                        rts9 = mp.polyroots(list(reversed(q9)), maxsteps=200, extraprec=30)
                        im9 = max(abs(r.imag) for r in rts9)
                        if im9 > mp.mpf('1e-10'):
                            print(f"  BOOM!! UNIVERSAL REAL-ROOTEDNESS BROKEN AT k=9! w9={w9}, max|Im|={mp.nstr(im9,6)}", flush=True)
                            broke_rr = True
                            break
                    break
        if broke_int: break
    if not broke_int:
        print("  Interlacing survived all tested phi8 at k=8!", flush=True)

print("\n=== (B) CHECK GLZ TURAN ON (-1, 0) FOR k = 35, 40, 45, 50 (ESP. NEAR t = -1/(k+1)!) ===", flush=True)
# Compute ht up to 48
ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 49):
    ht[k] = {}
    for i in range(k + 1):
        p1 = apply_T(ht[k-1][i], 1, k + 1) if i <= k - 1 else [0]
        p2 = ht[k-1][i-1] if i >= 1 else [0]
        ht[k][i] = apply_T(poly_add(p1, p2), 1, k + 1)

for k in [35, 40, 44, 48]:
    tgrid = [-mp.mpf(1)/(k+1), -mp.mpf(1)/(2*k), -mp.mpf(2)/(k+1), -mp.mpf(1)/(k*k),
             -mp.mpf('1e-6'), -mp.mpf('1e-12'), -mp.mpf('0.01'), -mp.mpf('0.05'), -mp.mpf('0.1'),
             -mp.mpf('0.25'), -mp.mpf('0.5'), -mp.mpf('0.75'), -mp.mpf('0.9'), -mp.mpf('0.99')]
    min_r = mp.mpf('1e100')
    neg_cnt = 0
    worst = None
    for i in range(1, k):
        for tv in tgrid:
            v0 = poly_eval(ht[k][i], tv)
            vm = poly_eval(ht[k][i-1], tv)
            vp = poly_eval(ht[k][i+1], tv)
            if v0*v0 - vm*vp < 0:
                neg_cnt += 1
            if vm*vp > 0:
                r = (v0*v0)/(vm*vp)
                if r < min_r:
                    min_r = r
                    worst = (i, float(tv))
    print(f"  GLZ k={k:2d}: Turan(-1,0) neg={neg_cnt}, min h_i^2/(h- h+)={mp.nstr(min_r, 8)} at {worst}", flush=True)
