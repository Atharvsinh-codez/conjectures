"""
Fast 4-way probe:
1. q-stability Q_k(t,y) at q=1/2, 1/3, 1/4: t=0 interlacing Q'_k(0,y) << Q_k(0,y) for k=10..40 (150-digit mpmath)
   + full UHP stability for k=15, 20.
2. Pair (3,2): D_{k+2} D_{k+1} htilde_{k-1,i} << htilde_{k-1,i-3} for k=4..16.
3. Family (I) prod (D_{b_j}+y)(1) and Family (II) prod D_{b_j}(D_{b_j}+y)(1) with wild b_j in [1.01, 1e5]:
   check interlacing, Turan on (-1,0), and Quad stability!
4. Combinatorial refinements of JSP_{k,i} at k=4,5: test last_letter (sorted by bar1 < 1 < ... < k),
   pos_of_last_k, and number_of_final_run for ALL-PAIRS interlacing!
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
        res[i] += c
        res[i+1] += j * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c
        res[i+2] += c
    return poly_trim(res)
def poly_eval(p, x):
    acc = mp.mpc(0) if isinstance(x, mp.mpc) else mp.mpf(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc

print("=== 1. q-STABILITY AT t=0 FOR k=10..40 (150-digit mpmath) ===", flush=True)
d_val = {1: {0: 2, 1: 0}}
for k in range(2, 41):
    d_val[k] = {}
    for i in range(k+1):
        c_im1 = math.comb(k-1, i-1) if i >= 1 else 0
        c_i = math.comb(k-1, i) if i <= k-1 else 0
        d_im1 = d_val[k-1][i-1] if i >= 1 else 0
        d_i = d_val[k-1][i] if i <= k-1 else 0
        u_prime = (2*d_i + (k+1)*c_i) + d_im1
        d_val[k][i] = 2*u_prime + (k+1)*math.comb(k, i)

for q_str, q_val in [("1/2", mp.mpf('0.5')), ("0.45", mp.mpf('0.45')), ("1/3", mp.mpf(1)/3), ("1/4", mp.mpf('0.25'))]:
    status = "PASS all"
    min_gap_log = mp.mpf('1e100')
    for k in [10, 15, 20, 25, 30, 35, 40]:
        P0 = [mp.mpf(math.comb(k, i)) * (q_val**(i*(i-1)//2)) for i in range(k+1)]
        P1 = [mp.mpf(d_val[k][i]) * (q_val**(i*(i-1)//2)) for i in range(k+1)]
        rts0 = mp.polyroots(list(reversed(P0)), maxsteps=400, extraprec=50)
        rts1 = mp.polyroots(list(reversed(P1)), maxsteps=400, extraprec=50)
        if max(abs(r.imag) for r in rts0) > mp.mpf('1e-20') or max(abs(r.imag) for r in rts1) > mp.mpf('1e-20'):
            status = f"FAIL complex at k={k}"
            break
        r0 = sorted([r.real for r in rts0])
        r1 = sorted([r.real for r in rts1])
        # All roots negative; expect r1[m] < r0[m] < r1[m+1] (since d_val/C(k,i) decreases with i -> roots of P1 shift LEFT!)
        ok = all(r1[m] < r0[m] for m in range(k)) and all(r0[m] < r1[m+1] for m in range(k-1))
        if not ok:
            status = f"FAIL interlace at k={k}"
            break
        # Check tightest ratio r1[m+1]/r0[m] (wait: all negative, so r0[m] < r1[m+1] < 0 means |r0[m]| > |r1[m+1]|)
        for m in range(k-1):
            g = mp.log(abs(r0[m]) / abs(r1[m+1]))
            if g < min_gap_log: min_gap_log = g
    print(f"  q={q_str:4s}: {status} (min log|r0[m]/r1[m+1]| = {mp.nstr(min_gap_log, 6)})", flush=True)

# Also check UHP stability of Q_k(t, y) at q=1/4 for k=15, 20!
mp.mp.dps = 80
ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 21):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_Dj(ht[k-1][idx], k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_Dj(poly_add(p1, p2), k + 1)

for q_str, q_val in [("1/2", mp.mpf('0.5')), ("1/4", mp.mpf('0.25'))]:
    for k in [12, 16, 20]:
        max_im = mp.mpf('-1e100')
        for tv in [mp.mpc(-10, '1e-6'), mp.mpc(-1, '1e-6'), mp.mpc(-0.99, '1e-6'), mp.mpc(-0.5, '1e-6'),
                   mp.mpc(-0.01, '1e-6'), mp.mpc(0, '1e-6'), mp.mpc(0.01, '1e-6'), mp.mpc(1, '1e-6'), mp.mpc(0, 1)]:
            cy = [poly_eval(ht[k][i], tv) * (q_val**(i*(i-1)//2)) for i in range(k+1)]
            rts = mp.polyroots(list(reversed(cy)), maxsteps=400, extraprec=50)
            for r in rts:
                if r.imag > max_im: max_im = r.imag
        print(f"  Q_{k}(t,y) at q={q_str}: UHP t -> max Im(y) = {mp.nstr(max_im, 6)}", flush=True)

print("\n=== 2. CHECK PAIR (3,2): D_{k+2} D_{k+1} htilde_{k-1,i} << htilde_{k-1,i-3} ===", flush=True)
mp.mp.dps = 50
for k in [4, 6, 8, 10, 12, 16, 20]:
    fails = 0
    for i in range(3, k):
        lhs = apply_Dj(apply_Dj(ht[k-1][i], k+1), k+2)
        rhs = ht[k-1][i-3]
        rL = sorted([r.real for r in mp.polyroots(list(reversed(lhs)), maxsteps=300, extraprec=30)])
        rR = sorted([r.real for r in mp.polyroots(list(reversed(rhs)), maxsteps=300, extraprec=30)])
        # Expect deg(rhs) = deg(lhs)+1, rR[0] < rL[0] < rR[1] < ... < rL[-1] < rR[-1]
        ok = all(rR[m] < rL[m] < rR[m+1] for m in range(len(rL)))
        if not ok:
            fails += 1
    print(f"  k={k:2d}: Pair (3,2) fails = {fails}/{k-3}", flush=True)

print("\n=== 3. FAMILY (I) prod(D_{b_j}+y)(1) AND FAMILY (II) prod D_{b_j}(D_{b_j}+y)(1) WITH WILD b_j ===", flush=True)
random.seed(99)
for fam_name, doubled in [("Family I: prod(D_{b_j}+y)(1)", False), ("Family II: prod D_{b_j}(D_{b_j}+y)(1)", True)]:
    fail_int = 0
    fail_tur = 0
    fail_p21 = 0
    for tr in range(25):
        k = 8
        # Extreme alternating/wild b_j in [1.01, 1e5]!
        if tr == 0: bs = [2, 10000, 2, 10000, 2, 10000, 2, 10000]
        elif tr == 1: bs = [10000, 2, 10000, 2, 10000, 2, 10000, 2]
        elif tr == 2: bs = [10000, 10000, 10000, 10000, 1.01, 1.01, 1.01, 1.01]
        elif tr == 3: bs = [1.01, 1.01, 1.01, 1.01, 10000, 10000, 10000, 10000]
        else: bs = [1.0 + 10**random.uniform(-2, 5) for _ in range(k)]
        polys = {0: [mp.mpf(1)]}
        for j in range(1, k+1):
            b = mp.mpf(bs[j-1])
            new_p = {}
            for i in range(j+1):
                p1 = apply_Dj(polys[i], b) if i in polys else [mp.mpf(0)]
                p2 = polys[i-1] if (i-1) in polys else [mp.mpf(0)]
                u = poly_add(p1, p2)
                new_p[i] = apply_Dj(u, b) if doubled else u
            polys = new_p
        # Check real roots + interlacing
        roots = {}
        ok = True
        for i in range(k+1):
            p = polys[i]
            if len(p) <= 1:
                roots[i] = []
                continue
            # Divide out (1+t) factor if doubled? No, mp.polyroots handles -1 fine.
            rts = mp.polyroots(list(reversed(p)), maxsteps=300, extraprec=30)
            if max(abs(r.imag) for r in rts) > mp.mpf('1e-12'):
                ok = False; break
            roots[i] = sorted([r.real for r in rts])
        if ok:
            # For doubled, all polynomials share a single root at -1 (except k=0); strip -1 before checking strict interlacing
            for i in range(1, k):
                r_lo = [r for r in roots[i] if r > mp.mpf('-0.999999999999')]
                r_hi = [r for r in roots[i-1] if r > mp.mpf('-0.999999999999')]
                if len(r_hi) != len(r_lo) + 1:
                    ok = False; break
                if not all(r_hi[m] < r_lo[m] < r_hi[m+1] for m in range(len(r_lo))):
                    ok = False; break
        if not ok: fail_int += 1
        # Check Turan on (-1, 0)
        tok = True
        for i in range(1, k):
            for tv in [-mp.mpf(m)/30 for m in range(1, 30)] + [-mp.mpf('1e-5'), -mp.mpf('0.99999')]:
                v0 = poly_eval(polys[i], tv)
                vm = poly_eval(polys[i-1], tv)
                vp = poly_eval(polys[i+1], tv)
                if v0*v0 - vm*vp < 0:
                    tok = False; break
            if not tok: break
        if not tok: fail_tur += 1
    print(f"  {fam_name}: Int_fail={fail_int}/25, Turan(-1,0)_fail={fail_tur}/25", flush=True)
