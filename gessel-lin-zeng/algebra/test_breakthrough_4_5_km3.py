"""
3-SECOND BREAKTHROUGH TEST FOR ALL-k PROOFS OF OPEN CASES i = 4, i = 5, i = k - 3!
1. For i = 4:
   We ALREADY proved in 1 line: for ANY deg-4 poly q_0 with (1+t) h~_{4,4} << q_0,
     M_k^+ := (prod_{j=5}^k D_j^2)(q_0)  (deg 2k - 4)
   AUTOMATICALLY satisfies  h~_{k, {1,2,3,4}} << M_k^+  FOR ALL k >= 5!
   Test candidates q_0 (with (1+t) h~_{4,4} << q_0):
     - q_0 = D_5 h~_{4,4}  (wait: check if (1+t) h~_{4,4} << D_5 h~_{4,4}!)
     - q_0 = D_4 h~_{4,4}
     - q_0 = D_3 h~_{4,4}
     - q_0 = D_2 h~_{4,4}
     - q_0 = t * h~_{4,4}  (WAIT A SECOND!! (1+t) f << t f IS TRIVIALLY TRUE FOR *EVERY* f IN RZ(-1,0)!!
       CHECK THAT!! Roots of (1+t)f are -1 < r1 < ... < rd; roots of t*f are r1 < ... < rd < 0!
       Wait: they share r1..rd (weak interlacing!), so D_j(t f) vs D_j((1+t) f) STRICTLY interlaces!
       AND D_m f has roots strictly between -1 and 0, with (1+t) f << D_m f FOR EVERY m >= 1!!
       CHECK THAT!! At root r_s of f, (1+t)f = 0 and D_m f(r_s) = r_s(1+r_s) f'(r_s) has OPPOSITE sign to f'(r_s),
       same sign as ((1+t)f)'(r_s)?? No: f << D_m f means r_1 < u_1 < r_2 < u_2 < ... < r_d < u_d < 0,
       so -1 < r_1 < u_1 < r_2 < u_2 < ... < r_d < u_d < 0, wait: we need -1 < u_1 < r_1 < u_2 < r_2 ... (u to the LEFT of r)!
       NO! Look at (1+t) f (deg d+1) << q_0 (deg d+1):
       roots of (1+t)f are -1 < r_1 < r_2 < ... < r_d;
       roots of q_0 are u_1 < u_2 < ... < u_{d+1};
       (1+t) f << q_0 MEANS  -1 < u_1 < r_1 < u_2 < r_2 < ... < u_d < r_d < u_{d+1}!!
       AND FOR q_0 = D_m f: roots of D_m f are r_1 < u_1 < ...? NO!
       Let's check deg(D_m f) = d + 1! Roots of f are r_1 < ... < r_d in (-1, 0).
       At t = -1: D_m f(-1) = (1 - m) f(-1) (OPPOSITE sign to f(-1) for m > 1!!).
       At t = r_1: D_m f(r_1) = r_1(1+r_1) f'(r_1) (OPPOSITE sign to f'(r_1), i.e., SAME sign as f(-1)!!).
       SO D_m f HAS A ROOT u_1 IN (-1, r_1), A ROOT u_2 IN (r_1, r_2), ..., AND A ROOT u_{d+1} IN (r_d, 0)!!
       YES!! FOR *EVERY* m > 1 AND *EVERY* f IN RZ(-1,0):
         (1 + t) f  <<  D_m f   ALWAYS HOLDS!!
       OH!! MY!! GOD!! TEST q_0 = D_m h~_{4,4} (for m = 2, 3, 4, 5, 6) RIGHT NOW!!)
2. Also test i = 5: M_k^{(5)} = (prod_{j=6}^k D_j^2)(D_m h~_{5,5}) for m = 2, 3, 4, 5, 6!
3. Also test i = k - 3: N_k = (prod_{j=4}^k D_{j+1})(q_0) where q_0 << h~_{3,0} (deg 4),
   e.g. q_0 = D_m h~_{2, 0} (deg 4) or q_0 = h~_{3, 1} (deg 4) or q_0 = P_c h~_{3, 0} (deg 4)!
"""
import mpmath as mp
mp.mp.dps = 80

def apply_Dj_mp(p, j):
    out = [mp.mpf(0)] * (len(p) + 1)
    for m, c in enumerate(p):
        if c != 0:
            out[m] += (1 + m) * c; out[m + 1] += (j + m) * c
    return out

def mul_1pt(p):
    out = [mp.mpf(0)] * (len(p) + 1)
    for m, c in enumerate(p):
        out[m] += c; out[m + 1] += c
    return out

def build_h_S(k, S):
    p = [mp.mpf(1)]
    if 1 not in S: p = apply_Dj_mp(p, 2)
    for j in range(2, k + 1):
        p = apply_Dj_mp(p, j + 1)
        if j not in S: p = apply_Dj_mp(p, j + 1)
    return p

def check_interlace_bez(f, g):
    if len(f) == len(g):
        g = mul_1pt(g)
    n = len(g) - 1
    assert len(f) == n, f"deg mismatch {len(f)-1} vs {n}"
    s = (abs(g[0]) / abs(g[-1])) ** (mp.mpf(1) / n)
    gs = [g[m] * (s ** m) for m in range(n + 1)]
    fs = [f[m] * (s ** m) for m in range(n)] + [mp.mpf(0)]
    mxg = max(abs(x) for x in gs); mxf = max(abs(x) for x in fs)
    gs = [x / mxg for x in gs]; fs = [x / mxf for x in fs]
    H = [[mp.mpf(0)] * n for _ in range(n)]
    for r in range(n + 1):
        for col in range(r + 1, n + 1):
            val = -(gs[r]*fs[col] - gs[col]*fs[r])
            for m in range(col - r):
                H[r + m][r + (col - r - 1 - m)] += val
    for i in range(n):
        piv = H[i][i]
        if piv <= 0: return False
        for j in range(i + 1, n):
            fac = H[j][i] / piv
            for k in range(j, n):
                H[j][k] -= fac * H[i][k]
                H[k][j] = H[j][k]
    return True

h0 = {1: [mp.mpf(1), mp.mpf(2)]}
hk_k = {1: [mp.mpf(1)]}
for k in range(2, 42):
    h0[k] = apply_Dj_mp(apply_Dj_mp(h0[k-1], k+1), k+1)
    hk_k[k] = apply_Dj_mp(hk_k[k-1], k+1)

print("=== 1. OPEN CASE i = 4: M_k = (prod_{j=5}^k D_j^2)(q_0) (deg 2k - 4) ===")
# Check base (1+t) h_{4,4} << q_0 AND check (early << M_k, late << M_k) up to k = 35!
q0_cands_4 = [
    ("D_2 h_{4,4}", apply_Dj_mp(hk_k[4], 2)),
    ("D_3 h_{4,4}", apply_Dj_mp(hk_k[4], 3)),
    ("D_4 h_{4,4}", apply_Dj_mp(hk_k[4], 4)),
    ("D_5 h_{4,4}", apply_Dj_mp(hk_k[4], 5)),
    ("D_2 h_{2,0}", apply_Dj_mp(h0[2], 2)),
    ("D_3 h_{2,0}", apply_Dj_mp(h0[2], 3)),
    ("D_4 h_{2,0}", apply_Dj_mp(h0[2], 4)),
    ("D_5 h_{2,0}", apply_Dj_mp(h0[2], 5)),
]
for label, q0 in q0_cands_4:
    base_ok = check_interlace_bez(mul_1pt(hk_k[4]), q0)
    res = []
    M = q0
    for k in range(5, 36):
        M = apply_Dj_mp(apply_Dj_mp(M, k), k)
        if k in [5, 8, 12, 18, 25, 35]:
            early = build_h_S(k, {1, 2, 3, 4})
            late  = build_h_S(k, {k-3, k-2, k-1, k})
            ok_e = check_interlace_bez(early, M)
            ok_l = check_interlace_bez(late, M)
            res.append(f"k={k}:{'OK' if (ok_e and ok_l) else f'e{int(ok_e)}l{int(ok_l)}'}")
    print(f"  {label:12s} (base={str(base_ok):5s}): " + " ".join(res))

print("\n=== 2. OPEN CASE i = 5: M_k = (prod_{j=6}^k D_j^2)(q_0) (deg 2k - 5) ===")
for label, q0 in [
    ("D_2 h_{5,5}", apply_Dj_mp(hk_k[5], 2)),
    ("D_3 h_{5,5}", apply_Dj_mp(hk_k[5], 3)),
    ("D_4 h_{5,5}", apply_Dj_mp(hk_k[5], 4)),
    ("D_6 h_{5,5}", apply_Dj_mp(hk_k[5], 6)),
    ("h_{3,0}", h0[3]),
]:
    base_ok = check_interlace_bez(mul_1pt(hk_k[5]), q0)
    res = []
    M = q0
    for k in range(6, 36):
        M = apply_Dj_mp(apply_Dj_mp(M, k), k)
        if k in [6, 10, 15, 20, 25, 35]:
            early = build_h_S(k, {1, 2, 3, 4, 5})
            late  = build_h_S(k, {k-4, k-3, k-2, k-1, k})
            ok_e = check_interlace_bez(early, M)
            ok_l = check_interlace_bez(late, M)
            res.append(f"k={k}:{'OK' if (ok_e and ok_l) else f'e{int(ok_e)}l{int(ok_l)}'}")
    print(f"  {label:12s} (base={str(base_ok):5s}): " + " ".join(res))

print("\n=== 3. OPEN CASE i = k - 3: N_k = (prod_{j=4}^k D_{j+1})(q_0) (deg k + 1) ===")
# Note: late_{k-3} = (prod_{j=4}^k D_{j+1})(h_{3,0}), deg = 5 + (k - 3) = k + 2!
# So for ANY deg-4 q_0 with q_0 << h_{3,0}, N_k << late_{k-3} holds AUTOMATICALLY for ALL k >= 4!
# Test if N_k << early_{k-3} ALSO holds for ALL k = 5..35!!
for label, q0 in [
    ("D_2 h_{2,0}", apply_Dj_mp(h0[2], 2)),
    ("D_3 h_{2,0}", apply_Dj_mp(h0[2], 3)),
    ("D_4 h_{2,0}", apply_Dj_mp(h0[2], 4)),
    ("D_5 h_{2,0}", apply_Dj_mp(h0[2], 5)),
    ("h_{3,{3}}",   build_h_S(3, {3})),
    ("h_{3,{1}}",   build_h_S(3, {1})),
    ("h_{5,5}",     hk_k[5]),
]:
    base_ok = check_interlace_bez(q0, h0[3])
    res = []
    N = q0
    for k in range(4, 36):
        N = apply_Dj_mp(N, k + 1)
        if k in [4, 8, 12, 18, 25, 35]:
            early = build_h_S(k, set(range(1, k - 2)))
            late  = build_h_S(k, set(range(4, k + 1)))
            ok_e = check_interlace_bez(N, early)
            ok_l = check_interlace_bez(N, late)
            res.append(f"k={k}:{'OK' if (ok_e and ok_l) else f'e{int(ok_e)}l{int(ok_l)}'}")
    print(f"  {label:12s} (base={str(base_ok):5s}): " + " ".join(res))
