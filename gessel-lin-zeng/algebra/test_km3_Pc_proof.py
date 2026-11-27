"""
2-SECOND EXCITING CHECK FOR AN EXACT ALL-k PROOF OF OPEN CASE i = k - 3:
Remember Ma-Wang Prop 1: 2 A_{k+1, k} = P_{2k+2} A_{k+2, k+2} = S_{2k+1}(P_{2k+2} A_{k+1, k+1})!
And A_{k, k-3} = S_{2k+1}(B_{k, k-3}) where B_{k, k-3} = S_{2k} A_{k-1, k-3} + A_{k-1, k-4}!
BOTH have the EXACT SAME outer operator S_{2k+1}!!
Check for k = 4..25:
1. Does A_{k+1, k} << A_{k, k-3} hold?
2. Does P_{2k+2} A_{k+1, k+1} << S_{2k} A_{k-1, k-3} hold?
3. Does P_{2k+2} A_{k+1, k+1} << A_{k-1, k-4} hold?
4. Does P_{2k} A_{k, k} << A_{k-1, k-3} hold (for individual A_{k-1, S} too?!)?
5. Even better: DOES P_{2k+2} A_{k+1, k+1} << A_{k-1, S} (for ALL |S| = k-4) hold, or ONLY aggregate?!
   Wait: check BOTH individual S (late/early) AND aggregate!
   AND check general r: does A_{k+1, i+2} << A_{k, i} hold?!
"""
import mpmath as mp
mp.mp.dps = 80

def apply_SL_mp(p, L):
    out = [mp.mpf(0)] * (len(p) + 1)
    for m, c in enumerate(p):
        if c != 0:
            out[m] += (1 + m) * c; out[m + 1] += (L - m) * c
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out

def padd_mp(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for idx in range(len(a)): out[idx] += a[idx]
    for idx in range(len(b)): out[idx] += b[idx]
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out

def apply_Pc_mp(p, c):
    return [(c - m) * v for m, v in enumerate(p)]

def check_interlace_bez_x(f, g):
    """Exact 80-digit check of f << g in (-inf, 0) where deg g = deg f + 1, positive lc!
    Wait: for g (deg n) and f (deg n-1) with positive lc, f << g <=> (g(x)f(y)-g(y)f(x))/(x-y) is PD!
    """
    n = len(g) - 1
    assert len(f) == n, f"deg mismatch: len(f)={len(f)}, len(g)={len(g)}"
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
        if piv <= 0: return False, i
        for j in range(i + 1, n):
            fac = H[j][i] / piv
            for k in range(j, n):
                H[j][k] -= fac * H[i][k]
                H[k][j] = H[j][k]
    return True, n

def build_A_S(k, S):
    p = [mp.mpf(1)]
    L = 0
    for j in range(1, k + 1):
        if j in S:
            p = apply_SL_mp(p, L)
            L += 2
        else:
            p = apply_SL_mp(apply_SL_mp(p, L), L + 1)
            L += 3
    return p

# Build A_{k, i} up to k = 22
A = {(1, 0): [mp.mpf(1), mp.mpf(1)], (1, 1): [mp.mpf(1)]}
for k in range(2, 23):
    for i in range(k + 1):
        L = 3*k - 2 - i
        u = [mp.mpf(0)]
        if i <= k - 1: u = padd_mp(u, apply_SL_mp(A[(k-1, i)], L - 1))
        if i >= 1: u = padd_mp(u, A[(k-1, i-1)])
        A[(k, i)] = apply_SL_mp(u, L)

print("=== CHECK OPEN CASE i = k - 3 VIA P_{2k+2} A_{k+1, k+1} ===")
for k in [4, 5, 6, 8, 10, 12, 15, 18, 20]:
    # 1. A_{k+1, k} << A_{k, k-3}
    ok1 = check_interlace_bez_x(A[(k+1, k)], A[(k, k-3)])[0]
    # 2. P_{2k+2} A_{k+1, k+1} << S_{2k} A_{k-1, k-3}
    PA = apply_Pc_mp(A[(k+1, k+1)], 2*k + 2)
    SA = apply_SL_mp(A[(k-1, k-3)], 2*k)
    ok2 = check_interlace_bez_x(PA, SA)[0]
    # 3. P_{2k+2} A_{k+1, k+1} << A_{k-1, k-4}
    ok3 = check_interlace_bez_x(PA, A[(k-1, k-4)])[0]
    # 4. P_{2k} A_{k, k} << A_{k-1, k-3} (AND check individual early/late!)
    ok4_agg = check_interlace_bez_x(apply_Pc_mp(A[(k, k)], 2*k), A[(k-1, k-3)])[0]
    # Check individual late/early for A_{k, S} (|S|=k-3) vs A_{k+1, k}!!
    late = build_A_S(k, set(range(4, k + 1)))
    early = build_A_S(k, set(range(1, k - 2)))
    ok_ind_l = check_interlace_bez_x(A[(k+1, k)], late)[0]
    ok_ind_e = check_interlace_bez_x(A[(k+1, k)], early)[0]
    # Check P_{2k} A_{k, k} << A_{k-1, late} (|S| = k-3)
    late_km1 = build_A_S(k - 1, set(range(3, k)))
    early_km1 = build_A_S(k - 1, set(range(1, k - 2)))
    ok4_l = check_interlace_bez_x(apply_Pc_mp(A[(k, k)], 2*k), late_km1)[0]
    ok4_e = check_interlace_bez_x(apply_Pc_mp(A[(k, k)], 2*k), early_km1)[0]
    print(f"  k={k:2d}: A_{{k+1,k}}<<A_{{k,k-3}}:{ok1} | PA<<SA:{ok2} | PA<<A_{{k-1,k-4}}:{ok3} | P A_k<<A_{{k-1,k-3}}:(agg={ok4_agg},l={ok4_l},e={ok4_e}) | A_{{k+1,k}}<<A_{{k,S}}:(l={ok_ind_l},e={ok_ind_e})")
