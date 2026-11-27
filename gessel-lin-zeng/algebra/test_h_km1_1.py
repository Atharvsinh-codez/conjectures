"""
2-second check:
Does M_k = h~_{k-1, {1}} = D_k^2 D_{k-1}^2 ... D_3^2 (1) (deg 2k - 4)
interlace BOTH early AND late for i = 4 (ALL k = 5..35)?!
AND check h~_{k-1, {k-1}} = D_k D_{k-1}^2 ... D_3^2 D_2(1) (= H_{k-1}? Wait: H_{k-1} was D_k h_{k-2,0} = D_k D_{k-1}^2 ... D_2(1)!)
Wait: compare h~_{k-1, {1}} vs H_{k-1} = h~_{k-1, {k-1}}!
Because root monotonicity moves {k-1} -> {1} strictly RIGHT, h~_{k-1, {1}} IS SHIFTED RIGHT OF H_{k-1}!!
AND H_{k-1} ONLY FAILED `early` AT k=10 BECAUSE H_{k-1} WAS TOO FAR LEFT!!
WHAT ABOUT h~_{k-1, {1}} (WHICH IS RIGHTMOST AMONG |S'|=1)?!
TEST h~_{k-1, {1}} FOR i = 4, AND h~_{k-1, {1, 2}} FOR i = 5, AND h~_{k+1, {k+1}} / h~_{k+1, {1}} FOR i = k - 3!!
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
    if len(f) == len(g): g = mul_1pt(g)
    n = len(g) - 1
    assert len(f) == n
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

print("=== 1. OPEN CASE i = 4: TEST M_k = h~_{k-1, {1}} = D_k^2 ... D_3^2 (1) ===")
# Check base case (1+t) h~_{4,4} << D_4^2 D_3^2(1)!!
q0 = build_h_S(3, {1})  # D_4^2 D_3^2 (1), deg 4!
h44 = build_h_S(4, {1, 2, 3, 4})  # D_5 D_4 D_3 (1), deg 3!
print(f"  BASE CHECK (1+t) h~_{{4,4}} << h~_{{3,{1}}}? {check_interlace_bez(mul_1pt(h44), q0)}")
for k in [5, 8, 12, 18, 25, 35]:
    early = build_h_S(k, {1, 2, 3, 4})
    late  = build_h_S(k, {k-3, k-2, k-1, k})
    M = build_h_S(k - 1, {1})
    print(f"  k={k:2d}: early << h_{{k-1,{{1}}}}: {check_interlace_bez(early, M)} | late << h_{{k-1,{{1}}}}: {check_interlace_bez(late, M)}")

print("\n=== 2. OPEN CASE i = 5: TEST M_k = h~_{k-1, {1, 2}} = D_k^2 ... D_4^2 D_3 (1) ===")
for k in [6, 10, 15, 20, 25, 35]:
    early = build_h_S(k, {1, 2, 3, 4, 5})
    late  = build_h_S(k, {k-4, k-3, k-2, k-1, k})
    M = build_h_S(k - 1, {1, 2})
    print(f"  k={k:2d}: early << h_{{k-1,{{1,2}}}}: {check_interlace_bez(early, M)} | late << h_{{k-1,{{1,2}}}}: {check_interlace_bez(late, M)}")

print("\n=== 3. OPEN CASE i = k - 3: TEST N_k = h~_{k+1, {2..k+1}} (1 double at 1!) vs {1..k-3} and {4..k} ===")
for k in [5, 8, 12, 18, 25, 35]:
    early = build_h_S(k, set(range(1, k - 2)))
    late  = build_h_S(k, set(range(4, k + 1)))
    # N_k has deg k+1 = 2(k+1) - |S'| - 1 => |S'| = k (1 double step!)
    # Late-most |S'|=k at k+1 is S' = {2, 3, ..., k+1} (1 double step at j=1: D_{k+2}...D_3 D_2(1)!)
    N_late = build_h_S(k + 1, set(range(2, k + 2)))
    N_early = build_h_S(k + 1, set(range(1, k + 1)))  # 1 double step at j=k+1
    print(f"  k={k:2d}: h_{{k+1,{{2..k+1}}}} << (l,e): ({int(check_interlace_bez(N_late, late))},{int(check_interlace_bez(N_late, early))}) | h_{{k+1,{{1..k}}}} << (l,e): ({int(check_interlace_bez(N_early, late))},{int(check_interlace_bez(N_early, early))})")
