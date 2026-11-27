"""
3-SECOND BREAKTHROUGH TEST FOR OPEN CASES i = 4, i = 5, i = k - 3!
Using the EXACT conjugation identity (1+t) D_{j+1} = D_j (1+t):
  (1+t) h~_{k, {1..4}} = (prod_{j=5}^k D_j^2) [ (1+t) h~_{4,4} ]!
  h~_{k-2, 0}          = (prod_{j=5}^k D_{j-1}^2) [ h~_{2,0} ]!
TEST RIGHT NOW for k = 5..30:
1. For i = 4:
   (a) Does h~_{k-2, 0} << h~_{k, S} (equal deg 2k-5!) hold for BOTH late and early?!
       (Wait: earlier D_m h~_{k-2, 0} got better and better as m -> infty, which IS h~_{k-2, 0} << h~_{k, S}!!)
   (b) Does M_k := (prod_{j=5}^k D_j^2)(h~_{2, 0}) (or (prod_{j=5}^k D_j^2)(q_0) for q_0 = h~_{4,4} or D_4 h~_{3,3} or h~_{3,1}?)
       satisfy M_k << h~_{k, early} AND h~_{k, late} << M_k  (WAIT: common interlacer means BOTH on SAME side:
       either M_k << late AND M_k << early, OR late << M_k AND early << M_k!!)!
2. For i = 5:
   Does (prod_{j=6}^k D_j^2)(q_0) give a common interlacer for ALL |S| = 5?!
3. For i = k - 3:
   Remember h~_{k, {1..k-3}} (early) = D_{k+1}^2 D_k^2 D_{k-1}^2 h~_{k-3, k-3},
   while h~_{k, {4..k}} (late) = D_{k+1} D_k ... D_5 (h~_{3, 0})!
   Notice (1+t) h~_{k, {4..k}} = (prod_{j=5}^k D_{j-1}) [(1+t) h~_{3, 0}]!!
   Compare with h~_{k+2, k+2} = D_{k+3} D_{k+2} ... D_5 (h~_{3, 3})!
   Check (prod_{j=5}^k D_{j+1})(q_0) for i = k - 3!!
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
    """Check f << g in (-1, 0) where deg g = deg f + 1 OR deg g = deg f!
    Wait: for deg g = deg f = n, f << g in (-1, 0) <=> f << (1+t) g (deg n+1, since root -1 is left of all roots in (-1,0)!)
    Check! If roots of f are r1<...<rn and roots of g are s1<...<sn in (-1,0):
      roots of (1+t)g are -1 < s1 < ... < sn.
      f << (1+t)g <=> -1 < r1 < s1 < r2 < s2 < ... < rn < sn <=> f << g (equal degree)!
    AND check f << t*g? Roots of t*g are s1<...<sn<0, so f << t*g <=> s1<r1<s2<r2<...<sn<rn<0 <=> g << f!
    GENIUS!! 100% exact equal-degree interlacing via deg+1 Bezoutiant on (f, (1+t)g)!
    """
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

# Build h0[k] and hk_k[k]
h0 = {1: [mp.mpf(1), mp.mpf(2)]}
hk_k = {1: [mp.mpf(1)]}
for k in range(2, 35):
    h0[k] = apply_Dj_mp(apply_Dj_mp(h0[k-1], k+1), k+1)
    hk_k[k] = apply_Dj_mp(hk_k[k-1], k+1)

print("=== 1. OPEN CASE i = 4: CHECK EQUAL-DEGREE INTERLACING h~_{k-2, 0} << h~_{k, S} AND M_k << h~_{k, S} ===")
# Note: M_k^{(shift)} = prod_{j=5}^k D_j^2 (q_0) where q_0 has deg 3!
# Also check base case q_0 << (1+t) h~_{4,4} (which AUTOMATICALLY proves M_k << early for ALL k >= 5!!)!
print(f"  BASE CHECK: h~_{{2,0}} << (1+t) h~_{{4,4}}? {check_interlace_bez(h0[2], mul_1pt(hk_k[4]))}")
print(f"  BASE CHECK: h~_{{4,4}} << (1+t) h~_{{4,4}}? {check_interlace_bez(hk_k[4], mul_1pt(hk_k[4]))} (wait: same roots!)")
print(f"  BASE CHECK: h~_{{3,2}} (deg 3) << (1+t) h~_{{4,4}}? {check_interlace_bez(build_h_S(3, {1,2}), mul_1pt(hk_k[4]))}")
print(f"  BASE CHECK: h~_{{2,0}} << h~_{{4,4}}? {check_interlace_bez(h0[2], hk_k[4])}")
print(f"  BASE CHECK: h~_{{4,4}} << h~_{{2,0}}? {check_interlace_bez(hk_k[4], h0[2])}")

# Test 4 candidates for common lower interlacer M_k (deg 2k - 5) of {h~_{k, S} : |S| = 4}:
#   C1: h~_{k-2, 0} = prod_{j=5}^k D_{j-1}^2 (h~_{2, 0})
#   C2: prod_{j=5}^k D_j^2 (h~_{2, 0})
#   C3: prod_{j=5}^k D_j^2 (h~_{4, 4})
#   C4: h~_{k, late} << h~_{k, early} itself! (Check if early and late STILL interlace at k=30!)
for k in [5, 8, 12, 18, 25, 30]:
    early = build_h_S(k, {1, 2, 3, 4})
    late  = build_h_S(k, {k-3, k-2, k-1, k})
    ok_le = check_interlace_bez(late, early)
    # C1: h_{k-2, 0}
    c1_l = check_interlace_bez(h0[k-2], late)
    c1_e = check_interlace_bez(h0[k-2], early)
    # C2: prod_{j=5}^k D_j^2 (h_{2, 0})
    m2 = h0[2]
    for j in range(5, k + 1): m2 = apply_Dj_mp(apply_Dj_mp(m2, j), j)
    c2_l = check_interlace_bez(m2, late); c2_e = check_interlace_bez(m2, early)
    # C3: prod_{j=5}^k D_j^2 (h_{4, 4})
    m3 = hk_k[4]
    for j in range(5, k + 1): m3 = apply_Dj_mp(apply_Dj_mp(m3, j), j)
    c3_l = check_interlace_bez(m3, late); c3_e = check_interlace_bez(m3, early)
    # Also check UPPER interlacer: late << m2 and early << m2?
    c2_u_l = check_interlace_bez(late, m2); c2_u_e = check_interlace_bez(early, m2)
    print(f"  k={k:2d}: late<<early:{ok_le} | h_{{k-2,0}}<<(l,e):({int(c1_l)},{int(c1_e)}) | D_j^2(h20)<<(l,e):({int(c2_l)},{int(c2_e)}) | (l,e)<<D_j^2(h20):({int(c2_u_l)},{int(c2_u_e)}) | D_j^2(h44)<<(l,e):({int(c3_l)},{int(c3_e)})")

print("\n=== 2. CHECK late << early FOR i = 5 AND i = k - 3 UP TO k = 30! ===")
for k in [6, 10, 15, 20, 25, 30]:
    e5 = build_h_S(k, {1,2,3,4,5}); l5 = build_h_S(k, {k-4,k-3,k-2,k-1,k})
    ek3 = build_h_S(k, set(range(1, k-2))); lk3 = build_h_S(k, set(range(4, k+1)))
    # Also test C_k3: prod_{j=5}^k D_j (h_{3, 0}) vs lk3, ek3!
    mk3 = h0[3]
    for j in range(5, k + 1): mk3 = apply_Dj_mp(mk3, j)
    # deg h0[3] = 5; apply (k-4) times -> deg k+1 = deg(lk3) - 1!
    ok_mk3_l = check_interlace_bez(mk3, lk3)
    ok_mk3_e = check_interlace_bez(mk3, ek3)
    # Also prod_{j=5}^k D_{j+1} (h_{3, 0})? That's literally lk3 itself! (Wait: lk3 has D_5..D_{k+1} on h0[3]!)
    # Check prod_{j=5}^k D_{j+1} (q) where q << h_{3, 0} (deg 4): e.g., q = D_4 h_{2, 0} (deg 4)!!
    mk3b = apply_Dj_mp(h0[2], 4)
    for j in range(5, k + 1): mk3b = apply_Dj_mp(mk3b, j + 1)
    ok_mk3b_l = check_interlace_bez(mk3b, lk3)
    ok_mk3b_e = check_interlace_bez(mk3b, ek3)
    print(f"  k={k:2d}: i=5 late<<early:{check_interlace_bez(l5,e5)} | i=k-3 late<<early:{check_interlace_bez(lk3,ek3)} | D_j(h30)<<(l,e):({int(ok_mk3_l)},{int(ok_mk3_e)}) | D_{{j+1}}(D_4 h20)<<(l,e):({int(ok_mk3b_l)},{int(ok_mk3b_e)})")
