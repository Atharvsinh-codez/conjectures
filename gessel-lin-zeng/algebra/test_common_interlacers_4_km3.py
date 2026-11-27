"""
3-SECOND TEST OF EXACT COMMON INTERLACERS FOR OPEN CASES i = 4 AND i = k - 3!!
1. For i = 4 (where mu(early) - mu(late) -> 4/5 < 1, so common interlacer DOES EXIST for all k!!):
   Test M_k = D_m h~_{k-2, 0} (which is S_{2k - 6 + m} A_{k-2, 0}!) for m = k+1, k+2, 1.5*k, 2*k, 3*k, 5*k!
   (Wait: as m increases, D_m h~_{k-2, 0} shifts RIGHT! At m = k it was too far left for `early` at k=10!)
   Also check if t * h~_{k-2, 0} (m -> infty!) interlaces `early` AND `late`!
   Also check Ma-Wang's exact candidates: H_{k} = D_{k+1} h~_{k-1, 0} has deg 2k-2 (too big),
   we need deg 2k - 4 = deg(h~_{k, 4}) + 1 = deg(h~_{k-2, 0}) + 1!
   Other deg 2k - 4 candidates: P_c h~_{k-1, 0}? No, deg(h~_{k-1, 0}) = 2k - 3, so h'_{k-1, 0} and P_{2k-3} h~_{k-1, 0} have deg 2k - 4!!
   WHOA!! TEST BOTH D_m h~_{k-2, 0} AND A'_{k-1, 0} AND P_{2k-3} A_{k-1, 0} FOR i = 4!!
2. For i = k - 3:
   Test M_k = P_c A_{k+2, k+2} for c = 1.2k, 1.4k, 1.5k, 1.6k, 1.8k, 2k, 2k+2!
   AND test whether `early` and `late` for i = k - 3 EVEN HAVE a common interlacer at k = 20, 30!
"""
import mpmath as mp
mp.mp.dps = 80

def apply_Dj_mp(p, j):
    out = [mp.mpf(0)] * (len(p) + 1)
    for m, c in enumerate(p):
        if c != 0:
            out[m] += (1 + m) * c; out[m + 1] += (j + m) * c
    return out

def apply_SL_mp(p, L):
    out = [mp.mpf(0)] * (len(p) + 1)
    for m, c in enumerate(p):
        if c != 0:
            out[m] += (1 + m) * c; out[m + 1] += (L - m) * c
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out

def apply_Pc_mp(p, c):
    out = [(c - m) * v for m, v in enumerate(p)]
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out

def build_h_S(k, S):
    p = [mp.mpf(1)]
    if 1 not in S: p = apply_Dj_mp(p, 2)
    for j in range(2, k + 1):
        p = apply_Dj_mp(p, j + 1)
        if j not in S: p = apply_Dj_mp(p, j + 1)
    return p

def build_A_S(k, S):
    p = [mp.mpf(1)]
    L = 0
    for j in range(1, k + 1):
        if j in S:
            p = apply_SL_mp(p, L); L += 2
        else:
            p = apply_SL_mp(apply_SL_mp(p, L), L + 1); L += 3
    return p

def check_interlace_bez(f, g):
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

# Precompute h0[k] and A0[k] and E[k] = A_{k,k}
h0 = {1: [mp.mpf(1), mp.mpf(2)]}
A0 = {1: [mp.mpf(1), mp.mpf(1)]}
Ek = {1: [mp.mpf(1)]}
for k in range(2, 30):
    h0[k] = apply_Dj_mp(apply_Dj_mp(h0[k-1], k+1), k+1)
    A0[k] = apply_SL_mp(apply_SL_mp(A0[k-1], 3*k-3), 3*k-2)
    Ek[k] = apply_SL_mp(Ek[k-1], 2*k-2)

print("=== 1. OPEN CASE i = 4: TEST COMMON UPPER INTERLACERS FOR ALL {A_{k, S} : |S|=4} ===")
for label, cand_fn in [
    ("D_{k+1} h_{k-2,0}", lambda k: apply_Dj_mp(h0[k-2], k+1)),
    ("D_{k+2} h_{k-2,0}", lambda k: apply_Dj_mp(h0[k-2], k+2)),
    ("D_{1.5k} h_{k-2,0}", lambda k: apply_Dj_mp(h0[k-2], mp.mpf('1.5')*k)),
    ("D_{2k} h_{k-2,0}", lambda k: apply_Dj_mp(h0[k-2], 2*k)),
    ("D_{3k} h_{k-2,0}", lambda k: apply_Dj_mp(h0[k-2], 3*k)),
]:
    res = []
    for k in [6, 10, 15, 20, 25]:
        early = build_h_S(k, {1, 2, 3, 4})
        late  = build_h_S(k, {k-3, k-2, k-1, k})
        M = cand_fn(k)
        ok_e = check_interlace_bez(early, M)
        ok_l = check_interlace_bez(late, M)
        res.append(f"k={k}:({'OK' if (ok_e and ok_l) else f'e={int(ok_e)},l={int(ok_l)}'})")
    print(f"  {label:20s}: " + " ".join(res))

print("\n=== 1B. OPEN CASE i = 4 IN x-VARIABLE: TEST A'_{k-1, 0}, P_{2k-3} A_{k-1, 0}, S_L A_{k-2, 0} ===")
for label, cand_fn in [
    ("A'_{k-1, 0}", lambda k: [m*A0[k-1][m] for m in range(1, len(A0[k-1]))]),
    ("P_{2k-3} A_{k-1,0}", lambda k: apply_Pc_mp(A0[k-1], 2*k - 3)),
    ("(1+x)A_{k-2,0}?no", lambda k: None),
]:
    if label.endswith("no"): continue
    res = []
    for k in [6, 10, 15, 20, 25]:
        early = build_A_S(k, {1, 2, 3, 4})
        late  = build_A_S(k, {k-3, k-2, k-1, k})
        M = cand_fn(k)
        ok_e = check_interlace_bez(early, M)
        ok_l = check_interlace_bez(late, M)
        res.append(f"k={k}:({'OK' if (ok_e and ok_l) else f'e={int(ok_e)},l={int(ok_l)}'})")
    print(f"  {label:20s}: " + " ".join(res))

print("\n=== 2. OPEN CASE i = k - 3: TEST P_{alpha*k} A_{k+2, k+2} << A_{k, S} (|S| = k-3) ===")
for label, c_fn in [
    ("c = 2k+2 (=2A_{k+1,k})", lambda k: 2*k + 2),
    ("c = 1.8k", lambda k: mp.mpf('1.8')*k),
    ("c = 1.6k", lambda k: mp.mpf('1.6')*k),
    ("c = 1.5k", lambda k: mp.mpf('1.5')*k),
    ("c = 1.4k", lambda k: mp.mpf('1.4')*k),
    ("c = 1.2k", lambda k: mp.mpf('1.2')*k),
    ("c = k+2",  lambda k: k + 2),
]:
    res = []
    for k in [6, 10, 15, 20, 25]:
        early = build_A_S(k, set(range(1, k - 2)))
        late  = build_A_S(k, set(range(4, k + 1)))
        M = apply_Pc_mp(Ek[k + 2], c_fn(k))
        ok_e = check_interlace_bez(M, early)
        ok_l = check_interlace_bez(M, late)
        res.append(f"k={k}:({'OK' if (ok_e and ok_l) else f'e={int(ok_e)},l={int(ok_l)}'})")
    print(f"  {label:24s}: " + " ".join(res))
