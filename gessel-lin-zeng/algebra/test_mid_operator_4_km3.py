"""
3-second check:
1. For i = 4: test M_k^+ = prod_{j=5}^k (D_{j - alpha}^2)(q_0) for alpha = 0.25, 0.5, 0.75, 1.0,
   AND prod_{j=5}^k (D_j D_{j-1})(q_0) (exact integer operators!) up to k = 45!
   AND check `late << early` up to k = 45 for i = 4, i = 5, and i = k - 3!
2. Also check i = k - 3: test N_k = prod_{j=4}^k D_{j + alpha}(q_0) for alpha = 0, 0.25, 0.5, 0.75, 1.0, 1.5!
"""
import mpmath as mp
mp.mp.dps = 100

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

h0_2 = apply_Dj_mp(apply_Dj_mp([mp.mpf(1), mp.mpf(2)], 3), 3)
h0_3 = apply_Dj_mp(apply_Dj_mp(h0_2, 4), 4)
h4_4 = build_h_S(4, {1, 2, 3, 4})

print("=== 1. OPEN CASE i = 4: TEST prod_{j=5}^k D_{j-alpha}^2 (D_3 h_{2,0}) AND D_j D_{j-1} UP TO k = 45 ===")
ks = [8, 15, 25, 35, 45]
earlys = {k: build_h_S(k, {1, 2, 3, 4}) for k in ks}
lates  = {k: build_h_S(k, {k-3, k-2, k-1, k}) for k in ks}
print("  late << early (i=4): " + " ".join(f"k={k}:{check_interlace_bez(lates[k], earlys[k])}" for k in ks))

for label, step_fn in [
    ("D_{j-0.5}^2", lambda p, j: apply_Dj_mp(apply_Dj_mp(p, j - mp.mpf('0.5')), j - mp.mpf('0.5'))),
    ("D_j D_{j-1}", lambda p, j: apply_Dj_mp(apply_Dj_mp(p, j - 1), j)),
    ("D_{j-1} D_j", lambda p, j: apply_Dj_mp(apply_Dj_mp(p, j), j - 1)),
    ("D_{j-0.75}^2", lambda p, j: apply_Dj_mp(apply_Dj_mp(p, j - mp.mpf('0.75')), j - mp.mpf('0.75'))),
    ("D_{j-1}^2 (h_{k-2,0} style)", lambda p, j: apply_Dj_mp(apply_Dj_mp(p, j - 1), j - 1)),
]:
    M = apply_Dj_mp(h0_2, 3)  # deg 4
    res = []
    for j in range(5, 46):
        M = step_fn(M, j)
        if j in ks:
            oe = check_interlace_bez(earlys[j], M)
            ol = check_interlace_bez(lates[j], M)
            res.append(f"k={j}:{'OK' if (oe and ol) else f'e{int(oe)}l{int(ol)}'}")
    print(f"  {label:28s}: " + " ".join(res))

print("\n=== 2. CHECK late << early FOR i = 5 AND i = k - 3 UP TO k = 45 ===")
for k in [8, 15, 25, 35, 45]:
    e5 = build_h_S(k, {1, 2, 3, 4, 5}); l5 = build_h_S(k, {k-4, k-3, k-2, k-1, k})
    ek3 = build_h_S(k, set(range(1, k-2))); lk3 = build_h_S(k, set(range(4, k+1)))
    print(f"  k={k:2d}: i=5 late<<early: {check_interlace_bez(l5, e5)} | i=k-3 late<<early: {check_interlace_bez(lk3, ek3)}")

print("\n=== 3. OPEN CASE i = k - 3: TEST N_k = prod_{j=4}^k D_{j + alpha} (D_4 h_{2,0}) UP TO k = 45 ===")
earlys3 = {k: build_h_S(k, set(range(1, k-2))) for k in ks}
lates3  = {k: build_h_S(k, set(range(4, k+1))) for k in ks}
for alpha in [1.0, 1.25, 1.5, 1.75, 2.0]:
    N = apply_Dj_mp(h0_2, 4)  # deg 4
    res = []
    for j in range(4, 46):
        N = apply_Dj_mp(N, j + mp.mpf(str(alpha)))
        if j in ks:
            oe = check_interlace_bez(N, earlys3[j])
            ol = check_interlace_bez(N, lates3[j])
            res.append(f"k={j}:{'OK' if (oe and ol) else f'e{int(oe)}l{int(ol)}'}")
    print(f"  D_{{j+{alpha:.2f}}}: " + " ".join(res))
