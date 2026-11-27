import numpy as np
import sympy as sp

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0]*n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_S(p, L):
    d = len(p) - 1
    res = [0]*(d + 2)
    for n in range(d + 2):
        v = 0
        if n <= d: v += (n + 1) * p[n]
        if 1 <= n <= d + 1: v += (L - n + 1) * p[n - 1]
        res[n] = v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def roots_real(p):
    c = np.array(p[::-1], dtype=np.float64)
    r = np.roots(c)
    if np.max(np.abs(r.imag)) > 1e-6:
        return None
    return np.sort(r.real)

KMAX = 12
A = {(0, 0): [1]}
for k in range(1, KMAX + 1):
    for i in range(k + 1):
        u1 = apply_S(A[(k-1, i)], 3*k - 3 - i) if i <= k - 1 else [0]
        u2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S(poly_add(u1, u2), 3*k - 2 - i)

# 1. Check EQUAL-DEGREE interlacing along DIAGONALS: A_{k, i} vs A_{k+1, i+2}!!
#    WAIT A SECOND! They have the EXACT SAME DEGREE 2k - i - 1!!
#    And at i = k - 1, A_{k, k-1} = (1/2) P_{2k} A_{k+1, k+1}, so A_{k+1, k+1} << A_{k, k-1} (or vice versa)!
#    WAIT A SECOND!! Does A_{k, i} << A_{k+1, i+2} (or A_{k+1, i+2} << A_{k, i}) for ALL (k, i)?!
print("=== 1. EQUAL-DEGREE DIAGONAL INTERLACING: A_{k, i} vs A_{k+1, i+2} ===")
ok_LR = 0; ok_RL = 0; fail = 0
for k in range(1, KMAX):
    for i in range(k + 1):
        if i + 2 <= k + 1:
            f = A[(k, i)]
            g = A[(k+1, i+2)]
            rf = roots_real(f); rg = roots_real(g)
            assert len(rf) == len(rg)
            if len(rf) == 0: continue
            lr = all(rf[m] < rg[m] for m in range(len(rf))) and all(rg[m] < rf[m+1] for m in range(len(rf)-1))
            rl = all(rg[m] < rf[m] for m in range(len(rf))) and all(rf[m] < rg[m+1] for m in range(len(rf)-1))
            if lr: ok_LR += 1
            elif rl: ok_RL += 1
            else:
                fail += 1
                print(f"  FAIL at k={k}, i={i}: rf={np.round(rf,3)}, rg={np.round(rg,3)}")
print(f"  A_{{k,i}} left of A_{{k+1,i+2}}: {ok_LR} | right of: {ok_RL} | fail: {fail}")

# 2. Check Brändén-Solus SYMMETRIC DECOMPOSITION of A_{k,i}(x)!
#    Write A_{k,i}(x) = a(x) + x b(x) where a(x) is palindromic of degree d = 2k-i-1, b(x) is palindromic of degree d-1!
#    Recall: A(x) = a(x) + x b(x) uniquely with x^d a(1/x) = a(x), x^{d-1} b(1/x) = b(x):
#      a(x) = (A(x) - x^{d+1} A(1/x)) / (1 - x) ?? Wait:
#      a(x) + x b(x) = A(x),  a(x) + b(x) = x^d A(1/x) (check: x^d(a(1/x) + (1/x)b(1/x)) = a(x) + b(x)!)
#      => (1 - x) b(x) = x^d A(1/x) - A(x),  (1 - x) a(x) = A(x) - x^{d+1} A(1/x)!
#    Branden-Solus Theorem: If a(x), b(x) are real-rooted and b << a (or a << b), then A(x) is real-rooted!
#    AND we can check if S_L preserves the symmetric-decomposition interlacing (b << a)!!
print("\n=== 2. BRANDEN-SOLUS SYMMETRIC DECOMPOSITION A_{k,i}(x) = a(x) + x b(x) ===")
def sym_decomp(p):
    # p ascending of deg d
    d = len(p) - 1
    if d == 0: return p[:], [0]
    p_rev = p[::-1]
    # (1-x) b(x) = p_rev - p  -> b[n] = sum_{j=0}^n (p_rev[j] - p[j])
    b = []
    acc = 0
    for j in range(d):
        acc += p_rev[j] - p[j]
        b.append(acc)
    assert acc + (p_rev[d] - p[d]) == 0
    # a(x) = p(x) - x b(x)
    a = [p[0]] + [p[j] - b[j-1] for j in range(1, d+1)]
    return a, b

sd_pass = 0; sd_fail = 0
for k in range(2, 11):
    for i in range(k + 1):
        p = A[(k, i)]
        if len(p) <= 2: continue
        a, b = sym_decomp(p)
        ra = roots_real(a); rb = roots_real(b)
        if ra is None or rb is None:
            sd_fail += 1
            print(f"  SymDecomp COMPLEX at k={k}, i={i}: a_real={ra is not None}, b_real={rb is not None}")
            continue
        ok = all(ra[m] < rb[m] < ra[m+1] for m in range(len(rb)))
        if ok: sd_pass += 1
        else:
            sd_fail += 1
            print(f"  SymDecomp NO-INTERLACE at k={k}, i={i}: b<<a={ok}")
print(f"  Symmetric decomposition b << a (k=2..10): pass={sd_pass}, fail={sd_fail}")

# 3. Check INTERLACING ALONG THE OTHER DIAGONAL: A_{k,i} vs A_{k-1, i-1} (same diff k-i) and A_{k,i} vs A_{k-1, i}
print("\n=== 3. VERTICAL / DIAGONAL INTERLACING: A_{k-1,i-1} vs A_{k,i} (deg diff 1) ===")
v_pass = 0; v_fail = 0
for k in range(2, 12):
    for i in range(1, k + 1):
        f = A[(k-1, i-1)]
        g = A[(k, i)]
        rf = roots_real(f); rg = roots_real(g)
        assert len(rg) == len(rf) + 1
        ok = all(rg[m] < rf[m] < rg[m+1] for m in range(len(rf)))
        if ok: v_pass += 1
        else:
            v_fail += 1
            print(f"  FAIL A_{{{k-1},{i-1}}} << A_{{{k},{i}}}")
print(f"  A_{{k-1,i-1}} << A_{{k,i}}: pass={v_pass}, fail={v_fail}")

# 4. CHECK A BRAND-NEW IDEA: REFINE A_{k,i}(x) BY THE FIRST OR LAST LETTER (OR NUMBER OF ASCENTS/RUNS) -> INTERLACING MATRIX!
#    Wait: in GLZ, A_{k,i} = sum_{|S|=i} A_{k,S}.
#    We ALREADY know the {A_{k,S}} themselves do NOT have a common interlacer (SOURCES_AND_STATE Section 4).
#    WHAT ABOUT refining A_{k,i}(x) by the POSITION OF THE LAST ENTRY (or splitting by k vs k-1, i.e., the insertion vector)?
#    Wait: what IS the combinatorial recursion for A_{k,i}(x) by WORD POSITION / DESCENT / RUNS?
#    We have the recurrence a_{k,i,j} = j^2 a_{k-1,i,j} + ... (a 5-term recurrence in j)!
#    Wait: that's the SAME as S_{3k-2-i}(S_{3k-3-i} A_{k-1,i} + A_{k-1,i-1})!
#    WAIT A SECOND! Check if B_{k,i} = S_{3k-3-i} A_{k-1,i} + A_{k-1,i-1} satisfies symmetric decomposition b << a!
print("\n=== 4. CHECK SYMMETRIC DECOMPOSITION OF B_{k,i} AND OF S_L(f) IN GENERAL! ===")
bs_fail = 0; bs_pass = 0
for k in range(2, 11):
    for i in range(1, k):
        B = poly_add(apply_S(A[(k-1, i)], 3*k - 3 - i), A[(k-1, i-1)])
        if len(B) <= 2: continue
        a, b = sym_decomp(B)
        ra = roots_real(a); rb = roots_real(b)
        if ra is not None and rb is not None and all(ra[m] < rb[m] < ra[m+1] for m in range(len(rb))):
            bs_pass += 1
        else:
            bs_fail += 1
            print(f"  B_{{{k},{i}}} SymDecomp FAIL!")
print(f"  B_{{k,i}} SymDecomp b << a (k=2..10): pass={bs_pass}, fail={bs_fail}")

# 5. Check if A_{k,i}(x) has ALL-REAL ROOTS in the gamma-basis (Wait: A_{k,i} is NOT palindromic,
#    so check the (a(x), b(x)) gamma-positivity / real-rootedness of a_{k,i} << a_{k,i-1} and b_{k,i} << b_{k,i-1}!)
print("\n=== 5. INTERLACING OF SYMMETRIC DECOMPOSITION PARTS (a_{k,i}, b_{k,i}) ACROSS i! ===")
a_inter_fail = 0; b_inter_fail = 0; ab_cross_fail = 0; cnt5 = 0
for k in range(3, 11):
    for i in range(1, k):
        a1, b1 = sym_decomp(A[(k, i)])
        a0, b0 = sym_decomp(A[(k, i-1)])
        if len(b1) <= 1: continue
        cnt5 += 1
        ra1, rb1 = roots_real(a1), roots_real(b1)
        ra0, rb0 = roots_real(a0), roots_real(b0)
        if not all(ra0[m] < ra1[m] < ra0[m+1] for m in range(len(ra1))):
            a_inter_fail += 1
        if not all(rb0[m] < rb1[m] < rb0[m+1] for m in range(len(rb1))):
            b_inter_fail += 1
        # Check cross: b1 << a0? (deg b1 = d-1, deg a0 = d+1 -> deg diff 2)
        # Check a1 << b0 (EQUAL DEGREE d!!)!
        ok_eq_LR = all(ra1[m] < rb0[m] for m in range(len(ra1))) and all(rb0[m] < ra1[m+1] for m in range(len(ra1)-1))
        ok_eq_RL = all(rb0[m] < ra1[m] for m in range(len(ra1))) and all(ra1[m] < rb0[m+1] for m in range(len(ra1)-1))
        if not (ok_eq_LR or ok_eq_RL):
            ab_cross_fail += 1
print(f"  (k=3..10, n={cnt5}): a_{{k,i}}<<a_{{k,i-1}} fails={a_inter_fail}, b_{{k,i}}<<b_{{k,i-1}} fails={b_inter_fail}, a_{{k,i}}||b_{{k,i-1}} fails={ab_cross_fail}")
