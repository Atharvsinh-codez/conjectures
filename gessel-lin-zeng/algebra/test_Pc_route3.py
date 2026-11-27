"""
Investigate Route III: The P_c = c - x d/dx intertwining invariant!
Remember @verifier (#23) verified 100% through k=22:
  (I_k): A_{k,i} << P_{3k-1-i} A_{k,i-1}  (253/253 PASS)
  (P21): S_{3k-3-i} A_{k-1,i} << P_{3k-1-i} A_{k-1,i-2}  (210/210 PASS)
Wait: let's check:
1. Does f << P_c g (deg g = deg f + 1, c > deg g) ALWAYS imply S_{c+1} f << g?
   (Or S_L f << g for what L?)
2. What is the EXACT relation between P_c and S_L on A_{k,i}, and what is the
   strongest c = c*(k,i) for which A_{k,i} << P_c A_{k,i-1} holds?
   Wait: as c decreases toward deg(A_{k,i-1}), P_c shifts roots FURTHER left (STRONGER!).
   Wait: check that! P_c g = c g - x g'. At x < 0, larger c pushes P_c g closer to g!
   So SMALLER c is STRONGER!
   Let's check the smallest c for which A_{k,i} << P_c A_{k,i-1}, AND check whether
   A_{k-1,i} << P_c A_{k-1,i-1} closes a 1-step induction!
"""
import numpy as np

def poly_trim(p):
    p = list(p)
    while len(p) > 1 and abs(p[-1]) == 0:
        p.pop()
    return p

def apply_SL(p, L):
    # (1 + Lx) p(x) + x(1 - x) p'(x)
    # p[m] x^m -> (1+m) p[m] x^m + (L - m) p[m] x^{m+1}
    out = [0.0] * (len(p) + 1)
    for m, c in enumerate(p):
        out[m] += (1 + m) * c
        out[m + 1] += (L - m) * c
    return poly_trim(out)

def apply_Pc(p, c):
    # c p(x) - x p'(x) -> (c - m) p[m] x^m
    out = [(c - m) * val for m, val in enumerate(p)]
    return poly_trim(out)

def roots_neg(p):
    # returns sorted real roots (increasing) or None if non-real
    r = np.roots(list(reversed(p)))
    if np.max(np.abs(r.imag)) > 1e-7:
        return None
    return np.sort(r.real)

def interlaces(p, q):
    """
    Check p << q (Ma-Wang p prec q) for polynomials with negative real roots and positive lc:
    deg q in {deg p, deg p + 1};
    roots of p: r_1 < ... < r_m; roots of q: s_1 < ... < s_n (n in {m, m+1}).
    If n == m + 1: s_1 < r_1 < s_2 < r_2 < ... < r_m < s_{m+1}.
    If n == m:     r_1 < s_1 < r_2 < s_2 < ... < r_m < s_m.
    (In both cases, rightmost root of q is > rightmost root of p, and roots strictly alternate.)
    """
    rp = roots_neg(p)
    rq = roots_neg(q)
    if rp is None or rq is None:
        return False
    m, n = len(rp), len(rq)
    if n == m + 1:
        for k in range(m):
            if not (rq[k] < rp[k] < rq[k + 1]):
                return False
        return True
    elif n == m:
        for k in range(m):
            if not (rp[k] < rq[k]):
                return False
            if k + 1 < m and not (rq[k] < rp[k + 1]):
                return False
        return True
    return False

# 1. Test universal implication: if f << P_c g (deg g = deg f + 1 = d, c > d),
#    does S_L f << g ALWAYS hold? For what L?
print("=== TEST 1: Does f << P_c g ALWAYS imply S_L f << g? ===")
rng = np.random.default_rng(7)
for d in [4, 6, 8]:
    c = d + 3  # like GLZ: deg g = 2k-i-2, c = 3k-4-i, so c - d = k - 2
    for L_off in [-1, 0, 1, 2]:
        L = c + L_off
        fails = 0
        trials = 500
        for _ in range(trials):
            # generate random negative roots for g of deg d
            sg = -np.sort(np.exp(rng.normal(0, 1.5, d)))[::-1]  # increasing negative
            g = np.poly(sg)[::-1]  # ascending coeffs, wait: np.poly has leading 1 -> all positive!
            Pg = apply_Pc(g, c)
            rPg = roots_neg(Pg)
            # pick f of deg d-1 strictly interlacing Pg: s'_1 < r_1 < s'_2 < ... < r_{d-1} < s'_d
            rf = [rng.uniform(rPg[k], rPg[k+1]) for k in range(d - 1)]
            f = np.poly(rf)[::-1]
            Sf = apply_SL(f, L)
            if not interlaces(Sf, g):
                fails += 1
        print(f"  d={d}, c={c}, L={L} (L-c={L_off:+d}): fails = {fails}/{trials}")

# 2. Test whether P_c g << g and what S_L f << g vs f << P_c g really is!
# Wait: for A_{k,i}, L_step = 3k - 3 - i, while c_prev = 3(k-1) - 1 - i = 3k - 4 - i = L_step - 1!
# So L_step = c_prev + 1!
# Wait: check whether S_{c+1} f << g holds when c = d + (k-2)!
# Let's also check if f << P_c g ALWAYS implies S_{L-1} f << P_{L+1} S_L g?
print("\n=== TEST 2: Check ALL 4 pairs of (S_{L-1} f + g) vs (P_{L+1} S_L g + P_{L+1} h) ===")
# Given A_{k-1,i} = f (deg d-1), A_{k-1,i-1} = g (deg d), A_{k-1,i-2} = h (deg d+1),
# with L = 3k - 2 - i, c_prev = 3k - 4 - i = L - 2, and for (g, h): c_prev_gh = 3k - 3 - i = L - 1.
# Hypothesis at k-1:
#   f << P_{L-2} g   AND   g << P_{L-1} h.
# Wait: what is the relation between S_{L-1} f + g and P_{L+1}(S_L g + h) = P_{L+1} S_L g + P_{L+1} h?
# Let's check for random (f, g, h) satisfying f << P_{L-2} g and g << P_{L-1} h:
#   Pair (1,1): S_{L-1} f << P_{L+1} S_L g ?
#   Pair (1,2): S_{L-1} f << P_{L+1} h ?
#   Pair (2,1): g << P_{L+1} S_L g ?
#   Pair (2,2): g << P_{L+1} h ?  (AUTOMATIC since g << P_{L-1} h << P_{L+1} h!)
for k_sim, i_sim in [(6, 3), (8, 4), (10, 5)]:
    d = 2*k_sim - i_sim - 2  # deg g
    L = 3*k_sim - 2 - i_sim
    counts = {'11': 0, '12': 0, '21': 0, '22': 0, 'U_real': 0, 'U_ll_V': 0, 'Sf_ll_g': 0, 'Sg_ll_h': 0}
    trials = 400
    for _ in range(trials):
        sh = -np.sort(np.exp(rng.normal(0, 1.5, d + 1)))[::-1]
        h = np.poly(sh)[::-1]
        Ph = apply_Pc(h, L - 1)
        rPh = roots_neg(Ph)
        sg = [rng.uniform(rPh[m], rPh[m+1]) for m in range(d)]
        g = np.poly(sg)[::-1]
        Pg = apply_Pc(g, L - 2)
        rPg = roots_neg(Pg)
        sf = [rng.uniform(rPg[m], rPg[m+1]) for m in range(d - 1)]
        f = np.poly(sf)[::-1]

        T1 = apply_SL(f, L - 1)
        T2 = g
        T3 = apply_Pc(apply_SL(g, L), L + 1)
        T4 = apply_Pc(h, L + 1)
        U = [T1[m] + T2[m] for m in range(len(T2))]
        V = [T3[m] + T4[m] for m in range(len(T4))]

        if not interlaces(T1, T3): counts['11'] += 1
        if not interlaces(T1, T4): counts['12'] += 1
        if not interlaces(T2, T3): counts['21'] += 1
        if not interlaces(T2, T4): counts['22'] += 1
        if not interlaces(T1, T2): counts['Sf_ll_g'] += 1
        if not interlaces(apply_SL(g, L), h): counts['Sg_ll_h'] += 1
        if roots_neg(U) is None: counts['U_real'] += 1
        if not interlaces(U, V): counts['U_ll_V'] += 1
    print(f"  k={k_sim}, i={i_sim} (d={d}, L={L}): fails out of {trials}: {counts}")

# 3. Also check on ACTUAL GLZ polynomials: what is the SMALLEST c_min(k,i) such that A_{k,i} << P_c A_{k,i-1}?
print("\n=== TEST 3: Smallest c such that A_{k,i} << P_c A_{k,i-1} on GLZ ===")
# Build exact A_{k,i} (float is fine for k<=11)
A = {(1, 0): [1.0, 2.0], (1, 1): [1.0]}
for k in range(2, 12):
    for i in range(k + 1):
        L = 3*k - 2 - i
        u = [0.0]
        if i <= k - 1:
            t1 = apply_SL(A[(k-1, i)], L - 1)
            u = [0.0]*max(len(u), len(t1))
            for m, v in enumerate(t1): u[m] += v
        if i >= 1:
            t2 = A[(k-1, i-1)]
            if len(t2) > len(u): u += [0.0]*(len(t2)-len(u))
            for m, v in enumerate(t2): u[m] += v
        A[(k, i)] = apply_SL(u, L)

for k in range(3, 12):
    worst_gap = 1e9
    worst_i = -1
    for i in range(1, k + 1):
        d_prev = 2*k - i  # deg A_{k, i-1}
        c_glz = 3*k - 1 - i  # = d_prev + k - 1
        # bisect c in (d_prev, c_glz + 5]
        lo = d_prev + 1e-6
        hi = c_glz
        if not interlaces(A[(k, i)], apply_Pc(A[(k, i-1)], hi)):
            print(f"  FAIL at k={k}, i={i}!")
            continue
        for _ in range(40):
            mid = 0.5*(lo + hi)
            if interlaces(A[(k, i)], apply_Pc(A[(k, i-1)], mid)):
                hi = mid
            else:
                lo = mid
        gap = c_glz - hi
        if gap < worst_gap:
            worst_gap = gap
            worst_i = i
    print(f"  k={k:2d}: worst (c_glz - c_min) = {worst_gap:.6f} at i={worst_i} (c_glz={3*k-1-worst_i}, d_prev={2*k-worst_i}, c_min={3*k-1-worst_i-worst_gap:.4f})")
