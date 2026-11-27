import numpy as np

# Adversarial test of the Pc 4-pair / combined U1+U2 << V1+V2 induction step!
# Suppose at level k-1 we have THREE polynomials f0, f1, f2 (playing A_{k-1,i}, A_{k-1,i-1}, A_{k-1,i-2})
# of degrees d, d+1, d+2 with all negative roots, satisfying:
#   IH1: f0 << f1 and f1 << f2
#   IH2: f0 << P_{L-3} f1 and f1 << P_{L-2} f2  (where L = 3k - 1 - i)
#   IH3 (optional, stronger): S_{L-2} f0 << P_L f2  (Pair 1,2)
# Question: Does U1 + U2 << V1 + V2 ALWAYS hold, where
#   U1 = S_{L-2} f0, U2 = f1
#   V1 = P_L S_{L-1} f1, V2 = P_L f2  ?!
# Especially when f1 has two VERY CLOSE roots!

def poly_add(a, b):
    n = max(len(a), len(b))
    res = np.zeros(n)
    res[:len(a)] += a
    res[:len(b)] += b
    return res

def apply_S(p, L):
    d = len(p) - 1
    res = np.zeros(d + 2)
    for n in range(d + 2):
        v = 0.0
        if n <= d: v += (n + 1) * p[n]
        if 1 <= n <= d + 1: v += (L - n + 1) * p[n - 1]
        res[n] = v
    return res

def apply_P(p, c):
    return np.array([(c - n) * p[n] for n in range(len(p))], dtype=np.float64)

def from_roots_neg(roots, scale=1.0):
    # prod(1 - x/r) or prod(x - r) with positive coeffs (r < 0)
    p = np.array([scale], dtype=np.float64)
    for r in roots:
        p = np.convolve(p, [-r, 1.0])
    return p

def check_interlace_deg_diff1(g, f):
    # deg g = deg f - 1, all negative roots: r_1(f) < r_1(g) < r_2(f) < ... < r_d(f) < 0
    rg = np.roots(g[::-1])
    rf = np.roots(f[::-1])
    if np.max(np.abs(rg.imag)) > 1e-7 or np.max(np.abs(rf.imag)) > 1e-7:
        return False, "complex"
    rg = np.sort(rg.real)
    rf = np.sort(rf.real)
    if rg[-1] >= 0 or rf[-1] >= 0:
        return False, "nonneg"
    ok = all(rf[m] < rg[m] < rf[m+1] for m in range(len(rg)))
    return ok, (rg, rf)

# Test 1: f1 has two close roots at -1.001, -0.999, with f0, f2 obeying IH1, IH2, IH3!
L = 10  # e.g. k=4, i=1 -> deg f0 = 3? Let's take d=1 (deg f0=1, f1=2, f2=3) or d=2 (deg f0=2, f1=3, f2=4), L=12
rng = np.random.default_rng(123)
for d, L in [(1, 8), (2, 11), (2, 15), (3, 15)]:
    valid_ih = 0
    valid_ih3 = 0
    fail_comb_ih = 0
    fail_comb_ih3 = 0
    fail_P21_ih3 = 0
    for trial in range(5000):
        # Generate roots of f2, f1, f0 interlacing; sometimes cluster two roots of f1
        if trial < 2500:
            # cluster a pair of roots of f1!
            pts = np.sort(-np.exp(rng.uniform(-3, 2, size=d+2)))
            # shrink gap between two adjacent roots of f1
            r2 = pts
            # pick r1 inside intervals of r2, with r1[0] and r1[1] very close to r2[1]!
            r1 = np.array([rng.uniform(r2[m], r2[m+1]) for m in range(d+1)])
            eps = 10**rng.uniform(-5, -1)
            if d + 1 >= 2:
                r1[0] = r2[1] - eps * (r2[1] - r2[0])
                r1[1] = r2[1] + eps * (r2[2] - r2[1])
            r0 = np.array([rng.uniform(r1[m], r1[m+1]) for m in range(d)])
        else:
            r2 = np.sort(-np.exp(rng.uniform(-3, 2, size=d+2)))
            r1 = np.array([rng.uniform(r2[m], r2[m+1]) for m in range(d+1)])
            r0 = np.array([rng.uniform(r1[m], r1[m+1]) for m in range(d)])

        # Scale f0, f1, f2 randomly (positive)
        s0 = np.exp(rng.uniform(-3, 3))
        s1 = np.exp(rng.uniform(-3, 3))
        s2 = np.exp(rng.uniform(-3, 3))
        f0 = from_roots_neg(r0, s0)
        f1 = from_roots_neg(r1, s1)
        f2 = from_roots_neg(r2, s2)

        # Check IH2: f0 << P_{L-3} f1 and f1 << P_{L-2} f2
        ok1, _ = check_interlace_deg_diff1(f0, apply_P(f1, L - 3))
        if not ok1: continue
        ok2, _ = check_interlace_deg_diff1(f1, apply_P(f2, L - 2))
        if not ok2: continue
        valid_ih += 1

        U1 = apply_S(f0, L - 2)
        U2 = f1
        V1 = apply_P(apply_S(f1, L - 1), L)
        V2 = apply_P(f2, L)

        # Check IH3: U1 << V2
        ok_ih3, _ = check_interlace_deg_diff1(U1, V2)
        if ok_ih3: valid_ih3 += 1

        ok_comb, _ = check_interlace_deg_diff1(poly_add(U1, U2), poly_add(V1, V2))
        if not ok_comb:
            fail_comb_ih += 1
            if ok_ih3:
                fail_comb_ih3 += 1
        if ok_ih3:
            ok_p21, _ = check_interlace_deg_diff1(U2, V1)
            if not ok_p21:
                fail_P21_ih3 += 1
    print(f"d={d}, L={L:2d}: valid_IH={valid_ih:4d} (fail_comb={fail_comb_ih:4d}) | valid_IH3={valid_ih3:4d} (fail_P21={fail_P21_ih3:4d}, fail_comb={fail_comb_ih3:4d})")
