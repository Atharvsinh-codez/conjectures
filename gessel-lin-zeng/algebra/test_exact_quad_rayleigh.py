import math
import numpy as np

def poly_trim(c):
    while len(c) > 1 and c[-1] == 0: c.pop()
    return c

def apply_D(j, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return poly_trim(out)

def add_poly(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return poly_trim(out)

def sub_poly(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] -= b[i]
    return poly_trim(out)

def mul_poly(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b): out[i + j] += ai * bj
    return poly_trim(out)

def deriv(p):
    return poly_trim([m * p[m] for m in range(1, len(p))]) if len(p) > 1 else [0]

def wronskian(f, g):
    # W(f, g) = f*g' - f'*g
    return sub_poly(mul_poly(f, deriv(g)), mul_poly(deriv(f), g))

def to_x_poly(h_coeffs, d):
    # Convert h(t) of degree <= d to (1-x)^d h(x/(1-x))
    pows = [[1]]
    for r in range(1, d + 1):
        pows.append(mul_poly(pows[-1], [1, -1]))
    out = [0] * (d + 1)
    for m, cm in enumerate(h_coeffs):
        pw = pows[d - m]
        for j, pj in enumerate(pw):
            out[m + j] += cm * pj
    return poly_trim(out)

def count_real_roots_sturm(p):
    """Exact Sturm count of DISTINCT real roots of integer poly p (low->high) on (-inf, +inf)."""
    p = p[:]
    while len(p) > 1 and p[-1] == 0: p.pop()
    if len(p) <= 1: return 0
    if p[-1] < 0: p = [-v for v in p]
    # Make squarefree first? Or just run PRS and divide by last nonzero (gcd)!
    def prem_neg(a, b):
        da, db = len(a) - 1, len(b) - 1
        cur = a[:]
        lc_b = b[-1]
        for shift in range(da - db, -1, -1):
            # multiply cur by |lc_b|? Careful: must multiply by POSITIVE factor to preserve signs!
            # Use lc_b^2? Wait: if we do one step with |lc_b|, if lc_b < 0 we adjust!
            c = cur[shift + db]
            if c == 0: continue
            if lc_b > 0:
                cur = [lc_b * v for v in cur]
                for m in range(db + 1):
                    cur[shift + m] -= c * b[m]
            else:
                cur = [(-lc_b) * v for v in cur]
                for m in range(db + 1):
                    cur[shift + m] += c * b[m]
        while len(cur) > 1 and cur[-1] == 0: cur.pop()
        if len(cur) == 1 and cur[0] == 0: return None
        nxt = [-v for v in cur]
        g = math.gcd(*[abs(v) for v in nxt])
        if g > 1: nxt = [v // g for v in nxt]
        return nxt
    chain = [p, deriv(p)]
    g = math.gcd(*[abs(v) for v in chain[1]])
    if g > 1: chain[1] = [v // g for v in chain[1]]
    while len(chain[-1]) > 1:
        nxt = prem_neg(chain[-2], chain[-1])
        if nxt is None: break
        chain.append(nxt)
    # Sign changes at -inf and +inf
    def sign_changes(seq):
        nz = [s for s in seq if s != 0]
        return sum(1 for a, b in zip(nz, nz[1:]) if a * b < 0)
    s_pos = [1 if q[-1] > 0 else -1 for q in chain]
    s_neg = [(1 if q[-1] > 0 else -1) * (1 if (len(q)-1) % 2 == 0 else -1) for q in chain]
    return sign_changes(s_neg) - sign_changes(s_pos)

# Compute ht up to k=14
KMAX = 14
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, KMAX + 1):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

# === TEST 1: EXACT 100%-RIGOROUS REAL STABILITY CERTIFICATE FOR Quad Q_{k,i}(t, y) = X - 2My + Yy^2! ===
# Check: 4*W(X,M)*W(M,Y) - W(X,Y)^2 has ZERO real roots and POSITIVE leading coeff!
print("=== TEST 1: EXACT 100%-RIGOROUS REAL STABILITY OF Quad X - 2My + Yy^2 via Sturm on 4*W(X,M)*W(M,Y) - W(X,Y)^2 ===")
for k in range(3, 12):
    fails = []
    min_x_coeff = 10**99
    min_t_coeff = 10**99
    for i in range(2, k):
        X = ht[(k-1, i-2)]
        M = apply_D(k+1, ht[(k-1, i-1)])
        Y = apply_D(k+1, apply_D(k+1, ht[(k-1, i)]))
        wxm = wronskian(X, M)
        wmy = wronskian(M, Y)
        wxy = wronskian(X, Y)
        disc = sub_poly([4 * v for v in mul_poly(wxm, wmy)], mul_poly(wxy, wxy))
        nroots = count_real_roots_sturm(disc)
        lc_pos = disc[-1] > 0
        if nroots > 0 or not lc_pos:
            fails.append((i, nroots, lc_pos))
        # Also check coefficients of disc in t and in x!
        min_t_coeff = min(min_t_coeff, min(disc))
        # degree of disc <= 4*(deg X - 1)
        dx = to_x_poly(disc, 4*(len(X) - 2))
        min_x_coeff = min(min_x_coeff, min(dx))
    print(f"  k={k:2d}: EXACT REAL-STABLE for ALL i=2..{k-1}: {'YES (0 real roots)!' if not fails else f'FAIL {fails}'} | min_t_coeff>=0: {min_t_coeff>=0} | min_x_coeff>=0: {min_x_coeff>=0}")

# === TEST 2: Check EXACT 100%-rigorous real stability of Quad for larger k=12, 13, 14 (check middle i and boundary i=2) ===
print("\n=== TEST 2: Check Quad at k=12, 13, 14 (fast numeric + exact check on boundary/middle i) ===")
t_R = np.concatenate([-10.0**np.linspace(3, -5, 100), [0.0], 10.0**np.linspace(-5, 3, 100), np.linspace(-1.5, 0.5, 400)])
for k in range(12, 15):
    fails = []
    for i in range(2, k):
        X = ht[(k-1, i-2)]
        M = apply_D(k+1, ht[(k-1, i-1)])
        Y = apply_D(k+1, apply_D(k+1, ht[(k-1, i)]))
        wxm = wronskian(X, M)
        wmy = wronskian(M, Y)
        wxy = wronskian(X, Y)
        disc = sub_poly([4 * v for v in mul_poly(wxm, wmy)], mul_poly(wxy, wxy))
        # Check sign on t_R and exact Sturm for i=2 and i=k//2
        vals = np.zeros_like(t_R)
        for c in reversed(disc): vals = vals * t_R + float(c)
        if np.min(vals) <= 0:
            fails.append((i, float(np.min(vals))))
        if i in (2, k//2, k-1):
            nr = count_real_roots_sturm(disc)
            if nr > 0 or disc[-1] <= 0:
                fails.append((i, f"sturm={nr}"))
    print(f"  k={k:2d}: Quad real-stable check: {'PASS!' if not fails else f'FAIL {fails}'}")

# === TEST 3: Check my prediction that B_10(x, -x*y) FAILS at x = -0.001 + 1e-7j! ===
print("\n=== TEST 3: Verify prediction that B_10(x, -x*y) fails at x = -0.001 + 1e-7j ===")
def eval_c(p, z):
    acc = 0j
    for c in reversed(p): acc = acc * z + c
    return acc
def apply_S(L, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (L - m) * cm
    return poly_trim(out)
A = {(0, 0): [1]}
for k in range(1, 11):
    for i in range(k + 1):
        t1 = apply_S(3*k - 3 - i, A[(k-1, i)]) if i <= k - 1 else [0]
        t2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S(3*k - 2 - i, add_poly(t1, t2))
xv = complex(-0.001, 1e-7)
yc = [math.comb(10, i) * eval_c(A[(10, i)], xv) * ((-xv)**i) for i in range(10, -1, -1)]
yc = [c / max(abs(z) for z in yc) for c in yc]
print(f"  B_10(x, -x*y) at x={xv}: max Im(y) = {max(r.imag for r in np.roots(yc)):+.4e} (>0 confirms FAIL!)")

# === TEST 4: Check fine points near t=0 and t=-1 for binom-E(t, -y) at k=10..14! ===
print("\n=== TEST 4: Stress test binom-E(t, -y) at fine t near 0, -1, (-1,0) for k=10..14 ===")
fine_ts = [complex(re, 1e-7) for re in [-1.01, -1.001, -0.999, -0.99, -0.5, -0.1, -0.01, -0.001, -1e-4, 1e-4, 0.001, 0.01, 0.1]]
for k in range(10, 15):
    E = []
    for i in range(k):
        p = ht[(k-1, i)][:]
        for _ in range(i): p = apply_D(k + 1, p)
        E.append(p)
    min_im = 1e9
    worst = None
    for tv in fine_ts:
        yc = [math.comb(k-1, i) * eval_c(E[i], tv) for i in range(k-1, -1, -1)]
        yc = [c / max(abs(z) for z in yc) for c in yc]
        mi = min(r.imag for r in np.roots(yc))
        if mi < min_im:
            min_im = mi; worst = tv
    print(f"  k={k:2d}: binom-E(t, -y) fine minIm(y) = {min_im:+.4e} at t={worst}")
