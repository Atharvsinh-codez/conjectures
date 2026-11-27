"""Fast exact Cauchy-index & Sturm test of ALL closure candidates for Gessel-Lin-Zeng!"""
import sys
import time
import sympy as sp

t = sp.Symbol('t')

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for idx, v in enumerate(a): res[idx] += v
    for idx, v in enumerate(b): res[idx] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def poly_sub(a, b):
    return poly_add(a, [-v for v in b])

def poly_mul(a, b):
    res = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            res[i + j] += u * v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_Dtilde(j, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(poly):
        res[m] += c * (m + 1)
        res[m + 1] += c * (m + j + 1)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def build_p(kmax):
    p = [[None] * (k + 1) for k in range(kmax + 1)]
    p[1][0] = [1, 2]
    p[1][1] = [1]
    for k in range(2, kmax + 1):
        p[k][0] = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][0]))
        for i in range(1, k):
            p[k][i] = apply_Dtilde(k, poly_add(apply_Dtilde(k, p[k - 1][i]), p[k - 1][i - 1]))
        p[k][k] = apply_Dtilde(k, p[k - 1][k - 1])
    return p

def exact_cauchy_interlace(A, B):
    """Exact Cauchy index I(B/A) via integer primitive pseudo-remainder sequence."""
    p0 = sp.Poly(list(reversed(A)), t, domain='ZZ')
    p1 = sp.Poly(list(reversed(B)), t, domain='ZZ')
    da, db = p0.degree(), p1.degree()
    if da == 0 or db == 0:
        return True, "trivial", 0
    if p0.LC() < 0: p0 = -p0
    if p1.LC() < 0: p1 = -p1
    chain = [p0, p1]
    while chain[-1].degree() > 0:
        a, b = chain[-2], chain[-1]
        mdeg = a.degree() - b.degree() + 1
        r = sp.prem(a, b)
        if r.is_zero:
            return False, "common_root", 0
        if b.LC() < 0 and (mdeg % 2 == 1):
            r = -r
        r = -r
        g = abs(r.content())
        if g > 1:
            r = r.exquo(sp.Poly(g, t, domain='ZZ'))
        chain.append(r)
    def sign_changes(seq):
        nz = [s for s in seq if s != 0]
        return sum(1 for u, v in zip(nz, nz[1:]) if u * v < 0)
    v_pos = sign_changes([1 if poly.LC() > 0 else -1 for poly in chain])
    v_neg = sign_changes([(1 if poly.LC() > 0 else -1) * ((-1) ** poly.degree()) for poly in chain])
    ci = v_neg - v_pos
    if da == db:
        if ci == da: return True, "B<A", ci
        if ci == -da: return True, "A<B", ci
        return False, f"FAIL({ci}/{da})", ci
    elif da == db + 1:
        if ci == da: return True, "A_outer", ci
        return False, f"FAIL({ci}/{da})", ci
    elif db == da + 1:
        ok, d, c2 = exact_cauchy_interlace(B, A)
        return ok, ("B_outer" if d == "A_outer" else d), -c2
    return False, f"deg({da},{db})", ci

t0 = time.time()
KMAX = 15
p = build_p(KMAX)

print(f"=== TEST A: Pair (2,1): Dtilde_{{k+1}} p_{{k,i}} << p_{{k,i-2}} (k=2..{KMAX}) ===", flush=True)
fails_21 = []
for k in range(2, KMAX + 1):
    for i in range(2, k + 1):
        A = p[k][i - 2]
        B = apply_Dtilde(k + 1, p[k][i])
        ok, d, ci = exact_cauchy_interlace(A, B)
        if not ok or d != "A_outer":
            fails_21.append((k, i, d))
print(f"  fails ({time.time()-t0:.1f}s): {fails_21[:12]} (total {len(fails_21)})", flush=True)

t1 = time.time()
print(f"\n=== TEST B: X < Y: p_{{k-1,i-2}} < Dtilde_k^2 p_{{k-1,i}} (k=3..{KMAX}) ===", flush=True)
fails_XY = []
for k in range(3, KMAX + 1):
    for i in range(2, k):
        X = p[k - 1][i - 2]
        Y = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][i]))
        ok, d, ci = exact_cauchy_interlace(Y, X)
        if not ok or d != "B<A":
            fails_XY.append((k, i, d))
print(f"  fails ({time.time()-t1:.1f}s): {fails_XY[:12]} (total {len(fails_XY)})", flush=True)

t2 = time.time()
print(f"\n=== TEST B2: EVEN STRONGER: p_{{k-1,i-2}} < Dtilde_{{k-1}}^2 p_{{k-1,i}}? (k=3..12) ===", flush=True)
fails_XY2 = []
for k in range(3, 13):
    for i in range(2, k):
        X = p[k - 1][i - 2]
        Y = apply_Dtilde(k - 1, apply_Dtilde(k - 1, p[k - 1][i]))
        ok, d, ci = exact_cauchy_interlace(Y, X)
        if not ok or d != "B<A":
            fails_XY2.append((k, i, d))
print(f"  fails ({time.time()-t2:.1f}s): {fails_XY2[:12]} (total {len(fails_XY2)})", flush=True)

t3 = time.time()
print(f"\n=== TEST C: Turan M^2 - X*Y > 0 on (-1,0) & ALL R (k=3..10) ===", flush=True)
fails_tm = []; fails_tR = []
for k in range(3, 11):
    for i in range(2, k):
        X = p[k - 1][i - 2]
        M = apply_Dtilde(k, p[k - 1][i - 1])
        Y = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][i]))
        Dpoly = poly_sub(poly_mul(M, M), poly_mul(X, Y))
        P = sp.Poly(list(reversed(Dpoly)), t, domain='ZZ')
        c_m10 = sp.count_roots(P, -1, 0)
        c_all = sp.count_roots(P)
        val_pos = sum(c * ((-1)**m) * (2**(len(Dpoly)-1-m)) for m, c in enumerate(Dpoly)) > 0
        if c_m10 > 0 or not val_pos: fails_tm.append((k, i, c_m10, val_pos))
        if c_all > 0 or not val_pos: fails_tR.append((k, i, c_all, val_pos))
print(f"  fails (-1,0): {fails_tm[:10]}, fails R: {fails_tR[:10]} ({time.time()-t3:.1f}s)", flush=True)

t4 = time.time()
print(f"\n=== TEST D: Row log-concavity p_{{k,i}}^2 - p_{{k,i-1}}*p_{{k,i+1}} (k=2..10) ===", flush=True)
fails_rm = []; fails_rR = []
for k in range(2, 11):
    for i in range(1, k):
        Dpoly = poly_sub(poly_mul(p[k][i], p[k][i]), poly_mul(p[k][i-1], p[k][i+1]))
        P = sp.Poly(list(reversed(Dpoly)), t, domain='ZZ')
        c_m10 = sp.count_roots(P, -1, 0)
        c_all = sp.count_roots(P)
        val_pos = sum(c * ((-1)**m) * (2**(len(Dpoly)-1-m)) for m, c in enumerate(Dpoly)) > 0
        if c_m10 > 0 or not val_pos: fails_rm.append((k, i, c_m10, val_pos))
        if c_all > 0 or not val_pos: fails_rR.append((k, i, c_all, val_pos))
print(f"  fails (-1,0): {fails_rm[:10]}, fails R: {fails_rR[:10]} ({time.time()-t4:.1f}s)", flush=True)
