"""Deep exact Cauchy-index check of:
1. A_{k,i} << A_{k,i-1} (same k)
2. Pair (2,1): S_{3k-3-i} A_{k-1,i} << A_{k-1,i-2}
3. X < Y: A_{k-1,i-2} < S_{3k-2-i} S_{3k-3-i} A_{k-1,i}
up to k = 25 in the x-variable (using fast primitive integer PRS)!
"""
import time
import sympy as sp

x = sp.Symbol('x')

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_S(L, p):
    d = len(p) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(p):
        res[m] += c * (m + 1)
        res[m + 1] += c * (L - m)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def build_A(kmax):
    A = [[None] * (k + 1) for k in range(kmax + 1)]
    A[0][0] = [1]
    for k in range(1, kmax + 1):
        A[k][0] = apply_S(3*k - 2, apply_S(3*k - 3, A[k-1][0]))
        for i in range(1, k):
            L = 3*k - 2 - i
            A[k][i] = apply_S(L, poly_add(apply_S(L - 1, A[k-1][i]), A[k-1][i-1]))
        A[k][k] = apply_S(2*k - 2, A[k-1][k-1])
    return A

def exact_cauchy_interlace(A, B):
    p0 = sp.Poly(list(reversed(A)), x, domain='ZZ')
    p1 = sp.Poly(list(reversed(B)), x, domain='ZZ')
    da, db = p0.degree(), p1.degree()
    if da == 0 or db == 0:
        return True, "trivial"
    if p0.LC() < 0: p0 = -p0
    if p1.LC() < 0: p1 = -p1
    chain = [p0, p1]
    while chain[-1].degree() > 0:
        a, b = chain[-2], chain[-1]
        mdeg = a.degree() - b.degree() + 1
        r = sp.prem(a, b)
        if r.is_zero:
            return False, "common_root"
        if b.LC() < 0 and (mdeg % 2 == 1):
            r = -r
        r = -r
        g = abs(r.content())
        if g > 1:
            r = r.exquo(sp.Poly(g, x, domain='ZZ'))
        chain.append(r)
    def sign_changes(seq):
        nz = [s for s in seq if s != 0]
        return sum(1 for u, v in zip(nz, nz[1:]) if u * v < 0)
    v_pos = sign_changes([1 if poly.LC() > 0 else -1 for poly in chain])
    v_neg = sign_changes([(1 if poly.LC() > 0 else -1) * ((-1) ** poly.degree()) for poly in chain])
    ci = v_neg - v_pos
    if da == db:
        if ci == da: return True, "B<A"
        if ci == -da: return True, "A<B"
        return False, f"FAIL({ci}/{da})"
    elif da == db + 1:
        if ci == da: return True, "A_outer"
        return False, f"FAIL({ci}/{da})"
    return False, f"deg({da},{db})"

KMAX = 22
A = build_A(KMAX)
t0 = time.time()
for k in range(16, KMAX + 1):
    # Check deep interior + full row!
    f_row = []; f_21 = []; f_xy = []
    for i in range(1, k + 1):
        ok, d = exact_cauchy_interlace(A[k][i-1], A[k][i])
        if not ok or d != "A_outer": f_row.append((i, d))
    for i in range(2, k):
        L = 3*k - 2 - i
        SA = apply_S(L - 1, A[k-1][i])
        ok1, d1 = exact_cauchy_interlace(A[k-1][i-2], SA)
        if not ok1 or d1 != "A_outer": f_21.append((i, d1))
        SSA = apply_S(L, SA)
        ok2, d2 = exact_cauchy_interlace(SSA, A[k-1][i-2])
        if not ok2 or d2 != "B<A": f_xy.append((i, d2))
    print(f"k={k:2d} ({time.time()-t0:5.1f}s): row_fails={f_row}, pair21_fails={f_21}, XY_fails={f_xy}", flush=True)
