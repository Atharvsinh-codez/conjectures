"""FAST + EXACT test of ALL candidate closure invariants using exact Cauchy index & exact Sturm (NO polyroots convergence issues!)."""
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
    """Exact Cauchy index I(B/A) via integer primitive pseudo-remainder sequence (FAST, no QQ blowup!)."""
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
        # pseudo-remainder: prem(a, b) = (LC(b))^{deg(a)-deg(b)+1} * a mod b
        # Sign of multiplier (LC(b))^{deg(a)-deg(b)+1}:
        mdeg = a.degree() - b.degree() + 1
        r = sp.prem(a, b)
        if r.is_zero:
            return False, "common_root", 0
        # If LC(b) < 0 and mdeg is odd, prem flipped sign! Compensate so prem has same sign as true rem!
        if b.LC() < 0 and (mdeg % 2 == 1):
            r = -r
        # Sturm step: next = -rem
        r = -r
        # Remove positive integer content!
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

def count_roots_m1_0(coeffs):
    P = sp.Poly(list(reversed(coeffs)), t, domain='ZZ')
    return sp.count_roots(P, -1, 0), sp.count_roots(P)

KMAX = 14
p = build_p(KMAX)

print(f"=== TEST A: Pair (2,1): Dtilde_{{k+1}} p_{{k,i}} << p_{{k,i-2}} (k=2..{KMAX}) ===")
fails_21 = []
for k in range(2, KMAX + 1):
    for i in range(2, k + 1):
        A = p[k][i - 2]
        B = apply_Dtilde(k + 1, p[k][i])
        ok, d, ci = exact_cauchy_interlace(A, B)
        if not ok or d != "A_outer":
            fails_21.append((k, i, d))
print(f"  fails: {fails_21[:12]} (total {len(fails_21)})")

print(f"\n=== TEST B: X < Y: p_{{k-1,i-2}} < Dtilde_k^2 p_{{k-1,i}} (k=3..{KMAX}) ===")
fails_XY = []
for k in range(3, KMAX + 1):
    for i in range(2, k):
        X = p[k - 1][i - 2]
        Y = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][i]))
        ok, d, ci = exact_cauchy_interlace(Y, X) # check X < Y => B<A with A=Y, B=X
        if not ok or d != "B<A":
            fails_XY.append((k, i, d))
print(f"  fails: {fails_XY[:12]} (total {len(fails_XY)})")

print(f"\n=== TEST C: Turan M^2 - X*Y > 0 on (-1,0) AND ALL R (k=3..{KMAX}) ===")
fails_turan_m10 = []
fails_turan_R = []
for k in range(3, KMAX + 1):
    for i in range(2, k):
        X = p[k - 1][i - 2]
        M = apply_Dtilde(k, p[k - 1][i - 1])
        Y = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][i]))
        Dpoly = poly_sub(poly_mul(M, M), poly_mul(X, Y))
        c_m10, c_all = count_roots_m1_0(Dpoly)
        # evaluate at t = -1/2 (in rational: sum c_m (-1/2)^m)
        val_pos = sum(c * ((-1)**m) * (2**(len(Dpoly)-1-m)) for m, c in enumerate(Dpoly)) > 0
        if c_m10 > 0 or not val_pos:
            fails_turan_m10.append((k, i, c_m10, val_pos))
        if c_all > 0 or not val_pos:
            fails_turan_R.append((k, i, c_all, val_pos))
print(f"  M^2 - X*Y fails on (-1,0): {fails_turan_m10[:10]} (total {len(fails_turan_m10)})")
print(f"  M^2 - X*Y fails on ALL R:  {fails_turan_R[:10]} (total {len(fails_turan_R)})")

print(f"\n=== TEST D: Row log-concavity p_{{k,i}}^2 - p_{{k,i-1}}*p_{{k,i+1}} on (-1,0) & ALL R (k=2..{KMAX}) ===")
fails_rlc_m10 = []
fails_rlc_R = []
for k in range(2, KMAX + 1):
    for i in range(1, k):
        Dpoly = poly_sub(poly_mul(p[k][i], p[k][i]), poly_mul(p[k][i-1], p[k][i+1]))
        c_m10, c_all = count_roots_m1_0(Dpoly)
        val_pos = sum(c * ((-1)**m) * (2**(len(Dpoly)-1-m)) for m, c in enumerate(Dpoly)) > 0
        if c_m10 > 0 or not val_pos:
            fails_rlc_m10.append((k, i, c_m10, val_pos))
        if c_all > 0 or not val_pos:
            fails_rlc_R.append((k, i, c_all, val_pos))
print(f"  p_{{k,i}}^2 - p_{{k,i-1}}p_{{k,i+1}} fails on (-1,0): {fails_rlc_m10[:10]} (total {len(fails_rlc_m10)})")
print(f"  p_{{k,i}}^2 - p_{{k,i-1}}p_{{k,i+1}} fails on ALL R:  {fails_rlc_R[:10]} (total {len(fails_rlc_R)})")

print(f"\n=== TEST E: Check p_{{k,i-1}} < Dtilde_k p_{{k,i}} up to k={KMAX} ===")
fails_dk = []
for k in range(2, KMAX + 1):
    for i in range(1, k + 1):
        ok, d, ci = exact_cauchy_interlace(apply_Dtilde(k, p[k][i]), p[k][i - 1])
        if not ok or d != "B<A":
            fails_dk.append((k, i, d))
print(f"  p_{{k,i-1}} < Dtilde_k p_{{k,i}} fails (k=2..{KMAX}): {fails_dk[:10]}")
