"""Exact Sturm + exact interlacing + high-precision root analysis for Gessel-Lin-Zeng h_{k,i}(t) and A_{k,i}(x)."""
import math
from fractions import Fraction
import sympy as sp
import mpmath as mp

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for idx, v in enumerate(a):
        res[idx] += v
    for idx, v in enumerate(b):
        res[idx] += v
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def apply_D(j, poly):
    """D_j = t(1+t) d/dt + 1 + j t on low-to-high coeff list."""
    d = len(poly) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(poly):
        res[m] += c * (m + 1)
        res[m + 1] += c * (m + j)
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def apply_Dtilde(j, poly):
    """Dtilde_j(g) = (1+t)^{-1} D_j((1+t) g) = t(1+t) g' + (1 + (j+1) t) g."""
    d = len(poly) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(poly):
        res[m] += c * (m + 1)
        res[m + 1] += c * (m + j + 1)
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def build_p(kmax):
    """p_{k,i}(t) = h_{k,i}(t) / (1+t) for k >= 1, 0 <= i <= k.
    Recurrence: p_{1,0} = [1, 2], p_{1,1} = [1].
    For k >= 2:
      p_{k,0} = Dtilde_k^2 p_{k-1,0}
      p_{k,i} = Dtilde_k( Dtilde_k p_{k-1,i} + p_{k-1,i-1} )
      p_{k,k} = Dtilde_k p_{k-1,k-1}
    """
    p = [[None] * (k + 1) for k in range(kmax + 1)]
    p[1][0] = [1, 2]
    p[1][1] = [1]
    for k in range(2, kmax + 1):
        p[k][0] = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][0]))
        for i in range(1, k):
            p[k][i] = apply_Dtilde(k, poly_add(apply_Dtilde(k, p[k - 1][i]), p[k - 1][i - 1]))
        p[k][k] = apply_Dtilde(k, p[k - 1][k - 1])
    return p

def exact_real_rooted_in_m1_0(coeffs):
    """Use sympy exact Sturm count to verify all roots are simple and in (-1, 0)."""
    t = sp.Symbol('t')
    d = len(coeffs) - 1
    if d == 0:
        return True
    poly = sp.Poly(list(reversed(coeffs)), t, domain='ZZ')
    # Check squarefree and count real roots in (-1, 0)
    if poly.eval(0) == 0 or poly.eval(-1) == 0:
        return False
    cnt = sp.count_roots(poly, -1, 0)
    # Also verify squarefree
    if sp.gcd(poly, poly.diff()).degree() > 0:
        return False
    return cnt == d

def exact_interlace_same_deg(A, B):
    """Exact Hermite-Biehler / sign-at-roots check for equal-degree polynomials A, B
    with all simple roots in (-1, 0) and positive leading coefficients (Wait: let's check sign at 0!).
    Wait: use Wronskian or Euclidean / Cauchy index!
    Even cleaner: A and B (same degree d, both real-rooted, same leading sign) strictly interlace
    IFF for EVERY real lambda (or equivalently, the Hermite-Biehler condition: A + i*B has all roots in one open half-plane!
    Wait: even simpler: compute Cauchy index of B/A on (-inf, inf) via exact Sturm chain (A, B, rem...)!
    Let (p_0, p_1, ..., p_r) be the Euclidean remainder sequence with p_0 = A, p_1 = B, p_{j+1} = -rem(p_{j-1}, p_j).
    Let V(-inf) - V(+inf) be the Cauchy index I(B/A).
    If A has d simple real roots:
      I(B/A) = sum_{r: A(r)=0} sgn(B(r)/A'(r)).
    Strict interlacing of equal-degree A, B (positive leading coeffs) holds IFF |I(B/A)| == d!
      - If I(B/A) == +d: at every root r_m of A, B(r_m)/A'(r_m) > 0.
        At largest root r_d of A, A'(r_d) > 0 => B(r_d) > 0, so B has its largest root < r_d: b_1 < a_1 < ... < b_d < a_d (B < A)!
      - If I(B/A) == -d: at every root r_m of A, B(r_m)/A'(r_m) < 0.
        At largest root r_d of A, A'(r_d) > 0 => B(r_d) < 0, so B has a root > r_d: a_1 < b_1 < ... < a_d < b_d (A < B)!
    Let's double check for deg A = d, deg B = d-1 (positive leading coeffs):
      - If I(B/A) == +d: at largest root r_d, B(r_d) > 0, at r_{d-1}, B(r_{d-1}) < 0, ..., so B's d-1 roots strictly interlace A's d roots (A_outer)!
    HOLD ON! That's 100% EXACT using sympy's exact polynomial remainder sequence over QQ!
    """
    t = sp.Symbol('t')
    p0 = sp.Poly(list(reversed(A)), t, domain='QQ')
    p1 = sp.Poly(list(reversed(B)), t, domain='QQ')
    if p0.degree() == 0 or p1.degree() == 0:
        return True, "trivial", 0
    # Ensure positive leading coeff
    if p0.LC() < 0:
        p0 = -p0
    if p1.LC() < 0:
        p1 = -p1
    chain = [p0, p1]
    while chain[-1].degree() > 0:
        _, r = sp.div(chain[-2], chain[-1], domain='QQ')
        if r.is_zero:
            return False, "common_root", 0
        # Normalize by positive denominator/content to keep numbers small
        r = -r
        # divide by abs(LC)
        r = r.exquo(sp.Poly(abs(r.LC()), t, domain='QQ'))
        chain.append(r)
    def sign_changes(seq):
        nz = [s for s in seq if s != 0]
        return sum(1 for a, b in zip(nz, nz[1:]) if a * b < 0)
    v_pos = sign_changes([1 if poly.LC() > 0 else -1 for poly in chain])
    v_neg = sign_changes([ (1 if poly.LC() > 0 else -1) * ((-1) ** poly.degree()) for poly in chain])
    cauchy_idx = v_neg - v_pos
    da, db = p0.degree(), p1.degree()
    if da == db:
        if cauchy_idx == da:
            return True, "B<A", cauchy_idx
        elif cauchy_idx == -da:
            return True, "A<B", cauchy_idx
        else:
            return False, f"FAIL(I={cauchy_idx}/{da})", cauchy_idx
    elif da == db + 1:
        if cauchy_idx == da:
            return True, "A_outer", cauchy_idx
        else:
            return False, f"FAIL(I={cauchy_idx}/{da})", cauchy_idx
    elif db == da + 1:
        # swap
        ok, d, ci = exact_interlace_same_deg(B, A)
        return ok, ("B_outer" if d == "A_outer" else d), -ci
    return False, f"deg({da},{db})", cauchy_idx

KMAX = 16
p = build_p(KMAX)

print(f"=== EXACT STURM REAL-ROOTEDNESS CHECK FOR k=1..{KMAX} ===")
rr_fails = []
for k in range(1, KMAX + 1):
    for i in range(k + 1):
        if not exact_real_rooted_in_m1_0(p[k][i]):
            rr_fails.append((k, i))
print(f"Real-rootedness fails in k=1..{KMAX}: {rr_fails}")

print(f"\n=== EXACT CAUCHY-INDEX INTERLACING CHECKS FOR k=1..{KMAX} ===")
# 1. p_{k, i-1} (deg 2k-i) vs p_{k, i} (deg 2k-i-1): does p_{k,i-1} outer-interlace p_{k,i}?
dirs_row = set(); fails_row = []
for k in range(2, KMAX + 1):
    for i in range(1, k + 1):
        ok, d, ci = exact_interlace_same_deg(p[k][i-1], p[k][i])
        dirs_row.add(d)
        if not ok:
            fails_row.append((k, i, d))
print(f"1. p_{{k,i-1}} vs p_{{k,i}} (same k): dirs={dirs_row}, fails={fails_row[:10]}")

# 2. Dtilde_k(p_{k-1, i}) vs p_{k-1, i-1} (THE EXACT TWO SUMMANDS IN THE RECURRENCE!)
dirs_sum = set(); fails_sum = []
for k in range(2, KMAX + 1):
    for i in range(1, k):
        A = apply_Dtilde(k, p[k-1][i])
        B = p[k-1][i-1]
        ok, d, ci = exact_interlace_same_deg(A, B)
        dirs_sum.add(d)
        if not ok:
            fails_sum.append((k, i, d))
print(f"2. Dtilde_k p_{{k-1,i}} vs p_{{k-1,i-1}} (RECURRENCE SUMMANDS!): dirs={dirs_sum}, fails={fails_sum[:10]}")

# 3. Dtilde_k(p_{k-1, i}) vs p_{k-1, i-1} WITH SHIFTED INDEX: what about Dtilde_{k+1} p_{k,i} vs p_{k,i-1}?
dirs_next = set(); fails_next = []
for k in range(1, KMAX + 1):
    for i in range(1, k + 1):
        A = apply_Dtilde(k + 1, p[k][i])
        B = p[k][i-1]
        ok, d, ci = exact_interlace_same_deg(A, B)
        dirs_next.add(d)
        if not ok:
            fails_next.append((k, i, d))
print(f"3. Dtilde_{{k+1}} p_{{k,i}} vs p_{{k,i-1}}: dirs={dirs_next}, fails={fails_next[:10]}")

# 4. What about Dtilde_k(p_{k,i}) vs p_{k,i-1} (same k in Dtilde)?
dirs_sameD = set(); fails_sameD = []
for k in range(1, KMAX + 1):
    for i in range(1, k + 1):
        A = apply_Dtilde(k, p[k][i])
        B = p[k][i-1]
        ok, d, ci = exact_interlace_same_deg(A, B)
        dirs_sameD.add(d)
        if not ok:
            fails_sameD.append((k, i, d))
print(f"4. Dtilde_k p_{{k,i}} vs p_{{k,i-1}}: dirs={dirs_sameD}, fails={fails_sameD[:10]}")
