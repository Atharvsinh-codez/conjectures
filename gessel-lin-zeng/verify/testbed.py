"""
Independent exact-arithmetic verification testbed for Gessel-Lin-Zeng Conjecture 15
(Jacobi-Stirling descent polynomials A_{k,i}(x) and transformed polynomials h_{k,i}(t)).

Zero floating-point handwaving: all polynomial arithmetic, Sturm counts, and Cauchy
index interlacing certificates use exact Python integers.
"""

from __future__ import annotations
from itertools import combinations
import math
import time
from typing import Dict, List, Tuple


# ---------------------------------------------------------------------------
# 1. Exact integer polynomial & Sturm / Cauchy-index verification library
# ---------------------------------------------------------------------------

def poly_trim(p: List[int]) -> List[int]:
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def poly_add(a: List[int], b: List[int]) -> List[int]:
    n = max(len(a), len(b))
    res = [0] * n
    for i, c in enumerate(a):
        res[i] += c
    for i, c in enumerate(b):
        res[i] += c
    return poly_trim(res)


def poly_sub(a: List[int], b: List[int]) -> List[int]:
    n = max(len(a), len(b))
    res = [0] * n
    for i, c in enumerate(a):
        res[i] += c
    for i, c in enumerate(b):
        res[i] -= c
    return poly_trim(res)


def poly_mul(a: List[int], b: List[int]) -> List[int]:
    if a == [0] or b == [0]:
        return [0]
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                res[i + j] += ai * bj
    return poly_trim(res)


def poly_deriv(p: List[int]) -> List[int]:
    if len(p) <= 1:
        return [0]
    return poly_trim([i * p[i] for i in range(1, len(p))])


def poly_div_exact_linear(p: List[int], root: int) -> List[int]:
    """Divide p(t) exactly by (t - root). Asserts zero remainder."""
    p = poly_trim(p)
    d = len(p) - 1
    assert d >= 1
    q = [0] * d
    rem = p[-1]
    q[d - 1] = rem
    for i in range(d - 2, -1, -1):
        rem = p[i + 1] + root * rem
        q[i] = rem
    final_rem = p[0] + root * rem
    assert final_rem == 0, f"Polynomial not divisible by (t - {root}), rem={final_rem}"
    return poly_trim(q)


def poly_primitive_pos(p: List[int]) -> List[int]:
    p = poly_trim(p)
    if p == [0]:
        return [0]
    g = 0
    for c in p:
        g = math.gcd(g, abs(c))
    if g > 1:
        return [c // g for c in p]
    return p


def poly_pseudo_rem_sign_preserving(a: List[int], b: List[int]) -> List[int]:
    """Return c * rem(a, b) with strictly positive scalar c > 0."""
    a = poly_trim(a)
    b = poly_trim(b)
    assert b != [0]
    da = len(a) - 1
    db = len(b) - 1
    if da < db:
        return poly_primitive_pos(a)
    r = list(a)
    lb = b[-1]
    abs_lb = abs(lb)
    sgn_lb = 1 if lb > 0 else -1
    for k in range(da, db - 1, -1):
        if k < len(r) and r[k] != 0:
            rk = r[k]
            factor = rk * sgn_lb
            for i in range(len(r)):
                r[i] *= abs_lb
            shift = k - db
            for j in range(db + 1):
                r[shift + j] -= factor * b[j]
            r = poly_primitive_pos(r)
    return poly_primitive_pos(r)


def poly_eval_int(p: List[int], x: int) -> int:
    acc = 0
    for c in reversed(p):
        acc = acc * x + c
    return acc


def sign_int(x: int) -> int:
    return 1 if x > 0 else (-1 if x < 0 else 0)


def sign_variations(signs: List[int]) -> int:
    nz = [s for s in signs if s != 0]
    return sum(1 for i in range(len(nz) - 1) if nz[i] * nz[i + 1] < 0)


def sturm_chain(p0: List[int], p1: List[int]) -> List[List[int]]:
    r0 = poly_primitive_pos(p0)
    r1 = poly_primitive_pos(p1)
    chain = [r0]
    if r1 == [0]:
        return chain
    chain.append(r1)
    while len(chain[-1]) > 1:
        rem = poly_pseudo_rem_sign_preserving(chain[-2], chain[-1])
        if rem == [0]:
            break
        chain.append([-c for c in rem])
    return chain


def chain_signs_at(chain: List[List[int]], pt: str | int) -> List[int]:
    if pt == "-inf":
        return [sign_int(p[-1]) if (len(p) - 1) % 2 == 0 else -sign_int(p[-1]) for p in chain]
    if pt == "+inf":
        return [sign_int(p[-1]) for p in chain]
    assert isinstance(pt, int)
    return [sign_int(poly_eval_int(p, pt)) for p in chain]


def verify_simple_roots_in_interval(p: List[int], lo: str | int, hi: str | int) -> Tuple[bool, str]:
    """
    Exact Sturm certificate that ALL deg(p) roots of p in Z[x] are real, simple,
    and lie strictly inside (lo, hi).
    """
    p = poly_trim(p)
    d = len(p) - 1
    if d == 0:
        return (p[0] != 0, f"const {p[0]}")
    if isinstance(lo, int) and poly_eval_int(p, lo) == 0:
        return (False, f"root at lo={lo}")
    if isinstance(hi, int) and poly_eval_int(p, hi) == 0:
        return (False, f"root at hi={hi}")
    chain = sturm_chain(p, poly_deriv(p))
    if len(chain[-1]) != 1 or chain[-1][0] == 0:
        return (False, f"not squarefree (gcd deg {len(chain[-1]) - 1})")
    v_lo = sign_variations(chain_signs_at(chain, lo))
    v_hi = sign_variations(chain_signs_at(chain, hi))
    cnt = v_lo - v_hi
    if cnt != d:
        return (False, f"deg={d}, roots in ({lo},{hi})={cnt}")
    return (True, f"deg={d}, all {d} simple roots in ({lo},{hi})")


def verify_strict_interlacing(p: List[int], q: List[int], lo: str | int = "-inf", hi: str | int = 0) -> Tuple[bool, str]:
    """
    Exact Cauchy-index certificate for strict directed interlacing P << Q on (lo, hi):
      lc(P) > 0, lc(Q) > 0, deg(Q) in {deg(P), deg(P) + 1},
      all roots of P and Q are real, simple, strictly inside (lo, hi), and
        - deg(Q) == deg(P) = m:     alpha_m < beta_m < ... < alpha_1 < beta_1
        - deg(Q) == deg(P) + 1 = m: beta_m < alpha_{m-1} < ... < alpha_1 < beta_1
      (i.e. rightmost root of Q is strictly to the right of rightmost root of P).

    Uses Cauchy index I_{-inf}^{+inf}(P/Q) = V(-inf) - V(+inf) on sturm_chain(Q, P)
    together with boundary checks ensuring all roots lie in (lo, hi).
    """
    p = poly_trim(p)
    q = poly_trim(q)
    dp = len(p) - 1
    dq = len(q) - 1
    if p[-1] <= 0 or q[-1] <= 0:
        return (False, "leading coefficients must be > 0")
    if dq not in (dp, dp + 1):
        return (False, f"degree mismatch: deg(P)={dp}, deg(Q)={dq}")
    if dq == 0:
        return (True, "trivial deg 0")
    if isinstance(lo, int) and (poly_eval_int(p, lo) == 0 or poly_eval_int(q, lo) == 0):
        return (False, f"boundary root at lo={lo}")
    if isinstance(hi, int) and (poly_eval_int(p, hi) == 0 or poly_eval_int(q, hi) == 0):
        return (False, f"boundary root at hi={hi}")
    chain = sturm_chain(q, p)
    if len(chain[-1]) != 1 or chain[-1][0] == 0:
        return (False, f"shared factor of degree {len(chain[-1]) - 1}")
    v_neg = sign_variations(chain_signs_at(chain, "-inf"))
    v_pos = sign_variations(chain_signs_at(chain, "+inf"))
    v_lo = sign_variations(chain_signs_at(chain, lo))
    v_hi = sign_variations(chain_signs_at(chain, hi))
    if v_neg - v_pos != dq:
        return (False, f"Cauchy index on R is {v_neg - v_pos}, expected {dq}")
    if v_lo - v_hi != dq:
        return (False, f"Cauchy index on ({lo},{hi}) is {v_lo - v_hi}, expected {dq}")
    # Also ensure P has all its roots inside (lo, hi) (in equal-degree case, P's leftmost root
    # is < beta_m; we must verify it's still > lo when lo is finite).
    if isinstance(lo, int):
        ok_p, msg_p = verify_simple_roots_in_interval(p, lo, hi)
        if not ok_p:
            return (False, f"P root escaped ({lo},{hi}): {msg_p}")
    return (True, f"P(deg {dp}) << Q(deg {dq}) certified on ({lo},{hi})")


# ---------------------------------------------------------------------------
# 2. Direct brute-force enumeration of Jacobi-Stirling permutations A_{k,S}(x)
# ---------------------------------------------------------------------------

def enum_jacobi_stirling_word_poly(k: int, deleted_S: Tuple[int, ...]) -> List[int]:
    r"""
    Directly enumerate all Jacobi-Stirling permutations on
      M_{k,S} = {bar 1, 1, 1, ..., bar k, k, k} \ {bar j : j in S},
    with linear order bar 1 < 1 < bar 2 < 2 < ... < bar k < k,
    encoded as integers: bar j -> 2*j - 1, unbarred j -> 2*j.

    Condition: for every j in 1..k, every entry strictly between the two copies
    of (2*j) must be strictly greater than (2*j).
    Equivalent stack/nesting condition while scanning left-to-right:
      if copy #1 of (2*j) has been placed and copy #2 has NOT yet been placed
      (so 2*j is currently "open"), then any symbol placed now must be >= 2*j
      (in particular, the currently open unbarred values form a strictly increasing
      stack, and the next symbol placed must be >= top of the open stack!).

    Proof of equivalence:
      Suppose a value v is placed between the two copies of 2*j.
      Then 2*j is open when v is placed, so v > 2*j (since v cannot be 2*j unless it closes 2*j).
      Conversely, if every symbol placed while 2*j is open satisfies v >= 2*j (with v == 2*j closing it),
      then every symbol strictly between the two copies of 2*j is > 2*j.
      Note: since any open 2*j requires subsequent symbols >= 2*j, the max open value is simply
      the top of the open stack!

    Statistic: internal descents #{1 <= r <= L-1 : w_r > w_{r+1}}.
    """
    S_set = set(deleted_S)
    # Build available symbols:
    # bars: 2*j - 1 with count 1 (for j not in S)
    # unbarred: 2*j with count 2 (state: 0=unplaced, 1=open, 2=closed)
    rem_bars = [2 * j - 1 for j in range(1, k + 1) if j not in S_set]
    L = 3 * k - len(S_set)
    counts: Dict[int, int] = {}

    bar_used = [False] * len(rem_bars)
    pair_state = [0] * (k + 1)  # 1..k: 0=0 placed, 1=1 placed (open), 2=2 placed (closed)
    open_stack: List[int] = []  # stores unbarred values 2*j currently open (strictly increasing)

    def dfs(pos: int, last_val: int, des: int) -> None:
        if pos == L:
            counts[des] = counts.get(des, 0) + 1
            return
        min_allowed = open_stack[-1] if open_stack else 0
        # Option 1: close the top of open_stack (value == min_allowed)
        if open_stack:
            v = open_stack.pop()
            j = v // 2
            pair_state[j] = 2
            add_des = 1 if (pos > 0 and last_val > v) else 0
            dfs(pos + 1, v, des + add_des)
            pair_state[j] = 1
            open_stack.append(v)
        # Option 2: place an unused bar > min_allowed
        for idx, bv in enumerate(rem_bars):
            if not bar_used[idx] and bv > min_allowed:
                bar_used[idx] = True
                add_des = 1 if (pos > 0 and last_val > bv) else 0
                dfs(pos + 1, bv, des + add_des)
                bar_used[idx] = False
        # Option 3: open an unplaced unbarred pair (2*j) > min_allowed
        start_j = (min_allowed // 2) + 1
        for j in range(start_j, k + 1):
            if pair_state[j] == 0:
                v = 2 * j
                pair_state[j] = 1
                open_stack.append(v)
                add_des = 1 if (pos > 0 and last_val > v) else 0
                dfs(pos + 1, v, des + add_des)
                open_stack.pop()
                pair_state[j] = 0

    dfs(0, -1, 0)
    max_d = max(counts) if counts else 0
    return poly_trim([counts.get(d, 0) for d in range(max_d + 1)])


# ---------------------------------------------------------------------------
# 3. Insertion operators S_L and transformed operator D_j
# ---------------------------------------------------------------------------

def op_S(L: int, p: List[int]) -> List[int]:
    """S_L p(x) = (1 + L x) p(x) + x(1 - x) p'(x)."""
    dp = poly_deriv(p)
    term1 = poly_mul([1, L], p)
    term2 = poly_mul([0, 1, -1], dp)
    return poly_add(term1, term2)


def A_k_S_via_operator(k: int, deleted_S: Tuple[int, ...]) -> List[int]:
    """Compute A_{k,S}(x) by sequential insertion of symbols j=1..k."""
    S_set = set(deleted_S)
    p = [1]
    L = 0
    for j in range(1, k + 1):
        if j not in S_set:
            p = op_S(L, p)
            L += 1
        p = op_S(L, p)
        L += 2
    return p


def op_D(j: int, f: List[int]) -> List[int]:
    """D_j f(t) = t(1 + t) f'(t) + (1 + j t) f(t)."""
    df = poly_deriv(f)
    term1 = poly_mul([0, 1, 1], df)
    term2 = poly_mul([1, j], f)
    return poly_add(term1, term2)


def compute_A_table(max_k: int) -> Dict[Tuple[int, int], List[int]]:
    A: Dict[Tuple[int, int], List[int]] = {(0, 0): [1]}
    for k in range(1, max_k + 1):
        for i in range(k + 1):
            inner = [0]
            if i <= k - 1:
                inner = poly_add(inner, op_S(3 * k - 3 - i, A[(k - 1, i)]))
            if i >= 1:
                inner = poly_add(inner, A[(k - 1, i - 1)])
            A[(k, i)] = op_S(3 * k - 2 - i, inner)
    return A


def compute_h_table(max_k: int) -> Dict[Tuple[int, int], List[int]]:
    h: Dict[Tuple[int, int], List[int]] = {(0, 0): [1]}
    for k in range(1, max_k + 1):
        for i in range(k + 1):
            inner = [0]
            if i <= k - 1:
                inner = poly_add(inner, op_D(k, h[(k - 1, i)]))
            if i >= 1:
                inner = poly_add(inner, h[(k - 1, i - 1)])
            h[(k, i)] = op_D(k, inner)
    return h


def transform_A_to_h(A_poly: List[int], m: int) -> List[int]:
    """Compute (1 + t)^m A(t / (1 + t)) = sum_r a_r t^r (1 + t)^{m - r}."""
    res = [0] * (m + 1)
    for r, ar in enumerate(A_poly):
        if ar == 0:
            continue
        assert r <= m
        exp = m - r
        for s in range(exp + 1):
            res[r + s] += ar * math.comb(exp, s)
    return poly_trim(res)


def main() -> None:
    t0 = time.perf_counter()
    print("=== 1. Direct Jacobi-Stirling permutation enumeration vs S_L operator ===")
    # Full enumeration for k=1,2,3 (all subsets) and selected subsets for k=4
    for k in range(1, 4):
        for i in range(k + 1):
            sum_bf = [0]
            for S in combinations(range(1, k + 1), i):
                p_bf = enum_jacobi_stirling_word_poly(k, S)
                p_op = A_k_S_via_operator(k, S)
                assert p_bf == p_op, f"Mismatch at k={k}, S={S}: {p_bf} != {p_op}"
                sum_bf = poly_add(sum_bf, p_bf)
            print(f"  k={k}, i={i}: A_{{{k},{i}}}(x) = {sum_bf} (all {math.comb(k, i)} subsets verified by brute-force!)")

    # k=4: verify ALL subsets with i >= 1 (and i=0 has 246400 words, let's check i=0..4 too!)
    t_k4 = time.perf_counter()
    for i in range(5):
        for S in combinations(range(1, 5), i):
            p_bf = enum_jacobi_stirling_word_poly(4, S)
            p_op = A_k_S_via_operator(4, S)
            assert p_bf == p_op, f"Mismatch at k=4, S={S}"
    print(f"  k=4: ALL 16 subsets S subseteq [4] (including 246,400 words at S=empty) verified identically in {time.perf_counter()-t_k4:.2f}s!")

    # k=5: verify all subsets for i in {3, 4, 5} (16 subsets, up to ~20k words each)
    t_k5 = time.perf_counter()
    k5_subsets = 0
    for i in (3, 4, 5):
        for S in combinations(range(1, 6), i):
            p_bf = enum_jacobi_stirling_word_poly(5, S)
            p_op = A_k_S_via_operator(5, S)
            assert p_bf == p_op, f"Mismatch at k=5, S={S}"
            k5_subsets += 1
    print(f"  k=5: all {k5_subsets} subsets for i in {{3,4,5}} verified identically in {time.perf_counter()-t_k5:.2f}s!")

    print("\n=== 2. Aggregate A_{k,i} recurrence vs transformed h_{k,i} = D_k(D_k h_{k-1,i} + h_{k-1,i-1}) ===")
    MAX_K = 24
    A = compute_A_table(MAX_K)
    h = compute_h_table(MAX_K)
    for k in range(1, MAX_K + 1):
        for i in range(k + 1):
            h_from_A = transform_A_to_h(A[(k, i)], 2 * k - i)
            assert h[(k, i)] == h_from_A, f"Transform mismatch at ({k},{i})"
            assert len(A[(k, i)]) - 1 == 2 * k - i - 1, f"Degree mismatch at ({k},{i})"
            htilde = poly_div_exact_linear(h[(k, i)], -1)
            assert poly_eval_int(htilde, -1) != 0
    print(f"  Exact identity h_{{k,i}}(t) = (1+t)^{{2k-i}} A_{{k,i}}(t/(1+t)) verified for all 1 <= k <= {MAX_K}, 0 <= i <= k.")

    print(f"\n=== 3. Exact Sturm real-rootedness & aggregate interlacing diagnostics (k <= {MAX_K}) ===")
    t1 = time.perf_counter()
    for k in range(1, MAX_K + 1):
        for i in range(k + 1):
            ok, msg = verify_simple_roots_in_interval(A[(k, i)], "-inf", 0)
            assert ok, f"A_{{{k},{i}}} NOT simple negative real-rooted: {msg}"
    print(f"  ALL A_{{k,i}}(x) for 1 <= k <= {MAX_K}, 0 <= i <= k certified simple negative real-rooted!")

    row_pass = 0
    row_fail = []
    sum_G_ll_F = 0
    sum_fail = []
    U_pass = 0
    U_fail = []
    # Also check vertical interlacing htilde_{k-1, i} << htilde_{k, i} and htilde_{k-1, i-1} << htilde_{k, i}
    vert_same_pass = 0
    vert_same_fail = []
    vert_diag_pass = 0
    vert_diag_fail = []

    htilde: Dict[Tuple[int, int], List[int]] = {}
    for k in range(1, MAX_K + 1):
        for i in range(k + 1):
            htilde[(k, i)] = poly_div_exact_linear(h[(k, i)], -1)

    for k in range(2, MAX_K + 1):
        # Inner sums U[(k,i)] = (D_k h_{k-1,i} + h_{k-1,i-1}) / (1+t)
        U: Dict[int, List[int]] = {}
        for i in range(k + 1):
            inner = [0]
            if i <= k - 1:
                inner = poly_add(inner, op_D(k, h[(k - 1, i)]))
            if i >= 1:
                inner = poly_add(inner, h[(k - 1, i - 1)])
            U[i] = poly_div_exact_linear(inner, -1)

        for i in range(1, k + 1):
            ok, _ = verify_strict_interlacing(htilde[(k, i)], htilde[(k, i - 1)], -1, 0)
            if ok:
                row_pass += 1
            else:
                row_fail.append((k, i))

            ok_u, _ = verify_strict_interlacing(U[i], U[i - 1], -1, 0)
            if ok_u:
                U_pass += 1
            else:
                U_fail.append((k, i))

        for i in range(1, k):
            F = poly_div_exact_linear(op_D(k, h[(k - 1, i)]), -1)
            G = htilde[(k - 1, i - 1)]
            ok_GF, _ = verify_strict_interlacing(G, F, -1, 0)
            if ok_GF:
                sum_G_ll_F += 1
            else:
                sum_fail.append((k, i))

        # Check diagonal vertical: htilde_{k-1, i-1} (deg 2k-i-2) vs htilde_{k, i} (deg 2k-i-1)
        for i in range(1, k + 1):
            ok_d, _ = verify_strict_interlacing(htilde[(k - 1, i - 1)], htilde[(k, i)], -1, 0)
            if ok_d:
                vert_diag_pass += 1
            else:
                vert_diag_fail.append((k, i))

    print(f"  Row interlacing htilde_{{k,i}} << htilde_{{k,i-1}} (k=2..{MAX_K}): PASS={row_pass}, FAIL={len(row_fail)} {row_fail[:10]}")
    print(f"  Inner sum interlacing htilde_{{k-1,i-1}} << D_{{k+1}} htilde_{{k-1,i}}: PASS={sum_G_ll_F}, FAIL={len(sum_fail)} {sum_fail[:10]}")
    print(f"  Inner sum row U_{{k,i}} << U_{{k,i-1}}: PASS={U_pass}, FAIL={len(U_fail)} {U_fail[:10]}")
    print(f"  Diagonal interlacing htilde_{{k-1,i-1}} << htilde_{{k,i}}: PASS={vert_diag_pass}, FAIL={len(vert_diag_fail)} {vert_diag_fail[:10]}")
    print(f"  Total elapsed: {time.perf_counter()-t0:.2f}s.")


if __name__ == "__main__":
    main()
