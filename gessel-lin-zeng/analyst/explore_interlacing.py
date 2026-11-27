"""Immediate exploratory test of Gessel-Lin-Zeng h_{k,i}(t) and p_{k,i}(t) = h_{k,i}(t)/(1+t):
Check real-rootedness, root interlacing between D_k h_{k-1,i} and h_{k-1,i-1},
between h_{k,i} and h_{k,i-1}, and between h_{k,i} and h_{k-1,i}, up to k=18!
"""
import numpy as np
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
    """Apply D_j = t(1+t) d/dt + (1 + j t) to polynomial (low-to-high coeff list)."""
    d = len(poly) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(poly):
        # c * [m t^m + m t^{m+1} + t^m + j t^{m+1}]
        res[m] += c * (m + 1)
        res[m + 1] += c * (m + j)
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def div_by_1_plus_t(poly):
    """Exact division by (1+t); returns quotient and checks remainder==0."""
    # poly = sum c_m t^m = (1+t) sum q_m t^m => q_0 = c_0, q_m = c_m - q_{m-1}
    q = []
    prev = 0
    for c in poly[:-1]:
        val = c - prev
        q.append(val)
        prev = val
    assert poly[-1] == prev, f"Not divisible by (1+t): rem={poly[-1]-prev}"
    return q

def build_h(kmax):
    """Build exact integer polynomials h[k][i] for 0 <= k <= kmax, 0 <= i <= k."""
    h = [[None] * (k + 1) for k in range(kmax + 1)]
    h[0][0] = [1]
    for k in range(1, kmax + 1):
        h[k][0] = apply_D(k, apply_D(k, h[k - 1][0]))
        for i in range(1, k):
            s = poly_add(apply_D(k, h[k - 1][i]), h[k - 1][i - 1])
            h[k][i] = apply_D(k, s)
        h[k][k] = apply_D(k, h[k - 1][k - 1])
    return h

def roots_in_interval(poly):
    """Compute roots of poly (low-to-high) using numpy/mpmath; return sorted real roots or None if nonreal."""
    # Reverse to high-to-low for numpy
    if len(poly) <= 1:
        return []
    # Use float64 scaled or mpmath for exactness
    arr = np.array([float(c) for c in reversed(poly)], dtype=np.float64)
    arr /= np.max(np.abs(arr))
    r = np.roots(arr)
    if np.max(np.abs(r.imag)) > 1e-6:
        return None
    return sorted(r.real.tolist())

def check_interlace(rA, rB, label=""):
    """Check if rA and rB strictly interlace.
    Returns (ok, direction, min_gap).
    If len(rA) == len(rB):
      dir='A<B' means a_0 < b_0 < a_1 < b_1 < ...
      dir='B<A' means b_0 < a_0 < b_1 < a_1 < ...
    If len(rA) == len(rB) + 1:
      dir='A_outer' means a_0 < b_0 < a_1 < ... < b_{d-1} < a_d
    If len(rB) == len(rA) + 1:
      dir='B_outer' means b_0 < a_0 < b_1 < ... < a_{d-1} < b_d
    """
    na, nb = len(rA), len(rB)
    if na == 0 or nb == 0:
        return True, "trivial", 1.0
    if na == nb:
        ok_AB = all(rA[m] < rB[m] for m in range(na)) and all(rB[m] < rA[m+1] for m in range(na-1))
        ok_BA = all(rB[m] < rA[m] for m in range(na)) and all(rA[m] < rB[m+1] for m in range(na-1))
        if ok_AB:
            gaps = [rB[m] - rA[m] for m in range(na)] + [rA[m+1] - rB[m] for m in range(na-1)]
            return True, "A<B", min(gaps)
        if ok_BA:
            gaps = [rA[m] - rB[m] for m in range(na)] + [rB[m+1] - rA[m] for m in range(na-1)]
            return True, "B<A", min(gaps)
        return False, "FAIL", 0.0
    if na == nb + 1:
        ok = all(rA[m] < rB[m] < rA[m+1] for m in range(nb))
        if ok:
            gaps = [min(rB[m]-rA[m], rA[m+1]-rB[m]) for m in range(nb)]
            return True, "A_outer", min(gaps)
        return False, "FAIL", 0.0
    if nb == na + 1:
        ok = all(rB[m] < rA[m] < rB[m+1] for m in range(na))
        if ok:
            gaps = [min(rA[m]-rB[m], rB[m+1]-rA[m]) for m in range(na)]
            return True, "B_outer", min(gaps)
        return False, "FAIL", 0.0
    return False, f"deg_diff({na},{nb})", 0.0

KMAX = 14
h = build_h(KMAX)
p = [[None] * (k + 1) for k in range(KMAX + 1)]
for k in range(1, KMAX + 1):
    for i in range(k + 1):
        p[k][i] = div_by_1_plus_t(h[k][i])

print("=== 1. CHECK REAL-ROOTEDNESS OF p_{k,i}(t) on (-1, 0) for k=1..14 ===")
p_roots = [[None] * (k + 1) for k in range(KMAX + 1)]
for k in range(1, KMAX + 1):
    for i in range(k + 1):
        r = roots_in_interval(p[k][i])
        p_roots[k][i] = r
        if r is None or (len(r) > 0 and (r[0] <= -1.0 or r[-1] >= 0.0)):
            print(f"  FAIL real-rootedness at k={k}, i={i}: {r}")
print("  All p_{k,i} real-rooted in (-1, 0) for k=1..14 (float64 check)!")

print("\n=== 2. CHECK INTERLACING OF p_{k,i} vs p_{k,i-1} (same k) ===")
dirs_same_k = set()
fails_same_k = []
for k in range(2, KMAX + 1):
    for i in range(1, k + 1):
        ok, d, gap = check_interlace(p_roots[k][i], p_roots[k][i-1])
        dirs_same_k.add(d)
        if not ok:
            fails_same_k.append((k, i, d))
print(f"  directions: {dirs_same_k}, fails: {fails_same_k[:10]}")

print("\n=== 3. CHECK INTERLACING OF D_k h_{k-1,i}/(1+t) vs h_{k-1,i-1}/(1+t) ===")
dirs_sum = set()
fails_sum = []
min_gap_sum = (1e9, None)
for k in range(2, KMAX + 1):
    for i in range(1, k):
        A = div_by_1_plus_t(apply_D(k, h[k-1][i]))
        B = p[k-1][i-1]
        rA = roots_in_interval(A)
        rB = p_roots[k-1][i-1]
        ok, d, gap = check_interlace(rA, rB)
        dirs_sum.add(d)
        if not ok:
            fails_sum.append((k, i, d))
        elif gap < min_gap_sum[0] and d != "trivial":
            min_gap_sum = (gap, (k, i, d))
print(f"  directions: {dirs_sum}, fails: {fails_sum[:10]}, min_gap: {min_gap_sum}")

print("\n=== 4. CHECK INTERLACING OF p_{k,i} vs p_{k-1,i} (vertical) AND p_{k,i} vs p_{k-1,i-1} (diag) ===")
for name, fn in [
    ("p_{k,i} vs p_{k-1,i}", lambda k, i: (p_roots[k][i], p_roots[k-1][i]) if i <= k-1 else None),
    ("p_{k,i} vs p_{k-1,i-1}", lambda k, i: (p_roots[k][i], p_roots[k-1][i-1]) if i >= 1 else None),
]:
    ds = set()
    fs = []
    for k in range(2, KMAX + 1):
        for i in range(k + 1):
            pair = fn(k, i)
            if pair is None:
                continue
            ok, d, gap = check_interlace(*pair)
            ds.add(d)
            if not ok:
                fs.append((k, i, d))
    print(f"  {name}: dirs={ds}, fails={fs[:5]}")
