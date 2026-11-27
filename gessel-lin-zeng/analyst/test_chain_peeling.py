import numpy as np

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
    if len(p) <= 1:
        return np.array([], dtype=np.float64)
    c = np.array(p[::-1], dtype=np.float64)
    r = np.roots(c)
    if np.max(np.abs(r.imag)) > 1e-6:
        return None
    return np.sort(r.real)

def check_prec(f, g):
    """Check f << g (deg g = deg f + 1, roots of g surround roots of f: r1(g) < r1(f) < r2(g) < ... < r_d(g) < 0)."""
    rf = roots_real(f)
    rg = roots_real(g)
    if rf is None or rg is None:
        return False, "COMPLEX"
    if len(rg) != len(rf) + 1:
        return False, f"DEG({len(rf)},{len(rg)})"
    if len(rf) == 0:
        return (rg[0] < 0), "ok0"
    ok = all(rg[m] < rf[m] < rg[m+1] for m in range(len(rf)))
    if not ok:
        # measure worst violation
        v1 = min(rf[m] - rg[m] for m in range(len(rf)))
        v2 = min(rg[m+1] - rf[m] for m in range(len(rf)))
        return False, f"min_gap=({v1:.3e},{v2:.3e})"
    return True, "OK"

KMAX = 16
A = {(0, 0): [1]}
B = {}
for k in range(1, KMAX + 1):
    for i in range(k + 1):
        u1 = apply_S(A[(k-1, i)], 3*k - 3 - i) if i <= k - 1 else [0]
        u2 = A[(k-1, i-1)] if i >= 1 else [0]
        B[(k, i)] = poly_add(u1, u2)
        A[(k, i)] = apply_S(B[(k, i)], 3*k - 2 - i)

tests = [
    ("B_{k,i} << B_{k,i-1} (main target)", lambda k, i: (B[(k,i)], B[(k,i-1)]) if 1 <= i <= k else None),
    ("S_{3k-3-i} A_{k-1,i} << A_{k-1,i-2} (4th pair)", lambda k, i: (apply_S(A[(k-1,i)], 3*k-3-i), A[(k-1,i-2)]) if 2 <= i <= k-1 else None),
    ("A_{k-1,i} << B_{k-1,i-2} (peel 1)", lambda k, i: (A[(k-1,i)], B[(k-1,i-2)]) if 2 <= i <= k-1 else None),
    ("B_{k,i} << A_{k-1,i-2} (level k)", lambda k, i: (B[(k,i)], A[(k-1,i-2)]) if 2 <= i <= k else None),
    ("A_{k-1,i} << A_{k-2,i-3} (peel 2a)", lambda k, i: (A[(k-1,i)], A[(k-2,i-3)]) if 3 <= i <= k-1 and k >= 3 else None),
    ("B_{k,i} << B_{k-1,i-3} (peel 2b)", lambda k, i: (B[(k,i)], B[(k-1,i-3)]) if 3 <= i <= k and k >= 3 else None),
    ("S_{3k-3-i} A_{k-1,i} << A_{k-2,i-4} (peel 3!)", lambda k, i: (apply_S(A[(k-1,i)], 3*k-3-i), A[(k-2,i-4)]) if 4 <= i <= k-1 and k >= 4 else None),
    ("A_{k,i} << A_{k-1,i-3} (general jump!)", lambda k, i: (A[(k,i)], A[(k-1,i-3)]) if 3 <= i <= k and k >= 3 else None),
    ("A_{k,i} << A_{k-2,i-5} (jump 2!)", lambda k, i: (A[(k,i)], A[(k-2,i-5)]) if 5 <= i <= k and k >= 5 else None),
    ("A_{k,i} << A_{k-r, i-2r-1} for ALL r>=0!", None),
]

for name, fn in tests[:-1]:
    ok_cnt = 0; fail_cnt = 0; first_fail = None
    for k in range(2, KMAX + 1):
        for i in range(k + 1):
            pair = fn(k, i)
            if pair is None: continue
            ok, msg = check_prec(pair[0], pair[1])
            if ok: ok_cnt += 1
            else:
                fail_cnt += 1
                if first_fail is None: first_fail = (k, i, msg)
    print(f"  {name:48s}: pass={ok_cnt:4d}, fail={fail_cnt:4d}  first_fail={first_fail}")

# Check GENERAL jump family: A_{k, i} << A_{k-r, i - 2r - 1} AND B_{k, i} << B_{k-r, i - 2r - 1}
# Note deg A_{k,i} = 2k - i - 1; deg A_{k-r, i-2r-1} = 2(k-r) - (i-2r-1) - 1 = 2k - i = deg A_{k,i} + 1!
# AND check B_{k, i} << A_{k-r, i - 2r} (deg B_{k,i} = 2k-i-2, deg A_{k-r, i-2r} = 2k-i-1 = deg B + 1)!
# AND check A_{k, i} << B_{k-r, i - 2r - 2} (deg B_{k-r, i-2r-2} = 2(k-r)-(i-2r-2)-2 = 2k-i = deg A_{k,i} + 1)!
print("\n=== GENERAL r-JUMP FAMILIES (k <= 16) ===")
for r in range(0, 6):
    # Family 1: A_{k,i} << A_{k-r, i-2r-1}
    p1 = f1 = 0; ff1 = None
    # Family 2: B_{k,i} << B_{k-r, i-2r-1}
    p2 = f2 = 0; ff2 = None
    # Family 3: B_{k,i} << A_{k-r, i-2r}
    p3 = f3 = 0; ff3 = None
    # Family 4: A_{k,i} << B_{k-r, i-2r-2}
    p4 = f4 = 0; ff4 = None
    # Family 5: EQUAL-DEGREE A_{k-r, i-2r} left of A_{k, i}
    for k in range(r + 1, KMAX + 1):
        for i in range(k + 1):
            if i - 2*r - 1 >= 0 and i - 2*r - 1 <= k - r:
                ok, msg = check_prec(A[(k,i)], A[(k-r, i-2*r-1)])
                if ok: p1 += 1
                else:
                    f1 += 1
                    if ff1 is None: ff1 = (k, i, msg)
                if k - r >= 1:
                    ok, msg = check_prec(B[(k,i)], B[(k-r, i-2*r-1)])
                    if ok: p2 += 1
                    else:
                        f2 += 1
                        if ff2 is None: ff2 = (k, i, msg)
            if r >= 1 and i - 2*r >= 0 and i - 2*r <= k - r:
                ok, msg = check_prec(B[(k,i)], A[(k-r, i-2*r)])
                if ok: p3 += 1
                else:
                    f3 += 1
                    if ff3 is None: ff3 = (k, i, msg)
            if i - 2*r - 2 >= 0 and i - 2*r - 2 <= k - r and k - r >= 1:
                ok, msg = check_prec(A[(k,i)], B[(k-r, i-2*r-2)])
                if ok: p4 += 1
                else:
                    f4 += 1
                    if ff4 is None: ff4 = (k, i, msg)
    print(f"  r={r}: A<<A: {p1:3d}/{f1:<3d} {ff1} | B<<B: {p2:3d}/{f2:<3d} {ff2} | B<<A: {p3:3d}/{f3:<3d} {ff3} | A<<B: {p4:3d}/{f4:<3d} {ff4}")
