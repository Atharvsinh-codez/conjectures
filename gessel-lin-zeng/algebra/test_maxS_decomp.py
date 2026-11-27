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

def exact_strict_interlaces(p, q):
    dp, dq = len(p) - 1, len(q) - 1
    if p[-1] <= 0 or q[-1] <= 0: return False, "nonpos_lc"
    if dq == 0: return dp == 0, "deg0"
    if dq == dp:
        lcp, lcq = p[-1], q[-1]
        r1 = [lcq * p[m] - lcp * q[m] for m in range(dq)]
        while len(r1) > 1 and r1[-1] == 0: r1.pop()
        if len(r1) - 1 != dq - 1 or r1[-1] <= 0: return False, "eq_lc_fail"
        r0 = q[:]
    elif dq == dp + 1:
        r0, r1 = q[:], p[:]
    else:
        return False, f"deg_jump({dp}->{dq})"
    cur0, cur1 = r0[:], r1[:]
    while len(cur1) > 1:
        d0, d1 = len(cur0) - 1, len(cur1) - 1
        if d0 != d1 + 1: return False, f"deg_drop({d0},{d1})"
        a0, a1 = cur0[-1], cur1[-1]
        tmp = [a1 * cur0[m] for m in range(d0 + 1)]
        for m in range(d1 + 1): tmp[m + 1] -= a0 * cur1[m]
        b = tmp[d1]
        nxt = [-(a1 * tmp[m] - b * cur1[m]) for m in range(d1)]
        while len(nxt) > 1 and nxt[-1] == 0: nxt.pop()
        if len(nxt) - 1 != d1 - 1 or nxt[-1] <= 0: return False, f"sturm_fail_deg{d1-1}"
        g = math.gcd(*[abs(v) for v in nxt])
        if g > 1: nxt = [v // g for v in nxt]
        cur0, cur1 = cur1, nxt
    return True, "ok"

# Compute htilde_{k, i} for k=0..16 (with htilde_{0,0} = [1] careful: h_{0,0}=1, h_{1,0}=(1+t)(1+2t), h_{1,1}=(1+t)*1)
# Using exact recurrence for k >= 2: ht_{k,i} = D_{k+1}(D_{k+1} ht_{k-1,i} + ht_{k-1,i-1})
# Note for m=1, i=1: M_{k, 1, 1} = (D_{k+1}^2 ... D_3^2) ht_{1,1} = (D_{k+1}^2 ... D_3^2)(1).
# For m >= 2, i >= 1: M_{k, i, m} = (D_{k+1}^2 ... D_{m+2}^2 D_{m+1}) ht_{m-1, i-1}!
# Wait: check at k=2, i=1:
#   ht_{2,1} = D_3^2 ht_{1,1} + D_3 ht_{1,0} = M_{2,1,1} + M_{2,1,2}!
#   Matches: M_{2,1,1} = D_3^2 ht_{1,1}, M_{2,1,2} = D_3 ht_{1,0}!
KMAX = 15
ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, KMAX + 1):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

# Build M_{k, i, m} for 1 <= i <= k, i <= m <= k!
M = {}
for k in range(1, KMAX + 1):
    for i in range(1, k + 1):
        # m = k:
        if k == 1:
            M[(1, 1, 1)] = [1]
        else:
            M[(k, i, k)] = apply_D(k + 1, ht[(k-1, i-1)])
        # m < k: M_{k, i, m} = D_{k+1}^2 M_{k-1, i, m}!
        for m in range(i, k):
            M[(k, i, m)] = apply_D(k + 1, apply_D(k + 1, M[(k-1, i, m)]))
        # Sanity check sum_m M_{k, i, m} == ht_{k, i}!
        s = [0]
        for m in range(i, k + 1):
            s = add_poly(s, M[(k, i, m)])
        assert s == ht[(k, i)], f"Mismatch at ({k},{i})"

print("=== TEST 1: Check M_{k, i, m} (grouped by max(S) = m) interlacing for k=2..15 (incl. (13,7) & (14,8)!) ===")
for k in range(3, KMAX + 1):
    adj_fails = []
    strong_fails = []
    for i in range(1, k):
        # Check M_{k,i,m} << M_{k,i,m-1} (Wait: later deletion m -> roots shifted LEFT! So M_{k,i,k} << M_{k,i,k-1} << ... << M_{k,i,i}!)
        for m in range(i + 1, k + 1):
            ok, msg = exact_strict_interlaces(M[(k, i, m)], M[(k, i, m - 1)])
            if not ok:
                adj_fails.append((i, m, msg))
        # Check extreme pair M_{k, i, k} << M_{k, i, i} (strong interlacing!)
        for r in range(i, k + 1):
            for m in range(r + 1, k + 1):
                ok, msg = exact_strict_interlaces(M[(k, i, m)], M[(k, i, r)])
                if not ok:
                    strong_fails.append((i, m, r, msg))
                    break
            if strong_fails and strong_fails[-1][0] == i:
                break
    print(f"  k={k:2d}: adj (M_m << M_{{m-1}}) fails={len(adj_fails):2d} {adj_fails[:2]} | strong (all m>r) fails={len(strong_fails):2d} {strong_fails[:3]}")

# === TEST 2: Compute EXACT mean separation mu(M_{k,i,i}) - mu(M_{k,i,k}) up to k=30! ===
# Note: in x, mu(A) = A'(1)/A(1). In t: x = t/(1+t) -> x=1 is t -> infinity!
# So A(1) = lc(ht), and A'(1) = deg*lc(ht) - sub_lc(ht)!
# Wait: check! ht(t) = (1+t)^d A(t/(1+t)). Let u = 1/t -> t/(1+t) = 1/(1+u) = 1 - u + u^2 - ...
# u^d ht(1/u) = (1+u)^d A(1 - u + ...) = (1 + d*u + ...)(A(1) - A'(1)*u + ...) = A(1) + (d*A(1) - A'(1))*u + ...
# So [t^{d-1}] ht / [t^d] ht = d - A'(1)/A(1) = d - mu(A)!!
# Therefore mu(A) = d - sub_lc(ht) / lc(ht)!!
print("\n=== TEST 2: Mean separation mu(M_{k,i,i}) - mu(M_{k,i,k}) for k up to 40 (exact rational/float!) ===")
# Track (lc, sub_lc, deg) for ht_{k, i} and M_{k, i, i}, M_{k, i, k}!
# Under D_j: p = c_d t^d + c_{d-1} t^{d-1} + ...
# D_j(p) = (d+j) c_d t^{d+1} + ((d+1) c_d + (d-1+j) c_{d-1}) t^d + ...
# Wait: check D_j(t^d) = (d+1) t^d + (d+j) t^{d+1}. YES!
# Let mu = d - c_{d-1}/c_d. Under D_j:
# new_deg = d + 1.
# new_c_{d+1} = (d+j) c_d.
# new_c_d = (d+1) c_d + (d-1+j) c_{d-1}.
# new_mu = (d+1) - [(d+1) + (d-1+j)(d - mu)] / (d+j) = ...
# Let's just track exact rational (lc, sub_lc) scaled by lc! Or float since normalized!
def apply_D_stat(j, d, mu):
    # Returns (d+1, new_mu)
    # c_{d-1}/c_d = d - mu
    ratio = ((d + 1) + (d - 1 + j) * (d - mu)) / (d + j)
    return d + 1, (d + 1) - ratio

# Actually track (log_lc, mu) or exact Fraction? Float is exact to 1e-15 when we only track ratio lc1/lc2!
# Wait: for M_{k, i, m} we need lc too for sums! Use Python ints for (lc, sub_lc) — only 2 ints per (k, i)!
lc_h = {(1, 0): 2, (1, 1): 1}
sub_h = {(1, 0): 1, (1, 1): 0}
for k in range(2, 41):
    for i in range(k + 1):
        # D_{k+1} ht_{k-1, i}
        if i <= k - 1:
            d1 = 2*k - i - 3
            a_lc = (d1 + k + 1) * lc_h[(k-1, i)]
            a_sub = (d1 + 1) * lc_h[(k-1, i)] + (d1 + k) * sub_h[(k-1, i)]
        else:
            a_lc, a_sub = 0, 0
        s_lc = a_lc + (lc_h[(k-1, i-1)] if i >= 1 else 0)
        s_sub = a_sub + (sub_h[(k-1, i-1)] if i >= 1 else 0)
        ds = 2*k - i - 2
        lc_h[(k, i)] = (ds + k + 1) * s_lc
        sub_h[(k, i)] = (ds + 1) * s_lc + (ds + k) * s_sub

# Compute mu(M_{k, i, i}) and mu(M_{k, i, k}) for k = 10, 20, 30, 40!
for k in [10, 15, 20, 30, 40]:
    max_sep = 0.0
    max_at = None
    for i in range(1, k):
        # M_{k, i, i} starts at ht_{i, i} (deg i-1) and applies D_{i+2}^2 ... D_{k+1}^2!
        # Wait: M_{i, i, i} = ht_{i, i} (check: at k=i, sum has only m=i, so M_{i,i,i} = ht_{i,i}!).
        d_lo = i - 1
        mu_lo = d_lo - sub_h[(i, i)] / lc_h[(i, i)]
        for r in range(i + 1, k + 1):
            d_lo, mu_lo = apply_D_stat(r + 1, d_lo, mu_lo)
            d_lo, mu_lo = apply_D_stat(r + 1, d_lo, mu_lo)
        # M_{k, i, k} = D_{k+1} ht_{k-1, i-1}!
        d_hi = 2*(k-1) - (i-1) - 1
        mu_hi = d_hi - sub_h[(k-1, i-1)] / lc_h[(k-1, i-1)]
        d_hi, mu_hi = apply_D_stat(k + 1, d_hi, mu_hi)
        sep = mu_lo - mu_hi
        if sep > max_sep:
            max_sep = sep
            max_at = (i, mu_lo, mu_hi)
    print(f"  k={k:2d}: max mean sep mu(M_{{k,i,i}}) - mu(M_{{k,i,k}}) = {max_sep:.6f} at i={max_at[0]} (mu_lo={max_at[1]:.3f}, mu_hi={max_at[2]:.3f})")
