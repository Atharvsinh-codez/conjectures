"""
2-second test of the EXACT algebraic recurrence for P_{k, i, r}(x) (first_base refinement):
  For 1 <= r <= k-1:  P_{k, i, r} = S_{L-1}(S_{L-2} P_{k-1, i, r} + P_{k-1, i-1, r})
  For r = k:          P_{k, i, k} = S_L(x A_{k-1, i}) + x S_{L-2} A_{k-1, i} + x A_{k-1, i-1}
where L = 3k - 2 - i, and sum_{r=1}^k P_{k, i, r} = A_{k, i}!
Wait: ALSO check in the t-variable:
  P~_{k, i, r} = D_k(D_k P~_{k-1, i, r} + P~_{k-1, i-1, r})   (1 <= r <= k-1)
  P~_{k, i, k} = D_{k+1}(t h~_{k-1, i}) + t(D_k h~_{k-1, i} + h~_{k-1, i-1})
1. Verify against exact JSP_{5, 3} and JSP_{5, 4} word sums (just check sum_r = A_{k,i}!).
2. Test ALL-PAIRS interlacing of {P_{k, i, 1}, ..., P_{k, i, k}} for k = 2..12!
3. Also test whether P_{k, i, r} << P_{k, i-1, r} or what interlacing order they satisfy!
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 50

def apply_SL_mp(p, L):
    out = [mp.mpf(0)] * (len(p) + 1)
    for m, c in enumerate(p):
        if c != 0:
            out[m] += (1 + m) * c
            out[m + 1] += (L - m) * c
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out

def poly_add_mp(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out

def mul_x_mp(p):
    return [mp.mpf(0)] + list(p)

def roots_mp_fast(p):
    # Use float64 with rescaling for fast screening, mpmath only if needed
    # Wait: normalize coefficients so np.roots is well-conditioned!
    d = len(p) - 1
    if d <= 0:
        return np.array([])
    # strip trailing zeros (roots at x=0!)
    z0 = 0
    while z0 < len(p) and p[z0] == 0:
        z0 += 1
    q = p[z0:]
    dq = len(q) - 1
    if dq == 0:
        return np.zeros(z0)
    # geometric scale s = (q[0]/q[-1])^(1/dq)
    s = (abs(q[0]) / abs(q[-1])) ** (mp.mpf(1) / dq)
    sc = [float(q[m] * (s ** m)) for m in range(dq + 1)]
    mx = max(abs(v) for v in sc)
    sc = [v / mx for v in sc]
    r = np.roots(sc[::-1])
    if np.max(np.abs(r.imag)) > 1e-5:
        return None
    r_real = np.sort(r.real * float(s))
    if z0 > 0:
        r_real = np.concatenate([r_real, np.zeros(z0)])
    return r_real

def interlaces_weak(r1, r2):
    # Check r1 << r2: len(r2) in {len(r1), len(r1)+1}, rightmost of r2 >= rightmost of r1
    m, n = len(r1), len(r2)
    tol = 1e-7
    if n == m + 1:
        return all(r2[q] <= r1[q] + tol and r1[q] <= r2[q+1] + tol for q in range(m))
    if n == m:
        return all(r1[q] <= r2[q] + tol for q in range(m)) and all(r2[q] <= r1[q+1] + tol for q in range(m-1))
    return False

# Build P[(k, i)][r] for r = 1..k, and A[(k, i)]
P = {
    (1, 0): {1: [mp.mpf(1), mp.mpf(1)]},
    (1, 1): {1: [mp.mpf(1)]},
}
A = {
    (1, 0): [mp.mpf(1), mp.mpf(1)],
    (1, 1): [mp.mpf(1)],
}

for k in range(2, 13):
    for i in range(k + 1):
        L = 3*k - 2 - i
        # Build A_{k,i} directly for verification
        uA = [mp.mpf(0)]
        if i <= k - 1: uA = poly_add_mp(uA, apply_SL_mp(A[(k-1, i)], L - 1))
        if i >= 1: uA = poly_add_mp(uA, A[(k-1, i-1)])
        A[(k, i)] = apply_SL_mp(uA, L)

        P[(k, i)] = {}
        for r in range(1, k):
            u = [mp.mpf(0)]
            if i <= k - 1: u = poly_add_mp(u, apply_SL_mp(P[(k-1, i)][r], L - 2))
            if i >= 1: u = poly_add_mp(u, P[(k-1, i-1)][r])
            P[(k, i)][r] = apply_SL_mp(u, L - 1)
        # r = k
        pk = [mp.mpf(0)]
        if i <= k - 1:
            pk = poly_add_mp(pk, apply_SL_mp(mul_x_mp(A[(k-1, i)]), L))
            pk = poly_add_mp(pk, mul_x_mp(apply_SL_mp(A[(k-1, i)], L - 2)))
        if i >= 1:
            pk = poly_add_mp(pk, mul_x_mp(A[(k-1, i-1)]))
        P[(k, i)][k] = pk

        # Verify sum_r P[(k,i)][r] == A[(k,i)]!
        tot = [mp.mpf(0)]
        for r in range(1, k + 1):
            tot = poly_add_mp(tot, P[(k, i)][r])
        assert tot == A[(k, i)], f"Sum mismatch at ({k},{i})!"

    # Check all-pairs interlacing for this k!
    nonreal = 0
    pair_fails = 0
    order_ok = 0
    order_tot = 0
    first_fail = None
    for i in range(k + 1):
        rts = {r: roots_mp_fast(P[(k, i)][r]) for r in range(1, k + 1)}
        for r in range(1, k + 1):
            if rts[r] is None:
                nonreal += 1
        for a in range(1, k + 1):
            for b in range(a + 1, k + 1):
                if rts[a] is None or rts[b] is None:
                    continue
                order_tot += 1
                ok_ab = interlaces_weak(rts[a], rts[b])
                ok_ba = interlaces_weak(rts[b], rts[a])
                if ok_ab:
                    order_ok += 1
                if not (ok_ab or ok_ba):
                    pair_fails += 1
                    if first_fail is None:
                        first_fail = (k, i, a, b)
    print(f"  k={k:2d}: nonreal={nonreal}, pair_fails={pair_fails}/{order_tot}, (P_a << P_b for a<b: {order_ok}/{order_tot}), first_fail={first_fail}")
