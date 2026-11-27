"""
Quick check (2s):
1. Fix direction bug in Test 1: does f << P_c g imply g << S_L f?!
   (In GLZ: h~_{k-1,i-1} << D_{k+1} h~_{k-1,i}, i.e. g << S_L f, NOT S_L f << g!)
2. Check if A_{k,i} << P_{d_prev} A_{k,i-1} (exact degree d_prev - 1!) holds on GLZ!
3. Check the EXACT strongest c = c_min(k,i) for A_{k,i} << P_c A_{k,i-1}, AND check
   whether f << P_{d} g ALWAYS implies g << S_L f for ANY f, g!
"""
import numpy as np

def poly_trim(p):
    p = list(p)
    while len(p) > 1 and abs(p[-1]) < 1e-13 * max(1.0, max(abs(x) for x in p)):
        p.pop()
    return p

def apply_SL(p, L):
    out = [0.0] * (len(p) + 1)
    for m, c in enumerate(p):
        out[m] += (1 + m) * c
        out[m + 1] += (L - m) * c
    return poly_trim(out)

def apply_Pc(p, c):
    out = [(c - m) * val for m, val in enumerate(p)]
    return poly_trim(out)

def roots_neg(p):
    r = np.roots(list(reversed(p)))
    if np.max(np.abs(r.imag)) > 1e-7:
        return None
    return np.sort(r.real)

def interlaces(p, q):
    """p << q: deg q in {deg p, deg p + 1}, rightmost root of q > rightmost root of p."""
    rp = roots_neg(p)
    rq = roots_neg(q)
    if rp is None or rq is None:
        return False
    m, n = len(rp), len(rq)
    if n == m + 1:
        return all(rq[k] < rp[k] < rq[k + 1] for k in range(m))
    elif n == m:
        return all(rp[k] < rq[k] for k in range(m)) and all(rq[k] < rp[k + 1] for k in range(m - 1))
    return False

rng = np.random.default_rng(7)
print("=== FIXED DIRECTION: Does f << P_c g ALWAYS imply g << S_L f? ===")
for d in [4, 6, 8]:
    for c_mode, c_val in [("c=d (deg d-1)", d), ("c=d+1", d+1), ("c=d+3", d+3)]:
        for L in [d + 1, d + 2, d + 4]:
            fails = 0
            trials = 400
            for _ in range(trials):
                sg = -np.sort(np.exp(rng.normal(0, 1.5, d)))[::-1]
                g = np.poly(sg)[::-1]
                Pg = apply_Pc(g, c_val)
                rPg = roots_neg(Pg)
                if c_val == d:
                    # deg Pg = d - 1; pick f of deg d-1 with f << Pg: rf_1 < rPg_1 < rf_2 < rPg_2 ...
                    rf = [rPg[0] - np.exp(rng.normal(0, 1.0))] + [rng.uniform(rPg[k], rPg[k+1]) for k in range(d - 2)]
                else:
                    # deg Pg = d; pick f of deg d-1 with f << Pg: rPg_0 < rf_0 < rPg_1 < ... < rf_{d-2} < rPg_{d-1}
                    rf = [rng.uniform(rPg[k], rPg[k+1]) for k in range(d - 1)]
                f = np.poly(rf)[::-1]
                Sf = apply_SL(f, L)
                if not interlaces(g, Sf):
                    fails += 1
            print(f"  d={d}, {c_mode:14s}, L={L:2d}: fails(g << S_L f) = {fails}/{trials}")

# Also test: on GLZ (exact integer via floats for k<=11), does A_{k,i} << P_{2k-i} A_{k,i-1} hold?!
# AND what about g << S_{L-1} f vs f << P_d g?
print("\n=== GLZ EXACT CHECK: A_{k,i} << P_{d_prev} A_{k,i-1} (c = 2k-i) ===")
A = {(1, 0): [1.0, 2.0], (1, 1): [1.0]}
for k in range(2, 12):
    for i in range(k + 1):
        L = 3*k - 2 - i
        u = [0.0]
        if i <= k - 1:
            t1 = apply_SL(A[(k-1, i)], L - 1)
            u = [0.0]*max(len(u), len(t1))
            for m, v in enumerate(t1): u[m] += v
        if i >= 1:
            t2 = A[(k-1, i-1)]
            if len(t2) > len(u): u += [0.0]*(len(t2)-len(u))
            for m, v in enumerate(t2): u[m] += v
        A[(k, i)] = apply_SL(u, L)

tot = 0
fails = 0
for k in range(2, 12):
    for i in range(1, k + 1):
        tot += 1
        d_prev = 2*k - i
        Pg = apply_Pc(A[(k, i-1)], d_prev)
        if not interlaces(A[(k, i)], Pg):
            fails += 1
            print(f"  FAIL A_{{{k},{i}}} << P_{{{d_prev}}} A_{{{k},{i-1}}}")
print(f"  Result for A_{{k,i}} << P_{{2k-i}} A_{{k,i-1}} (k=2..11): {fails}/{tot} fails!")
# Also check A_{k,i} << A'_{k,i-1} (ordinary derivative)!
fails_d = 0
for k in range(2, 12):
    for i in range(1, k + 1):
        g = A[(k, i-1)]
        gp = [m * g[m] for m in range(1, len(g))]
        if not interlaces(A[(k, i)], gp):
            fails_d += 1
            if fails_d <= 5:
                print(f"  (note: A_{{{k},{i}}} << A'_{{{k},{i-1}}} fails)")
print(f"  Result for A_{{k,i}} << A'_{{{k},{i-1}}} (k=2..11): {fails_d}/{tot} fails!")
# Also check reversed order A'_{k,i-1} << A_{k,i}
fails_dr = 0
for k in range(2, 12):
    for i in range(1, k + 1):
        g = A[(k, i-1)]
        gp = [m * g[m] for m in range(1, len(g))]
        if not interlaces(gp, A[(k, i)]):
            fails_dr += 1
print(f"  Result for A'_{{{k},{i-1}}} << A_{{k,i}} (k=2..11): {fails_dr}/{tot} fails!")
