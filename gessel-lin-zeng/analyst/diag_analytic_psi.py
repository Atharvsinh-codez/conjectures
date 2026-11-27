"""DIAGNOSTIC FOR ROUTE 1 (ANALYTIC / TWO-POLE + ALTERNATING-SERIES SIGN GRID ON psi_{k,i}(s)):
Recall:
  psi_{k, i}(s) = e^s / (1 + e^s)^{3k - i + 1} * A_{k, i}(-e^s),  s in R (x = -e^s < 0).
  Recurrence: psi_{k, i}(s) = 1/(1 + e^s) * (psi''_{k-1, i}(s) + psi'_{k-1, i-1}(s)),  psi_{0,0}(s) = e^s/(1+e^s).
  For s < 0: psi_{k,i}(s) = sum_{n>=1} (-1)^{n-1} P_{k,i}(n) e^{n s}.
  For s > 0: psi_{k,i}(s) = (-1)^{3k-i} sum_{m>=k+1} (-1)^m |P_{k,i}(-m)| e^{-m s}  (verify sign!).

Check for k = 6, 10, 15, 20 and slices i = 0, 4, k//2, k-3, k:
1. Exact roots s_1 < s_2 < ... < s_d of A_{k,i}(-e^s) (d = 2k - i - 1).
2. Count how many consecutive sign-alternating sample points are certified by:
   - 1-term dominance of s < 0 series (at s = -0.5*log(P(n+1)/P(n-1)), wait: check peak ratio!)
   - 1-term dominance of s > 0 series
   - Two-pole principal part 2 Re Q_{k,i}(1 / (s - pi*i))!
3. Inspect the EXACT Laurent polynomial Q_{k,i}(w) at s = pi*i:
   What are its coefficients c_{k, i, r}? Does 2 Re Q_{k,i}(1/(s-pi*i)) match psi_{k,i}(s)
   with tiny relative error across the ENTIRE central gap between the s<0 and s>0 series?!
"""
import math
import mpmath as mp

mp.mp.dps = 50

# Exact A_{k, i}(x)
def apply_S(L, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r+1); res[r+1] += c*(L-r)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def poly_add(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

KMAX = 20
A = [[None]*(k+1) for k in range(KMAX+1)]
A[0][0] = [1]
for k in range(1, KMAX+1):
    for i in range(k+1):
        u = [0]
        if i <= k-1: u = poly_add(u, apply_S(3*k-3-i, A[k-1][i]))
        if i >= 1: u = poly_add(u, A[k-1][i-1])
        A[k][i] = apply_S(3*k-2-i, u)

def eval_poly_mp(p, x):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

def eval_psi(k, i, s):
    es = mp.exp(s)
    return es / ((1 + es)**(3*k - i + 1)) * eval_poly_mp(A[k][i], -es)

def roots_in_s(k, i):
    """Exact high-precision roots of A_{k,i}(-e^s) = 0 in s = log(-x)."""
    p = A[k][i]
    d = len(p) - 1
    if d == 0: return []
    # Scale x = g * u where g = (p[0]/p[-1])**(1/d)
    g = (mp.mpf(p[0]) / mp.mpf(p[-1])) ** (mp.mpf(1) / d)
    cs = [mp.mpf(p[r]) * (g ** r) for r in range(d + 1)]
    rts = mp.polyroots(cs[::-1], maxsteps=400, extraprec=50)
    assert max(abs(mp.im(r)) for r in rts) < mp.mpf('1e-25')
    xs = sorted([mp.re(r) * g for r in rts])  # all negative, most negative first
    # x = -e^s => s = log(-x), increasing as x goes from 0^- to -inf!
    ss = sorted([mp.log(-x) for x in xs])
    return ss

# Compute P_{k, i}(n) for integer n (both positive and negative!) from A_{k, i}(x) directly:
# Recall: sum_{n>=1} P_{k,i}(n) x^n = x A_{k,i}(x) / (1 - x)^{D+1}, where D = 3k - i!
# Therefore: P_{k,i}(n) = sum_{r=0}^{d} A_{k,i,r} * binom(n - 1 - r + D, D)  (polynomial in n of deg D!)
# Wait: binom(n - 1 - r + D, D) = prod_{j=1}^D (n - 1 - r + j) / D! — VALID FOR ALL n in Z (positive AND negative)!!
def eval_P(k, i, n):
    D = 3*k - i
    p = A[k][i]
    # Use exact rational or mpmath
    Dfact = mp.factorial(D)
    tot = mp.mpf(0)
    for r, coef in enumerate(p):
        # prod_{j=1}^D (n - 1 - r + j)
        prod = mp.mpf(1)
        base = n - 1 - r
        for j in range(1, D + 1):
            prod *= (base + j)
        tot += coef * prod
    return tot / Dfact

# Verify s>0 series at k=4, i=2, s=2.0!
s_chk = mp.mpf("2.0")
v_exact = eval_psi(4, 2, s_chk)
v_ser = sum((-1)**(m - 1) * (-eval_P(4, 2, -m)) * mp.exp(-m * s_chk) for m in range(1, 30))
# Check sign!
print(f"[Sanity k=4,i=2,s=2: psi={mp.nstr(v_exact,10)}, sum (-1)^m P(-m) e^(-ms)={mp.nstr(sum((-1)**m * eval_P(4,2,-m)*mp.exp(-m*s_chk) for m in range(1,30)),10)}]")
# Verify P(4, 2, -m) == 0 for m = 0, 1, ..., 4!
print(f"[Sanity P(4,2,-m) for m=0..7: {[mp.nstr(eval_P(4,2,-m), 4) for m in range(8)]}]")

# Compute EXACT principal part Q_{k,i}(w) = sum_{r=1}^{D+1} c_r w^r at s = pi*i via recurrence!
# psi_{k,i}(pi*i + u) = sum_{r=1}^{D+1} c_r u^{-r} + O(1)
#Wait: psi_{k,i}(pi*i + u) = e^{pi*i + u} / (1 + e^{pi*i + u})^{D+1} * A_{k,i}(-e^{pi*i + u})
#   = -e^u / (1 - e^u)^{D+1} * A_{k,i}(e^u)!!
# OH!! LOOK HOW SIMPLE THAT IS!!
#   f(u) := -e^u / (1 - e^u)^{D+1} * A_{k,i}(e^u) = sum_{r=1}^{D+1} c_r u^{-r} + O(1)!!
# We can compute c_1, ..., c_{D+1} directly by multiplying truncated Taylor series of u/(1-e^u) and e^u A_{k,i}(e^u) up to u^D!!
def laurent_principal_part(k, i):
    """Return list c[1..D+1] (0-indexed: c[r] is coeff of u^{-r}, c[0]=0) for psi_{k,i}(pi*i + u)."""
    D = 3*k - i
    M = D + 1  # need u^0 .. u^D in [u/(1-e^u)]^{D+1} * (-e^u A_{k,i}(e^u))
    # First: g(u) = u / (1 - e^u) = -1 + u/2 - u^2/12 + ... = -sum_{m>=0} B_m / m! u^m
    # Easier: (e^u - 1)/u = sum_{m=0}^D u^m / (m+1)!, then invert!
    h = [mp.mpf(1) / mp.factorial(m + 1) for m in range(M)]
    # Invert h(u) mod u^M: h * g = -1 (since u/(1-e^u) = -1 / ((e^u-1)/u))
    g = [mp.mpf(0)] * M
    g[0] = mp.mpf(-1)
    for m in range(1, M):
        s = sum(h[j] * g[m - j] for j in range(1, m + 1))
        g[m] = -s
    # Power g(u)^{D+1} mod u^M by binary exponentiation
    def mul_trunc(a, b):
        res = [mp.mpf(0)] * M
        for p in range(M):
            if a[p] == 0: continue
            for q in range(M - p):
                res[p + q] += a[p] * b[q]
        return res
    gp = [mp.mpf(0)] * M; gp[0] = mp.mpf(1)
    base = list(g); expn = D + 1
    while expn:
        if expn & 1: gp = mul_trunc(gp, base)
        expn >>= 1
        if expn: base = mul_trunc(base, base)
    # Compute -e^u A_{k,i}(e^u) = -sum_{r=0}^d A_r e^{(r+1)u} mod u^M
    # [u^m] (-sum A_r e^{(r+1)u}) = -(1/m!) sum_{r=0}^d A_r (r+1)^m
    p = A[k][i]
    num = [mp.mpf(0)] * M
    for m in range(M):
        s = sum(mp.mpf(coef) * (mp.mpf(r + 1) ** m) for r, coef in enumerate(p))
        num[m] = -s / mp.factorial(m)
    prod = mul_trunc(gp, num)
    # prod[m] is coeff of u^{m - (D+1)} = u^{-(D+1-m)}!
    c = [mp.mpf(0)] * (D + 2)
    for m in range(M):
        c[D + 1 - m] = prod[m]
    return c

def eval_two_pole(c, s):
    """Evaluate 2 Re sum_{r=1}^{D+1} c[r] / (s - pi*i)^r."""
    w = mp.mpc(1) / mp.mpc(s, -mp.pi)
    val = mp.mpc(0)
    wp = w
    for r in range(1, len(c)):
        val += c[r] * wp
        wp *= w
    return 2 * mp.re(val)

print("\n=== COVERAGE ANALYSIS: 1-TERM SERIES (s<0, s>0) + TWO-POLE APPROX ON ALL ROOTS ===")
for k in [6, 10, 15, 20]:
    slices = sorted(set([0, min(4, k), k//2, max(0, k-3), k]))
    print(f"\n--- k = {k} ---")
    for i in slices:
        ss = roots_in_s(k, i)
        d = len(ss)
        if d == 0: continue
        c = laurent_principal_part(k, i)
        # Check midpoints (and outer points) between consecutive roots!
        pts = [ss[0] - 0.5] + [0.5*(ss[j] + ss[j+1]) for j in range(d - 1)] + [ss[-1] + 0.5]
        two_pole_ok = 0; two_pole_max_rel_err = mp.mpf(0)
        two_pole_bad_s = []
        for idx, sv in enumerate(pts):
            ex = eval_psi(k, i, sv)
            tp = eval_two_pole(c, sv)
            err = abs(ex - tp) / (abs(ex) + mp.mpf('1e-100'))
            if (ex > 0 and tp > 0) or (ex < 0 and tp < 0):
                two_pole_ok += 1
            else:
                two_pole_bad_s.append((idx, float(sv), float(err)))
            if abs(sv) < 2.5 and err > two_pole_max_rel_err:
                two_pole_max_rel_err = err
        # Check 1-term series dominance on s < 0: for n = 1, 2, ..., when is P(n) e^{n s} > sum_{m != n} P(m) e^{m s}?
        # More cleanly: for an alternating series sum (-1)^{n-1} T_n, at a point s where T_n is the strict peak
        # (T_n > T_{n-1} and T_n > T_{n+1}) AND T_n is unimodal in n, the alternating series ALWAYS HAS THE SIGN OF (-1)^{n-1}
        # PROVIDED T_n > T_{n-1} + T_{n+1} (WAIT: EVEN BETTER! For ANY unimodal alternating series with peak at n,
        # the left tail sum_{m<n} has sign (-1)^n and magnitude <= T_{n-1}, and the right tail sum_{m>n} has sign (-1)^n
        # and magnitude <= T_{n+1}! SO T_n > T_{n-1} + T_{n+1} IS 100% SUFFICIENT FOR SIGN (-1)^{n-1}!!)
        # WAIT A SECOND!! CHECK HOW MANY n HAVE AN s WHERE T_n(s) > T_{n-1}(s) + T_{n+1}(s)!!
        # Maximize T_n / (T_{n-1} + T_{n+1}) over s < 0!
        # Let u = e^s. Ratio is 1 / ( (P(n-1)/P(n))/u + (P(n+1)/P(n))*u ).
        # By AM-GM, the minimum of the denominator over u > 0 is 2 * sqrt(P(n-1) P(n+1)) / P(n)!!
        # SO A SAMPLE POINT s_n WITH STRICT 1-TERM DOMINANCE EXISTS <==> P(n)^2 > 4 * P(n-1) * P(n+1)!!
        # WAIT A SECOND!! STOP!! CHECK THAT!!
        n_left = 0; last_s_left = -999.0
        for n in range(1, 50):
            Pn = eval_P(k, i, n)
            Pm = eval_P(k, i, n - 1) if n > 1 else mp.mpf(0)
            Pp = eval_P(k, i, n + 1)
            if Pn * Pn > 4 * Pm * Pp:
                n_left += 1
                if n > 1:
                    last_s_left = float(0.5 * mp.log(Pm / Pp))
            else:
                break
        # Same for s > 0: m = k+1, k+2, ...
        n_right = 0; first_s_right = 999.0
        for m in range(k + 1, k + 60):
            Pm = abs(eval_P(k, i, -m))
            Pl = abs(eval_P(k, i, -(m - 1))) if m > k + 1 else mp.mpf(0)
            Pr = abs(eval_P(k, i, -(m + 1)))
            if Pm * Pm > 4 * Pl * Pr:
                n_right += 1
                if m > k + 1:
                    first_s_right = float(-0.5 * mp.log(Pl / Pr))
            else:
                break
        # Count how many sample points in pts lie inside [last_s_left, first_s_right] and whether two_pole certifies ALL of them!
        bad_in_mid = [b for b in two_pole_bad_s if last_s_left <= b[1] <= first_s_right]
        print(f"  i={i:2d} (d={d:2d}): roots in [{float(ss[0]):+.2f}, {float(ss[-1]):+.2f}] | "
              f"2pole ok={two_pole_ok}/{d+1} (bad_in_mid={len(bad_in_mid)}, max_err(|s|<2.5)={float(two_pole_max_rel_err):.2e}) | "
              f"1term L={n_left} (to s={last_s_left:+.2f}), R={n_right} (from s={first_s_right:+.2f})", flush=True)
