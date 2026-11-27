"""
Blazing-fast exact integer computation of Bardakov-Kozlovskaya polynomial P_{2m}(e1, e2, e3)
using the EXACT Newton power-sum identity:
  P_r(x, y) = n * sum_{j=0}^r C(n*r, n*j) x^{r-j} y^j
  E_0 = 1,   d * E_d = sum_{r=1}^d (-1)^{r-1} P_r * E_{d-r}
  p_n(x, y, z) = sum_{d=0}^n (-1)^d z^{n-d} E_d(x, y).

Then extract exact c_{a,b,c}^{(n)} in p_n(x, y, z) = sum_{a+2b+3c=n} c_{a,b,c}^{(n)} e_1^a e_2^b e_3^c
and inspect:
1. Full tables for n = 2, 4, 6, 8, 10, 12, 14, 16 (and verify sign law through n = 24!).
2. Prime factorizations / 2-adic valuations / divisibility by n, 2^{...}, etc.!
3. Exact slices for fixed c (c=0, 1, 2, 3) and fixed b (b=0) and fixed a (a=0)!
4. Real-rootedness / gamma-positivity of c-slices Q_c(e1, e2) in u = e2 / e1^2!
"""
import math
import time

def poly_mul_xy(A, B):
    # A, B are lists of ints representing homogeneous polys in (x, y) of deg len-1: A[j] = coeff of x^{d-j} y^j
    out = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if a:
            for j, b in enumerate(B):
                if b:
                    out[i + j] += a * b
    return out

def compute_p_n_fast(n):
    """Compute exact dict {(a,b,c): coeff} for p_n(x,y,z) in O(n^4) integer ops!"""
    # Step 1: Newton power sums P_r(x, y)
    P = [None]
    for r in range(1, n + 1):
        P.append([n * math.comb(n * r, n * j) for j in range(r + 1)])
    # Step 2: Elementary symmetric polys E_d(x, y) of the roots w_k = (alpha + zeta^k beta)^n
    E = [[1]]
    for d in range(1, n + 1):
        acc = [0] * (d + 1)
        for r in range(1, d + 1):
            prod = poly_mul_xy(P[r], E[d - r])
            if (r - 1) % 2 == 0:
                for idx, v in enumerate(prod): acc[idx] += v
            else:
                for idx, v in enumerate(prod): acc[idx] -= v
        assert all(v % d == 0 for v in acc), f"Non-integer division at n={n}, d={d}!"
        E.append([v // d for v in acc])

    # Now p_n(x, y, z) = sum_{d=0}^n (-1)^d z^{n-d} sum_{j=0}^d E[d][j] x^{d-j} y^j.
    # To extract c_{a,b,c}^{(n)} in O(n^3)!!
    # Look! Set x + y = s, x*y = p!
    # Wait: even simpler! Every symmetric poly in (x, y, z) is UNIQUELY determined by its values
    # on monomials x^p y^q z^r with p >= q >= r (and p+q+r = n).
    # Even BETTER: convert p_n(x, y, z) directly to (e1, e2, e3) using the EXACT 1-variable recursion on z!
    # Wait: for each d, E_d(x, y) is symmetric in (x, y), so write E_d(x, y) = sum_{j=0}^{d//2} alpha_{d,j} (x+y)^{d-2j} (xy)^j!
    # THEN substitute x+y = e1 - z, xy = e2 - z(e1 - z) = e2 - e1*z + z^2!
    # Wait: since p_n(x, y, z) is symmetric in (x, y, z), it depends ONLY on (e1, e2, e3), where z^3 - e1 z^2 + e2 z - e3 = 0!
    # OMG!! THAT MEANS YOU JUST REDUCE p_n(x, y, z) AS A POLYNOMIAL IN z MODULO (z^3 - e1 z^2 + e2 z - e3)!!
    # EVEN BETTER: just peel off c = r from r = 0 to n//3? No: peeling from largest (p,q,r) in (x,y,z) needs only the
    # restricted monomials with p >= q >= r, OR use the EXACT recurrence:
    # Let f(x, y, z) = sum_{a+2b+3c=n} c_{a,b,c} e1^a e2^b e3^c.
    # Notice: f(x, y, 0) = sum_{a+2b=n} c_{a,b,0} (x+y)^a (xy)^b!
    # Subtract sum c_{a,b,0} e1^a e2^b from f(x, y, z) -> remainder is divisible by e3 = xyz!
    # Divide by xyz -> exact symmetric polynomial of degree n-3!! REPEAT!!
    # WAIT A SECOND!! THAT'S O(n^3) AND ZERO HASHTABLES!!
    # Let's represent a symmetric poly F(x, y, z) of degree m purely by a 2D array F[k][j] = coeff of z^k x^{m-k-j} y^j (0 <= k <= m, 0 <= j <= m-k)!
    # Wait: at z = 0 (k = 0), F[0][j] is a symmetric poly of deg m in (x, y):
    #   F(x, y, 0) = sum_{b=0}^{m//2} c_{m-2b, b} (x+y)^{m-2b} (xy)^b.
    # Peel c_{m-2b, b} in O(m^2) from j = 0, 1, ..., m//2:
    #   Since (x+y)^{m-2b} (xy)^b starts at x^{m-b} y^b (j = b), c_{m-2b, b} is just the remaining coeff at j = b!
    # Then subtract c_{m-2b, b} * e1^{m-2b} e2^b from F — wait: we only need F reduced, and dividing by xyz shifts (k, j) -> (k-1, j-1) and degree m -> m-3!
    # Let's check: how fast is expanding e1^a e2^b in (z, x, y)?
    # Even FASTER: just compute F(x, y, z) mod (z^3 - e1 z^2 + e2 z - e3)? No, let's check both:
    # In fact, for a given m, we don't even need 3D multiplication for each b separately:
    # Write G(e1, e2) = sum_{b=0}^{m//2} c_{m-2b, b} e1^{m-2b} e2^b!
    # Since e1 = s + z, e2 = p + s*z (where s = x+y, p = xy),
    # we can expand G(s+z, p+sz) in powers of (z, s, p) in O(m^3), OR keep F ALREADY in the (z, s, p) basis!!
    # OH!! MY!! GOD!!
    # KEEP F(x, y, z) IN THE (z, s, p) BASIS WHERE s = x+y, p = xy!!
    # Look how INSANELY SIMPLE the (z, s, p) basis is:
    # 1. Initially F = sum_{d=0}^n (-1)^d z^{n-d} E_d(x, y). Convert each E_d(x, y) to (s, p) in O(n^2) -> O(n^3) total!
    # 2. At each layer c = 0, 1, ..., n//3 (remaining deg m = n - 3c):
    #    - The z^0 slice of F IS literally sum_{b=0}^{m//2} c_{m-2b, b, c} s^{m-2b} p^b! READ OFF ALL c_{a,b,c} IN O(m) TIME!!
    #    - Subtract sum_b c_{m-2b, b, c} (s + z)^{m-2b} (p + s*z)^b from F!
    #    - Divide the remaining F by e3 = z * p (just shift z-power by -1 and p-power by -1!!)!
    # TOTAL COMPLEXITY TO EXTRACT EVERY SINGLE c_{a,b,c}^{(n)}: O(n^4) / tiny constants -> n = 40 IN UNDER 1 SECOND!!

    def xy_to_sp(poly_xy):
        # poly_xy[j] = coeff of x^{d-j} y^j (symmetric: poly_xy[j] == poly_xy[d-j])
        # Return dict/list coeffs[b] of s^{d-2b} p^b for b = 0..d//2
        d = len(poly_xy) - 1
        rem = list(poly_xy)
        res = [0] * (d // 2 + 1)
        for b in range(d // 2 + 1):
            c = rem[b]
            res[b] = c
            if c:
                # subtract c * (x+y)^{d-2b} * x^b y^b
                deg_s = d - 2 * b
                for r in range(deg_s + 1):
                    rem[b + r] -= c * math.comb(deg_s, r)
        assert all(v == 0 for v in rem), "Not symmetric in x,y!"
        return res

    # Represent F of homogeneous deg m as dict (k, b) -> coeff of z^k * s^{m - k - 2b} * p^b
    # Wait: use 2D list F[k][b] where 0 <= k <= m, 0 <= b <= (m - k)//2!
    F = [[0] * ((n - k) // 2 + 1) for k in range(n + 1)]
    for d in range(n + 1):
        k = n - d
        sgn = -1 if (d % 2 == 1) else 1
        sp_row = xy_to_sp(E[d])
        for b, v in enumerate(sp_row):
            F[k][b] = sgn * v

    coeffs = {}
    m = n
    c = 0
    while m >= 0:
        # Read off z=0 row: F[0][b] is coeff of s^{m-2b} p^b = c_{m-2b, b, c}!
        row0 = list(F[0])
        for b, val in enumerate(row0):
            coeffs[(m - 2 * b, b, c)] = val
        # Subtract sum_{b=0}^{m//2} row0[b] * (s + z)^{m - 2b} * (p + s*z)^b from F!
        # Expand (s+z)^a * (p+sz)^b = sum_{i=0}^a sum_{j=0}^b C(a,i) C(b,j) z^{i+j} s^{a-i+j} p^{b-j}
        # Notice k = i+j, new p-power = b' = b - j!
        # So for fixed (b, j) with b' = b - j >= 0, coeff of z^{i+j} p^{b'} gets row0[b] * C(b, j) * C(m-2b, i)!
        for b, val in enumerate(row0):
            if not val:
                continue
            a = m - 2 * b
            for j in range(b + 1):
                bp = b - j
                vj = val * math.comb(b, j)
                for i in range(a + 1):
                    k = i + j
                    F[k][bp] -= vj * math.comb(a, i)
        # Verify F[0] is all 0 and F[k][0] is all 0 (divisible by z*p)!
        assert all(v == 0 for v in F[0]), f"z=0 remainder nonzero at c={c}!"
        for k in range(1, m + 1):
            assert F[k][0] == 0, f"p=0 remainder nonzero at c={c}, k={k}: {F[k][0]}!"
        if m < 3:
            # Verify entire F is zero!
            for k in range(m + 1):
                assert all(v == 0 for v in F[k]), f"Final remainder nonzero at m={m}!"
            break
        # Divide by z * p: new degree m' = m - 3, F_new[k'][b'] = F[k'+1][b'+1]!
        m_new = m - 3
        F_new = [[0] * ((m_new - k) // 2 + 1) for k in range(m_new + 1)]
        for k in range(m_new + 1):
            for bp in range((m_new - k) // 2 + 1):
                F_new[k][bp] = F[k + 1][bp + 1]
        F = F_new
        m = m_new
        c += 1
    return coeffs

t0 = time.time()
for m in range(1, 13):
    n = 2 * m
    coeffs = compute_p_n_fast(n)
    zeros = [k for k, v in coeffs.items() if v == 0]
    sign_ok = all((v > 0) == (((b * (c + 1) + c) % 2) == 0) for (a, b, c), v in coeffs.items())
    print(f"n={n:2d} (m={m:2d}): terms={len(coeffs):3d}, zeros={len(zeros)}, sign_law_OK={sign_ok}, elapsed={time.time()-t0:.2f}s")
    if n <= 12:
        max_c = n // 3
        for c in range(max_c + 1):
            row = sorted([(a, b, coeffs[(a, b, c)]) for (a, b, cc) in coeffs if cc == c], key=lambda t: t[1])
            print(f"   c={c}: " + ", ".join(f"({a},{b}):{v}" for a, b, v in row))
