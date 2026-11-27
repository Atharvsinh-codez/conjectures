"""
BLAZING-FAST exact integer c_{a,b,c}^{(n)} extractor + TEST THE GAMMA/ROOT HYPOTHESIS!
For each even n = 2m and each c = 0..floor(n/3):
  f_{n, c}(u) := Q_{n, c}(1+u, u) = sum_b c_{n-2b-3c, b, c}^{(n)} u^b (1+u)^{n-3c-2b}.
Check:
1. Exact c_{a,b,c}^{(2m)} for n = 2, 4, ..., 24 (and check sign law (-1)^{b(c+1)+c}).
2. For EVERY (n, c), check the ROOTS of f_{n, c}(u)!!
   - For c ODD: are ALL roots of f_{n,c}(u) negative real?! (which IMMEDIATELY implies all c_{a,b,c} < 0!)
   - For c EVEN: are ALL roots of f_{n,c}(u) positive real?! (which IMMEDIATELY implies (-1)^b c_{a,b,c} > 0!)
   - Even if roots are complex, check WHERE they lie in C!
3. Also check the roots of the z-coefficient polynomials [z^r] p_n(1, u, z) and the
   exact generating function P_n(1+u, u, e3) = sum_c f_{n,c}(u) e3^c!!
"""
from math import comb
import numpy as np

def p_coeffs_table(n):
    """Return 3D table P[j][k] = [x^(j-k) y^k z^(n-j)] p_n(x,y,z) via Newton identities."""
    # power sums s_r(x,y) = n * (-1)^{n*r} * sum_{k=0}^r comb(n*r, n*k) x^{r-k} y^k
    powers = [[0]] + [
        [n * ((-1)**(n*r)) * comb(n*r, n*k) for k in range(r + 1)]
        for r in range(1, n + 1)
    ]
    a = [[1]]
    for j in range(1, n + 1):
        row = [0] * (j + 1)
        for r in range(1, j + 1):
            prev = a[j - r]
            pr = powers[r]
            for u, au in enumerate(prev):
                if au == 0: continue
                for v, bv in enumerate(pr):
                    row[u + v] -= au * bv
        a.append([v // j for v in row])
    return a  # a[j][k] = coeff of x^{j-k} y^k z^{n-j}

def extract_e_coeffs_fast(n):
    """
    Fast pure-integer extraction of c_{a,b,c}^{(n)} from p_n(x,y,z)!
    Method: work with the polynomial p(x,y,z) represented as dict (i,j,k) -> val (i>=j>=k not needed, full table).
    Actually, EVEN FASTER: work ONLY on monomials (i,j,k) with i >= j >= k? No, subtracting e1^a e2^b e3^c
    takes time if done naively.
    Better O(n^3) peeling:
      Store p(x,y,z) as 2D array in (i, j) where k = n - i - j (homogeneity!).
      At step c = 0, 1, ..., n//3:
        1. Read p(x, y, 0) (i.e. k = 0, so i + j = n - 3c after dividing by (xyz)^c).
           Since p(x, y, 0) = Q_{n,c}(x+y, xy) is palindromic in (x,y) of degree d = n - 3c:
           Read off gamma coeffs c_{d-2b, b, c} directly from Q(1+u, u) in O(d^2)!
        2. Subtract Q_{n,c}(x+y+z, xy+yz+zx) from the current polynomial (of degree d),
           and divide by xyz (shift (i, j, k) -> (i-1, j-1, k-1))!
      Wait: computing Q_{n,c}(e1, e2) in (x,y,z) takes O(d^3) if we precompute powers of e1, e2, OR
      just use multiplication by e1 = x+y+z and e2 = xy+yz+zx in Horner/DP!
      Let's compute Q_{n,c}(e1, e2) by DP on b:
        Start from e2^0, multiply by e1 repeatedly! Even better: Horner in e2 and e1!
    """
    a_tab = p_coeffs_table(n)
    # Build poly[i][j] = coeff of x^i y^j z^{d - i - j} for current degree d = n
    d = n
    poly = [[0]*(d+1) for _ in range(d+1)]
    for j_tot in range(n+1):
        for ky in range(j_tot+1):
            ix = j_tot - ky
            poly[ix][ky] = a_tab[j_tot][ky]

    def mul_e1(P, deg):
        # multiply homogeneous poly P(x,y,z) of degree deg by (x + y + z) -> degree deg+1
        out = [[0]*(deg+2) for _ in range(deg+2)]
        for i in range(deg+1):
            for j in range(deg+1-i):
                v = P[i][j]
                if v:
                    out[i+1][j] += v
                    out[i][j+1] += v
                    out[i][j] += v
        return out

    def mul_e2(P, deg):
        # multiply by (xy + xz + yz) -> degree deg+2
        out = [[0]*(deg+3) for _ in range(deg+3)]
        for i in range(deg+1):
            for j in range(deg+1-i):
                v = P[i][j]
                if v:
                    out[i+1][j+1] += v
                    out[i+1][j] += v
                    out[i][j+1] += v
        return out

    coeffs = {} # (a, b, c) -> int
    for c in range(n // 3 + 1):
        d = n - 3 * c
        # Read z=0 slice: i + j = d -> slice[j] = coeff of x^{d-j} y^j
        sl = [poly[d - j][j] for j in range(d + 1)]
        # Verify palindrome
        assert sl == sl[::-1], f"Not symmetric at n={n}, c={c}"
        # Extract gamma coeffs: sum_b g_b u^b (1+u)^{d-2b} = sum_j sl[j] u^j
        rem = sl[:]
        gammas = []
        for b in range(d // 2 + 1):
            gb = rem[b]
            gammas.append(gb)
            coeffs[(d - 2*b, b, c)] = gb
            # subtract gb * u^b * (1+u)^{d-2b}
            m_pow = d - 2*b
            for r in range(m_pow + 1):
                rem[b + r] -= gb * comb(m_pow, r)
        assert all(v == 0 for v in rem)
        if c == n // 3:
            break
        # Build Q_{n,c}(e1, e2) in (x,y,z) via Horner in e2!
        # Q = sum_{b=0}^{d//2} gammas[b] * e1^{d-2b} * e2^b
        # Wait: build powers of e1 up to d? Easier: Horner from b = d//2 down to 0!
        # Actually deg_e1 varies by 2 each step:
        # Let R = gammas[b_max] * (e1 if d%2 else 1)
        # For b from b_max-1 down to 0:
        #   R = R * e2 + gammas[b] * e1^{d - 2b} -> NO! That builds e1^{d-2b} each time!
        # Better: build e1_pow[k] for k=0..d (or just multiply R by e2 and track e1^2!):
        # Wait: Q = e1^{d%2} * sum_{b=0}^{B} gammas[b] * (e1^2)^{B-b} * e2^b!
        # Even simpler: compute P_cur = Q by accumulating e2_pow * e1_pow, or just:
        # Compute e1_pows[r] for r = 0..d? Too much memory? For d<=30, d^3 = 27000 ints, TINY!
        # Even better: Horner in (e2 / e1^2): no, polynomial Horner with TWO vars:
        #   S_0 = gammas[0] (deg 0)
        #   No: if we know e1_pows, wait: e1^r has simple multinomial coeffs!
        #   Just precompute e1_pows[r] for r=0..d (only O(d^3) total)!
        # Actually just compute Q(x,y,z) by multiplying by e1 and e2:
        # Let B = d // 2.
        # Compute E2_cur = [[1]] (deg 0). No, wait:
        # Q = sum_{b=0}^B gammas[b] * e1^{d-2b} * e2^b.
        # Compute from b=0 to B? No, from b=B down to 0:
        # Start with cur_e1 = [[1]] (if d%2==0) else e1 (deg d%2).
        # Wait: we can accumulate Q_poly of degree d directly by tracking E1_pow of deg d-2b?
        # Simplest O(d^3) total per c:
        #   Precompute E1[r] for r=0..d (takes O(d^3) ONCE per c! Wait: reuse across c!).
        # Let's just compute E1_pows once outside the loop!
        break

    # Clean version with precomputed E1_pows!
    E1 = [[[1]]] # E1[r] has degree r
    for r in range(1, n + 1):
        E1.append(mul_e1(E1[-1], r - 1))

    # Re-run peeling loop cleanly!
    poly = [[0]*(n+1) for _ in range(n+1)]
    for j_tot in range(n+1):
        for ky in range(j_tot+1):
            poly[j_tot - ky][ky] = a_tab[j_tot][ky]

    coeffs = {}
    f_polys = {} # (n, c) -> list of coeffs of f_{n,c}(u) = Q_{n,c}(1+u, u)
    for c in range(n // 3 + 1):
        d = n - 3 * c
        sl = [poly[d - j][j] for j in range(d + 1)]
        assert sl == sl[::-1]
        f_polys[c] = sl[:]
        rem = sl[:]
        gammas = []
        for b in range(d // 2 + 1):
            gb = rem[b]
            gammas.append(gb)
            coeffs[(d - 2*b, b, c)] = gb
            m_pow = d - 2*b
            for r in range(m_pow + 1):
                rem[b + r] -= gb * comb(m_pow, r)
        assert all(v == 0 for v in rem)
        if c == n // 3:
            break
        # Compute Q(x,y,z) = sum_{b=0}^{d//2} gammas[b] * E1[d-2b] * e2^b
        # Use Horner in e2:
        # Wait: Q = sum_{b=0}^B gammas[b] * e1^{d-2b} * e2^b
        # Start from b = B: acc = gammas[B] * E1[d-2B] (deg d-2B)
        # Wait: multiplying acc by e2 raises degree by 2, matching E1[d-2(B-1)]? NO, that lowers b by 1 so e1 deg increases by 2!
        # Reverse Horner: start from b = B down to 0, no:
        # Write Q = gammas[0]*e1^d + e2*(gammas[1]*e1^{d-2} + e2*(...))!
        # Let B = d // 2.
        # acc = gammas[B] * E1[d - 2*B] (deg d - 2*B)
        # For b = B-1 down to 0:
        #   acc = mul_e2(acc, cur_deg) -> deg cur_deg + 2 = d - 2*b!
        #   acc += gammas[b] * E1[d - 2*b]!!
        B = d // 2
        cur_deg = d - 2 * B
        acc = [[gammas[B] * E1[cur_deg][i][j] for j in range(cur_deg + 1)] for i in range(cur_deg + 1)]
        for b in range(B - 1, -1, -1):
            acc = mul_e2(acc, cur_deg)
            cur_deg += 2
            gb = gammas[b]
            e1p = E1[cur_deg]
            for i in range(cur_deg + 1):
                for j in range(cur_deg + 1 - i):
                    acc[i][j] += gb * e1p[i][j]
        # Subtract acc from poly and divide by xyz: new_poly[i-1][j-1] = poly[i][j] - acc[i][j]
        # Verify boundary terms (i=0, j=0, i+j=d) all vanish!
        for i in range(d + 1):
            for j in range(d + 1 - i):
                poly[i][j] -= acc[i][j]
                if i == 0 or j == 0 or i + j == d:
                    assert poly[i][j] == 0, f"Divisibility failed at n={n}, c={c}, i={i}, j={j}"
        new_d = d - 3
        new_poly = [[0]*(new_d + 1) for _ in range(new_d + 1)]
        for i in range(new_d + 1):
            for j in range(new_d + 1 - i):
                new_poly[i][j] = poly[i + 1][j + 1]
        poly = new_poly
    return coeffs, f_polys

print("=== Fast exact extraction & root analysis of f_{n,c}(u) = Q_{n,c}(1+u, u) ===")
for n in range(2, 22, 2):
    coeffs, f_polys = extract_e_coeffs_fast(n)
    # Check sign law
    bad_signs = []
    for (a, b, c), val in coeffs.items():
        if val == 0 or (1 if val > 0 else -1) != (-1)**(b*(c+1) + c):
            bad_signs.append((a, b, c, val))
    # Analyze roots of f_{n,c}(u) for each c!
    root_summary = []
    for c in range(n // 3 + 1):
        sl = f_polys[c] # low-to-high (palindromic)
        if len(sl) <= 1:
            continue
        r = np.roots(sl[::-1])
        max_imag = float(np.max(np.abs(r.imag)))
        all_real = max_imag < 1e-6
        if all_real:
            all_neg = bool(np.all(r.real < 0))
            all_pos = bool(np.all(r.real > 0))
            tag = "ALL_NEG" if all_neg else ("ALL_POS" if all_pos else "MIXED_REAL")
        else:
            # Count complex, check real parts
            n_cplx = int(np.sum(np.abs(r.imag) >= 1e-6))
            tag = f"CPLX({n_cplx}/{len(r)}, Re:[{r.real.min():.2f},{r.real.max():.2f}])"
        root_summary.append(f"c={c}:{tag}")
    print(f"n={n:2d}: terms={len(coeffs):3d}, bad_signs={len(bad_signs)} | " + "  ".join(root_summary))
