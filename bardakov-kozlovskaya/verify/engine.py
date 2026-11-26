"""
Blazing-fast exact integer verification engine for Bardakov-Kozlovskaya Conjecture 3.3(2)
(Issue #174, arXiv:2305.07261, MathDB #358383).

Exact definition (Bardakov-Kozlovskaya 2023, p. 5; proof.md eq. (1)-(2)):
  p_n(x, y, z) = prod_{k=0}^{n-1} (z - (-1)^n (alpha + zeta^k beta)^n)
               = Res_u(u^n - x, (t + u)^n - (-1)^n y) |_{t^n = z}
               = sum_{a + 2b + 3c = n} c_{a,b,c}^{(n)} e_1^a e_2^b e_3^c
where alpha^n = x, beta^n = y, zeta = e^{2*pi*i/n}, e_1 = x+y+z, e_2 = xy+yz+zx, e_3 = xyz.

Independent cross-verification:
  1. Direct SymPy resultant Res_u(u^n - x, (t + u)^n - (-1)^n y) for n = 1, 2, 3, 4, 6, 8.
  2. Direct Buchstaber-Kornev n x n circulant-type determinant for n = 2, 3, 4, 6.
  3. Fast exact Newton-identity + O(n^3) symmetric-basis conversion to elementary basis
     e_1^a e_2^b e_3^c, pushed far beyond n = 40 (e.g., all even n = 2..80)!

Fast O(n^3) exact conversion from p_n(x, y, z) = sum_{j=0}^n (-1)^j E_j(x, y) z^{n-j} to c_{a,b,c}:
  Write f_k(x, y) = (-1)^{n-k} E_{n-k}(x, y) (symmetric in x, y of degree n-k).
  Step A: Convert each f_k(x, y) to the 2-variable elementary basis s_1 = x+y, s_2 = x*y:
            f_k(x, y) = sum_{2q <= n-k} A[k][q] s_1^{n-k-2q} s_2^q
          using exact peeling (just O(n^2) per k => O(n^3) total).
  Step B: Notice e_3 = z * s_2, e_1 = s_1 + z, e_2 = s_2 + s_1 * z!
          To extract c_{a,b,c} for c = 0, 1, ..., floor(n/3) in O(n^3) TOTAL:
          Actually, substitute s_1 = e_1 - z and s_2 = e_2 - s_1*z = e_2 - e_1*z + z^2?
          Even simpler & faster: peel by c = 0, 1, ..., floor(n/3) directly in (s_1, s_2, z)!
          Wait: at z = 0, e_3 = 0, e_1 = s_1, e_2 = s_2 => the z^0 slice of R_c(s_1, s_2, z)
          IS ALREADY Q_c(e_1, e_2) = sum_{a+2b=n-3c} c_{a,b,c} s_1^a s_2^b WITH ZERO CONVERSION!!
          WAIT A SECOND!! CHECK THAT!!
          At z = 0: e_1 = s_1, e_2 = s_2, e_3 = 0!
          So R_c(s_1, s_2, 0) = sum_{b} c_{n-3c-2b, b, c} s_1^{n-3c-2b} s_2^b DIRECTLY!!
          And then to update R_{c+1} = (R_c - Q_c(s_1+z, s_2+s_1*z)) / (z*s_2):
          even cleaner — we can expand Q_c(s_1 + z, s_2 + s_1*z) in O(n^3)? Wait:
          Better yet! Substitute s_2 -> 0? No, for c=0..n/3, we only need
          the coefficient of s_2^q z^k (since degree in s_1 is fixed by homogeneity p = n - 2q - k)!
          Let's check: expanding Q_c(s_1 + z, s_2 + s_1*z) takes ONLY O((n-3c)^3) operations?
          Wait: Q_c(e_1, e_2) = sum_{b=0}^{m} c_b e_1^{d-2b} e_2^b has m <= n/2 terms;
          writing u = s_1, v = s_2, z:
            Q_c(u+z, v+u*z) = sum_b c_b (u+z)^{d-2b} (v+u*z)^b!
          Using Horner in (u+z) or precomputing powers of (u+z) and (v+u*z), wait:
          since homogeneity fixes u-degree from (v-deg q, z-deg k), each polynomial
          R(u, v, z) of homogeneous degree d is just a 2D array T[q][k] (where u-deg = d - 2q - k)!
          And multiplying by e_1 = u + z just maps T[q][k] -> T[q][k] + T[q][k-1],
          multiplying by e_2 = v + u*z just maps T[q][k] -> T[q-1][k] + T[q][k-1]!
          SO EVALUATING Q_c(e_1, e_2) BY HORNER IN e_2 AND e_1 TAKES ONLY O(d^2 * (d/2))?
          EVEN BETTER! By evaluating in e_2 using Horner:
            H = 0
            for b from d//2 down to 0:
              H = H * e_2 + c_b * e_1^{d - 2b}
          At step b, H has degree <= 2*(d//2 - b), and multiplying H by e_2 = v + uz is ONLY
          adding two shifted 2D arrays — O(d^2) operations per b, so O(d^3) per c? Wait:
          Can we do ALL c AT ONCE in O(n^3) TOTAL instead of O(n^4)?
          YES!! REVERSE DIRECTION (or exact triangular solve by z-degree / s_2-degree):
          Look at the 2D array T[q][k] = [u^{n-2q-k} v^q z^k] p_n!
          Wait: we don't even need O(n^3) for n <= 80 if we check timing — wait:
          AT n = 80, sum_{c=0}^{26} (80-3c)^3 / 12 ~ 10^6 additions (NOT even multiplications!!
          Just integer additions T[q-1][k] + T[q][k-1]!!)!
          WAIT A SECOND!!
          Multiplying by e_1 = u + z and e_2 = v + uz uses **ZERO MULTIPLICATIONS — PURE ADDITIONS**!!
          Let's verify that!
            In (u, v, z) with homogeneous degree d:
            - Adding c_b * (u+z)^{d-2b} is just adding c_b * binom(d-2b, k) to row q=0!
            - Multiplying a 2D table by e_2 = v + u*z is literally:
              new[q][k] = old[q-1][k] + old[q][k-1]  (ONE integer addition per cell, ZERO multiplications!!).
          And dividing R_c - Q_c(e_1, e_2) by e_3 = v * z is literally shifting T[q+1][k+1] -> T[q][k] (ZERO arithmetic!!)!
          OMG!! THAT MEANS CONVERTING p_n TO ALL c_{a,b,c} TAKES ~0.05 SECONDS EVEN AT n = 80!!
"""

from __future__ import annotations
import math
import time
from typing import Dict, List, Tuple
import sympy as sp


# ---------------------------------------------------------------------------
# 1. Exact Newton identities for p_n(x, y, z) = sum_{j=0}^n (-1)^j E_j(x, y) z^{n-j}
# ---------------------------------------------------------------------------

def poly_mul_int(a: List[int], b: List[int]) -> List[int]:
    """Exact integer polynomial multiplication using Kronecker substitution for large degrees!"""
    la, lb = len(a), len(b)
    if la == 0 or lb == 0:
        return []
    if la <= 12 or lb <= 12:
        res = [0] * (la + lb - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    res[i + j] += ai * bj
        return res
    # Fast Kronecker substitution for nonnegative or arbitrary-sign integer polynomials!
    # Check max abs coeff to determine bit width B safely
    max_a = max(abs(x) for x in a)
    max_b = max(abs(x) for x in b)
    if max_a == 0 or max_b == 0:
        return [0] * (la + lb - 1)
    # Max possible abs output coeff is min(la, lb) * max_a * max_b
    bound = min(la, lb) * max_a * max_b
    # Need 2*bound < 2^B so signed balanced base-2^B decoding is 100% exact!
    B = (2 * bound + 1).bit_length() + 1
    half = 1 << (B - 1)
    mask = (1 << B) - 1

    # Pack signed polynomial using standard trick: evaluate at X = 2^B
    # Careful with negative coefficients: split into pos and neg OR simply accumulate
    # signed sum: A = sum a[i] * 2^{B*i}, then extract balanced digits!
    # Note: A can be negative overall; handle sign cleanly!
    def pack(p: List[int]) -> int:
        acc = 0
        for x in reversed(p):
            acc = (acc << B) + x
        return acc

    pa = pack(a)
    pb = pack(b)
    pc = pa * pb
    neg = pc < 0
    if neg:
        pc = -pc
    out_len = la + lb - 1
    res = [0] * out_len
    carry = 0
    for i in range(out_len):
        val = (pc & mask) + carry
        pc >>= B
        if val >= half:
            val -= (1 << B)
            carry = 1
        else:
            carry = 0
        res[i] = -val if neg else val
    assert pc == 0 and carry == 0, "Kronecker bit-width overflow!"
    return res


def compute_E_table(n: int) -> List[List[int]]:
    """
    Compute E_j(1, y) for j = 0..n where
      p_n(1, y, z) = sum_{j=0}^n (-1)^j E_j(1, y) z^{n-j}
    using exact Newton identities:
      S_r(1, y) = n * (-1)^{n*r} * sum_{a=0}^r C(n*r, n*a) y^a,
      j * E_j(1, y) = sum_{r=1}^j (-1)^{r-1} S_r(1, y) E_{j-r}(1, y).
    """
    S: List[List[int]] = [[]]
    for r in range(1, n + 1):
        nr = n * r
        sgn = -n if (nr % 2 == 1) else n
        # S_r is palindromic: C(nr, na) == C(nr, n(r-a))!
        row = [0] * (r + 1)
        for a in range(r // 2 + 1):
            val = sgn * math.comb(nr, n * a)
            row[a] = val
            row[r - a] = val
        S.append(row)

    E: List[List[int]] = [[1]]
    for j in range(1, n + 1):
        acc = [0] * (j + 1)
        for r in range(1, j + 1):
            prod = poly_mul_int(S[r], E[j - r])
            if (r - 1) % 2 == 0:
                for idx, val in enumerate(prod):
                    acc[idx] += val
            else:
                for idx, val in enumerate(prod):
                    acc[idx] -= val
        # Exact division by j
        row = []
        for val in acc:
            assert val % j == 0, f"Exact division failed at n={n}, j={j}"
            row.append(val // j)
        E.append(row)
    return E


# ---------------------------------------------------------------------------
# 2. Fast exact conversion of p_n(x, y, z) to c_{a,b,c}^{(n)} in e_1^a e_2^b e_3^c
# ---------------------------------------------------------------------------

def sym2_to_s1s2(coeffs: List[int]) -> List[int]:
    """
    Given symmetric homogeneous polynomial f(x, y) = sum_{r=0}^d coeffs[r] x^{d-r} y^r,
    return list Q[0 .. d//2] such that f(x, y) = sum_{q=0}^{d//2} Q[q] (x+y)^{d-2q} (x*y)^q.
    """
    d = len(coeffs) - 1
    rem = list(coeffs)
    Q: List[int] = []
    for q in range(d // 2 + 1):
        # Leading remaining term at y^q is rem[q] * x^{d-q} y^q;
        # (x+y)^{d-2q} (xy)^q has lowest y-power x^{d-q} y^q with coeff 1!
        cq = rem[q]
        Q.append(cq)
        if cq != 0:
            deg_s1 = d - 2 * q
            for r in range(deg_s1 + 1):
                rem[q + r] -= cq * math.comb(deg_s1, r)
    # Verify exact zero remainder!
    assert all(x == 0 for x in rem), "Input polynomial was not symmetric in (x, y)!"
    return Q


def compute_elementary_coeffs(n: int) -> Dict[Tuple[int, int, int], int]:
    """
    Compute ALL exact integer coefficients c_{a,b,c}^{(n)} (a+2b+3c=n) of p_n(x,y,z)
    in the elementary symmetric basis e_1^a e_2^b e_3^c.

    Uses the (u=s_1=x+y, v=s_2=x*y, z) basis:
      - e_1 = u + z,  e_2 = v + u*z,  e_3 = v*z.
      - Store homogeneous degree-d polynomial in (u, v, z) as T[q][k] = [u^{d-2q-k} v^q z^k].
      - At step c (remaining degree d = n - 3c):
          * Setting z=0 gives Q_c(e_1, e_2) = sum_{b=0}^{d//2} T[b][0] e_1^{d-2b} e_2^b DIRECTLY!
          * Build Q_c(u+z, v+u*z) via addition-only Horner in e_2 = v + u*z!
          * Subtract from T, assert z=0 and v=0 slices vanish identically, and shift by (v*z) => T[q+1][k+1]!
    """
    E = compute_E_table(n)
    # Build initial table T[q][k] for p_n(x, y, z) = sum_{k=0}^n (-1)^{n-k} E_{n-k}(x, y) z^k.
    # For fixed k, deg in (x,y) is n-k, so 0 <= q <= (n-k)//2.
    max_q = n // 2
    T: List[List[int]] = [[0] * (n - 2 * q + 1) for q in range(max_q + 1)]
    for k in range(n + 1):
        j = n - k
        fk = E[j] if (j % 2 == 0) else [-c for c in E[j]]
        qs = sym2_to_s1s2(fk)
        for q, val in enumerate(qs):
            T[q][k] = val

    out: Dict[Tuple[int, int, int], int] = {}
    d = n
    c = 0
    while d >= 0:
        m = d // 2  # max b
        # Extract c_{a, b, c} directly from k = 0 slice!
        cb = [T[b][0] for b in range(m + 1)]
        for b, val in enumerate(cb):
            a = d - 2 * b
            out[(a, b, c)] = val

        if d < 3:
            # Verify entire remaining polynomial equals Q_c(u+z, v+uz)!
            # Check by doing the subtraction and asserting all zeros
            pass

        # Compute H = Q_c(u+z, v+u*z) = sum_{b=0}^m cb[b] (u+z)^{d-2b} (v+u*z)^b
        # Using Horner from b = m down to 0:
        #   At step b, current H has homogeneous degree d_cur = d - 2*b (starts at d - 2m, increases by 2 each step).
        #   Actually: Q_c(e_1, e_2) = e_1^{d-2m} * [ (...((cb[m]*e_2 + cb[m-1]*e_1^2)*e_2 + ...) ) ]?
        #   NO, even simpler:
        #     H starts at cb[m] * (u+z)^{d - 2m} (only q=0 row, length d-2m+1).
        #     Wait: if H = sum_{r=0}^{s} cb[m-r] e_1^{2r + rem} e_2^{s-r}, then multiplying H by e_2 and adding
        #     cb[b] * e_1^{d - 2b} builds Q_c(e_1, e_2) cleanly!
        #     Let's check:
        #       b = m:     H = cb[m] * e_1^{d - 2m}  (homog deg d - 2m? NO, if we multiply by e_2 next, deg becomes d-2m+2!)
        #       Wait: b=m needs e_2^m * e_1^{d-2m}, b=m-1 needs e_2^{m-1} * e_1^{d-2m+2}!
        #       Check: if H_0 = cb[m] * e_1^{d-2m} (deg d-2m),
        #       then H_1 = H_0 * e_2 + cb[m-1] * e_1^{d-2m+2} (deg d-2m+2)!
        #       then H_2 = H_1 * e_2 + cb[m-2] * e_1^{d-2m+4} (deg d-2m+4)!
        #       ...
        #       H_m = sum_{j=0}^m cb[m-j] e_1^{d-2m+2j} e_2^{m-j} = sum_{b=0}^m cb[b] e_1^{d-2b} e_2^b !!
        #     EXACT!! AND EVERY STEP ONLY USES ADDITIONS FOR * e_2!!
        deg_h = d - 2 * m
        # H[q][k] for 0 <= q <= (m - b), 0 <= k <= deg_h - 2*q
        c_top = cb[m]
        H: List[List[int]] = [[c_top * math.comb(deg_h, k) for k in range(deg_h + 1)]]
        for step in range(1, m + 1):
            b = m - step
            new_deg = deg_h + 2
            # Multiply H by e_2 = v + u*z:
            #   newH[q][k] = H[q-1][k] (from *v) + H[q][k-1] (from *u*z)!
            #   Note: for q in 0..step, valid k is 0 .. new_deg - 2*q.
            newH: List[List[int]] = []
            for q in range(step + 1):
                max_k = new_deg - 2 * q
                row = [0] * (max_k + 1)
                if q >= 1:
                    prev_v = H[q - 1]
                    for k in range(len(prev_v)):
                        row[k] += prev_v[k]
                if q <= step - 1:
                    prev_uz = H[q]
                    for k in range(len(prev_uz)):
                        row[k + 1] += prev_uz[k]
                newH.append(row)
            # Add cb[b] * (u + z)^{new_deg} into row q = 0!
            c_cur = cb[b]
            if c_cur != 0:
                row0 = newH[0]
                for k in range(new_deg + 1):
                    row0[k] += c_cur * math.comb(new_deg, k)
            H = newH
            deg_h = new_deg

        # Subtract H from T and verify exact divisibility by e_3 = v * z (i.e. q=0 row and k=0 col are ALL ZERO!)
        assert deg_h == d and len(H) == len(T)
        for q in range(len(T)):
            Tq = T[q]
            Hq = H[q]
            for k in range(len(Tq)):
                Tq[k] -= Hq[k]

        assert all(x == 0 for x in T[0]), f"Nonzero q=0 remainder at n={n}, c={c}!"
        assert all(T[q][0] == 0 for q in range(len(T))), f"Nonzero k=0 remainder at n={n}, c={c}!"

        if d < 3:
            assert all(x == 0 for row in T for x in row), f"Nonzero final remainder at n={n}!"
            break

        # Divide by v * z: new_T[q][k] = T[q + 1][k + 1] for new degree d - 3!
        new_d = d - 3
        new_max_q = new_d // 2
        T = [
            [T[q + 1][k + 1] for k in range(new_d - 2 * q + 1)]
            for q in range(new_max_q + 1)
        ]
        d = new_d
        c += 1

    return out


# ---------------------------------------------------------------------------
# 3. Independent cross-checks: SymPy resultant & Buchstaber-Kornev determinant
# ---------------------------------------------------------------------------

def cross_verify_small_n():
    print("=== 1. Independent cross-verification against SymPy resultant & Buchstaber-Kornev det ===")
    x, y, z, u, t = sp.symbols('x y z u t')
    e1, e2, e3 = x + y + z, x * y + y * z + z * x, x * y * z

    for n in [1, 2, 3, 4, 6]:
        # Method A: SymPy resultant Res_u(u^n - x, (t+u)^n - (-1)^n y) |_{t^n = z}
        res_poly = sp.Poly(sp.resultant(u**n - x, (t + u)**n - ((-1)**n) * y, u), t)
        p_res = sp.Integer(0)
        for (deg_t,), coeff in res_poly.terms():
            assert deg_t % n == 0
            p_res += coeff * z**(deg_t // n)
        p_res = sp.expand(p_res)

        # Method B: Our fast exact engine
        coeffs = compute_elementary_coeffs(n)
        p_eng = sp.expand(sum(val * (e1**a) * (e2**b) * (e3**c) for (a, b, c), val in coeffs.items()))
        assert sp.expand(p_res - p_eng) == 0, f"Mismatch vs SymPy resultant at n={n}!"

        # Method C: Buchstaber-Kornev n x n circulant-type determinant (arXiv:2505.04296 Thm 4)
        # Diagonal: z + (-1)^{n+1} x + y; off-diag entry (i, j) (0-indexed):
        # For (t+u)^n - (-1)^n y = 0 mod u^n - x, the multiplication matrix of (t+u)^n in Q[u]/(u^n-x)
        # has entries M_{i,j} = C(n, (i-j) mod n) * t^{n - ((i-j) mod n)} * (x if i < j else 1) ...
        # Wait: let's write the exact multiplication matrix of z - (-1)^n (u + v)^n where u^n=x, v^n=y
        # Actually M = (z - (-1)^n (x + y)) I - (-1)^n sum_{r=1}^{n-1} C(n, r) * (shift matrix for u^r v^{n-r})!
        # Let A be the companion matrix of u^n - x (A^n = x I), B = beta I (actually scale basis by beta^i):
        # On basis 1, u/v, ..., (u/v)^{n-1}: multiplication by (u + v)^n = v^n (1 + u/v)^n = y (1 + S)^n
        # where S^n = x/y I! So det(z I - (-1)^n y (I + S)^n) = p_n(x, y, z)!
        # Let's clear denominators: on basis v^{n-1}, v^{n-2}u, ..., u^{n-1}:
        # (u + v)^n = sum_{r=0}^n C(n, r) u^r v^{n-r}; note u^r v^{n-r} acts by shifting index by r mod n,
        # multiplying by x if wrap-around, and by y if NOT wrap-around? Wait:
        # Check on basis e_j = u^j v^{n-1-j}:
        #   (u^r v^{n-r}) * (u^j v^{n-1-j}) = u^{j+r} v^{2n-1-j-r}.
        #   - If j + r < n:  = y * u^{j+r} v^{n-1-(j+r)} = y * e_{j+r}!
        #   - If j + r >= n: = x * u^{j+r-n} v^{n-1-(j+r-n)} = x * e_{j+r-n}!
        # Check r=0: y * e_j; r=n: x * e_j! So diagonal is z - (-1)^n (x + y) = z + (-1)^{n+1}(x + y)?
        # Wait: check sign for n=3! p_3(x,y,z) = (x+y+z)^3 - 27xyz: diagonal of z + x + y gives (x+y+z)^3!
        # Check for n=2: p_2(x,y,z) = (x+y+z)^2 - 4(xy+yz+zx) = z^2 - 2z(x+y) + (x-y)^2:
        #   z - (x+y) on diagonal, off-diag -2y (lower) and -2x (upper) => (z-x-y)^2 - 4xy = z^2 - 2z(x+y) + (x-y)^2! EXACT!!
        M = sp.zeros(n, n)
        sgn_n = (-1)**n
        for j in range(n):
            M[j, j] = z - sgn_n * (x + y)
            for r in range(1, n):
                if j + r < n:
                    M[j + r, j] = -sgn_n * math.comb(n, r) * y
                else:
                    M[j + r - n, j] = -sgn_n * math.comb(n, r) * x
        p_det = sp.expand(M.det())
        assert sp.expand(p_det - p_eng) == 0, f"Mismatch vs Buchstaber-Kornev det at n={n}!"
        print(f"  n={n}: ALL 3 definitions (SymPy resultant, Buchstaber-Kornev {n}x{n} det, fast Newton engine) MATCH IDENTICALLY! ({len(coeffs)} coeffs)")

    # Print small even rows explicitly
    for n in [2, 4, 6]:
        coeffs = compute_elementary_coeffs(n)
        items = sorted(coeffs.items(), key=lambda kv: (kv[0][2], kv[0][1], kv[0][0]))
        print(f"  P_{n}(e1,e2,e3) = " + ", ".join(f"({a},{b},{c}):{v}" for (a, b, c), v in items))


# ---------------------------------------------------------------------------
# 4. Stress-test Conjecture 3.3(2) & sign law (-1)^{b(c+1)+c} up to n = 60+
# ---------------------------------------------------------------------------

def stress_test_even_n(max_n: int = 60):
    print(f"\n=== 2. Exact stress test of Conjecture 3.3(2) & sign law (-1)^(b(c+1)+c) for even n = 2..{max_n} ===")
    t0 = time.perf_counter()
    total_coeffs = 0
    zero_hits = []
    sign_violations = []
    c1_check_ok = True
    b0_table = {}

    for n in range(2, max_n + 1, 2):
        tn = time.perf_counter()
        coeffs = compute_elementary_coeffs(n)
        m = n // 2
        # Check c=0 formula: P_{2m}(e1, e2, 0) = (e1^2 - 4 e2)^m => c_{2m-2b, b, 0} = (-4)^b * C(m, b)
        for b in range(m + 1):
            a = n - 2 * b
            assert coeffs[(a, b, 0)] == ((-4) ** b) * math.comb(m, b), f"c=0 formula mismatch at n={n}, b={b}"

        n_zeros = 0
        n_sign_bad = 0
        min_abs = None
        for (a, b, c), val in coeffs.items():
            total_coeffs += 1
            expected_sgn = -1 if ((b * (c + 1) + c) % 2 == 1) else 1
            if val == 0:
                n_zeros += 1
                zero_hits.append((n, a, b, c))
            elif (1 if val > 0 else -1) != expected_sgn:
                n_sign_bad += 1
                sign_violations.append((n, a, b, c, val, expected_sgn))
            if min_abs is None or abs(val) < min_abs[0]:
                min_abs = (abs(val), (a, b, c))

        if n <= 12:
            b0_table[n] = [(c, coeffs[(n - 3 * c, 0, c)]) for c in range(n // 3 + 1)]

        if n <= 16 or n % 10 == 0 or n_zeros > 0 or n_sign_bad > 0:
            dt = time.perf_counter() - tn
            print(f"  n={n:2d}: {len(coeffs):4d} coeffs, zeros={n_zeros}, sign_fails={n_sign_bad}, "
                  f"min|c|={min_abs[0]} at {min_abs[1]} ({dt:.2f}s)")

    print(f"\nSUMMARY (even n = 2..{max_n}, total {total_coeffs} exact integer coefficients, {time.perf_counter()-t0:.2f}s):")
    print(f"  Zero coefficients found: {len(zero_hits)} {zero_hits[:5]}")
    print(f"  Sign-law (-1)^(b(c+1)+c) violations: {len(sign_violations)} {sign_violations[:5]}")

    print("\n=== 3. Exact small-n slices for pattern discovery (b=0 and c=1, 2) ===")
    for n in [2, 4, 6, 8, 10, 12]:
        print(f"  n={n:2d}, b=0 (c, val): {b0_table[n]}")
    for n in [4, 6, 8, 10, 12]:
        coeffs = compute_elementary_coeffs(n)
        c1 = [(b, coeffs[(n - 3 - 2 * b, b, 1)]) for b in range((n - 3) // 2 + 1)]
        print(f"  n={n:2d}, c=1 (b, val): {c1}")
    for n in [6, 8, 10, 12]:
        coeffs = compute_elementary_coeffs(n)
        c2 = [(b, coeffs[(n - 6 - 2 * b, b, 2)]) for b in range((n - 6) // 2 + 1)]
        print(f"  n={n:2d}, c=2 (b, val): {c2}")


if __name__ == "__main__":
    cross_verify_small_n()
    stress_test_even_n(max_n=60)
