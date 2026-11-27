"""Stress test:
1. Real stability of B_k(x, y) = sum binom(k,i) A_{k,i}(x) y^i (AND binom(k,i) h_{k,i}(t) y^i)
   on general positive-direction lines (a + s, b + c*s), c > 0!
2. Also check what recurrence B_k(x, y) satisfies, AND test other weights w_{k,i} = i!, 1/i!, binom(k,i)^alpha.
"""
import math
import numpy as np

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_S(L, p):
    d = len(p) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(p):
        res[m] += c * (m + 1)
        res[m + 1] += c * (L - m)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def build_A(kmax):
    A = [[None] * (k + 1) for k in range(kmax + 1)]
    A[0][0] = [1]
    for k in range(1, kmax + 1):
        A[k][0] = apply_S(3*k - 2, apply_S(3*k - 3, A[k-1][0]))
        for i in range(1, k):
            L = 3*k - 2 - i
            A[k][i] = apply_S(L, poly_add(apply_S(L - 1, A[k-1][i]), A[k-1][i-1]))
        A[k][k] = apply_S(2*k - 2, A[k-1][k-1])
    return A

def test_bivariate_stability(polys_list, weights, name, n_trials=400):
    """Test real stability of P(x, y) = sum_i weights[i] * polys_list[i](x) * y^i
    by restricting to lines x = a + s, y = b + c*s (c > 0) AND checking Hurwitz/root reality.
    Also tests complex (x, y) with Im(x) > 0, Im(y) > 0 directly by fixing y in H+ and checking roots in x!
    WAIT: Fixing y = u + i*v with v > 0 and checking if ALL roots x of P(x, y) have Im(x) <= 0
    is an EXACT characterization of real stability (for P with deg_x(P(x, y)) constant for y in H+)!
    Wait: for B_k(x, y), the leading term in x is A_{k,0}(x) (y^0), so deg_x(B_k(x, y)) = 2k-1 is CONSTANT for ALL y!
    Therefore B_k(x, y) is real-stable IFF for EVERY y with Im(y) > 0, ALL roots x of B_k(x, y) = 0 have Im(x) <= 0!!
    THAT'S A 100x STRONGER, FASTER, DIRECT TEST!!
    """
    rng = np.random.default_rng(42)
    max_im_x = -1e9
    worst_y = None
    dmax = max(len(p) - 1 for p in polys_list)
    # Test grid + random y in upper half plane (including very small Im(y) near real axis!)
    test_ys = []
    for re_y in [-100, -20, -5, -2, -1, -0.5, -0.2, -0.05, 0, 0.05, 0.2, 0.5, 1, 2, 5, 20, 100]:
        for im_y in [1e-4, 1e-2, 0.1, 0.5, 2.0, 10.0, 100.0]:
            test_ys.append(complex(re_y, im_y))
    for _ in range(n_trials):
        re_y = rng.normal(0, 1) * (10 ** rng.uniform(-2, 2))
        im_y = 10 ** rng.uniform(-4, 2)
        test_ys.append(complex(re_y, im_y))

    for y in test_ys:
        coeffs = np.zeros(dmax + 1, dtype=np.complex128)
        for i, p in enumerate(polys_list):
            w = weights[i] * (y ** i)
            for m, c in enumerate(p):
                coeffs[m] += w * c
        # Trim near-zero leading? No, y^0 term A_{k,0} has nonzero leading coeff!
        r = np.roots(coeffs[::-1])
        m_im = float(np.max(r.imag))
        if m_im > max_im_x:
            max_im_x = m_im
            worst_y = y
    print(f"  {name}: max Im(x) over {len(test_ys)} UHP y's = {max_im_x:+.3e} (at y={worst_y})")
    return max_im_x

A = build_A(9)
print("=== DIRECT UHP STABILITY TEST FOR B_k(x, y) = sum binom(k,i) A_{k,i}(x) y^i ===")
for k in range(1, 9):
    w = [math.comb(k, i) for i in range(k + 1)]
    test_bivariate_stability(A[k], w, f"B_{k}(x,y) [binom(k,i)]")

print("\n=== ALSO TEST OTHER WEIGHTS AT k = 3, 4, 5 ===")
for k in [3, 4, 5]:
    test_bivariate_stability(A[k], [1]*(k+1), f"k={k} [w=1]")
    test_bivariate_stability(A[k], [math.factorial(i) for i in range(k+1)], f"k={k} [w=i!]")
    test_bivariate_stability(A[k], [1/math.comb(k, i) for i in range(k+1)], f"k={k} [w=1/binom]")
    test_bivariate_stability(A[k], [math.comb(k, i)**2 for i in range(k+1)], f"k={k} [w=binom^2]")
    # What about (1-x)^i * binom(k,i) * A_{k,i}(x) y^i? (i.e. in h_{k,i}(t)!)
    # Check h_{k,i}(t) with binom(k,i)
    # Wait: if B_k(x, y) is stable, let's check B_k(t/(1+t), y) * (1+t)^{2k-1} = sum binom(k,i) (1+t)^i \tilde{h}_{k,i}(t) y^i!
    # Also check sum binom(k,i) \tilde{h}_{k,i}(t) y^i directly!
