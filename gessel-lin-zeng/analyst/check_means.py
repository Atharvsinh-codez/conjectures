"""Super-fast exact root/mean/interlacing probe in x-variable up to k=22 (using exact A_{k,i}(x) and float/mpmath)."""
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
    """S_L p(x) = (1 + L x) p(x) + x(1-x) p'(x), low-to-high."""
    d = len(p) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(p):
        # c * [(1 + m) x^m + (L - m) x^{m+1}]
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

def mean_x(p):
    """mu(p) = p'(1)/p(1)."""
    s0 = sum(p)
    s1 = sum(m * c for m, c in enumerate(p))
    return s1 / s0

A = build_A(25)
print("=== CHECK MEAN DIFFERENCES AT k = 25 ===")
k = 25
for i in range(2, k + 1, 2):
    L = 3*k - 2 - i
    m_A = mean_x(A[k-1][i])
    m_SA = mean_x(apply_S(L - 1, A[k-1][i]))
    m_Am1 = mean_x(A[k-1][i-1])
    m_Am2 = mean_x(A[k-1][i-2])
    m_SSA = mean_x(apply_S(L, apply_S(L - 1, A[k-1][i])))
    print(f"k={k}, i={i:2d}: mu(A_{{i-2}})-mu(S A_i) = {m_Am2 - m_SA:+.5f}, "
          f"mu(S A_i)-mu(A_{{i-1}}) = {m_SA - m_Am1:+.5f}, "
          f"mu(S^2 A_i)-mu(A_{{i-2}}) = {m_SSA - m_Am2:+.5f}")
