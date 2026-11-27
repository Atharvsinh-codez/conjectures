import numpy as np
import sympy as sp

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

def apply_D(m, p):
    res = [0]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

KMAX = 8
A = {(0, 0): [1]}
for k in range(1, KMAX + 1):
    for i in range(k + 1):
        u1 = apply_S(A[(k-1, i)], 3*k - 3 - i) if i <= k - 1 else [0]
        u2 = A[(k-1, i-1)] if i >= 1 else [0]
        A[(k, i)] = apply_S(poly_add(u1, u2), 3*k - 2 - i)

H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    H[k][0] = apply_D(m, apply_D(m, H[k-1][0]))
    for i in range(1, k):
        H[k][i] = apply_D(m, poly_add(apply_D(m, H[k-1][i]), H[k-1][i-1]))
    H[k][k] = apply_D(m, H[k-1][k-1])

# 1. Check Total Nonnegativity (2x2 and 3x3 minors) of the row-coefficient matrix M_k(i, j) = [x^j] A_{k, i}(x)
#    AND of M_k^{(t)}(i, j) = [t^j] tilde_h_{k, i}(t) (and reversed i)!
print("=== 1. TOTAL NONNEGATIVITY OF ROW-COEFFICIENT MATRICES M_k(i, j) ===")
def check_minors(M, name):
    rows, cols = M.shape
    neg2 = 0; neg3 = 0
    for r1 in range(rows):
        for r2 in range(r1+1, rows):
            for c1 in range(cols):
                for c2 in range(c1+1, cols):
                    det = int(M[r1, c1])*int(M[r2, c2]) - int(M[r1, c2])*int(M[r2, c1])
                    if det < 0: neg2 += 1
    if neg2 == 0 and rows >= 3 and cols >= 3:
        for r1 in range(rows):
            for r2 in range(r1+1, rows):
                for r3 in range(r2+1, rows):
                    for c1 in range(cols):
                        for c2 in range(c1+1, cols):
                            for c3 in range(c2+1, cols):
                                sub = sp.Matrix([[int(M[r, c]) for c in (c1,c2,c3)] for r in (r1,r2,r3)])
                                if sub.det() < 0: neg3 += 1
    print(f"  {name:35s}: neg_2x2={neg2:4d}, neg_3x3={neg3:4d}")

for k in [3, 4, 5]:
    # Matrix of A_{k, k-r} (increasing degree from top to bottom)
    cols = 2*k
    MA = np.zeros((k+1, cols), dtype=object)
    MH = np.zeros((k+1, cols), dtype=object)
    for r in range(k+1):
        i = k - r
        for j, v in enumerate(A[(k, i)]): MA[r, j] = v
        for j, v in enumerate(H[k][i]): MH[r, j] = v
    check_minors(MA, f"A_{{k, k-r}} (k={k})")
    check_minors(MH, f"tilde_h_{{k, k-r}} (k={k})")

# 2. Check Stieltjes continued fraction of A_{k, i}(x) / A_{k, i-1}(x) for small k!
print("\n=== 2. CHECK STIELTJES CONTINUED FRACTION OF A_{k,i}(x) / A_{k,i-1}(x) ===")
x = sp.Symbol('x')
def poly_expr(p): return sum(c * x**r for r, c in enumerate(p))

for k in [2, 3, 4]:
    for i in [1, k//2 + 1 if k//2 + 1 <= k else k]:
        num = sp.Poly(poly_expr(A[(k, i)]), x, domain='QQ')
        den = sp.Poly(poly_expr(A[(k, i-1)]), x, domain='QQ')
        # Euclidean algorithm in x (from high degree)
        qs = []
        a, b = den, num
        while b != 0 and b.degree() >= 0:
            q, r = sp.div(a, b, domain='QQ')
            qs.append(q.as_expr())
            a, b = b, -r
            if b == 0: break
        print(f"  k={k}, i={i}: quotients = {qs[:4]}")

# 3. Check if A_{k, i}(x) satisfies a LOW-ORDER differential equation or simple recurrence in i at fixed k!
#    Also check: for fixed k, consider the bivariate polynomial F_k(x, y) = sum_{i=0}^k A_{k, i}(x) y^i (or with weights):
#    Check if A_{k, i}(x) is a REAL-ROOTEDNESS-PRESERVING LINEAR TRANSFORM of A_{k, 0}(x) or Eulerian polynomials!
print("\n=== 3. CHECK EXPANSION OF A_{k, i}(x) IN EULERIAN BASIS E_d(x) OR BINOMIAL (1+x)^j x^{d-j} (GAMMA-LIKE) ===")
# Expand A_{k, i}(x) in (x + 1)^m * x^r? No, A_{k,i} is not palindromic.
# Expand A_{k, i}(x) in A_{k, k}, A_{k, k-1}, ...? No, different degrees.
# Check the roots of the RATIO of coefficients a_{k, i, j} / a_{k, i-1, j}!!
for k in [4, 6]:
    for i in [1, k//2]:
        a1 = A[(k, i)]; a0 = A[(k, i-1)]
        rats = [round(a1[j]/a0[j], 4) for j in range(len(a1))]
        print(f"  k={k}, i={i}: coeff ratios A_{{k,i}}[j]/A_{{k,i-1}}[j] = {rats}")
