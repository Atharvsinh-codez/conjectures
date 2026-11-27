"""EXACT INTEGER / RATIONAL CHECK OF TOTAL NONNEGATIVITY (ZERO FLOAT ROUNDOFF!):
1. Is C_{N, i} = [y^i] prod_{j=1}^N (A_j + lam_j y I) EXACTLY TOTALLY NONNEGATIVE
   (all 2x2, 3x3, 4x4 minors >= 0) for:
   (a) GLZ k = 2, 3, 4, 5 (in t-space AND in x-space!)
   (b) Wild Laguerre L_{c_j, m_j} (N = 4, 5, 6, exact rational c_j, m_j, lam_j!)
   (c) Wild Atomic D_{c_j, m_j} (N = 4, 5, 6!)
   (d) Arbitrary positive lower bidiagonal A_j (just to see if ALL bidiag fail TN for i>=1!)
2. Is the coefficient matrix M_k with entries (M_k)_{r, d} = [x^d] A_{k, k-r}(x)
   (and [t^d] htilde_{k, k-r}(t)) EXACTLY TOTALLY NONNEGATIVE (all 2x2, 3x3, 4x4 minors >= 0)
   for GLZ k = 2..8?!
3. ALSO check the TOEPLITZ matrix of coefficients of A_{k, i}(x) (Wait: PF sequence <=> Toeplitz matrix is TN!!),
   and check the 2-row interlacing matrix [[A_{k,i}], [A_{k,i-1}]]!
"""
import itertools
import sympy as sp
from fractions import Fraction

def det_exact(mat):
    """Exact Bareiss determinant for small matrix of ints or Fractions."""
    n = len(mat)
    if n == 1: return mat[0][0]
    if n == 2: return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]
    if n == 3:
        a,b,c = mat[0]; d,e,f = mat[1]; g,h,i = mat[2]
        return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
    # Use sympy for 4x4+
    return sp.Matrix(mat).det()

def check_TN_exact(M, max_order=4):
    rows = len(M); cols = len(M[0])
    res = {}
    for order in range(2, min(max_order, rows, cols) + 1):
        neg = None
        for r in itertools.combinations(range(rows), order):
            for c in itertools.combinations(range(cols), order):
                sub = [[M[ri][ci] for ci in c] for ri in r]
                d = det_exact(sub)
                if d < 0:
                    neg = (r, c, d)
                    break
            if neg is not None:
                break
        res[order] = (neg is None, neg)
    return res

def mat_mul(A, B):
    n = len(A)
    C = [[0]*n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            aik = A[i][k]
            if aik == 0: continue
            for j in range(n):
                if B[k][j] != 0:
                    C[i][j] += aik * B[k][j]
    return C

def mat_prod_y_exact(As, lams):
    N = len(As); n = len(As[0])
    C = [[[1 if r == c else 0 for c in range(n)] for r in range(n)]]
    for j in range(N):
        Aj = As[j]; lj = lams[j]
        newC = [[[0]*n for _ in range(n)] for _ in range(len(C) + 1)]
        for i, M in enumerate(C):
            AM = mat_mul(Aj, M)
            for r in range(n):
                for c in range(n):
                    newC[i][r][c] += AM[r][c]
                    if lj != 0:
                        newC[i+1][r][c] += lj * M[r][c]
        C = newC
    return C

print("=== 1(a). EXACT TN CHECK OF C_{N, i} FOR GLZ k=4 (t-space) ===")
k = 4; N = 2*k; n = 7
As = []; lams = []
for j in range(1, k + 1):
    m = j + 1
    A = [[0]*n for _ in range(n)]
    for r in range(n): A[r][r] = 1 + r
    for r in range(n-1): A[r+1][r] = m + r
    As.append(A); lams.append(1)
    As.append(A); lams.append(0)
C_glz = mat_prod_y_exact(As, lams)
for i in range(k + 1):
    res = check_TN_exact(C_glz[i], max_order=4)
    print(f"  GLZ k=4, i={i}: " + ", ".join(f"{o}x{o}:{'OK' if ok else f'NEG({info[0]},{info[1]})'}" for o,(ok,info) in res.items()))

print("\n=== 1(b). EXACT TN CHECK OF C_{N, i} FOR WILD LAGUERRE N=5 ===")
params = [(1, 1), (2, 50), (1, 100), (10, 3), (3, 7)]
N = len(params); n = N + 1
As_L = []
for c, m in params:
    A = [[0]*n for _ in range(n)]
    for r in range(n): A[r][r] = c + r
    for r in range(n-1): A[r+1][r] = m
    As_L.append(A)
C_L = mat_prod_y_exact(As_L, [1, 5, 2, 10, 3])
for i in range(N + 1):
    res = check_TN_exact(C_L[i], max_order=4)
    print(f"  Wild Lag N=5, i={i}: " + ", ".join(f"{o}x{o}:{'OK' if ok else f'NEG({info[0]},{info[1]},val={info[2]})'}" for o,(ok,info) in res.items()))

print("\n=== 1(c). EXACT TN CHECK OF C_{N, i} FOR WILD ATOMIC D_{c,m} N=5 ===")
params_D = [(1, 2), (3, 100), (1, 50), (5, 7), (2, 1000)]
N = len(params_D); n = N + 1
As_D = []
for c, m in params_D:
    A = [[0]*n for _ in range(n)]
    for r in range(n): A[r][r] = c + r
    for r in range(n-1): A[r+1][r] = m + r
    As_D.append(A)
C_D = mat_prod_y_exact(As_D, [1, 7, 2, 11, 3])
for i in range(N + 1):
    res = check_TN_exact(C_D[i], max_order=4)
    print(f"  Wild Atomic D N=5, i={i}: " + ", ".join(f"{o}x{o}:{'OK' if ok else f'NEG({info[0]},{info[1]},val={info[2]})'}" for o,(ok,info) in res.items()))

print("\n=== 1(d). EXACT TN CHECK OF C_{N, i} FOR ARBITRARY POSITIVE BIDIAGONAL N=4 ===")
import random
random.seed(7)
arb_tn_fails = {2: 0, 3: 0}
for trial in range(20):
    N = 4; n = 5
    As_a = []
    for _ in range(N):
        A = [[0]*n for _ in range(n)]
        for r in range(n): A[r][r] = random.randint(1, 20)
        for r in range(n-1): A[r+1][r] = random.randint(1, 20)
        As_a.append(A)
    C_a = mat_prod_y_exact(As_a, [1]*N)
    for o in [2, 3]:
        for i in range(1, N):
            res = check_TN_exact(C_a[i], max_order=o)
            if not res[o][0]:
                arb_tn_fails[o] += 1; break
print(f"  Arbitrary positive bidiagonal TN fails (20 trials): {arb_tn_fails}")

print("\n=== 2. EXACT TN CHECK OF COEFFICIENT MATRIX M_k = ([x^d] A_{k, k-r}(x)) AND ([t^d] htilde_{k, k-r}(t)) ===")
def apply_S(L, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r+1); res[r+1] += c*(L-r)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_D(m, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

A = [[None]*(kk+1) for kk in range(9)]
A[0][0] = [1]
for kk in range(1, 9):
    for ii in range(kk+1):
        u = [0]
        if ii <= kk-1: u = poly_add(u, apply_S(3*kk-3-ii, A[kk-1][ii]))
        if ii >= 1: u = poly_add(u, A[kk-1][ii-1])
        A[kk][ii] = apply_S(3*kk-2-ii, u)

H = [[None]*(kk+1) for kk in range(9)]
H[1][0] = [1, 2]; H[1][1] = [1]
for kk in range(2, 9):
    m = kk + 1
    U = [None]*(kk+1)
    U[0] = apply_D(m, H[kk-1][0])
    for ii in range(1, kk): U[ii] = poly_add(apply_D(m, H[kk-1][ii]), H[kk-1][ii-1])
    U[kk] = H[kk-1][kk-1]
    for ii in range(kk+1): H[kk][ii] = apply_D(m, U[ii])

for kk in range(2, 8):
    Mx = [[A[kk][kk-r][d] if d < len(A[kk][kk-r]) else 0 for d in range(2*kk)] for r in range(kk+1)]
    Mt = [[H[kk][kk-r][d] if d < len(H[kk][kk-r]) else 0 for d in range(2*kk)] for r in range(kk+1)]
    rx = check_TN_exact(Mx, max_order=4 if kk <= 5 else 3)
    rt = check_TN_exact(Mt, max_order=4 if kk <= 5 else 3)
    print(f"  k={kk} A_{{k,i}}(x) matrix: " + ", ".join(f"{o}x{o}:{'OK' if ok else f'NEG({info[0]},{info[1]})'}" for o,(ok,info) in rx.items()))
    print(f"  k={kk} htilde(t)  matrix: " + ", ".join(f"{o}x{o}:{'OK' if ok else f'NEG({info[0]},{info[1]})'}" for o,(ok,info) in rt.items()))
