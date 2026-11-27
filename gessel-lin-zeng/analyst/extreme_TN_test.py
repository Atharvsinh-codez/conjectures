"""100-trial extreme adversarial test of:
1. Exact TN of C_{N, i} = [y^i] prod_{j=1}^N (A_j + lam_j y I) for EXTREME wild Laguerre & Atomic D!
2. Exact TN of M_N = ([t^d] F_{N, N-r}(t)) for EXTREME wild Laguerre & Atomic D, AND GLZ k = 2..8!
3. ALSO check the EXACT Toeplitz-concatenation / Hurwitz / Hankel / Wronskian positivity for Wild Laguerre & Atomic!
"""
import itertools
import random
import sympy as sp

def det_exact(mat):
    n = len(mat)
    if n == 1: return mat[0][0]
    if n == 2: return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]
    if n == 3:
        a,b,c = mat[0]; d,e,f = mat[1]; g,h,i = mat[2]
        return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
    return sp.Matrix(mat).det()

def check_TN_exact(M, max_order=4):
    rows = len(M); cols = len(M[0])
    for order in range(2, min(max_order, rows, cols) + 1):
        for r in itertools.combinations(range(rows), order):
            for c in itertools.combinations(range(cols), order):
                sub = [[M[ri][ci] for ci in c] for ri in r]
                d = det_exact(sub)
                if d < 0:
                    return False, (order, r, c, d)
    return True, None

def mat_mul(A, B):
    n = len(A)
    C = [[0]*n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            aik = A[i][k]
            if aik == 0: continue
            for j in range(n):
                if B[k][j] != 0: C[i][j] += aik * B[k][j]
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
                    if lj != 0: newC[i+1][r][c] += lj * M[r][c]
        C = newC
    return C

rng = random.Random(2026)

for kind in ["EXTREME Laguerre", "EXTREME Atomic D"]:
    C_fails = 0; M_fails = 0
    first_C_fail = None; first_M_fail = None
    for trial in range(50):
        N = 5; n = N + 1
        As = []; lams = []
        for j in range(N):
            c = rng.choice([1, 2, 5, 20, 100, 1000])
            m = (c + rng.choice([1, 2, 10, 100, 10000])) if "Atomic" in kind else rng.choice([1, 2, 10, 100, 10000])
            lam = rng.choice([1, 2, 10, 100, 1000])
            A = [[0]*n for _ in range(n)]
            for r in range(n): A[r][r] = c + r
            for r in range(n-1): A[r+1][r] = (m + r) if "Atomic" in kind else m
            As.append(A); lams.append(lam)
        C = mat_prod_y_exact(As, lams)
        # Check C_{N, i} TN up to 4x4
        for i in range(1, N):
            ok, info = check_TN_exact(C[i], max_order=4)
            if not ok:
                C_fails += 1
                if first_C_fail is None: first_C_fail = (trial, i, info[:3])
                break
        # Check M_N = ([t^d] F_{N, N-r}(t)) TN up to 4x4
        MN = [[C[N - r][d][0] for d in range(n)] for r in range(N + 1)]
        okM, infoM = check_TN_exact(MN, max_order=4)
        if not okM:
            M_fails += 1
            if first_M_fail is None: first_M_fail = (trial, infoM[:3])
    print(f"{kind:18s} (50 extreme trials, N=5, up to 4x4): C_TN_fails={C_fails} {first_C_fail}, M_TN_fails={M_fails} {first_M_fail}", flush=True)
