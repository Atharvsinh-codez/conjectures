"""PROBE THE LOWER BIDIAGONAL TRANSFER MATRIX REPRESENTATION:
  F_N(t, y) = e(t)^T [ prod_{j=1}^N (A_j + lam_j y I) ] e_0
where (A_j)_{r,r} = r + c_j, (A_j)_{r+1,r} = r + m_j (Atomic D) or m_j (Laguerre L).

Tests:
1. For ARBITRARY positive lower bidiagonal matrices A_j (random a_{j,r}, b_{j,r} > 0):
   Does [y^i] e(t)^T prod (A_j + y I) e_0 ALWAYS have all real roots in t?
   (If YES -> huge linear-algebra theorem! If NO -> tells us the EXACT condition on a_{j,r}, b_{j,r}!)
2. What if a_{j,r} = c_j + r, b_{j,r} = m_j + alpha * r (or arbitrary functions of r)? Which class preserves RZ?
3. Is the matrix C_{N, i} = [y^i] prod_{j=1}^N (A_j + lam_j y I) TOTALLY NONNEGATIVE (all 2x2 and 3x3 minors >= 0)
   for GLZ, Atomic D, and Laguerre L?!
4. Is the matrix M_{N} with entries (M_N)_{i, d} = [t^d] F_{N, i}(t) (or [x^d] A_{k,i}(x)) TOTALLY NONNEGATIVE?!
   (Check 2x2 and 3x3 minors!)
5. WAIT A SECOND: What about the EIGENVALUES of a matrix whose characteristic polynomial is F_{N,i}(t)?!
   Check the companion / tridiagonal / Jacobi representation? No, first run 1-4!
"""
import numpy as np
import itertools

def mat_prod_y(As, lams):
    """Return list C[0..N] where C[i] = [y^i] prod_{j=N..1} (A_j + lam_j y I)."""
    N = len(As)
    n = As[0].shape[0]
    C = [np.eye(n)]
    for j in range(N):
        Aj = As[j]; lj = lams[j]
        newC = [np.zeros((n, n)) for _ in range(len(C) + 1)]
        for i, M in enumerate(C):
            newC[i] += Aj @ M
            newC[i + 1] += lj * M
        C = newC
    return C

def check_all_real_roots(coeffs):
    # strip trailing zeros
    c = list(coeffs)
    while len(c) > 1 and abs(c[-1]) < 1e-14:
        c.pop()
    if len(c) <= 2: return True, 0.0
    rts = np.roots(c[::-1])
    mi = float(np.max(np.abs(rts.imag)))
    return (mi < 1e-7), mi

def min_minor(M, order):
    rows, cols = M.shape
    mn = 1e300
    for r in itertools.combinations(range(rows), order):
        for c in itertools.combinations(range(cols), order):
            val = float(np.linalg.det(M[np.ix_(r, c)]))
            if val < mn:
                mn = val
    return mn

print("=== 1. ARBITRARY POSITIVE LOWER BIDIAGONAL MATRICES A_j ===")
rng = np.random.default_rng(42)
arb_fails = 0
for trial in range(100):
    N = int(rng.integers(3, 8))
    n = N + 1
    As = []
    for j in range(N):
        A = np.zeros((n, n))
        for r in range(n):
            A[r, r] = 10 ** rng.uniform(-2, 2)
        for r in range(n - 1):
            A[r + 1, r] = 10 ** rng.uniform(-2, 2)
        As.append(A)
    lams = [10 ** rng.uniform(-1, 1) for _ in range(N)]
    C = mat_prod_y(As, lams)
    bad = False
    for i in range(N + 1):
        ok, mi = check_all_real_roots(C[i][:, 0])
        if not ok:
            arb_fails += 1; bad = True
            if arb_fails <= 5:
                print(f"  ARB BIDIAG FAIL trial {trial} (N={N}, i={i}): max|Im|={mi:.3e}")
            break
print(f"  Arbitrary positive lower bidiagonal fails: {arb_fails}/100")

print("\n=== 2. SPECIALIZED CLASSES OF LOWER BIDIAGONAL MATRICES A_j ===")
# Test:
# (a) Toeplitz bidiagonal: a_{j,r} = c_j, b_{j,r} = m_j (constant in r!) -> (c_j + m_j t + lam_j y) commutes!
# (b) Independent of j in r-slope: a_{j,r} = c_j + r, b_{j,r} = m_j + r (Atomic D!)
# (c) General linear in r: a_{j,r} = c_j + u_j * r, b_{j,r} = m_j + v_j * r (positive u_j, v_j)
# (d) Quadratic in r: a_{j,r} = c_j + r^2, b_{j,r} = m_j + r^2
# (e) Decreasing in r: a_{j,r} = c_j + 1/(r+1), b_{j,r} = m_j + 1/(r+1)
# (f) Random strictly increasing in r (same for all j): a_{j,r} = c_j + f(r), b_{j,r} = m_j + g(r)
# (g) Same A_j up to diagonal shift + subdiagonal scale: a_{j,r} = c_j + a(r), b_{j,r} = m_j * b(r)
for name, make_A in [
    ("linear u_j, v_j varying", lambda j, n, r: (10**r.uniform(-1,1) + (10**r.uniform(-1,1))*np.arange(n),
                                                 10**r.uniform(-1,1) + (10**r.uniform(-1,1))*np.arange(n-1))),
    ("quadratic r^2", lambda j, n, r: (10**r.uniform(-1,1) + np.arange(n)**2,
                                       10**r.uniform(-1,1) + np.arange(n-1)**2)),
    ("exp 2^r", lambda j, n, r: (10**r.uniform(-1,1) + 2.0**np.arange(n),
                                 10**r.uniform(-1,1) + 2.0**np.arange(n-1))),
    ("decreasing 1/(r+1)", lambda j, n, r: (10**r.uniform(-1,1) + 1.0/(np.arange(n)+1),
                                            10**r.uniform(-1,1) + 1.0/(np.arange(n-1)+1))),
]:
    if make_A is None: continue
    fails = 0
    for trial in range(60):
        N = int(rng.integers(3, 8)); n = N + 1
        As = []
        for j in range(N):
            diag, sub = make_A(j, n, rng)
            A = np.diag(diag) + np.diag(sub, -1)
            As.append(A)
        lams = [10 ** rng.uniform(-1, 1) for _ in range(N)]
        C = mat_prod_y(As, lams)
        for i in range(N + 1):
            ok, mi = check_all_real_roots(C[i][:, 0])
            if not ok:
                fails += 1; break
    print(f"  {name:30s}: fails = {fails}/60")

# Handle (f) and (g): fixed random sequences a(r), b(r) across j!
for name, mode in [("fixed rand a(r), b=m_j", 1), ("fixed rand inc a(r), b=m_j", 2),
                   ("a=c_j+r, fixed rand b(r)", 3), ("a=c_j*a(r), b=m_j*b(r) fixed", 4),
                   ("Laguerre a=c_j+u_j*r, b=m_j", 5), ("Jacobi-type a=c_j+r, b=m_j-r?", 6)]:
    fails = 0
    for trial in range(60):
        N = int(rng.integers(3, 8)); n = N + 1
        ar = 10**rng.uniform(-1, 1, n)
        ar_inc = np.cumsum(10**rng.uniform(-1, 1, n))
        br = 10**rng.uniform(-1, 1, n-1)
        As = []
        for j in range(N):
            cj = 10**rng.uniform(-1, 1); mj = 10**rng.uniform(-1, 1)
            if mode == 1: diag, sub = cj + ar, np.full(n-1, mj)
            elif mode == 2: diag, sub = cj + ar_inc, np.full(n-1, mj)
            elif mode == 3: diag, sub = cj + np.arange(n), mj * br
            elif mode == 4: diag, sub = cj * ar, mj * br
            elif mode == 5: diag, sub = cj + (10**rng.uniform(-1,1))*np.arange(n), np.full(n-1, mj)
            elif mode == 6: diag, sub = cj + np.arange(n), (n + 5 + mj) - np.arange(n-1)
            As.append(np.diag(diag) + np.diag(sub, -1))
        lams = [10**rng.uniform(-1, 1) for _ in range(N)]
        C = mat_prod_y(As, lams)
        for i in range(N + 1):
            ok, mi = check_all_real_roots(C[i][:, 0])
            if not ok:
                fails += 1; break
    print(f"  {name:30s}: fails = {fails}/60")

print("\n=== 3. IS THE MATRIX C_{N, i} = [y^i] prod (A_j + y I) TOTALLY NONNEGATIVE (GLZ & LAGUERRE)? ===")
# GLZ: k = 4 (N = 8 operators: D_2, D_2+y, D_3, D_3+y, D_4, D_4+y, D_5, D_5+y)
k = 4; N = 2*k; n = N + 1
As = []; lams = []
for j in range(1, k + 1):
    m = j + 1
    A = np.diag(1.0 + np.arange(n)) + np.diag(m + np.arange(n - 1), -1)
    As.append(A); lams.append(1.0)
    As.append(A); lams.append(0.0)
C_glz = mat_prod_y(As, lams)
# Wait: only i = 0..k are nonzero since half of lams are 0!
for i in range(k + 1):
    m2 = min_minor(C_glz[i][:7, :5], 2)
    m3 = min_minor(C_glz[i][:6, :4], 3)
    print(f"  GLZ k=4, i={i}: min 2x2 minor of C_{{N,i}} = {m2:+.3e}, min 3x3 = {m3:+.3e}")

# Also test Laguerre N=5 with wild c_j, m_j!
N = 5; n = N + 1
As_L = [np.diag(c + np.arange(n)) + np.diag(np.full(n-1, m), -1)
        for c, m in [(0.5, 1.0), (2.0, 50.0), (0.1, 0.2), (10.0, 3.0), (1.5, 100.0)]]
C_L = mat_prod_y(As_L, [1.0]*N)
for i in range(N + 1):
    m2 = min_minor(C_L[i], 2)
    m3 = min_minor(C_L[i], 3)
    print(f"  Wild Lag N=5, i={i}: min 2x2 minor of C_{{N,i}} = {m2:+.3e}, min 3x3 = {m3:+.3e}")

print("\n=== 4. IS THE COEFFICIENT MATRIX M_{i, d} = [t^d] htilde_{k,i}(t) (and [x^d] A_{k,i}(x)) TN? ===")
# Build GLZ A_{k, i}(x) for k = 2..8 and check min 2x2 and 3x3 minors (reversing i -> k-i so deg increases with row!)
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

A = [[None]*(kk+1) for kk in range(9)]
A[0][0] = [1]
for kk in range(1, 9):
    for ii in range(kk+1):
        u = [0]
        if ii <= kk-1: u = poly_add(u, apply_S(3*kk-3-ii, A[kk-1][ii]))
        if ii >= 1: u = poly_add(u, A[kk-1][ii-1])
        A[kk][ii] = apply_S(3*kk-2-ii, u)

for kk in [3, 4, 5, 6, 7]:
    # Rows i = k, k-1, ..., 0 (so degree goes k-1, k, ..., 2k-1)
    M = np.zeros((kk + 1, 2*kk))
    for row_idx, ii in enumerate(range(kk, -1, -1)):
        for d, val in enumerate(A[kk][ii]):
            M[row_idx, d] = float(val)
    m2 = min_minor(M, 2)
    m3 = min_minor(M, 3) if kk <= 5 else 0.0
    print(f"  GLZ A_{{k,i}}(x) coeff matrix k={kk}: min 2x2 = {m2:+.3e}, min 3x3 = {m3:+.3e}")
