"""
3-second check:
1. Check if Ma-Wang Prop 1 (A_{k, k-1} = (1/2) P_{2k} A_{k+1, k+1}) generalizes to:
   - A_{k, k-2} as a linear combination of P_c P_{c'} A_{k+2, k+2} and A_{k+1, k+1} (or x^r d^r A_{k+2, k+2})!
   - Also check in the t-variable: recall E_n(x) = A_{n, n}(x) <-> h~_{n, n}(t) = D_{n+1} D_n ... D_3 (1)!
   -Wait: what IS A_{k, k-1} = (1/2) P_{2k} A_{k+1, k+1} in the t-variable?!
    Check! deg A_{k+1, k+1} = k, c = 2k, so P_{2k} <-> (2k + k t - t(1+t) d/dt) = (2k+1)(1+t) - D_{k+1}!
    Wait: h~_{k, k-1} = (1/2) [ (2k+1)(1+t) h~_{k+1, k+1} - D_{k+1} h~_{k+1, k+1} ]??
    LET'S CHECK! Remember D_{k+2} h~_{k+1, k+1} = h~_{k+2, k+2}!
    AND D_{k+1} = D_{k+2} - t!!
    Wait: check if h~_{k, k-1} is a simple combo of h~_{k+1, k+1} and h~_{k+2, k+2}!
2. Check AGGREGATE A_{k+2, k+2} << A_{k, k-3} through k = 35!
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 80

def apply_SL_int(p, L):
    out = [0]*(len(p)+1)
    for m, c in enumerate(p):
        out[m] += (1+m)*c; out[m+1] += (L-m)*c
    while len(out)>1 and out[-1]==0: out.pop()
    return out

def padd_int(a, b):
    n = max(len(a), len(b))
    out = [0]*n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return out

def apply_Pc_int(p, c):
    return [(c - m)*v for m, v in enumerate(p)]

A = {(1,0): [1, 1], (1,1): [1]}
for k in range(2, 12):
    for i in range(k + 1):
        L = 3*k - 2 - i
        u = [0]
        if i <= k - 1: u = padd_int(u, apply_SL_int(A[(k-1, i)], L - 1))
        if i >= 1: u = padd_int(u, A[(k-1, i-1)])
        A[(k, i)] = apply_SL_int(u, L)

# 1. Verify Ma-Wang Prop 1: 2 * A_{k, k-1} == P_{2k} A_{k+1, k+1}
for k in range(1, 8):
    lhs = [2*c for c in A[(k, k-1)]]
    rhs = apply_Pc_int(A[(k+1, k+1)], 2*k)
    assert lhs == rhs, f"Mismatch at k={k}: {lhs} vs {rhs}"
print("Verified Ma-Wang Prop 1: 2 A_{k, k-1} = P_{2k} A_{k+1, k+1} for k=1..7!")

# Check if A_{k, k-2} is a linear combo of {A_{k+2,k+2}, x A'_{k+2,k+2}, x^2 A''_{k+2,k+2}, A_{k+1,k+1}, x A'_{k+1,k+1}}!
print("\n=== Test linear relation for A_{k, k-2} in terms of {theta^r E_{k+2}, theta^r E_{k+1}} ===")
for k in [3, 4, 5, 6]:
    E2 = A[(k+2, k+2)]  # deg k+1
    E1 = A[(k+1, k+1)] + [0]  # pad to deg k+1
    # basis: E2, th E2, th^2 E2, E1, th E1, x E1
    def th(p): return [m*v for m, v in enumerate(p)]
    xE1 = [0] + A[(k+1, k+1)]
    basis = [E2, th(E2), th(th(E2)), E1, th(E1), xE1]
    # Solve M * c = A[(k, k-2)] over Q!
    d = len(E2)
    M = sp.Matrix([[basis[j][m] for j in range(len(basis))] for m in range(d)])
    rhs = sp.Matrix([A[(k, k-2)][m] for m in range(d)])
    # Check if solvable!
    aug = M.row_join(rhs)
    rM = M.rank(); rA = aug.rank()
    if rM == rA:
        sol = M.LUsolve(rhs) if M.rows == M.cols else M.gauss_jordan_solve(rhs)
        print(f"  k={k}: EXACT SOLUTION FOUND! coeffs = {sol}")
    else:
        print(f"  k={k}: no exact combo in 6-dim basis (rank {rM} vs {rA})")

# Wait: what about in terms of E_{k+2}, E_{k+1}, E_k (3 consecutive Eulerian polys, with theta)?!
# Notice A_{k, k-1} = (1/2)(2k - theta) E_{k+1}!
# Compute the residual A_{k, k-2} - c * P_a P_b E_{k+2}:
# Recall earlier: P_{2k+1} P_{2k+2} E_{k+2} = S_{2k}(2 P_{2k+1} A_{k, k-1} + 2 P_{2k+1} E_{k+1})!
# Wait: what is 2 A_{k, k-1} = P_{2k} E_{k+1}, so P_{2k+1} E_{k+1} = 2 A_{k, k-1} + E_{k+1}!
# Let's check the EXACT formula for A_{k, k-2} in terms of {theta^r E_{k+2}, x^a theta^b E_{k+1}}!
for k in [4, 5, 6, 7]:
    E2 = A[(k+2, k+2)]  # deg k+1
    # Try ONLY {E2, th E2, th^2 E2} vs A_{k, k-2} AND {E_{k+2}, th E_{k+2}, th^2 E_{k+2}, E_{k+1}, x E_{k+1}, x^2 E_k}
    def th(p): return [m*v for m, v in enumerate(p)]
    def mulx(p, pad_len):
        q = [0] + list(p)
        return q + [0]*(pad_len - len(q))
    d = len(E2)
    def pad(p): return list(p) + [0]*(d - len(p))
    # Try polynomial in (x, theta) acting on E2: wait, S_{2k+2} E_{k+1} = E_{k+2}, so E_{k+1} is NOT needed if we include x? No, E_{k+1} is inverse S!
    basis = [E2, th(E2), th(th(E2))]
    M = sp.Matrix([[basis[j][m] for j in range(3)] for m in range(d)])
    rhs = sp.Matrix([A[(k, k-2)][m] for m in range(d)])
    aug = M.row_join(rhs)
    print(f"  k={k}: pure theta-poly of E_{{k+2}}? rank {M.rank()} vs {aug.rank()}")

# 2. Check AGGREGATE A_{k+2, k+2} << A_{k, k-3} through k = 30!
print("\n=== 2. AGGREGATE A_{k+2, k+2} << A_{k, k-3} THROUGH k = 30 ===")
def apply_Dj_mp(poly, j):
    out = [mp.mpf(0)] * (len(poly) + 1)
    for m, c in enumerate(poly):
        if c != 0:
            out[m] += (1 + m) * c; out[m + 1] += (j + m) * c
    return out

def padd_mp(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for idx in range(len(a)): out[idx] += a[idx]
    for idx in range(len(b)): out[idx] += b[idx]
    return out

def check_interlace_bez(f, g):
    n = len(g) - 1
    assert len(f) == n
    s = (abs(g[0]) / abs(g[-1])) ** (mp.mpf(1) / n)
    gs = [g[m] * (s ** m) for m in range(n + 1)]
    fs = [f[m] * (s ** m) for m in range(n)] + [mp.mpf(0)]
    mxg = max(abs(x) for x in gs); mxf = max(abs(x) for x in fs)
    gs = [x / mxg for x in gs]; fs = [x / mxf for x in fs]
    H = [[mp.mpf(0)] * n for _ in range(n)]
    for r in range(n + 1):
        for col in range(r + 1, n + 1):
            val = -(gs[r]*fs[col] - gs[col]*fs[r])
            for m in range(col - r):
                H[r + m][r + (col - r - 1 - m)] += val
    for i in range(n):
        piv = H[i][i]
        if piv <= 0: return False, i
        for j in range(i + 1, n):
            fac = H[j][i] / piv
            for k in range(j, n):
                H[j][k] -= fac * H[i][k]
                H[k][j] = H[j][k]
    return True, n

h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
h_store = {1: h}
for k in range(2, 33):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = padd_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = padd_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
    h_store[k] = h

for k in [4, 6, 8, 10, 15, 20, 25, 30]:
    ok, idx = check_interlace_bez(h_store[k+2][k+2], h_store[k][k-3])
    print(f"  k={k:2d}: agg h~_{{k+2, k+2}} << h~_{{k, k-3}}: {'PASS' if ok else f'FAIL@{idx}'}")
