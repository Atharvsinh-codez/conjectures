import numpy as np

# Test the EXACT GLZ Stage 2+3 maps on d=2 and d=3 polynomials!
# 1. Exact i=1 map: (P0, P1, P2) -> (D_j^2 P0, D_{j+1} D_j P1, D_{j+1}^2 P2)
#    Wait: recall P0, P1, P2 are equal-degree polynomials with P0 << P1 << P2, J(P0, P1, P2) >= 1 + eps!
#    Wait: does (P0, P1, P2) -> (D_j^2 P0, D_{j+1} D_j P1, D_{j+1}^2 P2) preserve J >= 1 (even WITHOUT pole gap)?!
# 2. Exact i=2 map WITH final D_{j+1}: (U0, D_j U1, D_j^2 U2) -> (D_{j+1} D_j U0, D_{j+1}^2 D_j U1, D_{j+1}^3 D_j U2)
# 3. Also test the double-operator map: does D_m(D_m + y) have ANY counterexample for wild m_j,
#    OR does T_j(T_j + w_j y) fail at k=8..10 like single T_j + w_j y did?!

def apply_D(p, m):
    # p in ascending powers of t
    d = len(p) - 1
    res = np.zeros(d + 2, dtype=np.float64)
    for n in range(d + 1):
        res[n] += (1 + n) * p[n]
        res[n + 1] += (m + n) * p[n]
    return res

def wronskian(f, g):
    # W(f, g) = f g' - f' g
    fp = np.polyder(f[::-1])[::-1] if len(f) > 1 else np.array([0.0])
    gp = np.polyder(g[::-1])[::-1] if len(g) > 1 else np.array([0.0])
    a = np.convolve(f, gp)
    b = np.convolve(fp, g)
    n = max(len(a), len(b))
    res = np.zeros(n)
    res[:len(a)] += a
    res[:len(b)] -= b
    return res

def eval_poly(p, t):
    # Horner on ascending coeffs
    v = np.zeros_like(t, dtype=np.float64)
    for c in p[::-1]:
        v = v * t + c
    return v

def min_J_on_R(X, M, Y, ts):
    wxm = eval_poly(wronskian(X, M), ts)
    wmy = eval_poly(wronskian(M, Y), ts)
    wxy = eval_poly(wronskian(X, Y), ts)
    if np.min(wxm) <= 0 or np.min(wmy) <= 0:
        return -1.0
    # J = 4 W_XM W_MY / W_XY^2
    # avoid 0/0: where wxy == 0, J = +inf
    with np.errstate(divide='ignore', invalid='ignore'):
        J = 4.0 * wxm * wmy / (wxy * wxy)
    J = np.where(np.abs(wxy) < 1e-14 * np.sqrt(wxm * wmy), 1e12, J)
    return float(np.min(J))

# Dense grid on R (especially (-1, 0) and near 0, -1, +-inf)
ts = np.concatenate([
    -np.geomspace(1e-4, 1e3, 120),
    np.linspace(-0.9999, -0.0001, 400),
    -1.0 - np.geomspace(1e-4, 1e3, 120),
    np.geomspace(1e-4, 1e3, 120),
])

def from_roots(roots):
    # monic with positive constant term? No: prod(t - r_i) for r_i in (-1, 0) has ALL positive coeffs!
    p = np.array([1.0])
    for r in roots:
        p = np.convolve(p, [-r, 1.0])
    return p

print("=== 1. EXACT i=1 MAP: (P0, P1, P2) -> (D_j^2 P0, D_{j+1} D_j P1, D_{j+1}^2 P2) ===")
rng = np.random.default_rng(42)
for d in [2, 3]:
    for j in [5, 10, 20, 40, 80]:
        valid = 0
        min_Jout = 1e9
        max_drop = -1e9
        worst = None
        for trial in range(4000):
            # Generate interlacing roots in (-1, 0):
            # P0 << P1 << P2 means r_k(P2) < r_k(P1) < r_k(P0) or vice versa?
            # Check W(t+a, t+b) = (t+a)*1 - 1*(t+b) = a - b.
            # So W(X, M) > 0 requires -r(X) > -r(M) iff r(X) < r(M) (roots of X to the LEFT of roots of M!).
            # Wait: check GLZ! h_{1,0} = 1+2t (root -1/2), D_3 h_{1,1} = 1+3t (root -1/3):
            # X = h_{1,0} (root -1/2), M = D_3 h_{1,1} (root -1/3): r(X) = -1/2 < r(M) = -1/3!
            # W(1+2t, 1+3t) = (1+2t)*3 - 2*(1+3t) = 1 > 0! YES!
            # And in GLZ: P0 = h_{k-1,0} has roots to the LEFT of P1, which has roots to the LEFT of P2!
            # So X = P0, M = P1, Y = P2 has r_1(P0) < r_1(P1) < r_1(P2) < r_2(P0) < ... < 0!
            pts = np.sort(-rng.uniform(0.01, 0.99, size=3*d))
            r0 = pts[0::3]
            r1 = pts[1::3]
            r2 = pts[2::3]
            P0 = from_roots(r0)
            P1 = from_roots(r1) * np.sqrt(1.0 + 4.0 / j)
            P2 = from_roots(r2)
            # Check J_in
            Jin = min_J_on_R(P0, P1, P2, ts)
            if Jin < 1.0001:
                continue
            valid += 1
            Q0 = apply_D(apply_D(P0, j), j)
            Q1 = apply_D(apply_D(P1, j), j + 1)
            Q2 = apply_D(apply_D(P2, j + 1), j + 1)
            Jout = min_J_on_R(Q0, Q1, Q2, ts)
            drop = Jin - Jout
            if Jout < min_Jout:
                min_Jout = Jout
            if drop > max_drop:
                max_drop = drop
                worst = (Jin, Jout, r0, r1, r2)
        print(f"  d={d}, j={j:2d} [n={valid:4d}]: min_Jout={min_Jout:8.5f}, max_drop={max_drop:+8.5f} (j*d={j*max_drop:+7.4f}, j^2*d={j*j*max_drop:+8.3f})")
        if worst and min_Jout < 1.0:
            print(f"     FAIL BELOW 1! Jin={worst[0]:.5f}, Jout={worst[1]:.5f}")

print("\n=== 2. CHECK REVERSE ORDER TOO: (D_{j+1}^2 P0, D_{j+1} D_j P1, D_j^2 P2) ===")
for j in [5, 10, 20, 40]:
    valid = 0
    min_Jout = 1e9
    max_drop = -1e9
    for trial in range(3000):
        pts = np.sort(-rng.uniform(0.01, 0.99, size=6))
        P0, P1, P2 = from_roots(pts[0::3]), from_roots(pts[1::3]) * np.sqrt(1.0 + 4.0 / j), from_roots(pts[2::3])
        Jin = min_J_on_R(P0, P1, P2, ts)
        if Jin < 1.0001: continue
        valid += 1
        Q0 = apply_D(apply_D(P0, j+1), j+1)
        Q1 = apply_D(apply_D(P1, j), j+1)
        Q2 = apply_D(apply_D(P2, j), j)
        Jout = min_J_on_R(Q0, Q1, Q2, ts)
        min_Jout = min(min_Jout, Jout)
        max_drop = max(max_drop, Jin - Jout)
    print(f"  REV d=2, j={j:2d} [n={valid:4d}]: min_Jout={min_Jout:8.5f}, max_drop={max_drop:+8.5f} (j*d={j*max_drop:+7.4f})")

print("\n=== 3. SINGLE STEP: (D_j P0, D_{j+1} P1, D_{j+2} P2) vs (D_j P0, D_j P1, D_{j+1} P2)? No: check (D_j, D_{j+1}, D_{j+2}) on same (P0, P1, P2) ===")
for j in [5, 10, 20, 40]:
    valid = 0
    min_Jout = 1e9
    max_drop = -1e9
    for trial in range(3000):
        pts = np.sort(-rng.uniform(0.01, 0.99, size=6))
        P0, P1, P2 = from_roots(pts[0::3]), from_roots(pts[1::3]) * np.sqrt(1.0 + 4.0 / j), from_roots(pts[2::3])
        Jin = min_J_on_R(P0, P1, P2, ts)
        if Jin < 1.0001: continue
        valid += 1
        Q0 = apply_D(P0, j)
        Q1 = apply_D(P1, j+1)
        Q2 = apply_D(P2, j+2)
        Jout = min_J_on_R(Q0, Q1, Q2, ts)
        min_Jout = min(min_Jout, Jout)
        max_drop = max(max_drop, Jin - Jout)
    print(f"  L_j d=2, j={j:2d} [n={valid:4d}]: min_Jout={min_Jout:8.5f}, max_drop={max_drop:+8.5f}")
