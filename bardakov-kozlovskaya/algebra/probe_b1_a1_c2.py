"""
Probe algebraic proofs for the next three boundary strips: c = 2, b = 1, a = 1!
1. c = 2: Check whether Q_2(1+t, t) at t = s_j^{2m} (where -2m * arg(1 - s e^{i pi/m}) = j*pi)
   has exact alternating signs (giving m-3 roots in (0, 1/4)!) AND check dominance margin!
   Wait: at s_j where -2m*arg(1 - s e^{i pi/m}) = j*pi:
   th = pi/m, psi = j*pi/(2m), s_j = sin(psi)/sin(psi + th)!
   Range of psi: as s: 0 -> 1, -arg(1 - s e^{i pi/m}): 0 -> pi/2 - pi/(2m), so -2m*arg: 0 -> (m-1)*pi!
   Thus j = 1, 2, ..., m-2! (Wait: m-2 points -> m-3 sign changes = EXACTLY deg Q_2 = m - 3!!)
2. b = 1: Check exact formula for h_1(V) = R_1(1, V) at V = p/(1-p)^3 with p = -r^{2m} in (-1, 0):
   Verify formula against exact coeffs at n = 6, 8, 10, 12 AND check if explicit sample points r_j give all floor((2m-2)/3) negative roots!
3. a = 1: Check exact formula for S_1(W) = d/de1 P_{2m}(0, e2, e3) from P_{2m}(eps, p-1+eps, p) = prod (1 + L_{2m}(w_k) - eps)!
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 50

def poly_mul_xy(A, B):
    out = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if a:
            for j, b in enumerate(B):
                if b: out[i + j] += a * b
    return out

def compute_p_n_fast(n):
    P = [None]
    for r in range(1, n + 1):
        P.append([n * math.comb(n * r, n * j) for j in range(r + 1)])
    E = [[1]]
    for d in range(1, n + 1):
        acc = [0] * (d + 1)
        for r in range(1, d + 1):
            prod = poly_mul_xy(P[r], E[d - r])
            if (r - 1) % 2 == 0:
                for idx, v in enumerate(prod): acc[idx] += v
            else:
                for idx, v in enumerate(prod): acc[idx] -= v
        E.append([v // d for v in acc])
    def xy_to_sp(poly_xy):
        d = len(poly_xy) - 1; rem = list(poly_xy); res = [0] * (d // 2 + 1)
        for b in range(d // 2 + 1):
            c = rem[b]; res[b] = c
            if c:
                for r in range(d - 2*b + 1): rem[b + r] -= c * math.comb(d - 2*b, r)
        return res
    F = [[(-1 if d%2 else 1)*v for v in xy_to_sp(E[d])] for d in range(n, -1, -1)]
    coeffs = {}; m = n; c = 0
    while m >= 0:
        row0 = list(F[0])
        for b, val in enumerate(row0): coeffs[(m - 2*b, b, c)] = val
        for b, val in enumerate(row0):
            if not val: continue
            a = m - 2*b
            for j in range(b + 1):
                bp = b - j; vj = val * math.comb(b, j)
                for i in range(a + 1): F[i + j][bp] -= vj * math.comb(a, i)
        if m < 3: break
        m_new = m - 3
        F = [[F[k+1][bp+1] for bp in range((m_new-k)//2+1)] for k in range(m_new+1)]
        m = m_new; c += 1
    return coeffs

print("=== 1. TEST c = 2 SAMPLE POINTS s_j = sin(j pi/(2m)) / sin((j+2) pi/(2m)) FOR m = 4..16 ===")
for m in [4, 5, 6, 7, 8, 10, 12, 16]:
    n = 2 * m
    coeffs = compute_p_n_fast(n)
    deg_q2 = m - 3
    q2_coeffs = [coeffs[(n - 6 - 2*b, b, 2)] for b in range(deg_q2 + 1)]
    # Evaluate Q_2(1, u) at u = t/(1+t)^2 where t = s_j^{2m}, s_j = sin(j pi/(2m)) / sin((j+2) pi/(2m)) for j = 1..m-2
    th = mp.pi / m
    signs = []
    # Also check j = 0 (u = 0 -> Q_2(1, 0) > 0) and j = m-1? Wait: for j = m-2, psi + th = m*pi/(2m) = pi/2, s = sin((m-2)pi/(2m)) < 1!
    # Check j = 0 (u=0) plus j = 1..m-2 (wait: how many points? m-1 points -> m-2 intervals, need m-3 sign changes!)
    # Let's print the sign at u = 0 and at j = 1..m-2!
    vals = []
    for j in range(1, m - 1):
        psi = j * mp.pi / (2 * m)
        s = mp.sin(psi) / mp.sin(psi + th)
        t = s**(2 * m)
        u_val = t / ((1 + t)**2)
        val = sum(q2_coeffs[b] * (u_val ** b) for b in range(deg_q2 + 1))
        vals.append((j, float(u_val), '+' if val > 0 else '-'))
    # Count sign changes in vals (plus u=0 which is '+')
    seq = ['+'] + [x[2] for x in vals]
    sc = sum(1 for i in range(len(seq) - 1) if seq[i] != seq[i + 1])
    print(f"  m={m:2d} (deg Q2={deg_q2:2d}): sign changes={sc}, seq={''.join(seq)}")

print("\n=== 2. CHECK b = 1 EXACT FORMULA & SIGN CHANGES ON p = -r^{2m} IN (-1, 0) ===")
def L2m_mp(m, w):
    # ((1 + sqrt(1+4w))/2)^{2m} + ((1 - sqrt(1+4w))/2)^{2m}
    sq = mp.sqrt(1 + 4 * w)
    return ((1 + sq) / 2)**(2 * m) + ((1 - sq) / 2)**(2 * m)

for m in [4, 5, 6, 7, 8, 10]:
    n = 2 * m
    coeffs = compute_p_n_fast(n)
    max_c1 = (n - 2) // 3
    h1_coeffs = [coeffs[(n - 2 - 3*c, 1, c)] for c in range(max_c1 + 1)]
    # Check roots of h1(V) in V = p/(1-p)^3! Wait: what r values in (0.25, 1) give the roots of h1?
    # Remember the roots of F_{2m}(v) = 0 lie near theta_j = (j + 1/2) pi / (2m) in (0, pi/3), i.e. r_j = 1 / (4 cos^2 theta_j)!
    # Check signs of h1(V) at r_j = 1 / (4 cos^2(j*pi/(2m))) for j = 1 .. floor(2m/3)!
    seq = []
    for j in range(0, (2*m)//3 + 1):
        th = (j + mp.mpf('0.25')) * mp.pi / (2 * m)
        if th >= mp.pi / 3: continue
        r = 1 / (4 * mp.cos(th)**2)
        p_val = -(r**(2 * m))
        V_val = p_val / ((1 - p_val)**3)
        val = sum(h1_coeffs[c] * (V_val ** c) for c in range(max_c1 + 1))
        seq.append('+' if val > 0 else '-')
    sc = sum(1 for i in range(len(seq) - 1) if seq[i] != seq[i + 1])
    # Also check exact roots of h1(V) mapped back to r = (-p)^{1/(2m)} where p/(1-p)^3 = V!
    rts_V = np.roots([float(x) for x in h1_coeffs[::-1]])
    rts_V = sorted([r.real for r in rts_V if abs(r.imag) < 1e-6], reverse=True)
    # Map each V_k < 0 to p_k in (-1, 0) and theta_k = arccos(1/(2 sqrt(r_k))) * 2m / pi!
    th_scaled = []
    for Vv in rts_V:
        # solve p / (1-p)^3 = Vv for p in (-1, 0)
        p_sol = mp.findroot(lambda pv: pv / (1 - pv)**3 - Vv, mp.mpf(-0.01) if abs(Vv) < 0.01 else mp.mpf(-0.3))
        r_sol = (-p_sol)**(mp.mpf(1) / (2 * m))
        th_sol = mp.acos(1 / (2 * mp.sqrt(r_sol))) * (2 * m) / mp.pi
        th_scaled.append(round(float(th_sol), 3))
    print(f"  m={m:2d} (deg h1={max_c1:2d}): roots in 2m*theta/pi = {th_scaled}")
