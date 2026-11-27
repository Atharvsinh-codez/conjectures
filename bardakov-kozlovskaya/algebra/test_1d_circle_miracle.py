"""
TEST THE 1-DIMENSIONAL CIRCLE INTEGRAL MIRACLE!
At our exact universal sample points:
  s_j^{(c)} = sin(j pi / (2m)) / sin((j + c) pi / (2m)),   j = c .. m - ceil(c/2)
  t_j = (-1)^c * (s_j^{(c)})^{2m},   u_j = t_j / (1 + t_j)^2
we ALREADY verified 120/120 slices (n=4..40) have EXACT sign (-1)^j for Q_c^{(2m)}(1, u_j)!

Now test the 1-VARIABLE Cauchy integral:
  Q_c^{(2m)}(1, u_j) = (1/(2pi)) int_{-pi}^{pi} P_{2m}(1, u_j, R_v e^{i phi}) / (R_v^c e^{i c phi}) dphi!
Evaluate P_{2m}(1, u_j, v) DIRECTLY from its definition on the circle |v| = R_v,
where R_v is the Maillot saddle radius:
  r1 = sin((j+c)pi/2m), r2 = sin(j pi/2m), r3 = sin(c pi/2m),
  R_y = (r2/r1)^{2m} = |t_j|,  R_z = (r3/r1)^{2m},
  R_v = R_y * R_z / |1 + t_j|^3! (Check also R_v = R_y * R_z / (1 + R_y)^3 etc.)
Check |avg| / avg|F| AND min_phi Re((-1)^j * integrand) across EVERY (c, j) at n = 12, 20, 30, 40!
"""
import math
import numpy as np

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

print("=== TEST 1D CIRCLE INTEGRAL FOR Q_c^{(2m)}(1, u_j) AT MAILLOT RADIUS R_v ===")
for n in [12, 20, 30, 40]:
    m = n // 2
    coeffs = compute_p_n_fast(n)
    max_c = n // 3
    # Build Q_c(1, u) polynomial evaluator for each c
    Q_polys = {}
    for c in range(max_c + 1):
        deg_u = (n - 3 * c) // 2
        Q_polys[c] = [float(coeffs[(n - 3*c - 2*b, b, c)]) for b in range(deg_u + 1)]

    phis = np.linspace(-np.pi, np.pi, 256, endpoint=False)
    e_iphi = np.exp(1j * phis)

    worst_ratio = 1e9
    worst_case = None
    all_pos_re = 0
    tot_cases = 0
    for c in range(1, max_c + 1):
        deg_qc = (n - 3 * c) // 2
        if deg_qc == 0: continue
        j_start = c
        j_end = m - (c + 1) // 2
        for j in range(j_start, j_end + 1):
            tot_cases += 1
            r1 = math.sin((j + c) * math.pi / n)
            r2 = math.sin(j * math.pi / n)
            r3 = math.sin(c * math.pi / n)
            s = r2 / r1
            Ry = s ** n
            Rz = (r3 / r1) ** n
            t = Ry if (c % 2 == 0) else -Ry
            u_val = t / ((1.0 + t) ** 2)
            # Maillot saddle radius for v = e3 / e1^3:
            # Wait: at x=1, y=t, z~Rz, e1 = 1 + t, e3 = t*z -> |v| = Ry * Rz / |1 + t|^3!
            Rv = Ry * Rz / (abs(1.0 + t) ** 3)
            # Evaluate P_{2m}(1, u_val, Rv * e^{i phi}) = sum_{k=0}^{max_c} Q_k(1, u_val) * (Rv e^{i phi})^k!
            # Use log-scale or direct eval (for n<=40 float64 handles up to 1e308, let's be careful at n=40!)
            # Actually evaluate Q_k(1, u_val) * Rv^k in high-precision mpmath or normalized float!
            # Check overflow: at n=40, |c_{a,b,c}| ~ 1e150 -> safe in float64! Wait: let's check max!
            q_at_u = []
            for k in range(max_c + 1):
                # Horner
                val = 0.0
                for coef in reversed(Q_polys[k]):
                    val = val * u_val + coef
                q_at_u.append(val)
            # integrand = sum_{k=0}^{max_c} q_at_u[k] * Rv^{k-c} * e^{i(k-c)phi}
            integrand = np.zeros_like(e_iphi, dtype=complex)
            for k in range(max_c + 1):
                term = q_at_u[k] * (Rv ** (k - c))
                integrand += term * np.exp(1j * (k - c) * phis)
            avg = np.mean(integrand)
            abs_avg = np.mean(np.abs(integrand))
            ratio = abs(avg) / abs_avg if abs_avg > 0 else 0.0
            signed_re = ((-1) ** j) * integrand.real
            if np.min(signed_re) > 0:
                all_pos_re += 1
            if ratio < worst_ratio:
                worst_ratio = ratio
                worst_case = (c, j, float(ratio), float(np.min(signed_re) / np.mean(signed_re)), float(Rv))
    print(f"  n={n:2d} (tot={tot_cases:3d} pts): pointwise_Re>0 EVERYWHERE={all_pos_re}/{tot_cases}, worst |avg|/avg|F| = {worst_ratio:.5f} at (c,j,r,minRe/avg,Rv)={worst_case}")
