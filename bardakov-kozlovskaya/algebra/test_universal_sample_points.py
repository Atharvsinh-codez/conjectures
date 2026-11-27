"""
MOMENT OF TRUTH: Test the UNIVERSAL explicit sample points
  s_j^{(c)} = sin(j pi / (2m)) / sin((j + c) pi / (2m)),   j = c, c+1, ..., m - ceil(c/2)
across ALL even n = 2m = 4..40 and ALL c = 1, ..., floor(2m/3)!
Notice the count of points is ALWAYS (m - ceil(c/2)) - c + 1 = m - ceil(3c/2) + 1 = deg Q_c + 1!
Check:
1. Does Q_c(1, u_j) strictly alternate in sign at ALL deg Q_c + 1 points for EVERY (2m, c)?!
2. Compare Q_c(1+t_j, t_j) against the leading term of (-1)^c e_c(alpha_0, ..., alpha_{2m-1}) (or 1-var Cauchy integral)!
"""
import math
import mpmath as mp

mp.mp.dps = 80

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

print("=== UNIVERSAL SAMPLE POINT TEST: s_j^{(c)} = sin(j pi/2m) / sin((j+c) pi/2m), j = c .. m - ceil(c/2) ===")
total_slices = 0
pass_slices = 0
fails = []
for m in range(2, 21):
    n = 2 * m
    coeffs = compute_p_n_fast(n)
    max_c = n // 3
    m_fails = []
    for c in range(1, max_c + 1):
        deg_qc = (n - 3 * c) // 2
        if deg_qc == 0:
            continue
        total_slices += 1
        qc = [coeffs[(n - 3*c - 2*b, b, c)] for b in range(deg_qc + 1)]
        j_start = c
        j_end = m - (c + 1) // 2
        assert (j_end - j_start) == deg_qc, f"Count mismatch at m={m}, c={c}: {j_end-j_start} vs {deg_qc}"
        th = c * mp.pi / (2 * m)
        seq = []
        for j in range(j_start, j_end + 1):
            psi = j * mp.pi / (2 * m)
            s = mp.sin(psi) / mp.sin(psi + th)
            t = s**(2 * m)
            if c % 2 == 0:
                u_val = t / ((1 + t)**2)
            else:
                u_val = -t / ((1 - t)**2)
            val = sum(qc[b] * (u_val ** b) for b in range(deg_qc + 1))
            seq.append('+' if val > 0 else ('-' if val < 0 else '0'))
        sc = sum(1 for i in range(len(seq) - 1) if seq[i] != seq[i + 1] and seq[i] != '0' and seq[i+1] != '0')
        if sc == deg_qc:
            pass_slices += 1
        else:
            m_fails.append((c, deg_qc, sc, ''.join(seq)))
            fails.append((n, c, deg_qc, sc, ''.join(seq)))
    print(f"  n={n:2d} (m={m:2d}): fails={len(m_fails)} {m_fails[:3]}")
print(f"\nTOTAL: {pass_slices} / {total_slices} slices strictly alternate at ALL deg+1 explicit sample points! ({len(fails)} fails)")
