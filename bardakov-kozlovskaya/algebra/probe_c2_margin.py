"""
Check c = 2, 3, 4 sign-change points AND check whether c = 2 admits a clean algebraic proof
(either via Maclaurin/gamma-expansion of Q_2(1+u, u) or via dominance at s_j = sin(j pi/2m)/sin((j+2) pi/2m))!
Also check b = 1 roots relative to 4^{-2m}!
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 60

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

print("=== 1. c = 2 AT s_j = sin(j pi/2m) / sin((j+2) pi/2m) FOR j = 2..m-1 (INCLUDING s_{m-1}=1, u=1/4!) ===")
for m in [4, 5, 6, 8, 10, 12, 16]:
    n = 2 * m
    coeffs = compute_p_n_fast(n)
    deg_q2 = m - 3
    q2_coeffs = [coeffs[(n - 6 - 2*b, b, 2)] for b in range(deg_q2 + 1)]
    th = mp.pi / m
    seq = []
    min_ratio = mp.mpf('1e99')
    worst_j = -1
    for j in range(2, m):
        psi = j * mp.pi / (2 * m)
        s = mp.sin(psi) / mp.sin(psi + th)
        t = s**(2 * m)
        u_val = t / ((1 + t)**2)
        val = sum(q2_coeffs[b] * (u_val ** b) for b in range(deg_q2 + 1))
        seq.append('+' if val > 0 else '-')
        # Compare t^2 * (1+t)^{2m-6} * Q2(1, u) with dominant term (1-t)^{2m} * 2 * (1-s)^{-2m} * Re((1 - s e^{i pi/m})^{-2m})!
        # Wait: t^2 Q_2(1+t, t) = t^2 (1+t)^{2m-6} Q_2(1, u)!
        lhs_val = (t**2) * ((1 + t)**(2*m - 6)) * val
        if j < m - 1:
            dom = ((1 - t)**(2*m)) * 2 * ((1 - s)**(-2*m)) * mp.re((1 - s * mp.exp(1j * th))**(-2*m))
            # Expected sign of dom is (-1)^j!
            r = lhs_val / dom
            if abs(r) < min_ratio or r < 0:
                if r < min_ratio: min_ratio = r; worst_j = j
        else:
            # At j = m-1 (s = 1, t = 1): limit as s -> 1 of ((1-t)/(1-s))^{2m} = (2m)^{2m}!
            dom = ((2*m)**(2*m)) * 2 * mp.re((1 - mp.exp(1j * th))**(-2*m))
            r = lhs_val / dom
            if r < min_ratio: min_ratio = r; worst_j = j
    sc = sum(1 for i in range(len(seq) - 1) if seq[i] != seq[i + 1])
    print(f"  m={m:2d} (deg={deg_q2:2d}): sc={sc} == {deg_q2} ({sc==deg_q2}), seq={''.join(seq)}, min(Q2/dom)={float(min_ratio):.5f} (j={worst_j})")

print("\n=== 2. ALSO TEST c = 3 (u < 0, s_j = sin(j pi/2m)/sin((j+3) pi/2m)) AND c = 4 (u in (0,1/4)) ===")
for c_test in [3, 4]:
    print(f"  --- c = {c_test} ---")
    for m in [6, 8, 10, 12, 16]:
        n = 2 * m
        deg_qc = (n - 3 * c_test) // 2
        if deg_qc < 1: continue
        coeffs = compute_p_n_fast(n)
        qc = [coeffs[(n - 3*c_test - 2*b, b, c_test)] for b in range(deg_qc + 1)]
        th = c_test * mp.pi / (2 * m)
        # Note: as s: 0 -> 1, -2m*arg(1 - s e^{i th}) goes from 0 to (m - c/2) pi!
        # So j goes up to floor(m - c/2)!
        j_max = m - (c_test + 1) // 2
        seq = []
        for j in range(1, j_max + 1):
            psi = j * mp.pi / (2 * m)
            if psi + th > mp.pi: break
            s = mp.sin(psi) / mp.sin(psi + th)
            if s > 1: s = mp.mpf(1)
            t = s**(2 * m)
            u_val = (-t / ((1 - t)**2)) if (c_test % 2 == 1) else (t / ((1 + t)**2))
            if c_test % 2 == 1 and t == 1: continue
            val = sum(qc[b] * (u_val ** b) for b in range(deg_qc + 1))
            seq.append('+' if val > 0 else '-')
        sc = sum(1 for i in range(len(seq) - 1) if seq[i] != seq[i + 1])
        print(f"    m={m:2d} (deg={deg_qc:2d}): sc={sc}, seq={''.join(seq)}")
