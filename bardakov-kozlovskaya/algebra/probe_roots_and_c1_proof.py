"""
1. Check exact root intervals of Q_c(1, u) (for even c vs odd c), R_b(1, v), and S_a(1, W) for n = 6..24!
   Specifically: do ALL roots of even-c Q_c(1, u) lie in (0, 1/4)?
   And do consecutive slices Q_c, Q_{c+2} or R_b, R_{b+2} or S_a, S_{a+2} have a clean structure?
2. Test our 10-line c = 1 real-rootedness proof via k=0 dominance on s in (0, 1) where -2m*phi(s_j) = j*pi:
   Check |term(k=0)| / (sum_{k=1}^{m-1} |term(k)| + |(1-t)/(1+t)^2|) at every s_j for m = 2..20!
3. Check the EXACT formula for c = 2 slice Q_2(e1, e2) in terms of Veronese polynomials H_1, H_2!
4. Check if P_{2m}(e1, e2, e3) has a simple relation between adjacent m or a low-order PDE!
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

def get_roots_scaled(row):
    d = len(row) - 1
    if d <= 0: return np.array([])
    s = (abs(row[0]) / abs(row[-1])) ** (1.0 / d)
    sc = [float(row[k] * (s ** k)) for k in range(d + 1)]
    mx = max(abs(x) for x in sc)
    sc = [x / mx for x in sc]
    rts = np.roots(sc[::-1])
    if np.max(np.abs(rts.imag)) > 1e-5:
        return None
    return np.sort(rts.real * s)

print("=== 1. EXACT ROOT RANGES OF Q_c(1, u) FOR EVEN c >= 2 AND ODD c >= 1 (n = 8..24) ===")
for n in [8, 10, 12, 14, 16, 18, 20, 24]:
    coeffs = compute_p_n_fast(n)
    even_max = 0.0; even_min = 1e9
    odd_max = -1e9; odd_min = 1e9
    for c in range(1, n // 3 + 1):
        deg_u = (n - 3 * c) // 2
        if deg_u < 1: continue
        row = [coeffs[(n - 3*c - 2*b, b, c)] for b in range(deg_u + 1)]
        r = get_roots_scaled(row)
        if c % 2 == 0:
            even_min = min(even_min, float(r[0]))
            even_max = max(even_max, float(r[-1]))
        else:
            odd_min = min(odd_min, float(r[0]))
            odd_max = max(odd_max, float(r[-1]))
    print(f"  n={n:2d}: even c>=2 roots in [{even_min:.6f}, {even_max:.6f}] (all < 0.25? {even_max < 0.25}) | odd c>=1 roots in [{odd_min:.5f}, {odd_max:.6f}]")

print("\n=== 2. TEST c=1 PHASE-DOMINANCE PROOF (k=0 vs rest + subtracted term) FOR m = 2..18 ===")
for m in range(2, 19):
    th0 = mp.pi / (2 * m)
    # Solve -2m * arg(1 - s e^{i th0}) = j * pi for j = 1..m-1!
    # arg(1 - s e^{i th0}) = -psi_j where psi_j = j*pi/(2m).
    # tan(psi_j) = s sin(th0) / (1 - s cos(th0))  =>  s = sin(psi_j) / sin(psi_j + th0)!!
    # OH!! WHAT A CLEAN EXACT FORMULA FOR s_j!!
    # Check: 1 - s e^{i th0} = 1 - sin(psi)/sin(psi+th0) * (cos th0 + i sin th0)
    #   = (sin(psi+th0) - sin(psi)cos(th0) - i sin(psi)sin(th0)) / sin(psi+th0)
    #   = sin(th0) * (cos(psi) - i sin(psi)) / sin(psi+th0) = (sin(th0)/sin(psi+th0)) * e^{-i psi}!!
    # OH!! MY!! GOD!! |1 - s_j e^{i th0}| = sin(th0) / sin((j+1)*pi/(2m))!!
    min_margin = mp.mpf('1e99')
    worst_j = -1
    signs_ok = True
    for j in range(1, m):
        psi = j * mp.pi / (2 * m)
        s = mp.sin(psi) / mp.sin(psi + th0)
        t = s**(2 * m)
        # Main term (k=0, paired with k=2m-1 -> factor 1/m, real part = (-1)^j * (sin(psi+th0)/sin(th0))^{2m}):
        main = (mp.sin(psi + th0) / mp.sin(th0))**(2 * m) / m
        # Remaining k = 1..m-1 terms:
        rem = mp.mpf(0)
        for k in range(1, m):
            thk = (2*k + 1) * mp.pi / (2 * m)
            rem += abs(1 - s * mp.exp(1j * thk))**(-2 * m) / m
        sub = (1 - t) / ((1 + t)**2)
        ratio = main / (rem + sub)
        if ratio < min_margin:
            min_margin = ratio
            worst_j = j
        # Check exact sign of actual value!
        val = sum(mp.re((1 - s * mp.exp(1j * (2*k+1)*mp.pi/(2*m)))**(-2*m)) for k in range(m)) / m - sub
        if (val > 0) != (j % 2 == 0):
            signs_ok = False
    print(f"  m={m:2d} (n={2*m:2d}): signs_ok={signs_ok}, min(main / (rem + sub)) = {float(min_margin):.4f} (at j={worst_j})")

print("\n=== 3. CHECK EXACT ROOT RANGES OF R_b(1, v) (b-slices) AND S_a(1, W) (a-slices) AT n = 12, 18 ===")
for n in [12, 18]:
    coeffs = compute_p_n_fast(n)
    print(f"  n={n} b-slices R_b(1, v):")
    for b in range(0, min(6, n // 2)):
        max_c = (n - 2 * b) // 3
        if max_c < 1: continue
        row = [coeffs[(n - 2*b - 3*c, b, c)] for c in range(max_c + 1)]
        r = get_roots_scaled(row)
        print(f"    b={b}: deg={max_c}, roots={[round(float(x), 8) for x in r]}")
    print(f"  n={n} a-slices S_a(1, W):")
    for a in range(0, min(6, n - 5)):
        c0 = a % 2
        row = [coeffs[(a, (n - a - 3*c)//2, c)] for c in range(c0, (n - a)//3 + 1, 2)]
        r = get_roots_scaled(row)
        print(f"    a={a}: deg={len(row)-1}, roots={[round(float(x), 8) for x in r]}")
