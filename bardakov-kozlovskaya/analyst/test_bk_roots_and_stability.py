"""
Deep diagnostic of Bardakov-Kozlovskaya P_{2m}(e1, e2, e3):
1. Push exact integer extraction up to n = 34 (fast peeling) and check:
   - Sign law (-1)^{b(c+1)+c}
   - Roots of f_{n,c}(u) = Q_{n,c}(1+u, u) (ALL_NEG for odd c, ALL_POS for even c?)
   - Roots of G_{n,c}(t) = sum_b c_{n-2b-3c, b, c} t^b (in (-inf,0) for odd c, in (0, 1/4] for even c?)
   - Also check INTERLACING between G_{n, c}(t) and G_{n, c+2}(t) or G_{n-2, c}(t)!
2. Check roots in e3 of P_n(1, e2, e3), P_n^{even}(1, -e2, e3), -P_n^{odd}(1, e2, e3) for fixed e2!
3. Check roots of b=0 slice P_n(1, 0, e3), a=0 slice P_n(0, e2, e3), and check 2-var real stability
   of P_n^{even}(1, -e2, e3) and -P_n^{odd}(1, e2, e3) (or in sqrt(e3) / e3^2)!
4. Check exact prime factorizations / ratios of c_{a,b,c}^{(2m)} for n = 2..12!
"""
from math import comb
import numpy as np
import mpmath as mp

def p_coeffs_table(n):
    powers = [[0]] + [
        [n * ((-1)**(n*r)) * comb(n*r, n*k) for k in range(r + 1)]
        for r in range(1, n + 1)
    ]
    a = [[1]]
    for j in range(1, n + 1):
        row = [0] * (j + 1)
        for r in range(1, j + 1):
            prev = a[j - r]
            pr = powers[r]
            for u, au in enumerate(prev):
                if au == 0: continue
                for v, bv in enumerate(pr):
                    row[u + v] -= au * bv
        a.append([v // j for v in row])
    return a

def extract_e_coeffs_fast(n):
    a_tab = p_coeffs_table(n)
    def mul_e1(P, deg):
        out = [[0]*(deg+2) for _ in range(deg+2)]
        for i in range(deg+1):
            for j in range(deg+1-i):
                v = P[i][j]
                if v:
                    out[i+1][j] += v
                    out[i][j+1] += v
                    out[i][j] += v
        return out
    def mul_e2(P, deg):
        out = [[0]*(deg+3) for _ in range(deg+3)]
        for i in range(deg+1):
            for j in range(deg+1-i):
                v = P[i][j]
                if v:
                    out[i+1][j+1] += v
                    out[i+1][j] += v
                    out[i][j+1] += v
        return out

    E1 = [[[1]]]
    for r in range(1, n + 1):
        E1.append(mul_e1(E1[-1], r - 1))

    poly = [[0]*(n+1) for _ in range(n+1)]
    for j_tot in range(n+1):
        for ky in range(j_tot+1):
            poly[j_tot - ky][ky] = a_tab[j_tot][ky]

    coeffs = {}
    G_polys = {} # c -> list [c_{d,0,c}, c_{d-2,1,c}, ..., c_{d-2B,B,c}]
    for c in range(n // 3 + 1):
        d = n - 3 * c
        sl = [poly[d - j][j] for j in range(d + 1)]
        rem = sl[:]
        gammas = []
        for b in range(d // 2 + 1):
            gb = rem[b]
            gammas.append(gb)
            coeffs[(d - 2*b, b, c)] = gb
            m_pow = d - 2*b
            for r in range(m_pow + 1):
                rem[b + r] -= gb * comb(m_pow, r)
        G_polys[c] = gammas
        if c == n // 3:
            break
        B = d // 2
        cur_deg = d - 2 * B
        acc = [[gammas[B] * E1[cur_deg][i][j] for j in range(cur_deg + 1)] for i in range(cur_deg + 1)]
        for b in range(B - 1, -1, -1):
            acc = mul_e2(acc, cur_deg)
            cur_deg += 2
            gb = gammas[b]
            e1p = E1[cur_deg]
            for i in range(cur_deg + 1):
                for j in range(cur_deg + 1 - i):
                    acc[i][j] += gb * e1p[i][j]
        for i in range(d + 1):
            for j in range(d + 1 - i):
                poly[i][j] -= acc[i][j]
        new_d = d - 3
        new_poly = [[0]*(new_d + 1) for _ in range(new_d + 1)]
        for i in range(new_d + 1):
            for j in range(new_d + 1 - i):
                new_poly[i][j] = poly[i + 1][j + 1]
        poly = new_poly
    return coeffs, G_polys

print("=== 1. Check G_{n,c}(t) roots for even n = 2..30 ===")
all_G = {}
all_C = {}
for n in range(2, 32, 2):
    coeffs, G_polys = extract_e_coeffs_fast(n)
    all_G[n] = G_polys
    all_C[n] = coeffs
    anomalies = []
    for c in range(1, n // 3 + 1):
        g = G_polys[c] # low-to-high in t
        if len(g) <= 1: continue
        # normalize to float safely by dividing by max abs
        mx = max(abs(v) for v in g)
        gf = [float(v / mx) for v in g]
        r = np.roots(gf[::-1])
        max_im = float(np.max(np.abs(r.imag)))
        if max_im > 1e-6:
            anomalies.append(f"c={c}:CPLX(im={max_im:.2e},deg={len(g)-1})")
        else:
            rr = np.sort(r.real)
            if c % 2 == 1:
                if not np.all(rr < 0):
                    anomalies.append(f"c={c}:ODD_NOT_NEG({rr.min():.3f},{rr.max():.3f})")
            else:
                if not (np.all(rr > 0) and np.all(rr <= 0.25 + 1e-6)):
                    anomalies.append(f"c={c}:EVEN_OUT_0_1/4({rr.min():.4f},{rr.max():.4f})")
    print(f"n={n:2d}: G_{{n,c}}(t) anomalies = {anomalies if anomalies else 'NONE! ALL REAL IN (-inf,0) / (0,1/4]!'}")

print("\n=== 2. Check G_{n,c}(t) root ranges for n=16, 20, 24 ===")
for n in [16, 20, 24]:
    print(f"n={n}:")
    for c in range(1, n // 3 + 1):
        g = all_G[n][c]
        if len(g) <= 1: continue
        mx = max(abs(v) for v in g)
        r = np.roots([float(v/mx) for v in g][::-1])
        if np.max(np.abs(r.imag)) < 1e-6:
            rr = np.sort(r.real)
            print(f"  c={c} (deg {len(g)-1:2d}): REAL in [{rr[0]:.6f}, {rr[-1]:.6f}]")
        else:
            print(f"  c={c} (deg {len(g)-1:2d}): COMPLEX max_im={np.max(np.abs(r.imag)):.4e}, Re in [{r.real.min():.5f}, {r.real.max():.5f}]")

print("\n=== 3. Check roots of P_n(1, 0, e3), P_n^{even}(1, 0, e3), P_n^{odd}(1, 0, e3) in e3 ===")
for n in [6, 10, 14, 18, 22, 26]:
    c_dict = all_C[n]
    # Full b=0 poly in e3
    p_all = [c_dict.get((n - 3*c, 0, c), 0) for c in range(n // 3 + 1)]
    mx = max(abs(v) for v in p_all)
    r_all = np.roots([float(v/mx) for v in p_all][::-1])
    # Even c poly in w = e3^2
    p_ev = [p_all[c] for c in range(0, len(p_all), 2)]
    p_od = [p_all[c] for c in range(1, len(p_all), 2)]
    r_ev = np.roots([float(v/max(abs(x) for x in p_ev)) for v in p_ev][::-1]) if len(p_ev) > 1 else np.array([])
    r_od = np.roots([float(v/max(abs(x) for x in p_od)) for v in p_od][::-1]) if len(p_od) > 1 else np.array([])
    print(f"n={n:2d} (b=0): full max_im={np.max(np.abs(r_all.imag)):.2e} (Re:[{r_all.real.min():.3e},{r_all.real.max():.3e}]) | even(e3^2) max_im={np.max(np.abs(r_ev.imag)) if len(r_ev) else 0:.2e} | odd(e3^2) max_im={np.max(np.abs(r_od.imag)) if len(r_od) else 0:.2e}")

print("\n=== 4. Check roots of P_n(0, e2, e3) (a=0 boundary!) ===")
for n in [6, 12, 18, 24, 30]:
    c_dict = all_C[n]
    # a = 0 => 2b + 3c = n => c must be even! c = 2*k, b = n/2 - 3*k
    # Poly in u = e3^2 / e2^3: coeff of u^k is c_{0, n/2 - 3k, 2k}
    ks = list(range(n // 6 + 1))
    vals = [c_dict[(0, n//2 - 3*k, 2*k)] for k in ks]
    # Note sign is (-1)^{b} = (-1)^{n/2 - 3k} = (-1)^{n/2} * (-1)^k!
    mx = max(abs(v) for v in vals)
    r = np.roots([float(v/mx) for v in vals][::-1]) if len(vals) > 1 else np.array([])
    print(f"n={n:2d} (a=0, deg {len(vals)-1}): max_im={np.max(np.abs(r.imag)) if len(r) else 0:.2e}, roots={np.round(r.real, 5) if len(r) and np.max(np.abs(r.imag))<1e-6 else 'CPLX'}")

print("\n=== 5. Check 2-variable real stability of P_n^{even}(1, -e2, e3^2->w) and -P_n^{odd}(1, e2, e3^2->w) ===")
# Wait: test on a few real lines e2 = t0, w = t1!
rng = np.random.default_rng(42)
for n in [6, 8, 10, 12, 14, 16]:
    c_dict = all_C[n]
    # P_even(1, -e2, sqrt(w)) = sum_{k>=0, b>=0, 2b+6k<=n} (-1)^b c_{n-2b-6k, b, 2k} e2^b w^k (all coeffs > 0!)
    # -P_odd(1, e2, sqrt(w))/e3 = sum_{k>=0, b>=0, 2b+6k+3<=n} (-c_{n-2b-6k-3, b, 2k+1}) e2^b w^k (all coeffs > 0!)
    ev_fail = 0
    od_fail = 0
    for _ in range(50):
        # test ray (e2, w) = (alpha*s + beta, gamma*s + delta) with alpha, gamma > 0
        alpha, gamma = rng.uniform(0.1, 3.0, 2)
        beta, delta = rng.uniform(-2.0, 2.0, 2)
        # evaluate as poly in s
        deg_ev = n // 2
        pe = np.zeros(deg_ev + 1)
        for (a, b, c), v in c_dict.items():
            if c % 2 == 0:
                k = c // 2
                # (-1)^b * v * (alpha*s + beta)^b * (gamma*s + delta)^k
                p1 = np.array([1.0])
                for _ in range(b): p1 = np.convolve(p1, [alpha, beta])
                for _ in range(k): p1 = np.convolve(p1, [gamma, delta])
                pe[len(pe) - len(p1):] += ((-1)**b * float(v)) * p1
        # trim leading zeros
        nz = np.nonzero(np.abs(pe) > 1e-30 * np.max(np.abs(pe)))[0]
        if len(nz):
            re = np.roots(pe[nz[0]:])
            if np.max(np.abs(re.imag)) > 1e-5:
                ev_fail += 1
        po = np.zeros(deg_ev + 1)
        for (a, b, c), v in c_dict.items():
            if c % 2 == 1:
                k = (c - 1) // 2
                p1 = np.array([1.0])
                for _ in range(b): p1 = np.convolve(p1, [alpha, beta])
                for _ in range(k): p1 = np.convolve(p1, [gamma, delta])
                po[len(po) - len(p1):] += (-float(v)) * p1
        nz = np.nonzero(np.abs(po) > 1e-30 * np.max(np.abs(po)))[0]
        if len(nz):
            ro = np.roots(po[nz[0]:])
            if np.max(np.abs(ro.imag)) > 1e-5:
                od_fail += 1
    print(f"n={n:2d}: 2-var real-stability test (50 rays) -> P_even(1,-e2,w) fails={ev_fail}/50, -P_odd(1,e2,w) fails={od_fail}/50")
