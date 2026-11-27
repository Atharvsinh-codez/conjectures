"""
1. Push G_{n,c}(t) root check up to n = 44 (using exact integer peeling + mpmath when needed).
2. Check interlacing of G_{n,c}(t) vs G_{n,c+2}(t) and G_{n-2,c}(t).
3. Check the growth/scaling of log|c_{a,b,c}^{(n)}| / n across the simplex a+2b+3c = n!
   (Wait: does log|c_{a,b,c}| grow like n, or like n^2?! CHECK c_{0,0,n/3}!)
"""
from math import comb, log
import numpy as np
import mpmath as mp

def p_coeffs_table(n):
    powers = [[0]] + [[n * ((-1)**(n*r)) * comb(n*r, n*k) for k in range(r+1)] for r in range(1, n+1)]
    a = [[1]]
    for j in range(1, n+1):
        row = [0]*(j+1)
        for r in range(1, j+1):
            prev = a[j-r]; pr = powers[r]
            for u, au in enumerate(prev):
                if au:
                    for v, bv in enumerate(pr): row[u+v] -= au*bv
        a.append([v//j for v in row])
    return a

def extract_G_fast(n):
    a_tab = p_coeffs_table(n)
    def mul_e1(P, deg):
        out = [[0]*(deg+2) for _ in range(deg+2)]
        for i in range(deg+1):
            for j in range(deg+1-i):
                v = P[i][j]
                if v: out[i+1][j]+=v; out[i][j+1]+=v; out[i][j]+=v
        return out
    def mul_e2(P, deg):
        out = [[0]*(deg+3) for _ in range(deg+3)]
        for i in range(deg+1):
            for j in range(deg+1-i):
                v = P[i][j]
                if v: out[i+1][j+1]+=v; out[i+1][j]+=v; out[i][j+1]+=v
        return out
    E1 = [[[1]]]
    for r in range(1, n+1): E1.append(mul_e1(E1[-1], r-1))
    poly = [[0]*(n+1) for _ in range(n+1)]
    for j_tot in range(n+1):
        for ky in range(j_tot+1): poly[j_tot-ky][ky] = a_tab[j_tot][ky]
    G = {}
    for c in range(n//3 + 1):
        d = n - 3*c
        sl = [poly[d-j][j] for j in range(d+1)]
        rem = sl[:]
        gammas = []
        for b in range(d//2 + 1):
            gb = rem[b]; gammas.append(gb)
            for r in range(d-2*b+1): rem[b+r] -= gb*comb(d-2*b, r)
        G[c] = gammas
        if c == n//3: break
        B = d//2; cur_deg = d - 2*B
        acc = [[gammas[B]*E1[cur_deg][i][j] for j in range(cur_deg+1)] for i in range(cur_deg+1)]
        for b in range(B-1, -1, -1):
            acc = mul_e2(acc, cur_deg); cur_deg += 2; gb = gammas[b]; e1p = E1[cur_deg]
            for i in range(cur_deg+1):
                for j in range(cur_deg+1-i): acc[i][j] += gb*e1p[i][j]
        for i in range(d+1):
            for j in range(d+1-i): poly[i][j] -= acc[i][j]
        poly = [[poly[i+1][j+1] for j in range(d-2)] for i in range(d-2)]
    return G

mp.mp.dps = 80

print("=== 1. Check G_{n,c}(t) real-rootedness for n = 32, 34, 36, 38, 40, 42, 44 ===")
G_store = {}
for n in [12, 18, 24, 30, 32, 34, 36, 38, 40, 42, 44]:
    G = extract_G_fast(n)
    G_store[n] = G
    if n < 32: continue
    anom = []
    for c in range(1, n//3 + 1):
        g = G[c]
        if len(g) <= 1: continue
        mx = max(abs(v) for v in g)
        r = np.roots([float(v/mx) for v in g][::-1])
        if np.max(np.abs(r.imag)) > 1e-6:
            rm = mp.polyroots([mp.mpf(v) for v in g[::-1]], maxsteps=300, extraprec=60)
            max_im = max(abs(z.imag) for z in rm)
            if max_im > mp.mpf('1e-20'):
                anom.append(f"c={c}(deg={len(g)-1}):CPLX({float(max_im):.2e})")
            else:
                rr = sorted(float(z.real) for z in rm)
                if c % 2 == 1 and not all(x < 0 for x in rr): anom.append(f"c={c}:NOT_NEG")
                if c % 2 == 0 and not all(0 < x <= 0.250001 for x in rr): anom.append(f"c={c}:OUT_0_1/4")
        else:
            rr = np.sort(r.real)
            if c % 2 == 1 and not np.all(rr < 0): anom.append(f"c={c}:NOT_NEG")
            if c % 2 == 0 and not np.all((rr > 0) & (rr <= 0.25+1e-6)): anom.append(f"c={c}:OUT_0_1/4")
    print(f"n={n:2d}: anomalies = {anom if anom else 'NONE! 100% REAL IN (-inf,0)/(0,1/4]!'}")

print("\n=== 2. Check asymptotic scaling of log|c_{a,b,c}^{(n)}| / n^2 ===")
# Look at c_{0, 0, n/3} (for 6 | n) and c_{n-3, 0, 1} and interior (n/3, n/6, n/9)!
for n in [12, 18, 24, 30, 36, 42]:
    G = G_store[n]
    c_max = n // 3
    val_corner = abs(G[c_max][0])
    val_c1 = abs(G[1][0])
    val_mid = abs(G[n//6][n//12])
    print(f"n={n:2d}: log|c(c=1,b=0)|/n={log(val_c1)/n:.5f} | log|c(0,0,n/3)|/n^2={log(val_corner)/(n*n):.6f} | log|c(mid)|/n^2={log(val_mid)/(n*n):.6f}")

print("\n=== 3. Check interlacing between G_{n,c} and G_{n-2,c} or G_{n,c+2} (n=18, 24) ===")
def get_roots_sorted(g):
    mx = max(abs(v) for v in g)
    r = np.roots([float(v/mx) for v in g][::-1])
    assert np.max(np.abs(r.imag)) < 1e-5
    return np.sort(r.real)

def check_interlace(r1, r2):
    # deg r2 = deg r1 + 1
    if len(r2) == len(r1) + 1:
        return all(r2[j] < r1[j] < r2[j+1] for j in range(len(r1)))
    return None

# Check G_{24, c} vs G_{22, c} (Wait: deg G_{n,c} = floor((n-3c)/2), so n -> n+2 raises deg by +1!)
G22 = extract_G_fast(22)
G24 = G_store[24]
print("Interlacing G_{22,c} << G_{24,c} (deg increases by +1):")
for c in range(1, 22//3 + 1):
    r1 = get_roots_sorted(G22[c])
    r2 = get_roots_sorted(G24[c])
    print(f"  c={c} (deg {len(r1)}->{len(r2)}): strict_inside={check_interlace(r1, r2)}, r1_range=[{r1[0]:.3e},{r1[-1]:.3e}], r2_range=[{r2[0]:.3e},{r2[-1]:.3e}]")
