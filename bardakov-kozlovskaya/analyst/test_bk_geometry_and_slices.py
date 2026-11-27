"""
EXPLOSIVE DIAGNOSTIC:
1. Check real-rootedness of P_{2m}(1, e2, e3) in e3 (fixed e2), in e2 (fixed e3),
   and in e1 (fixed e2, e3), AND check 2-var real stability of P_{2m}(e1, e2, e3) itself!
2. Check odd-n G_{n,c}(t) real-rootedness and sign law (-1)^{ab+a+b+n}.
3. Check the exact geometry of the real curve (E2(y), E3(y)) for w = -1/2 + i*y (0 < y < sqrt(3)/2)
   AND w in R (the OTHER real symmetry axis of p_{2m}(1, w^{2m}, (1+w)^{2m}) = 0!!).
   WAIT A SECOND!! THERE ARE TWO REAL AXES IN THE w-PLANE:
   - Axis 1: w in R  --> x=1, y=w^{2m} >= 0, z=(1+w)^{2m} >= 0 (ALL THREE x,y,z >= 0 -> e1,e2,e3 >= 0!!)
   - Axis 2: Re(w) = -1/2 (w = -1/2 + i*y) --> x=1, y=conj(z) -> e1, e2 in R, e3 > 0!!
   - Axis 3: |w| = 1 (w = e^{i*phi}) --> wait: 1/w symmetry gives |1+w| = 1 too!
"""
from math import comb
import numpy as np

def p_coeffs_table(n):
    powers = [[0]] + [[n * ((-1)**(n*r)) * comb(n*r, n*k) for k in range(r+1)] for r in range(1, n+1)]
    a = [[1]]
    for j in range(1, n+1):
        row = [0]*(j+1)
        for r in range(1, j+1):
            for u, au in enumerate(a[j-r]):
                if au:
                    for v, bv in enumerate(powers[r]): row[u+v] -= au*bv
        a.append([v//j for v in row])
    return a

def extract_e_coeffs(n):
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
    coeffs = {}
    G = {}
    for c in range(n//3 + 1):
        d = n - 3*c
        sl = [poly[d-j][j] for j in range(d+1)]
        rem = sl[:]
        gammas = []
        for b in range(d//2 + 1):
            gb = rem[b]; gammas.append(gb)
            coeffs[(d-2*b, b, c)] = gb
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
    return coeffs, G

print("=== 1. Check P_{2m}(1, e2, e3) real-rootedness in e3 for fixed real e2 ===")
for n in [6, 8, 10, 12, 14, 16]:
    C, G = extract_e_coeffs(n)
    cplx_e2 = []
    pos_e3_neg_e2 = []
    for e2v in [-10.0, -2.0, -0.5, -0.1, -0.01, 0.0, 0.01, 0.05, 0.1, 0.2, 0.24, 0.26, 0.3, 0.5, 1.0, 5.0]:
        # poly in e3: coeff of e3^c is sum_b c_{n-2b-3c, b, c} e2v^b
        p = [sum(float(C[(n-2*b-3*c, b, c)]) * (e2v**b) for b in range((n-3*c)//2 + 1)) for c in range(n//3 + 1)]
        mx = max(abs(v) for v in p)
        r = np.roots([v/mx for v in p][::-1])
        if np.max(np.abs(r.imag)) > 1e-6:
            cplx_e2.append((e2v, int(np.sum(np.abs(r.imag)>1e-6))))
        else:
            if e2v <= 0 and np.any(r.real <= 0):
                pos_e3_neg_e2.append((e2v, float(r.real.min())))
    print(f"n={n:2d}: e3-roots complex at e2={cplx_e2}; nonpos at e2<=0: {pos_e3_neg_e2}")

print("\n=== 2. Check P_{2m}(1, e2, e3) real-rootedness in e2 for fixed real e3 ===")
for n in [6, 8, 10, 12, 14]:
    C, _ = extract_e_coeffs(n)
    cplx_e3 = []
    for e3v in [-1.0, -0.01, -1e-4, -1e-6, 1e-6, 1e-4, 0.01, 1.0]:
        p = [sum(float(C[(n-2*b-3*c, b, c)]) * (e3v**c) for c in range((n-2*b)//3 + 1)) for b in range(n//2 + 1)]
        mx = max(abs(v) for v in p)
        r = np.roots([v/mx for v in p][::-1])
        if np.max(np.abs(r.imag)) > 1e-6:
            cplx_e3.append((e3v, int(np.sum(np.abs(r.imag)>1e-6))))
    print(f"n={n:2d}: e2-roots complex at e3={cplx_e3}")

print("\n=== 3. Check P_{2m}(e1, e2, e3) real-rootedness in e1 for fixed (e2, e3) ===")
for n in [6, 8, 10, 12]:
    C, _ = extract_e_coeffs(n)
    # Test (e2, e3) in all 4 quadrants!
    for (s2, s3) in [(1,1), (1,-1), (-1,1), (-1,-1)]:
        # Wait: at e2=-1, e3=+1 or small e3!
        cplx_cnt = 0
        for mag3 in [1e-4, 1e-2, 0.1, 1.0]:
            e2v, e3v = float(s2), float(s3) * mag3
            p = [0.0]*(n+1)
            for (a,b,c), v in C.items():
                p[a] += float(v) * (e2v**b) * (e3v**c)
            mx = max(abs(v) for v in p)
            r = np.roots([v/mx for v in p][::-1])
            if np.max(np.abs(r.imag)) > 1e-5:
                cplx_cnt += 1
        print(f"n={n:2d}, sign(e2,e3)=({s2:+d},{s3:+d}): e1-roots complex in {cplx_cnt}/4 tests")

print("\n=== 4. Check P_{2m}(e1, e2, e3) 2-variable stability in (e2, e3) for e1=1 ===")
# Wait: along e2 = a*s + b, e3 = c*s + d (linear in e2 and e3 with e1=1)!
rng = np.random.default_rng(1)
for n in [6, 8, 10, 12]:
    C, _ = extract_e_coeffs(n)
    # Test rays with a > 0, c > 0 vs a < 0, c > 0!
    for sgn_a in [+1, -1]:
        fails = 0
        for _ in range(40):
            alpha = sgn_a * rng.uniform(0.1, 3.0)
            gamma = rng.uniform(0.001, 0.5)
            beta = rng.uniform(-1.0, 1.0)
            delta = rng.uniform(-0.1, 0.1)
            # poly in s of deg <= n//2
            deg = n // 2
            ps = np.zeros(deg + 1)
            for (a,b,c), v in C.items():
                p1 = np.array([1.0])
                for _ in range(b): p1 = np.convolve(p1, [alpha, beta])
                for _ in range(c): p1 = np.convolve(p1, [gamma, delta])
                ps[len(ps)-len(p1):] += float(v) * p1
            nz = np.nonzero(np.abs(ps) > 1e-25*np.max(np.abs(ps)))[0]
            r = np.roots(ps[nz[0]:])
            if np.max(np.abs(r.imag)) > 1e-5:
                fails += 1
        print(f"n={n:2d}, dir(e2,e3)=({sgn_a:+d},+1): complex rays = {fails}/40")

print("\n=== 5. Check ODD n = 3..21 G_{n,c}(t) real-rootedness and unified sign law ===")
for n in range(3, 22, 2):
    C, G = extract_e_coeffs(n)
    bad = 0
    anom = []
    for c in range(1, n//3 + 1):
        d = n - 3*c
        g = G[c]
        for b, val in enumerate(g):
            a = d - 2*b
            if val == 0 or (1 if val > 0 else -1) != (-1)**(a*b + a + b + n):
                bad += 1
        if len(g) <= 1: continue
        mx = max(abs(v) for v in g)
        r = np.roots([float(v/mx) for v in g][::-1])
        if np.max(np.abs(r.imag)) > 1e-6:
            anom.append(f"c={c}(d={d}):CPLX(im={np.max(np.abs(r.imag)):.2e})")
        else:
            rr = np.sort(r.real)
            if d % 2 == 1 and not np.all(rr < 0): anom.append(f"c={c}(d={d}):NOT_NEG")
            if d % 2 == 0 and not (np.all(rr > 0) and np.all(rr <= 0.25+1e-6)): anom.append(f"c={c}(d={d}):OUT_0_1/4")
    print(f"n={n:2d}: bad_signs={bad}, anomalies={anom if anom else 'NONE!'}")
