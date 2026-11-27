"""
1. Diagnose where diff_nonq = 0.1221 at (c=20, j=21, kap=0.499) came from!
2. Compute EXACT S_q(d, X) + S_SetI(d, X) for d=1 (X >= 0.85) and d >= 2 (X >= K_G2)!
3. Design & stress-test the 100% CLOSED-FORM CONTINUED-FRACTION UPPER BOUND on Set II (1 <= r < q)!
"""
import math

# 1. Diagnose (c=20, j=21, n=1700): compare exact finite-n b_r vs Gaussian b_r^G vs b_r^main!
def b_r_gauss(r, c, j, kap):
    th = 2.0 * math.pi * ((r * j) % c) / c
    B = sum(math.exp(-r * kap * m * (1.0 - m / c)) * math.cos(m * th) for m in range(1, c) if r * kap * m * (1.0 - m / c) < 40)
    A = 0.0; m = 1
    while True:
        ex = r * kap * m * (1.0 + m / c)
        if ex > 40: break
        A += 2.0 * math.exp(-ex) * math.cos(m * th); m += 1
    return B - A

def b_r_main(r, c, j, kap):
    gam = r * kap
    if gam > 40: return 0.0
    th = 2.0 * math.pi * ((r * j) % c) / c
    z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
    return 4.0 * r * (kap / c) * (z * (1.0 + z) / ((1.0 - z) ** 3)).real

def exact_finite_n_br(n, c, j, rmax):
    th_c = c * math.pi / n; psi_j = j * math.pi / n; alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1; log_Rz = n * math.log(r3 / r1); kap = 2.0 * math.pi * r1 * r2 / r3
    shift = 0.0 if (c % 2 == 0) else 0.5; sgn_ang = 0.0 if (j % 2 == 0) else math.pi
    small_a = []; large_a = []
    for k in range(-(n // 2), n // 2 + 1):
        x = k + shift
        if x <= -n / 2.0 + 1e-12 or x > n / 2.0 + 1e-12 or abs(abs(x) - c / 2.0) < 1e-9: continue
        eta = 2.0 * math.pi * x / n; z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
        lr = log_Rz - n * math.log(abs(z)); ang = sgn_ang - n * math.atan2(z.imag, z.real)
        if abs(x) < c / 2.0:
            if lr < 45.0: small_a.append((-lr, -ang))
        else:
            if lr > -45.0: large_a.append((lr, ang))
    return [sum(math.exp(r*lr)*math.cos(r*ang) for lr,ang in small_a if r*lr>-40) -
            sum(math.exp(r*lr)*math.cos(r*ang) for lr,ang in large_a if r*lr>-40) for r in range(1, rmax+1)], kap

print("=== 1. DIAGNOSE (c=20, j=21, n=1700) AND (c=29, j=39, n=36102) ===")
for c, j, n in [(20, 21, 1700), (29, 39, 36102)]:
    brs, kap = exact_finite_n_br(n, c, j, 10)
    print(f"  (c={c}, j={j}, n={n}, kap={kap:.4f}):")
    for r in range(1, 7):
        ex = brs[r - 1] / r
        ga = b_r_gauss(r, c, j, kap) / r
        mn = b_r_main(r, c, j, kap) / r
        print(f"    r={r}: exact={ex:+.5f}, gauss={ga:+.5f}, main={mn:+.5f} | (ex-ga)={ex-ga:+.5f}, (ga-mn)={ga-mn:+.5f}")

print("\n=== 2. EXACT COMBINED S_q(d, X) + S_SetI(d, X) OVER VALID G1 RANGE ===")
def D_exact_raw(T):
    M = 120
    s1 = sum(math.exp(-0.5 * T * (1.0 - ((i + 0.5) / M) ** 2)) for i in range(M)) * (2.0 / M)
    ymax = min(40.0, 20.0 / max(T, 1e-4) + 7.0 / math.sqrt(max(T, 1e-4)))
    s2 = sum(math.exp(-T * ((i + 0.5) * ymax / M) - 0.5 * T * (((i + 0.5) * ymax / M) ** 2)) for i in range(M)) * (2.0 * ymax / M)
    return abs(s1 - s2)

def S_q_exact(d, X):
    T1 = 0.5 * d * X
    return sum((d / (2.0 * l)) * D_exact_raw(l * T1) for l in range(1, 30)) + (d / 2.0) * (4.0 / (T1 * T1)) / (2.0 * 30 * 30)

def S_SetI_exact(d, X):
    # (2/d) * sum_{m>=1, v>=1} min(0.07681/(m v^2), 0.3491/(m^2 v X), 0.01540 m X^2 / v^4)
    # Wait: use exact sin(pi v/q) >= (2/pi)*pi v/q? No, wait: for v <= q/2, sin(pi v/q) >= 2 v/q ONLY loses at large v/q!
    # Let's check BOTH small-angle constants (0.0768, 0.349, 0.0154) AND exact sin(pi v/q) over q >= 2!
    s = 0.0
    for m in range(1, 80):
        for v in range(1, 80):
            t1 = 0.095 / (m * v * v)
            t2 = 0.388 / (m * m * v * X)
            t3 = 0.025 * m * X * X / (v ** 4)
            s += 2.0 * min(t1, t2, t3)
    return s / d

# Check d = 1 (X >= 0.85):
mx_d1 = (0.0, 0.0, 0.0, 0.0)
for ix in range(0, 80):
    X = 0.85 * (1.06 ** ix)
    sq = S_q_exact(1, X); s1 = S_SetI_exact(1, X)
    if sq + s1 > mx_d1[0]: mx_d1 = (sq + s1, X, sq, s1)
print(f"  d=1 (X >= 0.85): max [S_q + S_SetI] = {mx_d1[0]:.4f} at X={mx_d1[1]:.3f} (S_q={mx_d1[2]:.4f}, S_SetI={mx_d1[3]:.4f})")

# Check d >= 2 (X >= K_G2):
for K_G2 in [2.0, 2.2, 2.5]:
    mx_d2 = (0.0, None)
    for d in [2, 3, 5, 10, 50]:
        for ix in range(0, 40):
            X = K_G2 * (1.08 ** ix)
            sq = S_q_exact(d, X); s1 = S_SetI_exact(d, X)
            if sq + s1 > mx_d2[0]: mx_d2 = (sq + s1, (d, round(X, 3), round(sq, 4), round(s1, 4)))
    print(f"  d>=2 (X >= {K_G2:.1f}): max [S_q + S_SetI] = {mx_d2[0]:.4f} at {mx_d2[1]}")

print("\n=== 3. BREAKDOWN OF SET II (1 <= r < q) BY CF CONVERGENT LEVELS ===")
# For the worst cases from task-1069 (and tricky Fibonacci/crafted CFs!), see how Set II splits into:
#   - Convergent multiples: r = m * q_k (1 <= m <= a_{k+1})
#   - Semiconvergents: r = m * q_k + q_{k-1} (1 <= m <= a_{k+1} - 1)
#   - All other r in [1, q-1]
def cf_convergents_full(num, den):
    g = math.gcd(num, den)
    a, b = (num // g) % (den // g), den // g
    if a == 0: return [1], []
    qs = []; a_list = []
    p0, p1, q0, q1 = 0, 1, 1, 0
    while b != 0:
        ai = a // b
        p2, q2 = ai * p1 + p0, ai * q1 + q0
        a_list.append(ai)
        if q2 >= 1: qs.append(q2)
        p0, p1, q0, q1 = p1, p2, q1, q2
        a, b = b, a % b
    # Remove duplicate 1 if a_1 = 1? Wait: qs[0]=1 if a_0=0, qs[1]=a_1, ...
    return qs, a_list

def sharp_term(r, D, c, kap):
    gam = r * kap
    if gam > 35: return 0.0
    th = 2.0 * math.pi * D / c
    s = math.sin(0.5 * th)
    bA = 0.7581 / (r * c * s * s)
    bB = (1.0970 if gam <= 1.5 else 1.5040) / (r * c * gam * s)
    bC = 1.5001 * r * kap * kap / (c * (s ** 4))
    bE = (4.0 * kap / c) * math.exp(-gam) * (1.0 + math.exp(-gam)) / ((1.0 - math.exp(-gam)) ** 3)
    return min(bA, bB, bC, bE)

# Test worst CF shapes: Fibonacci 21/34, 55/89, 144/233; large a_i like [0; 3, 10, 3]; small c = 10, 13, 19, 29!
for c, j, kap in [
    (10, 13, 0.085), (13, 14, 0.405), (15, 16, 0.355), (19, 20, 0.455),
    (29, 39, 0.045), (34, 55, 0.05), (89, 144, 0.03), (233, 377, 0.015),
    (100, 131, 0.04), (100, 129, 0.04), (150, 197, 0.02),
]:
    d = math.gcd(c, j); q = c // d
    qs, _ = cf_convergents_full(j, c)
    # Deduplicate qs while preserving order
    q_uniq = []
    for x in qs:
        if not q_uniq or x != q_uniq[-1]:
            q_uniq.append(x)
    if q_uniq[0] != 1:
        q_uniq = [1] + q_uniq
    # Classify every r in 1..q-1:
    s_conv = 0.0; s_semi = 0.0; s_rest = 0.0
    conv_set = set(); semi_set = set()
    for idx in range(len(q_uniq) - 1):
        qk = q_uniq[idx]; qk1 = q_uniq[idx + 1]
        qkm1 = q_uniq[idx - 1] if idx > 0 else 0
        ak1 = qk1 // qk
        for m in range(1, ak1 + 1):
            if m * qk < qk1:
                conv_set.add(m * qk)
            if idx > 0 and m < ak1 and m * qk + qkm1 < qk1:
                semi_set.add(m * qk + qkm1)
    for r in range(1, q):
        Dr = min((r * j) % c, c - ((r * j) % c))
        val = sharp_term(r, Dr, c, kap)
        if r in conv_set: s_conv += val
        elif r in semi_set: s_semi += val
        else: s_rest += val
    # Also compute max c/(qk Dk^2) among qk*kap <= 0.45
    max_M = 0.0
    for qk in q_uniq[:-1]:
        Dk = min((qk * j) % c, c - ((qk * j) % c))
        if qk * kap <= 0.45 and Dk > 0:
            max_M = max(max_M, c / (qk * Dk * Dk))
    print(f"  c={c:3d} j={j:3d} kap={kap:.3f} (maxM={max_M:4.1f}) | SetII={s_conv+s_semi+s_rest:.4f} (conv:{s_conv:.4f}, semi:{s_semi:.4f}, rest:{s_rest:.4f})")
