"""
Test:
1. Exact finite-n b_r vs b_r^main across G1 (including large c*kap^2 up to 41!),
   checking how small |b_r^{exact} - b_r^{main}| / r is!
2. Design and test a 100% CLOSED-FORM CONTINUED-FRACTION UPPER BOUND (ZERO loop over r!)
   for sum_{q \nmid r} |b_r^main|/r + sum_{q|r} |b_r|/r across ALL of G1!
"""
import math
import cmath

def exact_finite_n_br(n, c, j, rmax):
    th_c = c * math.pi / n; psi_j = j * math.pi / n; alp1 = (j + c) * math.pi / n
    r1, r2, r3 = math.sin(alp1), math.sin(psi_j), math.sin(th_c)
    s = r2 / r1; log_Rz = n * math.log(r3 / r1)
    kap = 2.0 * math.pi * r1 * r2 / r3
    shift = 0.0 if (c % 2 == 0) else 0.5
    sgn_ang = 0.0 if (j % 2 == 0) else math.pi
    small_a = []; large_a = []
    for k in range(-(n // 2), n // 2 + 1):
        x = k + shift
        if x <= -n / 2.0 + 1e-12 or x > n / 2.0 + 1e-12: continue
        if abs(abs(x) - c / 2.0) < 1e-9: continue
        eta = 2.0 * math.pi * x / n
        z = complex(1.0 - s * math.cos(eta), -s * math.sin(eta))
        lr = log_Rz - n * math.log(abs(z))
        ang = sgn_ang - n * math.atan2(z.imag, z.real)
        if abs(x) < c / 2.0:
            if lr < 45.0:
                small_a.append((-lr, -ang))  # a_k^{-1}
        else:
            if lr > -45.0:
                large_a.append((lr, ang))    # a_k
    brs = []
    for r in range(1, rmax + 1):
        B = sum(math.exp(r * lr) * math.cos(r * ang) for lr, ang in small_a if r * lr > -40.0)
        A = sum(math.exp(r * lr) * math.cos(r * ang) for lr, ang in large_a if r * lr > -40.0)
        brs.append(B - A)
    return brs, kap

def b_r_main(r, c, j, kap):
    gam = r * kap
    if gam > 40.0: return 0.0
    th = 2.0 * math.pi * ((r * j) % c) / c
    z = math.exp(-gam) * complex(math.cos(th), math.sin(th))
    return 4.0 * r * (kap / c) * (z * (1.0 + z) / ((1.0 - z) ** 3)).real

print("=== 1. EXACT FINITE-n sum |b_r|/r vs sum |b_r^main|/r ON G1 (M=10, eps=0.45, K_G2=2.2) ===")
for c, j, kt in [
    (6, 7, 1.50), (6, 7, 2.30),
    (10, 13, 0.02), (10, 13, 0.30), (10, 13, 1.80),
    (12, 17, 0.50), (12, 17, 1.50),
    (15, 16, 0.50), (15, 16, 1.40),
    (20, 21, 0.50), (20, 21, 1.20),
    (29, 39, 0.05), (29, 39, 0.45), (29, 39, 1.00),
    (50, 69, 0.15), (50, 69, 0.80),
    (100, 137, 0.10), (100, 137, 0.60),
    (200, 273, 0.08), (200, 273, 0.42),
]:
    if c * kt * kt >= 41.0 or kt >= 2.45: continue
    n = max(132, int(round(2.0 * (math.pi ** 2) * j * (j + c) / (c * kt))) // 2 * 2)
    d = math.gcd(c, j); q = c // d
    rmax = min(400, int(18.0 / kt) + 15)
    brs, kap = exact_finite_n_br(n, c, j, rmax)
    S_ex = sum(abs(brs[r - 1]) / r for r in range(1, rmax + 1))
    S_ex_nonq = sum(abs(brs[r - 1]) / r for r in range(1, rmax + 1) if r % q != 0)
    S_main_nonq = sum(abs(b_r_main(r, c, j, kap)) / r for r in range(1, rmax + 1) if r % q != 0)
    S_diff_nonq = sum(abs(brs[r - 1] - b_r_main(r, c, j, kap)) / r for r in range(1, rmax + 1) if r % q != 0)
    S_ex_q = sum(abs(brs[r - 1]) / r for r in range(1, rmax + 1) if r % q == 0)
    print(f"  c={c:3d} j={j:3d} n={n:6d} kap={kap:.3f} ck^2={c*kap*kap:5.2f} | S_ex={S_ex:.4f} (q:{S_ex_q:.4f}, nonq:{S_ex_nonq:.4f}) | main_nonq={S_main_nonq:.4f}, diff_nonq={S_diff_nonq:.4f}")

print("\n=== 2. CHECK EXACT CLOSED-FORM CF UPPER BOUND (ZERO LOOP OVER r!) ===")
# Let's build and test a PURELY ANALYTIC CF UPPER BOUND for sum_{q \nmid r} |b_r^main|/r!
# Idea:
#   Split r >= 1 (q \nmid r) into:
#   (A) Set I (r >= q): blocks r = m*q + t (m >= 1, 1 <= t <= q-1).
#       Every block has 2 copies of D_r = d*v (v = 1..floor(q/2)).
#       As we proved analytically above:
#         sum_{m>=1, v>=1} 2 * min(K_A/(m d v^2 pi^2 * sinc^2), K_B/(m^2 d v pi X), K_C m X^2/(d v^4 pi^4 * sinc^4))
#       <= C_SetI / d (where C_SetI ~ 0.45, INDEPENDENT of X = q^2*kap and q!! Wait: check exact max over X>0!).
#   (B) Set II (1 <= r < q):
#       Let 1 = q_0 <= q_1 < q_2 < ... < q_s = q be the CF convergents of j/c (with D_0 > D_1 > ... > D_{s-1} = d > D_s = 0).
#       Wait: EVERY integer r in [1, q-1] lies in a unique interval [q_k, q_{k+1} - 1] (0 <= k <= s-1)!
#       Inside [q_k, q_{k+1} - 1]:
#       - Multiples r = m * q_k (1 <= m <= a_{k+1}): D_r = m * D_k!
#         Sum over m >= 1 is <= zeta(3) * term_bound(q_k, D_k) (since t(m*q_k) = t(q_k) is CONSTANT in m!).
#       - Non-multiples r in [q_k, q_{k+1} - 1]:
#         Write r = m * q_k + r' (1 <= m <= a_{k+1}, 1 <= r' <= q_k - 1).
#         WAIT A SECOND! For fixed m, as r' runs over 1..q_k - 1:
#         Any two r'1 != r'2 have |r'1 - r'2| < q_k -> ||(r'1 - r'2) alpha|| >= D_{k-1} / c!
#         And EVERY r = m*q_k + r' has ||r*alpha|| >= D_k / c (since r < q_{k+1})!
#         Wait: at most ONE r' can have D_r in [D_k, D_{k-1}) (namely the semiconvergent r' = q_k - q_{k-1}? Actually r = m*q_k + q_{k-1} has D_r = D_{k-1} - m*D_k >= D_{k+1}, wait: if m = a_{k+1}, m*q_k + q_{k-1} = q_{k+1} which is NOT < q_{k+1}!! So for ALL 1 <= m <= a_{k+1} and 1 <= r' < q_k with m*q_k + r' < q_{k+1}, D_r >= D_k ALWAYS!!)!
#         Let's verify that! By Best Approximation Theorem of convergents, for ALL 1 <= r < q_{k+1}, D_r >= D_k!
#         AND among {m*q_k + r' : 1 <= r' < q_k}, because their spacing is >= D_{k-1}, sorting them by D_r gives:
#           1st closest: >= D_k (actually if r' = q_{k-1}, D_r = D_{k-1} - m*D_k >= D_k!)
#           2nd, 3rd: >= D_{k-1} - m*D_k + D_{k-1} ... actually AT MOST TWO points in [v*D_{k-1}, (v+1)*D_{k-1}) for each v >= 1, and AT MOST ONE point in [D_k, D_{k-1})!
#         WAIT A SECOND!! Check how big the sum over ALL semiconvergents r = m*q_k + q_{k-1} (1 <= m <= a_{k+1} - 1) is!
#         Wait: for m = 1..a_{k+1}-1, D_r = D_{k-1} - m*D_k = (a_{k+1} - m)*D_k + D_{k+1} >= (a_{k+1} - m)*D_k!
#         Let u = a_{k+1} - m in {1, ..., a_{k+1} - 1}! Then r = (a_{k+1} - u)*q_k + q_{k-1} = q_{k+1} - u*q_k, and D_r = u*D_k + D_{k+1}!
#         CHECK THE PRODUCT r * D_r^2 FOR THOSE SEMICONVERGENTS!!
#           For u <= a_{k+1}/2: r >= q_{k+1}/2 and D_r >= u*D_k -> r*D_r^2 >= (1/2) q_{k+1} u^2 D_k^2 >= (1/2) a_{k+1} u^2 * (q_k D_k^2)!!
#           For u > a_{k+1}/2: r >= q_k and D_r >= (a_{k+1}/2) D_k -> r*D_r^2 >= (a_{k+1}^2 / 4) * (q_k D_k^2)!!
#         SO THE SUM OF ALL SEMICONVERGENTS AT LEVEL k IS <= (const / a_{k+1}) * (c / (q_k D_k^2)) — EVEN SMALLER THAN THE CONVERGENT ITSELF!!

# Let's numerically check the exact Set I + Convergents + Semiconvergents + Rest breakdown!
def cf_convergents(num, den):
    g = math.gcd(num, den)
    a, b = (num // g) % (den // g), den // g
    if a == 0: return [(0, 1)]
    res = []
    p0, p1, q0, q1 = 0, 1, 1, 0
    while b != 0:
        ai = a // b
        p2, q2 = ai * p1 + p0, ai * q1 + q0
        if q2 >= 1: res.append((p2, q2))
        p0, p1, q0, q1 = p1, p2, q1, q2
        a, b = b, a % b
    return res

# Check max over X > 0 of Set I (r >= q) universal double sum (times d):
# Set I = (2/d) * sum_{m>=1} sum_{v>=1} min(0.07681/(m v^2), 0.3491/(m^2 v X), 0.01540 m X^2 / v^4)
# Wait: include exact sin(pi v / q) >= 2 v / q? Wait: for v/q <= 1/2, sin(pi v/q) >= 2 v / q loses (pi/2)^2=2.47x on v/q near 1/2!
# Instead keep exact K_A/(m d q * sin^2(pi v / q)) — wait: for v <= q/2, we can just split v <= q/4 (sin >= 2*sqrt(2)*v/q) vs v > q/4 (sin >= 1/sqrt(2))!
max_set1 = 0.0
for ix in range(-40, 60):
    X = 1.2 ** ix
    s = 0.0
    for m in range(1, 200):
        for v in range(1, 200):
            # Using exact small-angle constants (since for large v/q, b_A has 1/(r c sin^2) <= 1/(r c * 0.5) which is tiny for large q!):
            t1 = 0.095 / (m * v * v)
            t2 = 0.388 / (m * m * v * X)
            t3 = 0.025 * m * X * X / (v ** 4)
            s += 2.0 * min(t1, t2, t3)
    max_set1 = max(max_set1, s)
print(f"  Universal Set I (r >= q) double sum bound (times d, max over ALL X > 0): {max_set1:.4f} / d !")
