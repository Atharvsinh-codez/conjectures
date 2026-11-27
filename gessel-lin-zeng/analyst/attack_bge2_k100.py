import mpmath as mp

# Targeted attack on b_j >= 2, a_j = 1, w_j = 1 (BOTH single D_{b_j}+y and double D_{b_j}(D_{b_j}+y))!
# Take b = [2]*38 (so at k=38, J_raw(0) = (8/9)*(20/19)^2 = 0.984 < 1!).
# Now hit step 39 with b_39 = 2^m for m = 1..30 (targeting t0 = -1/b_39 where Phi* = 1 - b_39*eps is O(1) vs s_2*eps = O(1))!
# Check:
#   1. Does plain Turan h_{39, i}(t)^2 >= h_{39, i-1}(t) h_{39, i+1}(t) STILL hold on (-1, 0)?
#   2. Does Q_{39, q} (q = 0.5, 0.25) Rayleigh STILL hold?
#   3. What is min_Jraw after step 39?
# Wait: computing 38 steps of [2]*38 in single-operator is INSTANT because D_2 commutes with itself:
#   h_{38, i}(t) = C(38, i) * D_2^{38-i}(1) !!
# We only need 38 applications of D_2 to [1] (takes 0.01 seconds!!)!

mp.mp.dps = 100

def apply_D(b, p):
    res = [mp.mpf(0)]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+b)
    return res

def add_mp(a, b):
    res = [mp.mpf(0)]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v

def deriv_mp(p):
    return [r*p[r] for r in range(1, len(p))] if len(p) > 1 else [mp.mpf(0)]

def eval_W(f, g, t):
    return eval_mp(f, t)*eval_mp(deriv_mp(g), t) - eval_mp(deriv_mp(f), t)*eval_mp(g, t)

# Precompute D_2^m(1) for m = 0..60!
K = 50  # At K=50, i=25: J_raw(0) = (8/9)*(26/25)^2 = 0.9614 < 1!
D2pow = [[mp.mpf(1)]]
for m in range(1, K + 2):
    D2pow.append(apply_D(mp.mpf(2), D2pow[-1]))

# Build P50[i] = C(K, i) * D_2^{K-i}(1)
import math
P50 = [[mp.mpf(math.comb(K, i)) * c for c in D2pow[K - i]] for i in range(K + 1)]

i0 = K // 2  # 25
print(f"=== TARGETED ATTACK ON b = [2]*{K} + [b_next] (SINGLE-OPERATOR D_b + y) AT i = {i0} ===")
# First check exact J_raw(0) at K=50, i=25:
w1 = eval_W(P50[i0], P50[i0-1], mp.mpf(0))
w2 = eval_W(P50[i0+1], P50[i0], mp.mpf(0))
w12 = eval_W(P50[i0+1], P50[i0-1], mp.mpf(0))
print(f"  At k={K}, i={i0}, t=0: J_raw(0) = {float(4*w1*w2/(w12*w12)):.6f} (exact formula = {8/9*(i0+1)*(K-i0+1)/(i0*(K-i0)):.6f})")

# Find where Turan of (D_b P50[i]) can go negative, AND check whether (D_b P50[i] + P50[i-1]) goes negative!
# Recall at t = -eps: P50[i](-eps) ~ C(K, i) * (1 - s_{K-i} eps), where s_m = 2*(2^m - 1) ~ 2^{m+1}.
# For i0 = 25, K - i0 = 25, so s_{25} ~ 2^{26} ~ 6.7e7!
# Try eps = alpha / 2^{26} and b_next = beta * 2^{26}!
s_scale = mp.power(2, K - i0 + 1)
print(f"  Slope scale s_{{K-i0}} ~ {float(s_scale):.3e}")

min_turan_Dp = (mp.mpf('1e100'), None)
min_turan_full = (mp.mpf('1e100'), None)
min_Jraw_next = (mp.mpf('1e100'), None)

for beta in [0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 20.0, 100.0, 1000.0]:
    b_next = mp.mpf(str(beta)) * s_scale
    # Compute D_b P50[r] and full P51[r] = D_b P50[r] + P50[r-1] for r in {i0-1, i0, i0+1}
    Dp = {r: apply_D(b_next, P50[r]) for r in [i0-1, i0, i0+1]}
    F  = {r: add_mp(Dp[r], P50[r-1]) for r in [i0-1, i0, i0+1]}
    for alpha in [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 5.0, 10.0]:
        tv = -mp.mpf(str(alpha)) / b_next  # so b_next * tv = -alpha!
        if tv <= -1: continue
        # Check Dp Turan
        u0, u1, u2 = eval_mp(Dp[i0-1], tv), eval_mp(Dp[i0], tv), eval_mp(Dp[i0+1], tv)
        rel_Dp = (u1*u1 - u0*u2) / (u1*u1 + abs(u0*u2) + mp.mpf('1e-300'))
        if rel_Dp < min_turan_Dp[0]: min_turan_Dp = (rel_Dp, (beta, alpha, float(tv)))
        # Check full F Turan
        v0, v1, v2 = eval_mp(F[i0-1], tv), eval_mp(F[i0], tv), eval_mp(F[i0+1], tv)
        rel_F = (v1*v1 - v0*v2) / (v1*v1 + abs(v0*v2) + mp.mpf('1e-300'))
        if rel_F < min_turan_full[0]: min_turan_full = (rel_F, (beta, alpha, float(tv)))
        # Check J_raw of F
        wf1 = eval_W(F[i0], F[i0-1], tv)
        wf2 = eval_W(F[i0+1], F[i0], tv)
        wf12 = eval_W(F[i0+1], F[i0-1], tv)
        if wf1 <= 0 or wf2 <= 0:
            min_Jraw_next = (mp.mpf(-1), (beta, alpha, float(tv), "W<=0"))
        elif wf12 != 0:
            J = 4*wf1*wf2/(wf12*wf12)
            if J < min_Jraw_next[0]: min_Jraw_next = (J, (beta, alpha, float(tv)))

print(f"  Min rel Turan of D_b(P50) alone:      {float(min_turan_Dp[0]):+.6f} at (beta,alpha,t)={min_turan_Dp[1]}")
print(f"  Min rel Turan of D_b(P50) + P50 full: {float(min_turan_full[0]):+.6f} at (beta,alpha,t)={min_turan_full[1]}")
print(f"  Min J_raw of P51 = D_b(P50) + P50:    {float(min_Jraw_next[0]):+.6f} at (beta,alpha,t)={min_Jraw_next[1]}")

# Also test K=200! At K=200, i=100, J_raw(0) = (8/9)*(101/100)^2 = 0.9067 < 1!
# Wait: can we compute D_2^m(1) at t = -eps directly up to m = 100? Yes, deg 100 takes 0.05s!
K2 = 100
D2pow2 = [[mp.mpf(1)]]
for m in range(1, K2 + 2):
    D2pow2.append(apply_D(mp.mpf(2), D2pow2[-1]))
i2 = 50
P100 = {r: [mp.mpf(math.comb(K2, r)) * c for c in D2pow2[K2 - r]] for r in [i2-2, i2-1, i2, i2+1]}
s2 = mp.power(2, K2 - i2 + 1)
min_Dp2 = mp.mpf('1e100'); min_F2 = mp.mpf('1e100'); min_J2 = mp.mpf('1e100')
loc_F2 = None; loc_J2 = None
for beta in [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 64.0]:
    b_next = mp.mpf(str(beta)) * s2
    Dp = {r: apply_D(b_next, P100[r]) for r in [i2-1, i2, i2+1]}
    F  = {r: add_mp(Dp[r], P100[r-1]) for r in [i2-1, i2, i2+1]}
    for alpha in [0.05, 0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0]:
        tv = -mp.mpf(str(alpha)) / b_next
        u0, u1, u2 = eval_mp(Dp[i2-1], tv), eval_mp(Dp[i2], tv), eval_mp(Dp[i2+1], tv)
        rel_Dp = (u1*u1 - u0*u2) / (u1*u1 + abs(u0*u2) + mp.mpf('1e-300'))
        if rel_Dp < min_Dp2: min_Dp2 = rel_Dp
        v0, v1, v2 = eval_mp(F[i2-1], tv), eval_mp(F[i2], tv), eval_mp(F[i2+1], tv)
        rel_F = (v1*v1 - v0*v2) / (v1*v1 + abs(v0*v2) + mp.mpf('1e-300'))
        if rel_F < min_F2: min_F2 = rel_F; loc_F2 = (beta, alpha)
        wf1 = eval_W(F[i2], F[i2-1], tv)
        wf2 = eval_W(F[i2+1], F[i2], tv)
        wf12 = eval_W(F[i2+1], F[i2-1], tv)
        if wf1 <= 0 or wf2 <= 0:
            min_J2 = mp.mpf(-1); loc_J2 = (beta, alpha, "W<=0")
        elif wf12 != 0:
            J = 4*wf1*wf2/(wf12*wf12)
            if J < min_J2: min_J2 = J; loc_J2 = (beta, alpha)

print(f"\n=== TARGETED ATTACK AT K=100, i=50 (J_raw(0) = {8/9*(51/50)**2:.5f}) ===")
print(f"  Min rel Turan of D_b(P100) alone:      {float(min_Dp2):+.6f}")
print(f"  Min rel Turan of D_b(P100)+P100 full:  {float(min_F2):+.6f} at {loc_F2}")
print(f"  Min J_raw of P101 = D_b(P100)+P100:    {float(min_J2):+.6f} at {loc_J2}")
