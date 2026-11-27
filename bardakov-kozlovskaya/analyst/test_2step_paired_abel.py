"""
TEST PAIRED 2-STEP SYMMETRIC ABEL SUMMATION ON |B_r - A_r| / r!!
Recall exact 1-step identity:
  s2 * (B_r - A_r) = 2*Re(U(1) - V(1)) + Re sum_{m=1}^{c-1} d2U(m) w^m - Re sum_{m=1}^{n-c-1} d2V(m) w^m,
where s2 = 4*sin^2(pi*r*nu).
Apply the EXACT SAME identity a second time to d2U and d2V (setting boundary values d2U(0)=d2V(0)=h0,
d2U(c)=conj(h0)*w^{-c}, d2V(n-c)=conj(h0)*w^{-(n-c)} where h0 = (d2U(1)+d2V(1))/2 -> exact cancellation at m=0, M!):
  s2^2 * (B_r - A_r) = s2 * 2*Re(U(1) - V(1))
                     + 2*Re(d2U(1) - d2V(1))  [wait: check exact sign!]
                     + Re sum_{m=1}^{c-1} d4U(m) w^m - Re sum_{m=1}^{n-c-1} d4V(m) w^m!
Test BOTH:
1. Exact identity error (< 1e-12) on finite-n!
2. Effective paired 2-step bound ub2_pair across ALL hard cases (incl. c=4 j=5, c=10 j=11, c=50 j=51, c=200 j=323 kap=0.001)!
"""
import numpy as np
from math import pi, sin, cos, sqrt, gcd

def check_2step_paired_abel(c, j, kap, r):
    # Test on Gaussian-phase model first (exact finite-n is identical up to O(1/n)):
    alpha = kap / c
    nu = j / c
    w = np.exp(-2j * pi * r * nu)
    s2 = 4.0 * (sin(pi * r * nu) ** 2)
    if s2 < 1e-12:
        return None
    ms_e = np.arange(0, c + 1, dtype=np.float64)
    U = np.exp(-r * alpha * ms_e * (c - ms_e)) # real! (stripped e^{-2pi i r m j/c}!)
    M_i = int(c + np.ceil(sqrt(45.0 / max(r * alpha, 1e-12)))) + 20
    ms_i = np.arange(0, M_i + 1, dtype=np.float64)
    V = np.exp(-r * alpha * ms_i * (c + ms_i)) # real!

    A_r = float(np.sum((U[1:c] * (w ** ms_e[1:c])).real))
    B_r = 2.0 * float(np.sum((V[1:M_i] * (w ** ms_i[1:M_i])).real))
    diff_ex = B_r - A_r

    # Extend U symmetrically around c/2 (already symmetric!) and V on [0, M_i]:
    # 1-step on A_r = sum_{m=1}^{c-1} U(m) w^m  and  B_r = 2 * Re sum_{m=1}^{M_i-1} V(m) w^m:
    # For sum_{m=1}^{c-1} U(m) w^m:
    #   s2 * A_r = Re[ U(0)w - U(1) + U(c)w^{c-1} - U(c-1)w^c ] - Re sum_{m=1}^{c-1} d2U(m) w^m
    #            = 2*Re(w - U(1)) - Re sum_{m=1}^{c-1} d2U(m) w^m  (since U(0)=U(c)=1, U(c-1)=U(1), w^c = e^{-2pi i r j} = 1!! Wait: w^c = 1!!)
    # For B_r = 2 * Re sum_{m=1}^{inf} V(m) w^m:
    #   s2 * B_r = 2*Re(w - V(1)) - 2*Re sum_{m=1}^{inf} d2V(m) w^m!
    # Therefore:
    #   s2 * (B_r - A_r) = 2*(U(1) - V(1)) + Re sum_{m=1}^{c-1} d2U(m) w^m - 2*Re sum_{m=1}^{inf} d2V(m) w^m!!
    d2U = U[2:] - 2.0 * U[1:-1] + U[:-2] # indices 1..c-1 (len c-1)
    d2V = V[2:] - 2.0 * V[1:-1] + V[:-2] # indices 1..M_i-1
    rhs1 = 2.0 * (U[1] - V[1]) + float(np.sum((d2U * (w ** ms_e[1:c])).real)) - 2.0 * float(np.sum((d2V * (w ** ms_i[1:M_i])).real))
    err1 = abs(s2 * diff_ex - rhs1)

    # 2nd step! Let h0 = d2V(0)? Wait: extend U(-1)=exp(-r*alpha*(-1)*(c+1))=exp(+r*alpha*(c+1))? No,
    # just set h_U(0) = h_V(0) = h0 := (d2U(1) + d2V(1))/2 (and h_U(c) = h0)!
    # Check: for h_U on 0..c with h_U(0)=h_U(c)=h0 and h_U(m)=d2U(m) for 1..c-1:
    #   s2 * Re sum_{m=1}^{c-1} d2U(m) w^m = 2*Re(h0*w - d2U(1)) - Re sum_{m=1}^{c-1} d2(h_U)(m) w^m!
    # For h_V on 0..inf with h_V(0)=h0 and h_V(m)=d2V(m) for m>=1:
    #   s2 * 2*Re sum_{m=1}^{inf} d2V(m) w^m = 2*Re(h0*w - d2V(1)) - 2*Re sum_{m=1}^{inf} d2(h_V)(m) w^m!
    # SUBTRACT! Notice 2*Re(h0*w) CANCELS IDENTICALLY!!
    #   s2 * [Re sum d2U w^m - 2*Re sum d2V w^m] = 2*(d2V(1) - d2U(1)) - Re sum d2(h_U) w^m + 2*Re sum d2(h_V) w^m!!
    # WAIT A SECOND!! Check optimal h0! Setting h0 = (d2U(1)+d2V(1))/2 makes d2(h_U)(1) = d2U(2) - 2*d2U(1) + h0,
    # OR setting h0 = 2*d2V(1) - d2V(2) (natural linear extrapolation of d2V to m=0, making d2(h_V)(1) = 0!!)!
    # Let's use h0 = (d2U(1) + d2V(1)) / 2 first (completely symmetric!):
    h0 = 0.5 * (d2U[0] + d2V[0])
    # Better: linear extrapolation of average! h0 = 2*avg(1) - avg(2)!
    if len(d2U) >= 2:
        h0 = (2.0*d2U[0] - d2U[1] + 2.0*d2V[0] - d2V[1]) * 0.5
    hU = np.concatenate([[h0], d2U, [h0]])
    hV = np.concatenate([[h0], d2V, [0.0]])
    d4U = hU[2:] - 2.0 * hU[1:-1] + hU[:-2] # len c-1
    d4V = hV[2:] - 2.0 * hV[1:-1] + hV[:-2] # len M_i-1
    rhs2 = s2 * 2.0 * (U[1] - V[1]) + 2.0 * (d2V[0] - d2U[0]) - float(np.sum((d4U * (w ** ms_e[1:c])).real)) + 2.0 * float(np.sum((d4V * (w ** ms_i[1:M_i])).real))
    err2 = abs((s2 * s2) * diff_ex - rhs2)

    # Paired triangle bound from 2-step!
    mh = (c - 1) // 2
    pair_d4 = 2.0 * float(np.sum(np.abs(d4U[:mh] - d4V[:mh])))
    mid_d4 = float(abs(d4U[mh])) if (c % 2 == 0) else 0.0
    tail_d4 = 2.0 * float(np.sum(np.abs(d4V[mh:])))
    ub2_pair = (abs(2.0 * (U[1] - V[1])) / s2 + (abs(2.0 * (d2V[0] - d2U[0])) + pair_d4 + mid_d4 + tail_d4) / (s2 * s2)) / r

    # Compare with 1-step paired bound:
    pair_d2 = 2.0 * float(np.sum(np.abs(d2U[:mh] - d2V[:mh])))
    mid_d2 = float(abs(d2U[mh])) if (c % 2 == 0) else 0.0
    tail_d2 = 2.0 * float(np.sum(np.abs(d2V[mh:])))
    ub1_pair = (abs(2.0 * (U[1] - V[1])) + pair_d2 + mid_d2 + tail_d2) / (s2 * r)
    return abs(diff_ex) / r, err1, err2, ub1_pair, ub2_pair

print("=== TEST 2-STEP PAIRED SYMMETRIC ABEL BOUND ON HARD CREVICE CASES ===")
for c, j, kap in [
    (4, 5, 0.20), (4, 5, 0.05),
    (10, 11, 0.03), (10, 11, 0.10), (10, 11, 0.30),
    (50, 51, 0.03), (50, 51, 0.10),
    (200, 323, 0.003), (200, 323, 0.001),
    (200, 205, 0.01), (200, 202, 0.01), (200, 201, 0.01),
]:
    sum1 = 0.0; sum2 = 0.0; sum_ex = 0.0
    max_err2 = 0.0
    for r in range(1, 30):
        res = check_2step_paired_abel(c, j, kap, r)
        if res is None: continue
        ex, e1, e2, u1, u2 = res
        max_err2 = max(max_err2, e2)
        sum_ex += ex; sum1 += u1; sum2 += min(u1, u2)
    print(f"c={c:3d} j={j:3d} kap={kap:5.3f}: sum_ex(r<30)={sum_ex:.4f} | sum_1step={sum1:.4f} | sum_min(1,2)={sum2:.4f} | err2={max_err2:.1e}")
