"""Test the AMAZING Herglotz-sum mechanism for Q_{k,q}(t, y) = sum_i q^{i(i-1)/2} htilde_{k,i}(t) y^i:
Recall EXACT recurrence:
  Q_{k,q}(t, y) = D_{k+1}( D_{k+1} Q_{k-1,q}(t, y) + y * Q_{k-1,q}(t, q*y) ).
Check on H+ x H+ (Im(t) > 0, Im(y) > 0):
  Term 1: Im( D_{k+1} Q_{k-1,q}(t,y) / Q_{k-1,q}(t,y) )
  Term 2: Im( y * Q_{k-1,q}(t, q*y) / Q_{k-1,q}(t,y) )
  Sum:    Im( Term 1 + Term 2 )
AND ALSO test in the x-variable!
AND ALSO test if Q_{k,q}(t, y) EVER fails real stability up to k=14 for q = 1/2, 1/4, 1/10!
"""
import math
import numpy as np

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_Dtilde(j, poly):
    d = len(poly) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(poly):
        res[m] += c * (m + 1)
        res[m + 1] += c * (m + j + 1)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def build_htilde(kmax):
    p = [[None] * (k + 1) for k in range(kmax + 1)]
    p[1][0] = [1, 2]
    p[1][1] = [1]
    for k in range(2, kmax + 1):
        p[k][0] = apply_Dtilde(k, apply_Dtilde(k, p[k - 1][0]))
        for i in range(1, k):
            p[k][i] = apply_Dtilde(k, poly_add(apply_Dtilde(k, p[k - 1][i]), p[k - 1][i - 1]))
        p[k][k] = apply_Dtilde(k, p[k - 1][k - 1])
    return p

H = build_htilde(10)

def eval_Q_and_DQ(k_idx, q, t, y, j_op):
    """Evaluate Q_{k_idx, q}(t, y), Q_{k_idx, q}(t, q*y), and D_{j_op} Q_{k_idx, q}(t, y)."""
    # Compute htilde_{k, i}(t) and htilde'_{k, i}(t)
    Q = 0j; Qt = 0j; Q_qy = 0j
    for i in range(k_idx + 1):
        w = q ** (i * (i - 1) // 2)
        val = 0j; dval = 0j
        for c in reversed(H[k_idx][i]):
            dval = dval * t + val
            val = val * t + float(c)
        Q += w * val * (y ** i)
        Qt += w * dval * (y ** i)
        Q_qy += w * val * ((q * y) ** i)
    DQ = t * (1.0 + t) * Qt + (1.0 + j_op * t) * Q
    return Q, DQ, Q_qy

rng = np.random.default_rng(2026)
print("=== CHECK HERGLOTZ TERMS ON H+ x H+ FOR Q_{k-1, q}(t, y) ===")
for q in [0.5, 0.25, 0.1]:
    print(f"\n--- q = {q} ---")
    for k in range(2, 8):
        min_im1 = 1e300; min_im2 = 1e300; min_im_sum = 1e300
        worst_sum = None
        for _ in range(3000):
            # sample t in H+ (focusing on near (-1,0) and general H+)
            if rng.random() < 0.5:
                t = complex(rng.uniform(-1.5, 0.5), 10 ** rng.uniform(-4, 0.5))
            else:
                t = complex(rng.normal(0, 2), 10 ** rng.uniform(-3, 1.5))
            # sample y in H+ spanning scale 1 to q^{-k}
            exp_y = rng.uniform(-1, k + 1)
            mag_y = (1.0 / q) ** exp_y
            ang_y = rng.uniform(1e-4, math.pi - 1e-4)
            y = mag_y * complex(math.cos(ang_y), math.sin(ang_y))
            Q, DQ, Q_qy = eval_Q_and_DQ(k - 1, q, t, y, k + 1)
            if abs(Q) < 1e-30: continue
            r1 = DQ / Q
            r2 = y * Q_qy / Q
            s = r1 + r2
            if r1.imag < min_im1: min_im1 = r1.imag
            if r2.imag < min_im2: min_im2 = r2.imag
            if s.imag < min_im_sum:
                min_im_sum = s.imag
                worst_sum = (t, y, r1, r2)
        print(f"  k={k}: min Im(DQ/Q)={min_im1:+.3e}, min Im(y Q(qy)/Q)={min_im2:+.3e}, min Im(sum)={min_im_sum:+.3e}")
        if min_im_sum < 0:
            t, y, r1, r2 = worst_sum
            print(f"       NEG SUM at t={t}, y={y}: r1={r1}, r2={r2}")
