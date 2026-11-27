import random
import math
import numpy as np

def poly_mul(a, b):
    out = np.zeros(len(a) + len(b) - 1, dtype=float)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out

def add_np(a, b):
    n = max(len(a), len(b))
    out = np.zeros(n, dtype=float)
    out[:len(a)] += a
    out[:len(b)] += b
    return out

def sub_np(a, b):
    n = max(len(a), len(b))
    out = np.zeros(n, dtype=float)
    out[:len(a)] += a
    out[:len(b)] -= b
    return out

def apply_T(a, b, p):
    # T(p) = t(1+t)p' + (a + b*t)p
    d = len(p) - 1
    out = np.zeros(d + 2, dtype=float)
    for m, cm in enumerate(p):
        out[m] += (m + a) * cm
        out[m + 1] += (m + b) * cm
    return out

def eval_poly(p, x):
    acc = np.zeros_like(x, dtype=float)
    for c in reversed(p):
        acc = acc * x + c
    return acc

def roots_real(p):
    # p is low->high
    r = np.roots(p[::-1])
    if np.max(np.abs(r.imag)) > 1e-6:
        return None
    return np.sort(r.real)

def interlaces(p, q):
    rp = roots_real(p)
    rq = roots_real(q)
    if rp is None or rq is None:
        return False
    dp, dq = len(rp), len(rq)
    if dq == dp + 1:
        return all(rq[m] < rp[m] < rq[m+1] for m in range(dp))
    elif dq == dp:
        return all(rp[m] < rq[m] for m in range(dp)) and all(rq[m] < rp[m+1] for m in range(dp - 1))
    return False

def poly_from_roots(rs):
    p = np.array([1.0])
    for r in rs:
        p = poly_mul(p, np.array([-r, 1.0]))
    return p

# Test the 1-step conjecture on random triples p2 << p1 << p0 in (-1, 0) satisfying p1^2 > c * p0 * p2 on (-1, 0) (or all R)!
t_grid_int = np.linspace(-1.0, 0.0, 501)
t_grid_R = np.concatenate([np.linspace(-20.0, -1.0, 200), t_grid_int, np.linspace(0.0, 20.0, 200)])

random.seed(2026)

for domain_name, grid in [("on (-1, 0)", t_grid_int), ("on all R", t_grid_R)]:
    print(f"\n=== 1-STEP CONJECTURE TEST: p2 << p1 << p0 with p1^2 >= C * p0*p2 {domain_name} ===")
    for C_req in [1.0, 1.5, 2.0, 4.0]:
        tested = 0
        fail_q2_q1 = 0
        fail_m2_xy_int = 0
        fail_turan_out_int = 0
        fail_turan_out_R = 0
        ex_q2_q1 = None
        for trial in range(400):
            d = random.choice([0, 1, 2, 3])
            # Generate p2 (deg d) << p1 (deg d+1) << p0 (deg d+2) in (-1, 0)
            # First generate p0 (d+2 roots in (-1, 0))
            r0 = sorted([-random.uniform(0.01, 0.99) for _ in range(d + 2)])
            # p1 has d+1 roots, one in each gap of r0
            r1 = [random.uniform(r0[m] + 1e-4, r0[m+1] - 1e-4) for m in range(d + 1)]
            # p2 has d roots, one in each gap of r1
            r2 = [random.uniform(r1[m] + 1e-4, r1[m+1] - 1e-4) for m in range(d)]
            p0 = poly_from_roots(r0)
            p1 = poly_from_roots(r1)
            p2 = poly_from_roots(r2)
            # Compute max ratio (p0*p2)/(p1^2) over grid where p0*p2 > 0!
            v0 = eval_poly(p0, grid)
            v1 = eval_poly(p1, grid)
            v2 = eval_poly(p2, grid)
            # Find max alpha such that p1^2 - alpha * p0 * p2 >= 0 on grid:
            # Wait: p0*p2 can be positive where p1 is small? If p1(t)=0, then r1[m] lies strictly between r0[m] and r0[m+1],
            # and between r2[m-1] and r2[m], check sign of p0(r1[m])*p2(r1[m])!
            # At r1[d] (largest root of p1): p0 has 1 root > r1[d] (neg), p2 has 0 roots > r1[d] (pos) -> p0*p2 < 0!
            # At r1[d-1]: p0 has 2 roots > r1[d-1] (pos), p2 has 1 root > r1[d-1] (neg) -> p0*p2 < 0!
            # SO p0(r1[m]) * p2(r1[m]) < 0 AT EVERY ROOT OF p1!!
            pos_mask = (v0 * v2) > 0
            if not np.any(pos_mask):
                max_ratio = 1e-12
            else:
                max_ratio = np.max((v0[pos_mask] * v2[pos_mask]) / (v1[pos_mask] ** 2))
            if max_ratio <= 0:
                continue
            # Scale p2 by s2 = (1 / (C_req * max_ratio)) * random factor in (0.1, 0.99)
            # so that p1^2 >= C_req * p0 * p2 holds tightly!
            s2 = random.uniform(0.3, 0.999) / (C_req * max_ratio)
            p2 = p2 * s2
            # Also randomly rescale p1 -> lam*p1, p2 -> lam^2*p2 (invariant under p1^2/(p0*p2)!)
            lam = 10 ** random.uniform(-2, 2)
            p1 = p1 * lam
            p2 = p2 * (lam ** 2)

            # Pick random T = t(1+t)d/dt + (a + b t), 0 < a < b, and w > 0
            a = random.choice([1.0, 2.0, 10.0])
            b = a + random.choice([0.5, 1.0, 5.0, 50.0, 500.0])
            w = 10 ** random.uniform(-2, 2)

            q0 = apply_T(a, b, p0)
            q1 = add_np(apply_T(a, b, p1), w * p0)
            q2 = add_np(apply_T(a, b, p2), w * p1)
            q3 = w * p2

            tested += 1
            if not interlaces(q2, q1):
                fail_q2_q1 += 1
                if ex_q2_q1 is None:
                    ex_q2_q1 = (d, a, b, w, r0, r1, r2, s2, lam)
            # Check M^2 - XY = (T p1)^2 - p0 * (T^2 p2) on (-1, 0)
            M = apply_T(a, b, p1)
            Y = apply_T(a, b, apply_T(a, b, p2))
            m2xy = eval_poly(M, t_grid_int)**2 - eval_poly(p0, t_grid_int) * eval_poly(Y, t_grid_int)
            if np.min(m2xy) < -1e-9:
                fail_m2_xy_int += 1
            # Check output Turan q1^2 - q0*q2 on (-1,0) and on R
            tur_int = eval_poly(q1, t_grid_int)**2 - eval_poly(q0, t_grid_int) * eval_poly(q2, t_grid_int)
            if np.min(tur_int) < -1e-9:
                fail_turan_out_int += 1
            tur_R = eval_poly(q1, t_grid_R)**2 - eval_poly(q0, t_grid_R) * eval_poly(q2, t_grid_R)
            if np.min(tur_R) < -1e-9:
                fail_turan_out_R += 1
        print(f"  C={C_req:3.1f} (n={tested}): fail(q2<<q1)={fail_q2_q1:3d} | fail(M^2>XY on (-1,0))={fail_m2_xy_int:3d} | fail(q1^2>q0q2 on (-1,0))={fail_turan_out_int:3d} | fail(q1^2>q0q2 on R)={fail_turan_out_R:3d}")
        if ex_q2_q1 is not None:
            d, a, b, w, r0, r1, r2, s2, lam = ex_q2_q1
            print(f"    ex fail q2<<q1: d={d}, a={a}, b={b}, w={w:.3g}, r0={[round(z,3) for z in r0]}, r1={[round(z,3) for z in r1]}, r2={[round(z,3) for z in r2]}")
