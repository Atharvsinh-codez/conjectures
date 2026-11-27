"""
Stress-test our 1D CIRCLE INTEGRAL MIRACLE directly from the cubic roots (x(phi), y(phi), z(phi))
of T^3 - T^2 + u_j T - R_v e^{i phi} = 0 at n = 20, 40, 80, 200 across ALL corners of the (c, j) triangle!
Recall:
  s_j^{(c)} = sin(j pi / n) / sin((j + c) pi / n),   j = c .. m - ceil(c/2)
  t_j = (-1)^c * (s_j^{(c)})^n,   u_j = t_j / (1 + t_j)^2,   R_v = |t_j| * (r3/r1)^n / |1 + t_j|^3.
Wait: at odd c with j = m - (c+1)/2 when s_j -> 1, |1 + t_j| = 1 - s_j^n is small -> u_j is large negative,
and R_v = s^n (r3/r1)^n / (1 - s^n)^3 is well-defined!
To avoid any float64 underflow/overflow at n = 80, 200, compute log|P| + i*arg(P) carefully!
"""
import math
import numpy as np

def eval_1d_circle_direct(n, c, j, N=128):
    """Evaluate (-1)^j * P_n(1, u_j, R_v e^{i phi}) / (R_v^c e^{i c phi}) directly from cubic roots!"""
    m = n // 2
    r1 = math.sin((j + c) * math.pi / n)
    r2 = math.sin(j * math.pi / n)
    r3 = math.sin(c * math.pi / n)
    s = r2 / r1
    # Use log-scale for Ry = s^n, Rz = (r3/r1)^n so n=200 never underflows!
    log_Ry = n * math.log(s)
    log_Rz = n * math.log(r3 / r1)
    # Wait: if log_Ry < -600, t_j = 0 in float64 -> work in (x=1, y, z) rescaled!
    # ACTUALLY BETTER: work in the (x, y, z) scaling where e1 = 1 + t_j (so x0 = 1, y0 = t_j, z0 = 0)!
    # Cubic for T_tilde = (1 + t_j) * T:
    #   T_tilde^3 - (1 + t_j) T_tilde^2 + t_j T_tilde - t_j * R_z e^{i phi} = 0!!
    # OH!! THAT RESCALING COMPLETELY ELIMINATES (1 + t_j)^3 AND NEVER OVERFLOWS EVEN WHEN t_j -> -1!!
    # Check: if T_tilde^3 - (1+t) T_tilde^2 + t T_tilde - t*z0 = 0, then
    #   e1 = 1 + t,  e2 = t,  e3 = t * z0 -> u = e2/e1^2 = t/(1+t)^2, v = e3/e1^3 = t*z0/(1+t)^3!
    #   And P_n(1+t, t, t*z0) / (t*z0)^c = (1+t)^{n-3c} * P_n(1, u, v) / v^c!
    # Wait: for odd c, t < 0, so (t*z0)^c has sign (-1)^c, and (1+t) = 1 - s^n > 0 ALWAYS (since s < 1 for odd c!)!
    # For even c, t > 0, so (1+t) > 0 ALWAYS and t^c > 0 ALWAYS!
    # So the sign of P_n(1+t, t, t * R_z e^{i phi}) / (t * R_z e^{i phi})^c is IDENTICAL to P_n(1, u, v) / v^c!
    if log_Ry < -650 or log_Rz < -650:
        return None  # handle extreme tiny radii with mpmath if needed
    Ry = math.exp(log_Ry)
    Rz = math.exp(log_Rz)
    t = Ry if (c % 2 == 0) else -Ry

    phis = np.linspace(-np.pi, np.pi, N, endpoint=False)
    # Solve T^3 - (1+t) T^2 + t T - t * Rz * e^{i phi} = 0 for each phi!
    # Compute log_val = sum_{k=0}^{n-1} log(z - (alpha + zeta^k beta)^n) - c * log(t * Rz * e^{i phi})
    # Better: evaluate P_n(x, y, z) via prod_{k=0}^{n-1} (z - (alpha + zeta^k beta)^n) where |x| >= |y| >= |z|!
    log_vals = np.zeros(N, dtype=complex)
    for idx, phi in enumerate(phis):
        e3_val = t * Rz * np.exp(1j * phi)
        roots = np.roots([1.0, -(1.0 + t), t, -e3_val])
        # Sort by decreasing modulus so |x| >= |y| >= |z|
        roots = roots[np.argsort(-np.abs(roots))]
        xv, yv, zv = roots[0], roots[1], roots[2]
        # Compute prod_{k=0}^{n-1} (zv/xv - (1 + zeta^k (yv/xv)^{1/n})^n) * xv^n
        yr = yv / xv
        zr = zv / xv
        ry = abs(yr) ** (1.0 / n)
        ay = np.angle(yr) / n
        ks = np.arange(n)
        w = 1.0 + ry * np.exp(1j * (ay + 2.0 * np.pi * ks / n))
        # (1 + ...)^n
        wn = np.exp(n * np.log(w))
        factors = zr - wn
        # Sum of logs: n*log(xv) + sum log(factors) - c*log(t*Rz) - i*c*phi
        # Keep track of exact sign of t^c = (-1)^c * Ry^c!
        lv = n * np.log(xv) + np.sum(np.log(factors)) - c * (log_Ry + log_Rz) - 1j * c * phi
        if c % 2 == 1:
            lv += 1j * np.pi  # divide by (-1)^c = -1
        if j % 2 == 1:
            lv += 1j * np.pi  # multiply by (-1)^j
        log_vals[idx] = lv

    # Shift real parts by max(Re) to avoid overflow/underflow
    max_re = np.max(log_vals.real)
    F = np.exp(log_vals - max_re)
    avg = np.mean(F)
    abs_avg = np.mean(np.abs(F))
    ratio = avg.real / abs_avg
    min_re_ratio = np.min(F.real) / avg.real
    im_err = abs(avg.imag) / abs(avg.real)
    return ratio, min_re_ratio, im_err

for n in [20, 40, 80, 200]:
    m = n // 2
    max_c = (n - 2) // 3  # largest c with deg Q_c >= 1
    # Test ALL critical corners + bulk!
    test_pts = []
    # 1. Bottom edge: c = 2, 3 (j = c, mid, j_end)
    for c in [1, 2, 3, 4]:
        if c > max_c: continue
        je = m - (c + 1) // 2
        for j in sorted(set([c, (c + je)//2, je])):
            test_pts.append((f"bot(c={c})", c, j))
    # 2. Top corner: c = max_c, max_c - 1, max_c - 2
    for c in [max_c, max_c - 1, max_c - 2]:
        if c < 1: continue
        je = m - (c + 1) // 2
        for j in sorted(set([c, (c + je)//2, je])):
            test_pts.append((f"top(c={c})", c, j))
    # 3. Bulk & right edge (u = 1/4 for even c, large -u for odd c) & left edge (j = c)
    for c in [max(2, max_c // 2), max(2, max_c // 2 + 1), max(2, (3 * max_c) // 4)]:
        if c > max_c: continue
        je = m - (c + 1) // 2
        for tag, j in [("left", c), ("mid", (c + je)//2), ("right", je)]:
            test_pts.append((f"{tag}(c={c})", c, j))

    # Deduplicate
    seen = set()
    worst = (1e9, None)
    print(f"\n=== n = {n} (m = {m}, max_c = {max_c}) ===")
    for tag, c, j in test_pts:
        if (c, j) in seen or j > m - (c + 1) // 2 or j < c: continue
        seen.add((c, j))
        res = eval_1d_circle_direct(n, c, j, N=128)
        if res is None:
            print(f"  {tag:12s} (c={c:2d}, j={j:2d}): [Ry or Rz < 1e-280, extreme asymptotic!]")
            continue
        ratio, min_re, im_err = res
        if ratio < worst[0]:
            worst = (ratio, (tag, c, j, min_re, im_err))
        if ratio < 0.95 or n <= 40 or "top" in tag or "right" in tag:
            print(f"  {tag:12s} (c={c:2d}, j={j:2d}): Re(avg)/avg|F| = {ratio:+.6f}, minRe/avg = {min_re:+.4f}, |Im/Re|={im_err:.1e}")
    print(f"  --> WORST at n={n}: Re(avg)/avg|F| = {worst[0]:+.6f} at {worst[1]}")
