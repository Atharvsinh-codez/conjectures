"""
High-k (k=25, 35, 50, 70) check of Two-Regime Split at t = -1/(k+2):
- Exact J_raw(0) and min J_raw on (-1, -1/(k+2)] and [-1/(k+2), 0)
- Wherever J_raw(t) < 1, check:
  (a) sign of B(t) = t(1+t)*A'(t),
  (b) whether t > r_max (rightmost root of htilde_{k,i}),
  (c) signs of A(t), C_unscaled = f'^2 - g'h', and (A*Phi^2 + B*Phi + C).
- Also test if Two-Regime condition (J_raw >= 1 on (-1, -1/m] and A>0, B>=0 where J_raw<1)
  is preserved under D_m on general triples!
"""
import mpmath as mp

mp.mp.dps = 180

def poly_eval_deriv(c, t):
    val = mp.mpf(0)
    dval = mp.mpf(0)
    for coeff in reversed(c):
        dval = dval * t + val
        val = val * t + coeff
    return val, dval

def poly_D(m, c):
    d = len(c) - 1
    out = [mp.mpf(0)] * (d + 2)
    for j, cj in enumerate(c):
        out[j] += (j + 1) * cj
        out[j + 1] += (j + m) * cj
    return out

def poly_add(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for i, v in enumerate(a): out[i] += v
    for i, v in enumerate(b): out[i] += v
    return out

def build_htilde(K):
    ht = {(1, 0): [mp.mpf(1), mp.mpf(2)], (1, 1): [mp.mpf(1)]}
    for k in range(2, K + 1):
        m = k + 1
        for i in range(k + 1):
            t1 = poly_D(m, ht[(k - 1, i)]) if i <= k - 1 else [mp.mpf(0)]
            t2 = ht[(k - 1, i - 1)] if i >= 1 else [mp.mpf(0)]
            ht[(k, i)] = poly_D(m, poly_add(t1, t2))
    return ht

ht = build_htilde(70)

print("=== Check J_raw(0) across k = 10..70 ===")
for k in [10, 20, 30, 40, 50, 60, 70]:
    min_J0 = mp.inf
    worst_i = -1
    for i in range(1, k):
        g, f, h = ht[(k, i-1)], ht[(k, i)], ht[(k, i+1)]
        gv, gd = g[0], g[1] if len(g) > 1 else mp.mpf(0)
        fv, fd = f[0], f[1] if len(f) > 1 else mp.mpf(0)
        hv, hd = h[0], h[1] if len(h) > 1 else mp.mpf(0)
        Wgf = gv*fd - gd*fv
        Wfh = fv*hd - fd*hv
        Wgh = gv*hd - gd*hv
        J0 = 4 * Wgf * Wfh / (Wgh * Wgh)
        if J0 < min_J0:
            min_J0 = J0
            worst_i = i
    print(f"k={k:2d}: min J_raw(0) = {float(min_J0):.8f} at i={worst_i}")

print("\n=== Detailed Regime I & II scan for k = 30, 50, 70 ===")
for k in [30, 50, 70]:
    m = k + 2
    t_split = mp.mpf(-1) / mp.mpf(m)
    # Test representative i values including worst_i, i=1, i=k//4, i=k//2, i=3*k//4, i=k-1
    i_list = sorted(set([1, 2, k//4, k//2, 3*k//4, k-2, k-1]))
    # Build dense grid in RegI and RegII, especially near 0 (scales 1/s down to 1e-20)
    grid_I = [-1 + mp.mpf(10)**(-p) for p in [1, 2, 4, 8]] + \
             [-1 + (t_split + 1) * mp.mpf(j) / 40 for j in range(1, 41)]
    grid_II = [t_split * mp.mpf(40 - j) / 40 for j in range(40)] + \
              [-t_split * mp.mpf(-1) * mp.mpf(10)**(-p) for p in [1, 2, 3, 4, 6, 8, 12]]
    # Fix sign of second part of grid_II:
    grid_II = [t_split * mp.mpf(40 - j) / 40 for j in range(40)] + \
              [t_split * mp.mpf(10)**(-p/2) for p in range(1, 24)] + [mp.mpf(0)]
    min_J_I = mp.inf
    min_J_II = mp.inf
    min_B_where_Jlt1 = mp.inf
    max_B_where_Jlt1 = -mp.inf
    cnt_Jlt1 = 0
    cnt_Jlt1_Bneg = 0
    worst_Jlt1_info = None
    min_Cunsc_II = mp.inf
    for i in i_list:
        g, f, h = ht[(k, i-1)], ht[(k, i)], ht[(k, i+1)]
        # Also locate approximate rightmost root scale via -f(0)/f'(0)
        s_f = f[1]/f[0] if len(f) > 1 else mp.mpf(1)
        extra_pts = [-mp.mpf(alpha)/s_f for alpha in [0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0] if -mp.mpf(alpha)/s_f > t_split]
        for t in grid_I:
            gv, gd = poly_eval_deriv(g, t)
            fv, fd = poly_eval_deriv(f, t)
            hv, hd = poly_eval_deriv(h, t)
            Wgf = gv*fd - gd*fv
            Wfh = fv*hd - fd*hv
            Wgh = gv*hd - gd*hv
            Jraw = 4*Wgf*Wfh/(Wgh*Wgh) if Wgh != 0 else mp.inf
            if Jraw < min_J_I: min_J_I = Jraw
        for t in grid_II + extra_pts:
            gv, gd = poly_eval_deriv(g, t)
            fv, fd = poly_eval_deriv(f, t)
            hv, hd = poly_eval_deriv(h, t)
            Wgf = gv*fd - gd*fv
            Wfh = fv*hd - fd*hv
            Wgh = gv*hd - gd*hv
            Jraw = 4*Wgf*Wfh/(Wgh*Wgh) if Wgh != 0 else mp.inf
            if Jraw < min_J_II: min_J_II = Jraw
            A = fv*fv - gv*hv
            Ap = 2*fv*fd - gd*hv - gv*hd
            u = t*(1+t)
            B = u * Ap
            Cunsc = fd*fd - gd*hd
            sc = fv*fv + abs(gv*hv)
            if sc > 0:
                if Cunsc / (fd*fd + abs(gd*hd) + mp.mpf('1e-300')) < min_Cunsc_II:
                    min_Cunsc_II = Cunsc / (fd*fd + abs(gd*hd) + mp.mpf('1e-300'))
            if Jraw < 1:
                cnt_Jlt1 += 1
                Bn = B / sc if sc > 0 else mp.mpf(0)
                if Bn < min_B_where_Jlt1:
                    min_B_where_Jlt1 = Bn
                    worst_Jlt1_info = (i, float(t), float(t*s_f), float(Jraw), float(Bn), float(fv/f[0]))
                if Bn > max_B_where_Jlt1:
                    max_B_where_Jlt1 = Bn
                if B < 0:
                    cnt_Jlt1_Bneg += 1
    print(f"k={k:2d}: RegI minJraw={float(min_J_I):.6f} | RegII minJraw={float(min_J_II):.6f}, minC_unsc={float(min_Cunsc_II):.4e}")
    print(f"      J<1 pts={cnt_Jlt1}, J<1 & B<0 pts={cnt_Jlt1_Bneg}, minB(J<1)={float(min_B_where_Jlt1) if cnt_Jlt1 else 'none'}, worst={worst_Jlt1_info}")
