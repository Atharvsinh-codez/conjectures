"""
Test the Two-Regime Split at t = -1/(k+2) for GLZ plain Turan preservation:
  (D_m f)^2 - (D_m g)(D_m h) = A * Phi^2 + B * Phi + C
where Phi = 1 + m*t, u = t*(1+t),
  A = f^2 - g*h,
  B = u * A',
  C = u^2 * (f'^2 - g'*h'),
  Disc_Phi = B^2 - 4*A*C = u^2 * (W(g,h)^2 - 4*W(g,f)*W(f,h)).

Check BOTH:
1. For htilde_{k,i} with m = k+2 (the first D_{k+2} in the step k -> k+1):
   signs of A, B, C, J_raw = 4*W(g,f)*W(f,h)/W(g,h)^2 on:
   - Regime I:  (-1, -1/(k+2)]
   - Regime II: [-1/(k+2), 0)
   Also locate exact t_*(k,i) where J_raw(t) = 1, and where B(t) changes sign!
2. For F_{k,i} = D_{k+2} htilde_{k,i} + htilde_{k,i-1} with m = k+2 (the second D_{k+2}):
   signs of A, B, C, J_raw on Regime I and Regime II!
3. Also check if J_raw >= 1 on (-1, -1/m] is preserved by D_m on general triples
   (adversarial check!).
"""
import mpmath as mp
import numpy as np

mp.mp.dps = 100

def poly_eval_deriv(c, t):
    # c[0] + c[1]*t + ... + c[d]*t^d
    val = mp.mpf(0)
    dval = mp.mpf(0)
    for coeff in reversed(c):
        dval = dval * t + val
        val = val * t + coeff
    return val, dval

def poly_D(m, c):
    # D_m(sum c_j t^j) = sum [(j+1)*c_j t^j + (j+m)*c_j t^{j+1}]
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

def build_htilde_and_F(K):
    ht = {(1, 0): [mp.mpf(1), mp.mpf(2)], (1, 1): [mp.mpf(1)]}
    F = {}
    for k in range(2, K + 1):
        m = k + 1
        for i in range(k + 1):
            t1 = poly_D(m, ht[(k - 1, i)]) if i <= k - 1 else [mp.mpf(0)]
            t2 = ht[(k - 1, i - 1)] if i >= 1 else [mp.mpf(0)]
            F[(k, i)] = poly_add(t1, t2)
            ht[(k, i)] = poly_D(m, F[(k, i)])
    return ht, F

def analyze_triple(g, f, h, m, label=""):
    t_split = mp.mpf(-1) / mp.mpf(m)
    # Dense log-spaced grids in Regime I (-1, t_split] and Regime II [t_split, 0)
    # Near -1, bulk, near t_split, near 0
    eps_list = [mp.mpf(10)**(-p) for p in [1, 2, 3, 4, 6, 8, 12]]
    grid_I = []
    for e in eps_list:
        if -1 + e < t_split:
            grid_I.append(-1 + e)
    for j in range(1, 60):
        grid_I.append(mp.mpf(-1) + (t_split + 1) * mp.mpf(j) / 60)
    grid_I.append(t_split)

    grid_II = [t_split]
    for j in range(1, 60):
        grid_II.append(t_split * mp.mpf(60 - j) / 60)
    for e in eps_list:
        if -e > t_split:
            grid_II.append(-e)

    def eval_stats(grid):
        min_A = mp.inf
        min_B = mp.inf
        max_B = -mp.inf
        min_C = mp.inf
        min_Jraw = mp.inf
        min_DmTur = mp.inf
        min_B_where_Jlt1 = mp.inf
        min_APhi2_plus_C_where_Bneg = mp.inf
        min_ratio_where_Bneg = mp.inf # (A*Phi^2 + C) / (-B*Phi) when B*Phi < 0
        t_minJ = None
        t_minB = None
        for t in grid:
            gv, gd = poly_eval_deriv(g, t)
            fv, fd = poly_eval_deriv(f, t)
            hv, hd = poly_eval_deriv(h, t)
            Phi = 1 + m * t
            u = t * (1 + t)
            A = fv*fv - gv*hv
            Ap = 2*fv*fd - gd*hv - gv*hd
            B = u * Ap
            C = u*u * (fd*fd - gd*hd)
            DmTur = A*Phi*Phi + B*Phi + C
            Wgf = gv*fd - gd*fv
            Wfh = fv*hd - fd*hv
            Wgh = gv*hd - gd*hv
            # Scale A, B, C by A_scale = fv^2 + |gv*hv| + 1e-300 to avoid underflow noise
            sc = fv*fv + abs(gv*hv)
            if sc == 0:
                continue
            A_n = A / sc
            B_n = B / sc
            C_n = C / sc
            Dm_n = DmTur / (sc * (Phi*Phi + u*u + 1))
            if Wgh != 0:
                Jraw = 4 * Wgf * Wfh / (Wgh * Wgh)
            else:
                Jraw = mp.inf if Wgf * Wfh > 0 else -mp.inf
            if A_n < min_A: min_A = A_n
            if B_n < min_B: min_B = B_n; t_minB = t
            if B_n > max_B: max_B = B_n
            if C_n < min_C: min_C = C_n
            if Jraw < min_Jraw: min_Jraw = Jraw; t_minJ = t
            if Dm_n < min_DmTur: min_DmTur = Dm_n
            if Jraw < 1 and B_n < min_B_where_Jlt1:
                min_B_where_Jlt1 = B_n
            if B * Phi < 0:
                rat = (A*Phi*Phi + C) / (-B*Phi)
                if rat < min_ratio_where_Bneg:
                    min_ratio_where_Bneg = rat
        return dict(min_A=float(min_A), min_B=float(min_B), max_B=float(max_B),
                    min_C=float(min_C), min_Jraw=float(min_Jraw),
                    min_DmTur=float(min_DmTur),
                    min_B_Jlt1=float(min_B_where_Jlt1),
                    min_rat_Bneg=float(min_ratio_where_Bneg),
                    t_minJ=float(t_minJ) if t_minJ is not None else None,
                    t_minB=float(t_minB) if t_minB is not None else None)

    return eval_stats(grid_I), eval_stats(grid_II)

ht, F = build_htilde_and_F(20)

print("=== PART 1: htilde_{k,i} with m = k+2 ===")
for k in [4, 6, 8, 12, 16, 20]:
    m = k + 2
    worst_I = dict(min_A=1e9, min_B=1e9, min_C=1e9, min_Jraw=1e9)
    worst_II = dict(min_A=1e9, min_B=1e9, min_C=1e9, min_Jraw=1e9, min_B_Jlt1=1e9, min_rat_Bneg=1e9)
    for i in range(1, k):
        sI, sII = analyze_triple(ht[(k, i-1)], ht[(k, i)], ht[(k, i+1)], m)
        for key in worst_I:
            if sI[key] < worst_I[key]: worst_I[key] = sI[key]
        for key in worst_II:
            if sII[key] < worst_II[key]: worst_II[key] = sII[key]
    print(f"k={k:2d} (m={m:2d}) | RegI: minA={worst_I['min_A']:.3e} minB={worst_I['min_B']:.3e} minC={worst_I['min_C']:.3e} minJraw={worst_I['min_Jraw']:.5f}")
    print(f"            | RegII: minA={worst_II['min_A']:.3e} minB={worst_II['min_B']:.3e} minC={worst_II['min_C']:.3e} minJraw={worst_II['min_Jraw']:.5f} minB(J<1)={worst_II['min_B_Jlt1']:.3e} minRat(BPhi<0)={worst_II['min_rat_Bneg']:.5f}")

print("\n=== PART 2: F_{k,i} = D_{k+1} htilde_{k-1,i} + htilde_{k-1,i-1} with m = k+1 ===")
for k in [4, 6, 8, 12, 16, 20]:
    m = k + 1
    worst_I = dict(min_A=1e9, min_B=1e9, min_C=1e9, min_Jraw=1e9)
    worst_II = dict(min_A=1e9, min_B=1e9, min_C=1e9, min_Jraw=1e9, min_B_Jlt1=1e9, min_rat_Bneg=1e9)
    for i in range(1, k):
        sI, sII = analyze_triple(F[(k, i-1)], F[(k, i)], F[(k, i+1)], m)
        for key in worst_I:
            if sI[key] < worst_I[key]: worst_I[key] = sI[key]
        for key in worst_II:
            if sII[key] < worst_II[key]: worst_II[key] = sII[key]
    print(f"k={k:2d} (m={m:2d}) | RegI: minA={worst_I['min_A']:.3e} minB={worst_I['min_B']:.3e} minC={worst_I['min_C']:.3e} minJraw={worst_I['min_Jraw']:.5f}")
    print(f"            | RegII: minA={worst_II['min_A']:.3e} minB={worst_II['min_B']:.3e} minC={worst_II['min_C']:.3e} minJraw={worst_II['min_Jraw']:.5f} minB(J<1)={worst_II['min_B_Jlt1']:.3e} minRat(BPhi<0)={worst_II['min_rat_Bneg']:.5f}")
