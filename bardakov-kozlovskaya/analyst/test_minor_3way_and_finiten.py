"""
1. TEST THE 4TH-ORDER PAIRED ABEL IDENTITY ON *EXACT FINITE-n* (NOT JUST MODEL!)
   Recall: at exact finite n:
     a_{-m}^{-1} = ((1 + r_j e^{i(th* + 2m*pi/n)}) / (s_j e^{i j pi/n}))^n
     a_m         = ((s_j e^{i j pi/n}) / (1 + r_j e^{i(th* - 2m*pi/n)}))^n
   With nu_exact := j / c (EXACT RATIONAL! So cos(c*th) = 1 IDENTICALLY, and e^{i*c*th} = 1 IDENTICALLY!!):
     U_r(m) := a_{-m}^{-r} * e^{+2*pi*i*r*m*j/c}
     V_r(m) := a_m^r       * e^{+2*pi*i*r*m*j/c}
   CHECK:
     (a) U_r(0) = V_r(0) = 1!
     (b) U_r(c) = 1 (since a_{-c}^{-1} = 1 and e^{2*pi*i*r*j} = 1)!
     (c) U_r(c - m) = conj(U_r(m))!
     (d) Re(U_r(m)) is EVEN around c/2, and Im(U_r(m)) is ODD around c/2!
     (e) What is Re(B_r - A_r) with th = 2*pi*r*j/c?
         B_r - A_r = Re[ sum_{m=1}^{n-c-1} V_r(m) e^{-i m th} - sum_{m=1}^{c-1} U_r(m) e^{-i m th} ]!
     Wait: check BOTH nu_exact = j/c AND nu_bdy = sin(psi_j)*cos(alpha_1)/sin(theta_c)!
2. TEST THE 3-WAY MINIMUM B_minor_3way ACROSS c = 2..100 IN G* = {c*kap^2 < 42, c*q*kap > 0.85, kap < 2.5}!
   Map EVERY surviving (c, jm, kap) where B_minor_3way >= 1.45!
"""
import numpy as np
from math import pi, sin, cos, log, exp, sqrt, gcd

# =============================================================================
# 1. TEST 4TH-ORDER PAIRED ABEL ON EXACT FINITE-n!
# =============================================================================
def check_finiten_4th_abel(n, c, j, r):
    rj = sin(j * pi / n) / sin((j + c) * pi / n)
    sj = sin(c * pi / n) / sin((j + c) * pi / n)
    th_star = pi * (n - c) / n
    kap = 2.0 * pi * sin(j * pi / n) * sin((j + c) * pi / n) / sin(c * pi / n)

    # Exact rational frequency nu = j / c:
    nu = j / c
    th = 2.0 * pi * r * nu
    s2 = 4.0 * (sin(th / 2.0) ** 2)
    if s2 < 1e-12: return None

    ms_e = np.arange(0, c + 1, dtype=np.float64)
    th_e = th_star + 2.0 * ms_e * pi / n
    # a_{-m}^{-1} = ((1 + rj e^{i th_e}) / (sj e^{i j pi/n}))^n
    log_ae_inv = n * np.log((1.0 + rj * np.exp(1j * th_e)) / (sj * np.exp(1j * j * pi / n)))
    U = np.exp(r * log_ae_inv + 2j * pi * r * nu * ms_e) # U[0] = U[c] = 1, U[c-m] = conj(U[m])!

    M = n - c
    ms_i = np.arange(0, M + 1, dtype=np.float64)
    th_i = th_star - 2.0 * ms_i * pi / n
    log_ai = n * np.log((sj * np.exp(1j * j * pi / n)) / (1.0 + rj * np.exp(1j * th_i)))
    V = np.exp(r * log_ai + 2j * pi * r * nu * ms_i) # V[0] = V[M] = 1, V[M-m] = conj(V[m])!

    # Note: A_r = Re sum_{m=1}^{c-1} U[m] e^{-i m th} = sum_{m=1}^{c-1} (Re U[m] cos(m th) + Im U[m] sin(m th))!
    # Wait: because U[c-m] = conj(U[m]) and cos((c-m)th) = cos(m th), sin((c-m)th) = -sin(m th):
    #   Im U[c-m] * sin((c-m)th) = (-Im U[m]) * (-sin(m th)) = +Im U[m] * sin(m th) (DOES NOT CANCEL!)
    #   SO KEEP COMPLEX FORM!
    #   Let w = e^{-i th}. Note w^c = e^{-2*pi*i*r*j} = 1!
    #   Wait: for COMPLEX sequence X(m) (m = 0..N) with X(0) = X(N) = 1 and w^N = 1:
    #     (2 - w - w^{-1}) * sum_{m=1}^{N-1} X(m) w^m
    #       = (w + w^{-1}) - (X(1) + X(N-1)) - sum_{m=1}^{N-1} (d2X)(m) w^m  (EXACT!)
    #     Notice X(N-1) = conj(X(1)) => X(1) + X(N-1) = 2 Re(X(1))!
    #     For V (length M = n - c): w^M = e^{-2*pi*i*r*j*(n-c)/c} != 1!
    #     WAIT A SECOND! For V, m = 1..M/2 is on one side and m = M-1..M/2 is the conjugate side (a_{-c-m})!
    #     And at m' = M - m = n - c - m:
    #       a_{n-c-m} = conj(a_m), so sum_{m=1}^{n-c-1} a_m^r = 2 * Re sum_{m=1}^{floor((n-c-1)/2)} a_m^r + [mid]!!
    #     SO B_r = 2 * Re sum_{m=1}^{M_h} V[m] e^{-i m th}  (with V[M_h] exponentially tiny ~ exp(-n)!!)!
    #     AND A_r = Re sum_{m=1}^{c-1} U[m] e^{-i m th}!!
    exact_A = float(np.sum(U[1:c] * np.exp(-1j * ms_e[1:c] * th)).real)
    exact_B = float(np.sum(V[1:M] * np.exp(-1j * ms_i[1:M] * th)).real)
    exact = exact_B - exact_A

    # Take REAL parts of (2*V(m) e^{-i m th}) = 2*(Re V cos(m th) + Im V sin(m th)):
    # Even cleaner! Apply (2 - w - w^{-1}) TWICE to the COMPLEX sums (using Re[...] at the end)!
    # For U (m=0..c, w^c = 1, U(c-m)=conj(U(m))):
    #   Define U(-1) := conj(U(1)) = U(c-1) and U(c+1) := U(1)!  [Periodic extension mod c!]
    #   WAIT A SECOND!! STOP EVERYTHING!!
    #   IF YOU EXTEND U PERIODICALLY MOD c: U(m+c) := U(m), THEN U(-1) = U(c-1) = conj(U(1))!!
    #   AND h_U(0) := U(1) - 2*U(0) + U(-1) = 2*Re(U(1)) - 2  (PURELY REAL!!)!
    #   AND FOR V (for 0 <= m <= M/2): DEFINE V(-1) := conj(V(1))!
    #   THEN h_V(0) := V(1) - 2*V(0) + V(-1) = 2*Re(V(1)) - 2  (PURELY REAL!!)!
    #   Now check the EXACT 4th-order identity for Re sum_{m=1}^{c-1} U(m) w^m AND 2 Re sum_{m=1}^{Mh} V(m) w^m:
    #   Step 1: s2 * Re sum_{m=1}^{c-1} U(m) w^m = 2*(cos th - Re U(1)) - Re sum_{m=1}^{c-1} d2U(m) w^m!
    #   Wait: check if s2 * Re sum_{m=1}^{Mh} V(m) w^m ALSO has boundary (cos th - Re V(1))!
    #   Check: (2 - w - w^{-1}) sum_{m=1}^{Mh} V(m) w^m
    #          = V(0)*w - V(1) - sum_{m=1}^{Mh} d2V(m) w^m + [tiny tail at Mh]!
    #   Take 2*Re[...]: 2*(cos th - Re V(1)) - 2*Re sum_{m=1}^{Mh} d2V(m) w^m!
    #   EXACT SAME BOUNDARY!
    #   Step 2: Apply (2 - w - w^{-1}) AGAIN to d2U (with d2U(0) = 2*Re U(1) - 2) and d2V (with d2V(0) = 2*Re V(1) - 2):
    #     s2^2 * (B_r - A_r) = 2 * Re[ (V(2) - U(2)) - 4*(V(1) - U(1)) ]
    #                          + 2 * Re sum_{m=1}^{Mh-1} d4V(m) w^m - Re sum_{m=1}^{c-1} d4U(m) w^m!
    Mh = min(M // 2, int(c + sqrt(45.0 / max(r * kap / c, 1e-10))) + 20)
    U_ext = np.concatenate([[np.conj(U[1])], U, [U[1]]]) # indices -1..c+1 (len c+3)
    d2U = U_ext[2:] - 2*U_ext[1:-1] + U_ext[:-2] # indices 0..c (len c+1), d2U[0] = d2U[c] = 2*Re(U[1]) - 2
    d2U_ext = np.concatenate([[np.conj(d2U[1])], d2U, [d2U[1]]]) # indices -1..c+1
    d4U = d2U_ext[2:] - 2*d2U_ext[1:-1] + d2U_ext[:-2] # indices 0..c -> we use 1..c-1 (d4U[1:-1])!

    V_s = V[:Mh+2]
    V_ext = np.concatenate([[np.conj(V_s[1])], V_s]) # indices -1..Mh+1
    d2V = V_ext[2:] - 2*V_ext[1:-1] + V_ext[:-2] # indices 0..Mh
    d2V_ext = np.concatenate([[np.conj(d2V[1])], d2V]) # indices -1..Mh
    d4V = d2V_ext[2:] - 2*d2V_ext[1:-1] + d2V_ext[:-2] # indices 0..Mh-1 -> we use 1..Mh-1 (d4V[1:])!

    bdy = 2.0 * float(((V[2] - U[2]) - 4.0 * (V[1] - U[1])).real)
    ws_e = np.exp(-1j * ms_e[1:c] * th)
    ws_i = np.exp(-1j * np.arange(1, Mh, dtype=np.float64) * th)
    rhs = (bdy + 2.0 * float(np.sum(d4V[1:] * ws_i).real) - float(np.sum(d4U[1:-1] * ws_e).real)) / (s2 * s2)
    err = abs(exact - rhs)

    # Paired triangle upper bound:
    mh = (c - 1) // 2
    # Pair m = 1..mh: Re( (2*d4V[m] - d4U[m] - conj(d4U[m])*...) )
    # Note d4U[c-m] = conj(d4U[m]) and w^{c-m} = w^{-m} = conj(w^m), so:
    #   Re(d4U[m]*w^m + d4U[c-m]*w^{c-m}) = 2 * Re(d4U[m] * w^m)!!
    # So 2*Re(d4V[m]*w^m) - 2*Re(d4U[m]*w^m) = 2*Re((d4V[m] - d4U[m]) * w^m)!!
    pair_d4 = 2.0 * float(np.sum(np.abs(d4V[1:mh+1] - d4U[1:mh+1])))
    mid_d4 = abs(float(d4U[mh+1].real)) if (c % 2 == 0) else 0.0
    tail_d4 = 2.0 * float(np.sum(np.abs(d4V[mh+1:])))
    num_p = (abs(bdy) + pair_d4 + mid_d4 + tail_d4) / (16.0 * r)
    ub4 = (abs(bdy) + pair_d4 + mid_d4 + tail_d4) / (s2 * s2 * r)
    sc = r * (kap**2) / c
    return kap, abs(exact) / r, ub4, err, num_p / sc

print("=== 1. EXACT FINITE-n 4TH-ORDER PAIRED ABEL VERIFICATION ===")
for n, c, j, r in [
    (1000, 10, 11, 1), (10000, 10, 11, 1), (10000, 30, 31, 1), (10000, 30, 32, 1),
    (10000, 50, 52, 1), (10000, 50, 81, 1), (10000, 100, 103, 1), (10000, 100, 109, 11),
    (100000, 100, 103, 1), (100000, 200, 323, 1), (100000, 200, 201, 1),
    (500, 10, 11, 1), (500, 20, 23, 1), (200, 6, 7, 1), (200, 8, 9, 1), (132, 6, 7, 1),
]:
    kap = 2*pi*sin(j*pi/n)*sin((j+c)*pi/n)/sin(c*pi/n)
    if kap >= 2.5 or c*kap*kap >= 42: continue
    res = check_finiten_4th_abel(n, c, j, r)
    if res is None: continue
    kap, ex, ub4, err, K4 = res
    print(f"  n={n:6d} c={c:3d} j={j:3d} r={r:2d} (kap={kap:.4f}, ck2={c*kap*kap:4.2f}): ex={ex:.5f} | ub4={ub4:.5f} | err={err:.1e} | K4={K4:.3f}")

# =============================================================================
# 2. CONTRAPOSITIVE MAP OF 3-WAY MINOR-ARC BOUND B_minor_3way >= 1.45 IN G*!
# =============================================================================
print("\n=== 2. CONTRAPOSITIVE MAP OF 3-WAY MINOR-ARC BOUND IN G* (c*q*kap > 0.85, c*kap^2 < 42) ===")
def B_minor_3way(c, jm, kap, C_p1=4.8, C_p2=0.60, K2=1.35, K4=4.5):
    nu = jm / c
    tot = 0.0
    for r in range(1, 81):
        tau_r = r * c * kap
        # Bound 0 (tight pair triangle):
        bp = (c / r) * min(C_p1 / (tau_r * tau_r), C_p2 / sqrt(tau_r))
        s2 = sin(pi * r * nu) ** 2
        if s2 > 1e-12:
            b2 = (K2 / s2) * min(1.0 / (r * c), sqrt(kap / (r * c)))
            b4 = (K4 / (s2 * s2)) * (r * kap * kap / c)
            tot += min(bp, b2, b4)
        else:
            tot += bp
    tot += C_p1 / (2.0 * 6400.0 * c * kap * kap)
    return tot

kaps = np.geomspace(0.002, 2.49, 60)
for c in [2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50, 100]:
    surv_D0 = []
    surv_Dne0 = []
    for jm in range(0, c // 2 + 1):
        g = gcd(c, jm) if jm > 0 else c
        q = c // g
        for kap in kaps:
            if c * q * kap <= 0.85 or c * kap * kap >= 42.0: continue
            bm = B_minor_3way(c, jm, kap)
            if bm >= 1.45:
                # Find best rational p0/r0 with r0 <= sqrt(c)+2
                best_r0 = q; best_D = 0; best_sc = 1e9
                for r0 in range(1, min(q, int(sqrt(c)) + 2) + 1):
                    p0 = round(r0 * jm / c)
                    D = abs(r0 * jm - c * p0)
                    if D == 0:
                        best_r0 = r0; best_D = 0; break
                    sc = D * sqrt(r0) / sqrt(c)
                    if sc < best_sc:
                        best_sc = sc; best_r0 = r0; best_D = D
                if best_D == 0:
                    surv_D0.append((jm, q, g, kap, bm))
                else:
                    surv_Dne0.append((jm, q, best_r0, best_D, c / best_r0, kap, c*best_r0*kap, bm))
    d0_qs = sorted(set((x[1], x[2]) for x in surv_D0))
    dne0_min = min((x[4] / (x[3]**2) for x in surv_Dne0), default=999.0)
    dne0_s = sorted(set((x[2], x[3], round(x[4], 1)) for x in surv_Dne0))[:6]
    print(f"  c={c:3d}: D=0 surv={len(surv_D0):3d} (q,g)={d0_qs[:6]} | D!=0 surv={len(surv_Dne0):3d} min(ct/D^2)={dne0_min:6.1f} (r0,D,ct)={dne0_s}")
