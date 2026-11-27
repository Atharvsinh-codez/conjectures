"""
TEST THE EXACT 1-STEP AND 2-STEP SYMMETRIC ABEL IDENTITY:
  (w + w^{-1} - 2) * sum_{m=1}^{M-1} g(m) w^m
    = sum_{m=1}^{M-1} Delta^2 g(m) w^m - [g(0) w - g(1) + g(M) w^{M-1} - g(M-1) w^M]
where w = exp(2*pi*i*r*nu), so w + w^{-1} - 2 = -4*sin^2(pi*r*nu) < 0!

Apply to U_r(m) = a_{-m}^{-r} * w^{-m} (m=0..c) and Vb_r(m) = conj(a_m^r) * w^{-m} (m=0..n-c):
  4*sin^2(pi*r*nu) * (B_r - A_r)
    = 2*Re(U_r(1) - Vb_r(1)) + Re sum_{m=1}^{c-1} Delta^2 U_r(m) w^m - Re sum_{m=1}^{n-c-1} Delta^2 Vb_r(m) w^m!

Wait: what if we apply it TWICE?!
  16*sin^4(pi*r*nu) * (B_r - A_r) = ...
Let's verify BOTH exact identities to 1e-13 AND measure the exact triangle-inequality upper bound
on |B_r - A_r| / r from 1-step vs 2-step vs paired-1-step across our hard Regime-2 points!
"""
import numpy as np
from math import pi, sin, cos, gcd

def check_abel_exact(n, c, j, r):
    rj = sin(j * pi / n) / sin((j + c) * pi / n)
    sj = sin(c * pi / n) / sin((j + c) * pi / n)
    th_star = pi * (n - c) / n
    nu = sin(j * pi / n) * cos((j + c) * pi / n) / sin(c * pi / n)
    w = np.exp(-2j * pi * r * nu)
    s2 = 4.0 * (sin(pi * r * nu) ** 2)

    # Exact U_r(m) for m = 0..c: a_{-m}^{-r} * w^{-m}
    # theta_{-m} = th_star + 2*m*pi/n => log(a_{-m}^{-1}) = -kap*m - 2*pi*i*nu*m + ...
    ms_e = np.arange(0, c + 1, dtype=np.float64)
    th_e = th_star + 2.0 * ms_e * pi / n
    log_ae_inv = n * np.log((1.0 + rj * np.exp(1j * th_e)) / (sj * np.exp(1j * j * pi / n)))
    U = np.exp(r * log_ae_inv + 2j * pi * r * nu * ms_e)

    # Exact V_r(m) for m = 0..n-c: a_m^r * w^{-m} (NO conj! a_m ALREADY has phase -2*pi*nu*m!)
    ms_i = np.arange(0, n - c + 1, dtype=np.float64)
    th_i = th_star - 2.0 * ms_i * pi / n
    log_ai = n * np.log((sj * np.exp(1j * j * pi / n)) / (1.0 + rj * np.exp(1j * th_i)))
    Vb = np.exp(r * log_ai + 2j * pi * r * nu * ms_i)

    A_r = float(np.sum((U[1:c] * (w ** ms_e[1:c])).real))
    B_r = float(np.sum((Vb[1:n-c] * (w ** ms_i[1:n-c])).real))
    diff_exact = B_r - A_r

    # 1-step symmetric Abel identity check:
    d2U = U[2:] - 2.0 * U[1:-1] + U[:-2]
    d2Vb = Vb[2:] - 2.0 * Vb[1:-1] + Vb[:-2]
    bdy1 = 2.0 * float((U[1] - Vb[1]).real) # Wait: check exact sign of bdy!
    rhs1 = (
        float(((U[0]*w - U[1] + U[c]*(w**(c-1)) - U[c-1]*(w**c)) -
               (Vb[0]*w - Vb[1] + Vb[n-c]*(w**(n-c-1)) - Vb[n-c-1]*(w**(n-c)))).real)
        - float(np.sum((d2U * (w ** ms_e[1:c])).real))
        + float(np.sum((d2Vb * (w ** ms_i[1:n-c])).real))
    )
    # Wait: let's check identity sign:
    # sum d2g w^m = (w + w^{-1} - 2) sum g w^m + (g0 w - g1 + gM w^{M-1} - g_{M-1} w^M)
    # => 4 sin^2 * sum g w^m = (g0 w - g1 + gM w^{M-1} - g_{M-1} w^M) - sum d2g w^m!
    # So 4 sin^2 * (B_r - A_r) = [bdy(Vb) - bdy(U)] - [sum d2Vb w^m - sum d2U w^m]!
    lhs1 = s2 * diff_exact
    # Notice bdy(Vb) - bdy(U) = 2*Re(U[1] - Vb[1])! Check!
    bdy_Vb = (Vb[0]*w - Vb[1] + Vb[n-c]*(w**(n-c-1)) - Vb[n-c-1]*(w**(n-c))).real
    bdy_U = (U[0]*w - U[1] + U[c]*(w**(c-1)) - U[c-1]*(w**c)).real
    rhs1_true = (bdy_Vb - bdy_U) - float(np.sum((d2Vb * (w ** ms_i[1:n-c])).real)) + float(np.sum((d2U * (w ** ms_e[1:c])).real))

    # Triangle bound from 1-step (pairing d2U and d2Vb for m = 1..c//2!):
    mh = (c - 1) // 2
    pair_d2 = 2.0 * float(np.sum(np.abs(d2U[:mh] - d2Vb[:mh])))
    mid_d2 = float(abs(d2U[mh])) if (c % 2 == 0) else 0.0
    tail_d2 = 2.0 * float(np.sum(np.abs(d2Vb[mh:(n-c)//2])))
    ub1_paired = (abs(bdy_Vb - bdy_U) + pair_d2 + mid_d2 + tail_d2) / s2
    ub1_unpaired = (abs(bdy_Vb - bdy_U) + float(np.sum(np.abs(d2U))) + float(np.sum(np.abs(d2Vb)))) / s2

    # 2-step symmetric Abel identity (apply again to d2U and d2Vb with endpoints d2g(0), d2g(M)??
    # Just apply to sum_{m=1}^{M-1} h(m) w^m where h(m) = d2g(m) for m=1..M-1, extending h(0)=0, h(M)=0!
    # Wait: if h(0)=0, h(M)=0, then bdy = -h(1) - h(M-1)*w^M = -2*Re(h(1))!
    # So sum_{m=1}^{M-1} h(m) w^m = [bdy_h - sum_{m=1}^{M-1} d2h(m) w^m] / s2!
    hU = np.concatenate([[0.0], d2U, [0.0]])
    hV = np.concatenate([[0.0], d2Vb, [0.0]])
    d4U = hU[2:] - 2.0 * hU[1:-1] + hU[:-2]
    d4V = hV[2:] - 2.0 * hV[1:-1] + hV[:-2]
    bdy2_U = (-hU[1] - hU[c-1]*(w**c)).real
    bdy2_V = (-hV[1] - hV[n-c-1]*(w**(n-c))).real
    ub2_unpaired = (abs(bdy_Vb - bdy_U) / s2) + (abs(bdy2_V - bdy2_U) + float(np.sum(np.abs(d4U))) + float(np.sum(np.abs(d4V)))) / (s2 * s2)
    return diff_exact, abs(lhs1 - rhs1_true), abs((bdy_Vb - bdy_U) - bdy1), ub1_paired, ub1_unpaired, ub2_unpaired

print("=== EXACT 1-STEP & 2-STEP SYMMETRIC ABEL SUMMATION CHECK ===")
for n, c, j, r in [
    (1000, 10, 11, 1),
    (10000, 10, 11, 1),
    (10000, 50, 81, 1),
    (10000, 50, 81, 2),
    (10000, 200, 323, 1),
    (10000, 200, 323, 2),
    (10000, 200, 241, 1),
    (100000, 200, 323, 1),
]:
    kap = 2*pi*sin(j*pi/n)*sin((j+c)*pi/n)/sin(c*pi/n)
    de, err1, errb, ub1p, ub1u, ub2u = check_abel_exact(n, c, j, r)
    print(f"n={n:6d} c={c:3d} j={j:3d} r={r} (kap={kap:.4f}): exact={de:+.4e} | id_err={err1:.1e} bdy_err={errb:.1e} | ub1_pair={ub1p:.3e} ub1_unp={ub1u:.3e} ub2_unp={ub2u:.3e}")
