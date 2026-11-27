import mpmath as mp
import math

mp.mp.dps = 120

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

# Take K = 100 (or 150), b = [2]*K (commuting D_2!).
# Notice for DOUBLE-operator too, D_2(D_2 + y) commutes with itself!
#   prod_{j=1}^K D_2(D_2 + y)(1) = D_2^K (D_2 + y)^K (1) = sum_{i=0}^K C(K, i) D_2^{2K - i}(1) y^i !!
# Let's test BOTH single and double at K = 100, hitting step K+1 (to get negative Turan at t0)
# and THEN hitting step K+2 with b_{K+2} chosen to target u_i(t0) = D_{b} h_i(t0) + h_{i-1}(t0) = 0!!

for mode in ["SINGLE", "DOUBLE"]:
    K = 100
    i0 = 50
    max_pow = (2*K - (i0 - 3) + 1) if mode == "DOUBLE" else (K - (i0 - 3) + 1)
    D2 = [[mp.mpf(1)]]
    for m in range(1, max_pow + 1):
        D2.append(apply_D(mp.mpf(2), D2[-1]))
    # Need indices i0-3 .. i0+2
    P = {}
    for r in range(i0 - 3, i0 + 3):
        pw = (2*K - r) if mode == "DOUBLE" else (K - r)
        P[r] = [mp.mpf(math.comb(K, r)) * c for c in D2[pw]]

    # Find exact slope scale at i0
    s_scale = eval_mp(deriv_mp(P[i0]), mp.mpf(0)) / eval_mp(P[i0], mp.mpf(0))
    print(f"\n=== {mode} K={K}, i0={i0}, slope_scale={float(s_scale):.3e} ===")

    # Step K+1: hit with b1 = 4 * s_scale
    b1 = mp.mpf(4) * s_scale
    t0 = -mp.mpf(3) / b1  # where Turan was -0.018!
    P1 = {}
    for r in range(i0 - 2, i0 + 3):
        u = add_mp(apply_D(b1, P[r]), P[r-1])
        P1[r] = apply_D(b1, u) if mode == "DOUBLE" else u

    # Check Turan of P1 at t0 and scan around t0!
    best_tur = (mp.mpf('1e100'), None)
    for mult in [0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0]:
        tt = t0 * mp.mpf(str(mult))
        v0, v1, v2 = eval_mp(P1[i0-1], tt), eval_mp(P1[i0], tt), eval_mp(P1[i0+1], tt)
        rel = (v1*v1 - v0*v2) / (v1*v1 + abs(v0*v2))
        if rel < best_tur[0]: best_tur = (rel, tt)
    t_neg = best_tur[1]
    print(f"  Step K+1 (b1={float(b1):.3e}): worst rel Turan = {float(best_tur[0]):+.6f} at t_neg = {float(t_neg):.4e}")
    # Print signs of P1[i0-1], P1[i0], P1[i0+1] at t_neg!
    vals = [eval_mp(P1[r], t_neg) for r in [i0-1, i0, i0+1]]
    print(f"    P1(t_neg) = {[float(v) for v in vals]}")

    # NOW STEP K+2: Can we choose b2 >= 2 so that W(P2[i0+1], P2[i0]) <= 0 (BROKEN INTERLACING!!) near t_neg?!
    # Recall u_r(t) = (1 + b2 t) P1[r](t) + t(1+t) P1'[r](t) + P1[r-1](t).
    # For u_{i0}(t_neg) = 0, solve for b2:
    #   b2 = -(P1[i0](t_neg) + t_neg(1+t_neg) P1'[i0](t_neg) + P1[i0-1](t_neg)) / (t_neg * P1[i0](t_neg))
    d0 = eval_mp(deriv_mp(P1[i0]), t_neg)
    b2_star = -(vals[1] + t_neg*(1+t_neg)*d0 + vals[0]) / (t_neg * vals[1])
    print(f"    Exact b2* making u_{{i0}}(t_neg)=0: b2* = {float(b2_star):.4e}")

    # Test a wide sweep of b2 (around b2* and powers of 10 times s_scale) for W(P2[i0+1], P2[i0]) and J_raw!
    min_W_rel = (mp.mpf('1e100'), None)
    min_J2 = (mp.mpf('1e100'), None)
    b2_cands = [mp.mpf(2), mp.mpf(10), s_scale * mp.mpf('0.1'), s_scale, s_scale * mp.mpf(4), s_scale * mp.mpf(16), s_scale * mp.mpf(100)]
    if b2_star > 2:
        b2_cands.extend([b2_star * mp.mpf(str(m)) for m in [0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 2.0]])
    for b2 in b2_cands:
        P2 = {}
        for r in [i0 - 1, i0, i0 + 1]:
            u = add_mp(apply_D(b2, P1[r]), P1[r-1])
            P2[r] = apply_D(b2, u) if mode == "DOUBLE" else u
        for mult in [0.2, 0.5, 0.8, 0.95, 1.0, 1.05, 1.2, 1.5, 2.0, 5.0]:
            tt = t_neg * mp.mpf(str(mult))
            w1 = eval_W(P2[i0], P2[i0-1], tt)
            w2 = eval_W(P2[i0+1], P2[i0], tt)
            w12 = eval_W(P2[i0+1], P2[i0-1], tt)
            if w1 <= 0 or w2 <= 0:
                print(f"  *** BOOM! INTERLACING BROKEN IN {mode} (b_j >= 2, a_j=1, w_j=1) AT STEP K+2! b2={float(b2):.3e}, t={float(tt):.3e}: w1={float(w1):.3e}, w2={float(w2):.3e}")
                min_W_rel = (mp.mpf(-1), (float(b2), float(tt)))
                break
            J = 4*w1*w2/(w12*w12) if w12 != 0 else mp.mpf('1e10')
            if J < min_J2[0]: min_J2 = (J, (float(b2), float(tt)))
        if min_W_rel[0] < 0: break
    print(f"  Step K+2 summary: min_W_status={float(min_W_rel[0])}, min_Jraw={float(min_J2[0]):+.6f} at {min_J2[1]}")
