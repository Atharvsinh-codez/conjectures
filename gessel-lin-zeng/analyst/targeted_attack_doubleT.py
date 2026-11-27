import numpy as np
import mpmath as mp

# TARGETED ADVERSARIAL ATTACK ON DOUBLE-OPERATOR prod_{j=1}^k T_j(T_j + w_j y)(1) AND prod D_{m_j}(D_{m_j}+y)(1)!
# Key question: Is double-operator ALWAYS real-rooted for ALL 0 < a_j < b_j, w_j > 0,
# OR can a targeted attack (like @algebraist's k=9 attack on single-operator) break it?!
#
# Recall HOW @algebraist broke single-operator at k=9:
#   Step 1: Find parameters for j=1..k0 such that raw J_i(t0) < 0 at some t0 in (-1, 0).
#           (WAIT: Even for GLZ itself, raw J_i(t0) < 0 at k=35! Can we get raw J_i(t0) < 0 at k=4,5,6 for double_T?!)
#   Step 2: Once raw J_i(t0) < 0, by Exact T-Turan equivalence, there EXISTS Phi* = a + b*t0 such that
#           (T p_i)^2 - (T p_{i-1})(T p_{i+1}) < 0 at t0!
#           Wait: in double_T, q_i = T(T p_i + w p_{i-1})!
#           If we pick w VERY SMALL (w -> 0+), then q_i ~ T^2 p_i!
#           WAIT A SECOND!!
#           CHECK WHAT HAPPENS AS w_{k+1} -> 0+ IN DOUBLE_T:
#             q_i = T^2 p_i + O(w)!
#           Can (T^2 p_i)^2 - (T^2 p_{i-1})(T^2 p_{i+1}) < 0 at t0?!
#           YES, IF raw J_i(t0) < 0 for (p_{i-1}, p_i, p_{i+1}) OR for (T p_{i-1}, T p_i, T p_{i+1})!!
#           AND THEN, once we have Turan q_i(t1)^2 < q_{i-1}(t1) q_{i+1}(t1) at step k1:
#           At step k1 + 1, pick (a, b, w) so u_i = T q_i + w q_{i-1} and u_{i+1} = T q_{i+1} + w q_i FAIL TO INTERLACE!
#           Since T preserves non-interlacing? Wait: T CAN'T fix non-interlacing if W(u_{i+1}, u_i) changes sign?
#           Wait: actually T CAN sometimes heal mild non-interlacing, OR if u_{i+1} has roots colliding past u_i,
#           then T(u_{i+1}) and T(u_i) STILL don't interlace, AND at step k1 + 2, T(T r_{i+1} + w r_i) gets COMPLEX ROOTS!!
#
# LET'S RUN THIS EXACT 3-STAGE TARGETED ATTACK IN 100-DPS MPMATH!!

mp.mp.dps = 100

def apply_T(p, a, b):
    d = len(p) - 1
    res = [mp.mpf(0)] * (d + 2)
    for n in range(d + 1):
        res[n] += (a + n) * p[n]
        res[n + 1] += (b + n) * p[n]
    return res

def add_p(p, q):
    n = max(len(p), len(q))
    r = [mp.mpf(0)] * n
    for i, v in enumerate(p): r[i] += v
    for i, v in enumerate(q): r[i] += v
    return r

def step_double(P, a, b, w):
    k = len(P)
    newP = []
    for i in range(k + 1):
        t2 = apply_T(P[i], a, b) if i < k else [mp.mpf(0)]
        wy = [w * c for c in P[i - 1]] if i >= 1 else [mp.mpf(0)]
        newP.append(apply_T(add_p(t2, wy), a, b))
    return newP

def eval_p(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v

def deriv_p(p):
    return [n * p[n] for n in range(1, len(p))] if len(p) > 1 else [mp.mpf(0)]

def eval_W(f, g, t):
    return eval_p(f, t) * eval_p(deriv_p(g), t) - eval_p(deriv_p(f), t) * eval_p(g, t)

# Stage 1: Check raw J_i(t) = 4 W(p_{i-1}, p_i) W(p_i, p_{i+1}) / W(p_{i-1}, p_{i+1})^2 on (-1, 0)
# for GLZ double-operator AND for step-jump double-operator (at t near 0^-, where GLZ has J_i < 1 at k=35)!
print("=== STAGE 1: FIND RAW J_i(t) < 1 IN DOUBLE-OPERATOR AT SMALL k! ===")
# For GLZ, at t = 0: check 4 W01 W12 / W02^2 as a function of k!
# Wait: at t = 0, p_i(0) = C(k, i), p_i'(0) = ...
# For general double_T with w_j = 1, a_j = 1, b_j = m_j:
# Try extreme m_j!
configs = [
    ("GLZ k=5", [(1, j+1, 1) for j in range(1, 6)]),
    ("m=[2]*2+[100]*3", [(1, 2, 1)]*2 + [(1, 100, 1)]*3),
    ("m=[100]*2+[2]*3", [(1, 100, 1)]*2 + [(1, 2, 1)]*3),
    ("m=[1.01]*3+[1e4]*3", [(1, mp.mpf('1.01'), 1)]*3 + [(1, mp.mpf('1e4'), 1)]*3),
    ("m=[1e4]*3+[1.01]*3", [(1, mp.mpf('1e4'), 1)]*3 + [(1, mp.mpf('1.01'), 1)]*3),
    ("a=1,b=100,w=[1]*2+[1e-4]*2", [(1, 100, 1)]*2 + [(1, 100, mp.mpf('1e-4'))]*2),
    ("a=1,b=2,w=[1,10,100,1000,1e4]", [(1, 2, mp.power(10, j)) for j in range(5)]),
    ("a=[1e-3]*2+[1]*3, b=[10]*5, w=1", [(mp.mpf('1e-3'), 10, 1)]*2 + [(1, 10, 1)]*3),
]

best_J_cfg = None
for name, steps in configs:
    P = [[mp.mpf(1)]]
    for a, b, w in steps:
        P = step_double(P, mp.mpf(a), mp.mpf(b), mp.mpf(w))
    min_J = mp.mpf('1e100')
    min_info = None
    for i in range(1, len(P) - 1):
        for tv in [mp.mpf('-1e-12'), mp.mpf('-1e-6'), mp.mpf('-1e-3'), mp.mpf('-0.05'), mp.mpf('-0.5'), mp.mpf('-0.95')]:
            w1 = eval_W(P[i-1], P[i], tv)
            w2 = eval_W(P[i], P[i+1], tv)
            w12 = eval_W(P[i-1], P[i+1], tv)
            if w12 != 0:
                J = 4 * w1 * w2 / (w12 * w12)
                if J < min_J:
                    min_J = J
                    min_info = (i, float(tv), P, steps)
    print(f"  {name:35s}: min_J = {float(min_J):+.6f} at i={min_info[0]}, t={min_info[1]:.1e}")
    if float(min_J) < 1.0 and (best_J_cfg is None or float(min_J) < best_J_cfg[0]):
        best_J_cfg = (float(min_J), name, min_info)

if best_J_cfg is not None:
    Jval, name, (i0, t0_f, P0, steps0) = best_J_cfg
    print(f"\n=== STAGE 2: FOUND RAW J < 1 ({Jval:.6f}) IN '{name}' AT i={i0}, t0={t0_f:.1e}! ===")
    t0 = mp.mpf(str(t0_f))
    # Compute exact optimal Phi* = a + b*t0 that MINIMIZES (T p_i)^2 - (T p_{i-1})(T p_{i+1}) at t0!
    # Recall T p = Phi * p + t(1+t) p'!
    # Let u_r = p_r(t0), v_r = t0*(1+t0)*p_r'(t0) for r in {i0-1, i0, i0+1}.
    # (Phi u1 + v1)^2 - (Phi u0 + v0)(Phi u2 + v2) = A Phi^2 + B Phi + C
    u0, u1, u2 = eval_p(P0[i0-1], t0), eval_p(P0[i0], t0), eval_p(P0[i0+1], t0)
    v0 = t0*(1+t0)*eval_p(deriv_p(P0[i0-1]), t0)
    v1 = t0*(1+t0)*eval_p(deriv_p(P0[i0]), t0)
    v2 = t0*(1+t0)*eval_p(deriv_p(P0[i0+1]), t0)
    A_q = u1*u1 - u0*u2
    B_q = 2*u1*v1 - u0*v2 - u2*v0
    C_q = v1*v1 - v0*v2
    Phi_star = -B_q / (2 * A_q)
    min_val = C_q - B_q*B_q / (4 * A_q)
    print(f"  At t0={t0_f:.1e}: A={float(A_q):.4e}, Phi*={float(Phi_star):.6f}, min_quad={float(min_val):.4e} (<0!)")
    # Choose 0 < a < b so that a + b*t0 = Phi*!
    # Since t0 < 0: b = (a - Phi*) / (-t0). Need b > a <=> a - Phi* > -t0 * a <=> (1 + t0)*a > Phi*!
    # Pick a = max(1e-3, Phi*/(1+t0) + 1) if Phi* > 0 else 1e-3!
    if Phi_star > 0:
        a_hit = Phi_star / (1 + t0) + mp.mpf('1')
    else:
        a_hit = mp.mpf('1e-3')
    b_hit = (a_hit - Phi_star) / (-t0)
    print(f"  Hitting Phi* with a={float(a_hit):.4e}, b={float(b_hit):.4e} (check a+b*t0={float(a_hit+b_hit*t0):.6f})")
    # Check BOTH single T and double T(T + w y) with tiny w = 1e-10!
    P_single = [apply_T(p, a_hit, b_hit) for p in P0]
    ts1 = eval_p(P_single[i0], t0)**2 - eval_p(P_single[i0-1], t0)*eval_p(P_single[i0+1], t0)
    print(f"  After 1x T: Turan at t0 = {float(ts1):.4e}")
    P_double = step_double(P0, a_hit, b_hit, mp.mpf('1e-15'))
    td1 = eval_p(P_double[i0], t0)**2 - eval_p(P_double[i0-1], t0)*eval_p(P_double[i0+1], t0)
    print(f"  After T(T + 1e-15 y): Turan at t0 = {float(td1):.4e}")
    # Also scan t near t0 for P_double!
    best_t_d = (mp.mpf('1e100'), None)
    for mult in [0.1, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 5.0, 10.0]:
        tt = t0 * mp.mpf(str(mult))
        if tt <= -1: continue
        val = eval_p(P_double[i0], tt)**2 - eval_p(P_double[i0-1], tt)*eval_p(P_double[i0+1], tt)
        denom = eval_p(P_double[i0], tt)**2 + abs(eval_p(P_double[i0-1], tt)*eval_p(P_double[i0+1], tt))
        rel = val / denom
        if rel < best_t_d[0]: best_t_d = (rel, float(tt))
    print(f"  Best rel Turan after T(T+1e-15 y) near t0: {float(best_t_d[0]):+.6f} at t={best_t_d[1]:.3e}")

    # WAIT: What if we hit the EXACT optimal Phi* for T^2 directly?!
    # Check if (T^2 p_i)^2 - (T^2 p_{i-1})(T^2 p_{i+1}) CAN BE NEGATIVE at t0 for ANY (a, b) with 0 < a < b!
    print("\n=== CHECKING IF (T^2 p_i)^2 - (T^2 p_{i-1})(T^2 p_{i+1}) CAN BE NEGATIVE AT t0! ===")
    best_t2 = (mp.mpf('1e100'), None)
    for kpow in np.linspace(-2, 14, 65):
        b_try = mp.power(10, kpow)
        for a_try in [mp.mpf('1e-4'), mp.mpf('1e-2'), mp.mpf('1'), mp.mpf('10'), b_try * mp.mpf('0.5')]:
            if a_try >= b_try: continue
            q0 = apply_T(apply_T(P0[i0-1], a_try, b_try), a_try, b_try)
            q1 = apply_T(apply_T(P0[i0], a_try, b_try), a_try, b_try)
            q2 = apply_T(apply_T(P0[i0+1], a_try, b_try), a_try, b_try)
            v_0, v_1, v_2 = eval_p(q0, t0), eval_p(q1, t0), eval_p(q2, t0)
            rel = (v_1*v_1 - v_0*v_2) / (v_1*v_1 + abs(v_0*v_2))
            if rel < best_t2[0]:
                best_t2 = (rel, (float(a_try), float(b_try)))
    print(f"  Min rel Turan of T^2(P0) at t0={t0_f:.1e}: {float(best_t2[0]):+.6f} at (a,b)={best_t2[1]}")
    # If negative, try to break interlacing and RZ on next 2 steps!!
    if float(best_t2[0]) < 0:
        a_n, b_n = mp.mpf(str(best_t2[1][0])), mp.mpf(str(best_t2[1][1]))
        P_neg = step_double(P0, a_n, b_n, mp.mpf('1e-15'))
        print(f"  *** NEGATIVE TURAN ACHIEVED IN DOUBLE-OPERATOR AT k={len(P_neg)-1}! Now trying to break interlacing & RZ! ***")
        # At t0, P_neg[i0]^2 < P_neg[i0-1]*P_neg[i0+1]!
        # Choose (a, b, w) at next step so that u_{i0}(t0) = T P_neg[i0](t0) + w P_neg[i0-1](t0) = 0!
        # Wait: need w = -T P_neg[i0](t0) / P_neg[i0-1](t0) > 0, and T(P_neg[i0+1]) + w P_neg[i0] to have WRONG sign!
        # Test a grid of (a, b, w) targeting u_{i0}(t0) = 0!
        found_break = False
        for a_k in [mp.mpf('1e-4'), mp.mpf('1')]:
            for b_pow in [-1, 1, 3, 5, 7, 9, 11, 13]:
                b_k = a_k + mp.power(10, b_pow)
                tp_i = eval_p(apply_T(P_neg[i0], a_k, b_k), t0)
                p_im1 = eval_p(P_neg[i0-1], t0)
                if p_im1 == 0: continue
                w_k = -tp_i / p_im1
                if w_k <= 0: continue
                P_next = step_double(P_neg, a_k, b_k, w_k)
                # Check interlacing of P_next[i0+1] and P_next[i0]!
                r1 = mp.polyroots(P_next[i0][::-1], maxsteps=200, extraprec=50)
                r2 = mp.polyroots(P_next[i0+1][::-1], maxsteps=200, extraprec=50)
                im1 = max(abs(mp.im(z)) for z in r1); im2 = max(abs(mp.im(z)) for z in r2)
                if im1 > 1e-15 or im2 > 1e-15:
                    print(f"  *** COMPLEX ROOT ALREADY AT STEP k+1! im1={float(im1):.3e}, im2={float(im2):.3e}!")
                    found_break = True; break
                r1s = sorted(mp.re(z) for z in r1); r2s = sorted(mp.re(z) for z in r2)
                ok_int = all(r1s[m] < r2s[m] < r1s[m+1] for m in range(len(r2s)))
                if not ok_int:
                    print(f"  *** INTERLACING BROKEN IN DOUBLE-OPERATOR AT k={len(P_next)-1}! (a={float(a_k)}, b={float(b_k):.1e}, w={float(w_k):.3e}) ***")
                    # Now one more step to break RZ!
                    for w_last in [mp.mpf('1e-6'), mp.mpf('1e-2'), mp.mpf('1'), mp.mpf('1e2'), mp.mpf('1e6'), mp.mpf('1e10')]:
                        P_fin = step_double(P_next, mp.mpf('1'), mp.mpf('10'), w_last)
                        rf = mp.polyroots(P_fin[i0+1][::-1], maxsteps=200, extraprec=50)
                        imf = float(max(abs(mp.im(z)) for z in rf))
                        if imf > 1e-10:
                            print(f"  *** BOOM! REAL-ROOTEDNESS BROKEN IN DOUBLE-OPERATOR AT k={len(P_fin)-1}, i={i0+1}! max_im={imf:.4e}, w_last={float(w_last):.1e} ***")
                            found_break = True; break
                    if found_break: break
            if found_break: break
        if not found_break:
            print("  Interlacing/RZ survived all targeted attacks!")
