import mpmath as mp
import math

# TEST THE EXACT PREDICTED BREAK OF UNIVERSAL prod_{j=1}^k (D_{b_j} + y)(1) WITH a_j=1, w_j=1, b_j >= 2!!
# Our exact discriminant formula for M_i^{(w)} at t = -eps gave:
#   Delta(w) = (2 + 2w + 3w^2)^2 - 4(1+w)(1 + 4w + 2w^2 + 2w^3)
#            = w^4 - 8w^3 - 10w^2 - 12w ... wait: check exact coeff!
#   At w = 6: Delta = +72 > 0!
#   In prod (D_2 + y)^K (1), at index i, the effective ratio is w_eff = C(K, i-1)/C(K, i) = i / (K - i + 1)!
#   Wait: we ALSO need J_raw(0) = (8/9) * (i+1)(K-i+1)/(i(K-i)) < 1, which requires K - i >= 9!
#   So pick K - i = 20 and i = 200 (K = 220, i = 200):
#     w_eff = 200 / 21 = 9.52 > 5.8 !!
#     J_raw(0) = (8/9) * (201/200) * (21/20) = 0.938 < 1 !!
#   Let's test BOTH SINGLE and DOUBLE at K = 220, i = 200 (takes 0.2 seconds using commuting D_2!)!

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

# First check exact Delta(w) for w = 1..15:
for w in [1, 5, 6, 10, 20]:
    A = 1 + w; B = 2 + 2*w + 3*w*w; C = 1 + 4*w + 2*w*w + 2*w**3
    print(f"  w={w:2d}: Delta = {B*B - 4*A*C:+d}, Phi_min = {-B/(2*A):.3f}, min_val = {C - B*B/(4*A):+.3f}")

# Test SINGLE at K = 220, i0 = 200 (so K - i0 = 20, max degree of D_2^{23}(1) is ONLY 23!! INSTANT!!):
K = 220; i0 = 200  # w_eff = 200/21 = 9.52!
D2 = [[mp.mpf(1)]]
for m in range(1, 30):
    D2.append(apply_D(mp.mpf(2), D2[-1]))

# Wait: also try K = 450, i0 = 420 (K-i0 = 30, w_eff = 13.5, J_raw(0) = 8/9 * 421/420 * 31/30 = 0.9207 < 1!)
for K, i0 in [(220, 200), (450, 420), (1000, 960)]:
    w_eff = i0 / (K - i0 + 1)
    max_p = K - (i0 - 3) + 1
    while len(D2) <= max_p:
        D2.append(apply_D(mp.mpf(2), D2[-1]))
    P = {r: [mp.mpf(math.comb(K, r)) * c for c in D2[K - r]] for r in range(i0 - 3, i0 + 3)}
    # Normalize by C(K, i0) so numbers stay O(1)
    norm = mp.mpf(math.comb(K, i0))
    P = {r: [c / norm for c in p] for r, p in P.items()}
    # s2 = slope of P[i0+1] / P[i0+1](0)
    s2 = eval_mp(deriv_mp(P[i0+1]), mp.mpf(0)) / eval_mp(P[i0+1], mp.mpf(0))
    # Predict Phi_min = -(2 + 2w + 3w^2) / (2(1+w))!
    Phi_opt = -(2 + 2*w_eff + 3*w_eff**2) / (2*(1 + w_eff))
    print(f"\n=== SINGLE K={K}, i0={i0} (K-i0={K-i0}, w_eff={w_eff:.2f}, Phi_opt={Phi_opt:.2f}, s2={float(s2):.3e}) ===")
    # Pick eps = alpha / s2, and b1 = (1 - Phi_opt) / eps!
    found_neg_M = None
    for alpha in [0.001, 0.005, 0.01, 0.05, 0.1, 0.2]:
        eps = mp.mpf(str(alpha)) / s2
        for phi_mult in [0.7, 0.85, 1.0, 1.15, 1.3]:
            Phi = mp.mpf(str(Phi_opt * phi_mult))
            b1 = (1 - Phi) / eps
            if b1 < 2: continue
            P1 = {r: add_mp(apply_D(b1, P[r]), P[r-1]) for r in range(i0 - 2, i0 + 3)}
            tv = -eps
            v0, v1, v2 = eval_mp(P1[i0-1], tv), eval_mp(P1[i0], tv), eval_mp(P1[i0+1], tv)
            w21 = eval_W(P1[i0+1], P1[i0], tv)
            M_val = (-tv)*(1+tv)*w21 + v1*v1 - v0*v2
            rel_M = M_val / (abs((-tv)*(1+tv)*w21) + v1*v1 + abs(v0*v2))
            if M_val < 0 and (found_neg_M is None or rel_M < found_neg_M[0]):
                found_neg_M = (rel_M, alpha, float(Phi), b1, tv, P1)
    if found_neg_M is None:
        print("  M_i stayed > 0 across all probes!")
    else:
        rel_M, alpha, Phi_f, b1, t_neg, P1 = found_neg_M
        print(f"  *** BOOM! M_i < 0 ACHIEVED WITH a=1, w=1, b>=2! rel_M={float(rel_M):+.6f} at alpha={alpha}, Phi={Phi_f:.2f}, b1={float(b1):.3e} ***")
        # Now hit step K+2 with exact b2* making P2[i0](t_neg) = 0 -> MUST BREAK INTERLACING!!
        v0, v1 = eval_mp(P1[i0-1], t_neg), eval_mp(P1[i0], t_neg)
        d1 = eval_mp(deriv_mp(P1[i0]), t_neg)
        b2_star = -(v1 + t_neg*(1+t_neg)*d1 + v0) / (t_neg * v1)
        print(f"  Exact targeted b2* = {float(b2_star):.4e} (>= 2? {b2_star >= 2})")
        if b2_star >= 2:
            P2 = {r: add_mp(apply_D(b2_star, P1[r]), P1[r-1]) for r in [i0 - 1, i0, i0 + 1]}
            # Check W(P2[i0+1], P2[i0]) at t_neg!
            # Recall P2[i0](t_neg) == 0, and P2[i0+1](t_neg) has WRONG SIGN relative to P1[i0](t_neg)!
            w_check = eval_W(P2[i0+1], P2[i0], t_neg)
            print(f"  P2[i0](t_neg) = {float(eval_mp(P2[i0], t_neg)):.3e}, P2[i0+1](t_neg)*P1[i0](t_neg) = {float(eval_mp(P2[i0+1], t_neg)*v1):.3e} (<0!)")
            print(f"  W(P2[i0+1], P2[i0])(t_neg) = {float(w_check):.4e}")
            # Also check roots of P2[i0] and P2[i0+1] near t_neg, AND check if W < 0 somewhere near t_neg!
            min_w = mp.mpf('1e100')
            for mult in [0.9, 0.95, 0.99, 0.999, 1.0, 1.001, 1.01, 1.05, 1.1]:
                ww = eval_W(P2[i0+1], P2[i0], t_neg * mp.mpf(str(mult)))
                if ww < min_w: min_w = ww
            print(f"  min W(P2[i0+1], P2[i0]) near t_neg = {float(min_w):.4e}")
            if min_w < 0:
                print(f"  *** INTERLACING OFFICIALLY BROKEN FOR a_j=1, w_j=1, b_j>=2 AT K+2 = {K+2}!! ***")
                # Check if 1 more step breaks REAL-ROOTEDNESS!
                for b3_mult in [0.1, 0.5, 1.0, 2.0, 10.0, 100.0]:
                    b3 = max(mp.mpf(2), b2_star * mp.mpf(str(b3_mult)))
                    P3_i = add_mp(apply_D(b3, P2[i0+1]), P2[i0])
                    rts = mp.polyroots(P3_i[::-1], maxsteps=200, extraprec=50)
                    mx_im = float(max(abs(mp.im(z)) for z in rts))
                    if mx_im > 1e-15:
                        print(f"  *** NUCLEAR BOOM! COMPLEX ROOTS FOR a_j=1, w_j=1, b_j>=2 AT K+3 = {K+3}, b3={float(b3):.3e}, max_im={mx_im:.4e}!! ***")
                        break
