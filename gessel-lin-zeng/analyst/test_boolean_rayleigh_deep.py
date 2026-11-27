import mpmath as mp
import numpy as np

# 1. Test Boolean-lattice Rayleigh h_{R u {a}} h_{R u {b}} >= h_R h_{R u {a,b}} on (-1, 0)
#    at k = 14, 20 (including extreme pairs a=1, b=14 AND random wild D_{m_j}!!)!
#    Wait: remember how at step b it was (D_b^2 p)(D_b q) - (D_b^2 q)(D_b p) = (-t)(1+t) W(D_b p, D_b q) > 0,
#    and THEN we apply V = D_{m_r} ... D_{m_1} to ALL FOUR!
#    WAIT A SECOND! Test directly: if p << q in (-1, 0), and Q0 = D_b p, Q1 = D_b q, R0 = D_b^2 p, R1 = D_b^2 q,
#    does V(R0)*V(Q1) - V(R1)*V(Q0) >= 0 on (-1, 0) for ANY V = D_{m_1} ... D_{m_r}?!
# 2. Test NON-ADJACENT exchange h_{1,2} h_{3,4} - h_1 h_{2,3,4} on (-1, 0) at k = 4, 5, 6!

mp.mp.dps = 80

def apply_D(m, p):
    res = [mp.mpf(0)]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def eval_mp(p, tv):
    v = mp.mpf(0)
    for c in reversed(p): v = v * tv + c
    return v

def eval_subset_poly(k, S, b_list=None):
    p = [mp.mpf(1)]
    for j in range(1, k + 1):
        b = mp.mpf(j) if b_list is None else mp.mpf(b_list[j - 1])
        p = apply_D(b, p)
        if j not in S:
            p = apply_D(b, p)
    return p

ts = [mp.mpf(-m) / mp.mpf(100) for m in range(1, 100)] + [mp.mpf('-0.999'), mp.mpf('-1e-4'), mp.mpf('-1e-8'), mp.mpf('-1e-12')]

print("=== 1. TEST BOOLEAN-LATTICE RAYLEIGH AT k = 14, 20 AND WILD b_j! ===")
cases = [
    ("GLZ k=14, a=1, b=14, R=empty", 14, set(), 1, 14, None),
    ("GLZ k=14, a=1, b=2,  R=empty", 14, set(), 1, 2, None),
    ("GLZ k=14, a=1, b=7,  R={2,4,6,8,10,12}", 14, {2,4,6,8,10,12}, 1, 7, None),
    ("GLZ k=14, a=2, b=3,  R={1,4,5,6,7,8,9,10,11,12,13,14}", 14, {1,4,5,6,7,8,9,10,11,12,13,14}, 2, 3, None),
    ("GLZ k=20, a=1, b=2,  R=empty", 20, set(), 1, 2, None),
    ("GLZ k=20, a=1, b=10, R={2,4,6,8,12,14,16,18}", 20, {2,4,6,8,12,14,16,18}, 1, 10, None),
    ("WILD k=8, b=[2,2,1e5,2,1e5,2,1e5,2], a=1, b=2", 8, set(), 1, 2, [2,2,1e5,2,1e5,2,1e5,2]),
    ("WILD k=8, b=[1.01,1e5,1.01,1e5,1.01,1e5,1.01,1e5], a=1, b=2", 8, set(), 1, 2, [1.01,1e5,1.01,1e5,1.01,1e5,1.01,1e5]),
    ("WILD k=8, b=[1e5,2,2,2,2,2,2,1e5], a=1, b=2", 8, set(), 1, 2, [1e5,2,2,2,2,2,2,1e5]),
]
for name, k, R, a, b, blist in cases:
    pa = eval_subset_poly(k, R | {a}, blist)
    pb = eval_subset_poly(k, R | {b}, blist)
    p0 = eval_subset_poly(k, R, blist)
    pab = eval_subset_poly(k, R | {a, b}, blist)
    mn = mp.mpf('1e100'); loc = None
    for tv in ts:
        va, vb, v0, vab = eval_mp(pa, tv), eval_mp(pb, tv), eval_mp(p0, tv), eval_mp(pab, tv)
        d = va*vb - v0*vab
        rel = d / (abs(va*vb) + abs(v0*vab) + mp.mpf('1e-300'))
        if rel < mn: mn = rel; loc = float(tv)
    print(f"  {name:56s}: min_rel = {float(mn):+.6f} at t={loc}")

print("\n=== 2. TEST NON-ADJACENT EXCHANGE h_S h_T - h_U h_V ON (-1, 0) AT k = 4, 6 ===")
# k=4, i=2: U={1}, V={2,3,4}. Pairs (S, T) by moving x in V to U: ({1,2},{3,4}), ({1,3},{2,4}), ({1,4},{2,3})
for k, U, V in [
    (4, {1}, {2,3,4}),
    (4, {4}, {1,2,3}),
    (4, {2}, {1,3,4}),
    (6, {1,2}, {3,4,5,6}),
    (6, {5,6}, {1,2,3,4}),
    (6, {1,6}, {2,3,4,5}),
]:
    pU = eval_subset_poly(k, U)
    pV = eval_subset_poly(k, V)
    print(f"  k={k}, U={sorted(U)}, V={sorted(V)}:")
    for x in sorted(V - U):
        S = U | {x}; T = V - {x}
        pS = eval_subset_poly(k, S); pT = eval_subset_poly(k, T)
        mn = mp.mpf('1e100'); loc = None
        for tv in ts:
            d = eval_mp(pS, tv)*eval_mp(pT, tv) - eval_mp(pU, tv)*eval_mp(pV, tv)
            sc = abs(eval_mp(pS, tv)*eval_mp(pT, tv)) + abs(eval_mp(pU, tv)*eval_mp(pV, tv)) + mp.mpf('1e-300')
            if d/sc < mn: mn = d/sc; loc = float(tv)
        print(f"    move {x} -> S={sorted(S)}, T={sorted(T)}: min_rel = {float(mn):+.6f} at t={loc}")
    # Also check AVERAGE over all x in V\U: sum_{x in V\U} (h_{U u x} h_{V \ x}) - |V\U| h_U h_V!
    diffs = []
    for tv in ts:
        lhs = sum(eval_mp(eval_subset_poly(k, U | {x}), tv)*eval_mp(eval_subset_poly(k, V - {x}), tv) for x in V - U)
        rhs = len(V - U) * eval_mp(pU, tv) * eval_mp(pV, tv)
        diffs.append(float((lhs - rhs) / (abs(lhs) + abs(rhs) + mp.mpf('1e-300'))))
    print(f"    AVG over all x in V\\U: min_rel = {min(diffs):+.6f}")
