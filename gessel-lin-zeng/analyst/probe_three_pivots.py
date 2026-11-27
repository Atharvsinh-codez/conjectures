import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 60

def apply_D_int(m, p):
    res = [0]*(len(p) + 1)
    for r, c in enumerate(p):
        res[r] += c*(r+1)
        res[r+1] += c*(r+m)
    return res

def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv_int(p):
    return [r*p[r] for r in range(1, len(p))] if len(p) > 1 else [0]

def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v

def eval_W_mp(f, g, t):
    return eval_mp(f, t)*eval_mp(deriv_int(g), t) - eval_mp(deriv_int(f), t)*eval_mp(g, t)

# Build GLZ tilde_h[k][i] up to k=16
KMAX = 16
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    H[k][0] = apply_D_int(m, apply_D_int(m, H[k-1][0]))
    for i in range(1, k):
        H[k][i] = apply_D_int(m, poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1]))
    H[k][k] = apply_D_int(m, H[k-1][k-1])

print("=== 1. FULL-ROW PAIRWISE INTERLACING W(E_{k,i}, E_{k,j}) > 0 FOR ALL 0 <= i < j <= k-1 ===")
# Recall E_{k,i} = D_{k+1}^i tilde_h_{k-1, i} (all degree 2k-3)
ts_check = [mp.mpf(s) for s in ["-100", "-10", "-2", "-1.1", "-1.01", "-1", "-0.99", "-0.9", "-0.5", "-0.1", "-0.01", "-1e-5", "0", "1e-5", "0.01", "0.1", "1", "10", "100"]]
for k in range(3, 16):
    m = k + 1
    E = []
    for i in range(k):
        p = H[k-1][i]
        for _ in range(i):
            p = apply_D_int(m, p)
        E.append(p)
    neg_pairs = []
    min_w_ratio = mp.mpf('1e100')
    for i in range(k):
        for j in range(i+1, k):
            for tv in ts_check:
                w = eval_W_mp(E[i], E[j], tv)
                if w <= 0:
                    neg_pairs.append((i, j, float(tv), float(w)))
                    break
    print(f"  k={k:2d}: all-pairs W<=0 count = {len(neg_pairs):3d}/{k*(k-1)//2:3d}, max_passing_r = {min(j-i for i,j,_,_ in neg_pairs)-1 if neg_pairs else k-1}" + (f", first_fail={neg_pairs[0][:3]}" if neg_pairs else ""))

print("\n=== 2. EXACT STAGE 3 FORMULA FOR Sigma(t) = coeff of 1/j^2 in d/dj log J(j; t) ON GLZ ===")
# For triple (h0, h1, h2) = (tilde_h_{k, i-1}, tilde_h_{k, i}, tilde_h_{k, i+1}),
# M(j) = D_j h1, Y(j) = D_j^2 h2, J(j; t) = 4 W(h0, M(j)) W(M(j), Y(j)) / W(h0, Y(j))^2.
# Compute exact W_XM(j), W_XY(j), W_MY(j) as polynomials in j at sample t values,
# and compute Sigma(t) = r_XM + sum(r_MY) - 2*sum(r_XY), AND exact d/dj log J(j; t) at j = k+1!
def poly_mul_t(p, q):
    res = [0]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            res[i+j] += a*b
    return res

def shift_t(p): return [0] + p

for k in [5, 8, 11, 14]:
    j0 = mp.mpf(k + 1)
    worst_dlogJ = (mp.mpf('1e100'), None)
    worst_Sigma = (mp.mpf('1e100'), None)
    max_Sigma = (-mp.mpf('1e100'), None)
    for i in range(1, k):
        h0, h1, h2 = H[k][i-1], H[k][i], H[k][i+1]
        # M(j) = D_0 h1 + j * (t h1)
        D0h1 = apply_D_int(0, h1)
        th1 = shift_t(h1)
        # Y(j) = D_0^2 h2 + j * (2 t D_0 h2 + t(1+t) h2) + j^2 * (t^2 h2)
        D0h2 = apply_D_int(0, h2)
        Y0 = apply_D_int(0, D0h2)
        Y1 = poly_add_int([2*c for c in shift_t(D0h2)], poly_mul_t([0, 1, 1], h2))
        Y2 = shift_t(shift_t(h2))
        for tv in [mp.mpf(s) for s in ["-10", "-2", "-1.1", "-0.9", "-0.5", "-0.1", "-0.01", "0.01", "0.1", "1", "10"]]:
            # W_XM(j) = a0 + a1*j
            a0 = eval_W_mp(h0, D0h1, tv); a1 = eval_W_mp(h0, th1, tv)
            # W_XY(j) = b0 + b1*j + b2*j^2
            b0 = eval_W_mp(h0, Y0, tv); b1 = eval_W_mp(h0, Y1, tv); b2 = eval_W_mp(h0, Y2, tv)
            # W_MY(j) = c0 + c1*j + c2*j^2 + c3*j^3
            c0 = eval_W_mp(D0h1, Y0, tv)
            c1 = eval_W_mp(D0h1, Y1, tv) + eval_W_mp(th1, Y0, tv)
            c2 = eval_W_mp(D0h1, Y2, tv) + eval_W_mp(th1, Y1, tv)
            c3 = eval_W_mp(th1, Y2, tv)
            if a1 == 0 or b2 == 0 or c3 == 0: continue
            # sum of roots: -a0/a1, -b1/b2, -c2/c3
            Sigma = (-a0/a1) + (-c2/c3) - 2*(-b1/b2)
            # exact d/dj log J at j0
            w1 = a0 + a1*j0; w1p = a1
            w2 = c0 + c1*j0 + c2*j0**2 + c3*j0**3; w2p = c1 + 2*c2*j0 + 3*c3*j0**2
            w12 = b0 + b1*j0 + b2*j0**2; w12p = b1 + 2*b2*j0
            dlogJ = w1p/w1 + w2p/w2 - 2*w12p/w12
            if dlogJ < worst_dlogJ[0]:
                worst_dlogJ = (dlogJ, (i, float(tv), float(Sigma)))
            if Sigma < worst_Sigma[0]:
                worst_Sigma = (Sigma, (i, float(tv)))
            if Sigma > max_Sigma[0]:
                max_Sigma = (Sigma, (i, float(tv)))
    print(f"  k={k:2d} (j={k+1:2d}): min d/dj log J = {float(worst_dlogJ[0]):+.6f} (j^2*dlogJ = {float((k+1)**2*worst_dlogJ[0]):+.4f}) at i={worst_dlogJ[1][0]}, t={worst_dlogJ[1][1]} (Sigma={worst_dlogJ[1][2]:+.4f})")
    print(f"            Sigma range across all (i, t): [{float(worst_Sigma[0]):+.4f}, {float(max_Sigma[0]):+.4f}]")

print("\n=== 3. ATTACK Q_{k, q} FOR b_j >= 2, a_j = 1, w_j = 1 (DOUBLE AND SINGLE!) ===")
# Check if b = [2]*m + [1e6]*n breaks Q_{k, q} at q = 1/4 or 1/8!
def apply_D_mp(m, p):
    res = [mp.mpf(0)]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1)
        res[r+1] += c*(r+m)
    return res

def add_mp(a, b):
    res = [mp.mpf(0)]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def test_Q_bge2(name, b_list, doubled=True, qs=(0.5, 0.25, 0.125)):
    P = [[mp.mpf(1)]]
    for b in b_list:
        k = len(P)
        newP = []
        for i in range(k + 1):
            t2 = apply_D_mp(b, P[i]) if i < k else [mp.mpf(0)]
            t1 = P[i-1] if i >= 1 else [mp.mpf(0)]
            u = add_mp(t2, t1)
            newP.append(apply_D_mp(b, u) if doubled else u)
        P = newP
    k = len(b_list)
    # Test fixed-t real-rootedness in y of sum q^{i(i-1)/2} P[i](t) y^i at critical t points (incl near 0, -0.5, -1, 1, -2)!
    t_probes = [mp.mpf(s) for s in ["-10", "-1.5", "-0.99", "-0.5", "-0.1", "-1e-3", "-1e-6", "-1e-10", "0", "1e-6", "1"]]
    res_str = []
    for q in qs:
        qm = mp.mpf(str(q))
        bad = None
        for tv in t_probes:
            coeffs = [qm**(i*(i-1)//2) * eval_mp(P[i], tv) for i in range(k+1)]
            # check y-roots
            try:
                roots = mp.polyroots(coeffs[::-1], maxsteps=200, extraprec=40)
                mx = float(max(abs(mp.im(z)) for z in roots))
                if mx > 1e-12:
                    bad = (float(tv), mx)
                    break
            except Exception:
                pass
        res_str.append(f"q={q}: {'FAIL@t='+str(bad[0]) if bad else 'PASS'}")
    print(f"  {name:38s} (k={k:2d}, {'dbl' if doubled else 'sgl'}): " + " | ".join(res_str))

for doubled in [False, True]:
    test_Q_bge2("[2]*4 + [1e4]*4", [mp.mpf(2)]*4 + [mp.mpf('1e4')]*4, doubled=doubled)
    test_Q_bge2("[2]*6 + [1e6]*6", [mp.mpf(2)]*6 + [mp.mpf('1e6')]*6, doubled=doubled)
    test_Q_bge2("[1e6]*6 + [2]*6", [mp.mpf('1e6')]*6 + [mp.mpf(2)]*6, doubled=doubled)
    test_Q_bge2("[2]*3 + [100]*3 + [1e5]*3", [mp.mpf(2)]*3 + [mp.mpf(100)]*3 + [mp.mpf('1e5')]*3, doubled=doubled)
    test_Q_bge2("[2]*12 + [1e8]*8", [mp.mpf(2)]*12 + [mp.mpf('1e8')]*8, doubled=doubled)
