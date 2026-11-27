"""Inspect:
1. ALL 6 roots of W_XM(j), W_XY(j), W_MY(j) in the j-plane for GLZ (k = 6, 10, 14) across t in R!
   Check which roots are real, which are complex, their exact distances to j_0 = k + 1,
   AND check the sign of (alpha_1 + 3*alpha_3 - 4*alpha_2)!!
2. Check min_{t in R} J(h_{k, i-2}, D_j h_{k, i-1}, D_j^2 h_{k, i}) as a function of j in [0, k+2]!
   At what exact threshold r_J(k) does J >= 1 first hold on ALL t in R?!
3. Quick check on @algebraist's Point #4: P_L = L - x d/dx on GLZ (k = 2..12):
   Does ALL 4 cross-pairs of B_{k, i} vs P_L B_{k, i-1} hold on GLZ?!
   (Specifically X2 < Y1: A_{k-1, i-1} < P_L S_{L-1} A_{k-1, i-1}, AND X1 < Y2: S_{L-2} A_{k-1, i} < P_L A_{k-1, i-2}!)
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 50

def apply_D_int(m, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly): res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def apply_S_int(L, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly): res[r] += c*(r+1); res[r+1] += c*(L-r)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_P_int(c, poly):
    return [(c - r) * coef for r, coef in enumerate(poly)]

def poly_add(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv_int(p): return [r*p[r] for r in range(1, len(p))]

def eval_mp(p, tv):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * tv + c
    return acc

KMAX = 14
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D_int(m, H[k-1][0])
    for i in range(1, k): U[i] = poly_add(apply_D_int(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1): H[k][i] = apply_D_int(m, U[i])

print("=== 1. ALL 6 ROOTS IN j OF W_XM(j) [deg 1], W_XY(j) [deg 2], W_MY(j) [deg 3] ===")
for k, i in [(6, 4), (10, 6), (14, 8)]:
    f0, f1, f2 = H[k][i-2], H[k][i-1], H[k][i]
    # D_j f = D_0 f + j*(t f)
    D0f1 = apply_D_int(0, f1); tf1 = [0] + list(f1)
    # D_j^2 f2 = D_0^2 f2 + j*(2 t D_0 f2 + t(1+2t) f2) + j^2*(t^2 f2)
    D0f2 = apply_D_int(0, f2); D00f2 = apply_D_int(0, D0f2)
    # Verify: D_j^2 = (D_0 + j t)^2 = D_0^2 + j(D_0(t.) + t D_0(.)) + j^2 t^2
    tf2 = [0] + list(f2)
    D0_tf2 = apply_D_int(0, tf2)
    t_D0f2 = [0] + list(D0f2)
    c1_f2 = poly_add(D0_tf2, t_D0f2)
    c2_f2 = [0, 0] + list(f2)
    # Sanity check at j = 5!
    chk1 = apply_D_int(5, apply_D_int(5, f2))
    chk2 = poly_add(poly_add(D00f2, [5*v for v in c1_f2]), [25*v for v in c2_f2])
    assert chk1 == chk2

    print(f"\n--- k={k}, i={i} (GLZ uses j_0 = {k+1} -> {k+2}) ---")
    for tv_str in ["-2", "-0.5", "-0.05", "0", "0.05", "0.3", "2"]:
        tv = mp.mpf(tv_str)
        def ev(p): return eval_mp(p, tv), eval_mp(deriv_int(p), tv)
        def W(a, b): return a[1]*b[0] - a[0]*b[1]
        v0 = ev(f0)
        m0, m1 = ev(D0f1), ev(tf1)
        y0, y1, y2 = ev(D00f2), ev(c1_f2), ev(c2_f2)
        # W_XM(j) = W(v0, m0) + j W(v0, m1)
        wxm = [float(W(v0, m0)), float(W(v0, m1))]
        # W_XY(j) = W(v0, y0) + j W(v0, y1) + j^2 W(v0, y2)
        wxy = [float(W(v0, y0)), float(W(v0, y1)), float(W(v0, y2))]
        # W_MY(j) = W(m0 + j m1, y0 + j y1 + j^2 y2)
        wmy = [
            float(W(m0, y0)),
            float(W(m0, y1) + W(m1, y0)),
            float(W(m0, y2) + W(m1, y1)),
            float(W(m1, y2)),
        ]
        r_xm = np.roots(wxm[::-1])
        r_xy = np.roots(wxy[::-1])
        r_my = np.roots(wmy[::-1])
        def fmt(rs):
            return "[" + ", ".join(f"{r.real:+.2f}" if abs(r.imag)<1e-6 else f"{r.real:+.2f}{r.imag:+.2f}i" for r in sorted(rs, key=lambda z: -z.real)) + "]"
        # Compute alpha_1, 2*alpha_2, 3*alpha_3 at j0 = k+1!
        j0 = k + 1
        a1 = sum(1.0/(j0 - r) for r in r_xm).real
        a2_2 = sum(1.0/(j0 - r) for r in r_xy).real
        a3_3 = sum(1.0/(j0 - r) for r in r_my).real
        dlogJ = a1 + a3_3 - 2*a2_2
        print(f"  t={tv_str:5s}: r_XM={fmt(r_xm):9s} r_XY={fmt(r_xy):19s} r_MY={fmt(r_my):28s} | d(log J)/dj={dlogJ:+.5f}")

print("\n=== 2. CHECK @algebraist's P_L CROSS-PAIRS ON GLZ (k = 2..10) ===")
A = [[None]*(k+1) for k in range(11)]
A[0][0] = [1]
for k in range(1, 11):
    for i in range(k+1):
        u = [0]
        if i <= k-1: u = poly_add(u, apply_S_int(3*k-3-i, A[k-1][i]))
        if i >= 1: u = poly_add(u, A[k-1][i-1])
        A[k][i] = apply_S_int(3*k-2-i, u)

def roots_x(p):
    d = len(p) - 1
    if d == 0: return []
    g = (mp.mpf(p[0])/mp.mpf(p[-1]))**(mp.mpf(1)/d)
    cs = [mp.mpf(p[r])*(g**r) for r in range(d+1)]
    rts = mp.polyroots(cs[::-1], maxsteps=400, extraprec=50)
    assert max(abs(mp.im(r)) for r in rts) < mp.mpf('1e-20')
    return sorted([mp.re(r)*g for r in rts], reverse=True)  # 0 > r_1 > r_2 > ...

def interlaces_MW(rf, rg):
    """Ma-Wang f < g (deg g = deg f + 1): 0 > s1 > r1 > s2 > r2 > ... > s_d."""
    assert len(rg) == len(rf) + 1
    for j in range(len(rf)):
        if not (rg[j] > rf[j] > rg[j+1]):
            return False
    return True

p_ih_fails = []; x2_y1_fails = []; x1_y2_fails = []
for k in range(2, 11):
    for i in range(1, k + 1):
        L = 3*k - 1 - i
        # Check IH: A_{k,i} < P_L A_{k,i-1}
        rf = roots_x(A[k][i])
        rg = roots_x(apply_P_int(L, A[k][i-1]))
        if not interlaces_MW(rf, rg): p_ih_fails.append((k, i))
        # Check X2 < Y1: A_{k-1, i-1} < P_L S_{L-1} A_{k-1, i-1}
        rx2 = roots_x(A[k-1][i-1])
        ry1 = roots_x(apply_P_int(L, apply_S_int(L-1, A[k-1][i-1])))
        if not interlaces_MW(rx2, ry1): x2_y1_fails.append((k, i))
        # Check X1 < Y2 (for 2 <= i <= k-1): S_{L-2} A_{k-1, i} < P_L A_{k-1, i-2}
        if 2 <= i <= k - 1:
            rx1 = roots_x(apply_S_int(L-2, A[k-1][i]))
            ry2 = roots_x(apply_P_int(L, A[k-1][i-2]))
            if not interlaces_MW(rx1, ry2): x1_y2_fails.append((k, i))
print(f"  IH (A_{{k,i}} < P_L A_{{k,i-1}}) fails (k=2..10): {p_ih_fails}")
print(f"  X2 < Y1 (A_{{k-1,i-1}} < P_L S_{{L-1}} A_{{k-1,i-1}}) fails (k=2..10): {x2_y1_fails}")
print(f"  X1 < Y2 (S_{{L-2}} A_{{k-1,i}} < P_L A_{{k-1,i-2}}) fails (k=2..10): {x1_y2_fails}")
