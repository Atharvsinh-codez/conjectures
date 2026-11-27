import mpmath as mp

# Check GLZ plain Turan AND raw J_raw(t) up to k = 80!
# Specifically:
# 1. Where does raw J_raw(tilde_h_{k,i+1}, tilde_h_{k,i}, tilde_h_{k,i-1})(t) < 1 live on (-1, 0)?
#    Does it ONLY occur for t in (-1/(k+1), 0) (where Phi = 1 + (k+1)t > 0),
#    OR does it also occur for t < -1/(k+1) (where Phi < 0)?!
# 2. Check the EXACT signs of A(t), B(t), C(t) in (D_{k+1} h_i)^2 - (D_{k+1} h_{i-1})(D_{k+1} h_{i+1}) = A Phi^2 + B Phi + C!
# 3. Check whether plain GLZ Turan tilde_h_{k,i}(t)^2 >= tilde_h_{k,i-1}(t) tilde_h_{k,i+1}(t) holds through k = 70!
# 4. Also check if GLZ has ANY sign-flip in (D_{k+1} h_i)^2 - (D_{k+1} h_{i-1})(D_{k+1} h_{i+1}) on (-1, 0) up to k = 70!

mp.mp.dps = 80

def apply_D_int(m, p):
    res = [0]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def deriv_int(p): return [r*p[r] for r in range(1, len(p))] if len(p) > 1 else [0]

def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v

def eval_W(f, g, t):
    return eval_mp(f, t)*eval_mp(deriv_int(g), t) - eval_mp(deriv_int(f), t)*eval_mp(g, t)

# Build GLZ tilde_h row by row up to KMAX = 72! (Only keep current row to save memory!)
KMAX = 72
row = [[1, 2], [1]]  # k = 1

check_ks = [20, 35, 45, 55, 65, 72]

for k in range(2, KMAX + 1):
    m = k + 1
    new_row = [None]*(k + 1)
    new_row[0] = apply_D_int(m, apply_D_int(m, row[0]))
    for i in range(1, k):
        new_row[i] = apply_D_int(m, poly_add_int(apply_D_int(m, row[i]), row[i-1]))
    new_row[k] = apply_D_int(m, row[k-1])
    row = new_row

    if k in check_ks:
        m_next = k + 2  # next-step operator D_{k+2} acting on tilde_h_{k, i}!
        # Grid covering BOTH (-1, -1/m_next) and (-1/m_next, 0)!
        t_left = [mp.mpf(s) for s in ["-0.9999", "-0.99", "-0.95", "-0.9", "-0.8", "-0.6", "-0.4", "-0.2", "-0.1"]]
        t_left += [-mp.mpf(s)/m_next for s in ["5.0", "2.0", "1.5", "1.1", "1.01", "1.0"]]
        t_left = [t for t in t_left if t > -1 and t <= -mp.mpf(1)/m_next]
        t_right = [-mp.mpf(s)/m_next for s in ["0.99", "0.9", "0.5", "0.2", "0.05", "1e-3", "1e-6", "1e-12", "1e-20"]]
        t_right.append(mp.mpf(0))

        # Check sample i's: 1, 2, 3, k//4, k//2, 3*k//4, k-2, k-1
        is_chk = sorted(set([1, 2, 3, k//4, k//2, (k+1)//2, 3*k//4, k-2, k-1]))
        is_chk = [i for i in is_chk if 1 <= i <= k - 1]

        min_tur = (mp.mpf('1e100'), None)
        min_Dtur = (mp.mpf('1e100'), None)
        min_J_left = (mp.mpf('1e100'), None)
        min_J_right = (mp.mpf('1e100'), None)
        min_B_right = (mp.mpf('1e100'), None)
        min_C_right = (mp.mpf('1e100'), None)

        # Precompute D_{m_next} row[i] only for needed indices
        need = sorted(set(x for i in is_chk for x in (i-1, i, i+1)))
        Drow = {r: apply_D_int(m_next, row[r]) for r in need}
        drow = {r: deriv_int(row[r]) for r in need}

        for i in is_chk:
            for tag, tlist in [("L", t_left), ("R", t_right)]:
                for tv in tlist:
                    v0, v1, v2 = eval_mp(row[i-1], tv), eval_mp(row[i], tv), eval_mp(row[i+1], tv)
                    p0, p1, p2 = eval_mp(drow[i-1], tv), eval_mp(drow[i], tv), eval_mp(drow[i+1], tv)
                    # Plain Turan
                    A = v1*v1 - v0*v2
                    sc = v1*v1 + abs(v0*v2) + mp.mpf('1e-300')
                    if A/sc < min_tur[0]: min_tur = (A/sc, (i, float(tv)))
                    # D_{k+2} Turan
                    dv0, dv1, dv2 = eval_mp(Drow[i-1], tv), eval_mp(Drow[i], tv), eval_mp(Drow[i+1], tv)
                    dA = dv1*dv1 - dv0*dv2
                    dsc = dv1*dv1 + abs(dv0*dv2) + mp.mpf('1e-300')
                    if dA/dsc < min_Dtur[0]: min_Dtur = (dA/dsc, (i, float(tv)))
                    # J_raw
                    w1 = v1*p0 - p1*v0  # W(h_i, h_{i-1}) > 0
                    w2 = v2*p1 - p2*v1  # W(h_{i+1}, h_i) > 0
                    w12 = v2*p0 - p2*v0 # W(h_{i+1}, h_{i-1}) > 0
                    J = 4*w1*w2/(w12*w12) if w12 != 0 else mp.mpf('1e10')
                    if tag == "L" and J < min_J_left[0]: min_J_left = (J, (i, float(tv)))
                    if tag == "R":
                        if J < min_J_right[0]: min_J_right = (J, (i, float(tv)))
                        # Check B and C on right (where Phi > 0 and t < 0)!
                        if tv < 0:
                            B = tv*(1+tv)*(2*v1*p1 - v0*p2 - p0*v2)
                            C = (tv*(1+tv))**2 * (p1*p1 - p0*p2)
                            Brel = B / (abs(tv*(1+tv))*(abs(v1*p1)+abs(v0*p2)+abs(p0*v2)) + mp.mpf('1e-300'))
                            Crel = C / ((tv*(1+tv))**2*(p1*p1+abs(p0*p2)) + mp.mpf('1e-300'))
                            if Brel < min_B_right[0]: min_B_right = (Brel, (i, float(tv), float(J)))
                            if Crel < min_C_right[0]: min_C_right = (Crel, (i, float(tv)))

        print(f"k={k:2d} (-1/m={-1/m_next:.4f}):")
        print(f"  plain_Tur >= {float(min_tur[0]):+.5f}@{min_tur[1]} | D_{{k+2}}_Tur >= {float(min_Dtur[0]):+.5f}@{min_Dtur[1]}")
        print(f"  J_raw on (-1, -1/m] >= {float(min_J_left[0]):.5f}@{min_J_left[1]} | on (-1/m, 0] >= {float(min_J_right[0]):.5f}@{min_J_right[1]}")
        print(f"  On (-1/m, 0): min_B_rel = {float(min_B_right[0]):+.5f}@{min_B_right[1]}, min_C_rel = {float(min_C_right[0]):+.5f}@{min_C_right[1]}")
