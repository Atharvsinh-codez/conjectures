import math

def poly_trim(c):
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c

def apply_D(j, c):
    d = len(c) - 1
    out = [0] * (d + 2)
    for m, cm in enumerate(c):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return poly_trim(out)

def add_poly(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return poly_trim(out)

def sub_poly(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] -= b[i]
    return poly_trim(out)

def mul_poly(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return poly_trim(out)

def to_x_poly(h_coeffs):
    d = len(h_coeffs) - 1
    pows = [[1]]
    for r in range(1, d + 1):
        pows.append(mul_poly(pows[-1], [1, -1]))
    out = [0] * (d + 1)
    for m, cm in enumerate(h_coeffs):
        pw = pows[d - m]
        for j, pj in enumerate(pw):
            out[m + j] += cm * pj
    return poly_trim(out)

ht = {(1, 0): [1, 2], (1, 1): [1]}
for k in range(2, 10):
    for i in range(k + 1):
        t1 = apply_D(k + 1, ht[(k-1, i)]) if i <= k - 1 else [0]
        t2 = ht[(k-1, i-1)] if i >= 1 else [0]
        ht[(k, i)] = apply_D(k + 1, add_poly(t1, t2))

print("=== Fast TEST 3: Turan and M^2-XY in x-coeffs and on (-1,0) for WILD b_j ===")
for name, bs in [
    ("GLZ inner [2,3,4,5,6,7]", [2,3,4,5,6,7]),
    ("Const [2]*6", [2]*6),
    ("Wild [2,100,2,100,2,100]", [2,100,2,100,2,100]),
    ("Wild [100,2,100,2,100,2]", [100,2,100,2,100,2]),
    ("Drop [100,50,25,12,6,3]", [100,50,25,12,6,3]),
    ("Spike [2,2,1000,2,2,2]", [2,2,1000,2,2,2]),
]:
    row = [[1]]
    turan_x_pos = True
    turan_t_pos = True
    m2xy_x_pos = True
    m2xy_R_pos = True
    for step, bj in enumerate(bs, 1):
        if step >= 3:
            for i in range(2, step):
                M = apply_D(bj, row[i-1])
                X = row[i-2]
                Y = apply_D(bj, apply_D(bj, row[i]))
                diff_t = sub_poly(mul_poly(M, M), mul_poly(X, Y))
                diff_x = to_x_poly(diff_t)
                if min(diff_x) < 0:
                    m2xy_x_pos = False
                for tv in [-0.99, -0.9, -0.7, -0.5, -0.3, -0.1, -0.01]:
                    val = sum(c * (tv**m) for m, c in enumerate(diff_t))
                    if val < 0:
                        m2xy_R_pos = False
        new_row = []
        for i in range(step + 1):
            t1 = apply_D(bj, row[i]) if i <= step - 1 else [0]
            t2 = row[i - 1] if i >= 1 else [0]
            new_row.append(add_poly(t1, t2))
        row = new_row
        for i in range(1, step):
            tur_t = sub_poly(mul_poly(row[i], row[i]), mul_poly(row[i-1], row[i+1]))
            if min(tur_t) < 0: turan_t_pos = False
            tur_x = to_x_poly(tur_t)
            if min(tur_x) < 0: turan_x_pos = False
    print(f"  {name:25s}: tur_t>=0:{str(turan_t_pos):5s} tur_x>=0:{str(turan_x_pos):5s} | M^2-XY_x>=0:{str(m2xy_x_pos):5s} | M^2-XY(-1,0)>=0:{str(m2xy_R_pos):5s}")

print("\n=== Fast TEST 4: ALL 2x2 minors p_i(x) p_j(x) - p_{i-1}(x) p_{j+1}(x) >= 0 (coeffwise in x and in t) ===")
for name, bs in [
    ("GLZ full D_j(D_j+y) k=8", "GLZ"),
    ("Const [2]*6", [2]*6),
    ("Wild [2,100,2,100,2,100]", [2,100,2,100,2,100]),
    ("Wild [100,2,100,2,100,2]", [100,2,100,2,100,2]),
    ("Drop [100,50,25,12,6,3]", [100,50,25,12,6,3]),
]:
    if bs == "GLZ":
        row = [ht[(8, i)] for i in range(9)]
    else:
        row = [[1]]
        for bj in bs:
            k = len(row)
            row = [add_poly(apply_D(bj, row[i]) if i <= k-1 else [0], row[i-1] if i >= 1 else [0]) for i in range(k+1)]
    xrow = [to_x_poly(p) for p in row]
    ok_2x2_x = True
    ok_2x2_t = True
    for i in range(1, len(xrow)):
        for j in range(i, len(xrow) - 1):
            det_x = sub_poly(mul_poly(xrow[i], xrow[j]), mul_poly(xrow[i-1], xrow[j+1]))
            det_t = sub_poly(mul_poly(row[i], row[j]), mul_poly(row[i-1], row[j+1]))
            if min(det_x) < 0: ok_2x2_x = False
            if min(det_t) < 0: ok_2x2_t = False
    print(f"  {name:25s}: 2x2_x>=0: {str(ok_2x2_x):5s} | 2x2_t>=0: {str(ok_2x2_t):5s}")

print("\n=== Fast TEST 5: Check Wronskian W(p_{i-1}, p_i) and Hankel/Toeplitz minors in t & x for GLZ k=6 and Wild ===")
def deriv(p):
    return poly_trim([m * p[m] for m in range(1, len(p))]) if len(p) > 1 else [0]

for name, bs in [("GLZ k=6", "GLZ"), ("Wild [2,100,2,100,2]", [2,100,2,100,2])]:
    if bs == "GLZ":
        row = [ht[(6, i)] for i in range(7)]
    else:
        row = [[1]]
        for bj in bs:
            k = len(row)
            row = [add_poly(apply_D(bj, row[i]) if i <= k-1 else [0], row[i-1] if i >= 1 else [0]) for i in range(k+1)]
    xrow = [to_x_poly(p) for p in row]
    # Check W(A_{i-1}, A_i) = A_{i-1}' A_i - A_{i-1} A_i' in x!
    w_pos_x = True
    for i in range(1, len(xrow)):
        # A_{i-1}' A_i - A_{i-1} A_i'
        w = sub_poly(mul_poly(deriv(xrow[i-1]), xrow[i]), mul_poly(xrow[i-1], deriv(xrow[i])))
        if min(w) < 0: w_pos_x = False
    # Check 2x2 minors of COEFFICIENT MATRIX M_{i, m} = [x^m] A_i(x) (and [t^m] h_i(t))!
    def check_coeff_tp2(polys):
        # rows = polys[0], polys[1], ... (decreasing degree!)
        # wait: degree of polys[0] > polys[1] > ...
        R = len(polys)
        C = max(len(p) for p in polys)
        M = [[polys[r][c] if c < len(polys[r]) else 0 for c in range(C)] for r in range(R)]
        # Check orientation: reversed rows (increasing degree!)
        Mrev = list(reversed(M))
        neg_rev = 0
        for r1 in range(R):
            for r2 in range(r1+1, R):
                for c1 in range(C):
                    for c2 in range(c1+1, C):
                        if Mrev[r1][c1]*Mrev[r2][c2] - Mrev[r1][c2]*Mrev[r2][c1] < 0:
                            neg_rev += 1
        return neg_rev
    print(f"  {name:22s}: W_x>=0: {w_pos_x} | coeff_matrix_x TP2 negs: {check_coeff_tp2(xrow)} | coeff_matrix_t TP2 negs: {check_coeff_tp2(row)}")
