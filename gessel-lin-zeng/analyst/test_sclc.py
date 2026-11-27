"""EUREKA TEST OF SYMMETRIZED COEFFICIENTWISE LOG-CONCAVITY (SCLC):
For polynomials f, g, u, v, define:
  SCLC(f, g; u, v)  <=>  f_p g_q + f_q g_p >= u_p v_q + u_q v_p   for ALL p, q >= 0!

WE JUST PROVED ALGEBRAICALLY THAT:
  (1) S_L (in x-space, for L >= deg) and D_m (in t-space, for m > 0) PRESERVE SCLC EXACTLY!!
  (2) Pascal convolution u_i = v_i + v_{i-1} PRESERVES strong SCLC (all i <= j) EXACTLY!!
  (3) SCLC(M, M; X, Y) IMMEDIATELY IMPLIES M(x)^2 - X(x)Y(x) >=_x 0!!
      Wait: does M(x)^2 - X(x)Y(x) >=_x 0 imply M(x)^2 > X(x)Y(x) on x < 0?!
      CAREFUL! Positive x-coeffs imply > 0 on x >= 0, NOT automatically on x < 0!
      WAIT: check BOTH x-space AND t-space (and what SCLC holds for!)!

Test RIGHT NOW:
1. Does strong SCLC hold for g_{k, i} (aligned) in x-space AND in t-space for k = 2..14?
2. Does SCLC(M, M; X, Y) hold in x-space AND in t-space for k = 3..14?
3. Does strong x-LC (g_i g_j - g_{i-1} g_{j+1} >=_x 0) hold for k = 2..14?
4. ALSO check: in t-space, does SCLC hold for (htilde_{k,i}, htilde_{k,j}; htilde_{k,i-1}, htilde_{k,j+1}) (UNALIGNED)?
5. ALSO check GLZ paired multi-affine (1 variable u_r per step r) just to close out our previous question!
"""
import math

def apply_S(L, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r+1); res[r+1] += c*(L-r)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_D(m, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

def poly_mul(a, b):
    res = [0]*(len(a)+len(b)-1)
    for i, va in enumerate(a):
        if va == 0: continue
        for j, vb in enumerate(b): res[i+j] += va*vb
    return res

def check_sclc(f, g, u, v):
    """Return (ok, worst_pq_val)."""
    d = max(len(f), len(g), len(u), len(v))
    def get(p, idx): return p[idx] if idx < len(p) else 0
    for p in range(d):
        for q in range(p, d):
            val = get(f, p)*get(g, q) + get(f, q)*get(g, p) - get(u, p)*get(v, q) - get(u, q)*get(v, p)
            if val < 0:
                return False, (p, q, val)
    return True, None

def check_poly_ge(f, g, u, v):
    """Check f*g - u*v >= 0 coeffwise."""
    fg = poly_mul(f, g); uv = poly_mul(u, v)
    n = max(len(fg), len(uv))
    for r in range(n):
        val = (fg[r] if r < len(fg) else 0) - (uv[r] if r < len(uv) else 0)
        if val < 0:
            return False, (r, val)
    return True, None

KMAX = 12
A = [[None]*(k+1) for k in range(KMAX+1)]
A[0][0] = [1]
for k in range(1, KMAX+1):
    for i in range(k+1):
        u = [0]
        if i <= k-1: u = poly_add(u, apply_S(3*k-3-i, A[k-1][i]))
        if i >= 1: u = poly_add(u, A[k-1][i-1])
        A[k][i] = apply_S(3*k-2-i, u)

H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    U = [None]*(k+1)
    U[0] = apply_D(m, H[k-1][0])
    for i in range(1, k): U[i] = poly_add(apply_D(m, H[k-1][i]), H[k-1][i-1])
    U[k] = H[k-1][k-1]
    for i in range(k+1): H[k][i] = apply_D(m, U[i])

print("=== 1. STRONG SCLC & STRONG LC FOR ALIGNED g_{k,i} IN x-SPACE (k=2..12) ===")
for k in range(2, KMAX+1):
    Gx = []
    for i in range(k+1):
        p = A[k][i]
        for s in range(i): p = apply_S(3*k - i - 1 + s, p)
        Gx.append(p)
    sclc_fails = []; lc_fails = []
    for i in range(1, k):
        for j in range(i, k):
            ok, info = check_sclc(Gx[i], Gx[j], Gx[i-1], Gx[j+1])
            if not ok: sclc_fails.append((i, j, info[:2]))
            ok2, info2 = check_poly_ge(Gx[i], Gx[j], Gx[i-1], Gx[j+1])
            if not ok2: lc_fails.append((i, j, info2[0]))
    print(f"  k={k:2d} (x-space aligned): SCLC fails={len(sclc_fails)} {sclc_fails[:3]}, strong x-LC fails={len(lc_fails)} {lc_fails[:3]}")

print("\n=== 2. STRONG SCLC & STRONG LC FOR ALIGNED g_{k,i} IN t-SPACE (k=2..12) ===")
for k in range(2, KMAX+1):
    Gt = []
    for i in range(k+1):
        p = list(H[k][i])
        for _ in range(i): p = apply_D(k+1, p)
        Gt.append(p)
    sclc_fails = []; lc_fails = []
    for i in range(1, k):
        for j in range(i, k):
            ok, info = check_sclc(Gt[i], Gt[j], Gt[i-1], Gt[j+1])
            if not ok: sclc_fails.append((i, j, info[:2]))
            ok2, info2 = check_poly_ge(Gt[i], Gt[j], Gt[i-1], Gt[j+1])
            if not ok2: lc_fails.append((i, j, info2[0]))
    print(f"  k={k:2d} (t-space aligned): SCLC fails={len(sclc_fails)} {sclc_fails[:3]}, strong t-LC fails={len(lc_fails)} {lc_fails[:3]}")

print("\n=== 3. SCLC(M, M; X, Y) IN x-SPACE AND t-SPACE (k=3..12) ===")
for k in range(3, KMAX+1):
    m = k + 1
    fx = []; ft = []
    for i in range(2, k):
        Xx = A[k-1][i-2]
        Mx = apply_S(3*k - i - 2, A[k-1][i-1])
        Yx = apply_S(3*k - i - 2, apply_S(3*k - i - 3, A[k-1][i]))
        ok, info = check_sclc(Mx, Mx, Xx, Yx)
        if not ok: fx.append((i, info[:2]))

        Xt = H[k-1][i-2]
        Mt = apply_D(m, H[k-1][i-1])
        Yt = apply_D(m, apply_D(m, H[k-1][i]))
        okt, infot = check_sclc(Mt, Mt, Xt, Yt)
        if not okt: ft.append((i, infot[:2]))
    print(f"  k={k:2d}: x-SCLC(M,M;X,Y) fails={len(fx)} {fx[:3]}, t-SCLC(M,M;X,Y) fails={len(ft)} {ft[:3]}")

print("\n=== 4. CHECK STRONG SCLC IN t-SPACE FOR UNALIGNED htilde_{k,i} (k=2..10) ===")
for k in range(2, 11):
    sclc_fails = []; lc_fails = []
    for i in range(1, k):
        for j in range(i, k):
            ok, info = check_sclc(H[k][i], H[k][j], H[k][i-1], H[k][j+1])
            if not ok: sclc_fails.append((i, j, info[:2]))
            ok2, info2 = check_poly_ge(H[k][i], H[k][j], H[k][i-1], H[k][j+1])
            if not ok2: lc_fails.append((i, j, info2[0]))
    print(f"  k={k:2d} (t-space UNALIGNED): SCLC fails={len(sclc_fails)} {sclc_fails[:3]}, strong t-LC fails={len(lc_fails)} {lc_fails[:3]}")
