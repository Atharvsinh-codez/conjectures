"""HIGH-PRECISION (dps=200) CHECK OF TWO-POLE APPROXIMATION FOR psi_{k,i}(s) AT k = 6, 10, 15, 20!
Notice: psi_{k,i}(pi*i + u) = -e^u / (1 - e^u)^{D+1} * A_{k,i}(e^u) = sum_{r=1}^{D+1} c_r u^{-r} + O(1).
At dps=200 there is ZERO precision loss even at k=20 (D=60)!
Also inspect the coefficients c_r (do they alternate? are they smooth? is 2-pole accurate EVERYWHERE on R?!).
"""
import mpmath as mp

mp.mp.dps = 200

def apply_S(L, poly):
    res = [0]*(len(poly)+1)
    for r, c in enumerate(poly):
        res[r] += c*(r+1); res[r+1] += c*(L-r)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def poly_add(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

KMAX = 20
A = [[None]*(k+1) for k in range(KMAX+1)]
A[0][0] = [1]
for k in range(1, KMAX+1):
    for i in range(k+1):
        u = [0]
        if i <= k-1: u = poly_add(u, apply_S(3*k-3-i, A[k-1][i]))
        if i >= 1: u = poly_add(u, A[k-1][i-1])
        A[k][i] = apply_S(3*k-2-i, u)

def eval_poly_mp(p, x):
    acc = mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

def eval_psi(k, i, s):
    es = mp.exp(s)
    return es / ((1 + es)**(3*k - i + 1)) * eval_poly_mp(A[k][i], -es)

def roots_in_s(k, i):
    p = A[k][i]; d = len(p) - 1
    if d == 0: return []
    g = (mp.mpf(p[0]) / mp.mpf(p[-1])) ** (mp.mpf(1) / d)
    cs = [mp.mpf(p[r]) * (g ** r) for r in range(d + 1)]
    rts = mp.polyroots(cs[::-1], maxsteps=400, extraprec=50)
    xs = sorted([mp.re(r) * g for r in rts])
    return sorted([mp.log(-x) for x in xs])

def laurent_principal_part(k, i):
    D = 3*k - i; M = D + 1
    h = [mp.mpf(1) / mp.factorial(m + 1) for m in range(M)]
    g = [mp.mpf(0)] * M; g[0] = mp.mpf(-1)
    for m in range(1, M):
        g[m] = -sum(h[j] * g[m - j] for j in range(1, m + 1))
    def mul_trunc(a, b):
        res = [mp.mpf(0)] * M
        for p in range(M):
            if a[p] == 0: continue
            for q in range(M - p): res[p + q] += a[p] * b[q]
        return res
    gp = [mp.mpf(0)] * M; gp[0] = mp.mpf(1)
    base = list(g); expn = D + 1
    while expn:
        if expn & 1: gp = mul_trunc(gp, base)
        expn >>= 1
        if expn: base = mul_trunc(base, base)
    p = A[k][i]
    num = [mp.mpf(0)] * M
    for m in range(M):
        s = sum(mp.mpf(coef) * (mp.mpf(r + 1) ** m) for r, coef in enumerate(p))
        num[m] = -s / mp.factorial(m)
    prod = mul_trunc(gp, num)
    c = [mp.mpf(0)] * (D + 2)
    for m in range(M): c[D + 1 - m] = prod[m]
    return c

def eval_two_pole(c, s):
    w = mp.mpc(1) / mp.mpc(s, -mp.pi)
    val = mp.mpc(0); wp = w
    for r in range(1, len(c)):
        val += c[r] * wp; wp *= w
    return 2 * mp.re(val)

# Sanity check: Mittag-Leffler identity at k=4, i=2, s=0.3: sum_{|l|<=5} Q(1/(s - (2l+1)pi*i)) == psi(s)!
c42 = laurent_principal_part(4, 2)
s0 = mp.mpf("0.3")
ex42 = eval_psi(4, 2, s0)
ml42 = sum(sum(c42[r] / (mp.mpc(s0, -(2*l+1)*mp.pi)**r) for r in range(1, len(c42))) for l in range(-20, 21))
print(f"[Mittag-Leffler sanity k=4,i=2,s=0.3: exact={mp.nstr(ex42,15)}, ML={mp.nstr(mp.re(ml42),15)}, diff={mp.nstr(abs(ex42-ml42),4)}]")

for k in [6, 10, 15, 20]:
    print(f"\n--- k = {k} ---")
    for i in sorted(set([0, min(4, k), k//2, max(0, k-3), k])):
        ss = roots_in_s(k, i); d = len(ss)
        c = laurent_principal_part(k, i)
        pts = [ss[0] - 0.5] + [0.5*(ss[j] + ss[j+1]) for j in range(d - 1)] + [ss[-1] + 0.5]
        ok_cnt = 0; bad_pts = []
        max_err_mid = mp.mpf(0)
        for idx, sv in enumerate(pts):
            ex = eval_psi(k, i, sv)
            tp = eval_two_pole(c, sv)
            err = abs(ex - tp) / (abs(ex) + mp.mpf('1e-300'))
            if (ex > 0 and tp > 0) or (ex < 0 and tp < 0):
                ok_cnt += 1
            else:
                bad_pts.append((idx, round(float(sv), 2)))
            if -5 <= sv <= 5 and err > max_err_mid:
                max_err_mid = err
        print(f"  i={i:2d} (d={d:2d}): 2pole ok = {ok_cnt:2d}/{d+1:2d}, max_err(|s|<=5) = {mp.nstr(max_err_mid, 3)}, bad_pts = {bad_pts[:6]}", flush=True)
