import sympy as sp
import mpmath as mp

mp.mp.dps = 50
t = sp.Symbol('t')

def apply_Dtilde(j, p):
    return sp.expand(t*(1+t)*sp.diff(p, t) + (1 + (j+1)*t)*p)

KMAX = 10
g = {(1, 0): 1 + 2*t, (1, 1): sp.Integer(1)}
S = {}
for k in range(2, KMAX + 1):
    for i in range(k + 1):
        term1 = apply_Dtilde(k, g[(k-1, i)]) if i <= k - 1 else sp.Integer(0)
        term2 = g[(k-1, i-1)] if i >= 1 else sp.Integer(0)
        S[(k, i)] = sp.expand(term1 + term2)
        g[(k, i)] = apply_Dtilde(k, S[(k, i)])

def roots_mp(expr):
    p = sp.Poly(expr, t)
    if p.degree() == 0:
        return []
    coeffs = [mp.mpf(int(c)) for c in p.all_coeffs()]
    r = mp.polyroots(coeffs, maxsteps=200, extraprec=40)
    assert max(abs(z.imag) for z in r) < mp.mpf('1e-25'), f"Non-real root! {max(abs(z.imag) for z in r)}"
    return sorted([z.real for z in r])

def check_interlace(rp, rq):
    dp, dq = len(rp), len(rq)
    if dq == dp:
        gaps = [rq[m] - rp[m] for m in range(dp)] + [rp[m+1] - rq[m] for m in range(dp - 1)]
        return all(z > 0 for z in gaps), float(min(gaps)) if gaps else 999.0
    elif dq == dp + 1:
        gaps = [rp[m] - rq[m] for m in range(dp)] + [rq[m+1] - rp[m] for m in range(dp)]
        return all(z > 0 for z in gaps), float(min(gaps)) if gaps else 999.0
    return False, -999.0

Rg = {k_i: roots_mp(v) for k_i, v in g.items()}
RS = {k_i: roots_mp(v) for k_i, v in S.items()}
RD = {(k, i): roots_mp(apply_Dtilde(k, g[(k-1, i)])) for k in range(2, KMAX+1) for i in range(k)}

def test_rel(name, pairs):
    fails = []
    mg = 999.0
    mat = None
    for tag, rp, rq in pairs:
        ok, gap = check_interlace(rp, rq)
        if gap < mg:
            mg = gap
            mat = tag
        if not ok:
            fails.append((tag, f"{gap:.3e}"))
    print(f"{name:38s}: {len(pairs)-len(fails):3d}/{len(pairs):3d} PASS | min_gap={mg:+.3e} at {mat} | fails={fails[:4]}")

test_rel("g_{k,i} << g_{k,i-1}", [((k,i), Rg[(k,i)], Rg[(k,i-1)]) for k in range(1,KMAX+1) for i in range(1,k+1)])
test_rel("S_{k,i} << S_{k,i-1}", [((k,i), RS[(k,i)], RS[(k,i-1)]) for k in range(2,KMAX+1) for i in range(1,k+1)])
test_rel("D_k g_{k-1,i} << D_k g_{k-1,i-1} (A<<C)", [((k,i), RD[(k,i)], RD[(k,i-1)]) for k in range(2,KMAX+1) for i in range(1,k)])
test_rel("D_k g_{k-1,i} << g_{k-1,i-2} (A<<D)", [((k,i), RD[(k,i)], Rg[(k-1,i-2)]) for k in range(3,KMAX+1) for i in range(2,k)])
test_rel("g_{k-1,i-1} << D_k g_{k-1,i} (B<<A)", [((k,i), Rg[(k-1,i-1)], RD[(k,i)]) for k in range(2,KMAX+1) for i in range(1,k)])
test_rel("D_k g_{k-1,i} << g_{k-1,i-1} (A<<B)", [((k,i), RD[(k,i)], Rg[(k-1,i-1)]) for k in range(2,KMAX+1) for i in range(1,k)])
test_rel("g_{k-1,i-1} << S_{k,i-1}", [((k,i), Rg[(k-1,i-1)], RS[(k,i-1)]) for k in range(2,KMAX+1) for i in range(1,k+1)])
test_rel("D_k g_{k-1,i} << S_{k,i-1}", [((k,i), RD[(k,i)], RS[(k,i-1)]) for k in range(2,KMAX+1) for i in range(1,k)])
test_rel("g_{k,i} << g_{k+1,i+1}", [((k,i), Rg[(k,i)], Rg[(k+1,i+1)]) for k in range(1,KMAX) for i in range(k+1)])
test_rel("D_{k+1} g_{k,i} << D_k g_{k,i}?", [((k,i), roots_mp(apply_Dtilde(k+1, g[(k,i)])), roots_mp(apply_Dtilde(k, g[(k,i)]))) for k in range(1,KMAX) for i in range(k+1)])
test_rel("D_k g_{k,i} << D_{k+1} g_{k,i}?", [((k,i), roots_mp(apply_Dtilde(k, g[(k,i)])), roots_mp(apply_Dtilde(k+1, g[(k,i)]))) for k in range(1,KMAX) for i in range(k+1)])
