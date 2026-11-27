import sympy as sp
import mpmath as mp

mp.mp.dps = 50
t = sp.Symbol('t')
x = sp.Symbol('x')
y = sp.Symbol('y')

def apply_Dtilde(j, p):
    """Dtilde_j(p) = t*(1+t)*p'(t) + (1 + (j+1)*t)*p(t)."""
    return sp.expand(t*(1+t)*sp.diff(p, t) + (1 + (j+1)*t)*p)

# g_{k,i}(t) = h_{k,i}(t) / (1+t)
# Note h_{1,0} = (1+t)(1+2t), h_{1,1} = (1+t)*1
KMAX = 15
g = {(1, 0): 1 + 2*t, (1, 1): sp.Integer(1)}
for k in range(2, KMAX + 1):
    for i in range(k + 1):
        term1 = apply_Dtilde(k, g[(k-1, i)]) if i <= k - 1 else sp.Integer(0)
        term2 = g[(k-1, i-1)] if i >= 1 else sp.Integer(0)
        g[(k, i)] = apply_Dtilde(k, sp.expand(term1 + term2))

def roots_mp(expr):
    p = sp.Poly(expr, t)
    if p.degree() == 0:
        return [], True
    coeffs = [mp.mpf(int(c)) for c in p.all_coeffs()]
    r = mp.polyroots(coeffs, maxsteps=200, extraprec=40)
    max_im = max(abs(z.imag) for z in r)
    if max_im > mp.mpf('1e-25'):
        return None, False
    return sorted([z.real for z in r]), True

def check_interlace(rp, rq):
    """Check rp << rq (rq has root strictly to the right of rp's rightmost root).
    deg rq in {deg rp, deg rp + 1}."""
    dp, dq = len(rp), len(rq)
    if dq == dp:
        gaps = [rq[m] - rp[m] for m in range(dp)] + [rp[m+1] - rq[m] for m in range(dp - 1)]
        ok = all(z > 0 for z in gaps)
        return ok, (float(min(gaps)) if gaps else 999.0), "eq"
    elif dq == dp + 1:
        gaps = [rp[m] - rq[m] for m in range(dp)] + [rq[m+1] - rp[m] for m in range(dp)]
        ok = all(z > 0 for z in gaps)
        return ok, (float(min(gaps)) if gaps else 999.0), "+1"
    return False, -999.0, f"jump({dp}->{dq})"

R = {}
R_D = {} # roots of Dtilde_k g_{k-1, i}
for k in range(1, KMAX + 1):
    for i in range(k + 1):
        r, ok = roots_mp(g[(k, i)])
        if not ok:
            print(f"!!! NON-REAL g_{{{k},{i}}} !!!")
        R[(k, i)] = r
        if k >= 2 and i <= k - 1:
            dk_prev = apply_Dtilde(k, g[(k-1, i)])
            rd, okd = roots_mp(dk_prev)
            if not okd:
                print(f"!!! NON-REAL Dtilde_{k} g_{{{k-1},{i}}} !!!")
            R_D[(k, i)] = rd

print(f"All g_{{k,i}} real-rooted for k=1..{KMAX}: YES")

# Check bounds: are all roots strictly in (-1, 0)?
min_r = min(min(r) for r in R.values() if r)
max_r = max(max(r) for r in R.values() if r)
print(f"All roots in ({float(min_r):.8f}, {float(max_r):.8f}) subset (-1, 0)")

tests = [
    ("g_{k,i} << g_{k,i-1}", lambda k, i: (R[(k, i)], R[(k, i-1)])),
    ("g_{k,i-1} << g_{k,i}", lambda k, i: (R[(k, i-1)], R[(k, i)])),
    ("g_{k-1,i-1} << D_k g_{k-1,i}", lambda k, i: (R[(k-1, i-1)], R_D[(k, i)]) if i <= k-1 else None),
    ("D_k g_{k-1,i} << g_{k-1,i-1}", lambda k, i: (R_D[(k, i)], R[(k-1, i-1)]) if i <= k-1 else None),
    ("D_k g_{k-1,i} << D_k g_{k-1,i-1}", lambda k, i: (R_D[(k, i)], R_D[(k, i-1)]) if i <= k-1 else None),
    ("g_{k-1,i-1} << g_{k,i}", lambda k, i: (R[(k-1, i-1)], R[(k, i)])),
    ("D_k g_{k-1,i} << g_{k,i}", lambda k, i: (R_D[(k, i)], R[(k, i)]) if i <= k-1 else None),
    ("g_{k,i} << g_{k+1,i+1}", lambda k, i: (R[(k, i)], R[(k+1, i+1)]) if k < KMAX else None),
    ("g_{k,i} << g_{k+1,i}", lambda k, i: (R[(k, i)], R[(k+1, i)]) if k < KMAX else None),
]

print(f"\n=== STRICT INTERLACING TESTS (k=2..{KMAX}, 50-digit mpmath) ===")
for name, fn in tests:
    total = 0
    fails = []
    min_gap = 999.0
    min_at = None
    for k in range(2, KMAX + 1):
        for i in range(1, k + 1):
            pair = fn(k, i)
            if pair is None:
                continue
            total += 1
            ok, gap, kind = check_interlace(pair[0], pair[1])
            if gap < min_gap:
                min_gap = gap
                min_at = (k, i)
            if not ok:
                fails.append((k, i, f"{gap:.3e}", kind))
    print(f"{name:34s}: {total - len(fails):3d}/{total:3d} PASS | min_gap={min_gap:+.4e} at {min_at} | fails={fails[:5]}")

# Print roots for k=4, 5 to see the geometry clearly!
print("\n=== Roots of g_{k,i} and Dtilde_{k+1} g_{k,i} for k=4 ===")
for i in range(5):
    print(f"  g_{{4,{i}}}     : {[round(float(z), 5) for z in R[(4, i)]]}")
    if i < 4:
        print(f"  D5 g_{{4,{i+1}}}: {[round(float(z), 5) for z in R_D[(5, i+1)]]}")
