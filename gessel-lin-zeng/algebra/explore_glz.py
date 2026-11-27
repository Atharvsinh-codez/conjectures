import sympy as sp
import numpy as np

t = sp.Symbol('t')
x = sp.Symbol('x')
y = sp.Symbol('y')

def apply_D(j, p):
    """D_j(p) = t*(1+t)*p'(t) + (1 + j*t)*p(t)."""
    return sp.expand(t*(1+t)*sp.diff(p, t) + (1 + j*t)*p)

# Compute h_{k,i}(t) for k = 0..16
KMAX = 16
h = {(0, 0): sp.Integer(1)}
for k in range(1, KMAX + 1):
    for i in range(k + 1):
        term1 = apply_D(k, h[(k-1, i)]) if i <= k - 1 else sp.Integer(0)
        term2 = h[(k-1, i-1)] if i >= 1 else sp.Integer(0)
        h[(k, i)] = apply_D(k, sp.expand(term1 + term2))

print("=== Small values of h_{k,i}(t) and A_{k,i}(x) ===")
for k in range(1, 5):
    for i in range(k + 1):
        # h_{k,i}(t) = (1+t)^{2k-i} A_{k,i}(t/(1+t))
        # So A_{k,i}(x) = (1-x)^{2k-i} h_{k,i}(x/(1-x))
        Aki = sp.simplify((1 - x)**(2*k - i) * h[(k, i)].subs(t, x / (1 - x)))
        print(f"k={k}, i={i}: h = {sp.factor(h[(k,i)])} | A = {sp.expand(Aki)}")

def roots_in_open_minus1_0(expr):
    """Return sorted real roots (excluding -1), multiplicity of -1, and whether all roots are real."""
    p = sp.Poly(expr, t)
    # Factor out (1+t)
    m_neg1 = 0
    cur = p
    while cur.eval(-1) == 0:
        cur = sp.Poly(sp.quo(cur, sp.Poly(t + 1, t)), t)
        m_neg1 += 1
    if cur.degree() == 0:
        return [], m_neg1, True
    coeffs = [float(c) for c in cur.all_coeffs()]
    r = np.roots(coeffs)
    if np.max(np.abs(r.imag)) > 1e-6:
        return None, m_neg1, False
    return list(np.sort(r.real)), m_neg1, True

def check_interlace_roots(rp, rq, strict=False):
    """Check if rp << rq (rq has root to the right of rp's largest root, deg rq in {deg rp, deg rp + 1}).
    Also returns min gap."""
    dp, dq = len(rp), len(rq)
    eps = -1e-9 if strict else 1e-9
    if dq == dp:
        # rp[0] <= rq[0] <= rp[1] <= rq[1] <= ... <= rq[dp-1]
        gaps = [rq[m] - rp[m] for m in range(dp)] + [rp[m+1] - rq[m] for m in range(dp - 1)]
        ok = all(g >= -eps for g in gaps) if not strict else all(g > 0 for g in gaps)
        return ok, (min(gaps) if gaps else 999.0), "eq"
    elif dq == dp + 1:
        # rq[0] <= rp[0] <= rq[1] <= rp[1] <= ... <= rp[dp-1] <= rq[dp]
        gaps = [rp[m] - rq[m] for m in range(dp)] + [rq[m+1] - rp[m] for m in range(dp)]
        ok = all(g >= -eps for g in gaps) if not strict else all(g > 0 for g in gaps)
        return ok, (min(gaps) if gaps else 999.0), "+1"
    else:
        return False, -999.0, f"jump({dp}->{dq})"

print("\n=== Check roots of h_{k,i}(t) and ALL interlacing relations for k=1..14 ===")
R = {}
M1 = {}
for k in range(1, 15):
    for i in range(k + 1):
        r, m1, is_real = roots_in_open_minus1_0(h[(k, i)])
        if not is_real:
            print(f"NON-REAL h_{{{k},{i}}}!")
        R[(k, i)] = r
        M1[(k, i)] = m1

# Check degree and multiplicity at -1 pattern
print("Multiplicity of (t+1) in h_{k,i} for k=1..6:")
for k in range(1, 7):
    print(f"  k={k}: " + ", ".join(f"i={i}:deg={len(R[(k,i)])}+m(-1)={M1[(k,i)]}" for i in range(k+1)))

# Test key interlacing hypotheses!
hyps = {
    "h_{k,i} << h_{k,i-1}": [],
    "h_{k,i-1} << h_{k,i}": [],
    "h_{k-1,i-1} << D_k h_{k-1,i}": [],
    "D_k h_{k-1,i} << h_{k-1,i-1}": [],
    "h_{k-1,i} << h_{k,i}": [],
    "h_{k-1,i-1} << h_{k,i}": [],
    "h_{k,i} << D_{k+1} h_{k,i+1}? no, h_{k,i} << D_{k+1} h_{k,i}": [],
    "h_{k-1,i} << D_k h_{k-1,i}": [],
    "D_k h_{k-1,i} << D_k h_{k-1,i-1}": [],
    "h_{k,i} << h_{k+1,i+1}": [],
}

print("\n=== Checking pairwise interlacing across k=1..14 ===")
for k in range(2, 15):
    for i in range(1, k + 1):
        ok1, g1, _ = check_interlace_roots(R[(k, i)], R[(k, i-1)])
        ok2, g2, _ = check_interlace_roots(R[(k, i-1)], R[(k, i)])
        hyps["h_{k,i} << h_{k,i-1}"].append((k, i, ok1, g1))
        hyps["h_{k,i-1} << h_{k,i}"].append((k, i, ok2, g2))
        if i <= k - 1:
            dk_prev = apply_D(k, h[(k-1, i)])
            r_dk, m_dk, real_dk = roots_in_open_minus1_0(dk_prev)
            ok3, g3, _ = check_interlace_roots(R[(k-1, i-1)], r_dk)
            ok4, g4, _ = check_interlace_roots(r_dk, R[(k-1, i-1)])
            hyps["h_{k-1,i-1} << D_k h_{k-1,i}"].append((k, i, ok3, g3))
            hyps["D_k h_{k-1,i} << h_{k-1,i-1}"].append((k, i, ok4, g4))
            ok5, g5, _ = check_interlace_roots(R[(k-1, i)], R[(k, i)])
            hyps["h_{k-1,i} << h_{k,i}"].append((k, i, ok5, g5))
            ok6, g6, _ = check_interlace_roots(R[(k-1, i-1)], R[(k, i)])
            hyps["h_{k-1,i-1} << h_{k,i}"].append((k, i, ok6, g6))
            # Check D_k h_{k-1,i} vs D_k h_{k-1,i-1}
            dk_prev2 = apply_D(k, h[(k-1, i-1)])
            r_dk2, _, _ = roots_in_open_minus1_0(dk_prev2)
            ok9, g9, _ = check_interlace_roots(r_dk, r_dk2)
            hyps["D_k h_{k-1,i} << D_k h_{k-1,i-1}"].append((k, i, ok9, g9))

for name, res in hyps.items():
    if not res:
        continue
    fails = [(k, i, round(g, 6)) for k, i, ok, g in res if not ok]
    min_g = min(g for _, _, _, g in res)
    print(f"{name:34s}: total={len(res):3d}, fails={len(fails):3d}, min_gap={min_g:+.6e} | first fails={fails[:5]}")

print("\n=== Print roots of h_{k,i} for k=5, 6 ===")
for k in [5, 6]:
    print(f"--- k={k} ---")
    for i in range(k + 1):
        print(f"  i={i} (m_-1={M1[(k,i)]}): {[round(float(z), 5) for z in R[(k,i)]]}")
