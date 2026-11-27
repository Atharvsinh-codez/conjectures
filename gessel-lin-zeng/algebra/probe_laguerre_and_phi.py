"""
Deep probe of WHY prod_{j=1}^k (L_j + w_j y)(1) and prod_{j=1}^k (T_j + w_j y)(1)
satisfy interlacing q_{k,i} << q_{k,i-1} and Turan q_{k,i}^2 >= q_{k,i-1} q_{k,i+1} on the root domain!

Tests:
1. Non-PF_inf / non-classical differential operators in u:
   P_k(u, y) = prod_{j=1}^k (d/du + Phi_j(u) + w_j y)(1)
   Test Phi_j(u) = a_j - b_j u^3 (deg 3, W_j = exp(a u - b u^4/4), log-concave, NOT PF_inf!)
   Test Phi_j(u) = a_j - b_j (u + c u^3)
   Test Phi_j(u) = a_j - b_j (e^u + c e^{2u})
   Test Phi_j(u) = a_j - b_j (u^2 + 1) ... (not monotone, skip)
2. Generalized bidiagonal operators B_j(x^m) = alpha_{j,m} x^m + beta_{j,m} x^{m+1}
   - When does prod (B_j + w_j y)(1) have real roots + interlacing + Turan on (-inf, 0)?
   - Test arbitrary increasing alpha_{j,m} > 0, constant beta_{j,m} = b_j
   - Test alpha_{j,m} = m^2 + a_j (non-linear in m!)
   - Test alpha_{j,m} = c_j m + a_j (different leading slopes c_j!)
3. EXACT STRUCTURE OF LAGUERRE L_j = x d/dx + a_j + b_j x:
   - Is prod_{j=1}^k (L_j + w_j y)(1) real-stable in ANY transformed variables (e.g. (x, y/(-x)) or (x, -x y) or q-stable)?
   - Check B_k^L(x, -x y) and E_{k,i}^L = L^i q_{k-1, i} and q-stability for Laguerre!
"""
import math, random
import numpy as np
import mpmath as mp
mp.mp.dps = 40

def poly_add(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]

def poly_deriv(a):
    return [i * a[i] for i in range(1, len(a))] if len(a) > 1 else [0]

def poly_mul(a, b):
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            res[i+j] += ai * bj
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

# 1. Test Phi_j(u) = a_j - b_j u^3 (polynomial in u! So P_k(u,y) is a polynomial in u,y!)
# Wait: for W_j = exp(a_j u - b_j u^4/4), (d/du + Phi_j) f has degree deg(f) + 3!
# So deg q_{k,i}(u) = 3*(k-i)! Does it even stay real-rooted?! Let's check!
print("=== TEST A: Phi_j(u) = a_j - b_j u^3 (log-concave W, NOT PF_inf) ===")
def apply_Phi_poly(f, phi):
    # (d/du + phi(u)) f(u)
    return poly_add(poly_deriv(f), poly_mul(phi, f))

# Single operator (d/du - u^3)^n (1):
f = [1.0]
for n in range(1, 6):
    f = apply_Phi_poly(f, [0.0, 0.0, 0.0, -1.0])
    rts = np.roots(list(reversed(f)))
    nreal = sum(1 for r in rts if abs(r.imag) < 1e-6)
    print(f"  (d/du - u^3)^{n}(1): deg={len(f)-1}, real_roots={nreal}")

# What about T_j = c_j x d/dx + a_j + b_j x (DIFFERENT c_j > 0! Still PF_inf weights W_j(u) = exp(a_j u - b_j e^{u/c_j})?!)
# Wait: c_j x d/dx + a_j + b_j x = c_j (x d/dx + (a_j/c_j) + (b_j/c_j) x)!
# So dividing by c_j just scales w_j -> w_j/c_j! That's the SAME operator up to scalar!
# What about DIFFERENT exponents: B_j(x^m) = (m^2 + a_j) x^m + b_j x^{m+1}?
# Or DIFFERENT branch points: T_j = t(t + r_j) d/dt + a_j + b_j t? (Roots in (-r_j, 0)!)
print("\n=== TEST B: Generalized bidiagonal B_j(x^m) = (alpha(j,m) + w_j y) x^m + beta(j,m) x^{m+1} ===")
def test_bidiag(name, alpha_fn, beta_fn, trials=40, k=6):
    fail_rr = 0
    fail_int = 0
    fail_tur = 0
    for tr in range(trials):
        params = [(10**random.uniform(-1, 1), 10**random.uniform(-1, 1), 10**random.uniform(-1, 1)) for _ in range(k)]
        polys = {0: [1.0]}
        for j in range(1, k+1):
            a_j, b_j, w_j = params[j-1]
            new_p = {}
            for i in range(j+1):
                res = [0.0] * (j - i + 1)
                if i in polys:
                    p = polys[i]
                    for m, c in enumerate(p):
                        res[m] += alpha_fn(a_j, b_j, m) * c
                        res[m+1] += beta_fn(a_j, b_j, m) * c
                if (i-1) in polys:
                    p = polys[i-1]
                    for m, c in enumerate(p):
                        res[m] += w_j * c
                new_p[i] = res
            polys = new_p
        # Check real roots in (-inf, 0)
        roots = {}
        ok_rr = True
        for i in range(k+1):
            p = polys[i]
            if len(p) <= 1:
                roots[i] = []
                continue
            rts = np.roots(list(reversed(p)))
            if np.max(np.abs(np.imag(rts))) > 1e-6 or np.max(np.real(rts)) >= 0:
                ok_rr = False; break
            roots[i] = np.sort(np.real(rts))
        if not ok_rr:
            fail_rr += 1
            continue
        ok_int = True
        for i in range(1, k):
            r_lo, r_hi = roots[i], roots[i-1]
            for m in range(len(r_lo)):
                if not (r_hi[m] < r_lo[m] < r_hi[m+1]):
                    ok_int = False; break
        if not ok_int:
            fail_int += 1
        # Check Turan on (-inf, 0)
        xs = -np.logspace(-4, 3, 250)
        ok_tur = True
        for i in range(1, k):
            v0 = np.polyval(list(reversed(polys[i])), xs)
            vm = np.polyval(list(reversed(polys[i-1])), xs)
            vp = np.polyval(list(reversed(polys[i+1])), xs)
            if np.min(v0*v0 - vm*vp) < -1e-8 * np.max(np.abs(v0*v0)):
                ok_tur = False; break
        if not ok_tur:
            fail_tur += 1
    print(f"  {name}: RR_fail={fail_rr}/{trials}, Int_fail={fail_int}/{trials}, Turan_fail={fail_tur}/{trials}")

random.seed(123)
test_bidiag("Laguerre: alpha=m+a, beta=b", lambda a,b,m: m+a, lambda a,b,m: b)
test_bidiag("T_j: alpha=m+a, beta=m+a+b", lambda a,b,m: m+a, lambda a,b,m: m+a+b)
test_bidiag("Sub-T_j (b<a!): alpha=m+a, beta=m+b (random a,b)", lambda a,b,m: m+a, lambda a,b,m: m+b)
test_bidiag("Quad alpha: alpha=m^2+a, beta=b", lambda a,b,m: m*m+a, lambda a,b,m: b)
test_bidiag("Sqrt alpha: alpha=sqrt(m+a), beta=b", lambda a,b,m: math.sqrt(m+a), lambda a,b,m: b)
test_bidiag("Const alpha=a, beta=b (commuting!)", lambda a,b,m: a, lambda a,b,m: b)
test_bidiag("q-Laguerre (q=1.5): alpha=(1.5^m-1)+a, beta=b", lambda a,b,m: (1.5**m - 1) + a, lambda a,b,m: b)
test_bidiag("q-Laguerre (q=0.5): alpha=(1-0.5^m)+a, beta=b", lambda a,b,m: (1 - 0.5**m) + a, lambda a,b,m: b)
test_bidiag("2nd-order: alpha=(m+1)(m+a), beta=(m+1)(m+a+b)", lambda a,b,m: (m+1)*(m+a), lambda a,b,m: (m+1)*(m+a+b))
test_bidiag("Decr alpha: alpha=1/(m+a), beta=b", lambda a,b,m: 1.0/(m+a), lambda a,b,m: b)
test_bidiag("Pure beta=m+b, alpha=a", lambda a,b,m: a, lambda a,b,m: m+b)
test_bidiag("Negative slope beta=max(0.1, 10-m), alpha=m+a", lambda a,b,m: m+a, lambda a,b,m: max(0.1, 10-m))

# 3. Test DIFFERENT branch points: T_j = t(r_j + t) d/dt + a_j + b_j t (r_j > 0, roots of T_j(1) = -a_j/b_j in (-r_j, 0) iff a_j < r_j b_j)
# Notice T_j(t^m) = (m r_j + a_j) t^m + (m + b_j) t^{m+1}!
# Rescaling by 1/r_j: (m + a_j/r_j) t^m + ((m + b_j)/r_j) t^{m+1}!
test_bidiag("Diff branch r_j: alpha=m*r+a, beta=m+b (a<r*b)", lambda a,b,m: m*(a+1)+0.5*a*(a+1)*b, lambda a,b,m: m+b)

print("\n=== TEST C: Real-stability of Laguerre B_k^L(x, -x*y) and q-Laguerre ===")
# For Laguerre L_j = x d/dx + a_j + b_j x, check P_k^L(x, y) = prod (L_j + w_j y)(1)
# Check:
# (i) Is P_k^L(x, -x*y) real stable? (Wait: if y -> -x*y, at k=2 with a=b=w=1:
#     (1+3x+x^2) - 2x(1+x)y + x^2 y^2 = (1+x-xy)^2 + x ... at x=1: (2-y)^2 + 1 has complex roots!)
# (ii) What about P_k^L(x, y) with y -> y/x? (1+3x+x^2) + 2(1+x)(y/x) + y^2/x^2 -> x^2(1+3x+x^2) + 2x(1+x)y + y^2!
#     Discriminant in y: 4x^2(1+x)^2 - 4x^2(1+3x+x^2) = -4x^3 (still negative for x > 0!).
# (iii) Wait: what about Q_k^L(x, y) = sum q^{C(i,2)} q_{k,i}(x) y^i at q = 1/2?
for q_val in [mp.mpf('0.5'), mp.mpf('0.25')]:
    max_im = mp.mpf('-1e100')
    for trial in range(5):
        k = 6
        polys = {0: [mp.mpf(1)]}
        for j in range(1, k+1):
            a_j = mp.mpf(10 ** random.uniform(-1, 1))
            b_j = mp.mpf(10 ** random.uniform(-1, 1))
            w_j = mp.mpf(10 ** random.uniform(-1, 1))
            new_p = {}
            for i in range(j+1):
                res = [mp.mpf(0)] * (j - i + 1)
                if i in polys:
                    for m, c in enumerate(polys[i]):
                        res[m] += (m + a_j) * c
                        res[m+1] += b_j * c
                if (i-1) in polys:
                    for m, c in enumerate(polys[i-1]):
                        res[m] += w_j * c
                new_p[i] = res
            polys = new_p
        # Check UHP stability of sum q^{i(i-1)/2} polys[i](x) y^i
        for xval in [mp.mpc(-10, '1e-5'), mp.mpc(-1, '1e-5'), mp.mpc(-0.1, '1e-5'), mp.mpc(0.1, '1e-5'), mp.mpc(1, '1e-5'), mp.mpc(10, '1e-5'), mp.mpc(0, 1)]:
            cy = []
            for i in range(k+1):
                val = mp.mpc(0)
                for c in reversed(polys[i]):
                    val = val * xval + c
                cy.append(val * (q_val ** (i*(i-1)//2)))
            rts = mp.polyroots(list(reversed(cy)), maxsteps=200, extraprec=20)
            for r in rts:
                if r.imag > max_im:
                    max_im = r.imag
    print(f"  q-Laguerre (q={q_val}, 5 wild trials k=6): UHP x -> max Im(y) = {mp.nstr(max_im, 6)}")
