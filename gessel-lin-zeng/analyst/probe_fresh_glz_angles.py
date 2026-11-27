"""
Probe fresh GLZ-specific structural angles:
1. Check whether A_{k,i} interlaces ANY boundary/known polynomial A_{m,0}, H_m, A_{m,m}
   for ALL (k,i) up to k=10 (especially deep interior 4 <= i <= k-3!).
2. Check exact scaling of s_{k,i} = htilde'_{k,i}(0)/htilde_{k,i}(0), rightmost root r_max,
   and leftmost root r_min of htilde_{k,i}.
3. Test whether Two-Regime condition (J_raw >= 1 on (-1, -1/m], A>0 on (-1,0)) is preserved
   by D_m on general triples (or breaks even when b=m is small)!
4. Check why Nelder-Mead broke G_j at j=80 in probe_jacobi_and_Gj_opt.py: does G_j break
   at GLZ-realistic degrees/scales?
"""
import numpy as np
import sympy as sp

def compute_A_ht(K):
    t, x = sp.symbols('t x')
    def D(m, f): return sp.expand(t*(1+t)*sp.diff(f, t) + (1 + m*t)*f)
    ht = {(1,0): 1+2*t, (1,1): sp.Integer(1)}
    for k in range(2, K+1):
        for i in range(k+1):
            ht[(k,i)] = D(k+1, (D(k+1, ht[(k-1,i)]) if i <= k-1 else 0) + (ht[(k-1,i-1)] if i >= 1 else 0))
    # Convert to numpy root arrays in t in (-1, 0)
    roots = {}
    for (k,i), p in ht.items():
        coeffs = [float(c) for c in sp.Poly(p, t).all_coeffs()] # high to low
        if len(coeffs) <= 1:
            roots[(k,i)] = np.array([])
        else:
            r = np.roots(coeffs)
            assert np.max(np.abs(r.imag)) < 1e-8
            roots[(k,i)] = np.sort(r.real) # ascending in (-1, 0)
    return ht, roots, t

def interlaces(r_small, r_big):
    # Check if r_small << r_big or r_big << r_small (deg diff 0 or 1)
    n1, n2 = len(r_small), len(r_big)
    if n2 == n1 + 1:
        # r_big[0] < r_small[0] < r_big[1] < ... < r_big[n1]
        ok = all(r_big[j] < r_small[j] < r_big[j+1] for j in range(n1))
        return "strict_inside" if ok else None
    if n1 == n2:
        if all(r_small[j] < r_big[j] for j in range(n1)) and all(r_big[j] < r_small[j+1] for j in range(n1-1)):
            return "left_of"
        if all(r_big[j] < r_small[j] for j in range(n1)) and all(r_small[j] < r_big[j+1] for j in range(n1-1)):
            return "right_of"
    return None

ht, roots, t = compute_A_ht(10)

print("=== 1. Check interlacing of htilde_{k,i} against boundary families htilde_{m,0} and htilde_{m,m} ===")
for k in range(5, 11):
    for i in range(2, k-1):
        d = 2*k - i - 1
        r = roots[(k, i)]
        # Same degree m_diag = d + 1 has htilde_{m_diag, m_diag} of degree d
        hits = []
        for m in range(1, 11):
            for typ in [0, m]:
                if (m, typ) in roots:
                    rb = roots[(m, typ)]
                    if abs(len(rb) - len(r)) <= 1:
                        rel = interlaces(r, rb) if len(r) <= len(rb) else interlaces(rb, r)
                        if rel:
                            hits.append(f"({m},{typ}):{rel}")
        if k >= 7 and 3 <= i <= k-2:
            print(f"k={k:2d}, i={i} (deg {d:2d}): boundary interlacers = {hits}")

print("\n=== 2. Check s_{k,i}, rightmost root r_max, and -1/(k+2) ===")
for k in [3, 5, 7, 10]:
    print(f"k={k}, -1/(k+2) = {-1/(k+2):.6f}:")
    for i in range(k+1):
        p = sp.Poly(ht[(k,i)], t)
        c0 = float(p.nth(0))
        c1 = float(p.nth(1))
        s = c1 / c0
        r = roots[(k,i)]
        rmax = r[-1] if len(r) else -1.0
        rmin = r[0] if len(r) else -1.0
        n_regII = int(np.sum(r > -1/(k+2)))
        print(f"  i={i:2d} (deg {len(r):2d}): s={s:10.3f}, -1/s={-1/s:10.6f}, r_max={rmax:10.6f}, #roots in (-1/(k+2),0) = {n_regII}")

print("\n=== 3. Check 3-term or differential identities in i at fixed k (k=4,5) ===")
# Check if htilde_{k,i} satisfies any low-degree polynomial relation in i, or if D_m preserves J_raw>=1 on (-1,-1/m]
# Test D_m on J_raw >= 1 on (-1, -1/m] for random triples!
rng = np.random.default_rng(42)
def eval_Jraw_on_grid(g, f, h, grid):
    # g,f,h numpy poly1d (high to low)
    gp, fp, hp = np.polyder(g), np.polyder(f), np.polyder(h)
    gv, fv, hv = np.polyval(g, grid), np.polyval(f, grid), np.polyval(h, grid)
    gd, fd, hd = np.polyval(gp, grid), np.polyval(fp, grid), np.polyval(hp, grid)
    Wgf = gv*fd - gd*fv
    Wfh = fv*hd - fd*hv
    Wgh = gv*hd - gd*hv
    with np.errstate(divide='ignore', invalid='ignore'):
        J = 4*Wgf*Wfh / (Wgh**2)
    return np.nanmin(J)

def apply_Dm_np(m, p):
    # p is high-to-low coeffs
    c = p[::-1] # low to high
    out = np.zeros(len(c) + 1)
    for j, cj in enumerate(c):
        out[j] += (j + 1) * cj
        out[j + 1] += (j + m) * cj
    return out[::-1]

fails_Jraw_Dm = 0
tested_Jraw = 0
min_Jraw_out = 1e9
for trial in range(400):
    m = rng.integers(4, 25)
    d = rng.integers(3, 8)
    # Generate strictly interlacing roots in (-1, 0): g deg d+1, f deg d, h deg d-1
    pts = -np.sort(rng.uniform(0.01, 0.99, 2*d))[::-1] # wait, we want ascending in (-1,0)
    pts = np.sort(-rng.uniform(0.01, 0.99, 2*d))
    # wait need d+1 + d + d-1 = 3d roots! Let's make g(d+1) > f(d) > h(d-1) interlacing
    r_all = np.sort(-rng.uniform(0.005, 0.995, 3*d))
    # g takes 0, 3, 6, ... wait: g has d+1 roots, f has d, h has d-1
    # Standard interlacing: g_0 < f_0 < g_1 < ... and f_0 < h_0 < f_1 < ...
    rg = np.sort(-rng.uniform(0.005, 0.995, d+1))
    rf = np.array([rng.uniform(rg[j], rg[j+1]) for j in range(d)])
    rh = np.array([rng.uniform(rf[j], rf[j+1]) for j in range(d-1)])
    g = np.poly(rg)
    f = np.poly(rf)
    h = np.poly(rh)
    # Normalize g(0)=f(0)=h(0)=1 (positive at 0)
    g = g / np.polyval(g, 0)
    f = f / np.polyval(f, 0)
    h = h / np.polyval(h, 0)
    # Boost f slightly? Actually J_raw is invariant under scaling (g,f,h) -> (cg, cf, ch) EXCEPT f^2/(gh)!
    # Wait: Wgf*Wfh / Wgh^2 scales as f^2/(gh)! So scale f to match binom ratio!
    scale_f = rng.uniform(1.0, 1.5)
    f = f * scale_f
    grid_in = np.linspace(-0.999, -1.0/m, 200)
    jin = eval_Jraw_on_grid(g, f, h, grid_in)
    if jin >= 1.0:
        tested_Jraw += 1
        gm, fm, hm = apply_Dm_np(m, g), apply_Dm_np(m, f), apply_Dm_np(m, h)
        grid_out = np.linspace(-0.999, -1.0/(m+1), 200)
        jout = eval_Jraw_on_grid(gm, fm, hm, grid_out)
        if jout < min_Jraw_out: min_Jraw_out = jout
        if jout < 1.0:
            fails_Jraw_Dm += 1
print(f"Random triples with J_raw>=1 on (-1,-1/m]: tested={tested_Jraw}, D_m drops J_raw<1 on (-1,-1/(m+1)]: {fails_Jraw_Dm}, min_Jout={min_Jraw_out:.5f}")
