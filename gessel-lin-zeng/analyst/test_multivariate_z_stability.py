"""EUREKA TEST:
For Laguerre: replace m_j * t by independent variable z_j!
Then L_{c_j, m_j} = t d/dt + c_j + m_j t  becomes  E + c_j + z_j  where E = sum z_s d/dz_s (degree operator!).
Define the MULTI-AFFINE polynomial:
  F_{N, i}(z_1, ..., z_N) := [y^i] prod_{j=1}^N (E + c_j + z_j + lam_j y)(1).
Note that for any subset S subseteq [N], the coefficient of prod_{j in S} z_j is POSITIVE!
Test:
1. Is F_{N, i}(z_1, ..., z_N) REAL-STABLE in (z_1, ..., z_N) for ALL N = 2..7, ALL i, ALL c_j > 0, lam_j >= 0?!
   (Test BOTH UHP zeros AND Rayleigh Delta_{a,b} = F_a F_b - F F_{ab} >= 0 on R^N!)
2. What about the BIVARIATE/MULTIVARIATE operator E + c_j + z_j + lam_j y?
   (We know at z_1=...=z_N it's NOT stable in (t, y), so we test EACH [y^i] coefficient F_{N, i}(z_1,...,z_N)!)
3. NOW TEST GLZ / ATOMIC D_{c_j, m_j}:
   D_{c_j, m_j} = (1 + t)*E + c_j + z_j where t = z_0 (so deg operator E = z_0 d/dz_0 + sum z_s d/dz_s)!
   Is G_{N, i}(z_0, z_1, ..., z_N) := [y^i] prod_{j=1}^N ((1+z_0)E + c_j + z_j + lam_j y)(1)
   REAL-STABLE in (z_0, z_1, ..., z_N) — OR with z_0 = t, z_j = (m_j - c_j) t, i.e.
   D_{c_j, m_j} = t(1+t) d/dt + c_j(1+t) + (m_j - c_j) t?!
   WAIT A SECOND!! Look at that split:
   D_{c_j, m_j} = (1 + t)(t d/dt + c_j) + (m_j - c_j) t!!
   Let u_j = (m_j - c_j) t! Test BOTH splits for multivariate real stability!!
"""
import numpy as np
import itertools

def build_laguerre_multiaffine(cs, lams):
    """Build F_{N, i}(z_1,...,z_N) as dicts {bitmask: coeff} for i=0..N."""
    N = len(cs)
    # polys[i] maps bitmask (subset of 0..N-1) -> float coeff
    polys = [{0: 1.0}]
    for j in range(N):
        cj = cs[j]; lj = lams[j]; bit = 1 << j
        new_p = [{} for _ in range(len(polys) + 1)]
        for i, d in enumerate(polys):
            for mask, coef in d.items():
                deg = mask.bit_count()
                # (E + c_j + z_j) term:
                v0 = coef * (deg + cj)
                new_p[i][mask] = new_p[i].get(mask, 0.0) + v0
                new_p[i][mask | bit] = new_p[i].get(mask | bit, 0.0) + coef
                # lam_j * y term:
                if lj != 0.0:
                    new_p[i + 1][mask] = new_p[i + 1].get(mask, 0.0) + lj * coef
        polys = new_p
    return polys

def eval_multiaffine(d, zs):
    val = 0.0 + 0.0j
    N = len(zs)
    for mask, coef in d.items():
        term = coef
        m = mask; idx = 0
        while m:
            if m & 1: term *= zs[idx]
            m >>= 1; idx += 1
        val += term
    return val

def check_multiaffine_uhp(d, N, rng, trials=120):
    """Solve linear in z_k for random z_{-k} in H+^ (N-1); real-stable iff Im(z_k) <= 0 ALWAYS!"""
    max_im = -1e9
    for _ in range(trials):
        k = int(rng.integers(0, N))
        bit = 1 << k
        zs = [complex(rng.normal(0, 2) * (10**rng.uniform(-1, 1)), 10**rng.uniform(-3, 1)) for _ in range(N)]
        # f = A + B * z_k
        A = 0.0 + 0.0j; B = 0.0 + 0.0j
        for mask, coef in d.items():
            has_k = bool(mask & bit)
            m = mask & (~bit)
            term = coef
            idx = 0
            while m:
                if m & 1: term *= zs[idx]
                m >>= 1; idx += 1
            if has_k: B += term
            else: A += term
        if abs(B) < 1e-20: continue
        zk = -A / B
        if zk.imag > max_im:
            max_im = zk.imag
    return max_im

def check_multiaffine_rayleigh(d, N, rng, trials=100):
    """Check F_a F_b - F F_{ab} >= 0 for random x in R^N and all pairs a < b."""
    min_rel = 1e9
    for _ in range(trials):
        xs = [rng.normal(0, 2) * (10**rng.uniform(-1, 1)) for _ in range(N)]
        for a in range(N):
            for b in range(a + 1, N):
                ba = 1 << a; bb = 1 << b
                # f = c00 + c10*xa + c01*xb + c11*xa*xb
                c00 = c10 = c01 = c11 = 0.0
                for mask, coef in d.items():
                    ha = bool(mask & ba); hb = bool(mask & bb)
                    m = mask & (~(ba | bb))
                    term = coef; idx = 0
                    while m:
                        if m & 1: term *= xs[idx]
                        m >>= 1; idx += 1
                    if ha and hb: c11 += term
                    elif ha: c10 += term
                    elif hb: c01 += term
                    else: c00 += term
                # Delta_{a,b} = c10*c01 - c00*c11 (independent of xa, xb!!)
                det = c10 * c01 - c00 * c11
                scale = abs(c10 * c01) + abs(c00 * c11) + 1e-300
                rel = det / scale
                if rel < min_rel:
                    min_rel = rel
    return min_rel

print("=== 1. MULTIVARIATE MULTI-AFFINE LAGUERRE F_{N, i}(z_1, ..., z_N) REAL STABILITY ===")
rng = np.random.default_rng(2026)
for N in [3, 4, 5, 6]:
    # Test c_j = 1, lam_j = 1
    polys = build_laguerre_multiaffine([1.0]*N, [1.0]*N)
    worst_im = -1e9; worst_ray = 1e9
    for i in range(N):
        mi = check_multiaffine_uhp(polys[i], N, rng, trials=150)
        mr = check_multiaffine_rayleigh(polys[i], N, rng, trials=80)
        worst_im = max(worst_im, mi); worst_ray = min(worst_ray, mr)
    print(f"  N={N} (c=1, lam=1): max Im(z_k) = {worst_im:+.3e}, min rel Rayleigh = {worst_ray:+.3e}")

    # Test WILD c_j, lam_j (10 random trials)
    wild_uhp_fails = 0; wild_ray_fails = 0
    for trial in range(10):
        cs = [float(10**rng.uniform(-1, 1.5)) for _ in range(N)]
        lams = [0.0 if rng.random() < 0.2 else float(10**rng.uniform(-1, 1)) for _ in range(N)]
        polys = build_laguerre_multiaffine(cs, lams)
        for i in range(N):
            if not polys[i]: continue
            mi = check_multiaffine_uhp(polys[i], N, rng, trials=80)
            mr = check_multiaffine_rayleigh(polys[i], N, rng, trials=40)
            if mi > 1e-7: wild_uhp_fails += 1; break
            if mr < -1e-7: wild_ray_fails += 1; break
    print(f"  N={N} (WILD c, lam, 10 trials): uhp_fails = {wild_uhp_fails}/10, ray_fails = {wild_ray_fails}/10", flush=True)

print("\n=== 2. MULTIVARIATE GLZ / ATOMIC D_{c_j, m_j} WITH z_0 = t AND z_j = (m_j - c_j)*t ===")
# Notice: D_{c_j, m_j} = t(1+t) d/dt + c_j + m_j t = (1 + t)(t d/dt + c_j) + (m_j - c_j) t!
# In terms of z_0 = t and z_j = (m_j - c_j) t (multi-affine in z_1..z_N, polynomial in z_0):
#   E = z_0 d/dz_0 + sum z_s d/dz_s (total degree in all z's!)
#   Operator_j = (1 + z_0) * (E + c_j) + z_j  !
# Check at z_j = (m_j - c_j) z_0:
#   For a homogeneous deg-d term in z_0, z_1..z_N, E gives d, so (1+z_0)(d + c_j) + (m_j - c_j) z_0
#   = (d + c_j) + (d + m_j) z_0 = D_{c_j, m_j}!! EXACT!!
# ALSO test the split: D_{c_j, m_j} = (E + c_j) + z_0*E + z_j  (where z_j = m_j t, WITHOUT subtracting c_j)!
for split_name, subtract_c in [("Split A: (1+z0)(E+c_j) + z_j", True), ("Split B: (E+c_j) + z0*E + z_j", False)]:
    print(f"--- {split_name} ---")
    for N in [3, 4, 5]:
        # Represent poly as dict {(deg_z0, bitmask_z): coeff}
        cs = [1.0]*N
        lams = [1.0 if j % 2 == 0 else 0.0 for j in range(N)] if N == 4 else [1.0]*N
        polys = [{(0, 0): 1.0}]
        for j in range(N):
            cj = cs[j]; lj = lams[j]; bit = 1 << j
            new_p = [{} for _ in range(len(polys) + 1)]
            for i, d in enumerate(polys):
                for (d0, mask), coef in d.items():
                    tot_deg = d0 + mask.bit_count()
                    if subtract_c:
                        # (1 + z0)*(tot_deg + cj) + z_j
                        v = coef * (tot_deg + cj)
                        new_p[i][(d0, mask)] = new_p[i].get((d0, mask), 0.0) + v
                        new_p[i][(d0 + 1, mask)] = new_p[i].get((d0 + 1, mask), 0.0) + v
                        new_p[i][(d0, mask | bit)] = new_p[i].get((d0, mask | bit), 0.0) + coef
                    else:
                        # (tot_deg + cj) + z0*tot_deg + z_j
                        new_p[i][(d0, mask)] = new_p[i].get((d0, mask), 0.0) + coef * (tot_deg + cj)
                        if tot_deg > 0:
                            new_p[i][(d0 + 1, mask)] = new_p[i].get((d0 + 1, mask), 0.0) + coef * tot_deg
                        new_p[i][(d0, mask | bit)] = new_p[i].get((d0, mask | bit), 0.0) + coef
                    if lj != 0.0:
                        new_p[i + 1][(d0, mask)] = new_p[i + 1].get((d0, mask), 0.0) + lj * coef
            polys = new_p
        # Test UHP stability in (z0, z1, ..., zN): fix z0, z_{-k} in H+, solve linear for z_k!
        # AND fix z1..zN in H+, solve polynomial for z0!
        max_im_zk = -1e9; max_im_z0 = -1e9
        for i, d in enumerate(polys):
            if not d: continue
            # Solve for z_k
            for _ in range(80):
                z0 = complex(rng.normal(0, 1), 10**rng.uniform(-2, 1))
                zs = [complex(rng.normal(0, 1), 10**rng.uniform(-2, 1)) for _ in range(N)]
                k = int(rng.integers(0, N)); bit = 1 << k
                A = 0j; B = 0j
                for (d0, mask), coef in d.items():
                    has_k = bool(mask & bit)
                    m = mask & (~bit)
                    term = coef * (z0 ** d0)
                    idx = 0
                    while m:
                        if m & 1: term *= zs[idx]
                        m >>= 1; idx += 1
                    if has_k: B += term
                    else: A += term
                if abs(B) > 1e-20:
                    max_im_zk = max(max_im_zk, (-A / B).imag)
            # Solve for z0
            max_d0 = max(d0 for (d0, _) in d.keys())
            if max_d0 >= 1:
                for _ in range(60):
                    zs = [complex(rng.normal(0, 1), 10**rng.uniform(-2, 1)) for _ in range(N)]
                    c0 = np.zeros(max_d0 + 1, dtype=np.complex128)
                    for (d0, mask), coef in d.items():
                        term = coef; m = mask; idx = 0
                        while m:
                            if m & 1: term *= zs[idx]
                            m >>= 1; idx += 1
                        c0[d0] += term
                    if abs(c0[-1]) > 1e-20:
                        rts = np.roots(c0[::-1])
                        max_im_z0 = max(max_im_z0, float(np.max(rts.imag)))
        print(f"  N={N}: max Im(z_k) = {max_im_zk:+.3e}, max Im(z_0) = {max_im_z0:+.3e}", flush=True)
