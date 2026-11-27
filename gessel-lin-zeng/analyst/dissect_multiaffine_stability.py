"""Investigate:
1. Why did Laguerre F_{N,i}(z_1,...,z_N) pass 100% for c=1, lam=1 (N=3..6) while failing 1/10 for WILD c, lam?
   Test separately:
   (a) c_j = 1 (constant), WILD lam_j >= 0 (including GLZ pattern [1,0,1,0,...]!)
   (b) WILD c_j > 0, lam_j = 1 (constant)!
   (c) Equal c_j = c > 0, lam_j = 1!
   (d) Print the EXACT counterexample when it fails!
2. Test FULL MULTI-AFFINE GLZ:
   (a) Operator_j = (1 + w_j)*(E + c_j) + z_j  (in 2N variables w_1..w_N, z_1..z_N!)
   (b) Operator_j = (E + c_j) + w_j * E + z_j
   (c) What about the x-space operator S_L = 1 + L*x + x(1-x) d/dx = (1 + E) + (L*x - x*E) -> wait: has MINUS x*E!
       Try in t-space: t(1+t) d/dt + c + m*t = (1 + w_j) * E + c_j + z_j!
   (d) WAIT A SECOND: What about substituting t -> -1/(1+u) or y -> ...? First test (a)-(b) for N=2, 3, 4!
3. Test Laguerre F_{N, i}(z_1, ..., z_N) for N = 7, 8 with c_j = 1, lam_j = 1 AND GLZ lam = [1,0,1,0,...]!
"""
import numpy as np

def build_laguerre_ma(cs, lams):
    N = len(cs)
    polys = [{0: 1.0}]
    for j in range(N):
        cj = cs[j]; lj = lams[j]; bit = 1 << j
        new_p = [{} for _ in range(len(polys) + 1)]
        for i, d in enumerate(polys):
            for mask, coef in d.items():
                deg = mask.bit_count()
                new_p[i][mask] = new_p[i].get(mask, 0.0) + coef * (deg + cj)
                new_p[i][mask | bit] = new_p[i].get(mask | bit, 0.0) + coef
                if lj != 0.0:
                    new_p[i + 1][mask] = new_p[i + 1].get(mask, 0.0) + lj * coef
        polys = new_p
    return polys

def check_ma_uhp(d, N, rng, trials=150):
    max_im = -1e9
    for _ in range(trials):
        k = int(rng.integers(0, N)); bit = 1 << k
        zs = [complex(rng.normal(0, 2) * (10**rng.uniform(-1, 1)), 10**rng.uniform(-3, 1)) for _ in range(N)]
        A = 0j; B = 0j
        for mask, coef in d.items():
            has_k = bool(mask & bit)
            m = mask & (~bit)
            term = coef; idx = 0
            while m:
                if m & 1: term *= zs[idx]
                m >>= 1; idx += 1
            if has_k: B += term
            else: A += term
        if abs(B) < 1e-20: continue
        im = (-A / B).imag
        if im > max_im: max_im = im
    return max_im

rng = np.random.default_rng(12345)

print("=== 1. DISSECTING LAGUERRE MULTI-AFFINE STABILITY ===")
for mode in ["c=1, WILD lam>0", "c=1, GLZ lam=[1,0,...]", "WILD c>0, lam=1", "c_j=const c, lam=1", "c_j INCREASING, lam=1", "c_j DECREASING, lam=1"]:
    fails = 0; example = None
    for trial in range(20):
        N = 4
        if mode == "c=1, WILD lam>0":
            cs = [1.0]*N; lams = [float(10**rng.uniform(-2, 2)) for _ in range(N)]
        elif mode == "c=1, GLZ lam=[1,0,...]":
            N = 6; cs = [1.0]*N; lams = [1.0 if j%2==0 else 0.0 for j in range(N)]
        elif mode == "WILD c>0, lam=1":
            cs = [float(10**rng.uniform(-2, 2)) for _ in range(N)]; lams = [1.0]*N
        elif mode == "c_j=const c, lam=1":
            c0 = float(10**rng.uniform(-2, 2)); cs = [c0]*N; lams = [1.0]*N
        elif mode == "c_j INCREASING, lam=1":
            cs = sorted([float(10**rng.uniform(-2, 2)) for _ in range(N)]); lams = [1.0]*N
        elif mode == "c_j DECREASING, lam=1":
            cs = sorted([float(10**rng.uniform(-2, 2)) for _ in range(N)], reverse=True); lams = [1.0]*N
        polys = build_laguerre_ma(cs, lams)
        for i, d in enumerate(polys):
            if not d: continue
            mi = check_ma_uhp(d, N, rng, trials=120)
            if mi > 1e-7:
                fails += 1
                if example is None: example = (i, [round(x,3) for x in cs], [round(x,3) for x in lams], mi)
                break
        if mode == "c=1, GLZ lam=[1,0,...]":
            break  # deterministic!
    tot = 1 if mode == "c=1, GLZ lam=[1,0,...]" else 20
    print(f"  {mode:25s}: fails = {fails}/{tot}  example={example}", flush=True)

print("\n=== 2. CHECK GLZ PATTERN c=1, lam=[1,0,1,0,...] FOR N = 4, 6, 8 ===")
for N in [4, 6, 8]:
    cs = [1.0]*N; lams = [1.0 if j%2==0 else 0.0 for j in range(N)]
    polys = build_laguerre_ma(cs, lams)
    worst = -1e9
    for i, d in enumerate(polys):
        if not d: continue
        mi = check_ma_uhp(d, N, rng, trials=200)
        worst = max(worst, mi)
    print(f"  Laguerre GLZ-pattern N={N}: max Im(z_k) = {worst:+.3e}", flush=True)

print("\n=== 3. FULL MULTI-AFFINE GLZ (2N vars: w_1..w_N, z_1..z_N) ===")
# Operator_j:
#   Form A: (1 + w_j) * (E + c_j) + z_j   (wait: D_{c_j, m_j} = (1+t)(E + c_j) + (m_j - c_j)t)
#   Form B: (1 + w_j) * E + c_j + z_j     (wait: D_{c_j, m_j} = (1+t)E + c_j + m_j t)
#   Form C: E + c_j + w_j * (E + m_j)     (ONLY N variables w_1..w_N!! Wait: D_{c_j, m_j} = (E + c_j) + t*(E + m_j)!!)
#   WHOA!! STOP!! LOOK AT FORM C!!
#   D_{c_j, m_j} = (E + c_j) + w_j * (E + m_j)  HAS ONLY N VARIABLES w_1, ..., w_N,
#   AND SPECIALIZING w_1 = ... = w_N = t RECOVERS D_{c_j, m_j} EXACTLY!!
for form in ["Form C (N vars): (E+c_j) + w_j*(E+m_j) [GLZ m_j=j+1]",
             "Form C (N vars): (E+c_j) + w_j*(E+m_j) [const m=3]",
             "Form C (N vars): (E+c_j) + w_j*(E+m_j) [GLZ pairs]",
             "Form A (2N vars): (1+w_j)(E+1) + z_j",
             "Form B (2N vars): (1+w_j)E + 1 + z_j"]:
    if form.startswith("Form C"):
        for N in [3, 4, 5, 6]:
            if "GLZ pairs" in form and N % 2 == 1: continue
            if "m_j=j+1" in form:
                cs = [1.0]*N; ms = [float(j+2) for j in range(N)]; lams = [1.0]*N
            elif "const m=3" in form:
                cs = [1.0]*N; ms = [3.0]*N; lams = [1.0]*N
            else:
                # GLZ pairs: D_{k+1}(D_{k+1} + y), so m = [2,2,3,3,4,4], lam = [1,0,1,0,1,0]
                cs = [1.0]*N; ms = [float(j//2 + 2) for j in range(N)]; lams = [1.0 if j%2==0 else 0.0 for j in range(N)]
            polys = [{0: 1.0}]
            for j in range(N):
                cj = cs[j]; mj = ms[j]; lj = lams[j]; bit = 1 << j
                new_p = [{} for _ in range(len(polys) + 1)]
                for i, d in enumerate(polys):
                    for mask, coef in d.items():
                        deg = mask.bit_count()
                        new_p[i][mask] = new_p[i].get(mask, 0.0) + coef * (deg + cj)
                        new_p[i][mask | bit] = new_p[i].get(mask | bit, 0.0) + coef * (deg + mj)
                        if lj != 0.0:
                            new_p[i + 1][mask] = new_p[i + 1].get(mask, 0.0) + lj * coef
                polys = new_p
            worst = -1e9
            for i, d in enumerate(polys):
                if not d: continue
                mi = check_ma_uhp(d, N, rng, trials=150)
                worst = max(worst, mi)
            print(f"  {form} N={N}: max Im(w_k) = {worst:+.3e}", flush=True)
    else:
        N = 3; nv = 2*N
        polys = [{0: 1.0}]
        for j in range(N):
            bw = 1 << (2*j); bz = 1 << (2*j + 1)
            new_p = [{} for _ in range(len(polys) + 1)]
            for i, d in enumerate(polys):
                for mask, coef in d.items():
                    deg = mask.bit_count()
                    if "Form A" in form:
                        new_p[i][mask] = new_p[i].get(mask, 0.0) + coef * (deg + 1.0)
                        new_p[i][mask | bw] = new_p[i].get(mask | bw, 0.0) + coef * (deg + 1.0)
                        new_p[i][mask | bz] = new_p[i].get(mask | bz, 0.0) + coef
                    else:
                        new_p[i][mask] = new_p[i].get(mask, 0.0) + coef * (deg + 1.0)
                        if deg > 0:
                            new_p[i][mask | bw] = new_p[i].get(mask | bw, 0.0) + coef * deg
                        new_p[i][mask | bz] = new_p[i].get(mask | bz, 0.0) + coef
                    new_p[i + 1][mask] = new_p[i + 1].get(mask, 0.0) + coef
            polys = new_p
        worst = -1e9
        for i, d in enumerate(polys):
            if not d: continue
            mi = check_ma_uhp(d, nv, rng, trials=150)
            worst = max(worst, mi)
        print(f"  {form} N={N}: max Im = {worst:+.3e}", flush=True)
