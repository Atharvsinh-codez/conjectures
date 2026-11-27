import random
import math
import numpy as np

def poly_mul(a, b):
    out = np.zeros(len(a) + len(b) - 1, dtype=float)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out

def apply_D_np(b, p):
    # p[m] is coeff of t^m
    d = len(p) - 1
    out = np.zeros(d + 2, dtype=float)
    for m, cm in enumerate(p):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + b) * cm
    return out

def add_np(a, b):
    n = max(len(a), len(b))
    out = np.zeros(n, dtype=float)
    out[:len(a)] += a
    out[:len(b)] += b
    return out

def eval_poly_c(p, z):
    acc = 0j
    for c in reversed(p):
        acc = acc * z + c
    return acc

def test_uhp_stable(polys_in_y, test_ts):
    """polys_in_y[i] is coeff of y^i (as array of t-coeffs low->high).
    Returns max Im(y) over all t in test_ts."""
    k = len(polys_in_y) - 1
    max_im = -1e9
    worst = None
    for tv in test_ts:
        yc = [eval_poly_c(polys_in_y[i], tv) for i in range(k, -1, -1)]
        if abs(yc[0]) < 1e-15:
            continue
        mx = max(abs(c) for c in yc)
        yc = [c / mx for c in yc]
        r = np.roots(yc)
        mi = max(z.imag for z in r)
        if mi > max_im:
            max_im = mi
            worst = tv
    return max_im, worst

# Dense adversarial grid of t in UHP
test_ts = []
for re in [-20.0, -5.0, -1.5, -1.05, -0.95, -0.7, -0.5, -0.3, -0.1, -0.01, 0.0, 0.01, 0.1, 0.5, 1.0, 5.0, 20.0]:
    for im in [1e-5, 1e-3, 0.02, 0.2, 2.0]:
        test_ts.append(complex(re, im))

# === TEST 1: Does M_{b, q}(F) = D_b F(t, y) + y F(t, q y) preserve real stability on
#             random products F(t, y) = prod (a_r + b_r t + c_r y) with 0 < a_r < b_r, c_r > 0?! ===
print("=== TEST 1: One step of M_{b, q} on random real-stable F(t, y) = prod_{r=1}^m (a_r + b_r t + c_r y), 0<a_r<b_r ===")
random.seed(7)
for q in [1.0, 0.5, 0.25, 0.1, 0.01]:
    fails = 0
    worst_info = None
    for trial in range(60):
        m = random.choice([1, 2, 3, 4])
        # Build F(t, y) = prod (a_r + b_r t + c_r y)
        F = [np.array([1.0])]
        for _ in range(m):
            ar = random.uniform(0.1, 5.0)
            br = ar + random.uniform(0.1, 20.0) # so -ar/br in (-1, 0)
            cr = 10 ** random.uniform(-2, 2)
            # Multiply F by (ar + br t) + cr * y
            newF = [np.zeros(1) for _ in range(len(F) + 1)]
            lin = np.array([ar, br])
            for i in range(len(F)):
                newF[i] = add_np(newF[i], poly_mul(F[i], lin))
                newF[i + 1] = add_np(newF[i + 1], F[i] * cr)
            F = newF
        # Also allow extra t-factors with roots in (-1, 0)
        for _ in range(random.choice([0, 1, 2])):
            ar = random.uniform(0.1, 5.0)
            br = ar + random.uniform(0.1, 20.0)
            lin = np.array([ar, br])
            F = [poly_mul(p, lin) for p in F]

        b = random.choice([2.0, 3.0, 5.0, 20.0, 100.0])
        # Apply M_{b, q}: G[i] = D_b F[i] + q^{i-1} F[i-1]
        k = len(F) - 1
        G = []
        for i in range(k + 2):
            t1 = apply_D_np(b, F[i]) if i <= k else np.array([0.0])
            t2 = (q ** (i - 1)) * F[i - 1] if i >= 1 else np.array([0.0])
            G.append(add_np(t1, t2))
        mi, wt = test_uhp_stable(G, test_ts)
        if mi > 1e-6:
            fails += 1
            if worst_info is None:
                worst_info = (m, b, mi, wt)
    print(f"  q={q:4.2f}: fails={fails:2d}/60 | first fail={worst_info}")

# === TEST 2: Wait! What about F(t, y) = p_0(t) + p_1(t) y (general degree-1 in y with p_1 << p_0 in (-1,0))?! ===
print("\n=== TEST 2: One step of M_{b, q} on arbitrary interlacing pair p_1 << p_0 in (-1, 0) ===")
for q in [1.0, 0.5, 0.25, 0.1, 0.01]:
    fails = 0
    worst_info = None
    for trial in range(100):
        d = random.choice([1, 2, 3, 4])
        # Pick 2d+1 sorted random roots in (-1, 0)
        pts = sorted([-random.uniform(0.01, 0.99) for _ in range(2*d + 1)])
        r0 = pts[0::2] # d+1 roots for p_0
        r1 = pts[1::2] # d roots for p_1 -> p_1 << p_0!
        p0 = np.array([1.0])
        for r in r0: p0 = poly_mul(p0, np.array([-r, 1.0]))
        p1 = np.array([1.0])
        for r in r1: p1 = poly_mul(p1, np.array([-r, 1.0]))
        # Random positive scaling!
        p1 = p1 * (10 ** random.uniform(-3, 3))
        b = random.choice([1.0, 2.0, 3.0, 10.0, 100.0])
        G = [apply_D_np(b, p0), add_np(apply_D_np(b, p1), p0), q * p1]
        mi, wt = test_uhp_stable(G, test_ts)
        if mi > 1e-6:
            fails += 1
            if worst_info is None:
                worst_info = (d, b, mi, wt)
    print(f"  q={q:4.2f} (pair p_1<<p_0): fails={fails:3d}/100 | first fail={worst_info}")

# === TEST 3: What about Q_k(t, y) = sum q^{i(i-1)/2} q_{k,i}(t) y^i for REAL-ROOTED IN t for fixed y > 0?
# Wait: what about q=1/2 on RANDOM (p_0 + p_1 y + p_2 y^2) with deg_t dropping by 1 and q-stable?! ===
print("\n=== TEST 3: Test q-operator M_{b,q} on F(t, y) where F ITSELF has roots only in (-1,0) at y=0 and deg_t(F_i) = d - i ===")
# Specifically, take F generated by 1, 2, 3 steps of M_{b_j, q} with random b_j, AND also perturb F slightly!
for q in [0.5, 0.25, 0.1]:
    fails = 0
    for trial in range(50):
        # Generate arbitrary F(t, y) of degree 2 in y that is stable by taking prod_{r=1}^2 (u_r(t) + v_r(t) y)
        # where v_r << u_r in (-1, 0), deg u_r = deg v_r + 1!
        F = [np.array([1.0])]
        for _ in range(2):
            d = random.choice([1, 2])
            pts = sorted([-random.uniform(0.02, 0.98) for _ in range(2*d + 1)])
            u = np.array([1.0])
            for r in pts[0::2]: u = poly_mul(u, np.array([-r, 1.0]))
            v = np.array([1.0])
            for r in pts[1::2]: v = poly_mul(v, np.array([-r, 1.0]))
            v = v * (10 ** random.uniform(-2, 2))
            newF = [np.zeros(1) for _ in range(len(F) + 1)]
            for i in range(len(F)):
                newF[i] = add_np(newF[i], poly_mul(F[i], u))
                newF[i + 1] = add_np(newF[i + 1], poly_mul(F[i], v))
            F = newF
        # Notice F is real-stable with ALL t-roots at y=0 in (-1, 0)!
        # Apply M_{b, q}:
        b = random.choice([2.0, 5.0, 50.0])
        k = len(F) - 1
        G = [add_np(apply_D_np(b, F[i]) if i <= k else np.array([0.0]),
                    (q ** (i - 1)) * F[i - 1] if i >= 1 else np.array([0.0]))
             for i in range(k + 2)]
        mi, wt = test_uhp_stable(G, test_ts)
        if mi > 1e-6:
            fails += 1
    print(f"  q={q:4.2f} (product of 2 interlacing factors): fails={fails:2d}/50")
