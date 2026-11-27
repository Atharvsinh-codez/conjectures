import numpy as np
import sympy as sp
from scipy.optimize import minimize

# 1. Expand tilde_h_{k, i}(t) in Jacobi bases P_n^{(alpha, beta)}(1 + 2t) for various (alpha, beta)!
#    Also check real-rootedness of the Jacobi-coefficient generating polynomial sum c_n z^n!
def apply_D_int(m, p):
    res = [0]*(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

KMAX = 8
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    H[k][0] = apply_D_int(m, apply_D_int(m, H[k-1][0]))
    for i in range(1, k):
        H[k][i] = apply_D_int(m, poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1]))
    H[k][k] = apply_D_int(m, H[k-1][k-1])

print("=== 1. JACOBI BASIS EXPANSION OF tilde_h_{k, i}((u-1)/2) ===")
u = sp.Symbol('u')
for k in [3, 4, 5, 6]:
    for alpha, beta_fn, bname in [
        (0, lambda k, i: 0, "Legendre(0,0)"),
        (0, lambda k, i: k, "Jacobi(0,k)"),
        (0, lambda k, i: k-1, "Jacobi(0,k-1)"),
        (0, lambda k, i: k+1, "Jacobi(0,k+1)"),
        (0, lambda k, i: k-i, "Jacobi(0,k-i)"),
        (1, lambda k, i: k, "Jacobi(1,k)"),
        (sp.Rational(-1, 2), lambda k, i: sp.Rational(-1, 2), "ChebyshevT"),
    ]:
        all_pos = True
        all_rz = True
        sample_c = None
        for i in range(k + 1):
            p = H[k][i]
            d = len(p) - 1
            if d == 0: continue
            # expr in u: t = (u - 1)/2
            expr = sp.expand(sum(p[r] * ((u - 1)/2)**r for r in range(d + 1)))
            bv = beta_fn(k, i)
            # Expand in P_n^{(alpha, bv)}(u) using greedy top-degree peeling!
            rem = expr
            coeffs = [0.0]*(d + 1)
            for n in range(d, -1, -1):
                Pn = sp.Poly(sp.jacobi(n, alpha, bv, u), u)
                lead_P = float(Pn.LC())
                rem_poly = sp.Poly(rem, u)
                lead_r = float(rem_poly.nth(n))
                c = lead_r / lead_P
                coeffs[n] = c
                rem = sp.expand(rem - c * Pn.as_expr())
            if min(coeffs) < -1e-9:
                all_pos = False
            # Check if sum coeffs[n] z^n is real-rooted!
            if d >= 2:
                rts = np.roots(coeffs[::-1])
                if np.max(np.abs(rts.imag)) > 1e-5:
                    all_rz = False
            if i == k // 2:
                sample_c = [round(c, 3) for c in coeffs]
        if k in [4, 6]:
            print(f"  k={k}, {bname:14s}: all_pos={str(all_pos):5s}, all_coeff_poly_RZ={str(all_rz):5s} | sample(i={k//2}): {sample_c[:6]}")

print("\n=== 2. SCIPY ADVERSARIAL OPTIMIZATION ON FULL COMBINED GLZ STEP G_j (WITH BALANCED SCALE!) ===")
def apply_D_np(b, p):
    res = np.zeros(len(p)+1)
    for n, c in enumerate(p):
        res[n] += c*(n+1); res[n+1] += c*(n+b)
    return res

def add_np(a, b):
    n = max(len(a), len(b))
    r = np.zeros(n); r[:len(a)] += a; r[:len(b)] += b
    return r

def wronskian_np(f, g):
    fp = np.arange(1, len(f))*f[1:] if len(f) > 1 else np.array([0.0])
    gp = np.arange(1, len(g))*g[1:] if len(g) > 1 else np.array([0.0])
    a = np.convolve(f, gp); b = np.convolve(fp, g)
    n = max(len(a), len(b))
    res = np.zeros(n); res[:len(a)] += a; res[:len(b)] -= b
    return res

def eval_np(p, t):
    v = np.zeros_like(t)
    for c in p[::-1]: v = v*t + c
    return v

def from_roots(rs, scale=1.0):
    p = np.array([scale], dtype=np.float64)
    for r in rs: p = np.convolve(p, [-r, 1.0])
    return p

ts = np.concatenate([
    -np.geomspace(1e-4, 1e2, 40),
    np.linspace(-0.999, -0.001, 120),
    np.geomspace(1e-4, 1e2, 40),
])

def min_J_np(X, M, Y, ts):
    w1 = eval_np(wronskian_np(X, M), ts)
    w2 = eval_np(wronskian_np(M, Y), ts)
    w12 = eval_np(wronskian_np(X, Y), ts)
    if np.min(w1) <= 0 or np.min(w2) <= 0: return -1.0
    with np.errstate(divide='ignore', invalid='ignore'):
        J = 4.0 * w1 * w2 / (w12 * w12)
    J = np.where(np.abs(w12) < 1e-14 * np.sqrt(w1 * w2), 1e12, J)
    return float(np.min(J))

# Optimize at d=(3,2,1,0) for j = 20, 40, 80 with s1 constrained to [0.5, 2.0] (balanced) AND unconstrained!
# Parameterize r0 = (-a3, -a2, -a1) in (-1, 0), r1 in (r0[0], r0[1]) x (r0[1], r0[2]), r2 in (r1[0], r1[1]), log(s1)
# Tune s2, s3 so J1 == 1 + 2/(j-1) and J2 == 1 + 2/(j-1) EXACTLY!
rng = np.random.default_rng(42)
for j in [20, 40, 80]:
    eps_in = 2.0 / (j - 1)
    eps_out = 2.0 / j
    for bal_label, s1_lim in [("BALANCED s1 in [0.5,2]", 0.7), ("FREE s1 in [1e-4,1e4]", 9.0)]:
        best_margin = 1e9
        best_info = None
        for restart in range(25):
            def obj(x):
                # Map x[0..5] to ordered roots in (-1, 0) via sigmoid
                sig = 1.0 / (1.0 + np.exp(-np.clip(x[:6], -10, 10)))
                # r0: 3 roots in (-1, 0)
                a1 = sig[0] * 0.9 + 0.001
                a2 = a1 + (0.995 - a1) * sig[1] * 0.9 + 1e-4
                a3 = a2 + (0.999 - a2) * sig[2] * 0.9 + 1e-4
                r0 = np.array([-a3, -a2, -a1])
                r1 = np.array([r0[0] + sig[3]*(r0[1]-r0[0]), r0[1] + sig[4]*(r0[2]-r0[1])])
                r2 = np.array([r1[0] + sig[5]*(r1[1]-r1[0])])
                s1 = np.exp(np.clip(x[6], -s1_lim, s1_lim))
                h0 = from_roots(r0, 1.0); h1_u = from_roots(r1, 1.0); h2_u = from_roots(r2, 1.0); h3_u = np.array([1.0])
                J1_u = min_J_np(h0, apply_D_np(j, h1_u), apply_D_np(j, apply_D_np(j, h2_u)), ts)
                J2_u = min_J_np(h1_u, apply_D_np(j, h2_u), apply_D_np(j, apply_D_np(j, h3_u)), ts)
                if J1_u <= 1e-6 or J2_u <= 1e-6: return 10.0
                s2 = s1*s1 * J1_u / (1.0 + eps_in)
                s3 = s2*s2 * J2_u / (s1 * (1.0 + eps_in))
                h1, h2, h3 = h1_u * s1, h2_u * s2, h3_u * s3
                # Enforce W(h0, D_j^3 h3) >= 0 with penalty!
                E = apply_D_np(j, apply_D_np(j, apply_D_np(j, h3)))
                w03 = np.min(eval_np(wronskian_np(h0, E), ts))
                if w03 < 0: return 5.0 - w03
                out0 = apply_D_np(j, add_np(h0, apply_D_np(j, h1)))
                out1 = apply_D_np(j+1, apply_D_np(j, add_np(h1, apply_D_np(j, h2))))
                out2 = apply_D_np(j+1, apply_D_np(j+1, apply_D_np(j, add_np(h2, apply_D_np(j, h3)))))
                Jout = min_J_np(out0, out1, out2, ts)
                return Jout - (1.0 + eps_out)
            x0 = rng.normal(0, 2.0, size=7)
            res = minimize(obj, x0, method='Nelder-Mead', options={'maxiter': 300})
            if res.fun < best_margin:
                best_margin = res.fun
                best_info = res.x
        print(f"  j={j:2d} [{bal_label:22s}]: min(Jout-(1+2/j)) = {best_margin:+.6f} (j^2*m = {j*j*best_margin:+8.2f})")
