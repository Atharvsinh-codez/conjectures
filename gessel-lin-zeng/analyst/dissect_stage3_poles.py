import numpy as np
import mpmath as mp

# 1. Dissect the EXACT 6 poles of J(j; t) = 4 W_XM(j) W_MY(j) / W_XY(j)^2 on:
#    (A) The task-695 worst-case counterexample at j = 80, 160 (where P <= 3.5, L <= 3, yet j^2*drop = +9.57!)
#    (B) GLZ at k = 6, 10, 14 (where j^2*dlogJ >= -0.281 everywhere on R!)
#
#    Goal: See the EXACT structural difference between (A) and (B)!

def apply_D(b, p):
    res = np.zeros(len(p)+1)
    for n, c in enumerate(p):
        res[n] += c*(n+1); res[n+1] += c*(n+b)
    return res

def wronskian(f, g):
    fp = np.arange(1, len(f))*f[1:] if len(f) > 1 else np.array([0.0])
    gp = np.arange(1, len(g))*g[1:] if len(g) > 1 else np.array([0.0])
    a = np.convolve(f, gp); b = np.convolve(fp, g)
    n = max(len(a), len(b))
    res = np.zeros(n); res[:len(a)] += a; res[:len(b)] -= b
    return res

def eval_p(p, t):
    v = np.zeros_like(t, dtype=np.float64)
    for c in p[::-1]: v = v*t + c
    return v

def eval_scalar(p, t):
    v = 0.0
    for c in p[::-1]: v = v*t + c
    return v

def from_roots(rs, scale=1.0):
    p = np.array([scale], dtype=np.float64)
    for r in rs: p = np.convolve(p, [-r, 1.0])
    return p

def poles_of_J(h0, h1, h2, tv):
    D0h1 = apply_D(0.0, h1); th1 = np.concatenate([[0.0], h1])
    D0h2 = apply_D(0.0, h2)
    Y0 = apply_D(0.0, D0h2)
    tD0h2 = np.concatenate([[0.0], D0h2])
    t1th2 = np.convolve([0.0, 1.0, 1.0], h2)
    n1 = max(len(tD0h2), len(t1th2))
    Y1 = np.zeros(n1); Y1[:len(tD0h2)] += 2*tD0h2; Y1[:len(t1th2)] += t1th2
    Y2 = np.concatenate([[0.0, 0.0], h2])
    def W_at(f, g): return eval_scalar(wronskian(f, g), tv)
    a0, a1 = W_at(h0, D0h1), W_at(h0, th1)
    b0, b1, b2 = W_at(h0, Y0), W_at(h0, Y1), W_at(h0, Y2)
    c0 = W_at(D0h1, Y0)
    c1 = W_at(D0h1, Y1) + W_at(th1, Y0)
    c2 = W_at(D0h1, Y2) + W_at(th1, Y1)
    c3 = W_at(th1, Y2)
    p1 = np.roots([a1, a0])
    q12 = np.roots([b2, b1, b0])
    s123 = np.roots([c3, c2, c1, c0])
    return p1, q12, s123

ts = np.concatenate([
    -np.geomspace(1e-4, 1e2, 40),
    np.linspace(-0.999, -0.001, 140),
    np.geomspace(1e-4, 1e2, 40),
])

def delta_vec(h_hi, h_lo, j, ts):
    th = np.concatenate([[0.0], h_lo])
    Dh = apply_D(j, h_lo)
    w_num = eval_p(wronskian(h_hi, th), ts)
    w_den = eval_p(wronskian(h_hi, Dh), ts)
    if np.min(w_den) <= 0: return None
    return w_num / w_den

def min_J_with_loc(X, M, Y, ts):
    w1 = eval_p(wronskian(X, M), ts)
    w2 = eval_p(wronskian(M, Y), ts)
    w12 = eval_p(wronskian(X, Y), ts)
    if np.min(w1) <= 0 or np.min(w2) <= 0: return -1.0, 0.0
    with np.errstate(divide='ignore', invalid='ignore'):
        J = 4.0 * w1 * w2 / (w12 * w12)
    J = np.where(np.abs(w12) < 1e-14 * np.sqrt(w1 * w2), 1e12, J)
    idx = int(np.argmin(J))
    return float(J[idx]), float(ts[idx])

print("=== 1. FIND AND INSPECT WORST task-695 FAILURE AT j = 160 (P <= 3.5, L <= 3) ===")
rng = np.random.default_rng(777)
j = 160
worst = (-1e9, None)
for _ in range(15000):
    if rng.random() < 0.6:
        a = rng.uniform(0.02, 0.4)
        b = rng.uniform(a + 0.05, 0.95)
        r0 = np.sort([-b, -0.5*(a+b) + rng.uniform(-0.1, 0.1)*(b-a), -a])
        r1 = np.array([0.5*(r0[0]+r0[1]) + rng.uniform(-0.15, 0.15)*(r0[1]-r0[0]),
                       0.5*(r0[1]+r0[2]) + rng.uniform(-0.15, 0.15)*(r0[2]-r0[1])])
        r2 = np.array([0.5*(r1[0]+r1[1]) + rng.uniform(-0.15, 0.15)*(r1[1]-r1[0])])
    else:
        r0 = np.sort(-rng.uniform(0.01, 0.99, size=3))
        r1 = np.array([rng.uniform(r0[0], r0[1]), rng.uniform(r0[1], r0[2])])
        r2 = np.array([rng.uniform(r1[0], r1[1])])
    h0, h1, h2 = from_roots(r0), from_roots(r1), from_roots(r2)
    d1 = delta_vec(h0, h1, j, ts)
    d2 = delta_vec(h1, h2, j, ts)
    if d1 is None or d2 is None: continue
    P = j * max(float(np.max(d1)), float(np.max(d2)))
    if P > 3.5: continue
    L = (j**2) * float(np.max(np.abs(d2 - d1)))
    if L > 3.0: continue
    Mj = apply_D(j, h1); Yj = apply_D(j, apply_D(j, h2))
    Ju, t_in = min_J_with_loc(h0, Mj, Yj, ts)
    if Ju <= 0: continue
    c = np.sqrt((1.0 + 2.0/j) / Ju)
    Mjp = apply_D(j+1, h1)*c; Yjp = apply_D(j+1, apply_D(j+1, h2))
    Jout, t_out = min_J_with_loc(h0, Mjp, Yjp, ts)
    drop_j2 = (j**2) * ((1.0 + 2.0/j) - Jout)
    if drop_j2 > worst[0]:
        worst = (drop_j2, (r0, r1, r2, P, L, t_in, t_out, Jout, h0, h1, h2, c))

drop_j2, (r0, r1, r2, P, L, t_in, t_out, Jout, h0, h1, h2, c) = worst
print(f"  Worst j=160: j^2*drop = {drop_j2:+.3f}, P={P:.3f}, L={L:.3f}, t_in={t_in:.5f}, t_out={t_out:.5f}")
print(f"  r0={np.round(r0, 4)}, r1={np.round(r1, 4)}, r2={np.round(r2, 4)}")
p1, q12, s123 = poles_of_J(h0, h1, h2, t_out)
print(f"  Poles at t_out={t_out:.5f}:")
print(f"    W_XM (1 pole):  {np.round(p1, 3)}")
print(f"    W_XY (2 poles): {np.round(q12, 3)}")
print(f"    W_MY (3 poles): {np.round(s123, 3)}")
print(f"  Notice Re(q12) vs j=160! Check J(j; t_out) = {min_J_with_loc(h0, apply_D(j, h1)*c, apply_D(j, apply_D(j, h2)), np.array([t_out]))[0]:.6f} -> J(j+1; t_out) = {Jout:.6f}")

print("\n=== 2. COMPARE WITH GLZ POLES AT k = 10 (j = 11) AND k = 14 (j = 15) ===")
def apply_D_int(m, p):
    res = [0]*(len(p)+1)
    for r, cc in enumerate(p):
        res[r] += cc*(r+1); res[r+1] += cc*(r+m)
    return res
def poly_add_int(a, b):
    res = [0]*max(len(a), len(b))
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    return res

H = [[None]*(k+1) for k in range(15)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, 15):
    m = k + 1
    H[k][0] = apply_D_int(m, apply_D_int(m, H[k-1][0]))
    for i in range(1, k):
        H[k][i] = apply_D_int(m, poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1]))
    H[k][k] = apply_D_int(m, H[k-1][k-1])

for k in [8, 12]:
    for i in [k//2, k-1]:
        h0 = np.array(H[k][i-1], dtype=np.float64)
        h1 = np.array(H[k][i], dtype=np.float64)
        h2 = np.array(H[k][i+1], dtype=np.float64)
        for tv in [-0.5, -0.05, 0.5]:
            p1, q12, s123 = poles_of_J(h0, h1, h2, tv)
            # compute exact d/dj log J at j = k+1
            jv = k + 1
            dlog = np.sum(1.0/(jv - p1)) + np.sum(1.0/(jv - s123)) - 2*np.sum(1.0/(jv - q12))
            print(f"  GLZ k={k:2d}, i={i:2d}, t={tv:+5.2f} (j={jv}): j^2*dlogJ={jv*jv*dlog.real:+.4f} | p1={np.round(p1,2)} | q12={np.round(q12,2)} | s123={np.round(s123,2)}")
