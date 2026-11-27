import numpy as np

# Test Boolean-lattice Rayleigh inequality on (-1, 0):
#   tilde_h_{k, R u {a}}(t) * tilde_h_{k, R u {b}}(t) - tilde_h_{k, R}(t) * tilde_h_{k, R u {a,b}}(t) >= 0 on (-1, 0)
# Recall h_{k, S}(t) = prod_{j=1}^k (D_j if j in S else D_j^2)(1), and tilde_h_{k, S} = h_{k, S} / (1+t).
# Wait: since (1+t)^2 cancels, h_{k, S} gives the exact same sign on (-1, 0)!

def apply_D(m, p):
    res = np.zeros(len(p)+1)
    for r, c in enumerate(p):
        res[r] += c*(r+1); res[r+1] += c*(r+m)
    return res

def eval_p(p, ts):
    v = np.zeros_like(ts)
    for c in p[::-1]: v = v*ts + c
    return v

ts = np.linspace(-0.999, -0.001, 400)

for k in [3, 4, 5, 6, 7]:
    # Build all 2^k polynomials h_{k, S}
    H = {0: np.array([1.0])}
    for j in range(1, k + 1):
        newH = {}
        for mask, p in H.items():
            # j not in S -> D_j^2
            newH[mask] = apply_D(j, apply_D(j, p))
            # j in S -> D_j
            newH[mask | (1 << (j - 1))] = apply_D(j, p)
        H = newH
    # Check ALL (R, a, b) with a < b not in R!
    total = 0; fails = 0; min_rel = 1e9; first_fail = None
    for a in range(k):
        for b in range(a + 1, k):
            for R in range(1 << k):
                if (R & (1 << a)) or (R & (1 << b)): continue
                total += 1
                va = eval_p(H[R | (1 << a)], ts)
                vb = eval_p(H[R | (1 << b)], ts)
                v0 = eval_p(H[R], ts)
                vab = eval_p(H[R | (1 << a) | (1 << b)], ts)
                diff = va * vb - v0 * vab
                rel = float(np.min(diff / (np.abs(va*vb) + np.abs(v0*vab) + 1e-300)))
                if rel < min_rel: min_rel = rel
                if rel < -1e-7:
                    fails += 1
                    if first_fail is None:
                        first_fail = (bin(R), a+1, b+1, rel, float(ts[np.argmin(diff / (np.abs(va*vb)+np.abs(v0*vab)+1e-300))]))
    print(f"  k={k}: total={total:4d}, fails={fails:4d}, min_rel={min_rel:+.6f}, first_fail={first_fail}")
