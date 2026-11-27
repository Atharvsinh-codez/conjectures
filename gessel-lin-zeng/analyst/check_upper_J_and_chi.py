import mpmath as mp

# Quick check: on GLZ (k=3..18), what is the MAXIMUM (not minimum!) of J(h0, D_j h1, D_j^2 h2) over ALL t in R,
# and what is max R = W01*W23/W12^2 over ALL t in R?!
# Recall: W_{AE} > 0 <==> chi = W01*W23/(W02*W13) = (sqrt(J1*J2)/4) * sqrt(R) < 1 !!
# So if J <= 4 (or J1*J2*R < 16) EVERYWHERE on R, then W_{AE} > 0 FOLLOWS AUTOMATICALLY!!

mp.mp.dps = 60

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

def deriv_int(p): return [r*p[r] for r in range(1, len(p))] if len(p) > 1 else [0]
def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v
def eval_W(f, g, t):
    return eval_mp(f, t)*eval_mp(deriv_int(g), t) - eval_mp(deriv_int(f), t)*eval_mp(g, t)

KMAX = 16
H = [[None]*(k+1) for k in range(KMAX+1)]
H[1][0] = [1, 2]; H[1][1] = [1]
for k in range(2, KMAX+1):
    m = k + 1
    H[k][0] = apply_D_int(m, apply_D_int(m, H[k-1][0]))
    for i in range(1, k):
        H[k][i] = apply_D_int(m, poly_add_int(apply_D_int(m, H[k-1][i]), H[k-1][i-1]))
    H[k][k] = apply_D_int(m, H[k-1][k-1])

ts = [mp.mpf(s) for s in [
    "-100", "-10", "-2", "-1.1", "-1.01", "-1", "-0.999", "-0.99", "-0.95", "-0.9",
    "-0.8", "-0.7", "-0.6", "-0.5", "-0.4", "-0.3", "-0.2", "-0.1", "-0.05", "-0.01",
    "-1e-4", "0", "0.01", "0.1", "1", "10", "100"
]]

print("=== MAX J, MIN J, MAX R, AND MAX chi ON GLZ BASE QUADS (h_{i-1}, D_j h_i, D_j^2 h_{i+1}, D_j^3 h_{i+2}) ===")
for k in range(4, 16):
    j = k + 2  # next-step j!
    min_J = mp.mpf('1e100'); max_J = mp.mpf(0); max_R = mp.mpf(0); max_chi = mp.mpf(0)
    loc_maxJ = None; loc_maxchi = None
    for i in range(1, k):
        A = H[k][i-1]
        B = apply_D_int(j, H[k][i])
        C = apply_D_int(j, apply_D_int(j, H[k][i+1]))
        for tv in ts:
            w01 = eval_W(A, B, tv); w12 = eval_W(B, C, tv); w02 = eval_W(A, C, tv)
            J = 4*w01*w12/(w02*w02)
            if J < min_J: min_J = J
            if J > max_J: max_J = J; loc_maxJ = (i, float(tv))
        if i + 2 <= k:
            E = apply_D_int(j, apply_D_int(j, apply_D_int(j, H[k][i+2])))
            for tv in ts:
                w01 = eval_W(A, B, tv); w12 = eval_W(B, C, tv); w23 = eval_W(C, E, tv)
                w02 = eval_W(A, C, tv); w13 = eval_W(B, E, tv)
                R = (w01*w23)/(w12*w12)
                chi = (w01*w23)/(w02*w13)
                if R > max_R: max_R = R
                if chi > max_chi: max_chi = chi; loc_maxchi = (i, float(tv))
    print(f"  k={k:2d} (j={j:2d}): J in [{float(min_J):.5f}, {float(max_J):7.4f}@{loc_maxJ}] | max_R={float(max_R):.5f} | max_chi={float(max_chi):.5f}@{loc_maxchi}")
