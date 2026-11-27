import mpmath as mp

# Test the M = infinity (t-lifted) quad: (tilde_h_{k, i-1}, t * tilde_h_{k, i}, t^2 * tilde_h_{k, i+1})!
# Check on GLZ up to k = 40, including t = +-1e-5, +-1e-10, +-1e-16, (-1, 0), and all R!
# ALSO test Fixed-M (M = k+2 vs M = 2k vs M = 1000) at k = 20, 30, 40!
mp.mp.dps = 120

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
def shift_t(p): return [0] + p
def eval_mp(p, t):
    v = mp.mpf(0)
    for c in reversed(p): v = v * t + c
    return v
def eval_W(f, g, t):
    return eval_mp(f, t)*eval_mp(deriv_int(g), t) - eval_mp(deriv_int(f), t)*eval_mp(g, t)

ts = [mp.mpf(s) for s in [
    "-100", "-10", "-2", "-1.1", "-1.001", "-1", "-0.999", "-0.95", "-0.8", "-0.5", "-0.2",
    "-0.05", "-1e-3", "-1e-5", "-1e-8", "-1e-12", "-1e-16",
    "1e-16", "1e-12", "1e-8", "1e-5", "1e-4", "1e-3", "0.01", "0.1", "1", "10", "100"
]]

KMAX = 40
row = [[1, 2], [1]]
for k in range(2, KMAX + 1):
    m = k + 1
    new_row = [None]*(k + 1)
    new_row[0] = apply_D_int(m, apply_D_int(m, row[0]))
    for i in range(1, k):
        new_row[i] = apply_D_int(m, poly_add_int(apply_D_int(m, row[i]), row[i-1]))
    new_row[k] = apply_D_int(m, row[k-1])
    row = new_row

    if k in [5, 10, 15, 20, 25, 30, 35, 40]:
        # Check t-lifted quad (h_{i-1}, t h_i, t^2 h_{i+1}) across ALL i (or dense sample of i)
        mn_J = mp.mpf('1e100'); loc = None
        is_chk = list(range(1, k)) if k <= 20 else sorted(set([1, 2, 3, 4, 5, k//4, k//2-1, k//2, k//2+1, 3*k//4, k-3, k-2, k-1]))
        for i in is_chk:
            X = row[i-1]
            M = shift_t(row[i])
            Y = shift_t(shift_t(row[i+1]))
            for tv in ts:
                w1 = eval_W(X, M, tv); w2 = eval_W(M, Y, tv); w12 = eval_W(X, Y, tv)
                if w1 <= 0 or w2 <= 0:
                    mn_J = mp.mpf(-1); loc = (i, float(tv), "W<=0"); break
                J = 4*w1*w2/(w12*w12) if w12 != 0 else mp.mpf('1e10')
                if J < mn_J: mn_J = J; loc = (i, float(tv))
            if mn_J < 0: break
        # Also check at t = +alpha / s_{k,i} where B_k(x, -xy) failed!
        for i in [1, 2, k//2, k-1]:
            X = row[i-1]; M = shift_t(row[i]); Y = shift_t(shift_t(row[i+1]))
            s1 = mp.mpf(row[i][1]) / mp.mpf(row[i][0])
            for mult in ["0.01", "0.1", "0.5", "1", "2", "10", "100"]:
                for sgn in [-1, 1]:
                    tv = sgn * mp.mpf(mult) / s1
                    w1 = eval_W(X, M, tv); w2 = eval_W(M, Y, tv); w12 = eval_W(X, Y, tv)
                    if w1 <= 0 or w2 <= 0:
                        mn_J = mp.mpf(-1); loc = (i, float(tv), "W<=0_sc"); break
                    J = 4*w1*w2/(w12*w12) if w12 != 0 else mp.mpf('1e10')
                    if J < mn_J: mn_J = J; loc = (i, float(tv), "sc")
        print(f"  k={k:2d}: min J(h_{{i-1}}, t h_i, t^2 h_{{i+1}}) = {float(mn_J):+9.5f} (k*(J-1)={k*(float(mn_J)-1):+7.3f}) at {loc}")
