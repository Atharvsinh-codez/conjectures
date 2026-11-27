"""
Four fast, high-signal checks (< 15s total):
1. GLZ Stage 3 d/dj ln J_j(t) across FINE grid of R: (-inf, -1), (-1, 0) down to -1e-12, and (0, +inf)!
   Check for k = 6, 10, 14, 18 and ALL i = 2..k!
2. Universal Stage 3 drop on t >= 0 and t <= -1 for arbitrary d=2,3 real-stable quads with roots in (-1, 0)!
   (Remember at d=2, adversarial drop O(1/j) happened INSIDE (-1, 0) at t ~ -a! What about OUTSIDE (-1, 0)?!)
3. Brändén-Solus symmetric decomposition A_{k,i}(x) = a_{k,i}(x) + x b_{k,i}(x) for k = 2..12:
   Are a_{k,i} and b_{k,i} real-rooted and interlacing b_{k,i} << a_{k,i}?!
4. Combinatorial refinement `first_base` of JSP_{5, i} (k = 5, i = 1, 2, 3):
   Does all-pairs interlacing survive k = 5?
"""
import numpy as np
import mpmath as mp
from scipy.optimize import minimize

mp.mp.dps = 60

def apply_Dj_mp(p, j):
    out = [mp.mpf(0)] * (len(p) + 1)
    for m, c in enumerate(p):
        if c != 0:
            out[m] += (1 + m) * c
            out[m + 1] += (j + m) * c
    return out

def poly_add_mp(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return out

def poly_eval_deriv(p, t):
    v = mp.mpf(0); dv = mp.mpf(0)
    for c in reversed(p):
        dv = dv * t + v
        v = v * t + c
    return v, dv

def mul_t(p):
    return [mp.mpf(0)] + list(p)

# ============================================================
# 1. GLZ STAGE 3: EXACT d/dj ln J_j(t) ON THREE REGIONS:
#    (-inf, -1), (-1, 0) [down to -1e-12], and (0, +inf)
# ============================================================
print("=== 1. GLZ STAGE 3: min k^2 * (ln J_{k+2} - ln J_{k+1}) BY REGION ===")
h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
h_store = {1: h}
for k in range(2, 19):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = poly_add_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = poly_add_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
    if k in [6, 10, 14, 18]:
        h_store[k] = h

regions = {
    "t<-1": [mp.mpf(x) for x in [-100.0, -10.0, -3.0, -1.5, -1.1, -1.01, -1.001]],
    "(-1,0)": [mp.mpf(x) for x in [-0.999, -0.99, -0.9, -0.7, -0.5, -0.3, -0.1, -0.03, -0.01, -1e-3, -1e-4, -1e-6, -1e-9, -1e-12]],
    "t>=0": [mp.mpf(x) for x in [0.0, 1e-6, 1e-3, 0.02, 0.1, 0.2, 0.5, 1.0, 3.0, 10.0, 100.0]],
}

for k in [6, 10, 14, 18]:
    hk = h_store[k]
    reg_min = {name: (1e9, None) for name in regions}
    for i in range(2, k + 1):
        f0, f1, f2 = hk[i - 2], hk[i - 1], hk[i]
        M_a = apply_Dj_mp(f1, k + 1); M_b = apply_Dj_mp(f1, k + 2)
        Y_a = apply_Dj_mp(apply_Dj_mp(f2, k + 1), k + 1)
        Y_b = apply_Dj_mp(apply_Dj_mp(f2, k + 2), k + 2)
        for rname, tlist in regions.items():
            for t in tlist:
                vX = poly_eval_deriv(f0, t)
                vMa = poly_eval_deriv(M_a, t); vYa = poly_eval_deriv(Y_a, t)
                vMb = poly_eval_deriv(M_b, t); vYb = poly_eval_deriv(Y_b, t)
                def W(a, b): return a[0]*b[1] - a[1]*b[0]
                Ja = 4*W(vX, vMa)*W(vMa, vYa) / (W(vX, vYa)**2)
                Jb = 4*W(vX, vMb)*W(vMb, vYb) / (W(vX, vYb)**2)
                val = float((k**2) * mp.log(Jb / Ja))
                if val < reg_min[rname][0]:
                    reg_min[rname] = (val, (i, float(t)))
    print(f"  k={k:2d}: " + " | ".join(f"{rn}: min k^2*dlnJ={v[0]:+.4f} at {v[1]}" for rn, v in reg_min.items()))

# ============================================================
# 2. UNIVERSAL STAGE 3 DROP OUTSIDE (-1, 0) (t >= 0 and t <= -1) AT d = 2!
# ============================================================
print("\n=== 2. ADVERSARIAL d=2 STAGE 3 DROP ON t >= 0 AND t <= -1 ===")
# f2 = 1, f1 = t + a, f0 = (t + b1)(t + b2) with 0 < b1 < a < b2 < 1!
# Wait: require J_j(t) >= 1 on ALL of R (or check J_j >= 1 on R and minimize J_{j+1}/J_j on t >= 0 and t <= -1)!
def eval_J_d2(a, b1, b2, j, t):
    # f2 = 1
    # D_j(1) = 1 + j t
    # D_j^2(1) = D_j(1 + j t) = t(1+t)*j + (1+jt)^2 = 1 + 3 j t + (j + j^2) t^2
    # Wait: check D_j(1 + j t): (t+t^2)*j + (1+jt)(1+jt) = 1 + (2j+j)t + (j^2+j)t^2 = 1 + 3jt + j(j+1)t^2!
    # f1 = t + a -> D_j f1 = (t+t^2)*1 + (1+jt)(t+a) = a + (1 + a*j + 1) t + (j + 1) t^2
    #   = a + (2 + a*j) t + (j + 1) t^2!
    # f0 = b1*b2 + (b1+b2) t + t^2
    X  = b1*b2 + (b1+b2)*t + t*t
    Xp = (b1+b2) + 2*t
    M  = a + (2 + a*j)*t + (j + 1)*t*t
    Mp = (2 + a*j) + 2*(j + 1)*t
    Y  = 1 + 3*j*t + j*(j + 1)*t*t
    Yp = 3*j + 2*j*(j + 1)*t
    wXM = X*Mp - Xp*M
    wMY = M*Yp - Mp*Y
    wXY = X*Yp - Xp*Y
    return 4*wXM*wMY / (wXY**2), wXM, wMY, wXY

ts_check = np.concatenate([np.linspace(-5.0, -1.05, 25), np.linspace(-0.99, -0.01, 40), np.linspace(0.0, 5.0, 25), [20.0, 100.0]])
for j in [20, 40, 80, 160]:
    worst_pos = (1e9, None)
    worst_neg = (1e9, None)
    rng = np.random.default_rng(77)
    for _ in range(8000):
        # sample b1 < a < b2 in (0, 1) with log-uniform near 0 and 1 too!
        r = np.sort(np.exp(rng.uniform(-8, 0, 3)))
        if rng.random() < 0.5:
            r = 1.0 - r[::-1]
        b1, a, b2 = r[0], r[1], r[2]
        # Check J_j >= 1 on ts_check and W > 0
        ok = True
        for tc in ts_check:
            J, w1, w2, w3 = eval_J_d2(a, b1, b2, j, tc)
            if w1 <= 0 or w2 <= 0 or w3 <= 0 or J < 1.0:
                ok = False
                break
        if not ok:
            continue
        # Check t >= 0 and t <= -1
        for tc in [0.0, 1e-4, 1e-2, 0.1, 0.5, 2.0, 10.0, 100.0]:
            Ja = eval_J_d2(a, b1, b2, j, tc)[0]
            Jb = eval_J_d2(a, b1, b2, j + 1, tc)[0]
            val = (j**2) * np.log(Jb / Ja)
            if val < worst_pos[0]:
                worst_pos = (val, (a, b1, b2, tc, Ja, Jb))
        for tc in [-1.001, -1.05, -1.5, -5.0, -50.0]:
            Ja = eval_J_d2(a, b1, b2, j, tc)[0]
            Jb = eval_J_d2(a, b1, b2, j + 1, tc)[0]
            val = (j**2) * np.log(Jb / Ja)
            if val < worst_neg[0]:
                worst_neg = (val, (a, b1, b2, tc, Ja, Jb))
    print(f"  d=2, j={j:3d}: worst t>=0 j^2*dlnJ = {worst_pos[0]:+.4f} (at t={worst_pos[1][3] if worst_pos[1] else None}), worst t<=-1 j^2*dlnJ = {worst_neg[0]:+.4f} (at t={worst_neg[1][3] if worst_neg[1] else None})")

# ============================================================
# 3. BRÄNDÉN-SOLUS SYMMETRIC DECOMPOSITION A_{k,i} = a_{k,i} + x b_{k,i}
# ============================================================
print("\n=== 3. BRÄNDÉN-SOLUS SYMMETRIC DECOMPOSITION A_{k,i}(x) = a(x) + x b(x) ===")
def apply_SL(p, L):
    out = [0.0]*(len(p)+1)
    for m, c in enumerate(p):
        out[m] += (1+m)*c
        out[m+1] += (L-m)*c
    while len(out)>1 and abs(out[-1]) < 1e-12: out.pop()
    return out

def symm_decomp(p):
    # p(x) = a(x) + x b(x), where a(x) = x^d a(1/x) (deg d), b(x) = x^{d-1} b(1/x) (deg d-1)
    # Exact formula: a(x) = (p(x) - x p*(x)) / (1 - x) ?? Check!
    # p*(x) = x^d p(1/x) = a(x) + b(x). So p(x) - x p*(x) = (1-x) a(x) => a(x) = (p(x) - x p*(x)) / (1 - x)!
    # Wait: check! p*(x) = x^d (a(1/x) + (1/x) b(1/x)) = a(x) + b(x)!
    # Then p(x) - x p*(x) = a(x) + x b(x) - x a(x) - x b(x) = (1 - x) a(x)!
    # And b(x) = (p*(x) - p(x)) / (1 - x)!
    d = len(p) - 1
    if d == 0:
        return list(p), [0.0]
    pstar = list(reversed(p))
    # num_a = p(x) - x p*(x)
    num_a = [0.0]*(d + 2)
    for m in range(d + 1):
        num_a[m] += p[m]
        num_a[m + 1] -= pstar[m]
    # divide num_a by (1 - x): a[0] = num_a[0], a[m] = num_a[m] + a[m-1]
    a = [0.0]*(d + 1)
    a[0] = num_a[0]
    for m in range(1, d + 1):
        a[m] = num_a[m] + a[m - 1]
    # b[m] = pstar[m] - a[m] for m = 0..d-1
    b = [pstar[m] - a[m] for m in range(d)]
    return a, b

def is_real_rooted(p):
    while len(p) > 1 and abs(p[-1]) < 1e-10 * max(1.0, max(abs(x) for x in p)):
        p = p[:-1]
    if len(p) <= 1: return True, np.array([])
    r = np.roots(list(reversed(p)))
    if np.max(np.abs(r.imag)) > 1e-6:
        return False, None
    return True, np.sort(r.real)

A = {(1, 0): [1.0, 1.0], (1, 1): [1.0]}  # A_{1,0} = 1 + x! (Check: S_1(1) = 1 + x!)
for k in range(2, 11):
    for i in range(k + 1):
        L = 3*k - 2 - i
        u = [0.0]
        if i <= k - 1:
            t1 = apply_SL(A[(k-1, i)], L - 1)
            u = [0.0]*max(len(u), len(t1))
            for m, v in enumerate(t1): u[m] += v
        if i >= 1:
            t2 = A[(k-1, i-1)]
            if len(t2) > len(u): u += [0.0]*(len(t2)-len(u))
            for m, v in enumerate(t2): u[m] += v
        A[(k, i)] = apply_SL(u, L)

fails_a = 0; fails_b = 0; fails_int = 0; neg_coeff = 0; tot = 0
for k in range(2, 11):
    for i in range(k + 1):
        tot += 1
        a, b = symm_decomp(A[(k, i)])
        if min(a) < -1e-9 or min(b) < -1e-9:
            neg_coeff += 1
        oka, ra = is_real_rooted(a)
        okb, rb = is_real_rooted(b)
        if not oka: fails_a += 1
        if not okb: fails_b += 1
        if oka and okb and len(rb) > 0:
            # check b << a (ra[0] < rb[0] < ra[1] < ... < rb[d-2] < ra[d-1])
            if len(ra) == len(rb) + 1:
                if not all(ra[m] < rb[m] < ra[m+1] for m in range(len(rb))):
                    fails_int += 1
print(f"  A_{{k,i}} = a + x b (k=2..10, tot={tot}): neg_coeff={neg_coeff}, fails_a={fails_a}, fails_b={fails_b}, fails_interlace={fails_int}")

# Also check SYMMETRIC DECOMPOSITION OF h~_{k,i}(t)!!
# Wait: h~_{k,i}(t) = (1+t)^d A_{k,i}(t/(1+t))!
print("  (Also checking A_{k,i}(x) reversed / shifted symmetric decomposition in t:)")
h_f = {(1, 0): [1.0, 2.0], (1, 1): [1.0]}
for k in range(2, 10):
    for i in range(k + 1):
        u = [0.0]
        if i <= k - 1:
            t1 = [0.0]*(len(h_f[(k-1,i)])+1)
            for m, c in enumerate(h_f[(k-1,i)]):
                t1[m] += (1+m)*c; t1[m+1] += (k+1+m)*c
            u = t1
        if i >= 1:
            t2 = h_f[(k-1,i-1)]
            if len(t2) > len(u): u += [0.0]*(len(t2)-len(u))
            for m, v in enumerate(t2): u[m] += v
        res = [0.0]*(len(u)+1)
        for m, c in enumerate(u):
            res[m] += (1+m)*c; res[m+1] += (k+1+m)*c
        h_f[(k, i)] = res
fa_t = 0; fb_t = 0; nc_t = 0; tot_t = 0
for k in range(2, 10):
    for i in range(k + 1):
        tot_t += 1
        a, b = symm_decomp(h_f[(k, i)])
        if min(a) < -1e-9 or min(b) < -1e-9: nc_t += 1
        if not is_real_rooted(a)[0]: fa_t += 1
        if not is_real_rooted(b)[0]: fb_t += 1
print(f"  h~_{{k,i}} = a + t b (k=2..9, tot={tot_t}): neg_coeff={nc_t}, fails_a={fa_t}, fails_b={fb_t}")

# ============================================================
# 4. CHECK `first_base` AT k = 5 (finish combinatorial refinements!)
# ============================================================
print("\n=== 4. CHECK `first_base` REFINEMENT OF JSP_{5, i} (k = 5) ===")
# Build JSP_k recursively by inserting letter k (3 copies or 2 copies):
# In JSP_{k, S}: letter j has 2 copies if j in S (singleton? wait: deg A_{k,i} = 2k-i-1,
# for i=k, A_{k,k} has deg k-1 -> 2 copies per letter? Let's check k=1:
# JSP_{1, {1}}: word "11" (0 descents) -> A_{1,1} = 1!
# JSP_{1, {}}:  word "111"?? Wait: A_{1,0} = 1 + x (2 words, 1 with 0 descents, 1 with 1 descent)!
# Remember Jacobi-Stirling permutations: letters 1..k, with 2 copies of j (for j not in S -> 3 copies? No, 2 copies + 1 bar!):
# Wait: for k=1, i=0: 1_1 1_2 vs 1_2 1_1? Let's use the exact recursion on words:
# Start at k=1: i=1 -> {(1,1): 1}, i=0 -> {(1,2): 1, (2,1): 1}? No, check k=2, i=2: A_{2,2} = 1+2x (sum=3).
# Inserting "22" into "11" in 3 slots: 2211 (des 1), 1221 (des 1), 1122 (des 0) => 1 + 2x! EXACT!
# Inserting "22" with bar into "11":Wait: for A_{1,0} = 1+x (sum 2), A_{2,0} = S_4(S_3(1+x)) = (1+4x)(1+4x+2x^2)??
# Wait: S_3(1+x) = (1+3x)(1+x) + x(1-x) = 1 + 5x + 2x^2 (sum 8), S_4(1+5x+2x^2) = 1 + 14x + 28x^2 + 6x^3? No,
# A_{2,0}(1) = 1 + 14 + 28 + 6? Check sum of S_L(p) at x=1: (1+L) p(1)!
# So A_{1,0}(1) = 2, A_{2,0}(1) = 5 * 4 * 2 = 40? Wait: 1*2 * 4*5 = 40!
# Exact JSP_{k, S} construction (from Gessel-Lin-Zeng):
# Each letter j not in S has TWO copies {j, j'} (or three copies {j, j, j} with Stirling rule?):
# Wait: inserting TWO copies of j into a word of length M gives M+1 slots for the pair (jj);
# inserting THREE copies (or two distinct pairs): (M+1) slots for j_1 j_1, then (M+2) slots for j_2?
# For k=1, S={}: length 0 -> 1 slot for 1_1 1_1, then 2 slots for 1_2 (Wait: 1_1 1_1 1_2 vs 1_2 1_1 1_1 -> 1+x)!
# For k=2, S={}: length 3 -> 4 slots for 2_1 2_1 (gives length 5), then 5 slots? Wait: 4*5 = 20, times 2 = 40!
# YES!! That's the EXACT Jacobi-Stirling permutation rule:
#   - if j in S: insert the adjacent pair (j, j) into any of the (M+1) slots (Wait: 1 copy of j_1 j_1);
#   - if j not in S: first insert (j, j) (Wait: (j, j) as a pair, then ONE single j' > j? Or first single j', then pair (j, j)?)
# Check A_{2,1}(1) = 2+12+6 = 20:
#   S={1}: k=1 gives (1,1) [len 2], k=2 not in S: 3 slots * 4 slots? Wait: 3*4 = 12!
#   S={2}: k=1 gives len 3 (2 words), k=2 in S: 4 slots => 2*4 = 8! Total = 12 + 8 = 20!! EXACT!!
# Check internal descents: for k=1, S={}: 1_1 1_1 1_2 (0 des) vs 1_2 1_1 1_1 (1 des if 1_2 > 1_1!) -> 1 + x!
# Let's represent j_1 as (j, 0) and j_2 as (j, 1)!
# Step k: if k in S (so i increases by 1): insert pair [(k,0),(k,0)] at any slot 0..len(w) (keeping Stirling order!).
#         if k not in S: first insert [(k,0),(k,0)] at any slot, THEN insert [(k,1)] at any slot? Wait: check S={1}, k=2!
#         From [(1,0),(1,0)] (len 2): insert [(2,0),(2,0)] (3 slots) -> len 4, then insert [(2,1)] (5 slots? No, 3*5=15 != 12!)
#         Instead: first insert [(k,1)] (len+1 = 3 slots -> len 3), THEN insert [(k,0),(k,0)] (len+1 = 4 slots) -> 3*4 = 12!!
#         CHECK descents on S={1}, k=2 (should give S_3(S_2(1)) = S_3(1+2x) = 1+7x+4x^2, sum 12!):
#         Wait: S_2(1) = 1+2x (sum 3), S_3(1+2x) = (1+3x)(1+2x) + x(1-x)*2 = 1 + 7x + 4x^2!
# Let's verify that exact rule in 5 lines!
def count_des(w):
    return sum(1 for idx in range(len(w) - 1) if w[idx] > w[idx + 1])

def step_words(words_dict, k, in_S):
    # words_dict: word -> count
    out = {}
    if in_S:
        # insert pair ((k,0),(k,0)) at any of len(w)+1 slots
        for w, cnt in words_dict.items():
            for pos in range(len(w) + 1):
                nw = w[:pos] + ((k, 0), (k, 0)) + w[pos:]
                out[nw] = out.get(nw, 0) + cnt
    else:
        # first insert single (k, 1), then insert pair ((k,0),(k,0))
        tmp = {}
        for w, cnt in words_dict.items():
            for pos in range(len(w) + 1):
                nw = w[:pos] + ((k, 1),) + w[pos:]
                tmp[nw] = tmp.get(nw, 0) + cnt
        for w, cnt in tmp.items():
            for pos in range(len(w) + 1):
                nw = w[:pos] + ((k, 0), (k, 0)) + w[pos:]
                out[nw] = out.get(nw, 0) + cnt
    return out

# Verify k=2, i=1 matches A_{2,1} = [2, 12, 6]!
w1_in = {((1,0),(1,0)): 1}
w1_out = {((1,0),(1,0),(1,1)): 1, ((1,1),(1,0),(1,0)): 1}
w2_1 = {}
for nw, c in step_words(w1_in, 2, False).items(): w2_1[nw] = w2_1.get(nw, 0) + c
for nw, c in step_words(w1_out, 2, True).items(): w2_1[nw] = w2_1.get(nw, 0) + c
poly21 = [0]*4
for w, c in w2_1.items(): poly21[count_des(w)] += c
print(f"  Verify JSP rule at (2,1): {poly21} (expected [2, 12, 6])")
# Check k=5, i=3 (smallest non-trivial k=5 case!) by first_base = w[0][0] in {1..5}:
# Build JSP_{k, i} for k=5, i=3 (wait: how many words? A_{5,3}(1): let's check size first!)
print(f"  A_{{5,3}}(1) = {sum(A[(5,3)]):.0f}, A_{{5,4}}(1) = {sum(A[(5,4)]):.0f}")
# A_{5,4}(1) is small (~few thousand), A_{5,3}(1) is ~few ten thousand! Do DP tracking ONLY (first_base, length, des)!
# EVEN FASTER AND EXACT!! When inserting into a word with (first_base, length, des):
# Actually inserting at pos=0 changes first_base to k, inserting at pos>0 keeps first_base!
# Wait: inserting a single (k, 1) (larger than ALL existing letters!):
#   - at pos=0: new first_base = k, creates +1 descent (since (k,1) > old w[0])!
#   - at pos=len(w): keeps first_base, creates +0 descents!
#   - at an existing descent pos (w[pos-1] > w[pos]): keeps des (+0)!
#   - at an existing non-descent pos (w[pos-1] <= w[pos]): creates +1 descent!
# Inserting a pair ((k,0),(k,0)) (larger than all < k, smaller than (k,1)):
#   Wait: in `in_S = False`, (k, 1) IS ALREADY IN THE WORD! So (k,0) is smaller than (k,1)!
#   Check the exact DP for `first_base` at k=5, i=1,2,3,4 using exact (first_base, des) — wait:
#   to be 100% sure of (k,0) vs (k,1) interaction, let's check via exact word enumeration for (5, 4) and (4, 2) first, or just DP with flag "is w[0] == (k,1)"!
#   Actually for k=5, i=4: A_{5,4}(1) is only ~2000 words! Let's check (5,4) and (5,3) directly if small enough!
dp_words = {(1, 0): w1_out, (1, 1): w1_in}
for k in range(2, 6):
    new_dp = {}
    for i in range(k + 1):
        if k == 5 and i < 3:
            continue  # skip huge i=0,1,2 at k=5
        if k == 4 and i < 1:
            continue
        cur = {}
        if (k - 1, i) in dp_words:
            for nw, c in step_words(dp_words[(k-1, i)], k, False).items():
                cur[nw] = cur.get(nw, 0) + c
        if (k - 1, i - 1) in dp_words:
            for nw, c in step_words(dp_words[(k-1, i-1)], k, True).items():
                cur[nw] = cur.get(nw, 0) + c
        new_dp[(k, i)] = cur
    dp_words = new_dp

for (k, i) in [(5, 3), (5, 4)]:
    by_fb = {}
    for w, c in dp_words[(k, i)].items():
        fb = w[0][0]
        d = count_des(w)
        if fb not in by_fb: by_fb[fb] = {}
        by_fb[fb][d] = by_fb[fb].get(d, 0) + c
    polys = {}
    for fb in sorted(by_fb):
        mx = max(by_fb[fb])
        polys[fb] = [by_fb[fb].get(m, 0) for m in range(mx + 1)]
    # Check pairwise interlacing for all 1 <= a < b <= 5!
    fails_p = []
    for a in sorted(polys):
        for b in sorted(polys):
            if a >= b: continue
            oka, ra = is_real_rooted(polys[a])
            okb, rb = is_real_rooted(polys[b])
            if not (oka and okb):
                fails_p.append((a, b, "non-real"))
                continue
            # check if either ra << rb or rb << ra
            def chk(r1, r2):
                m, n = len(r1), len(r2)
                if n == m + 1: return all(r2[q] <= r1[q] <= r2[q+1] for q in range(m))
                if n == m: return all(r1[q] <= r2[q] for q in range(m)) and all(r2[q] <= r1[q+1] for q in range(m-1))
                return False
            if not (chk(ra, rb) or chk(rb, ra)):
                fails_p.append((a, b, "no-interlace"))
    print(f"  JSP_{{{k},{i}}} `first_base` (total words={sum(dp_words[(k,i)].values())}): pairwise fails = {fails_p}")
