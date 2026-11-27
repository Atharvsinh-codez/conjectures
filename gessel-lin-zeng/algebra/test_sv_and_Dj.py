import math
import random
import numpy as np
from collections import defaultdict

# === PART 1: Test combinatorial refinements of A_{k,i}(x) (Savage-Visontai style!) for k=3, 4 ===
def gen_jsp(k, S):
    """Generate all Jacobi-Stirling permutations for (k, S) where S subset {1..k} are deleted barred letters.
    Represent letters as integers: bar_j = 2*j - 1, j = 2*j (so bar_1 < 1 < bar_2 < 2 < ...)."""
    words = [()]
    for j in range(1, k + 1):
        # Insert bar_j (if j not in S) into every slot of each word, then insert (j, j) as adjacent pair into every slot!
        # Wait: check order! bar_j < j, so inserting bar_j first, then (j, j) as adjacent block!
        if j not in S:
            bj = 2*j - 1
            nxt = []
            for w in words:
                for pos in range(len(w) + 1):
                    nxt.append(w[:pos] + (bj,) + w[pos:])
            words = nxt
        uj = 2*j
        nxt = []
        for w in words:
            for pos in range(len(w) + 1):
                nxt.append(w[:pos] + (uj, uj) + w[pos:])
        words = nxt
    return words

def des(w):
    return sum(1 for r in range(len(w) - 1) if w[r] > w[r+1])

def roots_real_x(coeffs):
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    if len(coeffs) <= 1:
        return []
    r = np.roots(coeffs[::-1])
    if np.max(np.abs(r.imag)) > 1e-6:
        return None
    return list(np.sort(r.real))

def interlaces_x(p, q):
    rp, rq = roots_real_x(p[:]), roots_real_x(q[:])
    if rp is None or rq is None: return False
    dp, dq = len(rp), len(rq)
    eps = 1e-7
    if dq == dp:
        return all(rp[m] <= rq[m] + eps for m in range(dp)) and all(rq[m] <= rp[m+1] + eps for m in range(dp - 1))
    elif dq == dp + 1:
        return all(rq[m] <= rp[m] + eps for m in range(dp)) and all(rp[m] <= rq[m+1] + eps for m in range(dp))
    return False

def check_sv_seq(polys):
    nz = [p for p in polys if any(c != 0 for c in p)]
    for p in nz:
        if roots_real_x(p[:]) is None:
            return "NON-REAL"
    # Check strong interlacing fwd and bwd
    sf = all(interlaces_x(nz[a], nz[b]) for a in range(len(nz)) for b in range(a+1, len(nz)))
    sb = all(interlaces_x(nz[b], nz[a]) for a in range(len(nz)) for b in range(a+1, len(nz)))
    af = all(interlaces_x(nz[a], nz[a+1]) for a in range(len(nz)-1))
    ab = all(interlaces_x(nz[a+1], nz[a]) for a in range(len(nz)-1))
    return f"strong_fwd={sf}, strong_bwd={sb}, adj_fwd={af}, adj_bwd={ab}"

from itertools import combinations
print("=== PART 1: Combinatorial refinements of A_{k,i}(x) for (k,i) = (3,1), (3,2), (4,2) ===")
for k, i in [(3, 1), (3, 2), (4, 2)]:
    all_words = []
    for S in combinations(range(1, k + 1), i):
        for w in gen_jsp(k, set(S)):
            all_words.append((S, w))
    # Verify total matches A_{k,i}
    tot = defaultdict(int)
    for S, w in all_words: tot[des(w)] += 1
    print(f"\n--- (k,i)=({k},{i}), total words={len(all_words)}, deg={max(tot)} ---")
    stats = {
        "last_letter": lambda S, w: w[-1],
        "first_letter": lambda S, w: w[0],
        "pos_of_kk": lambda S, w: w.index(2*k),
        "pos_of_kk_from_end": lambda S, w: len(w) - 2 - w.index(2*k),
        "num_ascents": lambda S, w: sum(1 for r in range(len(w)-1) if w[r] < w[r+1]),
        "num_left_to_right_max": lambda S, w: sum(1 for r in range(len(w)) if all(w[j] < w[r] for j in range(r))),
        "val_before_kk": lambda S, w: w[w.index(2*k)-1] if w.index(2*k) > 0 else 0,
        "k_in_S": lambda S, w: 1 if k in S else 0,
        "last_letter_is_k": lambda S, w: 1 if w[-1] == 2*k else 0,
        "des_before_kk": lambda S, w: des(w[:w.index(2*k)]),
    }
    for sname, sfn in stats.items():
        buckets = defaultdict(lambda: defaultdict(int))
        for S, w in all_words:
            buckets[sfn(S, w)][des(w)] += 1
        keys = sorted(buckets.keys())
        polys = []
        for ky in keys:
            d = max(buckets[ky].keys())
            polys.append([buckets[ky][m] for m in range(d + 1)])
        print(f"  {sname:22s} ({len(keys):2d} buckets): {check_sv_seq(polys)}")

# === PART 2: Test D_j on Turan p1^2 >= p0*p2 (on all R) for j = 1, 2, 3, 4, 5, 10, 50! ===
print("\n=== PART 2: Does D_j preserve p1^2 >= p0*p2 (on all R) for small j = 1, 2, 3, 4, 5, 10, 50? ===")
def poly_mul(a, b):
    out = np.zeros(len(a) + len(b) - 1, dtype=float)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out
def apply_D_np(j, p):
    d = len(p) - 1
    out = np.zeros(d + 2, dtype=float)
    for m, cm in enumerate(p):
        out[m] += (m + 1) * cm
        out[m + 1] += (m + j) * cm
    return out
def eval_poly(p, x):
    acc = np.zeros_like(x, dtype=float)
    for c in reversed(p): acc = acc * x + c
    return acc
def poly_from_roots(rs):
    p = np.array([1.0])
    for r in rs: p = poly_mul(p, np.array([-r, 1.0]))
    return p

t_int = np.linspace(-1.0, 0.0, 401)
t_R = np.concatenate([-10.0**np.linspace(2, 0, 50) - 1.0, np.linspace(-2.0, -1.0, 50), t_int, np.linspace(0.0, 1.0, 50), 10.0**np.linspace(0, 2, 50)])

random.seed(303)
triples = []
for _ in range(300):
    d = random.choice([1, 2, 3])
    r0 = sorted([-random.uniform(0.02, 0.98) for _ in range(d + 2)])
    r1 = [random.uniform(r0[m]+1e-4, r0[m+1]-1e-4) for m in range(d + 1)]
    r2 = [random.uniform(r1[m]+1e-4, r1[m+1]-1e-4) for m in range(d)]
    p0, p1, p2 = poly_from_roots(r0), poly_from_roots(r1), poly_from_roots(r2)
    v0, v1, v2 = eval_poly(p0, t_R), eval_poly(p1, t_R), eval_poly(p2, t_R)
    pos = (v0 * v2) > 1e-30
    mr = max(1.0, float(np.max((v0[pos]*v2[pos])/(v1[pos]**2)))) if np.any(pos) else 1.0
    p2 = p2 / mr
    triples.append((p0, p1, p2))

for j in [1, 2, 3, 4, 5, 10, 50]:
    fails = 0
    worst = 999.0
    for (p0, p1, p2) in triples:
        q0, q1, q2 = apply_D_np(j, p0), apply_D_np(j, p1), apply_D_np(j, p2)
        u0, u1, u2 = eval_poly(q0, t_R), eval_poly(q1, t_R), eval_poly(q2, t_R)
        pos = (u0 * u2) > 1e-30
        if np.any(pos):
            m = float(np.min((u1[pos]**2)/(u0[pos]*u2[pos])))
            if m < worst: worst = m
            if m < 1.0 - 1e-5: fails += 1
    print(f"  D_{j:2d}: fails(<1) = {fails:3d}/300 | worst min ratio = {worst:.6f}")
