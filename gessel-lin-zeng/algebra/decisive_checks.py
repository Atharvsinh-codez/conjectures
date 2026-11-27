"""
Decisive 10-second tests on BOTH major fronts:
1. Check if GLZ at k=34 (all alpha=1, beta>=1, w=1) + adversarial step 35..37 with alpha>=1, beta>=1
   breaks Turan (k=35), interlacing (k=36), and real-rootedness (k=37) — BOTH with arbitrary w_j > 0
   AND with fixed a_j=1, w_j=1, large b_j >= 2!
2. Test if D_j preserves R = W_AB W_CE / W_BC^2 <= 1 (or <= 4) on random 4-tuples (A,B,C,E)!
3. Test strengthened CUBIC C_c(t,y) = A - 3*c*B*y + 3*c*C*y^2 - E*y^3 for GLZ (k=6, 10, 14, 18):
   - What is the smallest c_min(k) such that C_c(t,y) is real-stable? (Expect ~ 1 - 2/k!)
   - Under Stage 1 (adding 2 overlapping strengthened Cubics with c = 1 - alpha/k), does c_out stay <= 1 - alpha/(k+1)
     WITH A FREE O(1/k^2) SURPLUS, AND AUTOMATICALLY KEEP x_1, x_2 <= 1?!
"""
import math, random
import numpy as np
import mpmath as mp
mp.mp.dps = 60

def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p
def poly_add(a, b):
    n = max(len(a), len(b))
    return poly_trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def poly_deriv(a):
    return [i * a[i] for i in range(1, len(a))] if len(a) > 1 else [0]
def apply_T(p, a, b):
    dp = poly_deriv(p)
    res = [0] * (len(p) + 2)
    for i, c in enumerate(p):
        res[i] += a * c; res[i+1] += b * c
    for i, c in enumerate(dp):
        if c == 0: continue
        res[i+1] += c; res[i+2] += c
    return poly_trim(res)
def poly_eval(p, x):
    acc = mp.mpc(0) if isinstance(x, mp.mpc) else mp.mpf(0)
    for c in reversed(p): acc = acc * x + c
    return acc

print("=== 1. TEST GLZ k=34 + ADVERSARIAL STEP 35 (alpha>=1, beta>=1) ===", flush=True)
ht = {1: {0: [1, 2], 1: [1]}}
for k in range(2, 35):
    ht[k] = {}
    for idx in range(k + 1):
        p1 = apply_T(ht[k-1][idx], 1, k + 1) if idx <= k - 1 else [0]
        p2 = ht[k-1][idx-1] if idx >= 1 else [0]
        ht[k][idx] = apply_T(poly_add(p1, p2), 1, k + 1)

p34 = ht[34]
tv = -mp.mpf('1e-16')
i = 17
f0, df0 = poly_eval(p34[i-1], tv), poly_eval(poly_deriv(p34[i-1]), tv)
f1, df1 = poly_eval(p34[i], tv),   poly_eval(poly_deriv(p34[i]), tv)
f2, df2 = poly_eval(p34[i+1], tv), poly_eval(poly_deriv(p34[i+1]), tv)
w01 = f0*df1 - df0*f1; w12 = f1*df2 - df1*f2; w02 = f0*df2 - df0*f2
print(f"  GLZ k=34, i=17, t=-1e-16: 4W01W12/W02^2 = {mp.nstr(4*w01*w12/(w02*w02), 8)}", flush=True)
delta = f1*f1 - f0*f2
ddelta_u = tv*(1+tv)*(2*f1*df1 - df0*f2 - f0*df2)
phi_star = -ddelta_u / (2*delta)
print(f"  Optimal phi_star = {mp.nstr(phi_star, 6)}", flush=True)
# Check with a35 = 1, b35 = 1 + beta35 >= 2: Phi = 1 + b35*tv = 1 - b35*1e-16!
# So b35 = (1 - phi_star) / 1e-16!
b35 = (1 - phi_star) / (-tv)
print(f"  With a35=1: needed b35 = {mp.nstr(b35, 6)} (b35 >= 2: {b35 >= 2})", flush=True)
if b35 >= 2:
    # Check Turan of T_{1, b35} p34 + w * p34 at tv for w=0 and w=1!
    for w35 in [mp.mpf(0), mp.mpf('1e-10'), mp.mpf(1)]:
        v_im1 = tv*(1+tv)*df0 + phi_star*f0 + (w35*poly_eval(p34[i-2], tv))
        v_i   = tv*(1+tv)*df1 + phi_star*f1 + (w35*f0)
        v_ip1 = tv*(1+tv)*df2 + phi_star*f2 + (w35*f1)
        print(f"    w35={mp.nstr(w35,3)}: q35_i^2 - q35_{{i-1}}*q35_{{i+1}} = {mp.nstr(v_i*v_i - v_im1*v_ip1, 6)}", flush=True)

print("\n=== 2. CHECK IF D_j PRESERVES R = W_AB W_CE / W_BC^2 <= 1 (OR <= 4) ON RANDOM 4-TUPLES ===", flush=True)
def pfromr_np(rts):
    p = np.array([1.0])
    for r in rts: p = np.polymul(p, np.array([1.0, -r]))
    return p[::-1]
def apply_D_np(j, p):
    d = len(p) - 1
    out = np.zeros(d + 2)
    for m, cm in enumerate(p):
        out[m] += (m + 1)*cm; out[m+1] += (m + j)*cm
    return out
def wronsk_np(f, g):
    fh, gh = f[::-1], g[::-1]
    return np.polysub(np.polymul(fh, np.polyder(gh)), np.polymul(np.polyder(fh), gh))

ts = np.concatenate([np.linspace(-10, -1.02, 80), np.linspace(-0.99, -0.01, 200), np.linspace(0.01, 10, 80)])
random.seed(7)
tested = 0; max_R_out_when_in_le_1 = 0.0; max_ratio_D = 0.0
for _ in range(5000):
    if tested >= 300: break
    # 4 polynomials of degree 3 with A << B << C << E in (-1, 0)
    pts = sorted([random.uniform(-0.95, -0.05) for _ in range(12)])
    # A: 0,4,8; B: 1,5,9; C: 2,6,10; E: 3,7,11 (so A << B << C << E AND A << E!)
    A = pfromr_np([pts[0], pts[4], pts[8]])
    B = pfromr_np([pts[1], pts[5], pts[9]]) * (10**random.uniform(-1,1))
    C = pfromr_np([pts[2], pts[6], pts[10]]) * (10**random.uniform(-1,1))
    E = pfromr_np([pts[3], pts[7], pts[11]]) * (10**random.uniform(-1,1))
    wAB = np.polyval(wronsk_np(A, B), ts)
    wBC = np.polyval(wronsk_np(B, C), ts)
    wCE = np.polyval(wronsk_np(C, E), ts)
    wAC = np.polyval(wronsk_np(A, C), ts)
    wBE = np.polyval(wronsk_np(B, E), ts)
    if np.min(4*wAB*wBC - wAC**2) <= 0 or np.min(4*wBC*wCE - wBE**2) <= 0: continue
    R_in = float(np.max(wAB*wCE / (wBC**2)))
    if R_in > 1.0: continue
    tested += 1
    j = random.choice([2, 5, 10, 30])
    DA, DB, DC, DE = apply_D_np(j, A), apply_D_np(j, B), apply_D_np(j, C), apply_D_np(j, E)
    dwAB = np.polyval(wronsk_np(DA, DB), ts)
    dwBC = np.polyval(wronsk_np(DB, DC), ts)
    dwCE = np.polyval(wronsk_np(DC, DE), ts)
    R_out = float(np.max(dwAB*dwCE / (dwBC**2)))
    if R_out > max_R_out_when_in_le_1: max_R_out_when_in_le_1 = R_out
    if R_out / R_in > max_ratio_D: max_ratio_D = R_out / R_in
print(f"  Random 4-tuples with J>=1 & R_in <= 1 ({tested} tested): max R_out after D_j = {max_R_out_when_in_le_1:.4f}, max R_out/R_in = {max_ratio_D:.4f}", flush=True)

print("\n=== 3. CHECK STRENGTHENED CUBIC C_c = A - 3*c*B*y + 3*c*C*y^2 - E*y^3 FOR GLZ ===", flush=True)
# Find smallest c_min(k) where C_c is real-stable on UHP test pts for k=6, 10, 14, 18!
t_uhp = [mp.mpc(-10, '1e-6'), mp.mpc(-1, '1e-6'), mp.mpc(-0.5, '1e-6'), mp.mpc(-0.01, '1e-6'),
         mp.mpc(0, '1e-6'), mp.mpc(0.05, '1e-6'), mp.mpc(0.2, '1e-6'), mp.mpc(1, '1e-6'), mp.mpc(10, '1e-6')]
for k in [6, 10, 14, 18]:
    i = k // 2 + 1
    A = ht[k-1][i-3]
    B = apply_T(ht[k-1][i-2], 1, k+1)
    C = apply_T(apply_T(ht[k-1][i-1], 1, k+1), 1, k+1)
    E = apply_T(apply_T(apply_T(ht[k-1][i], 1, k+1), 1, k+1), 1, k+1)
    # Binary search smallest c in [0.5, 1.0] where max Im(y) <= 0!
    lo, hi = mp.mpf('0.5'), mp.mpf('1.0')
    for _ in range(25):
        mid = (lo + hi) / 2
        ok = True
        for tv in t_uhp:
            cy = [poly_eval(A, tv), -3*mid*poly_eval(B, tv), 3*mid*poly_eval(C, tv), -poly_eval(E, tv)]
            rts = mp.polyroots(list(reversed(cy)), maxsteps=200, extraprec=30)
            if max(r.imag for r in rts) > mp.mpf('1e-12'):
                ok = False; break
        if ok: hi = mid
        else: lo = mid
    print(f"  GLZ k={k:2d}, i={i}: c_min = {mp.nstr(hi, 6)}  (k*(1 - c_min) = {mp.nstr(k*(1-hi), 5)})", flush=True)
