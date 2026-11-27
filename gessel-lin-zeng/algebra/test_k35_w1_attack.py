"""
CRITICAL 5-SECOND TEST OF ROUTE II (Universal Q_{k,q} for a_j=1, w_j=1, b_j >= 2):
At GLZ k=35, we KNOW 4 W(h18,h17)W(h17,h16)/W(h18,h16)^2 = 0.997128 < 1 at t=0!
Wait: remember Exact T-Turan Equivalence:
  (D_b h17)^2 - (D_b h16)(D_b h18) = A0 + Phi * B0 + Phi^2 * C0, where Phi = 1 + b t!
At t = -eps (tiny negative eps > 0), Phi = 1 - b * eps can be tuned to ANY real value < 1
by choosing b = (1 - Phi*) / eps > 2 (as long as Phi* < 1)!
CHECK THE SIGN OF Phi* = -B0 / (2 C0) AT t = -eps!!
WAIT A SECOND!!
At t = 0:
  h_i(0) = C(35, i) > 0, so C0 = h17(0)^2 - h16(0)h18(0) > 0!
  And what is Phi* = -B0 / (2 C0) = -[2 h17 h17' - (h16 h18' + h16' h18)] * t(1+t) / (2 C0) ??
  AT t = -eps: t(1+t) ~ -eps!
  WAIT A SECOND!
  B0 = t(1+t) * Delta'(t) where Delta(t) = h17(t)^2 - h16(t) h18(t)!
  If Delta'(0) > 0, then at t = -eps < 0:
    B0 = (-eps) * Delta'(0) < 0  =>  Phi* = -B0 / (2 C0) > 0!
    AND since B0 ~ O(eps), Phi* ~ O(eps) << 1!!
    SO b* = (1 - Phi*) / eps ~ 1/eps >> 2 IS POSITIVE AND HUGE!!
LET'S COMPUTE THE EXACT b36* AND CHECK WHETHER:
1. (D_{b36*} h35_17)^2 < (D_{b36*} h35_16)(D_{b36*} h35_18) at t = -eps!
2. With w36 = 1 (NOT w36 = 1e-30!), does U36_17 = D_{b36*} h35_17 + h35_16 STILL violate Turan?!
3. With a further step b37 >= 2, w37 = 1, can we break interlacing or real-rootedness of h37, h38?!
"""
import mpmath as mp
mp.mp.dps = 80

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

h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
for k in range(2, 36):
    j = k + 1
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = poly_add_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = poly_add_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn

# Inspect i = 17 at t = -eps!
f0, f1, f2 = h[16], h[17], h[18]  # note deg f0 > deg f1 > deg f2
eps = mp.mpf('1e-16')
t = -eps
v0, d0 = poly_eval_deriv(f0, t)
v1, d1 = poly_eval_deriv(f1, t)
v2, d2 = poly_eval_deriv(f2, t)
# u_deriv = t(1+t) f'
u0, u1, u2 = t*(1+t)*d0, t*(1+t)*d1, t*(1+t)*d2
A0 = u1**2 - u0*u2
B0 = 2*u1*v1 - (u0*v2 + u2*v0)
C0 = v1**2 - v0*v2
Phi_star = -B0 / (2*C0)
b_star = (1 - Phi_star) / (-t)
min_val = A0 + Phi_star*B0 + (Phi_star**2)*C0
print(f"At t = -1e-16:")
print(f"  C0 = {float(C0):.6e}, B0 = {float(B0):.6e}, A0 = {float(A0):.6e}")
print(f"  Phi* = {float(Phi_star):.6e}, b* = {float(b_star):.6e}")
print(f"  min [(D_{{b*}} h17)^2 - (D_{{b*}} h16)(D_{{b*}} h18)] = {float(min_val):.6e} (< 0 confirms Turan break!)")
print(f"  (D_{{b*}} h17)(t) = {float(u1 + Phi_star*v1):.6e}, h16(t) = {float(v0):.6e}")
print(f"  Ratio |D_{{b*}} h17| / |h16| at t = {float(abs(u1 + Phi_star*v1) / abs(v0)):.6e}!!")

# NOW CHECK WITH w36 = 1!!
# U36_i = D_{b*} h35_i + h35_{i-1} (WITH w36 = 1!)
# Look at the size of D_{b*} h35_i(t) vs h35_{i-1}(t) at t = -eps!!
# D_{b*} h35_i(t) = u_i + Phi* v_i ~ O(eps * v_i') + O(eps * v_i) ~ O(eps * 2^35 * v_i)!
# WAIT A SECOND!
# v_i' / v_i ~ 10^10 at k=35, while eps = 1e-16!
# So |D_{b*} h_i(t)| ~ 1e-6 * |h_i(t)| << |h_{i-1}(t)|!!
# THAT MEANS WHEN w36 = 1, h35_{i-1}(t) COMPLETELY SWAMPS D_{b*} h35_i(t) AT t = -1e-16!!
# WHAT IF eps IS CHOSEN LARGER SO |D_{b*} h_i(t)| >> |h_{i-1}(t)|?!
# Check for t = -1e-14, -1e-12, -1e-10, -1e-8, -1e-6, -1e-4, -1e-2!!
print("\n=== CHECK U36_i = D_{b*} h35_i + h35_{i-1} (w36 = 1!) ACROSS t = -10^{-m} ===")
for m in [16, 13, 11, 10, 9, 8, 7, 6, 5, 4, 3]:
    t = -mp.mpf(10)**(-m)
    v0, d0 = poly_eval_deriv(f0, t)
    v1, d1 = poly_eval_deriv(f1, t)
    v2, d2 = poly_eval_deriv(f2, t)
    # Also need h[15] for U36_16 = D_{b*} h16 + h15!
    vm1, dm1 = poly_eval_deriv(h[15], t)
    u0, u1, u2 = t*(1+t)*d0, t*(1+t)*d1, t*(1+t)*d2
    A0 = u1**2 - u0*u2
    B0 = 2*u1*v1 - (u0*v2 + u2*v0)
    C0 = v1**2 - v0*v2
    if C0 <= 0:
        print(f"  t=-1e-{m}: C0 <= 0 already!")
        continue
    Phi_star = -B0 / (2*C0)
    b_star = (1 - Phi_star) / (-t)
    turan_D = A0 + Phi_star*B0 + (Phi_star**2)*C0
    # Now evaluate U36_16, U36_17, U36_18 at t with b = b_star, w36 = 1!
    U16 = (u0 + Phi_star*v0) + vm1
    U17 = (u1 + Phi_star*v1) + v0
    U18 = (u2 + Phi_star*v2) + v1
    turan_U = U17**2 - U16*U18
    # Also check with DOUBLE D_{b*} (as in GLZ!)
    print(f"  t=-1e-{m:2d}: b*={float(b_star):+.3e}, Phi*={float(Phi_star):+.3e}, D_Turan={float(turan_D):+.3e}, U_Turan(w=1)={float(turan_U):+.3e}, |D17/h16|={float(abs(u1+Phi_star*v1)/abs(v0)):.3e}")
