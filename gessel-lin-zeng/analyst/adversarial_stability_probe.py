"""Critical diagnostic script:
1. Test whether plain sum_{i=0}^k A_{k,i}(x) y^i is real-rooted in y for ALL x < 0 (k=2..14)!
2. Test B_k(x,y) = sum_{i=0}^k binom(k,i) A_{k,i}(x) y^i for real stability VERY adversarially
   (near real roots, high k up to 14, negative y, positive lines x=a+s, y=b+c*s with c>0)!
3. Check the EXACT Wronskian / Rayleigh determinant of B_k(x,y):
   Delta_B(x,y) := (d_x B)(d_y B) - B (d_x d_y B)  (must be >= 0 for ALL real (x,y) if real-stable!).
   Wait: check sign! For P(x,y) = (x+1)(y+1) - eps = xy + x + y + 1 - eps:
   P_x P_y - P P_{xy} = (y+1)(x+1) - (xy+x+y+1-eps)(1) = +eps. Wait: for real-stable P,
   as x increases, roots y(x) MUST DECREASE (dy/dx = -P_x/P_y <= 0), wait:
   If P(x,y) = x + y: P_x P_y - P P_{xy} = 1 > 0!
   If P(x,y) = (x+a)(y+b) - c: P_x P_y - P P_{xy} = +c. At y = -b + c/(x+a), dy/dx = -c/(x+a)^2 < 0 iff c > 0!
   YES!! By Branden (2007), for a bivariate polynomial P(x,y), real stability =>
   P_x(x,y) P_y(x,y) - P(x,y) P_{xy}(x,y) >= 0 for ALL real (x,y)!! (Wait: check sign when deg_y = 1:
   P(x,y) = f(x) + y g(x). P_x P_y - P P_{xy} = f' g - f g' = W(f, g). If g << f (f outer-interlaces g,
   roots of f < roots of g? Wait: at largest root r of g, f(r) < 0, g'(r) > 0 => -f(r) g'(r) > 0!
   Check f(x) = (x+1)(x+3) = x^2+4x+3, g(x) = x+2: f'g - fg' = (2x+4)(x+2) - (x^2+4x+3)(1) = x^2+4x+5 > 0!
   YES!! P_x P_y - P P_{xy} >= 0 EVERYWHERE ON R^2 is NECESSARY (and for multiaffine, SUFFICIENT) for real stability!)
   Let's check P_x P_y - P P_{xy} for B_k(x, y) AND for plain A_k(x, y) on R^2!
"""
import math
import numpy as np

def poly_add(a, b):
    n = max(len(a), len(b))
    res = [0] * n
    for i, v in enumerate(a): res[i] += v
    for i, v in enumerate(b): res[i] += v
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def apply_S(L, p):
    d = len(p) - 1
    res = [0] * (d + 2)
    for m, c in enumerate(p):
        res[m] += c * (m + 1)
        res[m + 1] += c * (L - m)
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def build_A(kmax):
    A = [[None] * (k + 1) for k in range(kmax + 1)]
    A[0][0] = [1]
    for k in range(1, kmax + 1):
        A[k][0] = apply_S(3*k - 2, apply_S(3*k - 3, A[k-1][0]))
        for i in range(1, k):
            L = 3*k - 2 - i
            A[k][i] = apply_S(L, poly_add(apply_S(L - 1, A[k-1][i]), A[k-1][i-1]))
        A[k][k] = apply_S(2*k - 2, A[k-1][k-1])
    return A

A = build_A(14)

print("=== 1. TEST PLAIN sum_i A_{k,i}(x) y^i IN y FOR x < 0 (k=2..14) ===")
for k in range(2, 15):
    worst = (0.0, None)
    xs = [-1000, -100, -20, -10, -5, -3, -2, -1.5, -1.0, -0.7, -0.5, -0.3, -0.2, -0.1, -0.05, -0.01, -1e-4]
    xs += list(np.linspace(-5.0, -0.001, 500))
    for x0 in xs:
        ycoeffs = [sum(float(c) * (x0**m) for m, c in enumerate(A[k][i])) for i in range(k + 1)]
        if abs(ycoeffs[-1]) < 1e-30: continue
        r = np.roots(list(reversed(ycoeffs)))
        mi = float(np.max(np.abs(r.imag)))
        if mi > worst[0]:
            worst = (mi, x0)
    print(f"  k={k:2d} plain A_k(x0, y) on x0<0: max |Im(y)| = {worst[0]:.3e} (at x0={worst[1]})")

print("\n=== 2. ADVERSARIAL RAYLEIGH TEST P_x P_y - P P_{xy} FOR B_k(x,y) ON R^2 (k=2..12) ===")
def eval_poly_and_deriv(coeffs, x0):
    val = 0.0
    dval = 0.0
    for c in reversed(coeffs):
        dval = dval * x0 + val
        val = val * x0 + float(c)
    return val, dval

for k in range(2, 13):
    w = [math.comb(k, i) for i in range(k + 1)]
    min_val = (1e300, None)
    # Test grid of (x, y) in R^2, focusing on x < 0 (where roots live) and x > 0!
    xs = list(np.linspace(-10, 10, 101)) + list(np.linspace(-2, 0, 81))
    ys = list(np.linspace(-100, 100, 101)) + list(np.linspace(-10, 10, 81)) + [-1000, 1000]
    for x0 in xs:
        vals = [eval_poly_and_deriv(A[k][i], x0) for i in range(k + 1)]
        # normalize to avoid overflow
        scale = max(abs(v[0]) for v in vals) + 1e-300
        f = [w[i] * vals[i][0] / scale for i in range(k + 1)]
        fx = [w[i] * vals[i][1] / scale for i in range(k + 1)]
        for y0 in ys:
            # P, Px, Py, Pxy
            P = Px = Py = Pxy = 0.0
            ypow = 1.0
            for i in range(k + 1):
                P += f[i] * ypow
                Px += fx[i] * ypow
                if i >= 1:
                    Py += i * f[i] * (ypow / y0 if y0 != 0 else (1.0 if i == 1 else 0.0))
                    Pxy += i * fx[i] * (ypow / y0 if y0 != 0 else (1.0 if i == 1 else 0.0))
                ypow *= y0
            # Wait: avoid 0/0 when y0==0
            if y0 == 0:
                P = f[0]; Px = fx[0]; Py = f[1]; Pxy = fx[1]
            ray = Px * Py - P * Pxy
            # normalize by |Px*Py| + |P*Pxy| + 1e-300
            denom = abs(Px * Py) + abs(P * Pxy) + 1e-300
            rel = ray / denom
            if rel < min_val[0]:
                min_val = (rel, (x0, y0, ray))
    print(f"  B_{k:2d}(x,y): min normalized (Px*Py - P*Pxy) = {min_val[0]:+.4e} at x={min_val[1][0]:.3f}, y={min_val[1][1]:.3f}")

print("\n=== 3. ADVERSARIAL POSITIVE-LINE TEST FOR B_k(x, y): x = a + s, y = b + c*s (c > 0) (k=2..10) ===")
# Evaluate B_k(a + s, b + c*s) as poly in s and check all roots are real!
rng = np.random.default_rng(777)
for k in range(2, 11):
    w = [math.comb(k, i) for i in range(k + 1)]
    max_im = 0.0
    worst = None
    for trial in range(200):
        a = rng.uniform(-5, 5)
        b = rng.uniform(-50, 50)
        c = 10 ** rng.uniform(-3, 3)
        # Build poly in s via numpy.polynomial.polynomial
        from numpy.polynomial import polynomial as Pnp
        p_s = np.array([0.0])
        x_lin = np.array([a, 1.0])
        y_lin = np.array([b, c])
        y_pow = np.array([1.0])
        for i in range(k + 1):
            # eval A[k][i](a + s)
            ai_s = Pnp.polyval(x_lin, [float(v) for v in A[k][i]]) # wait polyval doesn't compose!
            # Use Horner with polymul
            acc = np.array([0.0])
            for coeff in reversed(A[k][i]):
                acc = Pnp.polyadd(Pnp.polymul(acc, x_lin), np.array([float(coeff)]))
            term = Pnp.polymul(acc * w[i], y_pow)
            p_s = Pnp.polyadd(p_s, term)
            y_pow = Pnp.polymul(y_pow, y_lin)
        # Normalize and find roots
        p_s = p_s / np.max(np.abs(p_s))
        r = np.roots(p_s[::-1])
        mi = float(np.max(np.abs(r.imag)))
        if mi > max_im:
            max_im = mi
            worst = (a, b, c)
    print(f"  B_{k:2d}(a+s, b+c*s): max |Im(s)| over 200 random positive lines = {max_im:.3e} (worst={worst})")
