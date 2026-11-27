"""Check coefficient ratios [x^j] A_{k, k-r}(x) / [x^j] A_{k+r, k+r}(x) AND
check whether A_{k, k-r}(x) is a differential operator in (x, theta_x) applied to A_{k+r, k+r}(x)!!"""
import sympy as sp

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

A = build_A(12)

print("=== 1. RATIO [x^j] A_{k, k-r} / [x^j] A_{k+r, k+r} FOR r = 1, 2, 3 ===")
for r in [1, 2, 3]:
    print(f"\n--- r = {r} ---")
    for k in range(r, r + 4):
        p1 = A[k][k - r]
        p2 = A[k + r][k + r]
        assert len(p1) == len(p2)
        ratios = [sp.Rational(p1[j], p2[j]) for j in range(len(p1))]
        # Also check finite differences of ratios to see if it's a polynomial in j!
        d1 = [ratios[j+1] - ratios[j] for j in range(len(ratios)-1)]
        d2 = [d1[j+1] - d1[j] for j in range(len(d1)-1)]
        print(f"  k={k}: ratios = {ratios[:6]}")
        if r == 2:
            print(f"       2nd diffs = {d2[:5]}")

print("\n=== 2. CHECK IF A_{k, k-2}(x) = (c0 + c1*x*d/dx + c2*(x*d/dx)^2 + d1*x + d2*x^2*d/dx + ...) A_{k+2, k+2}(x) ===")
# Test if A_{k, k-2} is in the span of { x^a (d/dx)^b A_{k+2, k+2} } or { A_{k+2,k+2}, A'_{k+2,k+2}, A''_{k+2,k+2} } with poly coeffs in x!
x = sp.Symbol('x')
for k in [3, 4, 5, 6]:
    target = sum(c * x**j for j, c in enumerate(A[k][k - 2]))
    base = sum(c * x**j for j, c in enumerate(A[k + 2][k + 2]))
    b0 = base
    b1 = sp.diff(base, x)
    b2 = sp.diff(b1, x)
    # Try c0*b0 + (c10 + c11*x)*b1 + (c20 + c21*x + c22*x^2)*b2 = target
    # Wait: deg(target) = deg(base) = k+1!
    # So only (a0)*b0 + (a1*x)*b1 + (a2*x^2)*b2 + (b10)*b1 + (b20 + b21*x)*b2 keep deg <= k+1!
    syms = sp.symbols('a0 a1 a2 b10 b20 b21')
    a0, a1, a2, b10, b20, b21 = syms
    expr = sp.expand(a0*b0 + (b10 + a1*x)*b1 + (b20 + b21*x + a2*x**2)*b2 - target)
    eqs = sp.Poly(expr, x).all_coeffs()
    sol = sp.solve(eqs, syms)
    print(f"  k={k} (r=2): diff-op on A_{{k+2,k+2}}: sol = {sol}")

print("\n=== 3. ALSO CHECK DIFF-OP RELATING A_{k, i} TO A_{k+1, i+1} (SHIFT BY 1 IN BOTH k AND i!) ===")
# Note deg A_{k, i} = 2k-i-1, deg A_{k+1, i+1} = 2k-i = deg + 1!
# Wait: at i = k-1, A_{k, k-1} = (k - x/2 d/dx) A_{k+1, k+1} had SAME degree!
# In general: deg A_{k, i} = 2k-i-1, and deg A_{k+1, i+2} = 2(k+1)-(i+2)-1 = 2k-i-1 (EXACT SAME DEGREE!!)!
# OH!! CHECK A_{k, i} vs A_{k+1, i+2}!!
# Notice at i = k-1: A_{k+1, k+1} IS A_{k+1, (k-1)+2}!!
print("--- CHECKING RATIO [x^j] A_{k, i} / [x^j] A_{k+1, i+2} (EXACT SAME DEGREE FOR ALL k, i!!) ---")
for k in [3, 4, 5, 6]:
    for i in range(0, k):
        if i + 2 > k + 1: continue
        p1 = A[k][i]
        p2 = A[k + 1][i + 2]
        assert len(p1) == len(p2)
        ratios = [sp.Rational(p1[j], p2[j]) for j in range(len(p1))]
        d1 = [ratios[j+1] - ratios[j] for j in range(len(ratios)-1)]
        print(f"  k={k}, i={i}: A_{{k,i}}/A_{{k+1,i+2}} ratios={ratios[:5]}, d1={d1[:4]}")
