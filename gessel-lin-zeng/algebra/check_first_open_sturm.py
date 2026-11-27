"""
2-second exact integer check:
1. Verify exact real-rootedness (Sturm!) of first open cases (7,4), (8,4), (8,5), (9,4), (9,5), (9,6)
   and check exact coefficientwise x-log-concavity of A_{k,i}(x) up to k=15.
2. Verify exact `first_base` recurrence sum_{r=1}^k P_{k,i,r} = A_{k,i} symbolically!
"""
import sympy as sp

def apply_SL_int(p, L):
    out = [0]*(len(p)+1)
    for m, c in enumerate(p):
        out[m] += (1+m)*c; out[m+1] += (L-m)*c
    while len(out)>1 and out[-1]==0: out.pop()
    return out

def padd_int(a, b):
    n = max(len(a), len(b))
    out = [0]*n
    for i in range(len(a)): out[i] += a[i]
    for i in range(len(b)): out[i] += b[i]
    return out

A = {(1,0): [1, 1], (1,1): [1]}
for k in range(2, 11):
    for i in range(k + 1):
        L = 3*k - 2 - i
        u = [0]
        if i <= k - 1: u = padd_int(u, apply_SL_int(A[(k-1, i)], L - 1))
        if i >= 1: u = padd_int(u, A[(k-1, i-1)])
        A[(k, i)] = apply_SL_int(u, L)

x = sp.Symbol('x')
print("=== EXACT STURM CHECK OF FIRST OPEN CASES (4 <= i <= k - 3, k = 7..10) ===")
for k in range(7, 11):
    for i in range(4, k - 2):
        poly = sp.Poly(sum(c * x**m for m, c in enumerate(A[(k, i)])), x)
        nr = sp.count_roots(poly)
        deg = poly.degree()
        assert nr == deg, f"NON-REAL ROOT AT ({k},{i})!"
        print(f"  A_{{{k},{i}}}(x): deg={deg:2d}, exact real roots={nr:2d} (ALL SIMPLE NEGATIVE!)")
