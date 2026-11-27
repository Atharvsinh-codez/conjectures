"""
TEST THE FOLDED PAIRED ABEL IDENTITIES (Orders 0, 1, 2, 4) on BOTH Gaussian and EXACT FINITE-n b_r!
Because U(c-m) = conj(U(m)) and w_0^{c-m} = conj(w_0^m):
  Re[ sum_{m=1}^{c-1} (Delta^p U(m)) w_0^m ] = 2 Re[ sum_{m=1}^{floor((c-1)/2)} (Delta^p U(m)) w_0^m ]
                                            + 1_{c even} Re[ Delta^p U(c/2) (-1)^{rj} ] !
Pairing with 2 Re[ sum_{m>=1} (Delta^p V(m)) w_0^m ] cancels e^{-r*kap*m} IDENTICALLY with ZERO jump at m=c/2!
Let's measure the EXACT constants C_0, C_1, C_2, C_4!
"""
import math

def test_folded_abel_gauss():
    # For Gaussian: U(m) = exp(-gam*m*(1-m/c)) (periodic mod c), V(m) = exp(-gam*|m|*(1+|m|/c)) (even on Z)!
    # Check exact identities AND exact L1 triangle-inequality constants!
    c0_mx = (0.0, None); c1_mx = (0.0, None); c2_mx = (0.0, None); c4_mx = (0.0, None)
    for c in [6, 10, 15, 25, 50, 100]:
        half = (c - 1) // 2
        for gam in [0.005, 0.02, 0.06, 0.15, 0.30, 0.50, 0.80, 1.20, 1.50, 2.50]:
            def U(m):
                mm = m % c
                return math.exp(-gam * mm * (1.0 - mm / c))
            def V(m):
                am = abs(m)
                ex = gam * am * (1.0 + am / c)
                return math.exp(-ex) if ex < 50.0 else 0.0
            Mmax = max(c + 5, int(40.0 / gam) + 5)
            # 1) Order 0: b = 2 sum_{m=1}^{half} (U(m)-V(m)) cos(m th) + 1_{even} U(c/2) cos(c/2 th) - 2 sum_{m>half} V(m) cos(m th)
            L0 = 2.0 * sum(abs(U(m) - V(m)) for m in range(1, half + 1))
            if c % 2 == 0: L0 += abs(U(c // 2))
            L0 += 2.0 * sum(abs(V(m)) for m in range(half + 1, Mmax))
            # Note: L0 only bounded by C_0 / (c * gam^2) when T = c*gam/2 is moderate/large!
            if c * gam >= 4.0:
                v0 = L0 * c * gam * gam
                if v0 > c0_mx[0]: c0_mx = (v0, (c, gam))

            # 2) Order 1: multiply (1 + B) - (1 + A) by (1 - w_0)!
            # Wait: sum_{m=0}^{c-1} U(m) w_0^m * (1 - w_0) = sum_{m=0}^{c-1} (U(m) - U(m-1)) w_0^m
            # Let D1U(m) = U(m) - U(m-1), D1V(m) = V(m) - V(m-1).
            # Notice D1U(c-m+1) = U(c-m+1) - U(c-m) = U(m-1) - U(m) = -D1U(m)!
            # Check identity at a random theta!
            # Even better: SYMMETRIC Order 1! Multiply by (w_0^{-1/2} - w_0^{1/2}) = -2i sin(th/2)!
            #   -2i sin(th/2) * (B - A) = sum_{m=0}^{c-1} (U(m+1) - U(m)) w_0^{m+1/2} - sum_{m=-M}^{M} (V(m+1) - V(m)) w_0^{m+1/2}
            # Fold m and c-1-m! For U: (U(c-m) - U(c-1-m)) = -(U(m+1) - U(m)), and w_0^{c-m-1/2} = conj(w_0^{m+1/2})!
            # So Re[-2i sin(th/2) * (B-A)] = 0 (of course!), and Im[-2i sin(th/2) * (B-A)] = -2 sin(th/2) (B-A)!
            # Thus: 2 sin(th/2) * (B - A) = -2 sum_{m=0}^{half_1} [(U(m+1)-U(m)) - (V(m+1)-V(m))] sin((m+1/2)th) + tail!
            # WAIT A SECOND!! CHECK THAT AT m = 0: U(0) = V(0) = 1, so (U(1)-U(0)) - (V(1)-V(0)) = U(1) - V(1)!!
            # Wait: for even c = 2h: m=0..h-1 pairs with c-1-m = 2h-1..h (EXACT partition of 0..c-1 into h pairs!!)!
            # For odd c = 2h+1: m=0..h-1 pairs with 2h..h+1, and m=h has U(h+1)-U(h) = U(h)-U(h) = 0 IDENTICALLY!!
            # SO FOR BOTH EVEN AND ODD c:
            #   2 sin(th/2) (B - A) = -2 sum_{m=0}^{floor(c/2)-1} [(U(m+1)-U(m)) - (V(m+1)-V(m))] sin((m+1/2)th)
            #                         + 2 sum_{m=floor(c/2)}^{Mmax} (V(m+1) - V(m)) sin((m+1/2)th) !!
            h1 = c // 2
            L1 = 2.0 * sum(abs((U(m + 1) - U(m)) - (V(m + 1) - V(m))) for m in range(0, h1))
            L1 += 2.0 * sum(abs(V(m + 1) - V(m)) for m in range(h1, Mmax))
            # Verify identity at th = 2pi/c!
            if c == 10 and gam == 0.30:
                th = 2.0 * math.pi / c
                lhs = 2.0 * math.sin(0.5 * th) * (sum(U(m)*math.cos(m*th) for m in range(1, c)) - 2.0*sum(V(m)*math.cos(m*th) for m in range(1, Mmax)))
                rhs = (-2.0 * sum(((U(m+1)-U(m)) - (V(m+1)-V(m))) * math.sin((m+0.5)*th) for m in range(0, h1))
                       + 2.0 * sum((V(m+1)-V(m)) * math.sin((m+0.5)*th) for m in range(h1, Mmax)))
                print(f"  Order-1 folded identity check (c=10, gam=0.30): err = {abs(lhs - rhs):.2e}")
            if c == 15 and gam == 0.30:
                th = 2.0 * math.pi / c
                lhs = 2.0 * math.sin(0.5 * th) * (sum(U(m)*math.cos(m*th) for m in range(1, c)) - 2.0*sum(V(m)*math.cos(m*th) for m in range(1, Mmax)))
                rhs = (-2.0 * sum(((U(m+1)-U(m)) - (V(m+1)-V(m))) * math.sin((m+0.5)*th) for m in range(0, h1))
                       + 2.0 * sum((V(m+1)-V(m)) * math.sin((m+0.5)*th) for m in range(h1, Mmax)))
                print(f"  Order-1 folded identity check (c=15, gam=0.30): err = {abs(lhs - rhs):.2e}")
            if gam <= 1.5:
                v1 = (L1 / 2.0) * c * gam
                if v1 > c1_mx[0]: c1_mx = (v1, (c, gam))

            # 3) Order 2 folded:
            #   4 sin^2(th/2) (B - A) = 2(V(1) - U(1)) + 2 sum_{m=1}^{half} (D2V(m) - D2U(m)) cos(m th)
            #                           - 1_{even} D2U(c/2) cos(c/2 th) + 2 sum_{m=half+1}^{Mmax} D2V(m) cos(m th)
            def D2U(m): return U(m + 1) - 2.0 * U(m) + U(m - 1)
            def D2V(m): return V(m + 1) - 2.0 * V(m) + V(m - 1)
            L2 = 2.0 * abs(V(1) - U(1)) + 2.0 * sum(abs(D2V(m) - D2U(m)) for m in range(1, half + 1))
            if c % 2 == 0: L2 += abs(D2U(c // 2))
            L2 += 2.0 * sum(abs(D2V(m)) for m in range(half + 1, Mmax))
            if c == 10 and gam == 0.30:
                th = 2.0 * math.pi / c
                lhs = 4.0 * (math.sin(0.5 * th) ** 2) * (sum(U(m)*math.cos(m*th) for m in range(1, c)) - 2.0*sum(V(m)*math.cos(m*th) for m in range(1, Mmax)))
                rhs = (2.0*(V(1)-U(1)) + 2.0*sum((D2V(m)-D2U(m))*math.cos(m*th) for m in range(1, half+1))
                       - D2U(c//2)*math.cos((c//2)*th) + 2.0*sum(D2V(m)*math.cos(m*th) for m in range(half+1, Mmax)))
                print(f"  Order-2 folded identity check (c=10, gam=0.30): err = {abs(lhs - rhs):.2e}")
            v2 = (L2 / 4.0) * c
            if v2 > c2_mx[0]: c2_mx = (v2, (c, gam))

            # 4) Order 4 folded:
            #   16 sin^4(th/2) (B - A) = 2[(U(2)-V(2)) - 4(U(1)-V(1))] + 2 sum_{m=1}^{half} (D4U(m) - D4V(m)) cos(m th)
            #                            + 1_{even} D4U(c/2) cos(c/2 th) - 2 sum_{m=half+1}^{Mmax} D4V(m) cos(m th)
            def D4U(m): return U(m+2) - 4*U(m+1) + 6*U(m) - 4*U(m-1) + U(m-2)
            def D4V(m): return V(m+2) - 4*V(m+1) + 6*V(m) - 4*V(m-1) + V(m-2)
            L4 = 2.0 * abs((U(2) - V(2)) - 4.0 * (U(1) - V(1))) + 2.0 * sum(abs(D4U(m) - D4V(m)) for m in range(1, half + 1))
            if c % 2 == 0: L4 += abs(D4U(c // 2))
            L4 += 2.0 * sum(abs(D4V(m)) for m in range(half + 1, Mmax))
            if c == 10 and gam == 0.30:
                th = 2.0 * math.pi / c
                lhs = 16.0 * (math.sin(0.5 * th) ** 4) * (sum(U(m)*math.cos(m*th) for m in range(1, c)) - 2.0*sum(V(m)*math.cos(m*th) for m in range(1, Mmax)))
                rhs = (2.0*((U(2)-V(2)) - 4.0*(U(1)-V(1))) + 2.0*sum((D4U(m)-D4V(m))*math.cos(m*th) for m in range(1, half+1))
                       + D4U(c//2)*math.cos((c//2)*th) - 2.0*sum(D4V(m)*math.cos(m*th) for m in range(half+1, Mmax)))
                print(f"  Order-4 folded identity check (c=10, gam=0.30): err = {abs(lhs - rhs):.2e}")
            v4 = (L4 / 16.0) * c / (gam * gam)
            if v4 > c4_mx[0]: c4_mx = (v4, (c, gam))

    print(f"\n  FOLDED PAIRED ABEL CONSTANTS ON GAUSSIAN (compare exact: 8.0, 1.097, 0.758, 1.500):")
    print(f"    Order 0 (T>=2): C0 = L0 * c * gam^2      = {c0_mx[0]:.5f} at {c0_mx[1]}")
    print(f"    Order 1 (g<=1.5): C1 = (L1/2) * c * gam  = {c1_mx[0]:.5f} at {c1_mx[1]}")
    print(f"    Order 2 (all g):  C2 = (L2/4) * c        = {c2_mx[0]:.5f} at {c2_mx[1]}")
    print(f"    Order 4 (all g):  C4 = (L4/16)*c / gam^2 = {c4_mx[0]:.5f} at {c4_mx[1]}")

test_folded_abel_gauss()
