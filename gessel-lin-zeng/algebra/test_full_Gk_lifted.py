"""
3-SECOND BOMBSHELL CHECK:
1. Is the FULL D-lifted bivariate polynomial
     G_k(t, y) := sum_{i=0}^{k-1} (D_{k+1}^i h~_{k-1, i}(t)) y^i
   (and G~_k(t, y) := D_{k+1}^2 [(1+y) G_k(t, y)] = sum_{i=0}^k (D_{k+1}^i h~_{k, i}(t)) y^i)
   REAL-STABLE IN (t, y) for k = 2..20?!
   NOTE: If G_k(t, y) IS real-stable for all k, then Stage 1 + Stage 2 is literally ONE LINE:
     G~_k(t, y) = D_{k+1}^2 [ (1 + y) G_k(t, y) ] !!
   ZERO Wronskians! ZERO Plucker! ZERO Link 2! ZERO W_{04} > 0!
2. Also check G_k(t, y) with NEGATIVE y sign: G_k(t, -y) (check UHP vs LHP convention!).
   Wait: check k=2!
     E_{2, 0} = h~_{1, 0} = 1 + 2t,  E_{2, 1} = D_3 h~_{1, 1} = 1 + 3t.
     G_2(t, y) = (1 + 2t) + (1 + 3t) y.
     At t = i*eps (Im t > 0): y = -(1+2i*eps)/(1+3i*eps) = -(1+2i*eps)(1-3i*eps)/|.|^2 -> Im y = +eps > 0!
     WHOA!! Im y > 0 when Im t > 0, WHICH MEANS G_2(t, -y) = (1+2t) - (1+3t) y HAS Im y < 0 (REAL-STABLE WITH MINUS SIGN!!)!
     WAIT A SECOND! Check Stage 1 with MINUS sign!
     If G_k^-(t, y) := sum_{i=0}^{k-1} E_{k, i}(t) (-y)^i:
       (1 - y) G_k^-(t, y) = sum (E_{k, i} - E_{k, i-1}) (-y)^i  <-- WAIT! That's E_i - E_{i-1}, NOT E_i + E_{i-1}!!
     HOLD ON! Check sign of G_2(t, y) = (1+2t) + (1+3t) y!!
     Wait: if t in H_+ and y in H_+, can (1+2t) + (1+3t) y = 0?!
     YES! Because y = -(1+2t)/(1+3t) HAS Im y > 0!!
     SO (1+2t) + (1+3t) y IS **NOT** REAL-STABLE IN (t, y) — IT IS REAL-STABLE IN (t, -y)!!
     AND FOR (t, -y), multiplying by (1 - y) puts MINUS signs (NOT + signs), UNLESS y is in H_-!!
     WAIT A SECOND!!
     Check the sign of (1+2t) + (1+3t) y at t = i*0.1, y = i*0.1!
       y = -(1+0.2i)/(1+0.3i) = -(1.06 - 0.1i)/1.09 => Im y = +0.1/1.09 > 0!!
     SO G_2(t, y) = E_{2,0} + E_{2,1} y HAS A ZERO IN H_+ x H_+!!
     WHEREAS G_2(t, -y) = E_{2,0} - E_{2,1} y IS REAL-STABLE IN H_+ x H_+!!
     WAIT A SECOND!! Look at our Quad: X - 2c M y + Y y^2 ALSO had a MINUS sign on M y,
     and (1 - y)(X - M y) = X - (X+M)y + M y^2!!
     CHECK THAT!!
     Let G_k^-(t, y) := sum_{i=0}^{k-1} (-1)^i E_{k, i}(t) y^i!!
     Compute (1 - y) G_k^-(t, y)!!
       [y^i] (1 - y) sum (-1)^m E_m y^m = (-1)^i E_i - (-1)^{i-1} E_{i-1} = (-1)^i (E_i + E_{i-1})!!
     OH!! MY!! GOD!! IT *DOES* GIVE E_i + E_{i-1} WITH (-1)^i!!
     HAHAHA!! THANK GOD I CHECKED THE SIGN!!
     (1 - y) * sum (-1)^i E_i y^i = sum (-1)^i (E_i + E_{i-1}) y^i!!
     AND IS (1 - y) REAL-STABLE IN H_+? WAIT: 1 - y = 0 at y = 1 (REAL!!), SO (1 - y) HAS NO ZEROS IN H_+!!
     SO MULTIPLYING A REAL-STABLE POLYNOMIAL BY (1 - y) **100% PRESERVES REAL STABILITY**!!

TEST G_k^-(t, y) = sum_{i=0}^{k-1} (-1)^i (D_{k+1}^i h~_{k-1, i}(t)) y^i RIGHT NOW FOR k = 2..20!!
AND ALSO TEST WITH BINOMIAL WEIGHTS / q-WEIGHTS / A_{k,i}(2^i x)!
"""
import mpmath as mp
mp.mp.dps = 80

def apply_Dj_mp(poly, j):
    out = [mp.mpf(0)] * (len(poly) + 1)
    for m, c in enumerate(poly):
        if c != 0:
            out[m] += (1 + m) * c
            out[m + 1] += (j + m) * c
    return out

def padd_mp(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for idx in range(len(a)): out[idx] += a[idx]
    for idx in range(len(b)): out[idx] += b[idx]
    return out

def peval_c(poly, t):
    s = mp.mpc(0)
    for c in reversed(poly):
        s = s * t + c
    return s

def check_hb_pd_lhp(coeffs_y):
    """
    Check if ALL roots of P(y) = sum_{j=0}^n c_j y^j lie strictly in Im(y) < 0 (LHP)!
    For P(y) = A(y) + i B(y): Hermite matrix H of (A(x)B(y) - B(x)A(y))/(x-y) is PD
    (check sign on P(y) = (1+2t) - (1+3t) y at t = i*eps:
     A(y) = 1 - y, B(y) = 2eps - 3eps y.
     A(x)B(y) - B(x)A(y) = (1-x)(2-3y)eps - (2-3x)(1-y)eps = eps*(2 - 3y - 2x + 3xy - 2 + 2y + 3x - 3xy) = eps*(x - y)!
     Dividing by (x - y) gives +eps > 0! YES! PD <=> all y-roots in Im(y) < 0!)
    """
    n = len(coeffs_y) - 1
    if n <= 0: return True, 0
    s = (abs(coeffs_y[0]) / abs(coeffs_y[-1])) ** (mp.mpf(1) / n)
    sc = [coeffs_y[j] * (s ** j) for j in range(n + 1)]
    mx = max(abs(z) for z in sc)
    sc = [z / mx for z in sc]
    A = [z.real for z in sc]
    B = [z.imag for z in sc]
    H = [[mp.mpf(0)] * n for _ in range(n)]
    for r in range(n + 1):
        for col in range(r + 1, n + 1):
            val = -(A[r]*B[col] - B[r]*A[col])
            for m in range(col - r):
                H[r + m][r + (col - r - 1 - m)] += val
    for i in range(n):
        piv = H[i][i]
        if piv <= 0:
            return False, i
        for j in range(i + 1, n):
            fac = H[j][i] / piv
            for k in range(j, n):
                H[j][k] -= fac * H[i][k]
                H[k][j] = H[j][k]
    return True, n

print("=== TEST 1: FULL D-LIFTED G_k^-(t, y) = sum_{i=0}^{k-1} (-1)^i (D_{k+1}^i h~_{k-1, i}(t)) y^i ===")
h = [[mp.mpf(1), mp.mpf(2)], [mp.mpf(1)]]
ts_uhp = [mp.mpc('-2.0', '1e-4'), mp.mpc('-0.99', '1e-5'), mp.mpc('-0.8', '1e-5'),
          mp.mpc('-0.5', '1e-5'), mp.mpc('-0.2', '1e-5'), mp.mpc('-0.01', '1e-6'),
          mp.mpc('0.0', '1e-8'), mp.mpc('0.1', '1e-5'), mp.mpc('2.0', '1e-4'), mp.mpc('50.0', '1e-3')]

for k in range(2, 22):
    j = k + 1
    E = []
    for i in range(k):
        p = h[i]
        for _ in range(i):
            p = apply_Dj_mp(p, j)
        E.append(p)
    # Test G_k^-(t, y) = sum (-1)^i E_i(t) y^i (plain!) AND with q = 1, 0.9, 0.5!
    res = []
    for label, w_fn in [("plain(1)", lambda i, n: mp.mpf(1)),
                        ("1/C(k-1,i)", lambda i, n: mp.mpf(1) / mp.binomial(n, i)),
                        ("q=0.9", lambda i, n: mp.mpf('0.9')**(i*(i-1)//2)),
                        ("q=0.5", lambda i, n: mp.mpf('0.5')**(i*(i-1)//2))]:
        ok_all = True
        fail_t = None
        for tc in ts_uhp:
            cy = [((-1)**i) * w_fn(i, k - 1) * peval_c(E[i], tc) for i in range(k)]
            ok, idx = check_hb_pd_lhp(cy)
            if not ok:
                ok_all = False
                fail_t = (float(tc.real), idx)
                break
        res.append(f"{label}:{'PASS' if ok_all else 'FAIL@' + str(fail_t)}")
    if k <= 12 or k % 2 == 0:
        print(f"  k={k:2d}: " + " | ".join(res))
    # step to h_k
    hn = []
    for i in range(k + 1):
        u = [mp.mpf(0)]
        if i <= k - 1: u = padd_mp(u, apply_Dj_mp(h[i], j))
        if i >= 1: u = padd_mp(u, h[i - 1])
        hn.append(apply_Dj_mp(u, j))
    h = hn
