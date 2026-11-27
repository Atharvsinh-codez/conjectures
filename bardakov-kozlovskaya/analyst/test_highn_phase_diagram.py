"""
DIAGNOSTIC FOR COMPLETE EFFECTIVE ANALYTIC PROOF:
1. Check how |avg| / avg|Q| behaves for FIXED small k3 = 1, 2, 3, 5, 10 at HUGE n = 40, 80, 160, 320!
   Wait: earlier at k3=1 obtuse (k1=0.7n), max|arg| grew like n, and |avg|/avg|Q| dropped slowly
   (0.897 -> 0.855 -> 0.816 -> 0.778 -> 0.708). Does |avg|/avg|Q| go to 0 as n -> inf at fixed k3,
   OR does it stay bounded away from 0?!
   AND check what happens if we integrate over psi FIRST (exact [z^{k3}]) vs full 2D torus!
2. Also check the EXACT ratio |avg| / avg|F| for the ELEMENTARY Vandermonde-Cauchy integral
   Q(y, z) * W(y, z) at n = 40, 80, 160 across ALL regimes:
   - bulk (a, b, c proportional to n),
   - small b (b = 1, 2, 5),
   - small c (c = 2, 5, 10, sqrt(n)),
   - small a (a = 1, 2, 5, sqrt(n)) in BOTH Chart I and Chart II!
"""
import numpy as np

def test_Q_2d(n, k1, k2, k3, M=32):
    sgn = (-1)**(k1*k2 + k1*k3 + k2*k3)
    r2 = np.sin(np.pi * k2 / n) / np.sin(np.pi * k1 / n)
    r3 = np.sin(np.pi * k3 / n) / np.sin(np.pi * k1 / n)
    zetas = np.exp(2j * np.pi * np.arange(n) / n)
    g = np.linspace(-np.pi + 0.0123, np.pi + 0.0123, M, endpoint=False)
    P2, P3 = np.meshgrid(g, g, indexing='ij')
    p2f, p3f = P2.ravel(), P3.ravel()
    beta_v = r2 * np.exp(1j * p2f / n)
    w = 1.0 + beta_v[:, None] * zetas[None, :]
    log_w = n * np.log(w)
    log_z = n * np.log(r3) + 1j * p3f[:, None]
    use_z = log_z.real > log_w.real
    log_fac = np.where(use_z,
                       log_z + np.log(1.0 - np.exp(np.where(use_z, log_w - log_z, -100.0))),
                       log_w + 1j*np.pi + np.log(1.0 - np.exp(np.where(~use_z, log_z - log_w, -100.0))))
    log_Q = np.sum(log_fac, axis=1) - k2*(n*np.log(r2) + 1j*p2f) - k3*(n*np.log(r3) + 1j*p3f)
    if sgn < 0: log_Q += 1j * np.pi
    mx = np.max(log_Q.real)
    Q = np.exp(log_Q - mx)
    avg_Q = np.mean(Q)
    abs_Q = np.mean(np.abs(Q))
    # Also compute Vandermonde-Cauchy W(y, z) if |y|+|z| < 0.99 and |z/y|+|z| < 0.99!
    Ry, Rz = r2**n, r3**n
    a, b, c = k1 - k2, k2 - k3, k3
    if Ry + Rz < 0.98 and (Rz / Ry) + Rz < 0.98:
        y = Ry * np.exp(1j * p2f)
        z = Rz * np.exp(1j * p3f)
        W = (1 - y) * (1 - z) * (1 - z/y) / ((1 + y + z)**(a + 1) * (1 + z/y + z)**(b + 1))
        QW = Q * W
        avg_QW = np.mean(QW)
        abs_QW = np.mean(np.abs(QW))
        rat_QW = float(avg_QW.real / abs_QW)
        max_Werr = float(np.max(np.abs(W - 1.0)))
    else:
        rat_QW = float('nan')
        max_Werr = float('nan')
    return float(avg_Q.real / abs_Q), float(np.min(Q.real) / avg_Q.real), rat_QW, max_Werr, float(Ry), float(Rz/Ry)

print("=== 1. Behavior of |avg|/avg|Q| and Vandermonde-Cauchy |avg|/avg|QW| at n = 40, 80, 160, 320 ===")
for n in [40, 80, 160, 320]:
    print(f"\n--- n = {n} ---")
    cases = [
        ("c=1 obtuse",    int(0.7*n), n - int(0.7*n) - 1, 1),
        ("c=2 obtuse",    int(0.7*n), n - int(0.7*n) - 2, 2),
        ("c=5 obtuse",    int(0.7*n), n - int(0.7*n) - 5, 5),
        ("c=sqrt(n) obt", int(0.7*n), n - int(0.7*n) - int(np.sqrt(n)), int(np.sqrt(n))),
        ("c=2 right",     n//2, n//2 - 2, 2),
        ("b=1 bulk",      int(0.4*n), int(0.3*n)+1, int(0.3*n)-1), # adjust below
        ("b=2 bulk",      int(0.4*n), int(0.3*n)+2, int(0.3*n)-2),
        ("a=1 bulk",      int(0.4*n), int(0.4*n)-1, n - 2*int(0.4*n) + 1),
        ("a=2 bulk",      int(0.4*n), int(0.4*n)-2, n - 2*int(0.4*n) + 2),
        ("a=5 bulk",      int(0.4*n), int(0.4*n)-5, n - 2*int(0.4*n) + 5),
        ("bulk interior", n//2, n//3, n - n//2 - n//3),
    ]
    for lbl, k1, k2, k3 in cases:
        if lbl == "b=1 bulk":
            k3 = int(0.25*n); k2 = k3 + 1; k1 = n - k2 - k3
        elif lbl == "b=2 bulk":
            k3 = int(0.25*n); k2 = k3 + 2; k1 = n - k2 - k3
        if not (k1 >= k2 >= k3 >= 1 and k1 + k2 + k3 == n):
            continue
        rQ, minQ, rQW, werr, Ry, Rzy = test_Q_2d(n, k1, k2, k3, M=24)
        print(f"  {lbl:14s} (a={k1-k2:3d},b={k2-k3:3d},c={k3:3d}): Q_rat={rQ:+.4f} (minRe={minQ:+.4f}) | QW_rat={rQW:+.4f} (max|W-1|={werr:.3e}, Ry={Ry:.2e}, Rz/Ry={Rzy:.2e})")
