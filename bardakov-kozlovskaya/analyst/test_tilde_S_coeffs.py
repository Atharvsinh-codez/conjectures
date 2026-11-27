"""
Check the exact coefficients tilde_S_c(n, T) = [v^T] (-1)^c E_{n-c}(1, v) / (1 - v)^n
for n = 12, 20, 30 and c = 2, 3, 4, 5!
"""
import sys
sys.path.insert(0, r"d:\solving-conjectures\work\bardakov-kozlovskaya\verify")
import engine
import numpy as np
from math import comb

for n in [12, 20, 30]:
    E = engine.compute_E_table(n)
    print(f"\n================ n = {n} ================")
    for c in [2, 3, 4, 5]:
        # E_k(1, v) = E[k] (check sign & length!)
        poly = E[n - c] # coefficient of z^c in p_n(1, v, z) up to (-1)^{n-c} = (-1)^c
        # Check poly[0]: should be comb(n, c) (or (-1)^c comb(n, c))
        sgn_fix = 1 if poly[0] > 0 else -1
        # Compute tilde_S_c(n, T) = sum_{m=0}^{min(T, len(poly)-1)} (sgn_fix * poly[m]) * comb(n + T - m - 1, T - m)
        T_max = n // 2 + 5
        tS = []
        for T in range(T_max + 1):
            val = 0
            for m in range(min(T + 1, len(poly))):
                val += (sgn_fix * poly[m]) * comb(n + T - m - 1, T - m)
            tS.append(val)
        signs = "".join("+" if v > 0 else ("-" if v < 0 else "0") for v in tS[:min(22, T_max+1)])
        print(f"  c={c}: poly[0]={poly[0]} (comb={comb(n,c)}) | signs(T=0..21): {signs}")
        # Check dominance of T=j term at r_j for j = c+1 .. (n-c)//2:
        for j in range(c + 1, (n - c) // 2 + 1):
            rj = np.sin(j * np.pi / n) / np.sin((j + c) * np.pi / n)
            uj = ((-1)**c) * (rj ** n)
            terms = [tS[T] * (uj ** T) for T in range(T_max + 1)]
            tot = sum(terms)
            t_j = terms[j]
            ratio_dom = abs(t_j) / sum(abs(x) for x in terms)
            exp_sgn = 1 if (j % 2 == 0) else -1
            if j in (c + 1, c + 2, (n - c) // 2):
                kap = 2*np.pi*np.sin(j*np.pi/n)*np.sin((j+c)*np.pi/n)/np.sin(c*np.pi/n)
                print(f"     j={j:2d} (kap={kap:5.2f}): sgn(t_j)={'+' if t_j>0 else '-'} sgn(tot)={'+' if tot>0 else '-'} exp={'+' if exp_sgn>0 else '-'} | |t_j|/sum|t_T|={ratio_dom:.6f} | tot/t_j={tot/t_j:.6f}")
