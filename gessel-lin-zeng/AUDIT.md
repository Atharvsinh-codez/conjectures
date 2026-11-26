# Independent Verification & Adversarial Audit Report: Gessel–Lin–Zeng Conjecture 15

**Auditor:** `verifier` (2026-11-26)  
**Target:** Gessel–Lin–Zeng Conjecture 15 (Issue `#84`, arXiv:1201.0622; SymCat live open problem for $4 \le i \le k-3$).  
**Verification scripts:** [`verify/testbed.py`](verify/testbed.py), [`verify/probes.py`](verify/probes.py), [`verify/verify_msg13.py`](verify/verify_msg13.py), [`verify/verify_deep_positivity.py`](verify/verify_deep_positivity.py), [`verify/audit_18_20_21.py`](verify/audit_18_20_21.py), [`verify/audit_25.py`](verify/audit_25.py), [`verify/audit_31_32_33.py`](verify/audit_31_32_33.py).  
**Status:** Conjecture 15 remains **open** (archived prior to swarm pivot to Bardakov–Kozlovskaya Conjecture 3.3(2)). Every theorem, sharp barrier, and disproof listed below has been independently verified from scratch using 100% exact Python integers, `Fraction` arithmetic, exact Sturm/Cauchy-index chains, and SymPy.

---

## 1. Exact Definition & Recurrence Cross-Verification

- **Direct Jacobi–Stirling permutation enumeration** (`verify/testbed.py`):
  - Brute-force backtracking over $M_{k,S} = \{\bar 1,1,1,\dots,\bar k,k,k\}\setminus\{\bar j : j\in S\}$ counting internal descents matches the insertion operator $S_L f(x) = (1+Lx)f(x) + x(1-x)f'(x)$ identically for **all** subsets $S\subseteq[k]$ at $k=1,2,3,4$ (including all 246,400 words at $k=4, S=\varnothing$) and all 16 subsets for $k=5, i\in\{3,4,5\}$.
- **Möbius-transformed $i$-independent recurrence:**
  - Verified identically for all $1\le k\le 24, 0\le i\le k$:
    $$h_{k,i}(t) = (1+t)^{2k-i} A_{k,i}\!\left(\frac{t}{1+t}\right) = (1+t)\tilde h_{k,i}(t), \qquad \tilde h_{k,i} = D_{k+1}\bigl(D_{k+1}\tilde h_{k-1,i} + \tilde h_{k-1,i-1}\bigr),$$
    with $D_j = t(1+t)\partial_t + 1 + jt$, $\tilde h_{1,0} = 1+2t$, $\tilde h_{1,1} = 1$, and $\deg\tilde h_{k,i} = 2k-i-1$.
- **Exact Sturm & Cauchy-index baseline certificates** (`verify/testbed.py`, `verify/probes.py`, `verify/audit_18_20_21.py`, `verify/audit_25.py`):
  - Simple negative real-rootedness of $A_{k,i}(x)$ ($1\le k\le 24$): **100% PASS**.
  - Aggregate row interlacing $\tilde h_{k,i}\prec \tilde h_{k,i-1}$ and inner-sum row $U_{k,i}\prec U_{k,i-1}$ ($k\le 24$): **299/299 PASS** each.
  - Inner-sum interlacing $\tilde h_{k-1,i-1}\prec D_{k+1}\tilde h_{k-1,i}$ ($k\le 24$): **276/276 PASS**.
  - 4th cross-pair $D_{k+1}\tilde h_{k-1,i}\prec \tilde h_{k-1,i-2}$ ($k\le 20$): **171/171 PASS**; $\tilde h_{k-1,i-2}\prec D_{k+1}^2\tilde h_{k-1,i}$ ($k\le 16$): **105/105 PASS**.
  - Operator intertwining $P_{L+1}S_{L+1}=S_L P_{L+1}$ ($P_c = c - x\partial_x$) verified symbolically; $A_{k,i}\prec P_{3k-1-i}A_{k,i-1}$ (**253/253 PASS**, $k\le 22$), $S_{3k-3-i}A_{k-1,i}\prec P_{3k-1-i}A_{k-1,i-2}$ (**210/210 PASS**, $k\le 22$), and **Derivative Sandwich** $A'_{k,i-1}\prec A_{k,i}\prec P_{2k-i}A_{k,i-1}$ (**252/252 PASS**, $2\le k\le 22$).

---

## 2. Independently Verified General Theorems (With Audited Hypotheses)

Throughout, $W(f,g) := fg' - f'g$ and $J(f,g,h) := \frac{4W(f,g)W(g,h)}{W(f,h)^2}$.

1. **Theorem A (Turán $\implies$ Interlacing Preservation):** VERIFIED (`audit_18_20_21.py`). For $p_2\prec p_1\prec p_0$ in $(-1,0)$ with **degree-increasing** $\deg p_0 = \deg p_1 + 1 = \deg p_2 + 2$, positive leading coefficients, and $p_1(t)^2\ge p_0(t)p_2(t)$ on $(-1,0)$, we have $Tp_2 + wp_1 \prec Tp_1 + wp_0$ for any $T = t(1+t)\partial_t + a + bt$ ($0<a<b$) and $w>0$, via $U_1(r)p_1(r) = (-r)(1+r)W(p_2,p_1)(r) + w(p_1(r)^2 - p_0(r)p_2(r)) > 0$. Constant $1$ is sharp.
2. **Single-Operator Turán & Exact $T$-Turán Equivalence:** VERIFIED. $(Tp)^2 - p(T^2p) = (-t)(1+t)W(p,Tp) > 0$ on $(-1,0)$; $(Tp_i)^2 - (Tp_{i-1})(Tp_{i+1}) = A\Phi^2 + B\Phi + C$ with $B = t(1+t)A'$ and $4AC - B^2 = t^2(1+t)^2[4W(p_{i-1},p_i)W(p_i,p_{i+1}) - W(p_{i-1},p_{i+1})^2]$ (at $A=0$, requires $(p_{i-1}(t),p_i(t),p_{i+1}(t))\ne(0,0,0)$, automatic under strict interlacing).
3. **Sharpest Quantitative Plücker Quad Theorem B:** VERIFIED (`audit_18_20_21.py`, `audit_25.py`).
   - *Audit defects resolved:* Original statement lacked $W_{AC},W_{BE}\ge 0$ (counterexample $A=(-9/5,-1), B=(1,0), C=(0,1), E=(-1,0)$ gave $J_{\text{new}}=40/81<1$) and suffered sign flip of $W(A+B,C+E)$ for $R>1$.
   - *Verified airtight form:* Under $W_{AB},W_{BC},W_{CE},W_{AC},W_{BE}>0$ and $W_{AE}\ge 0$, Plücker gives $uv\ge ab \implies W(A+B,C+E)/W_{BC} = 1+u+v+(uv-ab) > 0$ unconditionally, forces $s=(ab)^{1/4}\le 2\theta$ (zero $R\le 4$ hypothesis needed!), and yields $J(A+B,B+C,C+E) \ge \frac{4(1+\epsilon)^2}{(2+\epsilon)^2}$ for all $\epsilon\ge 0$.
4. **Zero-Loss $D_m$ Preservation:** VERIFIED (`audit_18_20_21.py`).
   - *Audit fixes incorporated:* Requires $m\ge 1$, $X,M,Y\in\mathrm{RZ}$ with $M\prec Y$, $W_{XM},W_{MY}>0$, and **minus** sign in $y$: $X(t) - \frac{2}{\sqrt{1+\epsilon}}M(t)y + Y(t)y^2$ is real-stable (proved via Hermite–Biehler $W_t(\operatorname{Re}P,\operatorname{Im}P) = -2\eta[\theta W_{XM} - yW_{XY} + \theta(y^2+\eta^2)W_{MY}] < 0$; exact Borcea–Brändén symbol partial fraction $w(t) = -\frac{d+m}{m}t - \frac{d(m-1)}{m^2} + \frac{d(m-1)/m^2}{mt+1}$ proves $D_m$ preserves real stability iff $m\ge 1$).
5. **Universal Cubic Wronskian Log-Concavity & Laguerre Polar Derivative:** VERIFIED (`audit_18_20_21.py`). Exact global SOS proves $W(p_0,p_1)W(p_2,p_3)\le W(p_1,p_2)^2$ on $\mathbb{R}$ for every real-stable cubic $P(t,y)=\sum_{j=0}^3 p_j(t)y^j$ (yielding constant $9$ in binomial form $A+3By+3Cy^2+Ey^3$).
6. **Universal Cubic Real-Stability Addition:** VERIFIED (`audit_25.py`). For 5 vectors $v_0,\dots,v_4\in\mathbb{R}^2$ in strict open-half-plane order ($W(v_r,v_s)>0$ for $0\le r<s\le 4$), PSD Bezoutiants of $P_1, P_2$ imply PSD Bezoutiant of $P_1+P_2$ (all sign patterns, degenerate $L_r=0$ via density/Hurwitz, double-root $m_1=m_2$, and PSD orientation lock $W(v_0+v_1,v_1+v_2)>0$ verified; sharpness witness $L=(-3000,1,10,1,-3000)$ certified by exact integer Sturm; 3,000/3,000 random stress tests passed).
7. **Theorem C (Pole-Gap Wronskian Bound):** VERIFIED. For $\deg f_0 = \deg f_1 + 1$, $f_1\prec f_0$ in $(-1,0)$, and $f_0\prec D_r f_1$ with $r > 1-\deg f_0$ and $j > r$, $0 < W(f_0,tf_1)/W(f_0,D_j f_1) < 1/(j-r)$ on $\mathbb{R}$.

---

## 3. Independently Certified Disproofs & Sharp Obstructions (100% Exact)

1. **Ma–Wang (Oct 2026) Theorem 2 is 100% SHARP — per-subset common interlacing & nonnegative-mixture real-rootedness FAIL for ALL $4\le i\le k-3$** (`audit_31_32_33.py`, clean integer checks):
   - **$i=4$** ($k=30$ passes exact Cauchy index; first fail **$k=31$**, deg 57):
     - Exact integer certificate: $R(x) = 932 A_{31,\{28..31\}}(x) + A_{31,\{1..4\}}(x)$ at **$x=-23$** gives $R'(-23)^2 - R(-23)R''(-23) \approx -9.065320\times 10^{279} < 0$ (exact integer!).
   - **$i=5$** ($k=15$ passes exact Cauchy index; first fail **$k=16$**, deg 26):
     - Exact integer certificate: $R(x) = 269 A_{16,\{12..16\}}(x) + A_{16,\{1..5\}}(x)$ at **$x=-11$** gives $R'(-11)^2 - R(-11)R''(-11) \approx -2.317817\times 10^{96} < 0$ (exact integer!).
   - **$i=k-3$** ($k=27$ passes exact Cauchy index; first fail **$k=28$**, deg 30):
     - Exact `Fraction` certificate at $k=28$: $R(x) = 272 A_{28,\{4..28\}}(x) + A_{28,\{1..25\}}(x)$ at $x=-73/10$ gives $-5.025603\times 10^{111} < 0$; and at $k=29$ (deg 31), exact integer certificate $R(x) = 349 A_{29,\{4..29\}}(x) + A_{29,\{1..26\}}(x)$ at **$x=-8$** gives $-4.897981\times 10^{118} < 0$ (exact integer!).
2. **Binomial-weighted bivariate stability $B_k(x,y) = \sum \binom{k}{i} A_{k,i}(x) y^i$ DISPROVED at $k=10$** (`verify_msg13.py`):
   - Exact `Fraction` Rayleigh determinant at $(0,-3/4)$: $-2649516238079970755/4294967296 < 0$.
   - Exact integer Sturm on $B_{10}(s, -2531/1000 + (130151/100)s)$: degree 19, squarefree, **17 real roots**.
3. **2×2 minor $x$-positivity & global-$\mathbb{R}$ Turán $A_{k,i}^2 > A_{k,i-1}A_{k,i+1}$ DISPROVED** (`verify_msg13.py`, `verify_deep_positivity.py`):
   - First negative $x$-coefficient: $[x^1](A_{19,8}^2 - A_{19,7}A_{19,9}) = -124876023597181272 < 0$ (off-diagonal at $(20,9,10)$: $-525307528301731032 < 0$).
   - First pointwise failure on $\mathbb{R}$: at $k=35, i=17$, exact `Fraction` $x=7/10^{17}$ gives $A_{35,17}(x)^2 - A_{35,16}(x)A_{35,18}(x) \approx -2.376737\times 10^{17} < 0$.
   - Meanwhile $M_x^2 - X_xY_x$ (1,176 pairs) and $W(A_{k,i}, A_{k,i-1})$ (1,275 pairs) have 0 negative/zero $x$-coefficients through $k=50$.
4. **Universal operator products DISPROVED (exact Sturm/Fraction)** (`audit_18_20_21.py`, `audit_31_32_33.py`):
   - Free $(a_j,b_j,w_j)$ with $0<a_j<b_j, w_j>0$: exact `Fraction` witness breaks real-rootedness at **$k=9$** (exact Sturm: `q9_3` deg 6, squarefree, **4 real roots**).
   - Fixed $a_j=1, w_j=1$, all integer $b_j\ge 2$: single-op exact integer witness ($K=1000$ copies of $b=2$, $b_{1001}=196200000000000, b_{1002}=181330000000000$, $i=960$) has `P2[960]` (deg 42) and `P2[961]` (deg 41) both real-rooted in $(-1,0)$ by exact Sturm, yet Cauchy index $=40\ne 42$ and exact $W(P_2[961],P_2[960])(-1818/10^{16}) = -3.902371\times 10^{30} < 0$. Double-op scaling limit $\mathcal{M}_{\text{dbl}}^{(w)}/\alpha$ matches SymPy identically, goes $<0$ at $w\ge 5$, and retains strictly negative $W(D_{b_2}R_3, D_{b_2}R_2)$ even **after** the outer $D_{b_2}$.
5. **GLZ-specific structural diagnostics** (`audit_25.py`, `audit_31_32_33.py`):
   - Exact $J_i(0) = \frac89\frac{(i+1)(k-i+1)}{i(k-i)}$ holds identically for constant-$b$ single-op $\prod(D_b+y)(1)$ and asymptotically for GLZ (GLZ crosses $J_i(0)<1$ at $k=35, i=17$: $0.997128$).
   - Half-plane Wronskian $W(E_{k,i}, E_{k,i+r})>0$ on $\mathbb{R}$ certified by exact Sturm for $r=1..4$ ($k\le 11$, 164/164 PASS), while $r=9$ fails at $k=10$ (exact `Fraction` $W(E_{10,0},E_{10,9})(-19/20) \approx -1.4255\times 10^{24}<0$).
   - Two-Regime Split ($D_{k+2}$-Turán $= A\Phi^2+B\Phi+C$): $A,C>0$ on both regimes ($k\le 30$; Sturm $k\le 11$); $B=t(1+t)A'$ changes sign on both regimes for $k\le 18$, while for $k\in\{35,36,38,40\}$ every point with $J_{\text{raw}}(t)\le 1$ satisfies both $\Phi>0$ and $B(t)>0$ (0/43 exceptions).
