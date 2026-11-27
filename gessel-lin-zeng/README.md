# Gessel–Lin–Zeng Conjecture 15 — Sharp Obstructions, Real-Stability Theorems & Exact Verification Suite

> **Status: OPEN — Archived with Sharp Obstruction Theorems, Exact Disproofs & General Real-Stability Theorems**  
> **Authors:** **[Atharvsinh Jadav](https://athrix.me/)** ([X `@athrix_codes`](https://x.com/athrix_codes) · [GitHub `Atharvsinh-codez`](https://github.com/Atharvsinh-codez)), with the **Gemini 4 Neon Extra High multi-agent research swarm**  
> **Source Problem:** I. M. Gessel, Z. Lin, & J. Zeng, *Jacobi–Stirling polynomials and $P$-partitions*, *European J. Combin.* **33** (2012), 1987–2000 ([arXiv:1201.0622](https://arxiv.org/abs/1201.0622)), Conjecture 15  
> **Upstream Tracking:** [`1337isnot1337/solved-conjectures` Issue #84](https://github.com/1337isnot1337/solved-conjectures/issues/84) · SymCat (`realRootedOpenProblems.htm`)

> **Mathematical Honesty Notice:** **Gessel–Lin–Zeng Conjecture 15 itself remains OPEN** for $4 \le i \le k - 3$. This archive publishes exact-arithmetic theorems proving that all three standard proof engines in the field (per-subset common interlacers, binomial bivariate real stability, and universal differential-operator products) **provably fail** on the remaining open cases, alongside standalone general real-stability and Wronskian theorems and our complete **243-script (36,875-line)** Python research and verification suite.

---

## 1. Exact Conjecture & Structural Reformulations

For $k \ge 1$ and $0 \le i \le k$, let $M_k = \{\bar{1}, 1, 1, \dots, \bar{k}, k, k\}$ with $\bar{1} < 1 < \dots < \bar{k} < k$. For each subset $S \subseteq [k]$ of size $|S| = i$, delete $\{\bar{j} : j \in S\}$ and let $A_{k,S}(x)$ enumerate Jacobi–Stirling permutations of the resulting multiset by internal descents. Set $A_{k,i}(x) := \sum_{|S|=i} A_{k,S}(x)$ ($\deg A_{k,i} = 2k - i - 1$).

- **Gessel–Lin–Zeng Conjecture 15 (2012):** *$A_{k,i}(x)$ has only real zeros for every $k \ge 1$ and $0 \le i \le k$.*  
  (Known for $i \in \{0, 1, 2, 3, k-2, k-1, k\}$ via classical results and Ma–Wang, [arXiv:2610.03111](https://arxiv.org/abs/2610.03111), Oct 2026; **open for $4 \le i \le k-3$**.)

Under the Möbius substitution $x = t/(1+t)$ ($(-\infty, 0) \leftrightarrow (-1, 0)$), $\tilde{h}_{k,i}(t) := (1+t)^{2k-i-1} A_{k,i}(t/(1+t))$ obeys the **$i$-independent** differential-operator recurrence:
$$\tilde{h}_{k,i} = D_{k+1}\bigl(D_{k+1}\tilde{h}_{k-1,i} + \tilde{h}_{k-1,i-1}\bigr), \qquad D_j := t(1+t)\partial_t + 1 + jt, \qquad \tilde{h}_{1,0} = 1 + 2t,\;\; \tilde{h}_{1,1} = 1,$$
and in $s$-space ($x = -e^s$, $\psi_{k,i}(s) := (-t)(1+t)^{k+1}\tilde{h}_{k,i}(t)$), the changing parameter $j = k+1$ disappears completely into a pure $k$-th power of a single Sturm–Liouville operator:
$$\Psi_k(s, y) := \sum_{i=0}^k \psi_{k,i}(s) y^i = \mathcal{T}_y^k\!\left(\frac{e^s}{1 + e^s}\right), \qquad \mathcal{T}_y := \frac{1}{1 + e^s}\partial_s(\partial_s + y).$$

---

## 2. Main Verified Theorems & Exact Disproofs

Every result below has been independently audited in 100% exact integer / `fractions.Fraction` / rational Sturm & Cauchy-index / SymPy arithmetic (see [`AUDIT.md`](AUDIT.md)):

### A. Jadav's Ma–Wang Sharpness Theorem (100% Exact Boundary Obstruction across ALL $4 \le i \le k - 3$)

Ma and Wang ([arXiv:2610.03111](https://arxiv.org/abs/2610.03111), Thm. 2) proved real-rootedness for $i \in \{1, 2, 3, k-2, k-1\}$ by proving every nonnegative linear combination of $\{A_{k,S} : |S| = i\}$ is real-rooted (common interlacer). Using the rootwise extremal subsets $A_{k,\text{late}} = A_{k,\{k-i+1,\dots,k\}}$ and $A_{k,\text{early}} = A_{k,\{1,\dots,i\}}$, we prove **Ma–Wang Theorem 2 is 100% sharp** at every next boundary ($i = 4$, $i = 5$, and $i = k-3$), via exact integer/Fraction Laguerre violations $R'(x)^2 - R(x)R''(x) < 0$ (verified in [`verify/audit_31_32_33.py`](verify/audit_31_32_33.py) and [`verify/verify_clean_laguerre_certs.py`](verify/verify_clean_laguerre_certs.py)):

| Boundary Family | Exact First Fail $k$ ($k-1$ passes exact Cauchy index) | Degree | Exact Mixture $R(x)$ | Evaluation Point | Exact Laguerre Value $R'(x)^2 - R(x)R''(x)$ |
| :--- | :---: | :---: | :--- | :---: | :---: |
| **$i = 4$** | **$k = 31$** ($k = 30$ passes) | 57 | $932\,A_{31,\{28..31\}} + A_{31,\{1..4\}}$ | **$x = -23$** (integer) | $\approx -9.065320 \times 10^{279} < 0$ *(280-digit integer)* |
| **$i = 5$** | **$k = 16$** ($k = 15$ passes) | 26 | $269\,A_{16,\{12..16\}} + A_{16,\{1..5\}}$ | **$x = -11$** (integer) | $\approx -2.317817 \times 10^{96} < 0$ *(97-digit integer)* |
| **$i = k - 3$** | **$k = 28$** ($k = 27$ passes) | 30 | $272\,A_{28,\{4..28\}} + A_{28,\{1..25\}}$ | **$x = -73/10$** | $\approx -5.025603 \times 10^{111} < 0$ *(exact `Fraction`)* |
| **$i = k - 3$** *(int pt)* | **$k = 29$** | 31 | $349\,A_{29,\{4..29\}} + A_{29,\{1..26\}}$ | **$x = -8$** (integer) | $\approx -4.897981 \times 10^{118} < 0$ *(119-digit integer)* |

### B. Exact Disproofs of Bivariate Stability, Turán Positivity & Universal Operator Products

1. **Binomial-Weighted Bivariate Real Stability $B_k(x, y) = \sum_{i=0}^k \binom{k}{i} A_{k,i}(x) y^i$ DISPROVED at $k = 10$** ([`verify/verify_msg13.py`](verify/verify_msg13.py)):
   - Exact `Fraction` Rayleigh determinant at $(0, -3/4)$:
     $$B_{10,x}B_{10,y} - B_{10}B_{10,xy}\Big|_{(0, -3/4)} = -\frac{2649516238079970755}{4294967296} < 0.$$
   - Exact integer Sturm witness: $B_{10}\bigl(s, -\frac{2531}{1000} + \frac{130151}{100}s\bigr)$ is squarefree of degree $19$ with only **$17$ real roots**.
2. **Coefficientwise $x$-Positivity & Global-$\mathbb{R}$ Pointwise Turán DISPROVED** ([`verify/verify_deep_positivity.py`](verify/verify_deep_positivity.py)):
   - First negative $x$-coefficient: $[x^1](A_{19,8}^2 - A_{19,7}A_{19,9}) = -124{,}876{,}023{,}597{,}181{,}272 < 0$.
   - First pointwise failure on $\mathbb{R}$: at $k = 35, i = 17$, $x = 7/10^{17}$ gives $A_{35,17}(x)^2 - A_{35,16}(x)A_{35,18}(x) \approx -2.376737 \times 10^{17} < 0$.
3. **Universal Operator Products $\prod (D_{b_j} + y)(1)$ and $\prod D_{b_j}(D_{b_j} + y)(1)$ DISPROVED** ([`verify/audit_18_20_21.py`](verify/audit_18_20_21.py), [`verify/audit_31_32_33.py`](verify/audit_31_32_33.py)):
   - **Free parameters ($0 < a_j < b_j, w_j > 0$):** Exact Sturm counterexample breaks real-rootedness at **$k = 9$** (`q9_3`, deg 6, squarefree, 4 real roots).
   - **All-integer $b_j \ge 2$ with $a_j = w_j = 1$ (single-op):** Exact integer/Fraction certificate ($K = 1000$ copies of $b = 2$, $b_{1001} = 196{,}200{,}000{,}000{,}000$, $b_{1002} = 181{,}330{,}000{,}000{,}000$, $i = 960$) breaks interlacing with exact Cauchy index $40 \ne 42$ and exact Wronskian $-3.902371 \times 10^{30} < 0$.
   - **All-integer $b_j \ge 2$ (double-op):** Exact SymPy scaling limit $\mathcal{M}_{\text{dbl}}^{(w)}(\Phi)$ turns strictly negative for all $w \ge 5$, and the outer $D_{b_2}$ fails to heal the broken interlacing.
   - **No linear gauge eliminates Stage 3** ($D_{k+1}^i \mapsto D_{k+2}^i$): 5-line algebraic impossibility proof (`ALGEBRA_ROUTE.md` §3).

### C. Standalone Proved General Real-Stability & Wronskian Theorems

1. **Jadav's Cubic Real-Stability Addition Theorem** ([`ALGEBRA_ROUTE.md`](ALGEBRA_ROUTE.md), [`verify/audit_25.py`](verify/audit_25.py)): For 5 vectors $v_0, \dots, v_4 \in \mathbb{R}^2$ in strict open-half-plane order ($W(v_r, v_s) > 0$ for $0 \le r < s \le 4$), PSD $3 \times 3$ Bezoutiants of $P_1(y) = v_0 - 3v_1 y + 3v_2 y^2 - v_3 y^3$ and $P_2(y) = v_1 - 3v_2 y + 3v_3 y^2 - v_4 y^3$ unconditionally imply PSD Bezoutiant of $P_1(y) + P_2(y)$ (with sharpness witness $L = (-3000, 1, 10, 1, -3000)$ breaking the result without the Wronskian cone).
2. **Sharpest Quantitative Plücker Quad Theorem B** ([`verify/audit_18_20_21.py`](verify/audit_18_20_21.py), [`verify/audit_25.py`](verify/audit_25.py)): Under $W_{AB}, W_{BC}, W_{CE}, W_{AC}, W_{BE} > 0$, $W_{AE} \ge 0$, and $J(A, B, C), J(B, C, E) \ge 1 + \varepsilon$ ($\varepsilon \ge 0$), unconditionally $W(A+B, C+E) > 0$ and $J(A+B, B+C, C+E) \ge \frac{4(1+\varepsilon)^2}{(2+\varepsilon)^2} \ge 1$ with **zero upper bound** on the Plücker ratio $R = W_{AB}W_{CE}/W_{BC}^2$.
3. **Universal Cubic Wronskian Log-Concavity** ([`verify/audit_18_20_21.py`](verify/audit_18_20_21.py)): Global sum-of-squares proof that $W(p_0, p_1) W(p_2, p_3) \le W(p_1, p_2)^2$ on $\mathbb{R}$ for every real-stable cubic $P(t, y) = \sum_{j=0}^3 p_j(t) y^j$.
4. **Zero-Loss $D_m$ Preservation ($m \ge 1$)**, **Theorem A (Turán $\implies$ Interlacing)**, **Theorem C (Pole-Gap Wronskian Bound)**, and the **Two-Regime Turán Dichotomy** $(D_m f)^2 - (D_m g)(D_m h) = A\Phi^2 + B\Phi + C$ with $B^2 - 4AC = t^2(1+t)^2 W_{gh}^2(1 - J_{\text{raw}})$ ([`ANALYTIC_ROUTE.md`](ANALYTIC_ROUTE.md), [`AUDIT.md`](AUDIT.md)).

---

## 3. Manuscripts & Complete Python Suite Inventory

### Research & Audit Manuscripts (4 `.md` files)
- **[`AUDIT.md`](AUDIT.md)** — Complete independent `@verifier` audit of every positive theorem, hypothesis check, and exact disproof.
- **[`ALGEBRA_ROUTE.md`](ALGEBRA_ROUTE.md)** — Exact operator identities, 4 proved positive theorems, Ma–Wang sharpness table, and ruled-out algebraic routes.
- **[`ANALYTIC_ROUTE.md`](ANALYTIC_ROUTE.md)** — Pure-power Sturm–Liouville $\mathcal{T}_y^k$, Borcea–Brändén symbol calculus, Two-Regime Turán Dichotomy, universal operator disproofs, and high-precision ($k \le 72$) diagnostics.
- **[`SOURCES_AND_STATE.md`](SOURCES_AND_STATE.md)** — Literature & prior-work audit (`arXiv:1201.0622`, Ma–Wang `arXiv:2610.03111`, SymCat, Issue #84).

### Python Verification & Research Suite (**243 `.py` files, 36,875 lines**)

| Directory | `.py` Files | Total Lines | Contents & Key Entry Points |
|---|---:|---:|---|
| [`verify/`](verify/) | **8** | **2,224** | **Standalone portable exact-arithmetic verifiers:** [`testbed.py`](verify/testbed.py), [`probes.py`](verify/probes.py), [`verify_msg13.py`](verify/verify_msg13.py), [`verify_deep_positivity.py`](verify/verify_deep_positivity.py), [`audit_18_20_21.py`](verify/audit_18_20_21.py), [`audit_25.py`](verify/audit_25.py), [`audit_31_32_33.py`](verify/audit_31_32_33.py), [`verify_clean_laguerre_certs.py`](verify/verify_clean_laguerre_certs.py). |
| [`algebra/`](algebra/) | **104** | **17,594** | **Algebraic proof & obstruction scripts:** [`verify_all_algebra_theorems.py`](algebra/verify_all_algebra_theorems.py), [`exact_laguerre_barriers.py`](algebra/exact_laguerre_barriers.py), [`clean_laguerre_certs.py`](algebra/clean_laguerre_certs.py), [`final_go_pivot_check.py`](algebra/final_go_pivot_check.py), [`surgical_checkpoint.py`](algebra/surgical_checkpoint.py), plus 99 Wronskian, Plücker-quad, cubic Bezoutiant, and Sturm probe scripts. |
| [`analyst/`](analyst/) | **131** | **17,057** | **Analytic & adversarial probe scripts:** [`fast_certify_bge2_w1.py`](analyst/fast_certify_bge2_w1.py), [`check_double_M_asymptotic.py`](analyst/check_double_M_asymptotic.py), [`targeted_attack_doubleT.py`](analyst/targeted_attack_doubleT.py), [`verify_k72_and_turan_terms.py`](analyst/verify_k72_and_turan_terms.py), [`check_two_regime_highk.py`](analyst/check_two_regime_highk.py), [`verify_plucker_and_cubic_exact.py`](analyst/verify_plucker_and_cubic_exact.py), plus 125 operator-stability, Sturm–Liouville, and scaling-limit scripts. |
| **Total** | **243** | **36,875** | **100% reproducible exact integers, `Fraction`, Sturm/Cauchy chains, SymPy & `mpmath`.** |

---

## 4. Exact Reproduction Instructions

Requirements: **Python 3.10+** with `sympy` and `numpy` (and `mpmath` for optional high-precision `analyst/` scripts).

```bash
pip install sympy numpy mpmath
```

Run from repository root:
```bash
# 1. Direct Jacobi-Stirling enumeration + exact Sturm/Cauchy baseline (k <= 24)
python gessel-lin-zeng/verify/testbed.py

# 2. Supplementary aggregate & 4th-cross-pair probes
python gessel-lin-zeng/verify/probes.py

# 3. Exact disproof of B_10(x,y) real stability (Rayleigh det & Sturm line witness)
python gessel-lin-zeng/verify/verify_msg13.py

# 4. Exact first negative x-coeff (k=19, i=8) and pointwise R-failure (k=35, i=17)
python gessel-lin-zeng/verify/verify_deep_positivity.py

# 5. Audit of Theorem A, Plucker Quad B, Zero-Loss D_m, Cubic Wronskian SOS, k=9 disproof
python gessel-lin-zeng/verify/audit_18_20_21.py

# 6. Audit of Jadav's Cubic Real-Stability Addition, J_i(0), and r=9 Wronskian disproof
python gessel-lin-zeng/verify/audit_25.py

# 7. Exact Ma-Wang sharpness certificates (i=4,5,k-3) & all-integer b_j >= 2 disproofs
python gessel-lin-zeng/verify/audit_31_32_33.py

# 8. Standalone zero-dependency check of integer/Fraction Laguerre certificates (x=-23, -11, -73/10, -8)
python gessel-lin-zeng/verify/verify_clean_laguerre_certs.py
```

---

*See also the repository-wide catalogues: [`THEOREMS.md`](../THEOREMS.md) · [`EQUATIONS.md`](../EQUATIONS.md) · [`JOURNEY.md`](../JOURNEY.md) · [Root `README.md`](../README.md).*
