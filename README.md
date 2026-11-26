# Open Conjectures Research & Exact Verification Archive

> **Curated by [Atharsinh Jadav](https://athrix.me/)** — [Website](https://athrix.me/) · [X (`@athrix_codes`)](https://x.com/athrix_codes) · [GitHub (`Atharvsinh-codez`)](https://github.com/Atharvsinh-codez) *(with Gemini 4 Neon Extra High multi-agent research swarm)* · Licensed under the [MIT License](LICENSE).

Reproducible research archives, exact-arithmetic certificates, general algebraic/analytic theorems, and live openness audits for open problems in algebraic combinatorics and multivalued group theory.

Every mathematical theorem, counterexample, and barrier certificate in this repository is backed by standalone, exact-arithmetic Python scripts (exact integers, `fractions.Fraction`, exact rational Sturm/Cauchy-index chains, and SymPy — **zero** floating-point handwaving).

---

## Contents

| Directory / File | Target | Status | Highlights |
| :--- | :--- | :--- | :--- |
| [`gessel-lin-zeng/`](gessel-lin-zeng/) | **Gessel–Lin–Zeng Conjecture 15** ([arXiv:1201.0622](https://arxiv.org/abs/1201.0622), [Issue #84](https://github.com/1337isnot1337/solved-conjectures/issues/84)) | **Open** *(archived with sharp theorems & obstructions)* | Exact proof that Ma–Wang (Oct 2026) Thm. 2 is **100% sharp** across all $4 \le i \le k-3$; exact disproofs of $B_k(x,y)$ stability ($k=10$) and universal operator products; general real-stability theorems. |
| [`bardakov-kozlovskaya/`](bardakov-kozlovskaya/) | **Bardakov–Kozlovskaya Conjecture 3.3(2)** ([arXiv:2305.07261](https://arxiv.org/abs/2305.07261), [Issue #174](https://github.com/1337isnot1337/solved-conjectures/issues/174)) | **Active research** *(in progress)* | Exact integer verification of **11,069 coefficients** through $n=100$ (0 zeros, 0 violations of sign law $(-1)^{b(c+1)+c}$); exact $c=1, b=0$ closed formula; 5 exact representations catalogued. |
| [`UNLISTED_CONJECTURES.md`](UNLISTED_CONJECTURES.md) | Audit of [`1337isnot1337/solved-conjectures`](https://github.com/1337isnot1337/solved-conjectures) | Reference catalogue | Comprehensive deduplicated inventory of **308** still-unsolved conjectures across `main` partial results, `docs/*-research.md`, `origin/research`, open GitHub issues, and classic open problems. |

---

## 1. Gessel–Lin–Zeng Conjecture 15 (`gessel-lin-zeng/`)

**Conjecture 15** (Gessel–Lin–Zeng 2012, *European J. Combin.* **33**, 1987–2000) asks whether the Jacobi–Stirling descent polynomials $A_{k,i}(x) = \sum_{|S|=i} A_{k,S}(x)$ ($\deg = 2k-i-1$) have only real zeros for every $k \ge 1$ and $0 \le i \le k$. Ma and Wang ([arXiv:2610.03111](https://arxiv.org/abs/2610.03111), Oct 2026) proved real-rootedness for $i \in \{1,2,3,k-2,k-1\}$ by showing every nonzero nonnegative linear combination of $\{A_{k,S} : |S|=i\}$ is real-rooted (common interlacer); the range $4 \le i \le k-3$ remains **open**.

> **Honest status:** Gessel–Lin–Zeng Conjecture 15 itself remains **open**. This archive establishes sharp structural theorems, exact boundary obstructions ruling out every natural local/common-interlacer/universal-operator extension, and independent exact-arithmetic verification of every claim.

### Key Verified Results (see [`gessel-lin-zeng/AUDIT.md`](gessel-lin-zeng/AUDIT.md))

1. **Ma–Wang (Oct 2026) Theorem 2 is 100% sharp across ALL open cases $4 \le i \le k-3$:**
   Using rootwise extremal deletion schedules $A_{k,\text{late}} = A_{k,\{k-i+1,\dots,k\}}$ and $A_{k,\text{early}} = A_{k,\{1,\dots,i\}}$, non-real nonnegative mixtures $R(x) = \lambda A_{k,\text{late}}(x) + A_{k,\text{early}}(x)$ occur at every remaining boundary, certified by exact integer/rational Laguerre violations $R'(x)^2 - R(x)R''(x) < 0$:
   - **$i = 4$** ($k=30$ passes exact Cauchy index; first failure **$k=31$**, $\deg 57$): $R(x) = 932\,A_{31,\{28..31\}}(x) + A_{31,\{1..4\}}(x)$ at **$x = -23$** gives $R'(-23)^2 - R(-23)R''(-23) \approx -9.065320 \times 10^{279} < 0$ (exact integer).
   - **$i = 5$** ($k=15$ passes exact Cauchy index; first failure **$k=16$**, $\deg 26$): $R(x) = 269\,A_{16,\{12..16\}}(x) + A_{16,\{1..5\}}(x)$ at **$x = -11$** gives $R'(-11)^2 - R(-11)R''(-11) \approx -2.317817 \times 10^{96} < 0$ (exact integer).
   - **$i = k-3$** ($k=27$ passes exact Cauchy index; first failure **$k=28$**, $\deg 30$): $R(x) = 272\,A_{28,\{4..28\}}(x) + A_{28,\{1..25\}}(x)$ at **$x = -73/10$** gives $-5.025603 \times 10^{111} < 0$; and at **$k=29$** ($\deg 31$), $R(x) = 349\,A_{29,\{4..29\}}(x) + A_{29,\{1..26\}}(x)$ at **$x = -8$** gives $-4.897981 \times 10^{118} < 0$ (exact integer).
2. **Exact disproof of the surviving binomial-weighted bivariate stability candidate $B_k(x,y) = \sum_{i=0}^k \binom{k}{i} A_{k,i}(x) y^i$ at $k = 10$:**
   - Exact `Fraction` Rayleigh determinant at $(x,y) = (0, -3/4)$:
     $$B_{10,x}B_{10,y} - B_{10}B_{10,xy}\Big|_{(0,-3/4)} = -\frac{2649516238079970755}{4294967296} < 0.$$
   - Exact integer Sturm witness: $B_{10}\!\left(s,\, -\frac{2531}{1000} + \frac{130151}{100}s\right)$ is squarefree of degree $19$ with only **$17$ real roots**.
3. **Exact disproofs of universal operator products $\prod (D_{b_j} + y)(1)$ and $\prod D_{b_j}(D_{b_j} + y)(1)$:**
   - Free parameters $0 < a_j < b_j, w_j > 0$: exact Sturm counterexample breaks real-rootedness at **$k=9$** (`q9_3`, deg 6, squarefree, 4 real roots).
   - Fixed $a_j = 1, w_j = 1$ with all-integer $b_j \ge 2$: exact finite integer/Fraction witness ($K=1000$ copies of $b=2$, $b_{1001}=196200000000000$, $b_{1002}=181330000000000$, $i=960$) breaks interlacing (exact Cauchy index $40 \ne 42$; exact Wronskian $-3.902371 \times 10^{30} < 0$), and the exact SymPy double-operator scaling limit turns negative for all $w \ge 5$ even after the outer $D_{b_2}$.
4. **Proved general real-stability & Wronskian theorems (independently audited):**
   - **Universal Cubic Real-Stability Addition:** For $v_0,\dots,v_4 \in \mathbb{R}^2$ in strict open-half-plane order ($W(v_r,v_s)>0$ for $r<s$), PSD $3\times 3$ Bezoutiants of $P_1(y) = v_0 - 3v_1 y + 3v_2 y^2 - v_3 y^3$ and $P_2(y) = v_1 - 3v_2 y + 3v_3 y^2 - v_4 y^3$ imply PSD Bezoutiant of $P_1+P_2$.
   - **Sharpest Quantitative Plücker Quad Theorem:** Under $W_{AB}, W_{BC}, W_{CE}, W_{AC}, W_{BE} > 0$, $W_{AE} \ge 0$, and $J(A,B,C), J(B,C,E) \ge 1+\varepsilon$ ($\varepsilon \ge 0$), unconditionally $W(A+B,C+E)>0$ and $J(A+B, B+C, C+E) \ge \frac{4(1+\varepsilon)^2}{(2+\varepsilon)^2}$ with zero upper bound on the Plücker ratio.
   - **Universal Cubic Wronskian Log-Concavity:** Exact global sum-of-squares proof that $W(p_0,p_1)W(p_2,p_3) \le W(p_1,p_2)^2$ on $\mathbb{R}$ for every real-stable cubic $P(t,y) = \sum_{j=0}^3 p_j(t)y^j$.
   - **Zero-Loss $D_m$ Preservation ($m \ge 1$)** and **Two-Regime Turán Dichotomy** $(D_m f)^2 - (D_m g)(D_m h) = A\Phi^2 + B\Phi + C$.

**Documents:** [`SOURCES_AND_STATE.md`](gessel-lin-zeng/SOURCES_AND_STATE.md) · [`ALGEBRA_ROUTE.md`](gessel-lin-zeng/ALGEBRA_ROUTE.md) · [`ANALYTIC_ROUTE.md`](gessel-lin-zeng/ANALYTIC_ROUTE.md) · [`AUDIT.md`](gessel-lin-zeng/AUDIT.md)

---

## 2. Bardakov–Kozlovskaya Conjecture 3.3(2) (`bardakov-kozlovskaya/`) — Active Research

For the symmetric $n$-valued Buchstaber–Novikov group polynomial
$$p_n(x,y,z) = \prod_{k=0}^{n-1}\Bigl(z - (-1)^n(\alpha+\zeta^k\beta)^n\Bigr) = \sum_{a+2b+3c=n} c_{a,b,c}^{(n)}\, e_1^a e_2^b e_3^c \qquad (\alpha^n=x,\;\beta^n=y,\;\zeta=e^{2\pi i/n}),$$
**Conjecture 3.3(2)** (Bardakov–Kozlovskaya 2023, [arXiv:2305.07261](https://arxiv.org/abs/2305.07261)) asserts that for every even $n = 2m \ge 2$ and every admissible triple $(a,b,c) \in \mathbb{Z}_{\ge 0}^3$ with $a+2b+3c=2m$, $c_{a,b,c}^{(2m)} \ne 0$. (Part (1), prime-power divisibility, already follows from Buchstaber–Kornev 2025/2026 and a direct resultant Frobenius proof.)

### Current Verified State (see [`bardakov-kozlovskaya/SOURCES_AND_STATE.md`](bardakov-kozlovskaya/SOURCES_AND_STATE.md))
- **Exact three-way convention verification (`verify/engine.py`):** SymPy resultant $\operatorname{Res}_u(u^n-x, (t+u)^n-(-1)^n y)|_{t^n=z}$, Buchstaber–Kornev's $n\times n$ $y$-circulant Wendt determinant, and fast exact Newton-identity + Horner engine agree identically.
- **Exact integer verification of 11,069 coefficients (`verify/stress_to_84.py`):**
  - All even $n = 2, 4, \dots, 84$ (9,464 coefficients) plus $n = 90$ (721) and $n = 100$ (884): **0 zero coefficients** and **0 violations** of the exact sign law $\operatorname{sign}(c_{a,b,c}^{(2m)}) = (-1)^{b(c+1)+c}$ (and 0 violations of $(-1)^{c(b+1)}$ for all odd $n = 1, 3, \dots, 31$, $c \ge 1$).
- **Exact $c=1, b=0$ closed formula** (verified for all even $n=4..40$ and odd $n=3..21$):
  $$c_{n-3,0,1}^{(n)} = -n\left(\frac{1}{2}\binom{2n}{n} - 2 - (-1)^n\right).$$
- **Full catalogue of 5 exact representations** (double product, integer resultant, Wendt $y$-circulant determinant, Gaiur–Rubtsov–van Straten / Buchstaber–Kornev univariate discriminant & rational curve, Newton power sums).

---

## 3. Unlisted & Unsolved Conjectures Catalogue (`UNLISTED_CONJECTURES.md`)

A deduplicated audit of [`1337isnot1337/solved-conjectures`](https://github.com/1337isnot1337/solved-conjectures) indexing **308** distinct unsolved conjectures (15 partial results on `main`, 18 `docs/*-research.md` targets, 170 `origin/research` targets, 54 open GitHub issues without research folders, and 51 classic open problems), separated cleanly from 151 solved entries, 13 unpublished candidate full proofs, and 22 confirmed prior external resolutions.

---

## 4. Exact Reproduction Instructions

Requirements: **Python 3.10+** with `sympy` (and `numpy` for `gessel-lin-zeng/verify/audit_25.py`). Install once via:

```bash
pip install sympy numpy
```

### Gessel–Lin–Zeng (`gessel-lin-zeng/verify/`)
```bash
# 1. Direct Jacobi-Stirling permutation enumeration + exact Sturm/Cauchy baseline (k <= 24)
python gessel-lin-zeng/verify/testbed.py

# 2. Supplementary aggregate & 4th-cross-pair probes
python gessel-lin-zeng/verify/probes.py

# 3. Exact disproof of B_10(x,y) real stability (Rayleigh det & Sturm line witness) + minor checks
python gessel-lin-zeng/verify/verify_msg13.py

# 4. Exact first negative x-coeff (k=19, i=8) and pointwise R-failure (k=35, i=17) of Turan
python gessel-lin-zeng/verify/verify_deep_positivity.py

# 5. Audit of Theorem A, Quad Theorem B, Zero-Loss D_m, Cubic Wronskian SOS, and k=9 free-param disproof
python gessel-lin-zeng/verify/audit_18_20_21.py

# 6. Audit of Universal Cubic Real-Stability Addition, J_i(0), and r=9 Wronskian disproof
python gessel-lin-zeng/verify/audit_25.py

# 7. Exact integer/Fraction Ma-Wang sharpness certificates (i=4,5,k-3) & all-integer b_j >= 2 disproofs
python gessel-lin-zeng/verify/audit_31_32_33.py

# 8. Standalone zero-dependency check of the clean integer/Fraction Laguerre points (x=-23, -11, -73/10, -8)
python gessel-lin-zeng/verify/verify_clean_laguerre_certs.py
```

### Bardakov–Kozlovskaya (`bardakov-kozlovskaya/verify/`)
```bash
# 1. Three-way definition verification (resultant vs Wendt det vs Newton engine) & exact c=1, b=0 formula
python bardakov-kozlovskaya/verify/engine.py

# 2. Exact integer verification of all even n = 2..84 + n = 90, 100 (11,069 coefficients)
python bardakov-kozlovskaya/verify/stress_to_84.py
```

---

## 5. Attribution, Citation & License

- **Author / Curator:** **Atharsinh Jadav** — Website: https://athrix.me/ · X: https://x.com/athrix_codes · GitHub: https://github.com/Atharvsinh-codez *(with Gemini 4 Neon Extra High multi-agent research swarm)*.
- **Citation:** See [`CITATION.cff`](CITATION.cff) (or use GitHub's *Cite this repository* button).
- **License:** Released under the [MIT License](LICENSE) (`Copyright (c) 2026 Atharsinh Jadav`).

