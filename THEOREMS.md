# Catalogue of Proved Theorems, Lemmas & Exact Disproofs

> **Authors:** **[Atharsinh Jadav](https://athrix.me/)** ([X `@athrix_codes`](https://x.com/athrix_codes) · [GitHub `Atharvsinh-codez`](https://github.com/Atharvsinh-codez)), with the **Gemini 4 Neon Extra High multi-agent research swarm**  
> **Scope:** Complete, audited inventory of every theorem, corollary, proposition, lemma, and exact-arithmetic disproof established across both campaigns in this repository. Every entry links directly to its full proof manuscript and independent verification script(s).

---

## Executive Index

| Campaign | Upstream Target | Status | Flagship Results |
|---|---|---|---|
| **Part I: Multivalued Group Polynomials** | **Bardakov–Kozlovskaya Conjecture 3.3(2)** ([arXiv:2305.07261](https://arxiv.org/abs/2305.07261), [Issue #174](https://github.com/1337isnot1337/solved-conjectures/issues/174), [MathDB #358383](https://mathdb.com/p/358383/conjecture-on-coefficients-of-multivalued-group-polynomials)) | **✅ SOLVED & 100% VERIFIED** | **The Jadav Sign-Law & Real-Rootedness Theorem** (Thm 1.1) & **Corollary 1.2** (full resolution for all even $n = 2m \ge 2$), backed by 10/10 signed-off components. |
| **Part II: Jacobi–Stirling Polynomials & Real Stability** | **Gessel–Lin–Zeng Conjecture 15** ([arXiv:1201.0622](https://arxiv.org/abs/1201.0622), [Issue #84](https://github.com/1337isnot1337/solved-conjectures/issues/84)) | **OPEN** *(sharp obstructions & general theorems proved)* | **Jadav's Ma–Wang Sharpness Theorem**, **Jadav's Cubic Real-Stability Addition Theorem**, **Sharpest Quantitative Plücker Quad Theorem**, **Universal Cubic Wronskian Log-Concavity**, and exact disproofs of $B_{10}(x,y)$ stability & universal operator products. |

---

## Part I. Bardakov–Kozlovskaya Conjecture 3.3(2) — Complete Resolution

Full manuscript: [`bardakov-kozlovskaya/PROOF.md`](bardakov-kozlovskaya/PROOF.md) · Independent audit: [`bardakov-kozlovskaya/AUDIT.md`](bardakov-kozlovskaya/AUDIT.md)

### I.1 Flagship Main Theorem & Resolution Corollary

1. **Theorem BK-1.1 — The Jadav Sign-Law & Real-Rootedness Theorem** ([`PROOF.md` §1.2](bardakov-kozlovskaya/PROOF.md))  
   Let $n = 2m \ge 2$ be any even integer, and expand the $n$-valued Buchstaber–Novikov group polynomial as $p_{2m}(x, y, z) = \sum_{a+2b+3c=2m} c_{a,b,c}^{(2m)}\, e_1^a e_2^b e_3^c$.
   - **(Exact Universal Sign Law — ALL even $n = 2m \ge 2$):** Every admissible coefficient $c_{a,b,c}^{(2m)}$ is nonzero and satisfies
     $$\operatorname{sign}\!\bigl(c_{a,b,c}^{(2m)}\bigr) = (-1)^{b(c+1)+c} = (-1)^{b(a+1)+c} = \begin{cases} (-1)^b & \text{if } c \text{ is even,} \\ -1 & \text{if } c \text{ is odd.} \end{cases}$$
     Equivalently, every admissible coefficient of $P_{2m}^{\mathrm{even}}(e_1, -e_2, e_3)$ and of $-P_{2m}^{\mathrm{odd}}(e_1, e_2, e_3)$ is strictly positive.
   - **(Universal Slice Real-Rootedness — exact scope):**
     - *Interior slices ($n \ge 132$, $2 \le c \le \lfloor (2m-3)/3 \rfloor$):* $Q_c^{(2m)}(1, u) := \sum_{b=0}^{D_c} c_{2m-3c-2b,b,c}^{(2m)} u^b$ ($D_c = \lfloor(2m-3c)/2\rfloor$) has $D_c$ **simple real roots** — all strictly positive in $(0, +\infty)$ (with at least $D_c - 1$ in $(0, 1/4)$) when $c$ is even, and all strictly negative in $(-\infty, 0)$ when $c$ is odd — strictly separated by the universal trigonometric sample points $u_j$.
     - *All four boundary slices (all even $n \ge 2$):* $Q_0^{(2m)}(1, u) = (1 - 4u)^m$; $e_2^{-m} P_{2m}(0, e_2, e_3)$ has $\lfloor m/3 \rfloor$ simple positive real roots in $e_3^2/e_2^3$; $e_1^{-2m} P_{2m}(e_1, 0, e_3)$ has $\lfloor 2m/3 \rfloor$ simple positive real roots in $e_3/e_1^3$; and for $c = 1$, $c_{2m-3-2b, b, 1}^{(2m)} = 2m((-4)^{b+1}\binom{m-1}{b+1} - \gamma_{b+1}) < 0$ for all $0 \le b \le m - 2$. *(Note: $c = 1$ proves strict negativity without claiming real-rootedness of $Q_1^{(2m)}(1, u)$.)*

2. **Corollary BK-1.2 — Complete Resolution of Bardakov–Kozlovskaya Conjecture 3.3** ([`PROOF.md` §1.2](bardakov-kozlovskaya/PROOF.md))  
   For every even $n = 2m \ge 2$ and every $(a, b, c) \in \mathbb{Z}_{\ge 0}^3$ with $a + 2b + 3c = 2m$, $c_{a,b,c}^{(2m)} \ne 0$. Combined with Part (1) (prime-power divisibility; Buchstaber–Kornev 2025/2026), this completely resolves Conjecture 3.3 of Bardakov and Kozlovskaya (2023).

---

### I.2 Finite Base, Exact Structural Resultants & Boundary Theorems

3. **Theorem BK-2.1 — Exact Integer Base Certification ($n = 2, 4, \dots, 130$)** ([`PROOF.md` §2](bardakov-kozlovskaya/PROOF.md); [`verify/engine.py`](bardakov-kozlovskaya/verify/engine.py), [`verify/stress_to_84.py`](bardakov-kozlovskaya/verify/stress_to_84.py), [`verify/fast_certify_140.py`](bardakov-kozlovskaya/verify/fast_certify_140.py))  
   Across all $65$ even degrees $n \in \{2, 4, \dots, 130\}$, all **$33{,}417$ exact integer coefficients** ($9{,}464$ for $n \le 84$; $23{,}953$ for $86 \le n \le 130$) are nonzero and satisfy $\operatorname{sign}(c_{a,b,c}^{(n)}) = (-1)^{b(c+1)+c}$ with $0$ zeros and $0$ violations. Cross-checked identically against SymPy resultants ($n \le 8$) and the Buchstaber–Kornev Wendt determinant ($n \le 6$).

4. **Theorem BK-3.1 — Depressed-Cubic Resultant** ([`PROOF.md` §3.1](bardakov-kozlovskaya/PROOF.md); [`algebra/probe_depressed_cubic.py`](bardakov-kozlovskaya/algebra/probe_depressed_cubic.py), [`verify/audit_bk_40.py`](bardakov-kozlovskaya/verify/audit_bk_40.py))  
   For $3 \nmid m$, $P_{2m}(e_1, e_2, \eta^m) = \operatorname{Res}_u(\Pi_{2m}(u, \eta) - e_1,\, \Sigma_{2m}(u, \eta) - e_2)$ with exact normalization $+1$; for $m = 3K$, $P_{6K}(e_1, e_2, e_3)$ factors over cube roots of unity as $\prod_{\omega^3=1} R_K(e_1, e_2, \omega e_3^{1/3})$.

5. **Theorem BK-3.2 — Half-Size Lucas–Chebyshev Resultant** ([`PROOF.md` §3.1](bardakov-kozlovskaya/PROOF.md); [`verify/audit_bk_40.py`](bardakov-kozlovskaya/verify/audit_bk_40.py))  
   With $L_{2m}(t(1-t)) = t^{2m} + (1-t)^{2m}$, for every even $n = 2m \ge 2$ and all $p, q \in \mathbb{C}$,
   $$P_{2m}(1 + q - p,\, q,\, p) = \operatorname{Res}_u\bigl(L_{2m}(u) + u^{2m} - q,\; u^{2m} - p\bigr) = \prod_{k=0}^{2m-1}\bigl(L_{2m}(\zeta^k p^{1/(2m)}) + p - q\bigr).$$

6. **Theorem BK-3.3 — Boundary Family $c = 0$ (All Even $n \ge 2$)** ([`PROOF.md` §3.2](bardakov-kozlovskaya/PROOF.md))  
   $Q_0^{(2m)}(e_1, e_2) = (e_1^2 - 4e_2)^m$, hence $c_{2m-2b, b, 0}^{(2m)} = (-4)^b\binom{m}{b}$ with exact sign $(-1)^b$.

7. **Theorem BK-3.4 — Boundary Families $a = 0$ and $b = 0$ (All Even $n \ge 2$)** ([`PROOF.md` §3.2](bardakov-kozlovskaya/PROOF.md); [`verify/audit_bk_43.py`](bardakov-kozlovskaya/verify/audit_bk_43.py))  
   Exact trigonometric differentiation of $L_{2m}(1/(4\cos^2\theta)) = 2\cos(2m\theta)/(2\cos\theta)^{2m}$ isolates all $\lfloor m/3 \rfloor$ simple positive real roots of $f_a^{(2m)}(w)$ (at $\theta_r = \frac{r\pi}{2m-1}$) and all $\lfloor 2m/3 \rfloor$ simple positive real roots of $f_b^{(2m)}(w)$ (at $\tilde{\theta}_r = \frac{r\pi}{2m+1}$), proving $\operatorname{sign}(c_{0, m-3q, 2q}^{(2m)}) = (-1)^{m-q}$ and $\operatorname{sign}(c_{2m-3c, 0, c}^{(2m)}) = (-1)^c$.

8. **Theorem BK-3.5 — Boundary Family $c = 1$ (All Even $n \ge 4$)** ([`PROOF.md` §3.2](bardakov-kozlovskaya/PROOF.md); [`algebra/verify_c1_and_probe_c2.py`](bardakov-kozlovskaya/algebra/verify_c1_and_probe_c2.py), [`verify/audit_bk_43.py`](bardakov-kozlovskaya/verify/audit_bk_43.py), [`verify/audit_bk_50.py`](bardakov-kozlovskaya/verify/audit_bk_50.py))  
   Rolle real-rootedness of $F_r(u)$, palindromic factorization of $A_{2m-1}^{(2m)}(u)/(1+u)$, Maclaurin's inequality, and the exact alternating trigonometric sum $\gamma_{m-1} = \frac{1}{2m}\sum_{k=0}^{m-1}(-1)^k \frac{\cos\theta_k}{\sin^{2m+1}\theta_k} > 4^{m-1}$ prove $c_{2m-3-2b, b, 1}^{(2m)} = 2m((-4)^{b+1}\binom{m-1}{b+1} - \gamma_{b+1}) < 0$ for all $0 \le b \le m - 2$.

---

### I.3 Interior Analytic Machinery: 1D Circle Integral, Regime 1 & Zero-Crevice Regime 2

9. **Theorem BK-4.1 — Universal Sample Points & Root-Interleaving Reduction** ([`PROOF.md` §4.1](bardakov-kozlovskaya/PROOF.md))  
   For even $n \ge 6$ and $2 \le c \le \lfloor(n-3)/3\rfloor$, verifying $(-1)^j Q_c^{(2m)}(1, u_j) > 0$ for all strictly interior $j \in \mathcal{J}_{n,c} = \{c+1, \dots, m - \lceil c/2\rceil - \mathbf{1}_{c\text{ even}}\}$ forces $Q_c^{(2m)}(1, u)$ to have $D_c$ simple real roots of the required sign and forces $\operatorname{sign}(c_{n-3c-2b,b,c}^{(2m)}) = (-1)^{b(c+1)+c}$.

10. **Theorem BK-4.2 — Exact 4-Factor 1D Cauchy Circle Integral** ([`PROOF.md` §4.2](bardakov-kozlovskaya/PROOF.md); [`algebra/verify_exact_identities_mp.py`](bardakov-kozlovskaya/algebra/verify_exact_identities_mp.py))  
    Expresses $(-1)^j (1+t)^{n-3c} Q_c^{(2m)}(1, u_j)$ identically as $\frac{C_{m,c,j}}{2\pi}\int_{-\pi}^\pi F_1 F_2 F_3 B\,d\phi$ with exact phase cancellation into the nonnegative boundary kernel $B^{\text{pure}}(e^{i\phi}) = 4\sin^2((\phi + j\pi)/2) \ge 0$ for both parities of $c$.

11. **Lemmas BK-4.1–4.5 & Theorem BK-4.6 — Closed-Form Proof of Regime 1 ($\kappa \ge 2.5, n \ge 84$)** ([`PROOF.md` §4.3](bardakov-kozlovskaya/PROOF.md); [`algebra/verify_regime1_complete.py`](bardakov-kozlovskaya/algebra/verify_regime1_complete.py), [`verify/audit_bk_regime1.py`](bardakov-kozlovskaya/verify/audit_bk_regime1.py))  
    Zero-search 6-subcase elementary proof via Jensen midpoint separation (Lemma 4.1), Rouché isolation (Lemma 4.2), monotone Riemann sandwich for $L_{\text{sum}}$ (Lemma 4.3), off-boundary factor bounds (Lemma 4.4), and exact $B_{\text{avg}}$ (Lemma 4.5): worst closed-form margins **`+0.469913`** ($\kappa \ge 3.0$) and **`+0.074799`** ($\kappa \ge 2.5$), cross-checked across all $13{,}006$ finite-$n$ points (0 failures).

12. **Formal Zero-Crevice Partition of Regime 2 ($\kappa < 2.5, n \ge 132$)** ([`PROOF.md` §1.3](bardakov-kozlovskaya/PROOF.md), [`AUDIT.md` §2](bardakov-kozlovskaya/AUDIT.md))  
    Algebraically proves that Sub-regime 2A ($c\kappa^2 \ge 42$), Sub-regime 2C-Upgrade ($cq\kappa \le 0.90$), Piece G2 ($D=0$ Major Arcs), Piece G1 (Dirichlet Minor Arcs), and Piece G3 ($D \ne 0$ Major Arcs) cover every strictly interior Regime 2 point with strict seam overlaps and zero crevices (cross-checked across all $3{,}012$ Regime 2 points for $n \in \{132, \dots, 300\}$: $0$ uncovered points).

13. **Lemma BK-5.1 (R2.1) & Theorem BK-5.2 (R2.2) — Universal Regime 2 Perturbation & Sub-Regime 2A ($c\kappa^2 \ge 42$)** ([`PROOF.md` §5.1–5.2](bardakov-kozlovskaya/PROOF.md); [`analyst/certify_R2_1_and_2A_iv.py`](bardakov-kozlovskaya/analyst/certify_R2_1_and_2A_iv.py), [`verify/audit_bk_r2_1_2A_2C.py`](bardakov-kozlovskaya/verify/audit_bk_r2_1_2A_2C.py))  
    50-digit `iv` proof of universal Regime 2 stretched-exponential decay ($R_z/R_y \le e^{-37.8481/\kappa}$, $c_{\text{se}} \ge 3.16794$) and mirror-root pairing on $c\kappa^2 \ge 42$: $|\arg F_3^{\text{pure}}| \le 1.332945 < \pi/2$, $|\Delta_{123}| \le 8.246\times 10^{-6}$, `bdy_ratio <= 0.01499` $\implies$ **`>= +0.220613 > 0`**.

14. **Theorem BK-5.3 (R2.3) — Sub-Regime 2C-Upgrade ($c q \kappa \le 0.90$)** ([`PROOF.md` §5.3](bardakov-kozlovskaya/PROOF.md); [`analyst/certify_R2_3_iv.py`](bardakov-kozlovskaya/analyst/certify_R2_3_iv.py), [`verify/audit_bk_r2_1_2A_2C.py`](bardakov-kozlovskaya/verify/audit_bk_r2_1_2A_2C.py))  
    50-digit `iv` Pólya–Szegő positive-root ratio (`rat_main >= 0.473803`) + KL negative-binomial harmonics ($|R_{\text{non-}q}| \le 0.174420$) + conjugate-pair polylogarithms ($E_{\text{cubic},2C} \le 6.89\times 10^{-48}$) $\implies$ **`>= +0.441480`** ($q = 1$) and **`>= +0.266331`** ($q \ge 2$).

15. **Proposition BK-5.4, Lemmas G2.1–G2.7 & Theorem BK-5.5 — Piece G2 ($D = 0$ Major Arcs)** ([`PROOF.md` §5.4](bardakov-kozlovskaya/PROOF.md), [`ALGEBRA_ROUTE.md` §5](bardakov-kozlovskaya/ALGEBRA_ROUTE.md); [`algebra/certify_G2_D0.py`](bardakov-kozlovskaya/algebra/certify_G2_D0.py), [`verify/audit_bk_g2_checks.py`](bardakov-kozlovskaya/verify/audit_bk_g2_checks.py))  
    All-parity contour identity, explicit $Q_{\min} > 0.20$ cutoffs, 2nd-order Euler–Maclaurin curvature, $P^2$-monotone $E_{\text{phase}}$, $\pi$-line meromorphic residue + exact Cauchy–Riemann monotonicity ($\frac{d}{dy}V(x_*(y), y) = (U_x^2+U_y^2)/U_x < 0$), Chebyshev-mediant far-tail decay, and Pigeonhole-Annulus perturbation ($E_{\text{cubic}} \le 1.82\times 10^{-12}$ / $1.62\times 10^{-31}$) certify **`>= +0.216651`** ($q = 1, \kappa < 2.5$) and **`>= +0.00176348`** ($q \ge 2, q^2\kappa \le 2.0$); 28/28 finite-$n$ parity corners PASS (`min +0.97916`).

16. **Lemmas G1.1–G1.5 & Theorem BK-5.6 — Piece G1 (Dirichlet Minor Arcs)** ([`PROOF.md` §5.5](bardakov-kozlovskaya/PROOF.md), [`algebra/G1_NOTES.md`](bardakov-kozlovskaya/algebra/G1_NOTES.md); [`algebra/certify_G1_minor.py`](bardakov-kozlovskaya/algebra/certify_G1_minor.py), [`verify/audit_bk_g3_g1_checks.py`](bardakov-kozlovskaya/verify/audit_bk_g3_g1_checks.py))  
    Zero-cancellation Hubbard–Stratonovich duality kernel $H_D(t)$, continued-fraction Fibonacci ladder with 1-line $D_r^2$-cancellation, exact $m$-space difference ($S_{m\text{-space}} \le 0.277972$), and unit-circle perturbation certify $\sup|\Phi^{\text{ex}}| \le 1.536448 < \pi/2$ and **`I_full / L^1 >= +0.03355852 > 0`**; 27/27 finite-$n$ corners PASS (`min +0.88956`).

17. **Sub-Lemmas G3.1–G3.4 & Theorem BK-6.1 — Piece G3 ($D \ne 0$ Major Arcs)** ([`PROOF.md` §6](bardakov-kozlovskaya/PROOF.md), [`ANALYTIC_ROUTE.md` §4](bardakov-kozlovskaya/ANALYTIC_ROUTE.md); [`analyst/certify_G3_ABC_full_iv.py`](bardakov-kozlovskaya/analyst/certify_G3_ABC_full_iv.py), [`analyst/certify_G3_D_full_iv.py`](bardakov-kozlovskaya/analyst/certify_G3_D_full_iv.py), [`verify/audit_bk_g3_post106.py`](bardakov-kozlovskaya/verify/audit_bk_g3_post106.py))  
    Exact $(1-w)$-free IBP kernels $W_{1,\text{fin}}, W_{2,\text{fin}}, W_{3,\text{fin}}$, Harmonicity identity, zero-additive-constant shifted-line Poisson + 4-pole residue identity, horizontal-saddle shift, and Pigeonhole-Annulus perturbation ($E_{\text{cubic,G3}} \le 2.80\times 10^{-6}$) certify **`I_full / L^1 >= +0.0545 > 0`** globally across all 5 exhaustive sub-cases (**A1 `+0.08207`**, **A2 `+0.40597`**, **B `+0.13156`**, **C `+0.07706`**, **D1 `+0.1646`**, **D2 `+0.0545`**, **D3 `+0.0784`** across all 37 $\tau$-slices); 32/32 extreme finite-$n$ corners PASS (`min +0.54579`).

18. **Exact Structural Disproof BK-D1 — Failure of Naive Cross-Degree Interlacing** ([`verify/audit_bk_40.py`](bardakov-kozlovskaya/verify/audit_bk_40.py), board `#42`)  
    Naive cross-degree slice interlacing $G_{n-2, c} \ll G_{n, c}$ is **FALSE**: explicit counterexamples occur at $(n, c) = (16, 4)$, $(26, 7)$, and $(13, 3)$, ruling out naive degree-induction interlacing and motivating the universal sample-point 1D circle integral architecture.

---

## Part II. Gessel–Lin–Zeng Conjecture 15 — Sharp Obstructions & General Real-Stability Theorems

Full manuscripts: [`gessel-lin-zeng/AUDIT.md`](gessel-lin-zeng/AUDIT.md) · [`gessel-lin-zeng/ALGEBRA_ROUTE.md`](gessel-lin-zeng/ALGEBRA_ROUTE.md) · [`gessel-lin-zeng/ANALYTIC_ROUTE.md`](gessel-lin-zeng/ANALYTIC_ROUTE.md)  
*(Honest status: Conjecture 15 itself remains **OPEN** for $4 \le i \le k - 3$.)*

### II.1 Proved General Polynomial & Real-Stability Theorems

1. **Theorem GLZ-1 — Jadav's Ma–Wang Sharpness Theorem (100% Sharpness across ALL $4 \le i \le k - 3$)** ([`ALGEBRA_ROUTE.md` §2](gessel-lin-zeng/ALGEBRA_ROUTE.md), [`AUDIT.md` §3](gessel-lin-zeng/AUDIT.md); [`verify/audit_31_32_33.py`](gessel-lin-zeng/verify/audit_31_32_33.py), [`verify/verify_clean_laguerre_certs.py`](gessel-lin-zeng/verify/verify_clean_laguerre_certs.py))  
   Ma and Wang's Theorem 2 ([arXiv:2610.03111](https://arxiv.org/abs/2610.03111), proving common interlacers of $\{A_{k,S} : |S|=i\}$ for $i \in \{1,2,3,k-2,k-1\}$) is **100% sharp**: non-real nonnegative linear combinations $R(x) = \lambda A_{k,\text{late}}(x) + A_{k,\text{early}}(x)$ exist in every remaining boundary family, certified by exact integer/Fraction Laguerre violations $R'(x)^2 - R(x)R''(x) < 0$:
   - **$i = 4$** ($k = 30$ passes exact Cauchy index; first failure **$k = 31$**, deg $57$): $932\,A_{31,\{28..31\}} + A_{31,\{1..4\}}$ at **$x = -23$** $\implies -9.065320\times 10^{279} < 0$.
   - **$i = 5$** ($k = 15$ passes exact Cauchy index; first failure **$k = 16$**, deg $26$): $269\,A_{16,\{12..16\}} + A_{16,\{1..5\}}$ at **$x = -11$** $\implies -2.317817\times 10^{96} < 0$.
   - **$i = k - 3$** ($k = 27$ passes exact Cauchy index; first failure **$k = 28$**, deg $30$): $272\,A_{28,\{4..28\}} + A_{28,\{1..25\}}$ at **$x = -73/10$** $\implies -5.025603\times 10^{111} < 0$ (and integer point $x = -8$ at $k = 29$: $349\,A_{29,\{4..29\}} + A_{29,\{1..26\}} \implies -4.897981\times 10^{118} < 0$).

2. **Theorem GLZ-2 — Jadav's Cubic Real-Stability Addition Theorem** ([`ALGEBRA_ROUTE.md` §2](gessel-lin-zeng/ALGEBRA_ROUTE.md), [`AUDIT.md` §2](gessel-lin-zeng/AUDIT.md); [`verify/audit_25.py`](gessel-lin-zeng/verify/audit_25.py))  
   Let $v_0, v_1, v_2, v_3, v_4 \in \mathbb{R}^2$ lie in strict open-half-plane order ($W(v_r, v_s) > 0$ for all $0 \le r < s \le 4$). If the cubics $P_1(y) = v_0 - 3v_1 y + 3v_2 y^2 - v_3 y^3$ and $P_2(y) = v_1 - 3v_2 y + 3v_3 y^2 - v_4 y^3$ both have PSD $3 \times 3$ Bezoutiants (allowing mixed-sign $y$-roots), then $P_1(y) + P_2(y)$ also has PSD Bezoutiant. Moreover, $d = 3$ is the unique maximum degree with $\min(p, d + 1 - p) \le 2$, and the Wronskian cone hypothesis is sharp (witness $L = (-3000, 1, 10, 1, -3000)$ gives 1 real root by exact Sturm).

3. **Theorem GLZ-3 — Sharpest Quantitative Plücker Quad Theorem B** ([`ALGEBRA_ROUTE.md` §2](gessel-lin-zeng/ALGEBRA_ROUTE.md), [`AUDIT.md` §2](gessel-lin-zeng/AUDIT.md); [`verify/audit_18_20_21.py`](gessel-lin-zeng/verify/audit_18_20_21.py), [`verify/audit_25.py`](gessel-lin-zeng/verify/audit_25.py))  
   For differentiable real functions $A, B, C, E$ on $I$ with $W_{AB}, W_{BC}, W_{CE}, W_{AC}, W_{BE} > 0$, $W_{AE} \ge 0$, and $J(A, B, C), J(B, C, E) \ge 1 + \varepsilon$ ($\varepsilon \ge 0$), unconditionally $W(A+B, C+E) > 0$ and (with **zero upper bound** on the Plücker ratio $R = W_{AB}W_{CE}/W_{BC}^2$):
   $$J(A + B,\; B + C,\; C + E) \;\ge\; \frac{4(1 + \varepsilon)^2}{(2 + \varepsilon)^2} \;\ge\; 1.$$

4. **Theorem GLZ-4 — Universal Cubic Wronskian Log-Concavity** ([`AUDIT.md` §2](gessel-lin-zeng/AUDIT.md); [`verify/audit_18_20_21.py`](gessel-lin-zeng/verify/audit_18_20_21.py))  
   Exact global sum-of-squares proof that every real-stable cubic $P(t, y) = \sum_{j=0}^3 p_j(t) y^j$ satisfies $W(p_0, p_1) W(p_2, p_3) \le W(p_1, p_2)^2$ on all of $\mathbb{R}$ (constant $9$ in binomial form $A + 3By + 3Cy^2 + Ey^3$).

5. **Theorem GLZ-5 — Zero-Loss $D_m$ Preservation ($m \ge 1$)** ([`ANALYTIC_ROUTE.md` §1](gessel-lin-zeng/ANALYTIC_ROUTE.md), [`AUDIT.md` §2](gessel-lin-zeng/AUDIT.md); [`verify/audit_18_20_21.py`](gessel-lin-zeng/verify/audit_18_20_21.py))  
   By exact Borcea–Brändén symbol partial fraction $w(t) = -(1 + \frac{n}{m})t - \frac{n(m-1)}{m^2} + \frac{n(m-1)}{m^2(1+mt)}$, $D_m = t(1+t)\partial_t + 1 + mt$ preserves real-stability (and hence $q$-Turán $J \ge q$ on $\mathbb{R}$ and cubic real-stability) with zero loss **if and only if** $m \ge 1$.

6. **Theorem GLZ-6 — Theorem A (Turán $\implies$ Interlacing Preservation)** ([`AUDIT.md` §2](gessel-lin-zeng/AUDIT.md); [`verify/audit_18_20_21.py`](gessel-lin-zeng/verify/audit_18_20_21.py))  
   For $p_2 \prec p_1 \prec p_0$ in $(-1, 0)$ with $\deg p_0 = \deg p_1 + 1 = \deg p_2 + 2$, positive leading coefficients, and $p_1^2 \ge p_0 p_2$ on $(-1, 0)$, we have $T p_2 + w p_1 \prec T p_1 + w p_0$ for any $T = t(1+t)\partial_t + a + bt$ ($0 < a < b$) and $w > 0$ (constant $1$ is sharp).

7. **Theorem GLZ-7 — Theorem C (Pole-Gap Wronskian Bound)** ([`AUDIT.md` §2](gessel-lin-zeng/AUDIT.md))  
   For $\deg f_0 = \deg f_1 + 1$, $f_1 \prec f_0$ in $(-1, 0)$, and $f_0 \prec D_r f_1$ with $r > 1 - \deg f_0$ and $j > r$, $0 < W(f_0, t f_1) / W(f_0, D_j f_1) < 1/(j - r)$ on $\mathbb{R}$.

8. **Theorem GLZ-8 — Pure-Power Sturm–Liouville Representation & Two-Regime Turán Dichotomy** ([`ANALYTIC_ROUTE.md` §1](gessel-lin-zeng/ANALYTIC_ROUTE.md))  
   In $s$-space ($x = -e^s$), $\Psi_k(s, y) = \mathcal{T}_y^k(e^s/(1+e^s))$ with symmetric Sturm–Liouville operator $\mathcal{T}_y = \frac{1}{1+e^s}\partial_s(\partial_s + y)$; and $(D_m f)^2 - (D_m g)(D_m h) = A\Phi^2 + B\Phi + C$ has $\Phi$-discriminant $B^2 - 4AC = t^2(1+t)^2 W_{gh}^2(1 - J_{\text{raw}})$, forcing $A'(t) \ne 0$ across any connected interval where $A, C_0 > 0$ and $J_{\text{raw}} < 1$.

---

### II.2 Exact-Arithmetic Disproofs & Certified Obstructions (GLZ)

9. **Disproof GLZ-D1 — Failure of Binomial-Weighted Bivariate Real Stability $B_k(x, y)$ at $k = 10$** ([`verify/verify_msg13.py`](gessel-lin-zeng/verify/verify_msg13.py))  
   Exact `Fraction` Rayleigh determinant at $(0, -3/4)$ equals $-2649516238079970755/4294967296 < 0$; exact integer Sturm on $B_{10}(s, -2531/1000 + (130151/100)s)$ gives squarefree degree $19$ with only **$17$ real roots**.

10. **Disproof GLZ-D2 — Failure of $2\times 2$ Minor $x$-Positivity & Global-$\mathbb{R}$ Turán Positivity** ([`verify/verify_deep_positivity.py`](gessel-lin-zeng/verify/verify_deep_positivity.py))  
    - First negative $x$-coefficient: $[x^1](A_{19,8}^2 - A_{19,7}A_{19,9}) = -124{,}876{,}023{,}597{,}181{,}272 < 0$.
    - First pointwise failure on $\mathbb{R}$: $A_{35,17}(7/10^{17})^2 - A_{35,16}(7/10^{17})A_{35,18}(7/10^{17}) \approx -2.376737\times 10^{17} < 0$.

11. **Disproof GLZ-D3 — Failure of Universal Operator Products (Free Parameters & All-Integer $b_j \ge 2$)** ([`verify/audit_18_20_21.py`](gessel-lin-zeng/verify/audit_18_20_21.py), [`verify/audit_31_32_33.py`](gessel-lin-zeng/verify/audit_31_32_33.py))  
    - Free $(a_j, b_j, w_j)$ with $0 < a_j < b_j, w_j > 0$: exact Sturm failure at **$k = 9$** (`q9_3`, deg 6, squarefree, 4 real roots).
    - Fixed $a_j = w_j = 1$, all-integer $b_j \ge 2$ (single-op): exact finite witness ($K = 1000$ copies of $b = 2$, $b_{1001} = 196{,}200{,}000{,}000{,}000$, $b_{1002} = 181{,}330{,}000{,}000{,}000$, $i = 960$) breaks interlacing (exact Cauchy index $40 \ne 42$; exact `Fraction` Wronskian $-3.902371\times 10^{30} < 0$).
    - All-integer $b_j \ge 2$ (double-op): exact SymPy scaling limit $\mathcal{M}_{\text{dbl}}^{(w)}(\Phi) < 0$ for all $w \ge 5$, and the outer $D_{b_2}$ fails to repair the broken interlacing.

12. **Disproof GLZ-D4 — Impossibility of Linear Gauges & Auxiliary Interlacing Cones** ([`ALGEBRA_ROUTE.md` §3](gessel-lin-zeng/ALGEBRA_ROUTE.md), [`AUDIT.md` §3](gessel-lin-zeng/AUDIT.md); [`verify/audit_25.py`](gessel-lin-zeng/verify/audit_25.py))  
    5-line algebraic proof that no linear gauge eliminates the Stage 3 shift $D_{k+1}^i \mapsto D_{k+2}^i$; half-plane Wronskian $W(E_{k,0}, E_{k,9})(-19/20) \approx -1.4255\times 10^{24} < 0$ at $k = 10$; $Q_{k,q}$ proper position fails for all $q > 0$; and black-box 1-step Stage 3 fails even at GLZ's boundary scale $\lambda = 3/j$.

---

*See also: [`EQUATIONS.md`](EQUATIONS.md) · [`JOURNEY.md`](JOURNEY.md) · [`bardakov-kozlovskaya/README.md`](bardakov-kozlovskaya/README.md) · [`gessel-lin-zeng/README.md`](gessel-lin-zeng/README.md) · [Root `README.md`](README.md).*
