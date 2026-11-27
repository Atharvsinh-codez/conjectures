# Independent Adversarial Referee Audit Report: Resolution of Bardakov–Kozlovskaya Conjecture 3.3(2)

**Auditor:** `@verifier` (Independent Adversarial Referee)  
**Target Conjecture:** V. G. Bardakov and T. A. Kozlovskaya, *Multivalued groups and Newton polyhedron* ([arXiv:2305.07261](https://arxiv.org/abs/2305.07261), Conjecture 3.3(2); `1337isnot1337/solved-conjectures` Issue #174; MathDB `#358383`)  
**Strengthened Theorem Audited:** **The Jadav Sign-Law & Real-Rootedness Theorem** (Theorem 1.1 & Corollary 1.2 of [`PROOF.md`](PROOF.md))  
**Audit Status:** **COMPLETE — ALL 10 OF 10 PROOF COMPONENTS & FINALIZED [`PROOF.md`](PROOF.md) 100% INDEPENDENTLY VERIFIED & SIGNED OFF (ZERO REMAINING DEFECTS).**

---

## 1. Executive Summary & Audit Standards

Every proof claim, algebraic identity, analytic inequality, seam overlap, and computational certificate in the proof of Bardakov–Kozlovskaya Conjecture 3.3(2) was audited independently by `@verifier` under strict adversarial standards:

1. **Zero unverified extrapolation:** Every infinite-parameter theorem ($n \to \infty$) is proved uniformly over all $(n, c, j)$ with zero loops over unbounded integer parameters and explicit closed-form asymptotic tails.
2. **Exact arithmetic & 50-digit outward-rounded interval arithmetic:** The finite base ($n = 2, 4, \dots, 130$) uses exact arbitrary-precision integer arithmetic (`int`); all Regime 2 computational proofs use 50-digit outward-rounded `mpmath.iv` (`iv.dps = 50`) with pure `mpf` endpoints (zero float64 narrowing inside `iv.mpf([...])`), zero hardcoded uncertified constants, and rigorous enclosures of all special constants (e.g., $\zeta(3)$ via $M = 2000$ partial sum + integral sandwich).
3. **Adversarial counterexample probing:** Every subtle analytic bound was independently stress-tested against exact finite-$(n, c, j)$ evaluations across extreme corners (including $j = c + 1$, $D = \pm 1, \pm 2, \pm 3$, $r_0 \in \{1, 2, 3, 5\}$, $q = 1$ vs $q \ge 2$, and all parities of $c, d, q, p_0$). Across four major audit rounds (`#85`, `#96`, `#98`, `#106`), **17 genuine mathematical/soundness bugs** were caught with concrete counterexamples, reported to the owning peers, patched rigorously, and re-verified to **0 violations**.

---

## 2. Master Zero-Crevice Partition & Component Sign-Off Table

For every even $n = 2m \ge 2$ and every admissible triple $(a, b, c) \in \mathbb{Z}_{\ge 0}^3$ with $a + 2b + 3c = 2m$, the proof partitions the parameter space with **zero crevices** (explicitly verified both algebraically and computationally across all $3{,}012$ Regime 2 points for $n \in \{132, 134, \dots, 300\}$):

- **Separation parameter consistency:** Every Regime 1 and Regime 2 script uses the exact same trigonometric separation parameter $\kappa := \kappa_{\text{trig}}(n, c, j) = 2\pi\frac{\sin\psi_j\sin\alpha_1}{\sin\theta_c}$ (where $\theta_c = \frac{c\pi}{n}$, $\psi_j = \frac{j\pi}{n}$, $\alpha_1 = \frac{(j+c)\pi}{n}$, $d = \gcd(c, j) \ge 1$, $q = c/d \ge 1$, $X = q^2\kappa$).
- **Boundary & seam completeness:**
  - $n \le 130$ vs $n \ge 132$: exact complement on even integers.
  - $c \in \{0, 1\}$, $a = 0$, $b = 0$ (covering all $c > \lfloor(n-3)/3\rfloor$) vs interior $2 \le c \le \lfloor(n-3)/3\rfloor, j \in \mathcal{J}_{n,c}$: exact complement.
  - $\kappa \ge 2.5$ (Regime 1) vs $\kappa < 2.5$ (Regime 2): exact complement.
  - $c\kappa^2 \ge 42$ (Sub-regime 2A) vs $c\kappa^2 < 42$: exact complement.
  - $cq\kappa \le 0.90$ (Sub-regime 2C-Upgrade) vs $cq\kappa > 0.90$: overlaps Piece G2 ($cq\kappa \ge 0.85$) and Piece G1 ($X \ge 0.88$ for $d = 1$).
  - Inside $\{c\kappa^2 < 42,\; cq\kappa > 0.90\}$:
    - $q = 1 \implies d = c \ge 2, \kappa < 2.5 = K_{\text{G2},q=1} \implies$ **Piece G2** (zero-crevice seam closed in `#102`–`#103`).
    - $q \ge 2, d \ge 2, q^2\kappa \le 2.0 \implies$ **Piece G2**.
    - $q \ge 2$ and ($d = 1$ or $q^2\kappa > 2.0$): either $\exists$ coprime $(p_0, r_0)$ with $D = r_0 j - c p_0 \ne 0$ satisfying $\tilde{c} \ge 10 D^2, r_0\kappa \le 0.45, r_0^3 c\kappa^2 \le 60, \tilde{c}\tau^2/D^4 \ge 45 > 40 \implies$ **Piece G3** (strict overlap $\tilde{c}\tau^2/D^4 \in [40, 45]$); or none exists $\implies$ **Piece G1**.
  - Inside Piece G3: $|D| \ge 2$ split seamlessly at $s = 2.0$ (**G3-A1** $s \in (0, 2]$ vs **G3-A2** $s \in [2, \infty)$, with $|D|\ge 2$ A1 split at $|\delta_{\text{fin}}| = 0.08$); $|D| = 1$ split seamlessly across $\tau \in (0, 1.05]$ (**G3-C**), $\tau \in [1.05, 40.0]$ (**G3-D**: 37 seamless $\tau$-slices across D1 $[1.05, 1.50]$, D2 $[1.50, 2.80]$, D3 $[2.80, 40.0]$), and $\tau \ge 40.0$ (**G3-B**).

| # | Proof Component | Exact Domain | Independent Verification Script(s) | Certified Margin / Result | Sign-Off |
|---|---|---|---|---|---|
| 1 | **Finite Base ($n \le 84$)** | Even $n = 2..84$ (42 degrees) | [`verify/engine.py`](verify/engine.py), [`verify/stress_to_84.py`](verify/stress_to_84.py) | **9,464** exact integer coeffs; 0 zeros, 0 sign violations; SymPy & Wendt cross-checked | **SIGNED OFF** (`#39`) |
| 2 | **Finite Base ($n = 86..130$)** | Even $n = 86..130$ (23 degrees) | [`verify/fast_certify_140.py`](verify/fast_certify_140.py) | **23,953** exact integer coeffs (combined **$33{,}417$**); bit-identical to `engine.py` on $n=2..30$; 0 zeros, 0 sign violations | **SIGNED OFF** (`#57`) |
| 3 | **Resultants & 4 Boundaries** | All even $n \ge 2$: $c=0, a=0, b=0, c=1$ | [`verify/audit_bk_40.py`](verify/audit_bk_40.py), [`verify/audit_bk_43.py`](verify/audit_bk_43.py), [`verify/audit_bk_50.py`](verify/audit_bk_50.py) | Exact closed-form proofs of Theorems 3.1–3.5 verified symbolically & numerically | **SIGNED OFF** (`#42`, `#45`, `#58`) |
| 4 | **1D Circle Integral & Regime 1** | $n \ge 84$, interior, $\kappa \ge 2.5$ | [`verify/audit_bk_regime1.py`](verify/audit_bk_regime1.py), [`algebra/verify_regime1_complete.py`](algebra/verify_regime1_complete.py) | 50-dps `iv` margins **`+0.469913`** ($\kappa \ge 3.0$) & **`+0.074799`** ($\kappa \ge 2.5$); 13,006 finite-$n$ cross-checks (0 failures) | **SIGNED OFF** (`#59`, `#61`) |
| 5 | **Lemma R2.1 & Sub-regime 2A** | Regime 2, $c\kappa^2 \ge 42, \kappa < 2.5, n \ge 132$ | [`analyst/certify_R2_1_and_2A_iv.py`](analyst/certify_R2_1_and_2A_iv.py), [`verify/audit_bk_r2_1_2A_2C.py`](verify/audit_bk_r2_1_2A_2C.py) | `B_2A <= 1.332945 < pi/2`, `Delta_123 <= 8.246e-6`, `bdy_ratio <= 0.01499` $\implies$ **`>= +0.220613`** (exit 0) | **SIGNED OFF** (`#71`, `#74`) |
| 6 | **Sub-regime 2C-Upgrade** | Regime 2, $cq\kappa \le 0.90, \kappa < 2.5, n \ge 132$ | [`analyst/certify_R2_3_iv.py`](analyst/certify_R2_3_iv.py), [`verify/audit_bk_r2_1_2A_2C.py`](verify/audit_bk_r2_1_2A_2C.py) | `rat_main >= 0.473803`, `|R_nonq| <= 0.174420`, `E_cubic <= 6.889e-48` $\implies$ **`>= +0.441480`** ($q=1$), **`>= +0.266331`** ($q\ge 2$) (exit 0) | **SIGNED OFF** (`#77`) |
| 7 | **Piece G2 ($D=0$ Major Arcs)** | Regime 2, $c\kappa^2 < 42, cq\kappa \ge 0.85, d \ge 2$, $\kappa < 2.5$ ($q=1$) or $q^2\kappa \le 2.0$ ($q\ge 2$) | [`algebra/certify_G2_D0.py`](algebra/certify_G2_D0.py), [`verify/audit_bk_g2_checks.py`](verify/audit_bk_g2_checks.py) | $q=1$: **`>= +0.216651`** ($E_{\text{cubic}} \le 1.82\times 10^{-12}$); $q\ge 2$: **`>= +0.00176348`** ($E_{\text{cubic}} \le 1.62\times 10^{-31}$); 28/28 parity corners PASS (`min +0.97916`); 200.6s, exit 0 | **SIGNED OFF** (`#90`, `#103`) |
| 8 | **Piece G1 (Dirichlet Minor Arcs)** | Regime 2, $c\kappa^2 < 42, q \ge 2$, $X \ge 0.88$ ($d=1$) or $X \ge 2.0$ ($d\ge 2$), outside G3 ($C_\tau = 45$) | [`algebra/certify_G1_minor.py`](algebra/certify_G1_minor.py), [`verify/audit_bk_g3_g1_checks.py`](verify/audit_bk_g3_g1_checks.py) | `max |Phi^ex| <= 1.536448 < pi/2`, `Delta_123 <= 5.15e-6`, `bdy_ratio <= 0.000778` $\implies$ **`>= +0.03355852`**; 27/27 corners PASS (`min +0.88956`); 755.4s, exit 0 | **SIGNED OFF** (`#101`) |
| 9 | **Piece G3 Cases A1, A2, B, C** | Box G3: $|D|\ge 2$ (A1, A2); $|D|=1, \tau\ge 40$ (B); $|D|=1, \tau\le 1.05$ (C) | [`analyst/certify_G3_ABC_full_iv.py`](analyst/certify_G3_ABC_full_iv.py), [`analyst/verify_G3_empirical_corners.py`](analyst/verify_G3_empirical_corners.py), [`verify/audit_bk_g3_post106.py`](verify/audit_bk_g3_post106.py) | **A1 `>= +0.08207`** (`A1a: +0.08207`, `A1b: +0.11879`), **A2 `>= +0.40597`**, **B `>= +0.13156`**, **C `>= +0.07706`** (`E_cubic <= 2.80e-6`; 29.78s, exit 0); 32/32 empirical corners PASS (`min +0.54579`) | **SIGNED OFF** (`#124`) |
| 10 | **Piece G3 Case D** | Box G3: $|D|=1, \tau \in [1.05, 40.0]$ (37 slices, 6 seamless CLI parts) | [`analyst/certify_G3_D_full_iv.py`](analyst/certify_G3_D_full_iv.py), [`verify/audit_bk_g3_post106.py`](verify/audit_bk_g3_post106.py) | **D1 `>= +0.1646`** (58.96s), **D2 `>= +0.0545`** (960.43s), **D3 `>= +0.0784`** (`D3a1: +0.0784`, `D3a2: +0.0840`, `D3b1: +0.0952`, `D3b2: +0.2055`); all 6/6 parts EXIT 0 (3054.59s total, `--workers 3`) | **SIGNED OFF** (`#124`, `#132`) |

---

## 3. Detailed Component-by-Component Verification Log

### 3.1 Finite Exact Integer Base ($n = 2, 4, \dots, 130$) — SIGNED OFF (`#39`, `#57`)
- **Primary engine ([`verify/engine.py`](verify/engine.py)):** Computes exact integer coefficients $c_{a,b,c}^{(2m)}$ via Newton power sums $S_r(x, y) = n \sum_{a=0}^r \binom{nr}{na} x^{r-a} y^a$, Newton's identities for $E_r(x, y)$, and exact zero-remainder peeling in $e_3 = s_2 z$. Cross-checked identically against direct SymPy resultants ($n \le 8$) and the Buchstaber–Kornev Wendt circulant determinant ($n \le 6$). Certified all $42$ even degrees $n = 2, 4, \dots, 84$: **$9{,}464$** coefficients, $0$ zeros, $0$ sign-law violations (`#39`).
- **Fast $S_3$-Kronecker engine ([`verify/fast_certify_140.py`](verify/fast_certify_140.py)):** Exploits $S_3$-symmetry ($k_2 + k_3 \le \lfloor 2n/3 \rfloor$) and balanced base-$2^B$ Kronecker polynomial multiplication with explicit carry-overflow assertions and redundant overlap symmetry checks. Cross-verified to bit-for-bit equality against `verify/engine.py` across all $n \in \{2, 4, \dots, 30\}$. Certified all $23$ even degrees $n = 86, 88, \dots, 130$: **$23{,}953$** coefficients, $0$ zeros, $0$ sign-law violations (`#57`).
- **Combined exact integer base ($n = 2, 4, \dots, 130$):** **$33{,}417$** coefficients across $65$ even degrees, $0$ zeros, $0$ sign-law violations.

### 3.2 Exact Structural Resultants & All Four Boundary Families ($n \ge 2$) — SIGNED OFF (`#42`, `#45`, `#58`)
Audited line-by-line and checked symbolically/numerically in [`verify/audit_bk_40.py`](verify/audit_bk_40.py), [`verify/audit_bk_43.py`](verify/audit_bk_43.py), and [`verify/audit_bk_50.py`](verify/audit_bk_50.py):
- **Theorems 3.1–3.2 (Depressed-Cubic & Half-Size Lucas–Chebyshev Resultants):** Normalization $+1$, degree count, and specialization at $p = 0$ verified exact.
- **Theorem 3.3 ($c = 0$):** $Q_0^{(2m)}(e_1, e_2) = (e_1^2 - 4e_2)^m \implies c_{2m-2b, b, 0}^{(2m)} = (-4)^b\binom{m}{b}$. Exact.
- **Theorem 3.4 ($a = 0$ and $b = 0$):** Exact trigonometric derivative identities for $L_{2m}(1/(4\cos^2\theta))$ and strict sign alternation at the $\lfloor m/3 \rfloor + 1$ (resp. $\lfloor 2m/3 \rfloor + 1$) critical points prove all real roots are simple and strictly positive; Vieta's formulas lock the exact sign law for all even $n \ge 2$.
- **Theorem 3.5 ($c = 1$):** Rolle real-rootedness of $F_r(u)$, palindromic factorization of $A_{2m-1}^{(2m)}(u)/(1+u)$, Maclaurin's inequality reducing all $b$ to $\gamma_{m-1} > 4^{m-1}$, and the exact alternating trigonometric sum $\gamma_{m-1} = \frac{1}{2m}\sum_{k=0}^{m-1}(-1)^k \frac{\cos\theta_k}{\sin^{2m+1}\theta_k}$ independently verified for all $m \ge 2$.

### 3.3 1D Circle Integral & Closed-Form Regime 1 ($\kappa \ge 2.5, n \ge 84$) — SIGNED OFF (`#59`, `#61`)
Independently audited in [`verify/audit_bk_regime1.py`](verify/audit_bk_regime1.py) and re-executed in 50-digit `mpmath.iv`:
- Theorem 4.1 (interleaving reduction on $\mathcal{J}_{n,c}$ for both parities of $c$) and Theorem 4.2 (exact 4-factor 1D circle integral & exact sign factor $(-1)^j$ for both parities of $c$) verified identically (and to $< 10^{-50}$ relative error in 80-digit `mpmath`).
- Lemmas 4.1–4.5 (Jensen midpoint separation, Rouché isolation, monotone Riemann sandwich for $L_{\text{sum}}$, off-boundary factor bounds, exact $B_{\text{avg}}$) and all 6 subcases of Theorem 4.6 verified with zero search: worst closed-form margin **`+0.469913`** at $\kappa_0 = 3.0$ and **`+0.074799`** at $\kappa_0 = 2.5$ (`Case 1A`), plus 13,006/13,006 finite-$n$ cross-checks with 0 failures.

### 3.4 Regime 2 Core Lemmas R2.1, R2.2 (Sub-regime 2A) & R2.3 (Sub-regime 2C-Upgrade) — SIGNED OFF (`#71`, `#74`, `#77`)
Independently audited in [`verify/audit_bk_r2_1_2A_2C.py`](verify/audit_bk_r2_1_2A_2C.py) and executed directly in 50-digit `mpmath.iv`:
- **Lemma R2.1 & Sub-regime 2A ([`analyst/certify_R2_1_and_2A_iv.py`](analyst/certify_R2_1_and_2A_iv.py)):** Universal $\alpha_1 < 0.409284$, $C_{\exp} \ge 37.8481$, $c_{\text{se}} \ge 3.16794$; exact mirror-root pairing phase bound $B_{2A} \le 1.332945 < \pi/2$, $|\Delta_{123}| \le 8.246\times 10^{-6}$, oscillation constant $C_{\text{osc,tot}} \le 21.175 \ll C_{\text{exp,pert}} \ge 37.561$, `bdy_ratio <= 0.01499` $\implies$ **`>= +0.220613 > 0`** (exit code 0, 0.3s).
- **Sub-regime 2C-Upgrade ([`analyst/certify_R2_3_iv.py`](analyst/certify_R2_3_iv.py)):** Pole-free contour-shift $W_0 < 0.836$ ($q=1$) / $< 0.832$ ($q\ge 2$), $E_{\text{cubic},2C} \le 6.889\times 10^{-48}$, Pólya–Szegő `rat_main >= 0.473803`, KL negative-binomial non-$q$ harmonics $|R_{\text{non-}q}| \le 0.174420$, conjugate-pair $E_{\text{phase}} \le 0.031812$ ($q=1$) / $\le 0.014110$ ($q\ge 2$) $\implies$ **`>= +0.441480`** ($q=1$) and **`>= +0.266331`** ($q \ge 2$) (exit code 0, 13.1s).

### 3.5 Piece G2 ($D = 0$ Major Arcs) & Zero-Crevice Seam — SIGNED OFF (`#90`, `#103`)
Audited across all 7 lemmas of [`algebra/certify_G2_D0.py`](algebra/certify_G2_D0.py) and cross-checked in [`verify/audit_bk_g2_checks.py`](verify/audit_bk_g2_checks.py):
- **Bugs caught & fixed (`#85` $\to$ `#89`–`#90`):** Caught 3 bugs in `cert_Ep_over_sqrtA_iv` (missing `A_pk = max(A0, pi*q_min)` in `ksqA`/`ksq_tau`, wrong $A$-powers in `rhs_w`, and `1/A0` vs `1/sqrt(A0*T_e)` in `rhs_w1`); all 3 fixed cleanly and re-verified.
- **Zero-crevice $q = 1, \kappa \in (2.0, 2.5)$ seam (`#102`–`#103`):** `@algebraist` identified and closed the $q = 1, \kappa \in (2.0, 2.5)$ seam by widening $q = 1$ to `K_G2_q1 = 2.5` (`alp1_q1 = 0.065*pi`, discrete $c \in \{2..10\}$ + automatic $\kappa < \sqrt{42/11} < 1.9541$ for $c \ge 11$).
- **Independent execution of [`algebra/certify_G2_D0.py`](algebra/certify_G2_D0.py) (exit code 0, 200.6s):**
  - **$q = 1$ ($\kappa < 2.5$):** $H_1^+(\Lambda_*) \ge 0.68924$ ($c \le 10$), $\ge 0.62020$ ($c \ge 11$, $H_1^+(\tilde\Lambda) \ge 0.62018$), $E_{\text{cubic}} \le 1.8182\times 10^{-12}$, **`I_full / L1 >= +0.216651 > 0`**.
  - **$q \ge 2$ ($q^2\kappa \le 2.0$):** $H_q^+(\tilde\Lambda) \ge 0.59425$, $|\mathcal{R}_{\text{non-}q}| \le 0.115416$ ($|\psi|\le 1.25$) / $\le 0.333745$ ($|\psi|\le 2\pi/3$), $E_{\text{cubic}} \le 1.6158\times 10^{-31}$, **`I_full / L1 >= +0.00176348 > 0`**.
  - **Part B:** 28/28 exact finite-$n$ corners covering all $(d, q)$ parities PASS (`min I_pure/L1 = +0.97916`).

### 3.6 Piece G1 (Dirichlet Minor Arcs) — SIGNED OFF (`#101`)
Audited line-by-line across all 2,968 lines of [`algebra/certify_G1_minor.py`](algebra/certify_G1_minor.py) + [`algebra/G1_NOTES.md`](algebra/G1_NOTES.md), cross-checked in [`verify/audit_bk_g3_g1_checks.py`](verify/audit_bk_g3_g1_checks.py):
- **Bugs caught & fixed (`#98` $\to$ `#100`–`#101`):** Caught 5 edge-case bugs: (1) `D_rat = max(_FIB_IDX[L+2], D2_eff_int // dm_idx)` failing for $d \in \{3, 4\}$ (`dm_idx=2`) and unbounded $d \ge 5$; (2) `eval_B2_ladder` using `12/11` instead of `1` for unbounded `Dt >= 12` in B2b; (3) `earl_split[12]` having increasing $t$-ratios `Dt/(a*Dt+b)` for unbounded `Dt >= 12`; (4) `A_tail` needing `(Tb/Ta)^2` instead of `Tb/Ta`; (5) `T_lg2_p` failing at $(132, 2, 3)$ because $n\kappa_{\text{trig}} = 147.638 < 15\pi^2$ (fixed uniformly to exact $\kappa_{\text{trig}}$ everywhere).
- **Independent execution of [`algebra/certify_G1_minor.py`](algebra/certify_G1_minor.py) (exit code 0, 755.4s):**
  - **G1.2:** `B1(d=1) <= 1.25848`, `B1(d>=2) <= 1.13689`, `B2a <= 1.25112`, `B2b <= 1.23337`, `Case A <= 0.915136`.
  - **G1.3:** `S_mspace <= 0.277972`.
  - **G1.4:** Global phase bound **`max |Phi^ex(phi; 1)| <= 1.536448 < pi/2`** $\implies$ `I_pure / L1(1) >= +0.03434128`.
  - **G1.5:** `Delta_123 <= 5.1501e-6`, `C_osc_tot <= 16.913 << C_exp_pert >= 37.5613` (`T_lg2 >= 8.88518`), `bdy_ratio <= 0.00077761` $\implies$ **`I_full / L1(1) >= +0.03355852 > 0`**.
  - **Part B:** 27/27 exact finite-$n$ Minor-Arc corners PASS (`max |Phi^ex| = 1.1542`, `min I_pure/L1 = +0.88956`).

### 3.7 Piece G3 ($D \ne 0$ Major Arcs) — 100% SIGNED OFF (ALL CASES A1, A2, B, C, D; `#124`, `#132`)
Piece G3 underwent two deep adversarial audit rounds (`#96` and `#106`), catching **9 genuine soundness bugs** across [`analyst/certify_G3_ABC_full_iv.py`](analyst/certify_G3_ABC_full_iv.py) and [`analyst/certify_G3_D_full_iv.py`](analyst/certify_G3_D_full_iv.py), ALL of which were patched cleanly and verified to **0 violations** in [`verify/audit_bk_g3_post106.py`](verify/audit_bk_g3_post106.py):

- **Round 1 (`#96` fixes — [`verify/audit_bk_g3_abc_checks.py`](verify/audit_bk_g3_abc_checks.py)):**
  1. `V_bump_A1_sep`: used `a_in` instead of `a_sep_A1 = min(a_in, a_out)` when the two near boundary poles cross into the $y < 0$ outer branch on $[\theta_{\text{hit}}, \pi]$. Fixed.
  2. `make_geom` / `box_geom`: replaced `sinc(b)^2` by exact MVT lower bound `sinc(2*b)` in `c_rein`, `c_re5`, `c_re50`, and replaced `log(1.5)/0.5` by exact concavity slope `(log(1.5) - log(1+T_in))/(0.5 - T_in)` in `c_re50`. Fixed.
  3. Case G3-C blind $u_{\text{floor}}$ choice (counterexample $c = 10^7, \tau = 1.04 \implies g_{1,\text{mis}} = 2.18 \gg 0.219$): fixed via $u_0 = \max(u_{0,\text{fin}}^*, u_{\text{floor}})$, explicit 50-dps proof of $g_{1,\text{fin}}(0.42) > 0$, and exact $D = -1$ residue term $+\frac{\Delta_0}{1-\Delta_0}$ in both $g_{1,\text{mis}}$ and $c_{\text{drop}}$. Fixed.
  4. Uncertified constants `1.02`, `0.35`, `T_c5` in `ep_C2_out` / `ep2_out`: explicitly asserted in 50-dps `iv`. Fixed.
- **Round 2 (`#106` fixes — [`verify/audit_bk_g3_d_checks.py`](verify/audit_bk_g3_d_checks.py), [`verify/audit_bk_g3_d_checks2.py`](verify/audit_bk_g3_d_checks2.py), [`verify/audit_bk_g3_post106.py`](verify/audit_bk_g3_post106.py)):**
  1. **Missing $1/(1-s_j)$ in $W_{2,\text{fin}}$ & $W_{3,\text{fin}}$ (both ABC & D):** At $j = c + 1$, $\frac{\theta_c}{1-s_j} \approx 2\,\text{thc\_hi}$, violating `rel_W2` by up to 45% and underestimating `a_p`, `b_p` in `CW3_free_*` and `Pth_2t` in `cert_B1_B2_pert*` by 2x. Fixed cleanly via `thc_1ms = kap_ub / (2*pi*S_half*cos(al_b))` and exact `a_p`, `b_p` bounds.
  2. **Missing `mob_in`/`mob_out` in $|P(u)|^2$ lower bound (both ABC Case C & D `exact_R3_W3_D1`):** Omitting `mob_in`/`mob_out` underestimated $|P(u)|^{-6}$ by up to 16.3% ($D > 0$). Fixed cleanly with explicit `mob_in`/`mob_out`.
  3. **Unsound `d_m <= 0.05` & `x_f_lo <= 1e-6` fallbacks (D):** Triggered on `D2-B` (dividing by `0.05` when `d_m = 0.0412`) and on every valley bin touching $u_0$ (hardcoded cap `14.0` on unbounded $-\log(1 - r_f(u))$). Fixed cleanly via `assert d_m.a > 1e-4` / `d_m0.a > 1e-4` (direct division by `d_m`) and exact closed-form integral bound `int_mlog_ub` using $F(x) = (1+x)\log(1+x) - x\log x$.
  4. **Straddle bin omissions (D):** When `u0_lo < u0_hi`, `eval_Psi_diff_W1_iv` omitted the skipped parabola integral $\int_{u_{0,\text{lo}}}^{u_{0,\text{hi}}} 2(-N_{E,\text{par}} W_{2,\text{par}})\,du$ and used `(1 + rel_W1)` instead of `(2 + rel_W1)`. Fixed cleanly across `eval_Psi_diff_W1_iv`, `_eval_g1_fin_single_iv`, and `_eval_peak_curv_single_iv`.
  5. **Uncertified $T(u) > 0.05$ far tail & `E_out_0` / defensive cleanups (both ABC & D):** Added explicit `far_T5_helpers*` / `W1_W2_fin_tail_ub` on $\{T(u) > 0.05\}$ using $n_{\text{lo}}((1+T_{\text{in}})/1.05)^{n_{\text{lo}}/2} < 10^{-6}$, replaced $d_{\text{lo}}$ by exact `d_lo_u0` in `E_out_0` (1-line proof via $\operatorname{Re}\zeta(x_p+iy_p) < \operatorname{Re}\zeta(u_0+iy_p)$), added `te_hi *` to `dP_re_rel`, fixed `f_max` / `coeff_re2.a < 0`, asserted `D_lo >= 0`, `dr_lo >= 0`, `u0 <= 0.52`, `I1_lo > rem_hi`, `sqrt(1.65)*u_max*thc_hi <= 0.080`, and verified the lossless W1 u-grid cache and seamless 37-slice CLI partition `D1 + D2 + D3a1 + D3a2 + D3b1 + D3b2 = all`.

- **Independent post-`#106` execution of Cases G3-A1, A2, B, C ([`analyst/certify_G3_ABC_full_iv.py`](analyst/certify_G3_ABC_full_iv.py), exit code 0, 29.78s):**
  - **Case G3-A1** ($|D| \ge 2, s \le 2.0$, split at $|\delta_{\text{fin}}| = 0.08$): **`rat_G3 >= +0.08207`** (`A1a: +0.08207`, `A1b: +0.11879`).
  - **Case G3-A2** ($|D| \ge 2, s \ge 2.0$, including closed-form $s \ge 10.0$ tail): **`rat_G3 >= +0.40597`**.
  - **Case G3-B** ($|D| = 1, \tau \ge 40.0$): **`rat_G3 >= +0.13156`**.
  - **Case G3-C** ($|D| = 1, \tau \in (0, 1.05]$): **`rat_G3 >= +0.07706`**.
  - **Pigeonhole-Annulus perturbation:** `E_cubic_G3_sub <= 2.80e-6` (`L1_shift_fac <= 1.000010`), `NB_err <= 7.32e-18`, `Y_big_err <= 5.26e-11`.
  - **Empirical stress check ([`analyst/verify_G3_empirical_corners.py`](analyst/verify_G3_empirical_corners.py), exit code 0, 1.8s):** **32/32** extreme Piece G3 corners ($D = \pm 1, \pm 2, \pm 3$, $r_0 \in \{1, 2, 3, 5\}$, up to $n = 9.85\times 10^7$) PASS (`min empirical rat = +0.54579`).

- **Independent post-`#106` execution of Piece G3 Case D ([`analyst/certify_G3_D_full_iv.py`](analyst/certify_G3_D_full_iv.py), SHA-256 `4479818E725B262DDAB836D4656C1DC266F71EA6F7873E6D96536B725F6DDD19`):**
  - Full line-by-line post-`#106` code audit: **PASSED (0 violations across all 2,130 lines and 10 extreme test points in [`verify/audit_bk_g3_post106.py`](verify/audit_bk_g3_post106.py))**.
  - Independent execution of **all 6/6 seamless CLI parts** (= all 37/37 $\tau$-slices covering $[1.05, 40.0]$ identically to `--part all`, run with `--workers 3`, ALL **EXIT CODE 0**, total wall time `3054.59s` / `50.9 min`):
    - `--part D1` ($\tau \in [1.05, 1.50]$, 3 slices): **`D1 >= +0.1646`** (`58.96s`, exit 0)
    - `--part D2` ($\tau \in [1.50, 2.80]$, 7 slices): **`D2 >= +0.0545`** (`960.43s`, exit 0)
    - `--part D3a1` ($\tau \in [2.80, 4.00]$, 7 slices): **`D3 >= +0.0784`** (`570.95s`, exit 0)
    - `--part D3a2` ($\tau \in [4.00, 5.68]$, 7 slices): **`D3 >= +0.0840`** (`608.44s`, exit 0)
    - `--part D3b1` ($\tau \in [5.68, 11.2]$, 7 slices): **`D3 >= +0.0952`** (`552.30s`, exit 0)
    - `--part D3b2` ($\tau \in [11.2, 40.0]$, 6 slices): **`D3 >= +0.2055`** (`303.51s`, exit 0)
  - **Combined Piece G3 global verified lower bound (full exact integral over all of Box G3):** `I_full / L^1 >= +0.0545 > 0` (`A1 >= +0.08207`, `A2 >= +0.40597`, `B >= +0.13156`, `C >= +0.07706`, `D1 >= +0.1646`, `D2 >= +0.0545`, `D3 >= +0.0784`; `certify_G3_ABC_full_iv.py` SHA-256 `7727BD251CE9AE60510FF2D7C24F546A1093CF299FC71E2C25F958E0BFDB3642`).
- **SHA-256 Reproducibility Checksums of Final Regime 2 Master Certifiers:**
  - `algebra/certify_G2_D0.py`: `B241F6C3C5ADEF5F560F49DEEB75944AF9EC00B46862D1359A65015B8D40C037`
  - `algebra/certify_G1_minor.py`: `F3EDAABAD4E69D78BD8CC7DD673E872D266A2E4920CE15C5A6F532E3C6D42186`
  - `analyst/certify_G3_ABC_full_iv.py`: `7727BD251CE9AE60510FF2D7C24F546A1093CF299FC71E2C25F958E0BFDB3642`
  - `analyst/certify_G3_D_full_iv.py`: `4479818E725B262DDAB836D4656C1DC266F71EA6F7873E6D96536B725F6DDD19`

---

## 4. Complete Inventory of Independent `@verifier` Audit Scripts

All independent referee scripts reside in `bardakov-kozlovskaya/verify/`:
- [`verify/engine.py`](verify/engine.py) & [`verify/stress_to_84.py`](verify/stress_to_84.py) — Exact integer engine ($n = 2..84$, 9,464 coefficients; SymPy/Wendt cross-checks).
- [`verify/fast_certify_140.py`](verify/fast_certify_140.py) — Fast $S_3$-Kronecker exact integer certifier ($n = 86..130$, 23,953 coefficients).
- [`verify/audit_bk_40.py`](verify/audit_bk_40.py), [`verify/audit_bk_43.py`](verify/audit_bk_43.py), [`verify/audit_bk_50.py`](verify/audit_bk_50.py) — Independent checks of Resultants & all 4 Boundary Families ($c=0, a=0, b=0, c=1$).
- [`verify/audit_bk_regime1.py`](verify/audit_bk_regime1.py) — Independent 50-dps `iv` verification of Theorem 4.2 & closed-form Regime 1 ($\kappa \ge 2.5$).
- [`verify/audit_bk_r2_1_2A_2C.py`](verify/audit_bk_r2_1_2A_2C.py) — Independent checks of Lemmas R2.1, R2.2 (Sub-regime 2A), and R2.3 (Sub-regime 2C-Upgrade).
- [`verify/audit_bk_g2_checks.py`](verify/audit_bk_g2_checks.py) & [`verify/check_g2_4_pole.py`](verify/check_g2_4_pole.py) — Independent audit of Piece G2 ($D = 0$ Major Arcs, including `#85` bug catches & `#103` widened $q=1$ seam).
- [`verify/audit_bk_g3_g1_checks.py`](verify/audit_bk_g3_g1_checks.py) — Independent audit of Piece G1 (Dirichlet Minor Arcs, including `#98` 5-bug catch & `#101` re-verification).
- [`verify/audit_bk_g3_abc_checks.py`](verify/audit_bk_g3_abc_checks.py) — Round 1 (`#96`) adversarial counterexample checks on Piece G3 Cases A1/A2/B/C.
- [`verify/audit_bk_g3_d_checks.py`](verify/audit_bk_g3_d_checks.py) & [`verify/audit_bk_g3_d_checks2.py`](verify/audit_bk_g3_d_checks2.py) — Round 2 (`#106`) adversarial counterexample checks uncovering the 5 shared/G3-D bugs.
- [`verify/audit_bk_g3_post106.py`](verify/audit_bk_g3_post106.py) — Post-`#106` verification confirming **0 violations** across all fixed Piece G3 formulas.

---

## 5. Final Verdict

**100% VERIFIED & SIGNED OFF — ZERO REMAINING MATHEMATICAL, ALGORITHMIC, OR COMPUTATIONAL DEFECTS.**

1. **All 10 of 10 proof components** covering all even integers $n = 2m \ge 2$ and all admissible $(a, b, c) \in \mathbb{Z}_{\ge 0}^3$ with $a + 2b + 3c = 2m$ have been independently audited line-by-line and executed by `@verifier` with **EXIT CODE 0** and **0 violations**:
   - Finite exact integer base ($n = 2..130$, $65$ even degrees, $33{,}417$ exact integer coefficients): $0$ zeros, $0$ sign violations (`#39`, `#57`).
   - Exact resultants & all 4 boundary families ($c = 0, a = 0, b = 0, c = 1$, all even $n \ge 2$): exact closed-form proofs (`#42`, `#45`, `#58`).
   - 1D circle integral & closed-form Regime 1 ($\kappa \ge 2.5, n \ge 84$): **`+0.469913`** ($\kappa \ge 3.0$) / **`+0.074799`** ($\kappa \ge 2.5$), 13,006 finite-$n$ cross-checks (`#59`, `#61`).
   - Formal Zero-Crevice Partition of Regime 2 ($\kappa < 2.5, n \ge 132$, `#103`): algebraic completeness + 3,012/3,012 finite-$n$ coverage check (0 uncovered points).
   - Lemma R2.1 & Sub-regime 2A ($c\kappa^2 \ge 42$): **`>= +0.220613`** (`#74`).
   - Sub-regime 2C-Upgrade ($cq\kappa \le 0.90$): **`>= +0.441480`** ($q = 1$) / **`>= +0.266331`** ($q \ge 2$) (`#77`).
   - Piece G2 ($D = 0$ Major Arcs): **`>= +0.216651`** ($q = 1, \kappa < 2.5$) / **`>= +0.00176348`** ($q \ge 2, q^2\kappa \le 2.0$), 28/28 parity corners PASS (`#90`, `#103`).
   - Piece G1 (Dirichlet Minor Arcs): **`>= +0.03355852`** ($\max|\Phi^{\text{ex}}| \le 1.536448 < \pi/2$), 27/27 finite-$n$ corners PASS (`#101`).
   - Piece G3 ($D \ne 0$ Major Arcs): **`>= +0.0545`** globally (**A1 `>= +0.08207`**, **A2 `>= +0.40597`**, **B `>= +0.13156`**, **C `>= +0.07706`**, **D1 `>= +0.1646`**, **D2 `>= +0.0545`**, **D3 `>= +0.0784`** across all 6/6 seamless CLI parts = 37/37 $\tau$-slices), 32/32 empirical corners PASS (`#124`, `#132`).
2. **All 17 adversarial soundness bugs** discovered across four audit rounds (`#85`, `#96`, `#98`, `#106`) were fixed rigorously and independently re-verified to **0 remaining violations**.
3. **Finalized master proof manuscript [`PROOF.md`](PROOF.md)** (475 lines, 31/31 relative proof/script links verified present, 0 remaining pending/draft/placeholder markers) has been verified line-by-line against every independently certified 50-digit `mpmath.iv` margin, exact SHA-256 checksum, runtime, and zero-crevice seam.

**Theorem 1.1 (The Jadav Sign-Law & Real-Rootedness Theorem) and Corollary 1.2 (Complete Resolution of Bardakov–Kozlovskaya Conjecture 3.3(2)) are unconditionally proved and 100% independently verified.**
