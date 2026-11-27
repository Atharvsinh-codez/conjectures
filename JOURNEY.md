# How We Solved This: Atharsinh Jadav & Gemini 4 Neon Extra High

> **Authors:** **[Atharsinh Jadav](https://athrix.me/)** ([X `@athrix_codes`](https://x.com/athrix_codes) · [GitHub `Atharvsinh-codez`](https://github.com/Atharvsinh-codez)), orchestrating the **Gemini 4 Neon Extra High multi-agent research swarm** (`coordinator`, `scout`, `algebraist`, `analyst`, `verifier`)  
> **Timeline:** 26–27 November 2026  
> **Outcome:** Complete, unconditionally proved, and 100% independently verified resolution of **Bardakov–Kozlovskaya Conjecture 3.3(2)** (*The Jadav Sign-Law & Real-Rootedness Theorem*), alongside sharp boundary obstructions and general real-stability theorems on **Gessel–Lin–Zeng Conjecture 15**.

This document chronicles how the campaign actually unfolded — including the live external collision we dodged on Day 1, the disciplined scientific pivot that turned a blocked attack into reusable general theorems, the dead ends we disproved along the way, and the **17 genuine soundness bugs** caught and killed by our adversarial `@verifier` agent before final sign-off.

---

## 1. The Multi-Agent Architecture

Rather than asking a single agent to write a proof and grade its own homework, the research campaign ran as a tightly coupled 5-role swarm on a shared message board:

- **`@coordinator`:** Strategic direction, go/pivot gates, scope discipline, and publication standards.
- **`@scout`:** Live literature, arXiv, SymCat, MathDB, and GitHub prior-resolution audits (ensuring we never waste compute on an already-solved conjecture or overlook an existing result).
- **`@algebraist`:** Structural identities, resultants, differential operators, real-stability theory, closed-form bounds, and Dirichlet minor-arc / $D=0$ major-arc certifiers.
- **`@analyst`:** Complex-analytic contour representations, saddle-point & shifted-line Poisson analysis, IBP kernels, and high-precision outward-rounded interval certifiers.
- **`@verifier`:** Independent, zero-trust adversarial referee — writing from-scratch exact-arithmetic checkers, hunting for edge-case counterexamples to our own lemmas, and holding strict veto power over sign-off.

---

## 2. Reconnaissance: Auditing 308 Conjectures & Dodging a Live Collision

The campaign began by indexing the entire unsolved landscape around [`1337isnot1337/solved-conjectures`](https://github.com/1337isnot1337/solved-conjectures) into [`UNLISTED_CONJECTURES.md`](UNLISTED_CONJECTURES.md) — a deduplicated catalogue of **308** still-unsolved conjectures spanning partial results on `main`, `docs/*-research.md`, `origin/research`, open GitHub issues, and classic open problems.

Our very first target was **Alexandersson's Catalan-word peak conjecture** ([Issue #59](https://github.com/1337isnot1337/solved-conjectures/issues/59)). Before writing a single solver script, `@scout` ran a live prior-resolution audit across GitHub, Zenodo, and SymCat — and caught a **same-day external resolution** published hours earlier on 26 November 2026: `ProgramComputer/open-math-research` [PR #62](https://github.com/ProgramComputer/open-math-research/pull/62) (Lean 4 formalization archived at [DOI:10.5281/zenodo.24552076](https://doi.org/10.5281/zenodo.24552076), already promoted to a Theorem on SymCat). Stopping immediately saved hours of duplicated effort.

`@scout` then screened **17 candidate conjectures across 6 groups** ([`docs/TARGET_SCREEN.md`](docs/TARGET_SCREEN.md)), eliminating 4 already-solved internal issues (`#495`, `#497`, `#498`, `#499`) and 3 active external-collision targets (`#81`, `#108`, `#83`), and elevating two 100%-open, collision-free targets with deep structural leverage:
1. **Campaign 1:** Gessel–Lin–Zeng Conjecture 15 ([Issue #84](https://github.com/1337isnot1337/solved-conjectures/issues/84), [arXiv:1201.0622](https://arxiv.org/abs/1201.0622))
2. **Campaign 2:** Bardakov–Kozlovskaya Conjecture 3.3(2) ([Issue #174](https://github.com/1337isnot1337/solved-conjectures/issues/174), [arXiv:2305.07261](https://arxiv.org/abs/2305.07261), [MathDB #358383](https://mathdb.com/p/358383/conjecture-on-coefficients-of-multivalued-group-polynomials))

---

## 3. Campaign 1 — Gessel–Lin–Zeng Conjecture 15: Deep Insights, Sharp Barriers & a Disciplined Pivot

We attacked Gessel–Lin–Zeng Conjecture 15 (real-rootedness of Jacobi–Stirling descent polynomials $A_{k,i}(x)$ for $4 \le i \le k - 3$) by asking a structured cascade of mathematical questions and testing every single one in exact rational arithmetic:

1. **Can we decouple the $i$-dependence of the recurrence?**  
   *Answer: Yes.* Möbius transformation $x = t/(1+t)$ yielded the $i$-independent recurrence $\tilde{h}_{k,i} = D_{k+1}(D_{k+1}\tilde{h}_{k-1,i} + \tilde{h}_{k-1,i-1})$ ($D_j = t(1+t)\partial_t + 1 + jt$), and $s$-space transformation ($x = -e^s$) collapsed the entire bivariate generating function into a pure $k$-th power of one Sturm–Liouville operator: $\Psi_k(s,y) = \mathcal{T}_y^k(e^s/(1+e^s))$ ([`gessel-lin-zeng/ANALYTIC_ROUTE.md`](gessel-lin-zeng/ANALYTIC_ROUTE.md)).
2. **Does the surviving binomial-weighted generating function $B_k(x, y) = \sum \binom{k}{i} A_{k,i}(x) y^i$ stay real-stable?**  
   *Answer: No — disproved at $k = 10$!* Exact `Fraction` Rayleigh determinant at $(0, -3/4)$ is $-2649516238079970755/4294967296 < 0$, backed by a degree-19 squarefree Sturm line with only $17$ real roots ([`gessel-lin-zeng/verify/verify_msg13.py`](gessel-lin-zeng/verify/verify_msg13.py)).
3. **Can an equal-degree $D$-lift + Plücker Quad + Cubic Real-Stability Addition close the induction?**  
   *Answer: Partially — yielding three standalone general theorems!* Setting $E_{k,i} = D_{k+1}^i \tilde{h}_{k-1,i}$ proved that the Pascal + outer $D_{k+1}^2$ step is **zero-loss** via:
   - **Jadav's Cubic Real-Stability Addition Theorem** (15-line proof in an open half-plane Wronskian cone),
   - **Sharpest Quantitative Plücker Quad Theorem B** ($J \ge \frac{4(1+\varepsilon)^2}{(2+\varepsilon)^2}$ with zero upper bound on $R$),
   - **Universal Cubic Wronskian Log-Concavity** ($W_{01}W_{23} \le W_{12}^2$ global SOS), and
   - **Zero-Loss $D_m$ Preservation** ($m \ge 1$).  
   Every single difficulty condensed into one single step: the index shift **Stage 3** ($D_{k+1}^i \mapsto D_{k+2}^i$).
4. **Can any universal operator theorem, linear gauge, or extension of Ma–Wang (Oct 2026) bypass Stage 3?**  
   *Answer: Provably NO!* We proved **four exact obstructions**:
   - **Jadav's Ma–Wang Sharpness Theorem:** Ma–Wang's Theorem 2 ([arXiv:2610.03111](https://arxiv.org/abs/2610.03111)) is **100% sharp** across every remaining open boundary ($i = 4$ fails first at $k = 31, x = -23$ with a 280-digit negative Laguerre integer; $i = 5$ at $k = 16, x = -11$; $i = k-3$ at $k = 28, x = -73/10$ and $k = 29, x = -8$).
   - **Universal operator products fail even for all-integer $b_j \ge 2$:** Exact finite single-op counterexample ($K = 1000$, Cauchy index $40 \ne 42$) and double-op SymPy scaling limit ($w \ge 5$).
   - **Five-line linear-gauge impossibility** and **failure of black-box 1-step Stage 3** even at GLZ's actual boundary scale $\lambda = 3/j$.
   - **Failure of global-$\mathbb{R}$ Turán positivity** ($A_{19,8}^2 - A_{19,7}A_{19,9}$ has a negative $x^1$-coefficient; $A_{35,17}^2 - A_{35,16}A_{35,18} < 0$ at $x = 7/10^{17}$).

Rather than overclaiming or grinding heuristically past a genuine barrier, the swarm executed a disciplined **Go/Pivot gate**: we independently verified and archived all GLZ theorems and disproofs (`243` Python scripts, `36,875` lines), stated honestly that Conjecture 15 remains open, and pivoted full force to **Bardakov–Kozlovskaya Conjecture 3.3(2)**.

---

## 4. Campaign 2 — Cracking Bardakov–Kozlovskaya Conjecture 3.3(2)

### Phase 1: From Nonvanishing to the Universal Sign Law & Exact Boundaries
Bardakov and Kozlovskaya's 2023 conjecture only asked whether every admissible coefficient $c_{a,b,c}^{(2m)}$ of the even-degree $n$-valued Buchstaber–Novikov group polynomial $p_{2m}(x, y, z) = \sum c_{a,b,c}^{(2m)} e_1^a e_2^b e_3^c$ is **nonzero**. Prior research notes observed through $n \le 40$ that every coefficient actually seemed to obey a rigid sign law: $\operatorname{sign}(c_{a,b,c}^{(2m)}) = (-1)^{b(c+1)+c}$.

- `@verifier` built a fast exact-integer Newton + peeling engine ([`verify/engine.py`](bardakov-kozlovskaya/verify/engine.py)), cross-checked it identically against SymPy resultants and Buchstaber–Kornev's Wendt circulant determinant, and pushed verification to $n = 84$ (`9,464` coefficients), then built an $S_3$-symmetric balanced-base Kronecker engine ([`verify/fast_certify_140.py`](bardakov-kozlovskaya/verify/fast_certify_140.py)) certifying all even $n \le 130$ (**33,417 exact integer coefficients**, **0 zeros**, **0 sign violations**), plus all `3,640` raw monomial partitions through $n = 60$.
- `@algebraist` discovered two exact structural resultants — the **Depressed-Cubic Resultant** (Theorem 3.1) and the **Half-Size Lucas–Chebyshev Resultant** $P_{2m}(1+q-p, q, p) = \operatorname{Res}_u(L_{2m}(u) + u^{2m} - q, u^{2m} - p)$ (Theorem 3.2) — and used them to prove **all four boundary families for every even $n \ge 2$**:
  - $c = 0$: $Q_0^{(2m)}(e_1, e_2) = (e_1^2 - 4e_2)^m$ (Theorem 3.3).
  - $a = 0$ and $b = 0$: exact trigonometric derivatives of $L_{2m}(1/(4\cos^2\theta)) = 2\cos(2m\theta)/(2\cos\theta)^{2m}$ isolating every root on open half-periods (Theorem 3.4).
  - $c = 1$: Rolle real-rootedness of $F_r(u)$, palindromic factorization of $A_{2m-1}^{(2m)}(u)/(1+u)$, Maclaurin's inequality, and an exact alternating trigonometric sum proving $\gamma_{m-1} > 4^{m-1}$ (Theorem 3.5).

### Phase 2: An Early Disproof That Forced the True Breakthrough
Early probes showed that every interior slice polynomial $Q_c^{(2m)}(1, u)$ was simple real-rooted through $n \le 50$. An obvious temptation was to prove real-rootedness by degree-induction interlacing $G_{n-2, c} \ll G_{n, c}$.

`@verifier` immediately tested that conjecture in exact `Fraction` arithmetic — and **disproved it in both parities** (`#42`, [`verify/audit_bk_40.py`](bardakov-kozlovskaya/verify/audit_bk_40.py)): explicit counterexamples occurred at $(n, c) = (16, 4)$, $(26, 7)$, and $(13, 3)$!

That disproof was pure gold. Because degree-to-degree interlacing was false, `@algebraist` and `@analyst` went straight to the **direct complex-analytic geometry of the roots**, uncovering the centerpiece of the entire proof:
1. **Universal Trigonometric Sample Points (Theorem 4.1):** $s_j = \sin(j\pi/n)/\sin((j+c)\pi/n)$, $t_j = (-1)^c s_j^n$, $u_j = t_j/(1+t_j)^2$ ($j \in \mathcal{J}_{n,c}$), whose strict sign alternation $(-1)^j Q_c^{(2m)}(1, u_j) > 0$ simultaneously proves **simple real-rootedness** AND the **exact sign law** via Vieta's formulas!
2. **The Exact 4-Factor 1D Cauchy Circle Integral (Theorem 4.2):** Integrating $p_n/(xyz)^c$ around $e_1 = 1+t, e_2 = t, e_3 = t R_z w$ cancelled every oscillating phase miraculously into the **strictly nonnegative** boundary kernel $B^{\text{pure}}(e^{i\phi}) = 4\sin^2((\phi + j\pi)/2) \ge 0$ for **both** parities of $c$!

### Phase 3: Closed-Form Regime 1 & the Hardy–Littlewood Decomposition of Regime 2
With separation parameter $\kappa = 2\pi \frac{\sin\psi_j\sin\alpha_1}{\sin\theta_c}$, the problem split naturally into two worlds:
- **Regime 1 ($\kappa \ge 2.5, n \ge 84$):** `@algebraist` proved the entire Regime 1 in **pure closed-form elementary inequalities across 6 subcases** (zero computer search), achieving worst margins **`+0.469913`** ($\kappa \ge 3.0$) and **`+0.074799`** ($\kappa \ge 2.5$), cross-checked across $13{,}006$ finite-$n$ points.
- **Regime 2 ($\kappa < 2.5, n \ge 132$):** When $\kappa$ becomes small, roots crowd tightly near the boundary and rational resonances $j/c \approx p_0/r_0$ appear — an unmistakable **Hardy–Littlewood circle-method** landscape! The swarm designed a **5-piece Zero-Crevice Partition**:
  - **Sub-regime 2A ($c\kappa^2 \ge 42$, `@analyst`):** Mirror pairing of small & large roots keeps $|\arg F_3^{\text{pure}}| \le 1.332945 < \pi/2$ on $|w| = 1$ (`>= +0.220613`).
  - **Sub-regime 2C-Upgrade ($cq\kappa \le 0.90$, `@analyst`):** Rescaled zero-exterior-root circle + Pólya–Szegő positive-root ratio + KL negative-binomial harmonics (`>= +0.441480` / `+0.266331`).
  - **Piece G2 ($D = 0$ Major Arcs, `@algebraist`):** All-parity contour identity + 2nd-order Euler–Maclaurin + $\pi$-line meromorphic residues controlled by the 2-line **Cauchy–Riemann level-curve monotonicity law** $\frac{d}{dy} V(x_*(y), y) = (U_x^2 + U_y^2)/U_x < 0$ + Pigeonhole-Annulus perturbation (`>= +0.216651` for $q=1$; `>= +0.00176348` for $q \ge 2$).
  - **Piece G1 (Dirichlet Minor Arcs, `@algebraist`):** Zero-cancellation **Hubbard–Stratonovich duality** + **continued-fraction Fibonacci ladder** + 1-line **$D_r^2$-cancellation identity** + exact $m$-space comparison certifying $\sup|\Phi^{\text{ex}}| \le 1.536448 < \pi/2$ (`>= +0.03355852`).
  - **Piece G3 ($D \ne 0$ Major Arcs, `@analyst`):** Exact **$(1-w)$-free IBP kernels** $W_{1,\text{fin}}, W_{2,\text{fin}}, W_{3,\text{fin}}$, **Harmonicity identity** $\frac{d}{du_0} g_{1,\text{fin}} = (-\operatorname{Re}\Psi_{\text{fin}}''(0))(-\operatorname{Re}\zeta'(u_0))$, zero-additive-constant shifted-line Poisson + 4-pole residue formula, horizontal-saddle shift, and 37-slice 50-digit `mpmath.iv` certification (`>= +0.0545` globally across A1, A2, B, C, D1, D2, D3).

---

## 5. Adversarial Self-Auditing: 17 Soundness Bugs + 1 Zero-Crevice Seam Caught & Killed

What separates a *plausible numerical argument* from an *airtight mathematical theorem* is relentless adversarial auditing. `@verifier` inspected every single line of all four Regime 2 master certifiers (~7,600 lines of 50-digit interval arithmetic) and built dedicated counterexample scripts (`bardakov-kozlovskaya/verify/audit_bk_*.py`) to try to break every inequality.

Across **four major audit rounds**, `@verifier` caught **17 genuine mathematical/soundness bugs** — plus `@algebraist` caught and closed a subtle **zero-crevice seam** — and **every single one** was fixed cleanly and re-verified to **0 violations**:

| Audit Round | Target Component | Bugs Caught | Concrete Examples of Defects Killed | Verification Script |
|---|---|---:|---|---|
| **Round 1 (`#85`)** | Piece G2 (`certify_G2_D0.py`) | **3** | Non-monotone peak `A_pk = max(A0, pi*q_min)` in `ksqA`/`ksq_tau`, wrong $A$-powers in `rhs_w`, and `1/A0` vs `1/sqrt(A0*T_e)` in `rhs_w1` (`cert_Ep_over_sqrtA_iv`). | [`verify/audit_bk_g2_checks.py`](bardakov-kozlovskaya/verify/audit_bk_g2_checks.py) |
| **Round 2 (`#96`)** | Piece G3-ABC (`certify_G3_ABC_full_iv.py`) | **4** | `V_bump_A1_sep` using `a_in` instead of `min(a_in, a_out)` across $y<0$ crossing; `sinc(b)^2` vs exact MVT `sinc(2*b)` and concavity slope in `c_re50`; Case G3-C blind $u_{\text{floor}}$ failing at $c=10^7, \tau=1.04$ (`g1_mis = 2.18`) + missing $D=-1$ residue; uncertified outer-pole constants. | [`verify/audit_bk_g3_abc_checks.py`](bardakov-kozlovskaya/verify/audit_bk_g3_abc_checks.py) |
| **Round 3 (`#98`)** | Piece G1 (`certify_G1_minor.py`) | **5** | `D_rat` Fibonacci index edge cases for $d \in \{3,4\}$ and $d \ge 5$; `12/11` vs `1` in B2b for unbounded `Dt >= 12`; increasing $t$-ratio in `earl_split[12]`; `(Tb/Ta)^2` vs `Tb/Ta` in `A_tail`; `T_lg2_p` failing at $(132,2,3)$ fixed uniformly via exact $\kappa_{\text{trig}}$. | [`verify/audit_bk_g3_g1_checks.py`](bardakov-kozlovskaya/verify/audit_bk_g3_g1_checks.py) |
| **Seam Catch (`#102`–`#103`)** | Zero-Crevice Partition | **1 seam** | Caught uncovered band $q = 1, \kappa \in (2.0, 2.5)$ between Piece G2 ($q^2\kappa \le 2.0$) and Regime 1 ($\kappa \ge 2.5$); closed by widening $q = 1$ in `certify_G2_D0.py` to `K_G2_q1 = 2.5` (explicit discrete $c \in \{2..10\}$ + automatic $\kappa < \sqrt{42/11} < 1.9541$ for $c \ge 11$). | [`algebra/certify_G2_D0.py`](bardakov-kozlovskaya/algebra/certify_G2_D0.py) |
| **Round 4 (`#106`)** | Piece G3-ABC & G3-D (`certify_G3_D_full_iv.py`) | **5** | Missing $1/(1-s_j)$ factor in $W_{2,\text{fin}}/W_{3,\text{fin}}$ at $j=c+1$ (45% error); missing `mob_in`/`mob_out` in $|P(u)|^2$ (16.3% error); unsound `d_m <= 0.05` & hardcoded `14.0` valley-bin cap (replaced by exact $(1+x)\log(1+x)-x\log x$ bound `int_mlog_ub`); omitted straddle-bin parabola integral; uncertified $T(u) > 0.05$ far tail & `d_lo_u0`. | [`verify/audit_bk_g3_d_checks.py`](bardakov-kozlovskaya/verify/audit_bk_g3_d_checks.py), [`d_checks2.py`](bardakov-kozlovskaya/verify/audit_bk_g3_d_checks2.py), [`post106.py`](bardakov-kozlovskaya/verify/audit_bk_g3_post106.py) |

---

## 6. Parallel 50-Digit Certification & Final 10/10 Sign-Off

With all 17 bugs patched and verified to **0 violations** in [`verify/audit_bk_g3_post106.py`](bardakov-kozlovskaya/verify/audit_bk_g3_post106.py), the final hurdle was executing the complete 50-digit outward-rounded `mpmath.iv` certification of Piece G3 Case D (37 $\tau$-slices on $[1.05, 40.0]$, 420+ boxes). Adding lossless per-box `W1` $u$-grid caching and `--workers 3` multiprocessing cut wall time cleanly while guaranteeing identical bit-for-bit interval enclosures:

- `@verifier` independently executed all 6 seamless CLI parts (`D1`, `D2`, `D3a1`, `D3a2`, `D3b1`, `D3b2` = 37/37 slices, `3054.59s` total, all **EXIT CODE 0**), matching `@analyst`'s certified margins bit-for-bit (`D1 >= +0.1646`, `D2 >= +0.0545`, `D3 >= +0.0784`).
- SHA-256 checksums of all flagship certifiers were locked into [`bardakov-kozlovskaya/AUDIT.md`](bardakov-kozlovskaya/AUDIT.md) and preserved cross-platform via `.gitattributes`.
- On **27 November 2026**, `@verifier` posted final **10/10 sign-off** (`#134`) on [`bardakov-kozlovskaya/PROOF.md`](bardakov-kozlovskaya/PROOF.md) and [`bardakov-kozlovskaya/AUDIT.md`](bardakov-kozlovskaya/AUDIT.md): **Bardakov–Kozlovskaya Conjecture 3.3(2) is unconditionally proved.**

---

## 7. Summary Ledger: Questions Formulated, Theorems Proved & Repository Metrics

### Research Questions Formulated & Resolved

| # | Campaign | Core Question Formulated | Outcome |
|---:|---|---|---|
| 1 | Recon | Is Alexandersson Catalan peaks (`#59`) still genuinely open? | **Answered (No):** Caught live external Lean 4 resolution (`PR #62`, Zenodo `24552076`), saving the swarm from collision. |
| 2 | GLZ `#84` | Does binomial-weighted $B_k(x,y)$ stay real-stable? | **Disproved** at $k=10$ (exact Rayleigh `< 0` & Sturm deg 19 / 17 real roots). |
| 3 | GLZ `#84` | Does equal-degree $D$-lifting + cubic real-stability addition work? | **Proved** (*Jadav's Cubic Real-Stability Addition*, *Sharpest Plücker Quad B*, *Cubic Wronskian Log-Concavity*, *Zero-Loss $D_m$*), isolating Stage 3 as the sole barrier. |
| 4 | GLZ `#84` | Can Ma–Wang Thm 2 or universal operator products close $4 \le i \le k-3$? | **Disproved** (*Jadav's Ma–Wang Sharpness Theorem* across $i=4,5,k-3$; all-integer $b_j \ge 2$ single- & double-op counterexamples; linear-gauge impossibility). |
| 5 | BK `#174` | Does the universal sign law $(-1)^{b(c+1)+c}$ hold beyond $n=40$? | **Verified** across all **33,417** exact integer coefficients for $n=2..130$ (0 zeros, 0 violations) and **proved for ALL even $n \ge 2$**. |
| 6 | BK `#174` | Can all 4 boundaries ($c=0, a=0, b=0, c=1$) be closed algebraically for all even $n$? | **Proved** via Depressed-Cubic & Half-Size Lucas–Chebyshev resultants, exact trig derivatives, and Rolle/Maclaurin/alternating trig sum. |
| 7 | BK `#174` | Do interior slices satisfy cross-degree interlacing $G_{n-2,c} \ll G_{n,c}$? | **Disproved** at $(16,4), (26,7), (13,3)$ — pivoting directly to universal sample points $u_j$ and the 1D Cauchy circle integral. |
| 8 | BK `#174` | Can Regime 1 ($\kappa \ge 2.5$) be proved in closed form with zero search? | **Proved** (6 elementary subcases, margins `+0.469913` / `+0.074799`). |
| 9 | BK `#174` | Can all of Regime 2 ($\kappa < 2.5, n \ge 132$) be certified with zero crevices? | **Proved & 100% Verified** across 5 overlapping sub-regimes (2A, 2C, G2, G1, G3) in 50-digit outward-rounded `mpmath.iv`, after catching & fixing all 17 audit bugs + 1 seam. |

### Final Repository Metrics

| Metric | Bardakov–Kozlovskaya (`bardakov-kozlovskaya/`) | Gessel–Lin–Zeng (`gessel-lin-zeng/`) | Combined Repository |
|---|---:|---:|---:|
| **Conjecture Status** | **✅ SOLVED (10/10 signed off)** | **OPEN** *(sharp theorems & disproofs)* | 1 Solved · 1 Open (Archived) |
| **Python Scripts (`.py`)** | **182** (`verify/`: 15, `algebra/`: 69, `analyst/`: 98) | **243** (`verify/`: 8, `algebra/`: 104, `analyst/`: 131) | **425 `.py` scripts** |
| **Lines of Python Code** | **32,715** | **36,875** | **69,590 lines** |
| **Exact Integer Coefficients Certified** | **33,417** ($n = 2..130$, 0 errors) | $k \le 24$ exact Sturm/Cauchy + exact Laguerre certs | 100% exact & 50-dps `iv` |
| **Adversarial Audit Bugs Caught & Killed** | **17 bugs + 1 zero-crevice seam** (0 remaining) | **3 hypothesis defects fixed** (Plücker Quad, Zero-Loss $D_m$) | **Zero known defects** |

---

*See also: [`THEOREMS.md`](THEOREMS.md) · [`EQUATIONS.md`](EQUATIONS.md) · [`bardakov-kozlovskaya/PROOF.md`](bardakov-kozlovskaya/PROOF.md) · [`bardakov-kozlovskaya/AUDIT.md`](bardakov-kozlovskaya/AUDIT.md) · [`gessel-lin-zeng/AUDIT.md`](gessel-lin-zeng/AUDIT.md) · [Root `README.md`](README.md).*
