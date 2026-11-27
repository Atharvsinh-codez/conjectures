# Bardakov–Kozlovskaya Conjecture 3.3(2) — Complete Proof & Verification Suite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24601452.svg)](https://doi.org/10.5281/zenodo.24601452)

> **Status: ✅ SOLVED & 100% INDEPENDENTLY VERIFIED (2026-11-27)**  
> **Published Preprint (Zenodo DOI):** [*A universal sign law and slice real-rootedness for even-order multivalued group polynomials*](https://zenodo.org/records/24601452) (`doi:10.5281/zenodo.24601452`) · [LaTeX Source (`paper/main.tex`)](paper/main.tex)  
> **Authors:** **[Atharvsinh Jadav](https://athrix.me/)** ([X `@athrix_codes`](https://x.com/athrix_codes) · [GitHub `Atharvsinh-codez`](https://github.com/Atharvsinh-codez)), with the **Gemini 4 Neon Extra High multi-agent research swarm**  
> **Source Problem:** V. G. Bardakov & T. A. Kozlovskaya, *Multivalued groups and Newton polyhedron*, [arXiv:2305.07261](https://arxiv.org/abs/2305.07261) (May 2023), Conjecture 3.3(2); *Sib. Elektron. Mat. Izv.* **20**:2 (2023), 1590–1596  
> **Upstream Tracking:** [`1337isnot1337/solved-conjectures` Issue #174](https://github.com/1337isnot1337/solved-conjectures/issues/174) · [MathDB #358383](https://mathdb.com/p/358383/conjecture-on-coefficients-of-multivalued-group-polynomials)

This directory contains the complete, unconditionally proved, and 100% independently verified resolution of **Bardakov–Kozlovskaya Conjecture 3.3(2)** — proved via a strictly stronger structural result, **The Jadav Sign-Law & Real-Rootedness Theorem** — together with all **6 manuscripts** and all **182 Python verification and research scripts (32,715 lines)**.

---

## 1. Exact Conjecture & Main Theorem

For every positive integer $n$, the defining symmetric homogeneous degree-$n$ integer polynomial of the $n$-valued Buchstaber–Novikov group $\mathbb{G}_n$ on $\mathbb{C}$ is
$$p_n(x, y, z) := \prod_{k=0}^{n-1}\Bigl(z - (-1)^n(\alpha + \zeta^k \beta)^n\Bigr) = \sum_{\substack{a,b,c \ge 0 \\ a + 2b + 3c = n}} c_{a,b,c}^{(n)}\, e_1^a e_2^b e_3^c,$$
where $\alpha^n = x$, $\beta^n = y$, $\zeta = e^{2\pi i/n}$, and $e_1 = x+y+z$, $e_2 = xy+yz+zx$, $e_3 = xyz$.

- **Bardakov–Kozlovskaya Conjecture 3.3(2) (2023):** *For every even integer $n = 2m \ge 2$ and every admissible triple $(a, b, c) \in \mathbb{Z}_{\ge 0}^3$ with $a + 2b + 3c = 2m$, the coefficient $c_{a,b,c}^{(2m)}$ is nonzero.* (Part (1), prime-power divisibility, was already resolved via Buchstaber–Kornev 2025/2026 and resultant Frobenius reduction; resolving Part (2) completely closes Conjecture 3.3.)

For each $0 \le c \le \lfloor 2m/3 \rfloor$, define the univariate $c$-slice polynomial of degree $D_c := \lfloor (2m - 3c)/2 \rfloor$ by
$$Q_c^{(2m)}(1, u) := \sum_{b=0}^{D_c} c_{2m-3c-2b,\, b,\, c}^{(2m)}\, u^b.$$

### Theorem 1.1 (The Jadav Sign-Law & Real-Rootedness Theorem — [`PROOF.md`](PROOF.md))

For every even integer $n = 2m \ge 2$:

1. **Exact Universal Sign Law (ALL even $n = 2m \ge 2$):** Every admissible coefficient $c_{a,b,c}^{(2m)}$ ($a + 2b + 3c = 2m$) is nonzero and obeys
   $$\operatorname{sign}\!\bigl(c_{a,b,c}^{(2m)}\bigr) = (-1)^{b(c+1)+c} = (-1)^{b(a+1)+c} = \begin{cases} (-1)^b & \text{if } c \text{ is even,} \\ -1 & \text{if } c \text{ is odd.} \end{cases}$$
   Equivalently, splitting $P_{2m} = P_{2m}^{\mathrm{even}} + P_{2m}^{\mathrm{odd}}$ by parity of the $e_3$-exponent $c$, every admissible coefficient of $P_{2m}^{\mathrm{even}}(e_1, -e_2, e_3)$ and of $-P_{2m}^{\mathrm{odd}}(e_1, e_2, e_3)$ is strictly positive.
2. **Universal Slice Real-Rootedness (exact scope):**
   - **Interior slices ($n \ge 132$, $2 \le c \le \lfloor (2m-3)/3 \rfloor$):** $Q_c^{(2m)}(1, u)$ has $D_c$ **simple real roots** — all strictly positive in $(0, +\infty)$ (with at least $D_c - 1$ in $(0, 1/4)$) when $c$ is even, and all strictly negative in $(-\infty, 0)$ when $c$ is odd — strictly separated by the universal trigonometric sample points $u_j = t_j / (1 + t_j)^2$, $t_j = (-1)^c (\sin(j\pi/n)/\sin((j+c)\pi/n))^n$.
   - **All four boundary slices (all even $n = 2m \ge 2$):** $Q_0^{(2m)}(1, u) = (1 - 4u)^m$; the $a = 0$ slice $e_2^{-m} P_{2m}(0, e_2, e_3)$ has $\lfloor m/3 \rfloor$ simple strictly positive real roots in $e_3^2/e_2^3$; the $b = 0$ slice $e_1^{-2m} P_{2m}(e_1, 0, e_3)$ has $\lfloor 2m/3 \rfloor$ simple strictly positive real roots in $e_3/e_1^3$; and for $c = 1$, $c_{2m-3-2b, b, 1}^{(2m)} = 2m((-4)^{b+1}\binom{m-1}{b+1} - \gamma_{b+1}) < 0$ for all $0 \le b \le m - 2$, where $\gamma_r > 4^r\binom{m-1}{r}$ via Rolle real-rootedness of $A_{2m-1}^{(2m)}(u)$, Maclaurin's inequality, and an exact alternating trigonometric sum. *(Note: Theorem 1.1 proves $c_{a,b,1}^{(2m)} < 0$ for all even $n \ge 4$ without claiming real-rootedness of $Q_1^{(2m)}(1, u)$.)*

### Corollary 1.2 (Complete Resolution of Bardakov–Kozlovskaya Conjecture 3.3)
For every even $n = 2m \ge 2$ and every $(a, b, c) \in \mathbb{Z}_{\ge 0}^3$ with $a + 2b + 3c = 2m$, $c_{a,b,c}^{(2m)} \ne 0$.

---

## 2. Ten-Component Verified Proof Architecture & Margins

Every one of the **10 proof components** was independently audited line-by-line and executed by `@verifier` with exit code `0` and zero remaining defects (full audit ledger, including all **17 adversarial soundness bugs** caught and patched across 4 audit rounds, in [`AUDIT.md`](AUDIT.md)):

| # | Proof Component | Exact Domain | Method / Arithmetic | Worst Certified Margin | Master Script |
|---:|---|---|---|---|---|
| 1 | **Finite Exact Base** | Even $n = 2, 4, \dots, 130$ (65 degrees, **33,417 coeffs**) | Exact integer Newton + $S_3$-Kronecker | **0 zeros, 0 sign violations** | [`verify/engine.py`](verify/engine.py), [`verify/fast_certify_140.py`](verify/fast_certify_140.py) |
| 2 | **Resultants & 4 Boundaries** | All even $n \ge 2$; $c=0, a=0, b=0, c=1$ | Depressed-Cubic & Lucas–Chebyshev resultants, Rolle, Maclaurin | **Exact closed form** | [`PROOF.md` §3](PROOF.md), [`verify/audit_bk_40.py`](verify/audit_bk_40.py)–[`50.py`](verify/audit_bk_50.py) |
| 3 | **1D Circle Integral & Regime 1** | $n \ge 84$, $2 \le c \le \lfloor(n-3)/3\rfloor$, $\kappa \ge 2.5$ | Pure closed-form 6-subcase elementary bounds (0 search) | **`+0.469913`** ($\kappa \ge 3.0$) / **`+0.074799`** ($\kappa \ge 2.5$) | [`algebra/verify_regime1_complete.py`](algebra/verify_regime1_complete.py) |
| 4 | **Zero-Crevice Partition** | Regime 2 ($\kappa < 2.5, n \ge 132$) | Algebraic completeness + 3,012/3,012 finite-$n$ check | **0 uncovered points** | [`PROOF.md` §1.3](PROOF.md) |
| 5 | **Lemma R2.1 (Perturbation)** | All Regime 2 ($\kappa < 2.5, n \ge 132$) | 50-digit outward-rounded `mpmath.iv` | $R_z/R_y \le e^{-37.8481/\kappa}$, $c_{\text{se}} \ge 3.16794$ | [`analyst/certify_R2_1_and_2A_iv.py`](analyst/certify_R2_1_and_2A_iv.py) |
| 6 | **Sub-regime 2A** | Regime 2, $c\kappa^2 \ge 42$ | Mirror root pairing, 50-digit `mpmath.iv` | **`>= +0.220613`** (`B_2A <= 1.332945 < pi/2`) | [`analyst/certify_R2_1_and_2A_iv.py`](analyst/certify_R2_1_and_2A_iv.py) |
| 7 | **Sub-regime 2C-Upgrade** | Regime 2, $c q \kappa \le 0.90$ | Pólya–Szegő + KL negative-binomial harmonics, 50-digit `mpmath.iv` | **`>= +0.441480`** ($q=1$) / **`>= +0.266331`** ($q\ge 2$) | [`analyst/certify_R2_3_iv.py`](analyst/certify_R2_3_iv.py) |
| 8 | **Piece G2 ($D=0$ Major Arcs)** | $c\kappa^2<42, cq\kappa\ge 0.85, d\ge 2$, $\kappa<2.5$ ($q=1$) or $q^2\kappa\le 2.0$ ($q\ge 2$) | $\pi$-line meromorphic residue + Cauchy–Riemann + 2nd-order Euler–Maclaurin + Pigeonhole-Annulus (`iv`) | **`>= +0.216651`** ($q=1$) / **`>= +0.00176348`** ($q\ge 2$); 28/28 corners | [`algebra/certify_G2_D0.py`](algebra/certify_G2_D0.py) |
| 9 | **Piece G1 (Minor Arcs)** | $c\kappa^2<42, q\ge 2$, $X\ge 0.88$ ($d=1$) or $X\ge 2.0$ ($d\ge 2$), outside G3 | Hubbard–Stratonovich duality + CF Fibonacci ladder + $D_r^2$-cancellation + $m$-space (`iv`) | **`>= +0.03355852`** ($\max|\Phi^{\text{ex}}|\le 1.536448<\pi/2$); 27/27 corners | [`algebra/certify_G1_minor.py`](algebra/certify_G1_minor.py) |
| 10 | **Piece G3 ($D\ne 0$ Major Arcs)** | $c\kappa^2<42, r_0 j = cp_0+D$ ($D\ne 0$), $\tilde{c}\ge 10D^2, r_0\kappa\le 0.45, r_0^3 c\kappa^2\le 60, \tilde{c}\tau^2/D^4\ge 40$ | $(1-w)$-free IBP kernels + shifted-line Poisson + 4-pole residue + horizontal saddle + Pigeonhole-Annulus (`iv`) | **`>= +0.0545`** (A1 `+0.08207`, A2 `+0.40597`, B `+0.13156`, C `+0.07706`, D1 `+0.1646`, D2 `+0.0545`, D3 `+0.0784`); 32/32 corners | [`analyst/certify_G3_ABC_full_iv.py`](analyst/certify_G3_ABC_full_iv.py), [`analyst/certify_G3_D_full_iv.py`](analyst/certify_G3_D_full_iv.py) |

---

## 3. SHA-256 Reproducibility Checksums (Flagship Certifiers)

All certified `.py` scripts are stored byte-for-byte unchanged with `.gitattributes` enforcing LF line endings (`*.py text eol=lf`) across all operating systems. Every checksum below matches [`AUDIT.md`](AUDIT.md) identically:

| Certified Script | Role | Verified Runtime | SHA-256 Checksum |
|---|---|---:|---|
| [`algebra/certify_G2_D0.py`](algebra/certify_G2_D0.py) | Piece G2 ($D=0$ Major Arcs, 50-dps `iv`) | 200.6s | `B241F6C3C5ADEF5F560F49DEEB75944AF9EC00B46862D1359A65015B8D40C037` |
| [`algebra/certify_G1_minor.py`](algebra/certify_G1_minor.py) | Piece G1 (Dirichlet Minor Arcs, 50-dps `iv`) | 755.4s | `F3EDAABAD4E69D78BD8CC7DD673E872D266A2E4920CE15C5A6F532E3C6D42186` |
| [`analyst/certify_G3_ABC_full_iv.py`](analyst/certify_G3_ABC_full_iv.py) | Piece G3 Cases A1, A2, B, C (50-dps `iv`) | 29.8s | `7727BD251CE9AE60510FF2D7C24F546A1093CF299FC71E2C25F958E0BFDB3642` |
| [`analyst/certify_G3_D_full_iv.py`](analyst/certify_G3_D_full_iv.py) | Piece G3 Case D (D1/D2/D3, 37 $\tau$-slices, 50-dps `iv`) | 3054.6s (`--workers 3`) | `4479818E725B262DDAB836D4656C1DC266F71EA6F7873E6D96536B725F6DDD19` |
| [`analyst/certify_R2_1_and_2A_iv.py`](analyst/certify_R2_1_and_2A_iv.py) | Lemma R2.1 & Sub-regime 2A (50-dps `iv`) | 0.3s | `E370D06B1ED07E9ED9554056A2CB1E8E03B0DE0C3A161151B2EA3EE3C55BC4E3` |
| [`analyst/certify_R2_3_iv.py`](analyst/certify_R2_3_iv.py) | Sub-regime 2C-Upgrade (50-dps `iv`) | 13.1s | `646B958B5A04A70837F3511D75356F4BC03CF51BD3BE126181B102BE3B02BACB` |
| [`algebra/verify_regime1_complete.py`](algebra/verify_regime1_complete.py) | Closed-form Regime 1 + 13,006 checks | ~2s | `F2F3D6FAD49FFD2C0BF16A470E2AFE030C759C1CAC5902173BAAB7D1E7807785` |
| [`verify/engine.py`](verify/engine.py) | Exact integer base ($n=2..84$, SymPy/Wendt cross-checks) | ~2s | `B8FB9F46F68B786A1E79EA20F916868F2CB370CA0DAA3FD240A138684FE4B84D` |
| [`verify/fast_certify_140.py`](verify/fast_certify_140.py) | Exact integer $S_3$-Kronecker base ($n=86..130$) | ~120s | `5A136E4DE9D41B9A910FC80CF0909B8EBB4D9C28A80800A654A2413021439EA3` |

---

## 4. Manuscripts & Complete Python Suite Inventory

### Proof & Audit Manuscripts (LaTeX Paper + 6 `.md` files)
- **[`paper/main.tex`](paper/main.tex)** *(Zenodo preprint [`doi:10.5281/zenodo.24601452`](https://zenodo.org/records/24601452))* — Standalone Overleaf-ready LaTeX research paper with 5 TikZ figures, 4 tables, and Appendices A-C.
- **[`PROOF.md`](PROOF.md)** *(475 lines)* — Complete, self-contained master proof manuscript (Theorem 1.1 & Corollary 1.2).
- **[`AUDIT.md`](AUDIT.md)** — Independent `@verifier` sign-off report covering all 10 proof components, all 17 adversarial soundness bugs caught & fixed, and SHA-256 hashes.
- **[`ALGEBRA_ROUTE.md`](ALGEBRA_ROUTE.md)** — Detailed algebraic derivations of Resultants, all 4 Boundaries, Closed-Form Regime 1, and Piece G2 Lemmas G2.1–G2.7.
- **[`ANALYTIC_ROUTE.md`](ANALYTIC_ROUTE.md)** — Detailed analytic derivations of Lemma R2.1, Sub-regimes 2A & 2C-Upgrade, and Piece G3 ($D \ne 0$ Major Arcs) Sub-Lemmas G3.1–G3.4.
- **[`algebra/G1_NOTES.md`](algebra/G1_NOTES.md)** — Complete Lemma G1.1–G1.5 derivation for Piece G1 (Dirichlet Minor Arcs).
- **[`SOURCES_AND_STATE.md`](SOURCES_AND_STATE.md)** — Initial literature & prior-resolution audit, exact indexing bijection, and 5 classical representations of $p_n$.

### Python Verification & Research Suite (**182 `.py` files, 32,715 lines**)

| Directory | `.py` Files | Total Lines | Contents & Key Entry Points |
|---|---:|---:|---|
| [`verify/`](verify/) | **15** | **3,403** | **Independent `@verifier` exact engines & adversarial audit scripts:** [`engine.py`](verify/engine.py), [`stress_to_84.py`](verify/stress_to_84.py), [`fast_certify_140.py`](verify/fast_certify_140.py), [`audit_bk_40.py`](verify/audit_bk_40.py), [`audit_bk_43.py`](verify/audit_bk_43.py), [`audit_bk_50.py`](verify/audit_bk_50.py), [`audit_bk_regime1.py`](verify/audit_bk_regime1.py), [`audit_bk_r2_1_2A_2C.py`](verify/audit_bk_r2_1_2A_2C.py), [`audit_bk_g2_checks.py`](verify/audit_bk_g2_checks.py), [`check_g2_4_pole.py`](verify/check_g2_4_pole.py), [`audit_bk_g3_g1_checks.py`](verify/audit_bk_g3_g1_checks.py), [`audit_bk_g3_abc_checks.py`](verify/audit_bk_g3_abc_checks.py), [`audit_bk_g3_d_checks.py`](verify/audit_bk_g3_d_checks.py), [`audit_bk_g3_d_checks2.py`](verify/audit_bk_g3_d_checks2.py), [`audit_bk_g3_post106.py`](verify/audit_bk_g3_post106.py). |
| [`algebra/`](algebra/) | **69** | **14,496** | **Flagship certifiers & algebraic discovery scripts:** [`certify_G2_D0.py`](algebra/certify_G2_D0.py), [`certify_G1_minor.py`](algebra/certify_G1_minor.py), [`verify_regime1_complete.py`](algebra/verify_regime1_complete.py), [`verify_exact_identities_mp.py`](algebra/verify_exact_identities_mp.py), [`verify_c1_and_probe_c2.py`](algebra/verify_c1_and_probe_c2.py), [`probe_depressed_cubic.py`](algebra/probe_depressed_cubic.py), plus 63 exploration, calibration, Hubbard–Stratonovich, CF ladder, and intermediate verification scripts. |
| [`analyst/`](analyst/) | **98** | **14,816** | **Flagship certifiers & analytic discovery scripts:** [`certify_R2_1_and_2A_iv.py`](analyst/certify_R2_1_and_2A_iv.py), [`certify_R2_3_iv.py`](analyst/certify_R2_3_iv.py), [`certify_G3_ABC_full_iv.py`](analyst/certify_G3_ABC_full_iv.py), [`certify_G3_D_full_iv.py`](analyst/certify_G3_D_full_iv.py), [`verify_G3_empirical_corners.py`](analyst/verify_G3_empirical_corners.py), [`verify_regime2_analyst_iv.py`](analyst/verify_regime2_analyst_iv.py), [`map_exact_crevices.py`](analyst/map_exact_crevices.py), plus 91 saddle-point, contour-shift, Poisson-summation, IBP, and boundary-probe scripts. |
| **Total** | **182** | **32,715** | **100% pure Python (exact integers + 50-digit outward-rounded `mpmath.iv` + SymPy).** |

---

## 5. Exact Reproduction Instructions

Requirements: **Python 3.10+** with `mpmath`, `sympy`, and `numpy`.

```bash
pip install mpmath sympy numpy
```

### Quick Smoke Checks (< 30 seconds total, run from repository root)
```bash
# 1. Exact 3-way convention verification (resultant vs Wendt det vs Newton engine)
python bardakov-kozlovskaya/verify/engine.py

# 2. Closed-form Regime 1 (all 6 subcases + 13,006 finite-n checks)
python bardakov-kozlovskaya/algebra/verify_regime1_complete.py

# 3. Lemma R2.1 + Sub-regime 2A (50-digit mpmath.iv, ~0.3s)
python bardakov-kozlovskaya/analyst/certify_R2_1_and_2A_iv.py

# 4. Sub-regime 2C-Upgrade (50-digit mpmath.iv, ~13s)
python bardakov-kozlovskaya/analyst/certify_R2_3_iv.py

# 5. All 32 extreme finite-n Piece G3 stress corners (up to n = 9.85e7)
python bardakov-kozlovskaya/analyst/verify_G3_empirical_corners.py
```

### Full End-to-End 50-Digit `mpmath.iv` & Exact Integer Certification
```bash
# Finite exact integer base: n = 2..84 (9,464 coeffs) and n = 86..130 (23,953 coeffs)
python bardakov-kozlovskaya/verify/stress_to_84.py
python bardakov-kozlovskaya/verify/fast_certify_140.py

# Piece G2 (D = 0 Major Arcs, ~201s)
python bardakov-kozlovskaya/algebra/certify_G2_D0.py

# Piece G1 (Dirichlet Minor Arcs, ~755s)
python bardakov-kozlovskaya/algebra/certify_G1_minor.py

# Piece G3 Cases A1, A2, B, C (~30s)
python bardakov-kozlovskaya/analyst/certify_G3_ABC_full_iv.py

# Piece G3 Case D (all 37 tau-slices via 6 seamless parts or --part all, ~51 min on --workers 3)
python bardakov-kozlovskaya/analyst/certify_G3_D_full_iv.py --part D1 --workers 3
python bardakov-kozlovskaya/analyst/certify_G3_D_full_iv.py --part D2 --workers 3
python bardakov-kozlovskaya/analyst/certify_G3_D_full_iv.py --part D3a1 --workers 3
python bardakov-kozlovskaya/analyst/certify_G3_D_full_iv.py --part D3a2 --workers 3
python bardakov-kozlovskaya/analyst/certify_G3_D_full_iv.py --part D3b1 --workers 3
python bardakov-kozlovskaya/analyst/certify_G3_D_full_iv.py --part D3b2 --workers 3
```

### Verify Flagship SHA-256 Checksums
```bash
# Linux / macOS
sha256sum bardakov-kozlovskaya/algebra/certify_G2_D0.py \
          bardakov-kozlovskaya/algebra/certify_G1_minor.py \
          bardakov-kozlovskaya/analyst/certify_G3_ABC_full_iv.py \
          bardakov-kozlovskaya/analyst/certify_G3_D_full_iv.py
```
```powershell
# Windows PowerShell
Get-FileHash bardakov-kozlovskaya/algebra/certify_G2_D0.py,
             bardakov-kozlovskaya/algebra/certify_G1_minor.py,
             bardakov-kozlovskaya/analyst/certify_G3_ABC_full_iv.py,
             bardakov-kozlovskaya/analyst/certify_G3_D_full_iv.py -Algorithm SHA256
```

---

*See also the repository-wide catalogues: [`THEOREMS.md`](../THEOREMS.md) · [`EQUATIONS.md`](../EQUATIONS.md) · [`JOURNEY.md`](../JOURNEY.md) · [Root `README.md`](../README.md).*
