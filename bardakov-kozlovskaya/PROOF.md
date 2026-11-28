# The Jadav Sign-Law & Real-Rootedness Theorem (Resolution of Bardakov–Kozlovskaya Conjecture 3.3(2))

**Authors:** Atharvsinh Jadav, with the Gemini 4 Neon Extra High Multi-Agent Research Swarm  
**Target Problem:** V. G. Bardakov and T. A. Kozlovskaya, *Multivalued groups and Newton polyhedron*, [arXiv:2305.07261](https://arxiv.org/abs/2305.07261) (May 2023), Conjecture 3.3(2); *Sibirskie Elektronnye Matematicheskie Izvestiya* **20**:2 (2023), 1590–1596, DOI `10.33048/semi.2023.20.097`  
**Tracking:** `1337isnot1337/solved-conjectures` [Issue #174](https://github.com/1337isnot1337/solved-conjectures/issues/174) · [MathDB #358383](https://mathdb.com/p/358383/conjecture-on-coefficients-of-multivalued-group-polynomials)  
**Document Status:** **COMPLETE & 100% INDEPENDENTLY VERIFIED.** All 10/10 proof components — Finite Base ($n \le 130$, `#39`, `#57`), Exact Resultants & All Four Boundary Families ($c=0, a=0, b=0, c=1$, `#42`, `#45`), Closed-Form Regime 1 ($\kappa \ge 2.5$, `#59`, `#61`), the Formal Zero-Crevice Partition (`#103`), and all five Regime 2 components (Lemma R2.1 & Sub-regime 2A `#74`, Sub-regime 2C-Upgrade `#77`, Piece G2 `#90`/`#103`, Piece G1 `#101`, and Piece G3 `#124`/`#132`) — have been independently audited and executed with exit code `0` and zero remaining defects.

---

## Abstract

For each positive integer $n$, the defining polynomial of the $n$-valued Buchstaber–Novikov group $\mathbb{G}_n$ on $\mathbb{C}$ is the symmetric homogeneous degree-$n$ integer polynomial
$$p_n(x, y, z) := \prod_{k=0}^{n-1}\Bigl(z - (-1)^n(\alpha + \zeta^k \beta)^n\Bigr) = \sum_{\substack{a, b, c \ge 0 \\ a + 2b + 3c = n}} c_{a,b,c}^{(n)}\, e_1^a e_2^b e_3^c,$$
where $\alpha^n = x$, $\beta^n = y$, $\zeta = e^{2\pi i/n}$, and $e_1 = x+y+z$, $e_2 = xy+yz+zx$, $e_3 = xyz$ are the elementary symmetric polynomials. In 2023, Bardakov and Kozlovskaya conjectured (Conjecture 3.3(2)) that for every **even** integer $n = 2m \ge 2$, every admissible coefficient $c_{a,b,c}^{(2m)}$ is nonzero.

We prove a substantially stronger structural theorem: **the Jadav Sign-Law & Real-Rootedness Theorem**. For every even $n = 2m \ge 2$ and every admissible triple $(a, b, c)$ with $a + 2b + 3c = 2m$, the coefficient $c_{a,b,c}^{(2m)}$ is nonzero and obeys the exact universal sign law
$$\operatorname{sign}\!\bigl(c_{a,b,c}^{(2m)}\bigr) = (-1)^{b(c+1)+c} = \begin{cases} (-1)^b & \text{if } c \text{ is even,} \\ -1 & \text{if } c \text{ is odd.} \end{cases}$$
Behind this sign law lies a hidden analytic geometry: for every $n = 2m \ge 132$ and every interior slice $2 \le c \le \lfloor (2m-3)/3 \rfloor$, the univariate $c$-slice polynomial $Q_c^{(2m)}(1, u) := \sum_{b=0}^{\lfloor (2m-3c)/2 \rfloor} c_{2m-3c-2b, b, c}^{(2m)}\, u^b$ is **strictly simple real-rooted**, with all roots positive when $c$ is even and all roots negative when $c$ is odd, strictly separated by an explicit universal family of trigonometric sample points $u_j$.

The proof combines: (i) an exact integer certification of all $33{,}417$ coefficients across $n = 2, 4, \dots, 130$; (ii) two exact structural resultants (Depressed-Cubic and Half-Size Lucas–Chebyshev) yielding elementary all-$n$ proofs of all four boundary families ($c = 0, a = 0, b = 0, c = 1$); (iii) an exact 1D Cauchy circle-integral representation at the universal sample points $u_j$; (iv) a completely closed-form elementary proof of Regime 1 ($\kappa \ge 2.5$, all $n \ge 84$, zero computer search); and (v) a zero-crevice Hardy–Littlewood-style circle decomposition of Regime 2 ($\kappa < 2.5$, all $n \ge 132$) into five overlapping sub-regimes — large phase-coupling ($c\kappa^2 \ge 42$), ultra-small rational-ray regime ($c q \kappa \le 0.90$), $D = 0$ major arcs, $D \ne 0$ major arcs, and Dirichlet minor arcs — certified uniformly over all $(n, c, j)$ in 50-digit outward-rounded interval arithmetic (`mpmath.iv`) with zero loops over the infinite parameter space.

---

## 1. Main Theorem & Zero-Crevice Proof Architecture

### 1.1 Notation and Indexing Bijection

Throughout this paper, let $n = 2m \ge 2$ be an even integer, $\zeta = e^{2\pi i/n}$, and let $e_1 = x+y+z$, $e_2 = xy+yz+zx$, $e_3 = xyz$. Since $(-1)^{2m} = +1$,
$$p_{2m}(x, y, z) = \prod_{k=0}^{2m-1}\Bigl(z - (\alpha + \zeta^k \beta)^{2m}\Bigr) = \prod_{r=0}^{2m-1}\prod_{s=0}^{2m-1}(\gamma + \zeta^r \alpha + \zeta^s \beta),$$
where $\alpha^{2m} = x$, $\beta^{2m} = y$, $\gamma^{2m} = z$. By the fundamental theorem of symmetric polynomials, there exists a unique polynomial $P_{2m}(e_1, e_2, e_3) \in \mathbb{Z}[e_1, e_2, e_3]$ such that $p_{2m}(x, y, z) = P_{2m}(e_1, e_2, e_3)$, written in two bijectively equivalent indexings:
$$P_{2m}(e_1, e_2, e_3) = \sum_{\substack{k_1 \ge k_2 \ge k_3 \ge 0 \\ k_1 + k_2 + k_3 = 2m}} A_{k_1, k_2, k_3}\, e_1^{k_1-k_2} e_2^{k_2-k_3} e_3^{k_3} = \sum_{\substack{a, b, c \ge 0 \\ a + 2b + 3c = 2m}} c_{a,b,c}^{(2m)}\, e_1^a e_2^b e_3^c,$$
where $a = k_1 - k_2$, $b = k_2 - k_3$, $c = k_3$, $A_{k_1,k_2,k_3} = c_{a,b,c}^{(2m)}$, and $a \equiv c \pmod 2$.

For each $0 \le c \le \lfloor 2m/3 \rfloor$, set $D_c := \lfloor (2m - 3c)/2 \rfloor$ and define the homogeneous $c$-slice polynomial $Q_c^{(2m)}(e_1, e_2)$ and its univariate specialization $Q_c^{(2m)}(1, u)$ by
$$Q_c^{(2m)}(e_1, e_2) := \sum_{b=0}^{D_c} c_{2m-3c-2b,\, b,\, c}^{(2m)}\, e_1^{2m-3c-2b} e_2^b, \qquad Q_c^{(2m)}(1, u) := \sum_{b=0}^{D_c} c_{2m-3c-2b,\, b,\, c}^{(2m)}\, u^b,$$
so that
$$P_{2m}(e_1, e_2, e_3) = \sum_{c=0}^{\lfloor 2m/3 \rfloor} Q_c^{(2m)}(e_1, e_2)\, e_3^c = \sum_{c=0}^{\lfloor 2m/3 \rfloor} e_1^{2m-3c}\, Q_c^{(2m)}\!\Bigl(1, \frac{e_2}{e_1^2}\Bigr)\, e_3^c.$$

### 1.2 Statement of the Main Theorem

> **Theorem 1.1 (The Jadav Sign-Law & Real-Rootedness Theorem).**  
> Let $n = 2m \ge 2$ be any even integer.
>
> 1. **(Exact Universal Sign Law):** For every triple $(a, b, c) \in \mathbb{Z}_{\ge 0}^3$ with $a + 2b + 3c = 2m$, the coefficient $c_{a,b,c}^{(2m)}$ is nonzero and satisfies
>    $$\operatorname{sign}\!\bigl(c_{a,b,c}^{(2m)}\bigr) \;=\; (-1)^{b(c+1)+c} \;=\; (-1)^{b(a+1)+c} \;=\; (-1)^{k_1 k_2 + k_1 k_3 + k_2 k_3} \;=\; \begin{cases} (-1)^b & \text{if } c \text{ is even,} \\ -1 & \text{if } c \text{ is odd.} \end{cases}$$
>    Equivalently, decomposing $P_{2m} = P_{2m}^{\mathrm{even}} + P_{2m}^{\mathrm{odd}}$ according to the parity of the $e_3$-exponent $c$, every admissible coefficient of $P_{2m}^{\mathrm{even}}(e_1, -e_2, e_3)$ and of $-P_{2m}^{\mathrm{odd}}(e_1, e_2, e_3)$ is strictly positive.
>
> 2. **(Universal Slice Real-Rootedness):**
>    - *(Interior slices, $n \ge 132$):* For every even $n = 2m \ge 132$ and every $2 \le c \le \lfloor (2m-3)/3 \rfloor$, the degree-$D_c$ univariate polynomial $Q_c^{(2m)}(1, u)$ has $D_c$ **simple real roots**: all strictly positive in $(0, +\infty)$ (with at least $D_c - 1$ strictly inside $(0, 1/4)$) when $c$ is even, and all strictly negative in $(-\infty, 0)$ when $c$ is odd. These roots are strictly interleaved by the universal trigonometric sample points $u_j$ of Theorem 4.1.
>    - *(All four boundary slices, all even $n = 2m \ge 2$):* $Q_0^{(2m)}(1, u) = (1 - 4u)^m$; the $a = 0$ slice $e_2^{-m} P_{2m}(0, e_2, e_3)$ has $\lfloor m/3 \rfloor$ simple strictly positive real roots in $e_3^2/e_2^3$; the $b = 0$ slice $e_1^{-2m} P_{2m}(e_1, 0, e_3)$ has $\lfloor 2m/3 \rfloor$ simple strictly positive real roots in $e_3/e_1^3$; and for $c = 1$, $c_{2m-3-2b, b, 1}^{(2m)} = 2m\bigl((-4)^{b+1}\binom{m-1}{b+1} - \gamma_{b+1}\bigr) < 0$ for all $0 \le b \le m - 2$, where $\gamma_r > 4^r \binom{m-1}{r}$ are the elementary symmetric functions of the positive reciprocal-pair parameters of the simple negative-real-rooted palindromic polynomial $A_{2m-1}^{(2m)}(u)$.

> **Corollary 1.2 (Complete Resolution of Bardakov–Kozlovskaya Conjecture 3.3).**  
> For every even integer $n = 2m \ge 2$ and every $(a, b, c) \in \mathbb{Z}_{\ge 0}^3$ with $a + 2b + 3c = 2m$, we have $c_{a,b,c}^{(2m)} \ne 0$. Together with the already-resolved Part (1) (prime-power divisibility; Buchstaber–Kornev 2025/2026), this completely resolves Conjecture 3.3 of Bardakov and Kozlovskaya (2023).

### 1.3 Master Zero-Crevice Partition Table

For $n \ge 132$, every coefficient $c_{a,b,c}^{(2m)}$ belongs either to one of the four boundary families ($c = 0, a = 0, b = 0, c = 1$; note $n - 3c \in \{0, 1, 2\}$ contains only $a = 0$ and $b = 0$ coefficients) or to an interior slice $2 \le c \le \lfloor (n-3)/3 \rfloor$. By Theorem 4.1, sign-definiteness and simple real-rootedness of $Q_c^{(2m)}(1, u)$ follow from $(-1)^j Q_c^{(2m)}(1, u_j) > 0$ at every strictly interior sample index
$$j \in \mathcal{J}_{n, c} := \bigl\{\, c + 1,\; c + 2,\; \dots,\; m - \lceil c/2 \rceil - \mathbf{1}_{c\text{ even}} \,\bigr\}.$$
Define the separation parameter $\kappa(n, c, j) := 2\pi \frac{\sin\psi_j \sin\alpha_1}{\sin\theta_c}$, where $\theta_c = \frac{c\pi}{n}$, $\psi_j = \frac{j\pi}{n}$, $\alpha_1 = \frac{(j+c)\pi}{n}$, and set $d := \gcd(c, j) \ge 1$, $q := c/d \ge 1$. The infinite parameter space is partitioned without crevices as follows:

| Component | Exact Parameter Domain | Method / Arithmetic | Worst Certified Margin | Audit Status |
|---|---|---|---|---|
| **Finite Base** (Sec. 2) | All even $n = 2, 4, \dots, 130$; all admissible $(a,b,c)$ ($33{,}417$ coeffs) | Exact integer arithmetic (Newton + Kronecker + $S_3$ symmetry) | 0 zeros, 0 sign violations | **SIGNED OFF** (`#39`, `#57`) |
| **Boundaries** (Sec. 3) | All even $n \ge 2$; $c=0$, $a=0$, $b=0$, $c=1$ | Pure algebraic/analytic proofs (Resultants, Rolle, Maclaurin) | Exact closed form | **SIGNED OFF** (`#42`, `#45`) |
| **Regime 1** (Sec. 4) | $n \ge 84$, $2 \le c \le \lfloor\frac{n-3}{3}\rfloor$, $j \in \mathcal{J}_{n,c}$, $\kappa \ge 2.5$ | Pure closed-form 6-subcase elementary inequalities (0 search) | **`+0.4699`** ($\kappa \ge 3.0$), **`+0.0748`** ($\kappa \ge 2.5$) | **SIGNED OFF** (`#59`, `#61`) |
| **Lemma R2.1** (Sec. 5.1) | Regime 2 universal perturbation ($\kappa < 2.5$, $n \ge 132$) | 50-digit `mpmath.iv` (0 loops) | $R_z/R_y \le e^{-37.8481/\kappa}$, $c_{\text{se}} \ge 3.16794$ | **SIGNED OFF** (`#74`) |
| **Sub-regime 2A** (Sec. 5.2) | Regime 2, $c\kappa^2 \ge 42$ | Exact small/large root pairing, 50-digit `mpmath.iv` | **`+0.220613 > 0`** | **SIGNED OFF** (`#74`) |
| **Sub-regime 2C-Upgrade** (Sec. 5.3) | Regime 2, $c q \kappa \le 0.90$ | Rescaled circle + Pólya–Szegő + KL harmonics, 50-digit `mpmath.iv` | **`+0.441480`** ($q=1$), **`+0.266331`** ($q\ge 2$) | **SIGNED OFF** (`#77`) |
| **Piece G2 ($D=0$ Major)** (Sec. 5.4) | Regime 2, $c\kappa^2 < 42$, $c q \kappa \ge 0.85$, $d \ge 2$, $\kappa < 2.5$ ($q=1$) or $q^2\kappa \le 2.0$ ($q\ge 2$) | $\pi$-line meromorphic residue + Cauchy–Riemann + Euler–Maclaurin + Pigeonhole-Annulus, 50-digit `mpmath.iv` | **`+0.00176348`** ($q\ge 2$), **`+0.216651`** ($q=1$) | **SIGNED OFF** (`#90`, `#103`) |
| **Piece G1 (Minor Arcs)** (Sec. 5.5) | Regime 2, $c\kappa^2 < 42$, $q \ge 2$, $X \ge 0.88$ ($d=1$) or $X \ge 2.0$ ($d\ge 2$), outside G3 core ($C_\tau = 45$) | Exact HS duality + CF Fibonacci ladder + $D_r^2$-cancellation + $m$-space + unit-circle perturbation, 50-digit `mpmath.iv` | **`+0.03355852 > 0`** ($|\Phi^{\text{ex}}| \le 1.536448 < \pi/2$) | **SIGNED OFF** (`#101`) |
| **Piece G3 ($D\ne 0$ Major)** (Sec. 6) | Regime 2, $c\kappa^2 < 42$, $r_0 j = c p_0 + D$ ($D\ne 0$), $\tilde{c}\ge 10 D^2$, $r_0\kappa\le 0.45$, $r_0^3 c\kappa^2\le 60$, $\tilde{c}\tau^2/D^4 \ge 40$ | Exact finite-$n$ IBP + shifted-line Poisson + 4-pole residue + Pigeonhole-Annulus, 50-digit `mpmath.iv` | **`+0.0545 > 0`** (A1 `+0.08207`, A2 `+0.40597`, B `+0.13156`, C `+0.07706`, D1 `+0.1646`, D2 `+0.0545`, D3 `+0.0784`) | **SIGNED OFF** (`#124`, `#132`) |

**Formal Zero-Crevice Partition Verification (All Even $n \ge 2$, All Admissible $(a,b,c)$):**
1. *Degree & boundary bridge:* Finite Base covers all $n \le 130$; Boundary Theorems 3.3–3.5 cover $c = 0, a = 0, b = 0, c = 1$ (and all $c > \lfloor (n-3)/3 \rfloor$) for all $n \ge 2$.
2. *$\kappa$-seam & $c\kappa^2$-seam:* For $n \ge 132$ and $2 \le c \le \lfloor (n-3)/3 \rfloor$, $j \in \mathcal{J}_{n,c}$, Regime 1 ($\kappa \ge 2.5$) and Regime 2 ($\kappa < 2.5$) are exact complements; inside Regime 2, Sub-regime 2A ($c\kappa^2 \ge 42$) and $\{c\kappa^2 < 42\}$ are exact complements.
3. *$cq\kappa$-seam & complete trichotomy inside $\{c\kappa^2 < 42, cq\kappa > 0.90\}$:* Sub-regime 2C-Upgrade covers $cq\kappa \le 0.90$. For any remaining point ($cq\kappa > 0.90 > 0.88 > 0.85$, $c\kappa^2 < 42$, $\kappa < 2.5$), set $d = \gcd(c, j) \ge 1$, $q = c/d \ge 1$, $X = q^2\kappa$:
   - *Case (i) $q = 1$:* Then $d = c \ge 2$, $cq\kappa > 0.85$, and $\kappa < 2.5 = K_{\text{G2}, q=1} \implies$ covered by **Piece G2**!
   - *Case (ii) $q \ge 2, d \ge 2, X \le 2.0$:* Then $d \ge 2$, $q^2\kappa \le 2.0 = K_{\text{G2}}$, $cq\kappa > 0.85 \implies$ covered by **Piece G2**!
   - *Case (iii) $q \ge 2$, and either $d = 1$ (where $X = cq\kappa > 0.90 > 0.88$) or ($d \ge 2$ and $X > 2.0$):* Either there exists a coprime $(p_0, r_0)$ with $D = r_0 j - c p_0 \ne 0$ satisfying all four inner Piece G3 inequalities ($\tilde{c} \ge 10 D^2, r_0\kappa \le 0.45, r_0^3 c\kappa^2 \le 60, \tilde{c}\tau^2/D^4 \ge 45 > 40$) $\implies$ covered by **Piece G3** (with strict overlap $[40, 45]$); or no such $(p_0, r_0)$ exists $\implies$ covered by **Piece G1**!

*(Machine-verified across all $3{,}012$ Regime 2 points for $n \in \{132, 134, \dots, 300\}$: $0$ uncovered points.)*

---

## 2. Finite Exact Integer Base ($n = 2, 4, \dots, 130$)

> **Theorem 2.1 (Exact Integer Verification for $n \le 130$).**  
> For every even integer $n \in \{2, 4, \dots, 130\}$ ($65$ even degrees) and every admissible triple $(a, b, c) \in \mathbb{Z}_{\ge 0}^3$ with $a + 2b + 3c = n$, the exact integer coefficient $c_{a,b,c}^{(n)}$ is nonzero and satisfies $\operatorname{sign}(c_{a,b,c}^{(n)}) = (-1)^{b(c+1)+c}$.

*Proof and Verification Methodology.* Viewing $p_n(x, y, z) = \sum_{j=0}^{n} (-1)^j E_j(x, y) z^{n-j}$ as a monic degree-$n$ polynomial in $z$ with roots $r_k = (\alpha + \zeta^k\beta)^n$, the $r$-th Newton power sum of its roots is given identically by
$$S_r(x, y) := \sum_{k=0}^{n-1} r_k^r = n \sum_{a=0}^{r} \binom{nr}{na} x^{r-a} y^a \qquad (n \text{ even}).$$
From $S_1, \dots, S_n$, Newton's identities $r E_r = \sum_{i=1}^{r} (-1)^{i-1} E_{r-i} S_i$ compute the exact integer homogeneous polynomials $E_r(x, y)$ using exact arbitrary-precision integer arithmetic (`verify/engine.py`). Converting $p_n(x, y, z)$ from the $(s_1, s_2, z)$-basis ($s_1 = x+y, s_2 = xy$) to the $(e_1, e_2, e_3)$-basis uses exact peeling by $c = 0, 1, \dots, \lfloor n/3 \rfloor$: at $z = 0$ one reads off $Q_c^{(n)}(s_1, s_2)$ directly, subtracts $Q_c^{(n)}(s_1 + z, s_2 + s_1 z)$ using pure integer additions, and divides cleanly by $e_3 = s_2 z$ with zero remainder.

For $n \in \{86, 88, \dots, 130\}$, `verify/fast_certify_140.py` exploits $S_3$-symmetry (every partition $k_1 \ge k_2 \ge k_3$ with $k_1+k_2+k_3 = n$ satisfies $k_2 + k_3 \le \lfloor 2n/3 \rfloor$, so $E_0, \dots, E_{\lfloor 2n/3 \rfloor}$ completely determine $p_n$ with redundant symmetry self-checks on all overlapping monomials) together with exact divide-and-conquer Kronecker polynomial multiplication in balanced base-$2^B$ representation (asserting zero carry overflow).

**Independent Cross-Checks and Exact Counts (signed off in `#39`, `#57`):**
1. `verify/engine.py` cross-verified identically against direct SymPy resultants $\operatorname{Res}_u(u^n - x, (t+u)^n - y)\big|_{t^n = z}$ ($n \le 8$) and the Buchstaber–Kornev $n \times n$ $y$-circulant Wendt determinant ($n \le 6$).
2. `verify/fast_certify_140.py` cross-verified to exact bit-for-bit equality against `verify/engine.py` across all $n \in \{2, 4, \dots, 30\}$.
3. **Exact certified counts:**
   - $n = 2, 4, \dots, 84$ ($42$ even degrees, `verify/engine.py`): **$9{,}464$** exact integer coefficients; $0$ zeros; $0$ sign-law violations.
   - $n = 86, 88, \dots, 130$ ($23$ even degrees, `verify/fast_certify_140.py`): **$23{,}953$** exact integer coefficients; $0$ zeros; $0$ sign-law violations.
   - **Combined ($n = 2, 4, \dots, 130$, $65$ even degrees):** **$33{,}417$** exact integer coefficients; $0$ zeros; $0$ sign-law violations. $\square$

---

## 3. Exact Structural Resultants & Complete Proofs of All Four Boundary Families

### 3.1 Two Exact Structural Resultants

> **Theorem 3.1 (Depressed-Cubic Resultant).**  
> Let $t^3 - u t - v = 0$ have roots $t_1, t_2, t_3$ (so $t_1 + t_2 + t_3 = 0$), and set $\eta := v^2$, $x = t_1^{2m}$, $y = t_2^{2m}$, $z = t_3^{2m}$. Then $e_3(x, y, z) = \eta^m$, while $e_1(x, y, z) = \Pi_{2m}(u, \eta)$ and $e_2(x, y, z) = \Sigma_{2m}(u, \eta)$ are given explicitly by Waring's formulas:
> $$\Pi_{2m}(u, \eta) := \sum_{r=0}^{\lfloor m/3 \rfloor} \frac{2m}{m - r}\binom{m - r}{2r} u^{m - 3r}\eta^r, \qquad \Sigma_{2m}(u, \eta) := \sum_{j=0}^{\lfloor 2m/3 \rfloor} (-1)^j \frac{m}{m - j}\binom{2m - 2j}{j} u^{2m - 3j}\eta^j.$$
> Consequently:
> - When $3 \nmid m$: $P_{2m}(e_1, e_2, \eta^m) = \operatorname{Res}_u\bigl(\Pi_{2m}(u, \eta) - e_1,\; \Sigma_{2m}(u, \eta) - e_2\bigr)$ (exact normalization $+1$).
> - When $m = 3K$: $\Pi_{6K}$ and $\Sigma_{6K}$ are polynomials in $w = u^3$, and $P_{6K}(e_1, e_2, e_3) = \prod_{\omega^3 = 1} R_K(e_1, e_2, \omega e_3^{1/3})$, where $R_K(e_1, e_2, \eta) := \operatorname{Res}_w\bigl(\tilde{\Pi}_{6K}(w, \eta) - e_1,\; \tilde{\Sigma}_{6K}(w, \eta) - e_2\bigr)$.

*Proof.* The vanishing $p_{2m}(x, y, z) = 0$ is equivalent to $t_1 + t_2 + t_3 = 0$ for some $2m$-th roots $t_1, t_2, t_3$ of $x, y, z$ (since $-1 = \zeta^m$). Such $(t_1, t_2, t_3)$ are the roots of a depressed cubic $t^3 - u t - v = 0$ with $v^{2m} = xyz = e_3$ and $u = -(t_1 t_2 + t_2 t_3 + t_3 t_1)$. Waring's power-sum formulas for $t_1^{2m} + t_2^{2m} + t_3^{2m}$ and $(t_1 t_2)^{2m} + (t_2 t_3)^{2m} + (t_3 t_1)^{2m}$ give $\Pi_{2m}(u, v^2)$ and $\Sigma_{2m}(u, v^2)$. Eliminating $u$ gives an exact multiple of $P_{2m}(e_1, e_2, \eta^m)$; comparing total degrees and specializing to $\eta = 0$ ($\Pi_{2m}(u, 0) = 2u^m$, $\Sigma_{2m}(u, 0) = u^{2m}$, giving $\operatorname{Res}_u(2u^m - e_1, u^{2m} - e_2) = (e_1^2 - 4e_2)^m$ when $3 \nmid m$) fixes the normalization to $+1$. $\square$

> **Theorem 3.2 (Half-Size Lucas–Chebyshev Resultant).**  
> For any integer $n \ge 1$, define the Lucas–Chebyshev polynomial
> $$L_n(u) := \sum_{r=0}^{\lfloor n/2 \rfloor} (-1)^r \frac{n}{n - r}\binom{n - r}{r} u^r,$$
> which satisfies the algebraic identity $L_n(t(1 - t)) = t^n + (1 - t)^n$. Then for every even integer $n = 2m \ge 2$ and all $p, q \in \mathbb{C}$,
> $$P_{2m}(1 + q - p,\, q,\, p) \;=\; \operatorname{Res}_u\!\bigl(L_{2m}(u) + u^{2m} - q,\; u^{2m} - p\bigr) \;=\; \prod_{k=0}^{2m-1}\Bigl(L_{2m}(\zeta^k p^{1/2m}) + p - q\Bigr).$$

*Proof.* First verify $L_n(t(1 - t)) = t^n + (1 - t)^n$: setting $x_1 = t$, $x_2 = 1 - t$ gives $x_1 + x_2 = 1$ and $x_1 x_2 = t(1 - t) = u$, so Waring's formula yields $x_1^n + x_2^n = \sum_{r=0}^{\lfloor n/2 \rfloor} (-1)^r \frac{n}{n-r}\binom{n-r}{r} (x_1 + x_2)^{n-2r}(x_1 x_2)^r = L_n(u)$. Now fix $z = 1$ in $p_{2m}(x, y, z) = 0$ and write $x^{1/2m} = -t$, $y^{1/2m} = -(1 - t)$ with $u := t(1 - t)$. Then $x + y = t^{2m} + (1 - t)^{2m} = L_{2m}(u)$ and $xy = u^{2m}$, so
$$e_3 = xyz = u^{2m} =: p, \qquad e_2 = xy + x + y = L_{2m}(u) + u^{2m} =: q, \qquad e_1 = x + y + 1 = 1 + q - p.$$
Over the $2m$ roots $u_k = \zeta^k p^{1/2m}$ of $u^{2m} - p = 0$, the product $\prod_{k=0}^{2m-1}(L_{2m}(u_k) + p - q)$ equals $\operatorname{Res}_u(L_{2m}(u) + u^{2m} - q,\, u^{2m} - p)$ (since $u^{2m} - p$ is monic) and vanishes identically on the hypersurface $P_{2m}(1 + q - p, q, p) = 0$. Both sides have degree $2m$ in $q$, and at $p = 0$ ($L_{2m}(0) = 1$) the RHS equals $(1 - q)^{2m}$, matching $P_{2m}(1 + q, q, 0) = ((1 + q)^2 - 4q)^m = (1 - q)^{2m}$ (Theorem 3.3 below). Hence the two polynomials are identically equal. $\square$

### 3.2 Complete Proofs of All Four Boundary Families (All Even $n = 2m \ge 2$)

> **Theorem 3.3 (Boundary Family $c = 0$).**  
> For every even $n = 2m \ge 2$, $Q_0^{(2m)}(e_1, e_2) = (e_1^2 - 4e_2)^m$, hence $c_{2m-2b, b, 0}^{(2m)} = (-4)^b \binom{m}{b}$ has exact sign $(-1)^b$ for all $0 \le b \le m$.

*Proof.* Setting $z = 0$ ($e_3 = 0$, $e_1 = x + y$, $e_2 = xy$) gives $p_{2m}(x, y, 0) = \prod_{k=0}^{2m-1}\bigl(-(\alpha + \zeta^k \beta)^{2m}\bigr) = \bigl(\prod_{k=0}^{2m-1}(\alpha + \zeta^k \beta)\bigr)^{2m} = (\alpha^{2m} - \beta^{2m})^{2m} = (x - y)^{2m} = (e_1^2 - 4e_2)^m$. $\square$

> **Theorem 3.4 (Boundary Families $a = 0$ and $b = 0$).**  
> For every even integer $n = 2m \ge 2$:
> 1. The univariate polynomial $f_a^{(2m)}(w) := \sum_{q=0}^{\lfloor m/3 \rfloor} c_{0,\, m-3q,\, 2q}^{(2m)}\, w^q$ has $\lfloor m/3 \rfloor$ simple strictly positive real roots in $w = e_3^2/e_2^3 \in (0, +\infty)$, and consequently $\operatorname{sign}(c_{0, m-3q, 2q}^{(2m)}) = (-1)^{m-q}$ for all $0 \le q \le \lfloor m/3 \rfloor$.
> 2. The univariate polynomial $f_b^{(2m)}(w) := \sum_{c=0}^{\lfloor 2m/3 \rfloor} c_{2m-3c,\, 0,\, c}^{(2m)}\, w^c$ has $\lfloor 2m/3 \rfloor$ simple strictly positive real roots in $w = e_3/e_1^3 \in (0, +\infty)$, and consequently $\operatorname{sign}(c_{2m-3c, 0, c}^{(2m)}) = (-1)^c$ for all $0 \le c \le \lfloor 2m/3 \rfloor$.

*Proof.* Parameterize $u \in [1/4, +\infty)$ by $u = \frac{1}{4\cos^2\theta}$ for $\theta \in [0, \pi/2)$. Since $t = \frac{1 + i\tan\theta}{2} = \frac{e^{i\theta}}{2\cos\theta}$ satisfies $1 - t = \frac{e^{-i\theta}}{2\cos\theta}$ and $t(1 - t) = \frac{1}{4\cos^2\theta} = u$, the defining identity $L_{2m}(t(1 - t)) = t^{2m} + (1 - t)^{2m}$ gives
$$L_{2m}\!\Bigl(\frac{1}{4\cos^2\theta}\Bigr) = \frac{2\cos(2m\theta)}{(2\cos\theta)^{2m}}.$$
Differentiating directly yields the two exact trigonometric derivative identities:
$$\frac{d}{d\theta}\!\left[\frac{2\cos(2m\theta)}{(2\cos\theta)^{2m}}\right] = \frac{-8m\sin((2m - 1)\theta)}{(2\cos\theta)^{2m+1}}, \qquad \frac{d}{d\theta}\!\Bigl[2(2\cos\theta)^{2m}\cos(2m\theta) + 1\Bigr] = -8m(2\cos\theta)^{2m-1}\sin((2m + 1)\theta).$$

**Part 1 ($a = 0$, i.e. $e_1 = 0$):** Setting $q = p - 1$ in Theorem 3.2 gives $P_{2m}(0, p - 1, p) = \operatorname{Res}_u(L_{2m}(u) + 1,\, u^{2m} - p)$. On the other hand, $P_{2m}(0, p - 1, p) = (p - 1)^m f_a^{(2m)}\bigl(\frac{p^2}{(p - 1)^3}\bigr)$. Thus every root $u_* \in (1, +\infty)$ of $L_{2m}(u) + 1 = 0$ produces a root $p_* = u_*^{2m} \in (1, +\infty)$ of $P_{2m}(0, p - 1, p)$, hence a strictly positive real root $w_* = \frac{p_*^2}{(p_* - 1)^3} \in (0, +\infty)$ of $f_a^{(2m)}(w) = 0$. Since $\frac{d}{dp}\frac{p^2}{(p - 1)^3} = -\frac{p(p + 2)}{(p - 1)^4} < 0$ on $(1, +\infty)$, distinct roots $u_* \in (1, +\infty)$ map injectively to distinct roots $w_* \in (0, +\infty)$.

Define $h_a(\theta) := L_{2m}(\frac{1}{4\cos^2\theta}) + 1 = \frac{2\cos(2m\theta)}{(2\cos\theta)^{2m}} + 1$ on $\theta \in (\pi/3, \pi/2)$ (where $2\cos\theta \in (0, 1)$, so $u \in (1, +\infty)$). Set $r_{\min} := m - \lfloor m/3 \rfloor$ and evaluate $h_a$ at the critical points $\theta_r := \frac{r\pi}{2m - 1}$ for $r = r_{\min}, r_{\min}+1, \dots, m - 1$. Since $3 r_{\min} = 2m + (m - 3\lfloor m/3 \rfloor) \ge 2m > 2m - 1$, every $\theta_r$ lies strictly inside $(\pi/3, \pi/2)$, so $2\cos\theta_r \in (0, 1)$, hence $(2\cos\theta_r)^{1-2m} > 1$. Because $2m\theta_r = r\pi + \theta_r$, we have $\cos(2m\theta_r) = (-1)^r\cos\theta_r$, giving
$$h_a(\theta_r) = 1 + (-1)^r (2\cos\theta_r)^{1 - 2m}, \qquad \operatorname{sign}(h_a(\theta_r)) = (-1)^r \qquad (r = r_{\min}, \dots, m - 1).$$
Finally, as $\theta \to \frac{\pi}{2}^-$, $2\cos\theta \to 0^+$ and $\cos(2m\theta) \to \cos(m\pi) = (-1)^m$, so $h_a(\theta) \sim 2(-1)^m (2\cos\theta)^{-2m}$ has sign $(-1)^m$. Thus across the $\lfloor m/3 \rfloor + 1$ ordered points $\theta_{r_{\min}} < \dots < \theta_{m-1} < \frac{\pi}{2}^-$, $h_a(\theta)$ strictly alternates sign $\lfloor m/3 \rfloor$ times, yielding $\lfloor m/3 \rfloor$ distinct roots in $(\pi/3, \pi/2)$ and hence all $\lfloor m/3 \rfloor$ simple strictly positive real roots of $f_a^{(2m)}(w)$. Since $f_a^{(2m)}(0) = c_{0, m, 0}^{(2m)} = (-4)^m$ (Theorem 3.3, sign $(-1)^m$), Vieta's formulas force $\operatorname{sign}(c_{0, m-3q, 2q}^{(2m)}) = (-1)^m (-1)^q = (-1)^{m-q}$.

**Part 2 ($b = 0$, i.e. $e_2 = 0$):** Setting $q = 0$ in Theorem 3.2 gives $P_{2m}(1 - p, 0, p) = \operatorname{Res}_u(L_{2m}(u) + u^{2m},\; u^{2m} - p) = (1 - p)^{2m} f_b^{(2m)}\bigl(\frac{p}{(1 - p)^3}\bigr)$. Every root $u_* \in (1/4, 1)$ of $L_{2m}(u) + u^{2m} = 0$ yields $p_* = u_*^{2m} \in (4^{-2m}, 1)$ and a root $w_* = \frac{p_*}{(1 - p_*)^3} \in (0, +\infty)$ of $f_b^{(2m)}(w) = 0$; since $\frac{d}{dp}\frac{p}{(1 - p)^3} = \frac{1 + 2p}{(1 - p)^4} > 0$ on $(0, 1)$, distinct $u_* \in (1/4, 1)$ map injectively to distinct $w_* \in (0, +\infty)$.

Define $h_b(\theta) := (2\cos\theta)^{4m}(L_{2m}(\frac{1}{4\cos^2\theta}) + (\frac{1}{4\cos^2\theta})^{2m}) = 2(2\cos\theta)^{2m}\cos(2m\theta) + 1$ on $\theta \in [0, \pi/3)$ (where $2\cos\theta \in (1, 2]$, so $u \in [1/4, 1)$). At the $\lfloor 2m/3 \rfloor + 1$ critical points $\tilde{\theta}_r := \frac{r\pi}{2m + 1}$ for $r = 0, 1, \dots, \lfloor 2m/3 \rfloor$, we have $\tilde{\theta}_{\lfloor 2m/3 \rfloor} \le \frac{2m\pi}{3(2m + 1)} < \frac{\pi}{3}$, so $2\cos\tilde{\theta}_r > 1$ and $(2\cos\tilde{\theta}_r)^{2m+1} > 1$. Because $2m\tilde{\theta}_r = r\pi - \tilde{\theta}_r$, we have $\cos(2m\tilde{\theta}_r) = (-1)^r\cos\tilde{\theta}_r$, giving
$$h_b(\tilde{\theta}_r) = 1 + (-1)^r (2\cos\tilde{\theta}_r)^{2m + 1}, \qquad \operatorname{sign}(h_b(\tilde{\theta}_r)) = (-1)^r \qquad (r = 0, 1, \dots, \lfloor 2m/3 \rfloor).$$
Since $h_b(0) = 2^{2m+1} + 1 > 0$, this yields $\lfloor 2m/3 \rfloor$ distinct roots strictly inside $(0, \pi/3)$, hence all $\lfloor 2m/3 \rfloor$ simple strictly positive real roots of $f_b^{(2m)}(w)$. Since $f_b^{(2m)}(0) = c_{2m, 0, 0}^{(2m)} = +1 > 0$, Vieta's formulas give $\operatorname{sign}(c_{2m-3c, 0, c}^{(2m)}) = (-1)^c$. $\square$

> **Theorem 3.5 (Boundary Family $c = 1$).**  
> For every even integer $n = 2m \ge 4$ and every $0 \le b \le m - 2$, we have $c_{2m-3-2b,\, b,\, 1}^{(2m)} < 0$.

*Proof.* Set $x = 1$, $y = u$, and differentiate $p_n(1, u, z) = P_n(e_1, e_2, e_3) = (e_1^2 - 4e_2)^m + e_3 Q_1^{(2m)}(e_1, e_2) + O(e_3^2)$ with respect to $z$ at $z = 0$ (where $e_1 = 1+u$, $e_2 = u$, $e_3 = 0$, $\partial_z e_1 = 1$, $\partial_z e_2 = 1+u$, $\partial_z e_3 = u$). Comparing the two sides yields:
- **From $P_n(e_1, e_2, e_3)$:** $\bigl[\partial_z P_n\bigr]_{z=0} = -2m(1 + u)(1 - u)^{2m - 2} + u\, Q_1^{(2m)}(1 + u, u)$.
- **From $p_n(1, u, z) = \prod_{k=0}^{n-1}(z - (1 + \zeta^k u^{1/n})^n)$:**
  $$\Bigl[\frac{\partial p_n}{\partial z}\Bigr]_{z=0} = -(1 - u)^n \sum_{k=0}^{n-1} \frac{1}{(1 + \zeta^k u^{1/n})^n} = -2m\, A_{2m - 1}^{(2m)}(u),$$
  where for each $1 \le r \le n$ we define $S_r(u) := \frac{1}{n}\sum_{k=0}^{n-1} (1 + \zeta^k u^{1/n})^{-r}$, $F_r(u) := (1 - u)^r S_r(u)$, and $A_{2m - 1}^{(2m)}(u) := F_{2m}(u)$.

**Step 1 (Rolle Real-Rootedness of $F_r$):** Clearly $S_1(u) = (1 - u)^{-1}$, so $F_1(u) = 1$. Differentiating $S_{r-1}(u)$ term-by-term gives $S_r(u) = S_{r-1}(u) + \frac{n u}{r - 1} S'_{r-1}(u)$; multiplying by $(1 - u)^r$ produces the exact polynomial recurrence (for $2 \le r \le n$):
$$(r - 1)\, F_r(u) \;=\; n u(1 - u)\, F'_{r-1}(u) \;+\; (r - 1)\bigl(1 + (n - 1)u\bigr)\, F_{r-1}(u).$$
By induction on $r \in \{1, \dots, n\}$, $F_r(u)$ is a polynomial of exact degree $r - 1$ with constant term $F_r(0) = 1$, leading coefficient $\binom{n - 1}{r - 1} > 0$, and $r - 1$ simple strictly negative real roots: indeed, at the $r - 2$ simple negative roots $\xi_1 < \dots < \xi_{r-2} < 0$ of $F_{r-1}$, the factor $n\xi_j(1 - \xi_j)$ is strictly negative while $F'_{r-1}(\xi_j)$ alternates sign, giving strict sign alternation of $F_r(\xi_j)$; since $F_r(0) = 1 > 0$ and $F'_{r-1}(\xi_{r-2}) > 0 \implies F_r(\xi_{r-2}) < 0$, there is a root in $(\xi_{r-2}, 0)$, $r - 3$ roots between consecutive $\xi_j$, and one root in $(-\infty, \xi_1)$ (as the leading coefficient is positive).

**Step 2 (Palindromic Factorization & Exact Closed Form for $c_{2m-3-2b, b, 1}^{(2m)}$):** For $r = n = 2m$, the substitution $u \mapsto u^{-1}$ permutes the terms of $S_{2m}$ and shows that $A_{2m-1}^{(2m)}(u) = F_{2m}(u)$ is monic, palindromic ($A_{2m-1}^{(2m)}(u^{-1}) = u^{-(2m-1)} A_{2m-1}^{(2m)}(u)$), of odd degree $2m - 1$, with $2m - 1$ simple strictly negative real roots. Therefore $u = -1$ is a root, and the remaining $2m - 2$ roots pair into reciprocals $(r_j, r_j^{-1})$ ($r_j < -1$, $j = 1, \dots, m - 1$). Setting $c_j := -(r_j + r_j^{-1}) - 2 = (\sqrt{-r_j} - 1/\sqrt{-r_j})^2 > 0$ and $\gamma_r := e_r(c_1, \dots, c_{m-1}) > 0$ ($\gamma_0 = 1$), we obtain
$$\frac{A_{2m-1}^{(2m)}(u)}{1 + u} = \prod_{j=1}^{m-1}\bigl((1 + u)^2 + c_j u\bigr) = \sum_{r=0}^{m-1} \gamma_r\, (1 + u)^{2m - 2 - 2r}\, u^r.$$
Since $(1 - u)^{2m - 2} = ((1 + u)^2 - 4u)^{m - 1} = \sum_{r=0}^{m-1} (-4)^r \binom{m - 1}{r} (1 + u)^{2m - 2 - 2r} u^r$, the $r = 0$ terms cancel identically in $u\, Q_1^{(2m)}(1 + u, u) = 2m\bigl((1 + u)(1 - u)^{2m - 2} - A_{2m - 1}^{(2m)}(u)\bigr)$, yielding
$$Q_1^{(2m)}(1 + u, u) = 2m \sum_{b=0}^{m-2} \Bigl((-4)^{b+1}\binom{m - 1}{b + 1} - \gamma_{b+1}\Bigr)(1 + u)^{2m - 3 - 2b} u^b,$$
and therefore the exact closed-form coefficient identity
$$c_{2m - 3 - 2b,\, b,\, 1}^{(2m)} \;=\; 2m\Bigl((-4)^{b+1}\binom{m - 1}{b + 1} - \gamma_{b+1}\Bigr) \qquad (0 \le b \le m - 2).$$

**Step 3 (Maclaurin's Inequality & Alternating Trigonometric Sum):** Since $c_1, \dots, c_{m-1} > 0$, Maclaurin's inequality gives $\bigl(\gamma_r / \binom{m-1}{r}\bigr)^{1/r} \ge \gamma_{m-1}^{1/(m-1)}$ for every $1 \le r \le m - 1$. Thus $\gamma_r > 4^r \binom{m-1}{r}$ for all $1 \le r \le m - 1$ — and hence $c_{2m-3-2b, b, 1}^{(2m)} < 0$ for all $b$ — follows from the single inequality $\gamma_{m-1} > 4^{m-1}$. Evaluating $\frac{A_{2m-1}^{(2m)}(u)}{1 + u}$ at $u = -1$ via L'Hôpital's rule gives $\gamma_{m-1} = (-1)^{m-1} (A_{2m-1}^{(2m)})'(-1)$. Using $u^{1/2m}\big|_{u=-1} = e^{i\pi/(2m)}$ and $1 + e^{i(2k+1)\pi/(2m)} = 2\cos\frac{(2k+1)\pi}{4m}\, e^{i(2k+1)\pi/(4m)}$, differentiating $A_{2m-1}^{(2m)}(u)$ at $u = -1$ and pairing $k$ with $2m - 1 - k$ yields the exact alternating trigonometric sum (independently verified in `#45`):
$$\gamma_{m - 1} \;=\; \frac{1}{2m}\sum_{k=0}^{m-1} (-1)^k \frac{\cos\theta_k}{\sin^{2m + 1}\theta_k}, \qquad \theta_k := \frac{(2k + 1)\pi}{4m}.$$
Because $\theta \mapsto \frac{\cos\theta}{\sin^{2m+1}\theta}$ is strictly decreasing on $(0, \pi/2)$ and $\frac{\sin\theta_1}{\sin\theta_0} = 3 - 4\sin^2\frac{\pi}{4m} \ge 3 - 4\sin^2\frac{\pi}{8} = 1 + \sqrt{2} > 2$ for all $m \ge 2$, the alternating sum exceeds its first two terms:
$$\gamma_{m - 1} \;>\; \frac{\cos(\pi/8)}{2m}\Bigl(\frac{4m}{\pi}\Bigr)^{2m + 1}\Bigl(1 - \frac{1}{2^{2m + 1}}\Bigr) \;>\; 4^{m - 1} \qquad \text{for all } m \ge 2.$$
(Explicitly: at $m = 2$, $\gamma_1 = 28 > 4$; for $m \ge 3$, $\frac{4m}{\pi} > 3.81$, so $(\frac{4m}{\pi})^{2m} > 14.5^m \gg 4^{m-1}$.) This proves $c_{2m-3-2b, b, 1}^{(2m)} < 0$ for all $m \ge 2$ and all $0 \le b \le m - 2$. $\square$

---

## 4. Universal Sample Points, Exact 1D Circle Integral & Closed-Form Proof of Regime 1

### 4.1 Reduction to Strictly Interior Sample Points

> **Theorem 4.1 (Universal Sample Points & Root-Interleaving Reduction).**  
> Fix even $n = 2m \ge 6$ and $2 \le c \le \lfloor (n - 3)/3 \rfloor$ (noting that $n - 3c \in \{0, 1, 2\}$ contains only $a = 0$ and $b = 0$ boundary coefficients, already proved in Theorem 3.4). Set $D_c := \lfloor (n - 3c)/2 \rfloor$ and define the strictly interior sample index set
> $$\mathcal{J}_{n, c} := \bigl\{\, c + 1,\; c + 2,\; \dots,\; m - \lceil c/2 \rceil - \mathbf{1}_{c\text{ even}} \,\bigr\}.$$
> For each $j \in \mathcal{J}_{n, c}$, define
> $$\theta_c := \frac{c\pi}{n}, \quad \psi_j := \frac{j\pi}{n}, \quad \alpha_1 := \frac{(j + c)\pi}{n}, \quad r_1 := \sin\alpha_1, \quad r_2 := \sin\psi_j, \quad r_3 := \sin\theta_c,$$
> $$s_j := \frac{r_2}{r_1} \in (0, 1), \quad R_y := s_j^n, \quad R_z := \Bigl(\frac{r_3}{r_1}\Bigr)^n, \quad t_j := (-1)^c R_y, \quad u_j := \frac{t_j}{(1 + t_j)^2}.$$
> **If $(-1)^j Q_c^{(2m)}(1, u_j) > 0$ for every $j \in \mathcal{J}_{n, c}$, then $Q_c^{(2m)}(1, u)$ has $D_c$ simple real roots — all in $(0, +\infty)$ (with at least $D_c - 1$ in $(0, 1/4)$) when $c$ is even, and all in $(-\infty, 0)$ when $c$ is odd — and hence every coefficient $c_{n-3c-2b, b, c}^{(2m)}$ ($0 \le b \le D_c$) is nonzero with exact sign $(-1)^{b(c+1)+c}$.**

*Proof.* For $j \in \mathcal{J}_{n,c}$ we have $0 < j < \min(j + c, n - j - c) \le n/2$, so $s_j \in (0, 1)$, and $\frac{d}{dx}\frac{\sin x}{\sin(x + \theta_c)} = \frac{\sin\theta_c}{\sin^2(x + \theta_c)} > 0$ shows $s_j$ is strictly increasing in $j$.
- **Even $c = 2q$:** Here $D_c = m - 3q$ and $|\mathcal{J}_{n, c}| = (m - q - 1) - 2q = D_c - 1$. As $j$ increases, $t_j = s_j^n \in (0, 1)$ increases strictly, so $0 < u_{c+1} < u_{c+2} < \dots < u_{m-q-1} < \frac{1}{4}$. Combining the $b = 0$ boundary $Q_c^{(2m)}(1, 0) = c_{n-3c, 0, c}^{(2m)} > 0$ (sign $(-1)^c = +1$, Theorem 3.4), the $D_c - 1$ strictly alternating signs $(-1)^j Q_c^{(2m)}(1, u_j) > 0$, and the $a = 0$ leading-term limit $\lim_{u \to +\infty}\operatorname{sign}(Q_c^{(2m)}(1, u)) = \operatorname{sign}(c_{0, m-3q, 2q}^{(2m)}) = (-1)^{m-q}$ (opposite to $(-1)^{m-q-1}$, Theorem 3.4) produces $D_c$ sign changes on $(0, +\infty)$ — hence $D_c$ simple positive real roots (with $D_c - 1$ in $(0, 1/4)$).
- **Odd $c = 2q + 1$:** Here $D_c = m - 3q - 2$ and $|\mathcal{J}_{n, c}| = (m - q - 1) - (2q + 1) = D_c$. As $j$ increases, $t_j = -s_j^n \in (-1, 0)$ decreases strictly, and $\frac{d}{dt}\frac{t}{(1 + t)^2} = \frac{1 - t}{(1 + t)^3} > 0$ on $(-1, 0)$, so $-\infty < u_{m-q-1} < \dots < u_{c+1} < 0$. Combining these $D_c$ strictly alternating evaluations with the $b = 0$ boundary $Q_c^{(2m)}(1, 0) < 0$ (sign $(-1)^c = -1$, Theorem 3.4) produces $D_c$ sign changes strictly inside $(u_{m-q-1}, 0) \subset (-\infty, 0)$, hence $D_c$ simple negative real roots.

In both parities, Vieta's formulas and the sign of $Q_c^{(2m)}(1, 0)$ force $\operatorname{sign}(c_{n-3c-2b, b, c}^{(2m)}) = (-1)^{b(c+1)+c}$ for all $0 \le b \le D_c$. $\square$

### 4.2 Exact 1D Circle Integral & 4-Factor Decomposition

> **Theorem 4.2 (Exact 1D Circle Integral).**  
> Fix even $n = 2m \ge 6$, $2 \le c \le \lfloor (n - 3)/3 \rfloor$, and $j \in \mathcal{J}_{n, c}$. For $|w| = \lambda > 0$, consider the cubic equation
> $$T(T - 1)(T - t) \;=\; t R_z w \qquad \bigl(t := t_j = (-1)^c R_y\bigr).$$
> Suppose Rouché's theorem isolates unique simple roots $x(w) = 1 + \delta_x(w)$ and $y(w) = t(1 + \tilde{\delta}_y(w))$ holomorphic on a neighborhood of $|w| \le \lambda$ with $|\delta_x| < 1$ and $|\delta_y| < 1$, where
> $$1 + \delta_y(w) := \frac{y(w)}{t\, x(w)}, \qquad 1 + \delta_z(w) := \frac{z(w)}{x(w) R_z w} = \frac{1}{(1 + \delta_x(w))^3(1 + \delta_y(w))},$$
> and let $\varepsilon_s(w) := (1 + \delta_y(w))^{1/n} - 1$ (principal branch). Let the $n$ angles $\eta_k \in (-\pi, \pi]$ be $\frac{2k\pi}{n}$ (even $c$) or $\frac{(2k+1)\pi}{n}$ (odd $c$), partitioned into **small** ($|\eta_k| < \theta_c$, count $c - 1$), **boundary** ($|\eta_k| = \theta_c$, namely $\pm\theta_c$), and **large** ($|\eta_k| > \theta_c$, count $n - c - 1$). Set $w_k := (1 - s_j e^{i\eta_k})^n$, $a_k^{\mathrm{alg}} := \frac{-s_j e^{i\eta_k}}{1 - s_j e^{i\eta_k}}$, $\tilde{w}_k := w_k (1 + a_k^{\mathrm{alg}}\varepsilon_s)^n$, $1 + \varepsilon_\pm(w) := (1 + a^{\mathrm{alg}}(\pm\theta_c)\varepsilon_s(w))^n$, and $Z := (-1)^j w^{-1}$. Then
> $$(-1)^j (1 + t)^{n - 3c} Q_c^{(2m)}(1, u_j) \;=\; \frac{C_{m,c,j}}{2\pi}\int_{-\pi}^{\pi} \bigl(1 + \Delta_{\text{mult}}(\lambda e^{i\phi})\bigr)\, B(\lambda e^{i\phi})\, d\phi, \qquad C_{m,c,j} := \frac{R_z \prod_{\text{large}} |w_k|}{R_y^c} > 0,$$
> where $1 + \Delta_{\text{mult}}(w) := F_1(w)\, F_2(w)\, F_3(w)$ and
> $$F_1(w) := (1 + \delta_x)^{n - 3c}, \qquad F_2(w) := \frac{\prod_{\text{large}}(1 + a_k^{\mathrm{alg}} \varepsilon_s)^n}{(1 + \delta_y)^c}, \qquad F_3(w) := \prod_{\text{small}}\!\Bigl(1 - \frac{\tilde{w}_k}{z/x}\Bigr)\prod_{\text{large}}\!\Bigl(1 - \frac{z/x}{\tilde{w}_k}\Bigr),$$
> $$B(w) := \frac{(-Z^{-1})\bigl((1 + \delta_z) - Z(1 + \varepsilon_+)\bigr)\bigl((1 + \delta_z) - Z(1 + \varepsilon_-)\bigr)}{1 + \delta_z}$$
> $$= (2 - Z - Z^{-1}) - Z^{-1}\delta_z + (\varepsilon_+ + \varepsilon_-) - Z\Bigl((1 + \varepsilon_+)(1 + \varepsilon_-)(1 + \delta_x)^3(1 + \delta_y) - 1\Bigr).$$

*Proof.* On the circle $e_1 = 1 + t$, $e_2 = t$, $e_3 = t R_z w$ ($|w| = \lambda$), Cauchy's integral formula gives
$$(1 + t)^{n - 3c} Q_c^{(2m)}(1, u_j) = Q_c^{(2m)}(1 + t, t) = \frac{1}{2\pi}\int_{-\pi}^{\pi} \frac{p_n(x(\lambda e^{i\phi}), y(\lambda e^{i\phi}), z(\lambda e^{i\phi}))}{(t R_z \lambda e^{i\phi})^c}\, d\phi.$$
Since $y/x = (-1)^c s_j^n(1 + \varepsilon_s)^n$, choosing $\alpha = x^{1/n}$ and $\beta = \alpha s_j(1 + \varepsilon_s) e^{i\pi c/n}$ gives
$$\frac{p_n(x, y, z)}{(xyz)^c} = x^{n - 3c}\frac{\prod_{k=0}^{n-1}\bigl((z/x) - \tilde{w}_k\bigr)}{(y/x)^c (z/x)^c}.$$
Factor $z/x$ from the $c - 1$ small roots, $-\tilde{w}_k = (-w_k)(1 + a_k^{\mathrm{alg}}\varepsilon_s)^n$ from the $n - c - 1$ large roots, and pair the single remaining denominator $z/x = R_z w(1 + \delta_z)$ with the two boundary roots. Because $1 - s_j e^{\pm i\theta_c} = \frac{r_3}{r_1} e^{\mp i\psi_j}$ (Law of Sines on the triangle with angles $\theta_c, \psi_j, \pi - \alpha_1$), we have $w_{\pm\theta_c} = (-1)^j R_z$, $\tilde{w}_{\pm\theta_c} = (-1)^j R_z (1 + \varepsilon_\pm)$, and $\frac{((z/x) - \tilde{w}_+)((z/x) - \tilde{w}_-)}{z/x} = R_z (-1)^{j+1} B(w)$. Finally, $\frac{\prod_{\text{large}}(-w_k)}{t^c} = -\frac{\prod_{\text{large}}|w_k|}{R_y^c}$ for **both** parities of $c$: non-real large roots come in conjugate pairs contributing $+|w_k|^2 > 0$; for even $c$, the unique real large root $\eta = \pi$ has $-w_\pi = -(1+s_j)^n < 0$ while $t^c = +R_y^c > 0$; for odd $c$, there is no real large root while $t^c = -R_y^c < 0$. Multiplying $(-1)^{j+1}$ by $-1$ gives $(-1)^j$, and substituting $(1 + \delta_z)^{-1} = (1 + \delta_x)^3 (1 + \delta_y)$ yields the expanded form of $B(w)$. $\square$

*(Independently verified to $< 10^{-50}$ relative error in 80-digit `mpmath` against exact integer coefficients in [`algebra/verify_exact_identities_mp.py`](algebra/verify_exact_identities_mp.py).)*

### 4.3 Elementary Lemmas & Closed-Form Proof of Regime 1 ($\kappa \ge 2.5$, All Even $n \ge 84$)

Define the separation parameter $\kappa(n, c, j) := 2\pi \frac{r_1 r_2}{r_3} = 2\pi \frac{\sin\psi_j \sin\alpha_1}{\sin\theta_c}$, together with $a_1 := \min(j + c,\, n - j - c)$, $p := n - 2j - c$, $q_{\mathrm{gap}} := j - c$, $d_{\mathrm{deg}} := n - 3c$.

> **Lemma 4.1 (Jensen Midpoint Root-Separation Bounds).**  
> For every strictly interior $(n, c, j)$ ($2 \le c \le \lfloor (n-3)/3 \rfloor$, $j \in \mathcal{J}_{n, c}$):
> 1. $0 < c < j < a_1 \le n/2$, $p \ge 1 + \mathbf{1}_{c\text{ even}}$, $q_{\mathrm{gap}} \ge 1$, $a_1 - j = \min(c, p) \ge 1 + \mathbf{1}_{c\text{ even}}$, $a_1 - c = q_{\mathrm{gap}} + \min(c, p) \ge 2 + \mathbf{1}_{c\text{ even}}$, and $d_{\mathrm{deg}} = p + 2q_{\mathrm{gap}} \ge 3$.
> 2. By convexity of $\cot x$ on $(0, \pi/2)$, for any $0 < u < v \le \pi/2$,
>    $$\Bigl(\frac{\sin u}{\sin v}\Bigr)^n = \exp\!\Bigl(-n\int_u^v \cot x\,dx\Bigr) \;\le\; \exp\!\Bigl(-n(v - u)\cot\!\Bigl(\frac{u + v}{2}\Bigr)\Bigr).$$
> 3. Because $3(c + j) = 2(2j + c) - (j - c) \le 2(n - 1) - q_{\mathrm{gap}} < 2n$, we have $\frac{c + j}{2n}\pi < \frac{\pi}{3}$, hence
>    $$\frac{R_z}{R_y} = \Bigl(\frac{\sin(c\pi/n)}{\sin(j\pi/n)}\Bigr)^n \;<\; \exp\!\Bigl(-\frac{(j - c)\pi}{\sqrt{3}}\Bigr) \qquad \text{universally.}$$
> 4. For odd $c$, $R_z \le e^{-2\pi/\sqrt{3}} < 0.02660$; for even $c$ and $n \ge 84$, $3c \le n - 4$ gives $R_z \le \exp\bigl(-3\pi\cot((\frac{1}{3} + \frac{1}{504})\pi)\bigr) < 0.004691$.
> 5. $\frac{j + a_1}{2} \le \frac{n - c}{2}$ and $\sin\alpha_1 - \sin\psi_j = 2\sin\bigl(\frac{c\pi}{2n}\bigr)\sin\bigl(\frac{p\pi}{2n}\bigr)$. $\square$

> **Lemma 4.2 (Rouché Isolation & Holomorphic Bounds).**  
> If $\rho_x \in (0, \min(\frac{1}{2}, \frac{1}{2}|1 - t|))$ and $u_y \in (0, \frac{1}{2})$ with $u_y R_y < \frac{1}{2}|1 - t|$ satisfy $\rho_x(1 - \rho_x)(|1 - t| - \rho_x) > R_y R_z \lambda$ and $u_y(1 - u_y)(|1 - t| - u_y R_y) > \frac{R_z}{R_y}\lambda$, then on $|w| \le \lambda$ the cubic $T(T - 1)(T - t) = t R_z w$ has unique simple holomorphic roots $x(w) = 1 + \delta_x(w)$ in $|T - 1| \le \rho_x$ and $y(w) = t(1 + \tilde{\delta}_y(w))$ in $|T - t| \le u_y R_y$, satisfying
> $$|\delta_x| \le \rho_x, \qquad |\delta_y| \le \bar{\delta}_y := \frac{u_y + \rho_x}{1 - \rho_x}, \qquad |\delta_z| \le \bar{\delta}_z := \frac{1}{(1 - \rho_x)^3(1 - \bar{\delta}_y)} - 1, \qquad |\varepsilon_s| \le \frac{L_{\text{dy}}}{n},$$
> where $L_{\text{dy}} := -\log(1 - \bar{\delta}_y)$ (since $|\binom{1/n}{k}| = (-1)^{k-1}\binom{1/n}{k}$ and $1 - e^{-u} \le u$). $\square$

> **Lemma 4.3 (Monotone Riemann Sandwich & Zero-Endpoint $|a_k^{\mathrm{alg}}|$ Bounds).**  
> Let $g(\eta) := \operatorname{Re}\bigl(\frac{-s_j e^{i\eta}}{1 - s_j e^{i\eta}}\bigr) = \frac{s_j^2 - s_j\cos\eta}{1 + s_j^2 - 2s_j\cos\eta}$ and $L_{\text{sum}} := \sum_{\text{large}} a_k^{\mathrm{alg}} = \sum_{\text{large}} g(\eta_k) \in \mathbb{R}$. Then $g'(\eta) = \frac{s_j(1 - s_j^2)\sin\eta}{(1 + s_j^2 - 2s_j\cos\eta)^2} > 0$ on $(0, \pi)$, $g(\pi) = \frac{s_j}{1 + s_j} \in (0, \frac{1}{2})$, $g(\theta_c) = -\frac{r_2\cos\alpha_1}{r_3}$, and $\frac{n}{\pi}\int_{\theta_c}^\pi g(\eta)\,d\eta = \frac{n}{\pi}(\arg(1 + s_j) - \arg(1 - s_j e^{i\theta_c})) = j$. By uniform-spacing monotonicity (both parities of $c$),
> $$j - g(\pi) \;\le\; L_{\text{sum}} \;\le\; j + g(\pi) - 2 g(\theta_c).$$
> Moreover, $|a_k^{\mathrm{alg}}| \le \min\bigl(\frac{r_2}{r_3},\; \frac{1}{2\sin(|\eta_k|/2)}\bigr)$ for every large root, and right-Riemann comparison with zero endpoint correction gives
> $$\frac{1}{n}\sum_{\text{large}} |a_k^{\mathrm{alg}}| \;\le\; \frac{1}{\pi}\log\cot\!\Bigl(\frac{\theta_c}{4}\Bigr), \qquad \frac{1}{n}\sum_{\text{large}} |a_k^{\mathrm{alg}}|^2 \;\le\; \frac{1}{2\pi}\cot\!\Bigl(\frac{\theta_c}{2}\Bigr). \qquad \square$$

> **Lemma 4.4 (Off-Boundary Factor Bounds on $|w| = 1$).**  
> For every $\eta_k \ne \pm\theta_c$, $\frac{|w_k|^2}{R_z^{2/n}} = 1 + \frac{2\kappa}{\pi\sin\theta_c}\sin\bigl(\frac{|\eta_k|+\theta_c}{2}\bigr)\sin\bigl(\frac{|\eta_k|-\theta_c}{2}\bigr)$. With $N_0 = 84$, $c_J = 4\kappa_{\min}/\pi$, and $\eta_{\text{pert}} \ge \frac{r_2}{r_3} L_{\text{dy}}$:
> 1. **Large roots** ($|\eta_k| = \theta_c + \frac{2r\pi}{n}$, $r \ge 1$): $\sin\bigl(\frac{|\eta_k|+\theta_c}{2}\bigr) \ge \sin\theta_c$ and monotonicity of $(1 + A/n)^{-n/2}$ in $n \ge N_0$ yield $\sum_{\text{large}} \bigl|\frac{z/x}{\tilde{w}_k}\bigr| \le S_{\text{lg}}^\sim := S_{\text{lg}}(1 + \bar{\delta}_z)(1 - \frac{\eta_{\text{pert}}}{N_0})^{-N_0}$, where
>    $$S_{\text{lg}} := 2\sum_{r=1}^{2}\!\Bigl(1 + \frac{2\kappa_{\min} r(1 - \frac{r^2\pi^2}{6 N_0^2})}{N_0}\Bigr)^{-N_0/2} + 2\sum_{r=3}^{20}\!\Bigl(1 + \frac{c_J r}{N_0}\Bigr)^{-N_0/2} + \frac{2 N_0}{c_J(N_0/2 - 1)}\!\Bigl(1 + \frac{20 c_J}{N_0}\Bigr)^{-(N_0/2 - 1)}.$$
> 2. **Small roots** ($|\eta_k| = \theta_c - \frac{2r\pi}{n}$, $1 \le r \le \lfloor c/2 \rfloor$): for $\theta_c \le \theta_{\max}$ and $\lambda_s := \operatorname{sinc}(\theta_{\max}/2)\kappa_{\min}$, $|w_k|/R_z \le \exp(-\lambda_s r(1 - r/c))$, giving a closed-form upper bound $S_{\text{sm}}$ (explicitly maximized over $c \in \{2,3,4,5\}$ and $c \ge 6$ via $r(1-r/c) \ge r/2$ for $r \ge 3$). Splitting $\rho_k = |w_k/R_z|^{1/n}$ at $1/2$ yields
>    $$\sum_{\text{small}}\Bigl|\frac{\tilde{w}_k}{z/x}\Bigr| \;\le\; S_{\text{sm}}^\sim := \Bigl(S_{\text{sm}}\, e^{2\eta_{\text{pert}}} + \frac{N_0}{3}\Bigl(\frac{1}{2} + \frac{\eta_{\text{pert}}}{N_0}\Bigr)^{N_0}\Bigr)(1 + \rho_x)^3(1 + \bar{\delta}_y). \qquad \square$$

> **Lemma 4.5 (Exact Cauchy Integral of $B(w)$ & Master Lower Bound).**  
> Under Lemma 4.2 on $|w| \le 1$, differentiating $T(T-1)(T-t) = t R_z w$ at $w = 0$ gives $\delta'_y(0) = -\frac{(1+t^2)R_z}{t(1-t)}$, $\delta'_z(0) = \frac{(1-2t^2)R_z}{t(1-t)}$, and $\varepsilon'_+(0) + \varepsilon'_-(0) = 2 g(\theta_c)\delta'_y(0)$. Cauchy's integral formula on Theorem 4.2's expanded $B(w)$ gives identically
> $$B_{\text{avg}} := \frac{1}{2\pi}\int_{-\pi}^{\pi} B(e^{i\phi})\,d\phi \;=\; 2 - (-1)^j\bigl(2 g(\theta_c)\,\delta'_y(0) - \delta'_z(0)\bigr) \;\in\; \mathbb{R},$$
> together with $\frac{1}{2\pi}\int_{-\pi}^\pi |B(e^{i\phi})|\,d\phi \le 2 + R_{\text{bdy}}$, where $R_{\text{bdy}} := \bar{\delta}_z + 2\varepsilon_{\text{bdy}} + ((1 + \varepsilon_{\text{bdy}})^2(1 + \rho_x)^3(1 + \bar{\delta}_y) - 1)$ ($\varepsilon_{\text{bdy}} := e^{\eta_{\text{pert}}} - 1$). Setting $t_1 \ge d_{\mathrm{deg}}\log(1 + \rho_x)$,
> $$t_2 := |L_{\text{sum}} - c|\, L_{\text{dy}} + \frac{(\frac{1}{n}\sum_{\text{large}}|a_k^{\mathrm{alg}}|)\, L_{\text{dy}}^2}{2(1 - L_{\text{dy}}/n)} + \frac{(\frac{1}{n}\sum_{\text{large}}|a_k^{\mathrm{alg}}|^2)\, L_{\text{dy}}^2}{2(1 - \eta_{\text{pert}}/n)},$$
> and $M_1 := \exp(t_1 + t_2 + S_{\text{sm}}^\sim + S_{\text{lg}}^\sim) - 1$, we obtain $|\Delta_{\text{mult}}(e^{i\phi})| \le M_1$ and the **Master Lower Bound**:
> $$\frac{1}{2\pi}\int_{-\pi}^{\pi}\operatorname{Re}\bigl((1 + \Delta_{\text{mult}}(e^{i\phi}))\, B(e^{i\phi})\bigr)\,d\phi \;\ge\; B_{\text{avg}} - M_1\,(2 + R_{\text{bdy}}). \qquad \square$$

> **Theorem 4.6 (Complete Closed-Form Certification of Regime 1: $\kappa \ge 2.5$, All Even $n \ge 84$).**  
> For every even $n = 2m \ge 84$ and every strictly interior $(c, j)$ ($2 \le c \le \lfloor(n-3)/3\rfloor$, $j \in \mathcal{J}_{n,c}$) with $\kappa(n, c, j) \ge 2.5$, we have $B_{\text{avg}} - M_1(2 + R_{\text{bdy}}) > 0$, and hence $(-1)^j Q_c^{(2m)}(1, u_j) > 0$.

*Proof.* Partition Regime 1 into 6 mutually exhaustive closed-form subcases and apply Lemmas 4.1–4.5 (all elementary inequalities verified with zero search in [`algebra/verify_regime1_complete.py`](algebra/verify_regime1_complete.py) and independently signed off in `#59`, `#61`):
1. **Case 1A ($2 \le c \le 0.15n$, $\kappa \ge \kappa_0$):** If $R_y > 1/2$, then $\cot\psi_j < \frac{2^{1/84}-1}{\sin(2\pi/84)} + \tan(0.075\pi) < 0.3510 \implies \frac{R_z}{R_y} < 0.48114^n \le 2.045\times 10^{-27}$ while $|1 - t| > 4/n$ for both parities of $c$, so Rouché holds with $\rho_x = u_y = 10^{-20}/n$ and all perturbations are $< 10^{-19}$. If $R_y \le 1/2$, then $|1 - t| \ge 1/2$, $\frac{R_z}{R_y} \le 0.002411$, and splitting $\sin\alpha_1 \ge 0.70$ vs $< 0.70$ (where $\kappa \ge \kappa_0$ forces $j \ge c_j\sqrt{n}$) together with $q_{\mathrm{gap}} = 1$ vs $q_{\mathrm{gap}} \ge 2$ yields Margin **`+0.8540`** at $\kappa_0 = 3.0$ and **`+0.0748`** at $\kappa_0 = 2.5$ ($S_{\text{off}}^\sim = 0.6306$, $M_1 = 0.9329$).
2. **Case 1B-I ($c > 0.15n$, $j \ge c + 2$):** $\kappa \ge 2\pi\sin(0.30\pi) = 5.0832$, $\frac{R_z}{R_y} \le e^{-2\pi/\sqrt{3}} = 0.02660$, $R_y^2 \le 0.22125$ (both parities), $R_y R_z \le e^{-\frac{\pi}{2\sqrt{3}} d_{\mathrm{deg}}}$, yielding Margin **`+0.8435`**.
3. **Case 1B-IIa ($0.15n < c \le 0.25n$, $j = c + 1$):** $\kappa \ge 5.0832$, $\frac{R_z}{R_y} \le 0.04579$, $R_y < 7.5\times 10^{-5}$, $\frac{r_2}{r_3} \le 1.073$, yielding Margin **`+1.3934`**.
4. **Case 1B-IIb ($0.25n < c \le 0.31n$, $j = c + 1$):** $\alpha_1 \in (\frac{\pi}{2}, \frac{2\pi}{3}) \implies \kappa > \pi\sqrt{3} = 5.4414$ and $\cos\alpha_1 < 0 \implies g(\theta_c) > 0 \implies L_{\text{sum}} - c < 1.5$; $\frac{R_z}{R_y} \le 0.1160$, $R_y \le 0.00642$, yielding Margin **`+0.6359`**.
5. **Case 1B-III even ($0.31n < c \le \frac{n-3}{3}$, $j = c + 1$, even $c$):** $\kappa > \pi\sqrt{3}$, $L_{\text{sum}} - c \le 1.5 + 2\cos(0.62\pi) = 0.7638$, $d_{\mathrm{deg}} \ge 4 \implies R_y \le 0.03542$, $R_z < 0.004691$, $d_{\mathrm{deg}} R_y R_z < 0.00095$, yielding Margin **`+0.5017`**.
6. **Case 1B-III odd ($0.31n < c \le \frac{n-3}{3}$, $j = c + 1$, odd $c$):** Same $\kappa$ and $L_{\text{sum}} - c \le 0.7638$, $t = -R_y < 0 \implies |1 - t| \ge 1$ and $\frac{1 + t^2}{1 - t}, \frac{|1 - 2t^2|}{1 - t} \le 1$; $d_{\mathrm{deg}} \ge 3 \implies R_y \le 0.1882$, $R_z \le 0.02660$, $d_{\mathrm{deg}} R_y R_z < 0.0200$, yielding Margin **`+0.4699`**.

| Subcase ($\kappa_0 = 3.0$, $n \ge 84$) | $S_{\text{off}}^\sim$ | $t_1$ | $t_2$ | $M_1$ | $R_{\text{bdy}}$ | $B_{\min}$ | **Closed-Form Margin** |
|---|---:|---:|---:|---:|---:|---:|---:|
| **1A** ($c \le 0.15n$, $\kappa \ge 3.0$) | 0.4125 | 0.0011 | 0.0255 | 0.5513 | 0.0416 | 1.9796 | **`+0.8540`** *(at $\kappa_0=2.5$: **`+0.0748`**)* |
| **1B-I** ($c > 0.15n$, $j \ge c+2$) | 0.0523 | 0.0545 | 0.2292 | 0.3993 | 0.5785 | 1.8732 | **`+0.8435`** |
| **1B-IIa** ($0.15n < c \le 0.25n$, $j=c+1$) | 0.0454 | 0.0000 | 0.1467 | 0.2118 | 0.3460 | 1.8903 | **`+1.3934`** |
| **1B-IIb** ($0.25n < c \le 0.31n$, $j=c+1$) | 0.0410 | 0.0001 | 0.2617 | 0.3537 | 1.1844 | 1.7623 | **`+0.6359`** |
| **1B-III even** ($c > 0.31n$, $j=c+1$, even) | 0.0495 | 0.0010 | 0.2107 | 0.2984 | 1.8726 | 1.6574 | **`+0.5017`** |
| **1B-III odd** ($c > 0.31n$, $j=c+1$, odd) | 0.0492 | 0.0206 | 0.2023 | 0.3126 | 1.8382 | 1.6699 | **`+0.4699`** |

*(All 17 intermediate lemma inequalities were also cross-checked across all $13{,}006$ strictly interior points for $n \in \{84, 86, 88, 90, 100, 160, 320\}$ with zero failures.)* $\square$

---

## 5. Regime 2 ($\kappa < 2.5$, All Even $n \ge 132$): Certified & Signed-Off Components

Throughout Sections 5 and 6, fix even $n = 2m \ge 132$, strictly interior $2 \le c \le \lfloor (n - 3)/3 \rfloor$, $j \in \mathcal{J}_{n, c}$ ($j \ge c + 1$), and $\kappa(n, c, j) < 2.5$. All interval computations use **50-digit outward-rounded interval arithmetic (`mpmath.iv`, `iv.dps = 50`)** with pure `mpf` endpoints (zero float64 narrowing and zero loops over $(n, c, j)$). Denote the unperturbed integrand of Theorem 4.2 by
$$F_{\text{pure}}(w) := B^{\text{pure}}(w)\, F_3^{\text{pure}}(w), \qquad B^{\text{pure}}(w) := 2 - Z - Z^{-1}, \qquad F_3^{\text{pure}}(w) := \prod_{\text{small}}(1 - u_k Z)\prod_{\text{large}}(1 - u_k^{-1} Z^{-1}),$$
where $u_k := (-1)^j R_z^{-1}(1 - s_j e^{i\eta_k})^n$ and $a_k := u_k^{-1}$ (note $a_k$ here denotes the reciprocal scaled root, distinct from $a_k^{\mathrm{alg}}$ in Theorem 4.2).

### 5.1 Lemma R2.1: Universal Regime 2 Cubic-Root Perturbation & Stretched-Exponential Decay

> **Lemma 5.1 (Lemma R2.1 — Universal Regime 2 Perturbation; signed off `#74`).**  
> For every $(n, c, j)$ in Regime 2 ($n \ge 132$, $2 \le c \le \lfloor(n-3)/3\rfloor$, $j \in \mathcal{J}_{n,c}$, $\kappa < 2.5$):
> 1. **Universal angle & separation bounds:** $\alpha_1 < \alpha_{\max} := 0.409284$, $\theta_c < \alpha_{\max}/2$, $S_\alpha := \frac{\sin\alpha_{\max}}{\alpha_{\max}} \ge 0.972314$, $\cos(2\alpha_{\max}) \ge 0.683268$, and
>    $$R_z \le 2.971\times 10^{-39}, \qquad R_y \le 5.121\times 10^{-7}, \qquad \frac{R_z}{R_y} \le 2.668\times 10^{-7},$$
>    $$\frac{R_z}{R_y} \le \exp\!\Bigl(-\frac{C_{\exp}}{\kappa}\Bigr) \quad (C_{\exp} \ge 37.8481), \qquad \frac{1}{\theta_c}\frac{R_z}{R_y} \le 2.025\times 10^{-4}.$$
> 2. **Holomorphic perturbation on $|w| = 1$:** Lemma 4.2 holds with $\rho_x = 10^{-40}$ and $u_y = 1.001\, R_z/R_y$, certifying $t_1 \le 4.017\times 10^{-43}$, $t_2 \le 1.000\times 10^{-6}$, $\eta_{\text{pert}} \le 2.998\times 10^{-7}$, $R_{\text{bdy}} \le 1.733\times 10^{-6}$, $C_{\text{exp,pert}} := \frac{131}{132} C_{\exp} \ge 37.5614$, $R_{\text{bdy}} \le 6.0061\, (R_z/R_y)^{131/132}$, and uniform off-boundary coefficient $\eta_{3,\text{coeff}} \le 2.0996$.
> 3. **Stretched-exponential decay inside $c\kappa^2 < 42$:** Since $j \ge c + 1$ implies $j^2(j+c)^2/c > \frac{1}{2}(j+c)^3$, we have $\frac{\alpha_1}{2} \le K_{jc} n^{-1/3}$ ($K_{jc} := (\frac{42}{2\pi^4 S_\alpha^4})^{1/3}$), yielding
>    $$\frac{R_z}{R_y} \;\le\; \exp\bigl(-c_{\text{se}}\, n^{1/3}\bigr) \quad (c_{\text{se}} \ge 3.16794), \qquad n\,\frac{R_z}{R_y} \;\le\; 1.304\times 10^{-5} \quad (\forall n \ge 132).$$

*Proof.* Items (1)–(3) follow from Lemma 4.1, monotonicity of $x\cot x$ on $(0, \pi/2)$, and the algebraic inequality $2j^2 - c(j + c) = (2j + c)(j - c) > 0$. Certified in 50-digit `mpmath.iv` by [`analyst/certify_R2_1_and_2A_iv.py`](analyst/certify_R2_1_and_2A_iv.py). $\square$

### 5.2 Sub-Regime 2A ($c\kappa^2 \ge 42$, $\kappa < 2.5$, $n \ge 132$)

> **Theorem 5.2 (Lemma R2.2 — Sub-Regime 2A; signed off `#74`).**  
> For every strictly interior $(n, c, j)$ with $n \ge 132$, $\kappa < 2.5$, and $c\kappa^2 \ge 42$, the Theorem 4.2 integral on $|w| = 1$ satisfies
> $$\frac{\operatorname{Re}\int_{-\pi}^{\pi} (1 + \Delta_{\text{mult}}(e^{i\phi}))\, B(e^{i\phi})\,d\phi}{\int_{-\pi}^{\pi} (2 - 2\cos(\phi + j\pi))\, |F_3^{\text{pure}}(e^{i\phi})|\,d\phi} \;\ge\; \mathbf{+0.220613} \;>\; 0,$$
> and therefore $(-1)^j Q_c^{(2m)}(1, u_j) > 0$.

*Proof.* On $|w| = 1$, $B^{\text{pure}}(e^{i\phi}) = 4\sin^2((\phi + j\pi)/2) \ge 0$. Pairing each small root $\eta_k = \pm(\theta_c - \frac{2m\pi}{n})$ ($1 \le m \le \lfloor c/2 \rfloor$) with its mirror large root $\eta_{k'} = \pm(\theta_c + \frac{2m\pi}{n})$ cancels the leading linear phase terms in $\arg F_3^{\text{pure}}(e^{i\phi})$ up to $O(m^2/(c\kappa))$. Splitting the index at $m \le 0.26\, c$ vs $m \in (0.26\, c, \lfloor c/2 \rfloor]$ (plus any unpaired center root and unpaired large roots $m > c/2$) certifies
$$|\arg F_3^{\text{pure}}(e^{i\phi})| \;\le\; B_{2A} \;\le\; 1.332945 \;<\; \frac{\pi}{2} \qquad \implies \qquad \cos(B_{2A}) \;\ge\; 0.235615 \;>\; 0$$
uniformly on $|w| = 1$. Bounding $|\Delta_{123}| := |F_1 F_2 (F_3 / F_3^{\text{pure}}) - 1| \le 8.246\times 10^{-6}$ and controlling the oscillation ratio $\max_{|w|=1}|F_3^{\text{pure}}| / \min_{|w|=1}|F_3^{\text{pure}}| \le \exp(C_{\text{osc,tot}}/\kappa + c_{\text{osc}})$ via $\int_0^\infty \log\frac{1+e^{-t}}{1-e^{-t}}\,dt = \frac{\pi^2}{4}$ gives $C_{\text{osc,tot}} \le 21.175 \ll C_{\text{exp,pert}} \ge 37.561$ and $c_{\text{osc}} \le 1.2549$, yielding boundary ratio $\le 0.01499$. Subtracting gives $\cos(B_{2A}) - |\Delta_{123}| - (1 + |\Delta_{123}|)\cdot 0.01499 \ge +0.220613 > 0$. Certified in 50-digit `mpmath.iv` by [`analyst/certify_R2_1_and_2A_iv.py`](analyst/certify_R2_1_and_2A_iv.py). $\square$

### 5.3 Sub-Regime 2C-Upgrade ($c q \kappa \le 0.90$, $\kappa < 2.5$, $n \ge 132$)

> **Theorem 5.3 (Lemma R2.3 — Sub-Regime 2C-Upgrade; signed off `#77`).**  
> Set $d := \gcd(c, j) \ge 1$ and $q := c/d \ge 1$. For every strictly interior $(n, c, j)$ with $n \ge 132$, $c \ge 2$, and $c q \kappa \le 0.90$, shifting the contour of Theorem 4.2 to the exact zero-exterior-root radius $|w| = \lambda_* \in (0, 1)$ certifies
> $$\frac{\operatorname{Re}\int_{-\pi}^{\pi} F_{\text{full}}(\lambda_* e^{i\phi})\,d\phi}{\int_{-\pi}^{\pi} |F_{\text{main}}(\lambda_* e^{i\phi})|\,d\phi} \;\ge\; \begin{cases} \mathbf{+0.441480} \;>\; 0 & (q = 1), \\ \mathbf{+0.266331} \;>\; 0 & (q \ge 2), \end{cases}$$
> and therefore $(-1)^j Q_c^{(2m)}(1, u_j) > 0$.

*Proof.* Inside $c q \kappa \le 0.90$:
1. **Sharp automatic geometry & cubic perturbation:** $\kappa \le 0.45$ ($q = 1$) and $\kappa \le 0.225$ ($q \ge 2$), enforcing $n \ge 525$ ($q = 1$), $n \ge 657$ ($q \ge 2$), $c\kappa^2 \le 0.405$, and stretched-exponential constant $c_{\text{se},2C} \ge 15.6461$. Choosing $\lambda_*$ so that $W_0 := \lambda_*^q \max_k |a_k|^q < 1$ satisfies the zero-first-moment saddle equation $\frac{1}{q}\sum_k \frac{x_k}{1 + x_k} = d$ ($x_k := |a_k|^q \lambda_*^q$), a 3-line pole-free contour-shift Poisson evaluation at $y_{\text{sh}} = 0.80$ certifies $W_0 < 0.836$ ($q = 1$) and $W_0 < 0.832$ ($q \ge 2$), whence $1 - |a_k|\lambda_* \ge (-\log W_0) W_0 / n$ and $E_{\text{cubic},2C} \le 6.89\times 10^{-48}$, while the peak curvature satisfies $H_q := \frac{1}{q}\sum_k \frac{x_k}{(1+x_k)^2} \ge \frac{d}{1 + W_0} \ge 0.54466\, d$.
2. **Pólya–Szegő positive-root ratio (`rat_main >= 0.473803`):** Every factor of the $q$-periodic main envelope $F_{\text{main}}$ has positive real roots, yielding $\operatorname{Re}\log(F_{\text{main}}(e^{i\theta})/F_{\text{main}}(1)) \in [4 H_q \log\cos(\theta/2),\, -2 H_q \sin^2(\theta/2)]$ and $|\operatorname{Im}\log F_{\text{main}}(e^{i\theta})| \le 2 H_q(2\tan(\theta/2) - \theta)$. Closed-form Gaussian bounds ($H_q \ge 5.0$: $\ge 0.53777$) plus a 150-box 50-digit `iv` enclosure ($H_q \in [0.544, 5.0]$) certify $\text{rat}_{\text{main}}(H_q) \ge 0.473803$.
3. **Non-$q$ harmonics ($|R_{\text{non-}q}| \le 0.174420$ for $q \ge 2$, identically $0$ for $q = 1$) & conjugate-pair phase ($E_{\text{phase}}$):** Expanding $\log(F_{\text{pure}}/F_{\text{main}})$ in powers of $w^r$ ($q \nmid r$), exact negative-binomial Cauchy evaluation + Kullback–Leibler divergence from the Stirling peak over blocks $l = 0, \dots, 47$ plus geometric tails certifies $|R_{\text{non-}q}| \le 0.174420$. Pairing conjugate roots $\pm k$ cancels linear phase drift and converts the quadratic phase error into polylogarithms $\operatorname{Li}_{p-1/2}(W_0), \operatorname{Li}_{p-1}(W_0)$, certifying $E_{\text{phase}} \le 0.031812$ ($q = 1$) and $\le 0.014110$ ($q \ge 2$).
4. **Assembly:** $\text{rat}_{\text{main}} - (\exp(|R_{\text{non-}q}| + E_{\text{phase}} + E_{\text{cubic},2C}) - 1) \ge +0.441480$ ($q = 1$) and $\ge +0.266331$ ($q \ge 2$). Certified in 50-digit `mpmath.iv` by [`analyst/certify_R2_3_iv.py`](analyst/certify_R2_3_iv.py). $\square$

### 5.4 Piece G2: Complete 50-Digit Certified Proof of $D = 0$ Major Arcs (Signed Off `#90`, `#103`)

**Domain (Box G2):** $n \ge 132$, $\kappa < 2.5$, $c\kappa^2 < 42$, $d = \gcd(c, j) \ge 2$, $q = c/d \ge 1$, $p_0 = j/d$ ($\gcd(p_0, q) = 1$), $X := q^2\kappa \le K_{\text{G2}, q=1} := 2.5$ (for $q = 1$, covering ALL of Regime 2 for $q = 1$) and $X \le K_{\text{G2}} := 2.0$ (for $q \ge 2$), and $c q \kappa = d X \ge 0.85$ (so $T := \frac{1}{2} c q \kappa = \frac{1}{2} d X \ge 0.425$, strictly overlapping Theorem 5.3 on $cq\kappa \in [0.85, 0.90]$).

Index the $n$ roots by $x_k = k + \text{shift} \in (-n/2, n/2]$ ($\text{shift} = 0$ for even $c$, $\frac{1}{2}$ for odd $c$, so $x_k - c/2 \in \mathbb{Z}$ for both parities), with $\eta(x) = \frac{2\pi x}{n}$ and $a_k = b(x_k)\, e^{2\pi i p_0(x_k - c/2)/q}$, where
$$b(x) := \Bigl(\frac{r_3/r_1}{1 - s_j e^{i\eta(x)}}\Bigr)^n e^{-2\pi i\frac{j}{c}(x - c/2)}(-1)^j = \exp\bigl(\log|b(x)| + i\,\delta\psi(x)\bigr), \qquad b(\pm c/2) = 1.$$

> **Proposition 5.4 (Universal $D = 0$ Contour Identity — All Parities of $(d, q, p_0)$).**  
> Set $z_{\text{var}} := Z^{-1} = (-1)^j w$ along $|z_{\text{var}}| = \lambda$, partitioned into $q$ equal sectors centered at the roots of $z_{\text{var}}^q = -\lambda^q$:
> $$z_{\text{var}} = \lambda \exp\!\Bigl(i\,\frac{(2m + 1)\pi + \psi}{q}\Bigr) \qquad (m = 0, 1, \dots, q - 1,\quad \psi \in [-\pi, \pi]).$$
> Then for **every** even $n$, every $c = d q$, and every $j = d p_0$ with $\gcd(p_0, q) = 1$ (all parities of $d, q, p_0$):
> $$F_{\text{pure}}(w) \;=\; \Bigl(\prod_{\text{small}} |u_k|\Bigr)\,\lambda^{-c}\, e^{-i d\psi}\,\prod_{k=0}^{n-1}\bigl(1 - a_k z_{\text{var}}\bigr).$$
> Within each sector $m$, the $q$ residue classes of $x_k - c/2 \pmod q$ permute bijectively onto the $q$ rays $\theta_t(\psi) \equiv \frac{(2t+1)\pi + \psi}{q} \pmod{2\pi}$ ($t = 0, \dots, q-1$), satisfying $\prod_{t=0}^{q-1}(1 - b\lambda e^{i\theta_t(\psi)}) = 1 + (b\lambda)^q e^{i\psi}$.

*Proof.* Rewriting $1 - u_k Z = -u_k z_{\text{var}}^{-1}(1 - a_k z_{\text{var}})$ on all $c - 1$ small roots and $(1 - Z)^2 = z_{\text{var}}^{-2}(1 - z_{\text{var}})^2$ produces the exact prefactor $(-1)^c(\prod_{\text{small}} u_k) z_{\text{var}}^{-c}$. Checking both parities of $q$: when $q$ is even, $p_0$ is odd, $c$ is even, $j \equiv d \pmod 2$, and the unique unpaired small root $x = 0$ has $u_0 = (-1)^j|u_0| = (-1)^d|u_0|$, matching $z_{\text{var}}^{-c} = (-1)^d \lambda^{-c} e^{-i d\psi}$; when $q$ is odd, $c \equiv d \pmod 2$, giving $(-1)^c \prod_{\text{small}} u_k = (-1)^d \prod |u_k|$ for both parities of $d$, again cancelling the $(-1)^d$ from $z_{\text{var}}^{-c}$. Since $\gcd(p_0, q) = 1$, $p_0(x_k - c/2) \pmod q$ permutes all $q$ residues. $\square$

The complete certification of Piece G2 proceeds through seven self-contained lemmas, each executed in pure 50-digit `mpmath.iv` with zero hardcoded constants in [`algebra/certify_G2_D0.py`](algebra/certify_G2_D0.py):

- **Lemma G2.1 (Modulus, Explicit $Q_{\min} > 0.20$ Cutoffs, Global UHP Curvature, and Phase Laws):**  
  With $Q(x) := \frac{|1 - s_j e^{i\eta(x)}|^2}{(r_3/r_1)^2} - 1$ and $v_G(x) := \frac{\kappa}{c}|\frac{c^2}{4} - x^2|$, 50-digit `iv` certifies $Q(\rho n) > 0.20$ at $\rho = 0.100$ (all Regime 2: $Q \ge 0.2644$, dead tail $\sum_{Q \ge 0.20}|b|^q \le 132\cdot 1.20^{-66} < 0.000785$), $\rho = 0.050$ (Box G2, $q = 1$, $\kappa < 2.5$, $\alpha_1 < 0.065\pi$: $Q \ge 0.5548$, $c_{\text{in}} \in [0.999035, 1.01458]$, $c_{\text{out}} \in [0.903784, 1.00077]$, $C_{\psi, q1} \le 0.0473395$, $C_{u'', q1} \le 1.04309$), and $\rho = 0.060$ (Box G2, $q \ge 2$: $Q \ge 0.2775$, $c_{\text{in}} \in [0.999669, 1.00185]$, $c_{\text{out}} \in [0.900742, 1.00026]$, $C_{\psi, q2} \le 0.080404$, $C_{u'', q2} \le 1.00645$). Moreover, $\beta(z) := b'(z)/b(z)$ obeys the global upper-half-plane bound $|\beta'(x + i y)| \le C_{u''}\frac{2\kappa}{c}$ for all $x \in \mathbb{R}, y \ge 0$, and $0 \le -\delta\psi''(x) \le A_{\text{in/out}}\frac{\kappa^2}{c^2} x$ on $[0, n/2]$.

- **Lemma G2.2 (Saddle $\Lambda_* = \lambda_*^q$ & Curvature $H_q^+$ via 2nd-Order Euler–Maclaurin):**  
  Define $G_q^+(W) := \prod_{k=0}^{n-1}(1 + |b(x_k)|^q W)$ and let $\Lambda_* = \lambda_*^q > 0$ solve $M_q^+(\Lambda_*) := \frac{1}{q}\sum_{k=0}^{n-1}\sigma(u_k) = d$ with curvature $H_q^+(\Lambda) := \frac{1}{q}\sum_{k=0}^{n-1} h(u_k)$ ($u_k := q\log(\lambda|b(x_k)|)$, $\sigma(u) = (1+e^{-u})^{-1}$, $h(u) = e^u/(1+e^u)^2$). Using $|B_2(\{x\})/2| \le 1/12$, the global bound $|u''(x)| \le C_{u''}\frac{2q\kappa}{c}$, and 50-digit `iv` enclosures $J_{\sqrt{}}(\sigma) \le 0.49969$, $J_{\sqrt{}}(h) \le 0.41158$ (with explicit $\Lambda\cdot\text{tail\_sum\_Q20}$), direct discrete sums on $c \in \{2, \dots, 10\}$ (which cover every $q = 1$ point with $\kappa \ge \sqrt{42/11} \approx 1.9540$) plus 30 narrow $T$-bins on $[0.425, 14.0]$ and an exact FTC tail on $T \ge 14.0$ certify:
  - *Subcase 2a ($q = 1, c \in \{2, \dots, 10\}$):* $H_1^+(\Lambda_*) \ge 0.68924$.
  - *Subcase 2b ($q = 1, c \ge 11 \implies \kappa < \sqrt{42/11} < 1.9541 < 2.0$):* $\Lambda_* \in [0.5717, 2.162]$, $H_{\text{cont,end}} \ge 0.84276$, $\max\text{disc}_H \le 0.24292 \implies H_1^+(\Lambda_*) \ge 0.62020$.
  - *Subcase 2c ($q \ge 2, d \ge 2$):* $\Lambda_* \in [0.5519, 1.803]$, $H_{\text{cont,end}} \ge 0.85844$, $\max\text{disc}_H \le 0.095147 \implies H_q^+(\Lambda_*) \ge 0.59428$.

  Pre-widening every $\Lambda$-interval by $e^{\pm\delta_{\text{rad}}}$ ($\delta_{\text{rad}} = 0.35\cdot 10^{-4} \ge |\log(\tilde{\Lambda}/\Lambda_*)|$) transfers all bounds literally to the Pigeonhole-Annulus radius $\tilde{\lambda}$ (Lemma G2.7), with $H_1^+(\tilde{\Lambda}) \ge 0.62018$, $H_q^+(\tilde{\Lambda}) \ge 0.59425$, and $H_{\text{cont,end}}(\tilde{\Lambda}) \ge 0.85841$.

- **Lemma G2.3 ($P^2$-Monotone Conjugate-Pair Phase Factorization $E_{\text{phase}}$ at $\tilde{\Lambda}$):**  
  Pairing conjugate roots $(x_k, -x_k)$ of $G_q(W) := \prod_{k=0}^{n-1}(1 + b(x_k)^q W)$ yields $G_q(\tilde{\Lambda} e^{i\psi})/G_q^+(\tilde{\Lambda} e^{i\psi}) = E_{\text{phase}}(\tilde{\Lambda} e^{i\psi})$ with $|E_{\text{pair}, k} - 1| = \frac{\sin^2(q\delta\psi_k/2)}{\sinh^2(u_k/2) + \cos^2(\psi/2)}$. Applying $\sin^2(P) \le P^2$ (preserving exact linear $X$-monotonicity everywhere) and an exact $[0, 8.8] \cup [8.8, 0.40 T] \cup \{t \le \sqrt{0.20}\}$ tail split for $T \ge 22$ certifies:
  - *$q \ge 2$ ($|\psi| \le 2\pi/3$):* $\Delta_{\max} \le 0.011303$, $\bigl|\frac{1}{q}\log E_{\text{phase}}\bigr| \le 0.0210522$.
  - *$q = 1$ ($|\psi| \le 1.75$):* $\Delta_{\max} \le 0.039652$, $|\log E_{\text{phase}}| \le 0.121256$; AM–GM on conjugate pairs yields $S_{\text{AMGM}} \le 0.0887442$, giving $\varepsilon_{\text{mid}} \le 0.0443721$ on $|\psi| \in [1.75, 2\pi/3]$ and far-tail decay $c_{\text{far}} = 1.50$, $\varepsilon_{\text{far}} \le 0.0887442$ on $|\psi| \in [2\pi/3, \pi]$.

- **Lemma G2.4 ($\pi$-Line Meromorphic Residue & Cauchy–Riemann Level-Curve Bound on $\mathcal{R}_{\text{non-}q}$, $q \ge 2$):**  
  Write $\log\prod_{k=0}^{n-1}(1 - a_k z_{\text{var}}) = \frac{1}{q}\log G_q(\tilde{\Lambda} e^{i\psi}) + \mathcal{R}_{\text{non-}q}(z_{\text{var}})$ and apply step-$q$ Poisson summation + one integration by parts to each Fourier mode $\ell \ne 0$, replacing branch cuts by the meromorphic integrand $\frac{\beta(z) u(z)}{1 - u(z)} e^{2\pi i \ell z/q}$ ($u(z) = \tilde{\lambda} b(z) e^{i\theta_t}$). Writing $\log b(x+iy) = U(x,y) + i V(x,y)$, since $U_x(x,y) < 0$ strictly on $(0, n/2) \times [0, \infty)$, the level curve $U(x_*(y), y) \equiv -\log\tilde{\lambda}$ obeys the exact Cauchy–Riemann monotonicity law
  $$\frac{d}{dy}\, V(x_*(y), y) \;=\; \frac{U_x^2 + U_y^2}{U_x} \;<\; 0 \qquad \forall y \ge 0.$$
  Hence shifting to $H_\pi$ ($-V = \pi$) on Case 1 rays ($\theta_t \le 5\pi/6$ for $m = 0$, $\theta_t \le 0.72\pi$ for $m \ge 1$) encloses **exact zero count $= 1$**, while shifting to $H_{3\theta_t/4}$ ($-V = \frac{3}{4}\theta_t$) on Case 2 rays ($\theta_t \in (0.72\pi, \pi]$) encloses **exact zero count $= 0$**. Combining the global UHP upper bound $U(x, y) \le (1 - C_{\text{ratio}}\frac{\kappa}{c} y) U(x, 0) + C_\psi\kappa^2 y + C_{u''}\frac{\kappa}{c} y^2$ (split cleanly between $q = 2$, where $r_3/r_2 \le 0.668$, and $q \ge 3$, where $\kappa \le 2/9$), explicit $\sqrt{A}$-rescaling for all $A \ge A_0$, and zero-hardcoded top-line/side-line/higher-ray bounds certifies $E_H(\pi q \ge 2\pi) \ge 8.30$, $E_{H, C2}(\pi) \ge 5.28$, `line_C1 <= 0.0156888`, and
  $$|\mathcal{R}_{\text{non-}q}(z_{\text{var}})| \;\le\; 0.115416 \quad (|\psi| \le 1.25), \qquad |\mathcal{R}_{\text{non-}q}(z_{\text{var}})| \;\le\; 0.333745 \quad (|\psi| \le 2\pi/3).$$

- **Lemma G2.5 (Chebyshev-Mediant Far-Tail Decay on $|\psi| \in [2\pi/3, \pi]$, $q \ge 2$):**  
  Writing $|\psi| = \pi - \beta$ ($\beta \in [0, \pi/3]$) and replacing only the single near-zero ray $\theta_{t_0} = \beta/q$ by $\pi/(3q)$ restores $|\mathcal{R}_{\text{non-}q}| \le 0.333745$ while multiplying each root by $e^{f(u_k; q, \beta)}$, bounded via $\frac{\sinh(qt)}{\sinh t} \ge q \ge \frac{\sin(q\theta)}{\sin\theta}$ by $f_{\max}(u) := \frac{1}{2}\log\bigl((\frac{\tanh(u/2)}{u/2})^2(\frac{u^2}{4} + \frac{\pi^2}{36})\bigr)$. Splitting at $u_0 = 1.80$ ($c_0 = 1.52$, 360 `iv` subintervals) and wiring $H_{\text{cont,end}}(\tilde{\Lambda}) \ge 0.85841$ certifies $\max(\text{tail}/H) \le 0.66073$, $c_{\text{far}} \ge 0.85927$, and $\varepsilon_{\text{far,disc}} \le 0.11787$.

- **Lemma G2.7 (Exact Multiplicative Pigeonhole-Annulus Cubic Perturbation, ALL $n \ge 132$):**  
  Combining all $n$ roots into the purely multiplicative ratio $\frac{F_{\text{full}}(w)}{F_{\text{pure}}(w)} = F_1(w) F_2(w) F_{\text{pref}}(w) \prod_{k=0}^{n-1}\frac{1 - \tilde{a}_k z_{\text{var}}}{1 - a_k z_{\text{var}}}$ removes additive $\Delta_B$ and $\max/\min$ oscillation ratios completely. Since $d \ge 2$ forces $j - c \ge 2$, Lemma 4.1(2) improves stretched-exponential decay inside Box G2 to $R_z/R_y \le \exp(-c_{\text{se,G2}}\, n^{1/3})$ with $c_{\text{se,G2}} \ge 6.5902$. By the Pigeonhole Principle, $[\lambda_*(1 - \varepsilon_{\text{rad}}/n), \lambda_*(1 + \varepsilon_{\text{rad}}/n)]$ ($\varepsilon_{\text{rad}} = 10^{-4}$) contains a radius $\tilde{\lambda}$ separated from all $n$ root moduli $\{|a_k|^{-1}\}$ by at least $\frac{\varepsilon_{\text{rad}}\lambda_*}{n(n+1)}$, certifying
  $$\max_{|w| = \tilde{\lambda}}\Bigl|\frac{F_{\text{full}}(w)}{F_{\text{pure}}(w)} - 1\Bigr| \;\le\; E_{\text{cubic}} := \begin{cases} 1.8182\times 10^{-12} & (q = 1), \\ 1.6158\times 10^{-31} & (q \ge 2), \end{cases}$$
  and universal linear saddle-shift phase drift $\varepsilon_{\text{shift}} \le -\frac{\pi N_0}{4}\log(1 - \varepsilon_{\text{rad}}/N_0) \le 7.854\times 10^{-5}$.

> **Theorem 5.5 (Complete 50-Digit Certification of Piece G2 — Full Exact Integral; signed off `#90`, `#103`).**  
> Combining Lemmas G2.1–G2.7 with the exact positive-root modulus/phase laws for $1 + e^{u + i\psi}$ (using local per-bin $\mathcal{R}_{\text{non-}q}(\psi)$ on $I_{\text{pos}}$ and $I_{\text{mid}}$, subdivided up to $H = 60$ plus analytic Gaussian scaling for $H \ge 60$, and subtracting $E_{\text{cubic}}$ directly in [`algebra/certify_G2_D0.py`](algebra/certify_G2_D0.py)), the **full exact cubic-root integral** on $|w| = \tilde{\lambda}$ satisfies for EVERY $(n, c, j)$ in Box G2:
> $$\frac{I_{\text{full}}(\tilde{\lambda})}{L^1(\tilde{\lambda})} \;\ge\; \begin{cases} \mathbf{+0.216651} \;>\; 0 & (q = 1,\; \kappa < 2.5,\; H_1^+(\tilde{\Lambda}) \ge 0.62018), \\ \mathbf{+0.00176348} \;>\; 0 & (q \ge 2,\; q^2\kappa \le 2.0,\; H_q^+(\tilde{\Lambda}) \ge 0.59425), \end{cases}$$
> and therefore $(-1)^j Q_c^{(2m)}(1, u_j) > 0$. Furthermore, all 28 exact finite-$n$ stress corners covering all four parity combinations of $(d, q)$ (including the $q = 1, \kappa \in (2.0, 2.5)$ band) pass Part B of [`algebra/certify_G2_D0.py`](algebra/certify_G2_D0.py) with empirical $\min I_{\text{pure}}/L^1 = +0.97916$. $\square$

### 5.5 Piece G1: Complete 50-Digit Certified Proof of Dirichlet Minor Arcs (Signed Off `#101`)

**Domain (Box G1):** Every $(n, c, j)$ in Regime 2 ($n \ge 132$, $2 \le c \le \lfloor(n-3)/3\rfloor$, $j \in \mathcal{J}_{n,c}$, $\kappa = 2\pi\frac{r_1 r_2}{r_3} < 2.5$, $c\kappa^2 < 42$) with $q := c/\gcd(c, j) \ge 2$, $X := q^2\kappa$ satisfying $X \ge 0.88$ (when $d = 1$, overlapping Sub-regime 2C-Upgrade's $cq\kappa \le 0.90$) or $X \ge 2.0$ (when $d \ge 2$, overlapping Piece G2's $q^2\kappa \le 2.0$), such that **no** coprime pair $(p_0, r_0)$ with $D := r_0 j - c p_0 \ne 0$ satisfies all four inner Piece G3 inequalities simultaneously: (1) $\frac{c}{r_0 D^2} \ge 10$, (2) $r_0\kappa \le 0.45$, (3) $r_0^3 c\kappa^2 \le 60$, and (4) $\frac{c^3 r_0 \kappa^2}{D^4} \ge 45$ (strictly overlapping Piece G3's $\ge 40$).

On the unit circle $|w| = 1$, where $B^{\text{pure}}(e^{i\phi}) = 4\sin^2((\phi + j\pi)/2) \ge 0$, write $\Phi^{\text{ex}}(\phi; 1) := \arg F_3^{\text{pure}}(e^{i\phi})$. By pairing mirror small and large roots and expanding in Fourier harmonics $r \ge 1$, the Gaussian-model phase decomposes via the continued fraction $\alpha := \min(\{j/c\}, 1 - \{j/c\}) = [0; a_1, \dots, a_s] \in (0, 1/2]$ ($a_1 \ge 2$, $a_s \ge 2$). Five pure 50-digit `mpmath.iv` lemmas in [`algebra/certify_G1_minor.py`](algebra/certify_G1_minor.py) (documented in [`algebra/G1_NOTES.md`](algebra/G1_NOTES.md)) establish uniform positivity:

- **Lemma G1.1 (Zero-Cancellation Hubbard–Stratonovich Duality Kernel $H_D(t)$ & $D_0(T)$):** Pairing $+u$ and $-u$ cancels the singular leading terms identically into $2\pi^2 \nu^2 D(t\pi\nu, \nu) = \mathbb{E}_{Z \sim \mathcal{N}(0, 1)}[Z^2 F(t, \sqrt{t/(\pi\nu)}\, Z)]$ with $|F(t, u)| \le 14$ globally, certifying in 50-digit `iv`: $\sup_{t > 0} H_D(t) \le 3.14341$, $\sup_{t \ge 1.50} H_D(t) \le 0.553496$, $\sup_{t > 0} H_D(t)/t^2 \le 26.8497$ ($\le 24.4591$ on $t \le 0.05$), $\sup_{u > 0} u^2 H_D(u) \le 9.3423$ ($\le 8.45199$ on $t \ge 12$), and centered-form $D_0(T)$ tails $C_{D_0}(T \ge 30) \le 4.34617$, $C_{D_0}(T \ge 13) \le 7.72465$.
- **Lemma G1.2 (Continued-Fraction Minor-Arc Phase Bound $S^{\text{Gauss}} + S_{k \ne 0}$):** Using symmetry normalization $a_1 \ge 2$, the 1-line $D_r^2$-cancellation identity $\frac{1}{4\pi^2}\frac{c}{r D_r^2} H_D(t_r) = \frac{t_r^2 H_D(t_r)}{c\kappa^2 r^3} \le \frac{U_{2H}}{c\kappa^2 r^3}$ for all high levels $k \ge k_2 + 1$ (simultaneously across all four minor-arc contrapositive branches), exact Fibonacci-duality unconditional `REST(k)` bounds, paired $\theta_+ + \theta_- = 1$ `S_SetI`, and coupled B1/B2a/B2b/Case-A ladders, 50-digit `iv` certifies across all of Box G1:
  - *Subcase B1 ($d = 1, X \ge 0.88$):* $\le 1.25848$; *Subcase B1 ($d \ge 2, X \ge 2.0$):* $\le 1.13689$;
  - *Subcase B2a ($X > 6\pi$):* $\le 1.25112$; *Subcase B2b ($X > 6\pi$):* $\le 1.23337$; *Case A ($t_0 > 1.50, X > 6\pi$):* $\le 0.915136$.
- **Lemma G1.3 (Exact $m$-Space Finite-$n$ vs Gaussian Phase Difference):** Separating discrete $c \in \{2..19\}$ ($\kappa \ge 0.88/c^2$) and $c \in \{20..39\}$ ($T_1 \ge 7$) from zero-loop $c \ge 20$ ($T_1 \in [0.05, 7]$), $c \ge 40$ ($T_1 \in [7, 50]$), and analytic tails ($T_1 \le 0.05$: $\le 0.232621$; $T_1 \ge 50$: $\le 0.212014$; dead tail: $\le 0.00156907$) certifies
  $$\sup_{\text{Box G1}} \bigl|\Phi^{\text{ex}}(\phi; 1) - \Phi^G(\phi; 1)\bigr| \;\le\; S_{m\text{-space}} \;\le\; 0.277972.$$
- **Lemma G1.4 & Lemma G1.5 (Grand Phase Bound $< \pi/2$ & Unit-Circle Cubic-Root Perturbation):** Combining Lemmas G1.2 and G1.3 yields the uniform global phase bound across all of Box G1:
  $$\sup_{\text{Box G1}} \bigl|\Phi^{\text{ex}}(\phi; 1)\bigr| \;\le\; \mathbf{1.536448} \;<\; \frac{\pi}{2} \approx 1.57080 \qquad \implies \qquad \frac{I_{\text{pure}}(1)}{L^1(1)} \;\ge\; \cos(1.536448) \;\ge\; \mathbf{+0.03434128}.$$
  Finally, using exact $\kappa = \kappa_{\text{trig}} = 2\pi\frac{r_1 r_2}{r_3} < 2.5$ (so $n\kappa \ge 15\pi^2 S_\alpha^2 \implies T_{\text{lg2}} \ge 8.88518$) and the universal Regime 2 Rouché/oscillation bounds (with zero reliance on $c\kappa^2 \ge 42$), Lemma G1.5 certifies $\Delta_{123} \le 5.1501\times 10^{-6}$, $C_{\text{osc,tot}} \le 16.913 \ll C_{\text{exp,pert}} \ge 37.5613$, and $\text{bdy\_ratio} \le 0.00077761$.

> **Theorem 5.6 (Complete 50-Digit Certification of Piece G1 — Full Exact Integral; signed off `#101`).**  
> For every strictly interior $(n, c, j)$ in Box G1, the **full exact cubic-root integral** on $|w| = 1$ satisfies
> $$\frac{I_{\text{full}}(1)}{L^1(1)} \;\ge\; \cos(\Phi_{\max}) - \Delta_{123} - (1 + \Delta_{123})\,\text{bdy\_ratio} \;\ge\; \mathbf{+0.03355852} \;>\; 0,$$
> and therefore $(-1)^j Q_c^{(2m)}(1, u_j) > 0$. Independently audited across all 2,968 lines and executed by `@verifier` (exit code 0, `#101`); all 27/27 exact finite-$n$ Minor-Arc stress corners pass Part B ($\max|\Phi^{\text{ex}}| = 1.1542$, $\min I_{\text{pure}}/L^1 = +0.88956$). $\square$

---

## 6. Piece G3: Complete 50-Digit Certified Proof of $D \ne 0$ Major Arcs (Signed Off `#124`, `#132`)

**Domain (Box G3):** Every $(n, c, j)$ in Regime 2 ($n \ge 132$, $2 \le c \le \lfloor (n-3)/3 \rfloor$, $j \in \mathcal{J}_{n,c}$, $\kappa < 2.5$, $c\kappa^2 < 42$) admitting a coprime Dirichlet approximation $r_0 j = c p_0 + D$ ($\gcd(p_0, r_0) = 1$, $D \in \mathbb{Z}\setminus\{0\}$) with
$$\tilde{c} := \frac{c}{r_0} \ge 10 D^2, \qquad r_0\kappa \le 0.45, \qquad r_0^3 c\kappa^2 \le 60, \qquad \frac{\tilde{c}\,\tau^2}{D^4} \ge 40 \qquad (\tau := c r_0 \kappa).$$
Note that $j \ge c + 1$ and $\tilde{c} \ge 10 D^2$ enforce $p_0 \ge r_0 \ge 1$, $s_j > 1/2$, $\alpha_1 \le 0.07169$, $n \ge 1007$, $R_z/R_y \le 9.01\times 10^{-39}$, and $c_{\text{se,G3}} \ge 3.330$.

The complete 50-digit `mpmath.iv` proof (detailed in Section 4 of [`ANALYTIC_ROUTE.md`](ANALYTIC_ROUTE.md) and certified with zero hardcoded constants across [`analyst/certify_G3_ABC_full_iv.py`](analyst/certify_G3_ABC_full_iv.py) and [`analyst/certify_G3_D_full_iv.py`](analyst/certify_G3_D_full_iv.py)) rests on four components:

1. **Pigeonhole-Annulus Cubic Perturbation (Lemma R2.1(5)):** For any target contour radius $\lambda_* \in [0.05, 3.00]$ (certified in every sub-case below), partitioning $[\lambda_*(1 - \varepsilon_{\text{ph}}), \lambda_*(1 + \varepsilon_{\text{ph}})]$ ($\varepsilon_{\text{ph}} := 10^{-5}/n$) into $n$ equal subintervals yields a midpoint $\tilde{\lambda}$ separated from all $\le n - 1$ forbidden moduli $\{|a_k|^{-1}, 1\}$. Combined with $M_4 := \sup_{n \ge 1007} n^4(R_z/R_y)^{1-1/n} \le 6.14\times 10^{-20}$ and the Poisson-kernel $L^1$ inequality on the degree-$n$ polynomial $w^{c+1} F_{\text{pure}}(w)$, this certifies
   $$\frac{\operatorname{Re}\int_{-\pi}^{\pi} F_{\text{full}}(\tilde{\lambda} e^{i\phi})\,d\phi}{\int_{-\pi}^{\pi} |F_{\text{pure}}(\lambda_* e^{i\phi})|\,d\phi} \;\ge\; \frac{\operatorname{Re}\int_{-\pi}^{\pi} F_{\text{pure}}(\lambda_* e^{i\phi})\,d\phi}{\int_{-\pi}^{\pi} |F_{\text{pure}}(\lambda_* e^{i\phi})|\,d\phi} \;-\; 2.80\times 10^{-6}.$$

2. **Exact Finite-$n$ Phase-Unwrapped Geometry & $(1-w)$-Free IBP Identities (Sub-Lemmas G3.1–G3.2):** On $u \in [-L_n, L_n]$ ($L_n := \frac{n}{2c}$), define $\zeta(u) := r_0 n\log\bigl(\frac{r_3/r_1}{1 - s_j e^{2iu\theta_c}}\bigr) + i\pi r_0 j - 2\pi i c p_0(u + \frac{1}{2})$ (so $e^{\zeta(u_k)} = a_k^{r_0}$ identically for both parities of $c$, with $\operatorname{Im}\zeta(0) = \pi D$ and $\operatorname{Im}\zeta(\pm 1/2) = \pi D \pm \pi D$ with zero endpoint error). Setting $w(u) = s_j e^{2iu\theta_c}$, $A(u) = 2\pi i c((r_0 + p_0)w(u) - p_0)$, $B = 4\pi c r_0 \theta_c$, and using the exact uniform bound $\frac{\theta_c}{1 - s_j} \le \frac{\kappa}{2\pi S_{\text{half}}\cos\alpha_1}$, every factor of $(1 - w)$ cancels identically from the first three IBP kernels $W_{1,\text{fin}} = \frac{Bw}{A^2}$, $W_{2,\text{fin}}$, $W_{3,\text{fin}}$. Integrating the finite-$n$ continuum potential $\Psi_{\text{fin}}(\theta; u_0) := -\log\Lambda - i\theta + \int_{-L_n}^{L_n}\log(1 + \Lambda e^{\zeta(u)+i\theta})\,du$ ($\Lambda = e^{-\operatorname{Re}\zeta(u_0)}$) by parts 1, 2, and 3 times expresses $\Psi_{\text{fin}}(\theta) - \Psi_{\text{fin}}(0)$, $g_{1,\text{fin}}(u_0) := \Psi_{\text{fin}}'(0)/i$, and $\Psi_{\text{fin}}''(\theta)$ explicitly with exact harmonicity $\frac{d}{du_0} g_{1,\text{fin}}(u_0) = (-\operatorname{Re}\Psi_{\text{fin}}''(0; u_0))(-\operatorname{Re}\zeta'(u_0))$ and conjugation symmetry $\Psi_{\text{fin}}(-\theta) = \overline{\Psi_{\text{fin}}(\theta)}$.

3. **Exact Integrated Meromorphic Residue + Uniform-in-$r_0$ Shifted-Line Poisson Identity (Sub-Lemma G3.3):** Partitioning $|w| = \lambda$ into $r_0$ equal sectors $\phi = \psi_\nu + \theta/r_0$ ($\psi_\nu = \frac{(2\nu + j\bmod 2)\pi}{r_0}$, $\theta \in [-\pi, \pi]$) gives $\log F_{\text{pure}}(\lambda e^{i(\psi_\nu + \theta/r_0)}) - \tilde{c}\,\Psi_{\text{fin}}(\theta; u_0)$ with **zero additive constant** ($C_\nu \equiv 0$, verified at $\lambda = 1$). Shifting every Fourier harmonic $k \ne 0$ to height $\pm Y_{\text{big}} = \pm 0.50\, r_0 \ll y_{\text{branch}} \ge 6.96\, r_0$ certifies shifted-line error $\le 5.26\times 10^{-11}$, negative-binomial tail $\le 7.32\times 10^{-18}$, and $T(x_p) < 0.50$ ($c_{\text{re,50}} \ge 0.774$ via MVT $\operatorname{sinc}(2b)$ and concavity $(\log 1.5 - \log(1+T_{\text{in}}))/(0.5 - T_{\text{in}})$) for every pole inside $|\operatorname{Im} u| \le Y_{\text{big}}$. All 4 $m = 1$ boundary poles (2 inner in $y > 0$, 2 outer in $y < 0$), all $k \ge 2$ higher-pole branches (both $y > 0$ and $y < 0$, plus deep interior poles when $|D| \ge 2$), and the zero-hardcoded valley bump factor $V_{\text{bump}}$ (with $a_{\text{sep}} := \min(a_{\text{in}}, a_{\text{out}})$ accounting for $y < 0$ crossings) are explicitly bounded in 50-digit `iv`.

4. **Five Exhaustive Sub-Cases (Sub-Lemma G3.4):**

| Sub-Case | Parameter Slice inside Box G3 | Certified Method & Key Parameters | 50-Digit `iv` Margin (`I_full / L^1`) |
|---|---|---|---|
| **G3-A1** | $|D| \ge 2,\; s \le 2.0$ | $(1-w)$-free $|R_{3,\text{fin}}|/H_{\text{fin}0} \le 0.0995$; 2-pole split at $|\delta_{\text{fin}}| = 0.08$ (`A1a`: $V_{\text{sep}} \le 2.197$; `A1b`: $V \le 4.135$) | **`>= +0.08207`** (`A1a: +0.08207`, `A1b: +0.11879`) |
| **G3-A2** | $|D| \ge 2,\; s \ge 2.0$ (incl. closed-form tail $s \ge 10$) | $|\delta_*| \le \min(0.030, 0.40/s^2)$, drop $\ge 3.435 X$, $V_{\text{bump}} \le 4.105$ | **`>= +0.40597`** |
| **G3-B** | $|D| = 1,\; \tau \ge 40.0$ | Linear drift $r_{\text{drift}} \le 0.04280$, $|\delta_*| \le \min(0.020, 0.80/\tau)$, $V_{\text{bump}} \le 4.035$ | **`>= +0.13156`** |
| **G3-C** | $|D| = 1,\; \tau \in (0, 1.05]$ | $u_0 = \max(u_{0,\text{fin}}^*, u_{\text{floor}})$, horizontal saddle shift + certified $g_{1,\text{fin}}(0.42) > 0$, $D=-1$ residue, explicit $[\theta_{\text{use}}, \pi]$ buffer, $V_{\text{bump}} \le 4.007$ | **`>= +0.07706`** |
| **G3-D1** | $|D| = 1,\; \tau \in [1.05, 1.50]$ (3 $\tau$-slices) | Fixed-$u_0$ Tier A boxes + `certify_G3D_harmonic_tail` ($\tilde{c} \ge 105..135$) | **`>= +0.1646`** |
| **G3-D2** | $|D| = 1,\; \tau \in [1.50, 2.80]$ (7 $\tau$-slices) | Tier A + harmonicity tail (`[1.50, 1.65]`) / subdivided exact-saddle Tier B (`harmonic_u_star_bracket_iv` + `v0_valley_drop_floor_iv`) | **`>= +0.0545`** |
| **G3-D3** | $|D| = 1,\; \tau \in [2.80, 40.0]$ (27 $\tau$-slices) | Tier A + subdivided exact-saddle Tier B (`D3a1 [2.80,4.00]: +0.0784`, `D3a2 [4.00,5.68]: +0.0840`, `D3b1 [5.68,11.2]: +0.0952`, `D3b2 [11.2,40.0]: +0.2055`) | **`>= +0.0784`** |

> **Theorem 6.1 (Complete 50-Digit Certification of Piece G3 — Full Exact Integral; signed off `#124`, `#132`).**  
> For every strictly interior $(n, c, j)$ in Box G3, shifting the Theorem 4.2 contour to the Pigeonhole-Annulus radius $|w| = \tilde{\lambda}$ certifies for the **full exact cubic-root integral**:
> $$\frac{I_{\text{full}}(\tilde{\lambda})}{L^1(\tilde{\lambda})} \;\ge\; \mathbf{+0.0545} \;>\; 0,$$
> and therefore $(-1)^j Q_c^{(2m)}(1, u_j) > 0$. Independently audited line-by-line and executed by `@verifier` with exit code `0`: [`analyst/certify_G3_ABC_full_iv.py`](analyst/certify_G3_ABC_full_iv.py) (SHA-256 `7727BD251CE9AE60510FF2D7C24F546A1093CF299FC71E2C25F958E0BFDB3642`, 29.78s, `#124`) and all 6/6 seamless exhaustive CLI parts (37/37 $\tau$-slices, 420+ boxes) of [`analyst/certify_G3_D_full_iv.py`](analyst/certify_G3_D_full_iv.py) (SHA-256 `4479818E725B262DDAB836D4656C1DC266F71EA6F7873E6D96536B725F6DDD19`, 3054.59s on `--workers 3`, `#132`). All 32/32 extreme finite-$n$ Piece G3 corners (up to $n = 9.85\times 10^7$, $D = \pm 1, \pm 2, \pm 3$, $r_0 \in \{1, 2, 3, 5\}$) pass [`analyst/verify_G3_empirical_corners.py`](analyst/verify_G3_empirical_corners.py) with empirical $\min I_{\text{pure}}/L^1 = +0.54579$. $\square$

---

## 7. Proof Completion & Reproducibility Ledger

### 7.1 Assembly of Theorem 1.1 and Corollary 1.2

Fix any even integer $n = 2m \ge 2$ and any $(a, b, c) \in \mathbb{Z}_{\ge 0}^3$ with $a + 2b + 3c = 2m$.
- If $n \le 130$, Theorem 2.1 directly gives $c_{a,b,c}^{(2m)} \ne 0$ and $\operatorname{sign}(c_{a,b,c}^{(2m)}) = (-1)^{b(c+1)+c}$.
- If $n \ge 132$ and $c = 0$, $a = 0$, $b = 0$, or $c = 1$ (which includes all $c > \lfloor (n-3)/3 \rfloor$), Theorems 3.3–3.5 give the exact sign law (and simple real-rootedness of the $a = 0$ and $b = 0$ slices).
- If $n \ge 132$ and $2 \le c \le \lfloor (n - 3)/3 \rfloor$, every strictly interior sample index $j \in \mathcal{J}_{n, c}$ falls into either Regime 1 ($\kappa \ge 2.5$, Theorem 4.6) or one of the five overlapping Regime 2 sub-regimes established by the Zero-Crevice Partition (Section 1.3): Sub-regime 2A (Theorem 5.2), Sub-regime 2C-Upgrade (Theorem 5.3), Piece G2 (Theorem 5.5), Piece G1 (Theorem 5.6), or Piece G3 (Theorem 6.1). In every case, $(-1)^j Q_c^{(2m)}(1, u_j) > 0$. By Theorem 4.1, $Q_c^{(2m)}(1, u)$ has $D_c$ simple real roots of the stated sign, and every coefficient $c_{a,b,c}^{(2m)}$ is nonzero with exact sign $(-1)^{b(c+1)+c}$.

This completes the proof of Theorem 1.1 and Corollary 1.2 for all even integers $n = 2m \ge 2$. $\square$

### 7.2 Verification Script Manifest

| Proof Component | Script / Document | Verification Scope & Runtime | Audit Status |
|---|---|---|---|
| Finite Base ($n = 2..84$) | [`verify/engine.py`](verify/engine.py) | 9,464 exact integer coefficients; SymPy & Wendt cross-checks | **SIGNED OFF** (`#39`) |
| Finite Base ($n = 86..130$) | [`verify/fast_certify_140.py`](verify/fast_certify_140.py) | 23,953 exact integer coefficients; exact agreement with `engine.py` on $n=2..30$ | **SIGNED OFF** (`#57`) |
| Resultants & Boundaries | [`algebra/probe_depressed_cubic.py`](algebra/probe_depressed_cubic.py), [`algebra/verify_c1_and_probe_c2.py`](algebra/verify_c1_and_probe_c2.py) | Exact algebraic verification of Theorems 3.1–3.5 | **SIGNED OFF** (`#42`, `#45`) |
| 1D Circle Integral & $B_{\text{avg}}$ | [`algebra/verify_exact_identities_mp.py`](algebra/verify_exact_identities_mp.py) | 80-digit `mpmath` verification of Theorem 4.2 & Lemma 4.5 ($< 10^{-50}$ err) | **SIGNED OFF** (`#59`) |
| Regime 1 ($\kappa \ge 2.5, n \ge 84$) | [`algebra/verify_regime1_complete.py`](algebra/verify_regime1_complete.py) | 6 closed-form margins + 13,006 finite-$n$ cross-checks (0 failures) | **SIGNED OFF** (`#59`, `#61`) |
| R2.1 & Sub-regime 2A | [`analyst/certify_R2_1_and_2A_iv.py`](analyst/certify_R2_1_and_2A_iv.py) | 50-digit `mpmath.iv` certification of Lemma 5.1 & Theorem 5.2 (`+0.220613`) | **SIGNED OFF** (`#74`) |
| Sub-regime 2C-Upgrade | [`analyst/certify_R2_3_iv.py`](analyst/certify_R2_3_iv.py) | 50-digit `mpmath.iv` certification of Theorem 5.3 (`+0.441480` / `+0.266331`) | **SIGNED OFF** (`#77`) |
| Piece G2 ($D = 0$ Major Arcs) | [`algebra/certify_G2_D0.py`](algebra/certify_G2_D0.py) | 50-digit `mpmath.iv` certification of Theorem 5.5 (`+0.00176348` / `+0.216651`, widened $q=1$ to all $\kappa<2.5$) + 28 parity corners; 117.7s, exit 0 | **SIGNED OFF** (`#90`, `#103`) |
| Piece G1 (Dirichlet Minor Arcs) | [`algebra/certify_G1_minor.py`](algebra/certify_G1_minor.py), [`algebra/G1_NOTES.md`](algebra/G1_NOTES.md) | 50-digit `mpmath.iv` certification of Theorem 5.6 (`+0.03355852`, $|\Phi^{\text{ex}}| \le 1.536448$) + 27/27 finite-$n$ corners; exit 0 | **SIGNED OFF** (`#101`) |
| Piece G3 ($D \ne 0$ Major Arcs) | [`analyst/certify_G3_ABC_full_iv.py`](analyst/certify_G3_ABC_full_iv.py), [`analyst/certify_G3_D_full_iv.py`](analyst/certify_G3_D_full_iv.py), [`analyst/verify_G3_empirical_corners.py`](analyst/verify_G3_empirical_corners.py), [`analyst/verify_regime2_analyst_iv.py`](analyst/verify_regime2_analyst_iv.py) | 50-digit `mpmath.iv` certification of Theorem 6.1 (`>= +0.0545` global; A1 `+0.08207`, A2 `+0.40597`, B `+0.13156`, C `+0.07706`, D1 `+0.1646`, D2 `+0.0545`, D3 `+0.0784`; 37/37 $\tau$-slices, 32/32 empirical corners); `ABC` 29.78s (`7727BD25...`), `D` 3054.59s (`4479818E...`), exit 0 | **SIGNED OFF** (`#124`, `#132`) |

---

## References

1. V. G. Bardakov and T. A. Kozlovskaya, *Multivalued groups and Newton polyhedron*, [arXiv:2305.07261](https://arxiv.org/abs/2305.07261) (2023); *Sibirskie Elektronnye Matematicheskie Izvestiya* **20**:2 (2023), 1590–1596.
2. I. Gaiur, V. Rubtsov, and D. van Straten, *Higher Bessel functions and $n$-valued groups*, [arXiv:2405.03015](https://arxiv.org/abs/2405.03015) (2024).
3. V. M. Buchstaber and M. Kornev, *Multivalued groups and Wendt determinants*, [arXiv:2505.04296](https://arxiv.org/abs/2505.04296) (2025); *Sbornik: Mathematics* **217**:8 (2026).
