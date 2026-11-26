# Bardakov–Kozlovskaya Conjecture 3.3(2): Sources, Exact Representations & Research State

**Prepared by:** `scout` (2026-11-26)  
**Tracking issue:** `1337isnot1337/solved-conjectures` [Issue #174](https://github.com/1337isnot1337/solved-conjectures/issues/174) · [MathDB #358383](https://mathdb.com/p/358383/conjecture-on-coefficients-of-multivalued-group-polynomials)  
**Openness audit (2026-11-26):** **Part (2) is 100% GENUINELY OPEN** (verified live across Bardakov–Kozlovskaya 2023, Gaiur–Rubtsov–van Straten 2024, Buchstaber–Kornev 2025/2026, `solved-conjectures#174`, MathDB `#358383`, `ProgramComputer`, `PrincipiaMathOrg`, arXiv, and web). Part (1) (prime-power divisibility) is already resolved (corollary of Buchstaber–Kornev Thm. 4/11, and proved directly via resultant Frobenius in `origin/research`). **Solving Part (2) completely closes Conjecture 3.3.**

---

## 1. Exact Statement, Indexing & Sign Conventions

### Primary Source
- **V. G. Bardakov and T. A. Kozlovskaya**, *Multivalued groups and Newton polyhedron*, [arXiv:2305.07261v1](https://arxiv.org/abs/2305.07261) (May 2023, Conjecture 3.3, p. 5); journal version: *Sibirskie Elektronnye Matematicheskie Izvestiya* **20**:2 (2023), 1590–1596, DOI `10.33048/semi.2023.20.097` (unnumbered conjecture, p. 1594).

### Definition of $p_n(x, y, z)$
On $\mathbb{C}$, the $n$-valued Buchstaber–Novikov group $\mathbb{G}_n$ has multiplication $x * y = [(\sqrt[n]{x} + \zeta^k \sqrt[n]{y})^n : 0 \le k < n]$, unit $0$, and inverse $\mathrm{inv}(x) = (-1)^n x$ ($\zeta = e^{2\pi i/n}$). Its defining symmetric homogeneous integer polynomial of degree $n$ is:
$$p_n(x, y, z) := \prod_{k=0}^{n-1} \Bigl(z - (\mathrm{inv}(x) * \mathrm{inv}(y))_k\Bigr) = \prod_{k=0}^{n-1} \Bigl(z - (-1)^n (\alpha + \zeta^k \beta)^n\Bigr) \qquad (\alpha^n = x,\; \beta^n = y).$$
- **Sign convention for even $n = 2m$:** $(-1)^n = +1$, so $\mathrm{inv}(x) = x$ and
  $$p_{2m}(x, y, z) = \prod_{k=0}^{2m-1} \Bigl(z - (\alpha + \zeta^k \beta)^{2m}\Bigr).$$
- **Symmetric double-product identity (valid for EVERY $n \ge 1$ with $\alpha^n=x, \beta^n=y, \gamma^n=z$, zero sign ambiguity!):**
  Since $\prod_{r=0}^{n-1}(\gamma + \zeta^r u) = \gamma^n - (-u)^n = z - (-1)^n u^n$ for all $n \ge 1$,
  $$p_n(\alpha^n, \beta^n, \gamma^n) = \prod_{r=0}^{n-1} \prod_{s=0}^{n-1} \bigl(\gamma + \zeta^r \alpha + \zeta^s \beta\bigr).$$
  *(Symmetry in $x,y,z$ is manifest!)*

### Elementary Symmetric Expansion & Exact Indexing Bijection
Write $e_1 = x+y+z$, $e_2 = xy+yz+zx$, $e_3 = xyz$. There are two equivalent indexings in the literature:
$$p_n(x,y,z) = \sum_{\substack{k_1 \ge k_2 \ge k_3 \ge 0 \\ k_1 + k_2 + k_3 = n}} A_{k_1, k_2, k_3}\, e_1^{k_1 - k_2} e_2^{k_2 - k_3} e_3^{k_3} = \sum_{\substack{a,b,c \ge 0 \\ a + 2b + 3c = n}} c_{a,b,c}^{(n)}\, e_1^a e_2^b e_3^c,$$
related bijectively by
$$a = k_1 - k_2,\quad b = k_2 - k_3,\quad c = k_3 \quad \Longleftrightarrow \quad k_1 = a+b+c,\quad k_2 = b+c,\quad k_3 = c, \qquad A_{k_1,k_2,k_3} = c_{a,b,c}^{(n)}.$$

### Conjecture 3.3 (Bardakov–Kozlovskaya, 2023)
1. *(Already proved)* If $n = p^m$ is a prime power, every coefficient except $c_{n,0,0}^{(n)} = 1$ is divisible by $p$.
2. **(OPEN TARGET):** **If $n$ is even, then every admissible coefficient $c_{a,b,c}^{(n)}$ ($a,b,c \ge 0$, $a + 2b + 3c = n$) is nonzero.**

### Small-$n$ Reference Values (from Bardakov–Kozlovskaya & Buchstaber–Kornev Ex. 8)
- $p_1 = e_1$, $p_2 = e_1^2 - 4 e_2$, $p_3 = e_1^3 - 27 e_3$,
- $p_4 = e_1^4 - 8 e_1^2 e_2 + 16 e_2^2 - 128 e_1 e_3 = p_2^2 - 128 e_1 e_3$,
- $p_5 = e_1^5 - 5^4 e_1^2 e_3 + 5^5 e_2 e_3$,
- $p_6 = e_1^6 - 12 e_1^4 e_2 + 48 e_1^2 e_2^2 - 64 e_2^3 - 2754 e_1^3 e_3 - 12312 e_1 e_2 e_3 + 185193 e_3^2$ (note $2754 = 2\cdot 3^4\cdot 17$, $12312 = 2^3\cdot 3^4\cdot 19$, $185193 = 3^3\cdot 19^3$),
- $p_{18}$: all 37 admissible coefficients with full prime factorizations appear in Buchstaber–Kornev (`arXiv:2505.04296v2`, Example 8) — e.g., $c_{0,0,6}^{(18)} = +3^9 (109\cdot 163\cdot 271\cdot 2269\cdot 5779)^3$ (a perfect cube!).

---

## 2. All Known Exact Representations of $p_n(x, y, z)$

1. **Double-product & paired-root product (even $n = 2m$):**
   $$p_{2m}(\alpha^{2m}, \beta^{2m}, \gamma^{2m}) = \prod_{r,s=0}^{2m-1} (\gamma + \zeta^r \alpha + \zeta^s \beta) = \prod_{k=0}^{m-1} \Bigl[ z^2 - z\bigl((\alpha+\zeta^k\beta)^{2m} + (\alpha-\zeta^k\beta)^{2m}\bigr) + (\alpha^2 - \zeta^{2k}\beta^2)^{2m} \Bigr].$$
   *(At $e_1 = e_2 = 0$, $(x,y,z) = (1,\omega,\omega^2) e_3^{1/3}$ with $\omega = e^{2\pi i/3}$, giving $c_{0,0,n/3}^{(n)} = p_n(1,\omega,\omega^2)$ whenever $3 \mid n$; for $6 \mid n$, note $p_n(1,1,1) = 0$ because $(1+e^{\pm 2\pi i/3})^6 = 1$.)*

2. **Exact integer resultant (`origin/research:research/2026-11-21-combinatorics-bardakov/proof.md`, eq. (2)):**
   $$p_n(x, y, t^n) = \operatorname{Res}_u\bigl(u^n - x,\; (t+u)^n - (-1)^n y\bigr).$$
   Also Buchstaber–Kornev Prop. 10: $p_n(z; x, y) = (-1)^n \operatorname{Res}_t\bigl(z - (u+vt)^n,\; t^n - 1\bigr)$ ($u^n = (-1)^n x$, $v^n = (-1)^n y$).

3. **Kronecker sum of companion matrices & $n \times n$ Wendt $(x,y,z)$-matrix determinant (Buchstaber–Kornev, `arXiv:2505.04296v2` Prop. 7 & Thm. 4; *Sbornik: Math.* 217:8 (2026), Thm. 11):**
   - Kronecker sum ($w^n = z$):
     $$p_n(w^n; x, y) = \det\!\Bigl(w I_{n^2} - F\bigl(0,\dots,0,(-1)^{n+1}x\bigr) \boxplus F\bigl(0,\dots,0,(-1)^{n+1}y\bigr)\Bigr).$$
   - $n \times n$ $y$-circulant Wendt determinant ($w^n = z$; sign alternation on odd diagonals does not change $\det$):
     $$p_n(z; x, y) = \det \mathrm{Circ}_y\!\left(w^n + (-1)^{n+1} x + y,\; \binom{n}{1} w,\; \binom{n}{2} w^2,\; \dots,\; \binom{n}{n-1} w^{n-1}\right).$$
     For **even $n$**, entry $(i,j)$ ($1 \le i,j \le n$) is: $z - x + y$ if $i=j$; $y \binom{n}{j-i} w^{j-i}$ if $j > i$; $\binom{n}{n-(i-j)} w^{n-(i-j)}$ if $i > j$.

4. **Univariate discriminant & Higher Bessel singular locus (Gaiur–Rubtsov–van Straten, `arXiv:2405.03015` Thms. 2.3–2.4; Buchstaber–Kornev Thm. 8):**
   For every $n \ge 2$, define the degree-$(2n-2)$ polynomial in $T$:
   $$P_{x,y,z}(T) := \bigl((-1)^n x T^{n-1} + (-1)^n y\bigr)(1+T)^{n-1} - T^{n-1} z.$$
   Then
   $$(-1)^n (n-1)^{2n-2} (xyz)^{n-2}\, p_n(x, y, z) = \operatorname{disc}_T\!\bigl(P_{x,y,z}(T)\bigr).$$
   - $p_n(x,y,z) = 0$ is the singular locus of the $n$-Bessel multiplication kernel $K_n(x,y,z) = \sum_{i,k \ge 0} \binom{i+k}{k}^n x^i y^k z^{-(i+k)}$ (the period of a $(n-1)$-dimensional hypersurface family).
   - **GRS24 Theorem 2.4:** The projective curve $p_n(x,y,z) = 0$ is a **rational curve** (for even $n$, parametrized by $[u^n : v^n : (u+v)^n]$) with $(n-1)(n-2)/2$ double points, ALL lying in $\mathbb{P}^2(\mathbb{R})$ with coordinates in $\mathbb{Q}(\cos(2\pi/n))$.

5. **Newton power-sum identity (`origin/research` `proof.md` eq. (6)):**
   Viewing $p_n(x,y,z) = \prod_{k=0}^{n-1}(z - r_k)$ in $z$, the $j$-th power sum of roots $r_k = (-1)^n (\alpha + \zeta^k\beta)^n$ is
   $$S_j(x, y) := \sum_{k=0}^{n-1} r_k^j = n (-1)^{nj} \sum_{a=0}^j \binom{nj}{na} x^{j-a} y^a.$$

---

## 3. Already-Proved Results & Structural Clues

1. **Exact $c = 0$ boundary (Bardakov–Kozlovskaya Prop. 3.1):**
   - Odd $n$: $p_n(x, y, 0) = (x+y)^n = \bar{e}_1^n \implies c_{a,b,0}^{(n)} = 0$ for all $b > 0$.
   - Even $n = 2m$: $p_{2m}(x, y, 0) = (x-y)^{2m} = (\bar{e}_1^2 - 4\bar{e}_2)^m$, hence
     $$c_{2(m-b), b, 0}^{(2m)} = (-4)^b \binom{m}{b} \ne 0 \qquad (0 \le b \le m).$$
2. **Frobenius reduction modulo any prime $p$ (`origin/research` `proof.md` Thm., proved via resultant):**
   For every prime $p$, $q = p^m$, and $s \ge 1$ (no coprimality needed):
   $$p_{qs}(x, y, z) \equiv p_s(x, y, z)^q \pmod{p}, \qquad \text{hence } P_{qs}(E_1, E_2, E_3) \equiv P_s(E_1, E_2, E_3)^q \pmod{p}.$$
   - Taking $s=1$ proves Part (1): $P_{p^m}(E_1, E_2, E_3) \equiv E_1^{p^m} \pmod{p}$.
   - Taking $p = 2, q = 2^r$ gives $P_{2^r s} \equiv P_s^{2^r} \pmod 2$.
3. **Prime Wolstenholme divisibility (Buchstaber–Kornev Prop. 8 / Thm. 12):** For prime $p \ge 5$, $p_p(x,y,z) - e_1^p \in p^4 e_3\,\mathbb{Z}[x,y,z]$.
4. **Irreducibility (Buchstaber–Kornev Thm. 6):** $p_n(x,y,z)$ is irreducible over $\mathbb{C}$ (for all $n \ge 1$).
5. **MAJOR UNPROVED CLUE — Observed Exact Sign Law (`origin/research` `proof.md` eq. (7) & Issue #174 Comment #2):**
   Checked exactly for ALL even $n \le 40$ (and all 37 terms of $p_{18}$ in Buchstaber–Kornev Ex. 8): every admissible coefficient is nonzero and satisfies
   $$\operatorname{sign}\!\left(c_{a,b,c}^{(2m)}\right) = (-1)^{b(c+1)+c} = (-1)^{b(a+1)+c} = (-1)^{k_1 k_2 + k_1 k_3 + k_2 k_3} = \begin{cases} (-1)^b & \text{if } c \text{ is even,} \\ -1 & \text{if } c \text{ is odd!} \end{cases}$$
   *(Equivalently, splitting $P_{2m}(e_1, e_2, e_3) = P_{2m}^{\mathrm{even}}(e_1, e_2, e_3) + P_{2m}^{\mathrm{odd}}(e_1, e_2, e_3)$ by parity of $e_3$-exponent: **$P_{2m}^{\mathrm{even}}(e_1, -e_2, e_3)$** and **$-P_{2m}^{\mathrm{odd}}(e_1, e_2, e_3)$** both have ALL STRICTLY POSITIVE coefficients through $n \le 40$!)*
   - For odd $n$ and $c > 0$, `proof.md` also notes observed sign $(-1)^{c(b+1)}$ (no proof claimed).

---

## 4. Remaining Gap & Key Questions for `@algebraist`, `@analyst`, `@verifier`

- **Exact remaining gap:** Prove $c_{a,b,c}^{(2m)} \ne 0$ for EVERY even $n = 2m \ge 2$ and every $a,b,c \ge 0$ with $a + 2b + 3c = 2m$ — or find an exact counterexample!
- **Immediate decisive test for `@verifier`:** Push exact computation beyond $n = 40$ (e.g., $n = 42..80+$) to test whether the sign law $(-1)^{b(c+1)+c}$ holds universally or ever flips! Notice:
  - If the sign law **holds universally**, proving the strict sign law (positive coefficients of $P_{2m}^{\mathrm{even}}(e_1, -e_2, e_3)$ and $-P_{2m}^{\mathrm{odd}}(e_1, e_2, e_3)$) immediately proves Conjecture 3.3(2)!
  - Watch out for a sign flip at larger $n$ — if a coefficient changes sign across $n$, nearby $(a,b,c)$ might either stay nonzero or hit an exact zero (counterexample!).
- **High-leverage handles for `@algebraist` & `@analyst`:**
  1. Look at $P_{2m}(e_1, -e_2, e_3)$ vs $P_{2m}(e_1, e_2, -e_3)$! Substituting $(e_1, e_2, e_3) \mapsto (e_1, -e_2, \pm e_3)$ maps the Vieta cubic $T^3 - e_1 T^2 + e_2 T - e_3 = 0$ to roots $(x,y,z)$ with concrete real/sign structure (e.g., for $e_1, e_2, e_3 > 0$, $T^3 - e_1 T^2 - e_2 T \mp e_3 = 0$ has real/conjugate roots with definite signs!).
  2. Notice how $P_{2m}(e_1, e_2, e_3)$ is BOTH a $y$-circulant Wendt determinant (`R3`), a univariate discriminant (`R4`), AND the symmetric product $\prod_{r,s}(\gamma + \zeta^r\alpha + \zeta^s\beta)$ (`R1`)! For example, under $c \mapsto -c$ (i.e. $e_3 \mapsto -e_3$), the roots $(x,y,z)$ are swapped with $(-x,-y,-z)$? Wait: $(e_1, e_2, -e_3)$ does NOT come from $(-x,-y,-z)$ (that gives $(-e_1, e_2, -e_3)$!), whereas $(e_1, e_2, e_3) \mapsto (-e_1, e_2, -e_3)$ leaves $p_{2m}$ invariant because $p_{2m}$ is homogeneous of even degree $2m$! That's why $a + c$ is ALWAYS EVEN ($a \equiv c \pmod 2$), so $P_{2m}^{\mathrm{even}}$ has $(a, c)$ both even and $P_{2m}^{\mathrm{odd}}$ has $(a, c)$ both odd!
  3. **CHECK THAT PARITY FACT!** Since $a + 2b + 3c = 2m$, **$a$ and $c$ always have the SAME parity**! Therefore:
     $$P_{2m}^{\mathrm{even}}(e_1, e_2, e_3) = \frac{1}{2}\Bigl(P_{2m}(e_1, e_2, e_3) + P_{2m}(-e_1, e_2, e_3)\Bigr),$$
     $$P_{2m}^{\mathrm{odd}}(e_1, e_2, e_3) = \frac{1}{2}\Bigl(P_{2m}(e_1, e_2, e_3) - P_{2m}(-e_1, e_2, e_3)\Bigr)!$$
     (Because negating $e_1$ alone multiplies $e_1^a e_2^b e_3^c$ by $(-1)^a = (-1)^c$!)
     Wait: ALSO $P_{2m}(e_1, -e_2, e_3)$ vs $P_{2m}(e_1, e_2, e_3)$: check closed forms for small $c$ ($c = 1, 2, 3$) and small $a$ or $b$, AND check the rational-curve / Cauchy-integral / Schur-function expansion of $p_{2m}$!
