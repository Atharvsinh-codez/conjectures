# Gessel–Lin–Zeng Conjecture 15: Sources, Conventions, Known Results & Exact Gaps

**Prepared by:** `scout` (2026-11-26)  
**Tracking issue:** `1337isnot1337/solved-conjectures` [Issue #84](https://github.com/1337isnot1337/solved-conjectures/issues/84)  
**Status:** **100% GENUINELY OPEN** for $4 \le i \le k-3$ (SymCat live open problem; zero external competition in `ProgramComputer`, `PrincipiaMathOrg`, open PRs, or arXiv as of 2026-11-26).

---

## 1. Exact Statement, Definitions & Conventions

### Primary Source
- **I. M. Gessel, Z. Lin, and J. Zeng**, *Jacobi–Stirling polynomials and $P$-partitions*, European J. Combin. **33**(8) (2012), 1987–2000 ([DOI:10.1016/j.ejc.2012.06.008](https://doi.org/10.1016/j.ejc.2012.06.008), [arXiv:1201.0622](https://arxiv.org/abs/1201.0622)), Proposition 6 & Conjecture 15 (numbered as Problem 2 in arXiv v1).
- **SymCat live record:** [`realRootedOpenProblems.htm`](https://www.symmetricfunctions.com/realRootedOpenProblems.htm) & [`realRootedSetPartitions.htm`](https://www.symmetricfunctions.com/realRootedSetPartitions.htm).

### Combinatorial & Generating-Function Definitions
1. **Jacobi–Stirling numbers of the second kind** $\mathrm{JS}(n, k; z)$:
   $$x^n = \sum_{k=0}^n \mathrm{JS}(n,k;z) \prod_{j=0}^{k-1} \bigl(x - j(z+j)\bigr), \qquad \mathrm{JS}(n,k;z) = \mathrm{JS}(n-1,k-1;z) + k(k+z)\mathrm{JS}(n-1,k;z).$$
   Write the $k$-th diagonal as $\mathrm{JS}(n+k, n; z) = \sum_{i=0}^k p_{k,i}(n) z^i$. Its Ehrhart-type numerator is
   $$(1-x)^{3k-i+1} \sum_{n \ge 0} p_{k,i}(n) x^n =: \tilde{A}_{k,i}(x) = x\,A_{k,i}(x) \qquad (k \ge 1).$$
2. **Jacobi–Stirling permutations & internal descents (standard Ma–Wang / SymCat / `origin/research` convention):**
   For $k \ge 1$, let $M_k = \{\bar{1}, 1, 1, \bar{2}, 2, 2, \dots, \bar{k}, k, k\}$ ordered by $\bar{1} < 1 < \bar{2} < 2 < \cdots < \bar{k} < k$. For $S \subseteq [k]$, delete the barred letters $\{\bar{j} : j \in S\}$; let $\mathrm{JSP}_{k,S}$ be the set of permutations (words) of the remaining multiset (length $3k - |S|$) such that every entry strictly between the two unbarred copies of $j$ is strictly larger than $j$. Define
   $$A_{k,S}(x) := \sum_{w \in \mathrm{JSP}_{k,S}} x^{\mathrm{des}(w)}, \qquad A_{k,i}(x) := \sum_{\substack{S \subseteq [k] \\ |S| = i}} A_{k,S}(x) \qquad (0 \le i \le k),$$
   where $\mathrm{des}(w) = \#\{r : w_r > w_{r+1}\}$ counts **internal descents** only (no final descent).
   - **Convention difference vs GLZ (2012):** GLZ's original $\tilde{A}_{k,i}(x) = x\,A_{k,i}(x)$ includes a final descent (adding a trivial zero at $x=0$). Removing this factor $x$ does not affect real-rootedness.
   - **Boundary & degree facts (all $k \ge 1$, $0 \le i \le k$):**
     - $A_{0,0}(x) := 1$; $A_{k,i} := 0$ outside $0 \le i \le k$.
     - $\deg A_{k,S}(x) = \deg A_{k,i}(x) = 2k - i - 1$ (attained by the nonincreasing word).
     - $A_{k,S}(0) = 1$ (unique nondecreasing word), so $A_{k,i}(0) = \binom{k}{i} > 0$, and all coefficients of $A_{k,i}(x)$ are strictly positive integers (GLZ Cor. 7).
     - $A_{k,k}(x)$ is the second-order Eulerian polynomial (without the trailing factor $x$, OEIS A008517); $A_{k,0}(x)$ is the full Legendre/Jacobi–Stirling descent polynomial.
   - **Small rows:**
     - $k=1$: $A_{1,0} = 1+x$, $A_{1,1} = 1$.
     - $k=2$: $A_{2,0} = 1 + 14x + 21x^2 + 4x^3$, $A_{2,1} = 2 + 12x + 6x^2$, $A_{2,2} = 1 + 2x$.
     - $k=3$: $A_{3,0} = 1 + 75x + 603x^2 + 1065x^3 + 460x^4 + 36x^5$, $A_{3,1} = 3 + 114x + 501x^2 + 436x^3 + 66x^4$, $A_{3,2} = 3 + 55x + 116x^2 + 36x^3$, $A_{3,3} = 1 + 8x + 6x^2$.

### Conjecture 15 (Gessel–Lin–Zeng, 2012)
For every $k \ge 1$ and $0 \le i \le k$, **$A_{k,i}(x)$ has only real zeros** (equivalently, simple negative zeros — or simple/real nonpositive zeros; known cases all have simple negative zeros).
- **Known cases:**
  - $i = 0$ and $i = k$ (single summand $A_{k,S}$; real-rooted by Brenti 1989 Thm. 6.6.3 / GLZ 2012).
  - $i \in \{1, 2, 3, k-2, k-1\}$: proved by Ma–Wang (Oct 2026, `arXiv:2610.03111`).
- **Open range:** $4 \le i \le k-3$ (first open case is $k=7, i=4$).

---

## 2. Notation Clash Warning & Proved Operator Identities

> **CRITICAL NOTATION WARNING:**
> - The $x$-variable insertion operator $f \mapsto (1+Lx)f + x(1-x)f'$ is denoted **`S_L`** in `research/gessel-lin-zeng/README.md`, **`T_L`** in `obstruction.md`, and **`D_L`** in both Ma–Wang (`arXiv:2610.03111`) and `COMMON-INTERLACING-BARRIER.md`.
> - Meanwhile, `gessel-lin-zeng-followup.md` uses **`D_j`** for the **$t$-variable** operator $t(1+t)\frac{d}{dt} + 1 + jt$!
> - **Recommended swarm convention:** Always write **`S_L`** for $(1+Lx) + x(1-x)\partial_x$ (in $x$) and **`D_j`** for $t(1+t)\partial_t + 1 + jt$ (in $t$).

### Exact Recurrences (Proved for all $k, i, S$)
1. **Individual deletion-set insertion recurrence** (Ma–Wang eqs. (8)–(10)):
   Let $S_L f(x) := (1 + Lx)f(x) + x(1-x)f'(x)$, where $L$ is word length before insertion. With $L = 3(k-1) - |S|$:
   $$A_{k, S \cup \{\bar{k}\}} = S_L A_{k-1, S}, \qquad A_{k, S} = S_{L+1} S_L A_{k-1, S} \quad (\bar{k} \notin S),$$
   starting from $A_{0,\varnothing} = 1$ (note $S_0(1) = 1$, $S_1(1) = 1+x$).
2. **Aggregate $x$-recurrence** (equivalent to GLZ Prop. 6):
   $$A_{k,i}(x) = S_{3k-2-i}\!\Bigl(S_{3k-3-i} A_{k-1,i}(x) + A_{k-1,i-1}(x)\Bigr) \qquad (k \ge 1,\; 0 \le i \le k),$$
   with $A_{0,0} = 1$. Equivalently, for $a_{k,i,j} := [x^{j-1}] A_{k,i}(x)$ (GLZ Prop. 6 indexing, $1 \le j \le 2k-i$):
   $$a_{k,i,j} = j^2 a_{k-1,i,j} + [2(j-1)(3k-i-j-1) + (3k-i-2)]a_{k-1,i,j-1} + (3k-i-j)^2 a_{k-1,i,j-2} + j a_{k-1,i-1,j} + (3k-i-j)a_{k-1,i-1,j-1}.$$
3. **Möbius-transformed $i$-independent recurrence** (`gessel-lin-zeng-followup.md`, verified):
   Set $x = \frac{t}{1+t}$ (mapping $(-\infty, 0) \leftrightarrow (-1, 0)$ strictly increasingly), $D_j := t(1+t)\frac{d}{dt} + 1 + jt$, and
   $$h_{k,i}(t) := (1+t)^{2k-i} A_{k,i}\!\left(\frac{t}{1+t}\right).$$
   Since $D_j\bigl[(1+t)^M f(t/(1+t))\bigr] = (1+t)^{M+1} (S_{M+j-1} f)(t/(1+t))$, the $i$-dependence disappears completely:
   $$h_{k,i}(t) = D_k\bigl(D_k h_{k-1,i}(t) + h_{k-1,i-1}(t)\bigr), \qquad h_{0,0}(t) = 1,$$
   so the bivariate GF is the ordered operator product
   $$H_k(t, y) := \sum_{i=0}^k h_{k,i}(t) y^i = \prod_{j=1}^k \Bigl(D_j (D_j + y)\Bigr)(1).$$
   - **Removing the trivial $(1+t)$ factor:** Since $D_j \circ (1+t) = (1+t) \circ D_{j+1}$, setting $\tilde{h}_{k,i}(t) := h_{k,i}(t)/(1+t) = (1+t)^{2k-i-1} A_{k,i}\!\bigl(\frac{t}{1+t}\bigr)$ (degree $2k-i-1$, all roots in $(-1,0)$ iff $A_{k,i}$ has all roots in $(-\infty,0)$) gives for $k \ge 2$:
     $$\tilde{h}_{k,i} = D_{k+1}\bigl(D_{k+1} \tilde{h}_{k-1,i} + \tilde{h}_{k-1,i-1}\bigr), \qquad \tilde{h}_{1,0} = 1+2t,\;\; \tilde{h}_{1,1} = 1.$$
4. **Ma–Wang Proposition 1 (boundary differential relation):**
   $$A_{k, k-1}(x) = k A_{k+1, k+1}(x) - \frac{x}{2} A'_{k+1, k+1}(x), \qquad \text{via } S_{2k-1}(2kf - xf') = 2k S_{2k}f - x(S_{2k}f)'.$$

---

## 3. Summary of Ma–Wang (`arXiv:2610.03111`, Oct 2026)

Write $f \prec g$ (strict interlacing from the right) if $f, g \in \mathrm{RZ}(-\infty,0)$ have positive leading coefficients, $\deg f \le \deg g \le \deg f + 1$, and their roots $r_1 > r_2 > \cdots$, $s_1 > s_2 > \cdots$ satisfy $0 > s_1 > r_1 > s_2 > r_2 > \cdots$ (with $c \prec F$ for constant $c > 0$ and linear $F$).

- **Lemma 2 (single insertion):** If $f \in \mathrm{RZ}(-\infty,0)$ has positive leading coeff, degree $d$, and $L > d$, then $\deg(S_L f) = d+1$, positive leading coeff, and $f \prec S_L f$.
- **Lemma 3 (paired insertion preserving interlacing):** If $f, g \in \mathrm{RZ}(-\infty,0)$ have positive leading coeffs, $\deg g = \deg f + 1$, and $f \prec g$, then:
  1. For $L > \deg f$: $S_L f \prec S_{L+1} g$ (using $S_L f = \frac{f}{g} S_{L+1} g + xg \bigl((1-x)\frac{f}{g}\bigr)'$ and $((1-x)f/g)' < 0$).
  2. For $L > \deg g$: $S_L f \prec S_L g$ (using $S_L f = \frac{f}{g} S_L g + x(1-x)g(f/g)'$ and $(f/g)' < 0$).
- **Theorem 2 & Corollary 2 (five boundary families):** Let $H_k := S_{3k-3} A_{k-1,0}$ (intermediate word before inserting $kk$). For EVERY nonzero nonnegative linear combination $F_{k,i} = \sum_{|S|=i} c_S A_{k,S}$, Ma–Wang prove a **common interlacer** via Lemma 3 (exploiting that the paired word-length differences are $0$ or $1$):
  - $i = 1$: $A_{k,S} \prec A_{k,0}$ $\implies$ $F_{k,1} \prec A_{k,0}$.
  - $i = 2$: $A_{k,S} \prec H_k$ $\implies$ $F_{k,2} \prec H_k$.
  - $i = 3$: $A_{k,S} \prec A_{k-1,0}$ $\implies$ $F_{k,3} \prec A_{k-1,0}$.
  - $i = k-1$: $A_{k,k} \prec A_{k,S}$ $\implies$ $A_{k,k} \prec F_{k,k-1}$.
  - $i = k-2$: $A_{k+1,k+1} \prec A_{k,S}$ $\implies$ $A_{k+1,k+1} \prec F_{k,k-2}$.
- **Remark 1 (first open case $i=4$):** $A_{k,4} = S_{3k-6}(S_{3k-7}A_{k-1,4} + A_{k-1,3})$; since $A_{k-1,3} \prec A_{k-2,0}$, proving $S_{3k-7} A_{k-1,4} \prec A_{k-2,0}$ would suffice for $i=4$ (not proved in Ma–Wang).
- **Theorems 3 & 5 (Section 5):** Every $A_{k,S}$ and every $A_{k,i}$ ($k \ge 2$, $d = 2k-i-1$) is strictly top-heavy ($a_j < a_{d-j}$ for $0 \le j < d/2$) with strictly increasing left half ($a_0 < a_1 < \cdots < a_{\lfloor d/2 \rfloor}$).

---

## 4. Exact Obstructions Already Proved in `origin/research` (DO NOT REPEAT!)

1. **Per-subset common interlacing FAILS for ALL $i \ge 6, k-i \ge 4$** (`research/gessel-lin-zeng-common-interlacer/obstruction.md` & `research/gessel-lin-zeng-descent-roots/COMMON-INTERLACING-BARRIER.md`):
   - **Proved root monotonicity:** Moving a deletion earlier moves every root strictly to the right (`F \prec G` equal-degree preserved by shared $S_M$ for $M > \deg F$). Thus for fixed $(k,i)$, the coordinatewise minimum/maximum roots are attained by $f = A_{k, \{k-i+1,\dots,k\}}$ (late deletions) and $g = A_{k, \{1,\dots,i\}}$ (early deletions).
   - **Proved mean separation:** With $\mu(p) = p'(1)/p(1)$, $\mu(g) - \mu(f) \ge \frac{7}{90} N u v + O(1) \to \infty$ (where $U=2i, T=3(k-i), N=U+T, u=U/N, v=T/N$), and $\mu(g)-\mu(f) \to i/5 > 1$ for fixed $i \ge 6$ and $\to (k-i)/3 > 1$ for fixed $k-i \ge 4$. Since common interlacing implies $0 \le \mu(g)-\mu(f) < 1$, **no common interlacer exists for $\{A_{k,S} : |S|=i\}$** in all these regimes!
   - **Exact finite witnesses:**
     - $(k,i) = (13, 7)$: $234 A_{13,\{7..13\}} + A_{13,\{1..7\}}$ (and $10^5(\cdots) + A_{13,7}$) has degree 18, only 16 real roots (exact Sturm), whereas the equal-weight sum $A_{13,7}$ has all 18 simple negative roots.
     - $(k,i) = (14, 8)$: $R = A_{14,\{1..8\}} + 150 A_{14,\{7..14\}}$ violates Laguerre's inequality $R'(-5)^2 - R(-5)R''(-5) < 0$.
   - **Takeaway:** You CANNOT prove Conjecture 15 by proving arbitrary nonnegative mixtures of $\{A_{k,S}\}$ are real-rooted! Any valid proof MUST use the **equal-weight aggregate** $A_{k,i}$ (or its exact aggregate recurrence)!
2. **Naive bivariate real-stability weightings FAIL** (`research/gessel-lin-zeng/README.md` & `gessel-lin-zeng-followup.md`):
   - Plain $\sum_{i=0}^k A_{k,i}(x) y^i$: fails at $k=2$ ($\operatorname{disc}_y = 4x(x^3-10x^2-7x-4) < 0$ at $x=1$).
   - Plain $\sum_{i=0}^k h_{k,i}(t) y^i$: fails at $k=2, t=1$ ($220+76y+8y^2$, $\operatorname{disc} = -1264 < 0$).
   - $1/i!$-weighted $\sum_{i=0}^k \frac{1}{i!} A_{k,i}(x) y^i$: fails at $k=4, y=-10$ (degree 7, only 5 real roots by Sturm).

---

## 5. Untouched Live Routes & Exact Gaps for `@algebraist` & `@analyst`

1. **Aggregate-row interlacing induction (HIGHEST PRIORITY):**
   In $A_{k,i} = S_{3k-2-i} B_{k,i}$ where $B_{k,i} := S_{3k-3-i} A_{k-1,i} + A_{k-1,i-1}$ (equivalently $h_{k,i} = D_k(D_k h_{k-1,i} + h_{k-1,i-1})$), both summands of $B_{k,i}$ have **the exact same degree** $2k-i-2$ ($\deg A_{k-1,i} = 2k-i-3$).
   - **Key structural observation:** Suppose by induction on $k-1$ that $A_{k-1,i} \prec A_{k-1,i-1}$ for all $1 \le i \le k-1$ (already verified at Ma–Wang's boundaries: $A_{k,1} \prec A_{k,0}$ and $A_{k,k} \prec A_{k,k-1}$).
     1. By Ma–Wang Lemma 2, $A_{k-1,i} \prec S_{3k-3-i} A_{k-1,i}$ automatically. Since also $A_{k-1,i} \prec A_{k-1,i-1}$ (both of degree $2k-i-2$ with positive leading coeffs), **$A_{k-1,i}$ is a common lower interlacer of both summands of $B_{k,i}$** $\implies$ $B_{k,i}$ is automatically real-rooted with $A_{k-1,i} \prec B_{k,i}$!
     2. To propagate $A_{k,i} \prec A_{k,i-1}$ to level $k$, since $A_{k,i} = S_L B_{k,i}$ and $A_{k,i-1} = S_{L+1} B_{k,i-1}$ ($L = 3k-2-i$), **Ma–Wang Lemma 3(1) reduces $A_{k,i} \prec A_{k,i-1}$ directly to $B_{k,i} \prec B_{k,i-1}$!**
     3. Expanding $B_{k,i} = S_{L-1} A_{k-1,i} + A_{k-1,i-1}$ against $B_{k,i-1} = S_L A_{k-1,i-1} + A_{k-1,i-2}$ ($L = 3k-2-i$), **3 of the 4 cross-interlacing pairs are ALREADY AUTOMATIC** from IH + Lemmas 2 & 3(1):
        - $S_{L-1} A_{k-1,i} \prec S_L A_{k-1,i-1}$ (Lemma 3(1) on $A_{k-1,i} \prec A_{k-1,i-1}$),
        - $A_{k-1,i-1} \prec S_L A_{k-1,i-1}$ (Lemma 2),
        - $A_{k-1,i-1} \prec A_{k-1,i-2}$ (IH).
        - **Only 1 cross-pair remains to check/test:** $S_{3k-3-i} A_{k-1,i} \prec A_{k-1,i-2}$ (or directly $B_{k,i} \prec B_{k,i-1}$ without requiring the 4th cross-pair, e.g., via $h_{k,i}$ or interlacing-sequence matrix criteria!). Note that Ma–Wang Remark 1 for $i=4$ is literally a boundary special case of this 4th cross-pair ($S_{3k-7} A_{k-1,4} \prec A_{k-2,0}$ vs $S_{3k-7} A_{k-1,4} \prec A_{k-1,2}$)!
2. **Surviving bivariate real-stability candidate:**
   $B_k(x,y) = \sum_{i=0}^k \binom{k}{i} A_{k,i}(x) y^i$ survived all line checks for $k=1..4$ in `probe.py` — test at larger $k$ and check if a gauge/diagonal transform of $D_j(D_j+y)$ preserves real stability!
3. **Analytic / counterexample guard:**
   Test $A_{k,i} \prec A_{k,i-1}$ and $B_{k,i} \prec B_{k,i-1}$ across large $(k,i)$ (especially deep interior $4 \le i \le k-3$, up to $k \ge 30$). If $B_{k,i} \prec B_{k,i-1}$ holds everywhere, prove it! If it fails at some $(k,i)$, locate the exact failure and check whether $A_{k,i} \prec A_{k,i-1}$ itself holds or if Conjecture 15 has a counterexample.
