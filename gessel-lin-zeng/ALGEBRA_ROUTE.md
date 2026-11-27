# Gessel–Lin–Zeng Conjecture 15 — Algebra Route Archive

**Status:** Archived at swarm checkpoint (`#35–#36`). Every local 1-step operator, bivariate stability-preserver, linear-gauge, and per-subset common-interlacing route is provably blocked, while yielding several complete, independently verified algebraic theorems and a sharp obstruction theorem resolving the exact boundary of Ma–Wang (Oct 2026) Theorem 2.

---

## 1. Exact Setup & Verified Algebraic Identities

All identities below are verified in exact rational/polynomial arithmetic in [`algebra/verify_all_algebra_theorems.py`](algebra/verify_all_algebra_theorems.py) (exit code 0) and independently audited by `@verifier` (`#23`, `#28`, `#34`).

For $k \ge 1$ and $0 \le i \le k$, with $S_L f(x) := (1+Lx)f(x) + x(1-x)f'(x)$ and $P_c f(x) := cf(x) - xf'(x)$:
$$A_{k,i}(x) = S_{3k-2-i}\bigl(S_{3k-3-i} A_{k-1,i}(x) + A_{k-1,i-1}(x)\bigr), \qquad A_{1,0}(x)=1+x,\; A_{1,1}(x)=1.$$
Under the Möbius change of variable $x = \frac{t}{1+t}$ ($(-\infty,0) \leftrightarrow (-1,0)$):
$$\tilde{h}_{k,i}(t) := (1+t)^{2k-i-1} A_{k,i}\!\left(\frac{t}{1+t}\right), \qquad D_j := t(1+t)\partial_t + 1 + jt,$$
$$\tilde{h}_{k,i} = D_{k+1}\bigl(D_{k+1}\tilde{h}_{k-1,i} + \tilde{h}_{k-1,i-1}\bigr), \qquad \tilde{h}_{1,0}=1+2t,\; \tilde{h}_{1,1}=1.$$

### Verified Operator & Structural Identities
1. **Conjugation & commutation:**
   $$D_j f(t) = (1+t)^{2-j}\partial_t\bigl[t(1+t)^{j-1} f(t)\bigr], \qquad (1+t) D_{j+1}^r = D_j^r (1+t), \qquad [D_j, D_{j'}] = (j'-j)\,t(1+t).$$
2. **Polar intertwining:** $P_c S_L = S_{L-1} P_c + (c-L)x$ (hence $P_{L+1} S_{L+1} = S_L P_{L+1}$); $2 A_{k,k-1} = P_{2k} A_{k+1,k+1}$.
3. **Exact $q$-deformed $T$-Turán equivalence:** For $f_0,f_1,f_2$ and $D_j f_r = \Phi f_r + t(1+t)f_r'$ ($\Phi = 1+jt$),
   $$(D_j f_1)^2 - (D_j f_0)(D_j f_2) = A\Phi^2 + B\Phi + C, \qquad B = t(1+t)A',$$
   with discriminant $4AC - B^2 = t^2(1+t)^2\bigl(4 W(f_0,f_1)W(f_1,f_2) - W(f_0,f_2)^2\bigr)$ (generalizing with factor $q$ to $4A_q C_q - B_q^2 = q\,t^2(1+t)^2[4 W_{01}W_{12} - q W_{02}^2]$).
4. **Equal-degree $D$-lifted recurrence:** Setting $E_{k,i} := D_{k+1}^i \tilde{h}_{k-1,i}$ ($\deg = 2k-3$) and $\mathcal{B}_k^{(E)}(t,y) := \sum_{i=0}^{k-1} E_{k,i}(t) y^i$,
   $$D_{k+1}^i \tilde{h}_{k,i} = D_{k+1}^2\bigl(E_{k,i} + E_{k,i-1}\bigr), \qquad \sum_{i=0}^k D_{k+1}^i \tilde{h}_{k,i}(t)\, y^i = D_{k+1}^2\bigl((1+y)\,\mathcal{B}_k^{(E)}(t,y)\bigr).$$
5. **Ehrhart representation:** $\sum_{n\ge 0} p_{k,i}(n) x^n = A_{k,i}(x)/(1-x)^{3k-i}$ with $p_{k,i}(n) = \sum_{1\le m_1\le\dots\le m_k = n+1} (\prod m_j) e_{k-i}(m_1,\dots,m_k)$.
6. **Exact `first_base` refinement:** $P_{k,i,r} = S_{L-1}(S_{L-2} P_{k-1,i,r} + P_{k-1,i-1,r})$ ($1\le r<k$), $P_{k,i,k} = S_L(x A_{k-1,i}) + x(S_{L-2} A_{k-1,i} + A_{k-1,i-1})$, summing to $A_{k,i}$ with $P_{k,i,1} \prec \dots \prec P_{k,i,k}$ (3,003/3,003 pairs through $k=12$).
7. **Constant-$b$ $t=0$ Wronskian ratio:** For $D_b(D_b+y)$ with constant $b$, $J_i(0) = \frac{8}{9}\frac{(i+1)(k-i+1)}{i(k-i)}$ (exact asymptote matching GLZ's $J_{\text{raw}}(0)$ crossing below $1$ at $k=35, i=17$: `0.997128`).

---

## 2. Proved & Independently Verified Positive Theorems

### Theorem 1 (Sharpest Quad Theorem B — verified `#28`)
Let $A, B, C, E$ be differentiable real functions on $I$ with $W_{AB}, W_{BC}, W_{CE}, W_{AC}, W_{BE} > 0$ and $W_{AE} \ge 0$ on $I$, and suppose $J(A,B,C) \ge 1+\varepsilon$ and $J(B,C,E) \ge 1+\varepsilon$ ($\varepsilon \ge 0$), where $J(X,M,Y) := \frac{4 W(X,M) W(M,Y)}{W(X,Y)^2}$. Then $W(A+B, C+E) = W_{AC}+W_{AE}+W_{BC}+W_{BE} > 0$ unconditionally, and **with zero upper bound on the Plücker ratio $R = W_{AB}W_{CE}/W_{BC}^2$**,
$$J(A+B,\, B+C,\, C+E) \;\ge\; \frac{4(1+\varepsilon)^2}{(2+\varepsilon)^2} \;\ge\; 1.$$

### Theorem 2 (Universal Cubic Real-Stability Addition — 15-line proof, verified `#28`)
Let $v_0, v_1, v_2, v_3, v_4 \in \mathbb{R}^2$ satisfy $W(v_r, v_s) > 0$ for all $0 \le r < s \le 4$. If the cubics $P_1(y) = v_0 - 3v_1 y + 3v_2 y^2 - v_3 y^3$ and $P_2(y) = v_1 - 3v_2 y + 3v_3 y^2 - v_4 y^3$ both have positive semidefinite $3\times 3$ Bezoutiants (allowing mixed-sign $y$-roots), then $P_1(y) + P_2(y)$ also has positive semidefinite Bezoutiant.
- **Key mechanism:** Angular half-plane order forces at most one sign change in any linear functional $L(v_r)$, reducing all sign patterns to $(+,+,+,+,+)$, $(-,+,+,+,+)$, $(-,-,+,+,+)$ up to negation and $y\mapsto 1/y$, where Obreschkoff + IVT close every case; $d=3$ is the unique maximum degree with $\min(p, d+1-p) \le 2$.
- **Sharpness witness:** $L = (-3000, 1, 10, 1, -3000)$ breaks the statement without the Wronskian cone (exact Sturm: 1 real root).

### Theorem 3 (Zero-Loss $D_m$ Preservation — verified `#23`)
For every $m \ge 1$, $D_m$ preserves $q$-Turán ($J \ge q$ on $\mathbb{R}$) and cubic real-stability with zero loss.

### Theorem 4 (Exact Sharpness of Ma–Wang Oct 2026 Theorem 2 for ALL $4 \le i \le k-3$ — verified `#34`)
In `arXiv:2610.03111`, Ma and Wang proved real-rootedness of $A_{k,i}(x)$ for $i \in \{0,1,2,3,k-2,k-1,k\}$ via common interlacers of the subset summands $\{A_{k,S}(x) : |S|=i\}$. By root monotonicity across subsets, $\{A_{k,S} : |S|=i\}$ admits a common interlacer if and only if $A_{k,\text{late}} \prec A_{k,\text{early}}$, where $\text{early} = \{1,\dots,i\}$ and $\text{late} = \{k-i+1,\dots,k\}$.

We prove that **Ma–Wang Theorem 2 is 100% sharp across the entire open range $4 \le i \le k-3$**: non-real nonnegative mixtures $R(x) = \lambda A_{k,\text{late}}(x) + A_{k,\text{early}}(x)$ exist in every boundary family $i=4$, $i=5$, and $i=k-3$, certified via exact rational/integer Laguerre inequalities $R'(x)^2 - R(x)R''(x) < 0$ in [`algebra/exact_laguerre_barriers.py`](algebra/exact_laguerre_barriers.py) and [`algebra/clean_laguerre_certs.py`](algebra/clean_laguerre_certs.py) (independently verified in exact `Fraction`/Cauchy-index arithmetic by `@verifier` `#34`):

| Boundary Family | Exact First Fail $k$ ($k-1$ passes) | Degree | Exact Mixture $R(x)$ | Evaluation Point | $R'(x)^2 - R(x)R''(x)$ |
| :--- | :---: | :---: | :--- | :---: | :---: |
| $i = 4$ | $k = 31$ ($k=30$ passes) | 57 | $932\,A_{31,\{28..31\}} + A_{31,\{1..4\}}$ | $x = -23$ (integer) | $-9.065 \times 10^{279} < 0$ |
| $i = 5$ | $k = 16$ ($k=15$ passes) | 26 | $269\,A_{16,\{12..16\}} + A_{16,\{1..5\}}$ | $x = -11$ (integer) | $-2.318 \times 10^{96} < 0$ |
| $i = k-3$ | $k = 28$ ($k=27$ passes) | 30 | $272\,A_{28,\{4..28\}} + A_{28,\{1..25\}}$ | $x = -73/10$ | $-5.026 \times 10^{111} < 0$ |
| $i = k-3$ (int pt) | $k = 29$ | 31 | $349\,A_{29,\{4..29\}} + A_{29,\{1..26\}}$ | $x = -8$ (integer) | $-4.898 \times 10^{118} < 0$ |

---

## 3. Summary of Proved Obstructions & Remaining Irreducible Gap

From identity (4) above, the Pascal + outer-$D_{k+1}^2$ step preserves real-stability with zero loss. The **sole** remaining obstruction is the index shift $D_{k+1}^i \mapsto D_{k+2}^i$ (Stage 3), where every local/black-box closure mechanism is now ruled out:
1. **Universal operator stability disproved (`#31`, certified `#34`):** Even with $a_j = w_j = 1$ and integer $b_j \ge 2$, both $\prod (D_{b_j}+y)(1)$ and $\prod D_{b_j}(D_{b_j}+y)(1)$ break interlacing.
2. **No linear gauge eliminates Stage 3 (`#32`):** Five-line impossibility proof (any linear gauge forces $D_{k+2}D_{k+1} = D_{k+1}^2$); bivariate symbol of $D_{k+1}^i \mapsto D_{k+2}^i$ fails Borcea–Brändén at $k=4$.
3. **Black-box 1-step Stage 3 and combined GLZ map $\mathcal{G}_j$ fail (`#24`, `#30`, [`algebra/final_go_pivot_check.py`](algebra/final_go_pivot_check.py)):** Even at GLZ's actual boundary scale $\lambda = 3/j$ ($j=80$), inputs with $J_{\text{in}} \ge 1 + 2/j = 1.025$ on all of $\mathbb{R}$ drop to $J_{\text{out}} = 0.983221 < 1$.
4. **Two-Regime Split at $t = -1/(k+2)$ (`#34`, [`algebra/surgical_checkpoint.py`](algebra/surgical_checkpoint.py)):** Empirically holds through $k=45$ ($A,B,C > 0$ wherever $J_{\text{raw}} \le 1$, which occurs only at $|t| \lesssim 2^{-k}$), but has no known algebraic induction preserving $J_{\text{raw}} > 1$ on the proper subinterval $(-1,-1/(k+2)]$ or $A' < 0$ on $|t| \lesssim 2^{-k}$.
5. **Ruled-out auxiliary routes:** Bivariate $B_k(x,y)$ ($k=10$), $B_k(x,-xy)$ ($k=20$), $Q_k(t,y)$ ($k=35$), $Q_{k,q}$ proper position (all $q>0$), $W(E_{k,0},E_{k,r})>0$ ($r=8,9$), Quad $R \le 4$ closure ($R_{\text{out}} = 81/16$), Derivative Sandwich implying 4th cross-pair (486/500 fail), Brändén–Solus symmetric decomposition, Ma–Wang Prop 1 generalization ($k=5$), and fixed-$k$ degree-$\le 2$ differential syzygy ($k=6$, rank 12/18).
