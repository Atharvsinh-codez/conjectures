# Algebraic & Analytic Route to Bardakov–Kozlovskaya Conjecture 3.3(2)

**Target Conjecture ([arXiv:2305.07261](https://arxiv.org/abs/2305.07261), Conjecture 3.3(2)):**
For every even integer $n = 2m \ge 2$, define the symmetric homogeneous degree-$n$ integer polynomial
$$p_{2m}(x, y, z) := \prod_{k=0}^{2m-1}\Bigl(z - (\alpha + \zeta^k \beta)^{2m}\Bigr) = \sum_{\substack{a, b, c \ge 0 \\ a + 2b + 3c = 2m}} c_{a,b,c}^{(2m)}\, e_1^a e_2^b e_3^c \qquad (\alpha^{2m} = x,\; \beta^{2m} = y,\; \zeta = e^{2\pi i/(2m)}),$$
where $e_1 = x+y+z$, $e_2 = xy+yz+zx$, $e_3 = xyz$ (bijection with Bardakov–Kozlovskaya's indexing: $A_{k_1,k_2,k_3} = c_{a,b,c}^{(2m)}$ with $a = k_1 - k_2$, $b = k_2 - k_3$, $c = k_3$, $a \equiv c \pmod 2$).
Conjecture 3.3(2) asserts that $c_{a,b,c}^{(2m)} \ne 0$ for every even $n = 2m \ge 2$ and every $(a,b,c) \in \mathbb{Z}_{\ge 0}^3$ with $a + 2b + 3c = 2m$.

**Stronger Theorem (Exact Sign Law & Slice Real-Rootedness):**
For every even $n = 2m \ge 2$ and every $0 \le c \le \lfloor 2m/3 \rfloor$, set $D_c := \lfloor (2m - 3c)/2 \rfloor$ and define the univariate $c$-slice polynomial
$$Q_c^{(2m)}(1, u) := \sum_{b=0}^{D_c} c_{2m-3c-2b, b, c}^{(2m)}\, u^b, \qquad \text{so that } P_{2m}(e_1, e_2, e_3) = \sum_{c=0}^{\lfloor 2m/3 \rfloor} e_1^{2m-3c} Q_c^{(2m)}\!\Bigl(1, \frac{e_2}{e_1^2}\Bigr)\, e_3^c.$$
1. For every $1 \le c \le \lfloor (2m-3)/3 \rfloor$, $Q_c^{(2m)}(1, u)$ has $D_c$ simple real roots — all strictly positive in $(0, +\infty)$ (with at least $D_c - 1$ in $(0, 1/4)$) when $c$ is even, and all strictly negative in $(-\infty, 0)$ when $c$ is odd.
2. Consequently, by Vieta's formulas and the exact boundary signs at $c = 0$, $a = 0$, $b = 0$, every coefficient $c_{a,b,c}^{(2m)}$ is nonzero and satisfies the exact sign law
$$\operatorname{sign}\!\bigl(c_{a,b,c}^{(2m)}\bigr) = (-1)^{b(c+1)+c} = \begin{cases} (-1)^b & \text{if } c \text{ is even,} \\ -1 & \text{if } c \text{ is odd.} \end{cases}$$

**Status of Proof Components:**
- **Finite Bridge ($n \le 130$):** Certified in exact integer arithmetic by `@verifier` (`#39`, `#57`: all 65 even degrees $n = 2, 4, \dots, 130$, 33,417 exact integer coefficients, 0 zeros, 0 sign violations).
- **All 4 Boundary Families ($c=0, a=0, b=0, c=1$) for ALL even $n \ge 2$:** 100% proved below and independently audited/certified by `@verifier` (`#42`, `#45`).
- **Exact Resultant Representations (Depressed-Cubic & Half-Size Lucas–Chebyshev):** Proved in Section 1 and certified by `@verifier` (`#45`).
- **Regime 1 ($\kappa(n, c, j) \ge 2.5$, ALL even $n \ge 84$, ALL strictly interior $(c, j)$):** 100% closed-form elementary analytic proof complete in Section 4 below (global worst closed-form margin **`+0.4699`** at $\kappa_0 = 3.0$ and **`+0.0748`** at $\kappa_0 = 2.5$, zero computer search; independently referee-signed-off by `@verifier` `#59`, `#61`).
- **Regime 2 ($\kappa(n, c, j) < 2.5$, ALL even $n \ge 132$, strictly interior $2 \le c \le \lfloor(n-3)/3\rfloor$, $j \in \mathcal{J}_{n,c}$):** Complete zero-crevice 6-piece decomposition (overlapping seams locked in `#65`–`#69`):
  - **R2.1 (Universal Cubic Perturbation Reduction, $\lambda \le 3$), R2.2 ($c\kappa^2 \ge 42$), R2.3 ($c q \kappa \le 0.90$), and G3 ($D \ne 0$ Major Arcs):** Owned by `@analyst` in `ANALYTIC_ROUTE.md` (R2.1–R2.3 audited by `@verifier` `#71`).
  - **Piece G2 ($D = 0$ Major Arcs: $d = \gcd(c, j) \ge 2$, $q = c/d \ge 1$, $\kappa < 2.5$ for $q = 1$ and $q^2\kappa \le 2.0$ for $q \ge 2$, $c q \kappa \ge 0.85$, $c\kappa^2 < 42$):** 100% proved in closed form in Section 5.1 below and certified in 50-digit outward-rounded `mpmath.iv` ([`certify_G2_D0.py`](algebra/certify_G2_D0.py), zero loops over $(n, c, j)$): **$I_{\text{full}} / L^1(\tilde{\lambda}) \ge +0.00176348$ ($q \ge 2$) and $\ge +0.216651$ ($q = 1$)**.
  - **Piece G1 (Dirichlet Minor Arcs: $\mathcal{G}^* \setminus (\text{G2} \cup \text{G3})$):** Proved in closed form in Section 5.2 below and certified in 50-digit outward-rounded `mpmath.iv` ([`certify_G1_minor.py`](algebra/certify_G1_minor.py), zero loops over $(n, c, j, r)$).

---

## 1. Exact Structural Resultants

### Theorem 1.1 (Depressed-Cubic Resultant)
Let $t^3 - u t - v = 0$ have roots $t_1, t_2, t_3$ with $e_1(t) = 0, e_2(t) = -u, e_3(t) = v$, and set $\eta := v^2$, $x = t_1^{2m}$, $y = t_2^{2m}$, $z = t_3^{2m}$. Since $e_3(x,y,z) = (t_1 t_2 t_3)^{2m} = v^{2m} = \eta^m$, the elementary symmetric polynomials $e_1(x,y,z)$ and $e_2(x,y,z)$ are polynomials in $(u, \eta)$ given explicitly by Waring's formula:
$$\Pi_{2m}(u, \eta) := t_1^{2m} + t_2^{2m} + t_3^{2m} = \sum_{r=0}^{\lfloor m/3 \rfloor} \frac{2m}{m - r}\binom{m - r}{2r} u^{m - 3r}\eta^r,$$
$$\Sigma_{2m}(u, \eta) := (t_1 t_2)^{2m} + (t_2 t_3)^{2m} + (t_3 t_1)^{2m} = \sum_{j=0}^{\lfloor 2m/3 \rfloor} (-1)^j \frac{m}{m - j}\binom{2m - 2j}{j} u^{2m - 3j}\eta^j.$$
Then:
- When $3 \nmid m$: $P_{2m}(e_1, e_2, \eta^m) = \operatorname{Res}_u\bigl(\Pi_{2m}(u, \eta) - e_1,\; \Sigma_{2m}(u, \eta) - e_2\bigr)$ (exact sign $+1$).
- When $m = 3K$: $\Pi_{6K}$ and $\Sigma_{6K}$ are polynomials in $w = u^3$, and $P_{6K}(e_1, e_2, e_3) = \prod_{\omega^3 = 1} R_K(e_1, e_2, \omega e_3^{1/3})$, where $R_K(e_1, e_2, \eta) := \operatorname{Res}_w(\tilde\Pi_{6K}(w, \eta) - e_1, \tilde\Sigma_{6K}(w, \eta) - e_2)$.

*Proof.* The equation $z = (\alpha + \zeta^k\beta)^{2m}$ for some $k$ is equivalent to $t_1 + t_2 + t_3 = 0$ for some $2m$-th roots $t_1 = x^{1/2m}, t_2 = y^{1/2m}, t_3 = z^{1/2m}$ (since $-1 = \zeta^m$ is a $2m$-th root of unity). Such $(t_1, t_2, t_3)$ are precisely the roots of $t^3 - u t - v = 0$ with $v^{2m} = e_3(x,y,z)$ and $u = -e_2(t_1, t_2, t_3)$, which maps $(u, v^2) \mapsto (\Pi_{2m}(u, v^2), \Sigma_{2m}(u, v^2), v^{2m})$. Eliminating $u$ yields the resultant; comparing degrees in $e_1, e_2, e_3$ and leading terms at $\eta = 0$ ($\Pi_{2m}(u,0) = 2 u^m$, $\Sigma_{2m}(u,0) = u^{2m}$, giving $\operatorname{Res}_u(2u^m - e_1, u^{2m} - e_2) = (e_1^2 - 4e_2)^m$ for $3 \nmid m$) fixes the normalization. $\square$

### Theorem 1.2 (Half-Size Lucas–Chebyshev Resultant)
Let $L_n(u) := \sum_{r=0}^{\lfloor n/2 \rfloor} \frac{n}{n-r}\binom{n-r}{r} u^r$ be the Lucas–Chebyshev polynomial satisfying $L_n(t(1-t)) = t^n + (1-t)^n$. For every even $n = 2m \ge 2$ and all $p, q \in \mathbb{C}$,
$$P_{2m}(1 + q - p,\, q,\, p) = \operatorname{Res}_u\!\bigl(L_{2m}(u) + u^{2m} - q,\; u^{2m} - p\bigr) = \prod_{k=0}^{2m-1}\Bigl(L_{2m}(\zeta^k p^{1/2m}) + p - q\Bigr).$$

*Proof.* Set $z = 1$ in $p_{2m}(x, y, z)$ and parameterize $x^{1/2m} = -t$, $y^{1/2m} = -(1-t)$ with $u := t(1-t)$. Then $x + y = t^{2m} + (1-t)^{2m} = L_{2m}(u)$ and $xy = u^{2m}$, so $e_3 = p = u^{2m}$, $e_2 = q = xy + (x+y) = L_{2m}(u) + u^{2m}$, and $e_1 = x + y + 1 = 1 + q - p$. Eliminating $u$ over the $2m$ roots $u = \zeta^k p^{1/2m}$ yields the identity; at $p = 0$, $L_{2m}(0) = 1$ gives $(1 - q)^{2m}$, matching $P_{2m}(1 + q, q, 0) = ((1+q)^2 - 4q)^m = (1 - q)^{2m}$ and confirming unit normalization. $\square$

---

## 2. Complete Proofs of All 4 Boundary Families ($c=0, a=0, b=0, c=1$)

All four boundary families hold for **all even $n = 2m \ge 2$** (audited and certified by `@verifier` in `#42`, `#45`):

1. **$c = 0$ ($e_3 = 0$):** Setting $z = 0$ gives $p_{2m}(x, y, 0) = \prod_{k=0}^{2m-1} (\alpha + \zeta^k \beta)^{2m} = (x - y)^{2m} = ((x+y)^2 - 4xy)^m = (e_1^2 - 4e_2)^m$. Hence $c_{2m-2b, b, 0}^{(2m)} = (-4)^b \binom{m}{b}$, with exact sign $(-1)^b$.

2. **$a = 0$ ($e_1 = 0$) and $b = 0$ ($e_2 = 0$):** Setting $q = p - 1$ (for $e_1 = 0$) and $q = 0$ (for $e_2 = 0$) in Theorem 1.2 yields exact univariate products over the roots of $L_{2m}(u) + 1 = 0$ (resp. $L_{2m}(u) + u^{2m} = 0$). Under $u = -\frac{1}{4\cos^2\theta}$, exact differentiation gives
   $$\frac{d}{d\theta} L_{2m}\!\Bigl(-\frac{1}{4\cos^2\theta}\Bigr) = -\frac{4m\sin((2m-1)\theta)}{(2\cos\theta)^{2m+1}}, \qquad \frac{d}{d\theta}\Bigl[(2\cos\theta)^{4m}\bigl(L_{2m}(u) + u^{2m}\bigr)\Bigr] = -4m(2\cos\theta)^{2m-1}\sin((2m+1)\theta),$$
   which isolates all $\lfloor m/3 \rfloor$ roots of $P_{2m}(0, e_2, e_3)$ in $e_3^2/e_2^3$ and all $\lfloor 2m/3 \rfloor$ roots of $P_{2m}(e_1, 0, e_3)$ in $e_3/e_1^3$ strictly on $(0, +\infty)$ (proved with `@analyst` `#40`, certified by `@verifier` `#42`). Hence:
   - $c_{2m-3c, 0, c}^{(2m)} = (-1)^c |c_{2m-3c, 0, c}^{(2m)}| \ne 0$ for all $0 \le c \le \lfloor 2m/3 \rfloor$;
   - $c_{0, m-3q, 2q}^{(2m)} = (-1)^{m-q} |c_{0, m-3q, 2q}^{(2m)}| \ne 0$ for all $0 \le q \le \lfloor m/3 \rfloor$.

3. **$c = 1$ (Self-Contained Proof for ALL even $n = 2m \ge 4$):**
   Set $x = 1, y = u, z = 0$ (so $e_1 = 1+u$, $e_2 = u$, $e_3 = 0$) and apply $\frac{\partial}{\partial z}\big|_{z=0}$ to $P_n(e_1, e_2, e_3) = (e_1^2 - 4e_2)^m + e_3 Q_1^{(2m)}(e_1, e_2) + O(e_3^2)$.
   - On the LHS: $\frac{\partial e_1}{\partial z} = 1$, $\frac{\partial e_2}{\partial z} = 1+u$, $\frac{\partial e_3}{\partial z} = u$, so
     $$\Bigl[\frac{\partial P_n}{\partial z}\Bigr]_{z=0} = -2m(1+u)(1-u)^{2m-2} + u\, Q_1^{(2m)}(1+u, u).$$
   - On the RHS: $p_n(1, u, z) = \prod_{k=0}^{n-1}(z - (1 + \zeta^k u^{1/n})^n)$, so
     $$\Bigl[\frac{\partial p_n}{\partial z}\Bigr]_{z=0} = -(1 - u)^n \sum_{k=0}^{n-1} \frac{1}{(1 + \zeta^k u^{1/n})^n} = -2m\, A_{2m-1}^{(2m)}(u),$$
     where $F_r(u) := (1-u)^r S_r(u)$, $S_r(u) := \frac{1}{n}\sum_{k=0}^{n-1} (1 + \zeta^k u^{1/n})^{-r}$, and $A_{2m-1}^{(2m)}(u) := F_{2m}(u)$.

   Differentiating $S_{r-1}(u)$ gives $S_r(u) = S_{r-1}(u) + \frac{n u}{r - 1} S'_{r-1}(u)$, and multiplying by $(1 - u)^r$ yields the exact recurrence (for $r \ge 2$, with $F_1(u) = 1$):
   $$(r - 1)\, F_r(u) = n u(1 - u)\, F'_{r-1}(u) + (r - 1)\bigl(1 + (n - 1)u\bigr)\, F_{r-1}(u).$$
   By Rolle's theorem on the positive-coefficient degree-$(r-1)$ polynomial $F_r(u)$ (which satisfies $F_r(0) = 1$, $F_r(1) = n^{r-1} > 0$, while the derivative term alternates sign at the simple negative roots of $F_{r-1}$), **$F_r(u)$ has $r - 1$ simple strictly negative real roots for every $1 \le r \le n$**.
   In particular, $A_{2m-1}^{(2m)}(u) = F_{2m}(u)$ has $2m - 1$ simple strictly negative real roots and is palindromic of odd degree $2m - 1$ ($F_{2m}(1/u) = u^{-(2m-1)} F_{2m}(u)$). Pairing reciprocal roots $(r_j, r_j^{-1})$ and dividing by $1 + u$ gives
   $$\frac{A_{2m-1}^{(2m)}(u)}{1 + u} = \prod_{j=1}^{m-1}\bigl((1 + u)^2 + c_j u\bigr) = \sum_{r=0}^{m-1} \gamma_r\, (1 + u)^{2m - 2 - 2r}\, u^r, \qquad c_j = -(r_j + r_j^{-1}) - 2 > 0,$$
   where $\gamma_r = e_r(c_1, \dots, c_{m-1}) > 0$ ($\gamma_0 = 1$). Since $(1 - u)^{2m-2} = \sum_{r=0}^{m-1} (-4)^r \binom{m-1}{r} (1 + u)^{2m-2-2r} u^r$, subtracting cancels the $r = 0$ term and divides cleanly by $u$:
   $$Q_1^{(2m)}(1+u, u) = 2m \sum_{b=0}^{m-2} \Bigl((-4)^{b+1}\binom{m-1}{b+1} - \gamma_{b+1}\Bigr)(1+u)^{2m - 3 - 2b} u^{b},$$
   yielding the exact closed-form identity $c_{2m-3-2b,\, b,\, 1}^{(2m)} = 2m\bigl((-4)^{b+1}\binom{m-1}{b+1} - \gamma_{b+1}\bigr)$ for all $0 \le b \le m - 2$.
   By **Maclaurin's inequality** for $c_1, \dots, c_{m-1} > 0$, $(\gamma_r / \binom{m-1}{r})^{1/r} \ge \gamma_{m-1}^{1/(m-1)}$ for all $1 \le r \le m - 1$, so $c_{2m-3-2b, b, 1}^{(2m)} < 0$ for all $b$ follows from the single inequality $\gamma_{m-1} > 4^{m-1}$.
   Differentiating $A_{2m-1}^{(2m)}(u)$ at $u = -1$ (L'Hôpital's rule on $\frac{A_{2m-1}^{(2m)}(u)}{1+u}$, certified by `@verifier` `#45`) gives the exact alternating trigonometric sum
   $$\gamma_{m-1} = \frac{1}{2m}\sum_{k=0}^{m-1} (-1)^k \frac{\cos\theta_k}{\sin^{2m+1}\theta_k}, \qquad \theta_k = \frac{(2k+1)\pi}{4m}.$$
   Since $\frac{\cos\theta}{\sin^{2m+1}\theta}$ is strictly decreasing on $(0, \pi/2)$ and $\frac{\sin\theta_1}{\sin\theta_0} = 3 - 4\sin^2\frac{\pi}{4m} \ge 3 - 4\sin^2\frac{\pi}{8} = 1 + \sqrt{2} > 2$ for $m \ge 2$, the alternating sum exceeds its first two terms:
   $$\gamma_{m-1} > \frac{\cos(\pi/8)}{2m}\Bigl(\frac{4m}{\pi}\Bigr)^{2m+1}\Bigl(1 - \frac{1}{2^{2m+1}}\Bigr) > 4^{m-1} \qquad \text{for all } m \ge 2.$$
   (At $m = 2$, $\gamma_1 = 28 > 4$; for $m \ge 3$, $\frac{4m}{\pi} > 3.81$, so $(\frac{4m}{\pi})^{2m} > 14.5^m \gg 4^{m-1}$.)
   Hence $\gamma_r > 4^r \binom{m-1}{r}$ for all $1 \le r \le m - 1$, proving $c_{2m-3-2b, b, 1}^{(2m)} < 0$ for all $m \ge 2$ and all $0 \le b \le m - 2$. $\square$

---

## 3. Universal Sample Points & Exact 1D Circle Integral

### Theorem 3.1 (Reduction to Strictly Interior Sample Points)
Fix even $n = 2m \ge 4$ and $2 \le c \le \lfloor (n - 3)/3 \rfloor$ (note $n - 3c \in \{0, 1, 2\}$ contains only $a = 0$ and $b = 0$ boundary coefficients, already proved in Section 2). Set $D_c := \lfloor (n - 3c)/2 \rfloor$ and define the strictly interior index range (incorporating `@verifier`'s double-corner elimination `#51`):
$$j \in \mathcal{J}_{n, c} := \bigl\{\, c + 1,\; c + 2,\; \dots,\; m - \lceil c/2 \rceil - \mathbf{1}_{c\text{ even}} \,\bigr\}.$$
For each $j \in \mathcal{J}_{n, c}$, define
$$\theta_c := \frac{c\pi}{n}, \quad \psi_j := \frac{j\pi}{n}, \quad \alpha_1 := \frac{(j + c)\pi}{n}, \quad r_1 := \sin\alpha_1, \quad r_2 := \sin\psi_j, \quad r_3 := \sin\theta_c,$$
$$s_j := \frac{r_2}{r_1} \in (0, 1), \quad R_y := s_j^n, \quad R_z := \Bigl(\frac{r_3}{r_1}\Bigr)^n, \quad t_j := (-1)^c R_y, \quad u_j := \frac{t_j}{(1 + t_j)^2}.$$
**If $(-1)^j Q_c^{(2m)}(1, u_j) > 0$ for every $j \in \mathcal{J}_{n, c}$, then $Q_c^{(2m)}(1, u)$ has $D_c$ simple real roots in $(0, +\infty)$ (for even $c$) or $(-\infty, 0)$ (for odd $c$), and hence every coefficient $c_{n-3c-2b, b, c}^{(2m)}$ ($0 \le b \le D_c$) is nonzero with exact sign $(-1)^{b(c+1)+c}$.**

*Proof.* Since $\frac{d}{dx}\frac{\sin x}{\sin(x + \theta_c)} = \frac{\sin\theta_c}{\sin^2(x + \theta_c)} > 0$, $s_j$ is strictly increasing in $j \in \mathcal{J}_{n, c}$ and satisfies $0 < s_j < 1$ (since $0 < j < \min(j+c, n-j-c) \le n/2$).
- **Even $c = 2q$:** Here $D_c = m - 3q$ and $|\mathcal{J}_{n,c}| = (m - q - 1) - 2q = D_c - 1$. As $j$ increases, $t_j = s_j^n \in (0, 1)$ strictly increases, so $0 < u_{c+1} < u_{c+2} < \dots < u_{m-q-1} < \frac{1}{4}$. Combining the $b = 0$ boundary $Q_c^{(2m)}(1, 0) = c_{n-3c,0,c}^{(2m)} > 0$ (sign $(-1)^c$), the $D_c - 1$ strictly alternating evaluations $(-1)^j Q_c^{(2m)}(1, u_j) > 0$, and the $a = 0$ leading-term limit $\lim_{u \to +\infty} \operatorname{sign}(Q_c^{(2m)}(1, u)) = \operatorname{sign}(c_{0, m-3q, 2q}^{(2m)}) = (-1)^{m-q}$ (opposite to $(-1)^{m-q-1}$) yields $D_c$ sign changes on $(0, +\infty)$ — hence $D_c$ simple positive real roots (at least $D_c - 1$ in $(0, 1/4)$).
- **Odd $c = 2q + 1$:** Here $D_c = m - 3q - 2$ and $|\mathcal{J}_{n,c}| = (m - q - 1) - (2q + 1) = D_c$. As $j$ increases, $t_j = -s_j^n \in (-1, 0)$ strictly decreases, and $\frac{d}{dt}\frac{t}{(1+t)^2} = \frac{1-t}{(1+t)^3} > 0$ on $(-1, 0)$, so $-\infty < u_{m-q-1} < \dots < u_{c+1} < 0$. Combining these $D_c$ alternating evaluations with the $b = 0$ boundary $Q_c^{(2m)}(1, 0) < 0$ (sign $(-1)^c$) yields $D_c$ sign changes strictly inside $(u_{m-q-1}, 0) \subset (-\infty, 0)$, hence $D_c$ simple negative real roots.

In both parities, Vieta's formulas force every coefficient of $Q_c^{(2m)}(1, u)$ to be nonzero with exact sign $(-1)^{b(c+1)+c}$. $\square$

### Theorem 3.2 (Exact 1D Circle Integral & 4-Factor Decomposition)
Fix even $n = 2m \ge 4$, $2 \le c \le \lfloor (n-3)/3 \rfloor$, and $j \in \mathcal{J}_{n, c}$. For $w = e^{i\phi}$ ($|w| = 1$), consider the cubic equation
$$T(T - 1)(T - t) = t R_z w \qquad (t := t_j = (-1)^c R_y).$$
Assume Rouché's theorem isolates unique simple roots $x(w) = 1 + \delta_x(w)$ and $y(w) = t(1 + \tilde\delta_y(w))$ holomorphic on a neighborhood of $|w| \le 1$ with $|\delta_x| < 1$ and $|\delta_y| < 1$, where
$$1 + \delta_y(w) := \frac{y(w)}{t\,x(w)}, \qquad 1 + \delta_z(w) := \frac{z(w)}{x(w) R_z w} = \frac{1}{(1 + \delta_x(w))^3 (1 + \delta_y(w))}.$$
Let $\varepsilon_s(w) := (1 + \delta_y(w))^{1/n} - 1$ (principal branch). Let the $n$ angles $\eta_k \in (-\pi, \pi]$ be $\frac{2k\pi}{n}$ (even $c$) or $\frac{(2k+1)\pi}{n}$ (odd $c$), and partition them into:
- **Small** ($|\eta_k| < \theta_c$, count $c - 1$),
- **Boundary** ($|\eta_k| = \theta_c$, count $2$, namely $\eta = \pm\theta_c$),
- **Large** ($|\eta_k| > \theta_c$, count $n - c - 1$).

Set $w_k := (1 - s e^{i\eta_k})^n$, $a_k := \frac{-s e^{i\eta_k}}{1 - s e^{i\eta_k}}$, $\tilde{w}_k := w_k (1 + a_k \varepsilon_s)^n$, $1 + \varepsilon_\pm(w) := (1 + a(\pm\theta_c)\varepsilon_s(w))^n$, and $Z := (-1)^j w^{-1}$. Then
$$(-1)^j (1 + t)^{n - 3c} Q_c^{(2m)}(1, u_j) = \frac{C_{m,c,j}}{2\pi}\int_{-\pi}^{\pi} (1 + \Delta_{\text{mult}}(e^{i\phi}))\, B(e^{i\phi})\, d\phi, \qquad C_{m,c,j} := \frac{R_z \prod_{\text{large}} |w_k|}{R_y^c} > 0,$$
where $1 + \Delta_{\text{mult}}(w) := F_1(w)\, F_2(w)\, F_3(w)$ and
$$F_1(w) := (1 + \delta_x)^{n - 3c}, \qquad F_2(w) := \frac{\prod_{\text{large}} (1 + a_k \varepsilon_s)^n}{(1 + \delta_y)^c}, \qquad F_3(w) := \prod_{\text{small}}\!\Bigl(1 - \frac{\tilde{w}_k}{z/x}\Bigr)\prod_{\text{large}}\!\Bigl(1 - \frac{z/x}{\tilde{w}_k}\Bigr),$$
$$B(w) := \frac{(-Z^{-1})\bigl((1+\delta_z) - Z(1+\varepsilon_+)\bigr)\bigl((1+\delta_z) - Z(1+\varepsilon_-)\bigr)}{1+\delta_z}$$
$$= (2 - Z - Z^{-1}) - Z^{-1}\delta_z + (\varepsilon_+ + \varepsilon_-) - Z\Bigl((1+\varepsilon_+)(1+\varepsilon_-)(1+\delta_x)^3(1+\delta_y) - 1\Bigr).$$

*Proof.* Along the circle $e_1 = 1+t$, $e_2 = t$, $e_3 = t R_z w$, Cauchy's integral formula yields
$$(1 + t)^{n - 3c} Q_c^{(2m)}(1, u_j) = Q_c^{(2m)}(1+t, t) = \frac{1}{2\pi}\int_{-\pi}^{\pi} \frac{p_n(x(e^{i\phi}), y(e^{i\phi}), z(e^{i\phi}))}{(t R_z e^{i\phi})^c}\,d\phi.$$
Since $y/x = (-1)^c s^n (1 + \varepsilon_s)^n$, choosing $\alpha = x^{1/n}$ and $\beta = \alpha s(1+\varepsilon_s) e^{i\pi c/n}$ gives
$$\frac{p_n(x, y, z)}{(xyz)^c} = x^{n-3c}\frac{\prod_{k=0}^{n-1}\bigl((z/x) - \tilde{w}_k\bigr)}{(y/x)^c (z/x)^c}.$$
Factor $(z/x)$ from the $c - 1$ small roots, $(-\tilde{w}_k) = (-w_k)(1 + a_k\varepsilon_s)^n$ from the $n - c - 1$ large roots, and combine the remaining single $(z/x) = R_z w(1+\delta_z)$ in the denominator with the two boundary roots. Because $1 - s e^{\pm i\theta_c} = \frac{r_3}{r_1} e^{\mp i\psi_j}$, we have $w_{\pm\theta_c} = (\frac{r_3}{r_1})^n e^{\mp i j\pi} = (-1)^j R_z$, so $\tilde{w}_{\pm\theta_c} = (-1)^j R_z (1 + \varepsilon_\pm)$ and
$$\frac{((z/x) - \tilde{w}_+)((z/x) - \tilde{w}_-)}{z/x} = R_z (-1)^{j+1} B(w).$$
Finally, $\frac{\prod_{\text{large}} (-w_k)}{t^c} = -\frac{\prod_{\text{large}} |w_k|}{R_y^c}$ for **both** parities of $c$: conjugate pairs contribute $+|w_k|^2 > 0$; for even $c$, the unique real large root $\eta = \pi$ has $-w_\pi = -(1+s)^n < 0$ while $t^c = +R_y^c > 0$; for odd $c$, there is no real large root so $\prod_{\text{large}}(-w_k) > 0$ while $t^c = (-R_y)^c = -R_y^c < 0$. Multiplying $(-1)^{j+1}$ by $-1$ gives $(-1)^j$, and substituting $(1 + \delta_z)^{-1} = (1 + \delta_x)^3 (1 + \delta_y)$ yields the expanded form of $B(w)$. $\square$

*(Independently verified in 80-digit `mpmath` against exact integer coefficients of $Q_c^{(2m)}$ to $< 10^{-50}$ relative error in `algebra/verify_exact_identities_mp.py`.)*

---

## 4. Complete Closed-Form Proof of Regime 1 (All Even $n \ge 84$, $\kappa \ge 3.0$)

Define the separation parameter
$$\kappa(n, c, j) := 2\pi \frac{r_1 r_2}{r_3} = 2\pi \frac{\sin\psi_j \sin\alpha_1}{\sin\theta_c}.$$

### Lemma 4.1 (1-Line Jensen Midpoint Bounds)
Let $a_1 := \min(j + c,\, n - j - c)$, $p := n - 2j - c$, $q := j - c$, $d := n - 3c$. For every strictly interior $(n, c, j)$ ($2 \le c \le (n-3)/3$, $j \in \mathcal{J}_{n,c}$):
1. $0 < c < j < a_1 \le n/2$, $p \ge 1 + \mathbf{1}_{c\text{ even}}$, $q \ge 1$, $a_1 - j = \min(c, p) \ge 1 + \mathbf{1}_{c\text{ even}}$, $a_1 - c = q + \min(c, p) \ge 2 + \mathbf{1}_{c\text{ even}}$, and $d = p + 2q \ge 3$.
2. For any $0 < u < v \le \pi/2$, convexity of $\cot x$ on $(0, \pi/2)$ ($\frac{d^2}{dx^2}\cot x = \frac{2\cos x}{\sin^3 x} > 0$) gives
   $$\Bigl(\frac{\sin u}{\sin v}\Bigr)^n = \exp\!\Bigl(-n\int_u^v \cot x\,dx\Bigr) \le \exp\!\Bigl(-n(v - u)\cot\!\Bigl(\frac{u + v}{2}\Bigr)\Bigr).$$
3. **Universal $(c+j)/2 < n/3$ identity:** Since $3(c + j) = 2(2j + c) - (j - c) \le 2(n - 1) - q < 2n$, we have $\frac{c+j}{2n}\pi < \frac{\pi}{3}$, hence
   $$\frac{R_z}{R_y} = \Bigl(\frac{\sin(c\pi/n)}{\sin(j\pi/n)}\Bigr)^n < \exp\!\Bigl(-\frac{(j - c)\pi}{\sqrt{3}}\Bigr) \qquad \text{for ALL strictly interior } (n, c, j).$$
4. **Universal $R_z$ bounds:** For odd $c$, $c + 1 \le n/3$ so $R_z \le (\sin\theta_c / \sin((c+2)\pi/n))^n \le e^{-2\pi/\sqrt{3}} < 0.02660$. For even $c$ ($3c \le n - 4$ by parity) and $n \ge 84$, $\frac{c + (c+3)}{2n}\pi \le (\frac{1}{3} + \frac{1}{6n})\pi \le (\frac{1}{3} + \frac{1}{504})\pi$, so $R_z \le \exp(-3\pi\cot((\frac{1}{3}+\frac{1}{504})\pi)) < 0.004691$.
5. **Universal $R_y$ midpoint bound:** $\frac{j + a_1}{2} \le \frac{n - c}{2}$, and $\sin\alpha_1 - \sin\psi_j = 2\sin(\frac{c\pi}{2n})\sin(\frac{p\pi}{2n})$.

### Lemma 4.2 (Rouché Isolation & Holomorphic Bounds)
Suppose $\rho_x \in (0, \min(\frac{1}{2}, \frac{1}{2}|1-t|))$ and $u_y \in (0, \frac{1}{2})$ with $u_y R_y < \frac{1}{2}|1-t|$ satisfy the strict inequalities
$$\rho_x(1 - \rho_x)(|1 - t| - \rho_x) > R_y R_z, \qquad u_y(1 - u_y)(|1 - t| - u_y R_y) > \frac{R_z}{R_y}.$$
Then on $|w| \le 1$, the cubic $T(T-1)(T-t) = t R_z w$ has unique holomorphic simple roots $x(w) = 1 + \delta_x(w)$ in $|T - 1| \le \rho_x$ and $y(w) = t(1 + \tilde\delta_y(w))$ in $|T - t| \le u_y R_y$, and
$$|\delta_x| \le \rho_x, \qquad |\delta_y| \le \bar\delta_y := \frac{u_y + \rho_x}{1 - \rho_x}, \qquad |\delta_z| \le \bar\delta_z := \frac{1}{(1 - \rho_x)^3(1 - \bar\delta_y)} - 1, \qquad |\varepsilon_s| \le 1 - (1 - \bar\delta_y)^{1/n} \le \frac{L_{\text{dy}}}{n},$$
where $L_{\text{dy}} := -\log(1 - \bar\delta_y)$ (using $|\binom{1/n}{k}| = (-1)^{k-1}\binom{1/n}{k}$ and $1 - e^{-u} \le u$).

### Lemma 4.3 (4-Line Monotone Riemann Sandwich & Zero-Endpoint $|a_k|$ Bounds)
Let $g(\eta) := \operatorname{Re}\bigl(\frac{-s e^{i\eta}}{1 - s e^{i\eta}}\bigr) = \frac{s^2 - s\cos\eta}{1 + s^2 - 2s\cos\eta}$ and $L_{\text{sum}} := \sum_{\text{large}} a_k = \sum_{\text{large}} g(\eta_k) \in \mathbb{R}$. Then:
1. $g'(\eta) = \frac{s(1-s^2)\sin\eta}{(1+s^2-2s\cos\eta)^2} > 0$ on $(0, \pi)$, $g(\pi) = \frac{s}{1+s} \in (0, \frac{1}{2})$, and $g(\theta_c) = -\frac{r_2\cos\alpha_1}{r_3}$.
2. Since $\frac{d}{d\eta}\arg(1 - s e^{i\eta}) = \operatorname{Im}\frac{-i s e^{i\eta}}{1 - s e^{i\eta}} = g(\eta)$ and $1 - s e^{i\theta_c} = \frac{r_3}{r_1} e^{-i\psi_j}$ ($\psi_j \in (0, \pi/2)$),
   $$\frac{n}{\pi}\int_{\theta_c}^{\pi} g(\eta)\,d\eta = \frac{n}{\pi}\bigl(\arg(1+s) - \arg(1 - s e^{i\theta_c})\bigr) = \frac{n}{\pi}\,\psi_j = j.$$
3. By monotonicity of $g$ on $[\theta_c, \pi]$ with uniform spacing $h = 2\pi/n$ (exact for **both** even and odd $c$):
   $$j - g(\pi) \le L_{\text{sum}} \le j + g(\pi) - 2 g(\theta_c).$$
4. For every large root $|\eta_k| > \theta_c$: $|a_k| \le \min\bigl(\frac{r_2}{r_3},\; \frac{1}{2\sin(|\eta_k|/2)}\bigr)$. Since $\frac{1}{2\sin(\eta/2)}$ and $\frac{1}{4\sin^2(\eta/2)}$ are strictly decreasing on $(0, \pi]$, right-Riemann comparison with zero endpoint correction (both parities of $c$) yields
   $$\frac{1}{n}\sum_{\text{large}} |a_k| \le \frac{1}{\pi}\log\cot\!\Bigl(\frac{\theta_c}{4}\Bigr), \qquad \frac{1}{n}\sum_{\text{large}} |a_k|^2 \le \frac{1}{2\pi}\cot\!\Bigl(\frac{\theta_c}{2}\Bigr).$$

### Lemma 4.4 (Off-Boundary Factor Bounds)
For every $\eta_k \ne \pm\theta_c$, $\frac{|w_k|^2}{R_z^{2/n}} = 1 + \frac{2\kappa}{\pi\sin\theta_c}\sin\bigl(\frac{|\eta_k|+\theta_c}{2}\bigr)\sin\bigl(\frac{|\eta_k|-\theta_c}{2}\bigr)$.
1. **Large roots** ($|\eta_k| = \theta_c + \frac{2r\pi}{n}$, $r \ge 1$): since $\theta_c \le \pi/3$, $\frac{|\eta_k|+\theta_c}{2} \in [\theta_c, \pi - \theta_c]$, so $\sin\bigl(\frac{|\eta_k|+\theta_c}{2}\bigr) \ge \sin\theta_c$ and $\frac{|w_k|^2}{R_z^{2/n}} \ge 1 + \frac{2\kappa}{\pi}\sin\bigl(\frac{r\pi}{n}\bigr)$. Because $(1 + A/n)^{-n/2}$ is strictly decreasing in $n \ge N_0 = 84$ for every $A > 0$,
   $$S_{\text{lg}} := \sum_{\text{large}}\frac{R_z}{|w_k|} \le 2\sum_{r=1}^{2}\!\Bigl(1 + \frac{2\kappa_{\min} r(1 - \frac{r^2\pi^2}{6 N_0^2})}{N_0}\Bigr)^{-N_0/2} + 2\sum_{r=3}^{20}\!\Bigl(1 + \frac{c_J r}{N_0}\Bigr)^{-N_0/2} + \frac{2 N_0}{c_J(N_0/2 - 1)}\!\Bigl(1 + \frac{20 c_J}{N_0}\Bigr)^{-(N_0/2 - 1)},$$
   where $c_J = 4\kappa_{\min}/\pi$. Including perturbations ($|a_k|\le r_2/r_3$, $|\varepsilon_s|\le L_{\text{dy}}/n$, $(1 - \eta/n)^{-n}$ strictly decreasing in $n \ge N_0$):
   $$\sum_{\text{large}}\Bigl|\frac{z/x}{\tilde{w}_k}\Bigr| \le S_{\text{lg}}^{\sim} := S_{\text{lg}}\,(1 + \bar\delta_z)\Bigl(1 - \frac{\eta_{\text{pert}}}{N_0}\Bigr)^{-N_0}, \qquad \eta_{\text{pert}} \ge \frac{r_2}{r_3} L_{\text{dy}}.$$
2. **Small roots** ($|\eta_k| = \theta_c - \frac{2r\pi}{n}$, $1 \le r \le \lfloor c/2 \rfloor$): for $\theta_c \le \theta_{\max}$, $x = \frac{r\pi}{n} \le \frac{\theta_{\max}}{2}$ gives $\frac{\sin(\theta_c - x)}{\sin\theta_c} \ge 1 - \frac{r}{c}$ and $\frac{\sin x}{x} \ge \operatorname{sinc}(\theta_{\max}/2)$, hence
   $$\frac{|w_k|}{R_z} \le \exp\!\Bigl(-\lambda\, r\Bigl(1 - \frac{r}{c}\Bigr)\Bigr), \qquad \lambda := \operatorname{sinc}(\theta_{\max}/2)\,\kappa_{\min}.$$
   Maximizing the finite sum $S_{\text{sm}}$ over $c \in \{2, 3, 4, 5\}$ (each with $\le 2$ values of $r$) and $c \ge \max(c_{\min}, 6)$ (using $r(1-r/c) \ge r/2$ for $r \ge 3$) gives a closed-form upper bound $S_{\text{sm}}$. For perturbations, the triangle inequality $|1 - s(1+\varepsilon_s)e^{i\eta_k}| \le |1 - s e^{i\eta_k}| + s|\varepsilon_s|$ divided by $R_z^{1/n} = r_3/r_1$ gives $\frac{|\tilde{w}_k|^{1/n}}{R_z^{1/n}} \le \rho_k + \frac{r_2}{r_3}|\varepsilon_s| \le \rho_k + \frac{\eta_{\text{pert}}}{n}$ ($\rho_k = |w_k/R_z|^{1/n}$); splitting $\rho_k \ge 1/2$ vs $\rho_k < 1/2$ (at most $c - 1 < n/3$ small roots; $\frac{n}{3}(\frac{1}{2}+\frac{\eta}{n})^n$ strictly decreasing for $n \ge N_0$ when $\eta \le 1$) and using $|(1+\delta_z)^{-1}| \le (1+\rho_x)^3(1+\bar\delta_y)$ yields
   $$\sum_{\text{small}}\Bigl|\frac{\tilde{w}_k}{z/x}\Bigr| \le S_{\text{sm}}^{\sim} := \Bigl(S_{\text{sm}}\, e^{2\eta_{\text{pert}}} + \frac{N_0}{3}\Bigl(\frac{1}{2} + \frac{\eta_{\text{pert}}}{N_0}\Bigr)^{N_0}\Bigr)(1 + \rho_x)^3(1 + \bar\delta_y).$$

### Lemma 4.5 (Exact Cauchy Integral of $B(w)$ & Master Lower Bound)
Under Lemma 4.2, $\delta_x, \delta_y, \delta_z, \varepsilon_\pm$ are holomorphic on $|w| \le 1$ and vanish at $w = 0$. Differentiating $T(T-1)(T-t) = t R_z w$ at $w = 0$ gives $\delta'_x(0) = \frac{t R_z}{1 - t}$, $y'(0) = -\frac{R_z}{1-t}$, hence
$$\delta'_y(0) = -\frac{(1 + t^2) R_z}{t(1 - t)}, \qquad \delta'_z(0) = \frac{(1 - 2t^2) R_z}{t(1 - t)}, \qquad \varepsilon'_+(0) + \varepsilon'_-(0) = 2 g(\theta_c)\,\delta'_y(0).$$
In the expanded form of $B(w)$ (Theorem 3.2), the terms $-Z^{-1}\delta_z$ and $\varepsilon_+ + \varepsilon_-$ vanish at $w = 0$, while $-Z \cdot H(w) = -(-1)^j w^{-1}(1+\varepsilon_+)(1+\varepsilon_-)(1+\delta_x)^3(1+\delta_y)$ has residue $-(-1)^j H'(0) = -(-1)^j(2 g(\theta_c)\delta'_y(0) - \delta'_z(0))$ at $w = 0$. Thus **Cauchy's integral formula gives identically**
$$B_{\text{avg}} := \frac{1}{2\pi}\int_{-\pi}^{\pi} B(e^{i\phi})\,d\phi = 2 - (-1)^j\bigl(2 g(\theta_c)\,\delta'_y(0) - \delta'_z(0)\bigr) \in \mathbb{R}.$$
Furthermore, since $2 - Z - Z^{-1} = 4\sin^2((\phi + j\pi)/2) \ge 0$ on $[-\pi, \pi]$,
$$\frac{1}{2\pi}\int_{-\pi}^{\pi} |B(e^{i\phi})|\,d\phi \le 2 + R_{\text{bdy}}, \qquad R_{\text{bdy}} := \bar\delta_z + 2\varepsilon_{\text{bdy}} + \bigl((1 + \varepsilon_{\text{bdy}})^2(1 + \rho_x)^3(1 + \bar\delta_y) - 1\bigr),$$
with $\varepsilon_{\text{bdy}} := e^{\eta_{\text{pert}}} - 1 \ge |\varepsilon_\pm|$. Finally, writing $F_1 = 1 + E_1$ ($1 + |E_1| \le (1 + \rho_x)^d \le e^{t_1}$, $t_1 \ge d\log(1+\rho_x)$), $F_2 = e^{W_2}$ ($|W_2| \le t_2$), $F_3 = \prod_k(1 + z_k)$ ($1 + |E_3| \le e^{S_{\text{sm}}^\sim + S_{\text{lg}}^\sim}$), we obtain
$$|\Delta_{\text{mult}}(e^{i\phi})| \le M_1 := \exp\bigl(t_1 + t_2 + S_{\text{sm}}^{\sim} + S_{\text{lg}}^{\sim}\bigr) - 1,$$
where
$$t_2 := |L_{\text{sum}} - c|\, L_{\text{dy}} + \frac{(\frac{1}{n}\sum_{\text{large}}|a_k|)\, L_{\text{dy}}^2}{2(1 - L_{\text{dy}}/n)} + \frac{(\frac{1}{n}\sum_{\text{large}}|a_k|^2)\, L_{\text{dy}}^2}{2(1 - \eta_{\text{pert}}/n)},$$
and therefore the **Master Lower Bound** holds:
$$\frac{1}{2\pi}\int_{-\pi}^{\pi}\operatorname{Re}\bigl((1 + \Delta_{\text{mult}}(e^{i\phi}))\, B(e^{i\phi})\bigr)\,d\phi \;\ge\; B_{\text{avg}} - M_1\,(2 + R_{\text{bdy}}).$$

### Theorem 4.6 (Closed-Form Certification of Regime 1 for ALL Even $n \ge 84$)
For every even $n = 2m \ge 84$ and every strictly interior $(c, j)$ ($2 \le c \le \lfloor(n-3)/3\rfloor$, $j \in \mathcal{J}_{n,c}$) satisfying $\kappa(n, c, j) \ge 3.0$ (and indeed even $\kappa(n, c, j) \ge 2.5$), we have
$$B_{\text{avg}} - M_1\,(2 + R_{\text{bdy}}) > 0, \qquad \text{and hence } (-1)^j Q_c^{(2m)}(1, u_j) > 0.$$

*Proof.* Partition Regime 1 into 6 mutually exhaustive closed-form subcases and apply Lemmas 4.1–4.5 (every elementary formula is evaluated and cross-checked in [`verify_regime1_complete.py`](algebra/verify_regime1_complete.py)):

1. **Case 1A ($2 \le c \le 0.15 n$, $\kappa \ge \kappa_0$):**
   - If $R_y > 1/2$: then $s^{-1} = \frac{\sin(\psi_j + \theta_c)}{\sin\psi_j} = \cos\theta_c + \cot\psi_j\sin\theta_c < 2^{1/n}$. Since $\frac{2^{1/n}-1}{1/n} \le 84(2^{1/84}-1)$ and $\sin\theta_c \ge \sin\frac{2\pi}{n} \ge \frac{84\sin(2\pi/84)}{n}$ for all $n \ge 84$ and $c \ge 2$, the factors of $1/n$ cancel identically:
     $$\cot\psi_j < \frac{2^{1/n}-1}{\sin\theta_c} + \tan\frac{\theta_c}{2} \le \frac{2^{1/84}-1}{\sin(2\pi/84)} + \tan(0.075\pi) < 0.3510 \;\implies\; \sin\psi_j > 0.9435, \quad \frac{r_3}{r_2} < 0.48114,$$
     hence $\frac{R_z}{R_y} < 0.48114^n \le 2.045 \times 10^{-27}$ (and $n^2 \frac{R_z}{R_y} \le 1.45 \times 10^{-23}$ for all $n \ge 84$). Meanwhile, for odd $c$ we have $|1 - t| = 1 + R_y > 1 > \frac{4}{n}$; for even $c$, $c \ge 2$ and $p \ge 2$, so Lemma 4.1(5) gives $1 - s \ge 2\sin^2(\frac{\pi}{n}) \ge \frac{8}{n^2} \implies |1 - t| = 1 - s^n \ge 1 - e^{-8/n} > \frac{4}{n}$. Thus Rouché (Lemma 4.2) holds with $\rho_x = u_y = 10^{-20}/n$ (since $\frac{10^{-20}}{n^2} \gg 1.45\times 10^{-23}/n^2 \ge \frac{R_z}{R_y}$), making all perturbation terms ($t_1, t_2, \eta_{\text{pert}}, R_{\text{bdy}}, |B_{\text{avg}}-2|$) $< 10^{-19}$ and strictly dominated by Subcase 1A-i.
   - If $R_y \le 1/2$ (Subcase 1A-i): $|1-t| \ge 1/2$, $\frac{R_z}{R_y} \le \exp(-\pi\cot(0.15\pi + \frac{\pi}{168})) < 0.002411$. Also $\alpha_1 < 0.575\pi$, so either $\sin\alpha_1 \ge 0.70$ ($R_z \le 0.6486^n$) or $\alpha_1 < \arcsin(0.70) < 0.247\pi$ (whence $\kappa \ge \kappa_0$ forces $j \ge c_j\sqrt{n}$ with $c_j = \sqrt{\frac{\operatorname{sinc}(0.15\pi)\kappa_0}{2\pi^2}}$ and $R_z \le e^{-\pi\cot(\arcsin 0.70)\,j}$, giving $n R_z$ strictly decreasing in $n \ge 84$). Using $\frac{r_2}{r_3}\frac{R_z}{R_y} = (\frac{R_z}{R_y})^{1 - 1/n}$ and splitting $q = j-c = 1$ vs $q \ge 2$ ($j \le 0.25n$ via $\cot(0.20\pi)$ vs $j > 0.25n$ via $(\frac{\sin 0.15\pi}{\sin 0.25\pi})^n$), plus the trivial count $\#\text{large} \le n$ with $|a_k| \le r_2/r_3$ in $t_2$, yields **Margin $= +0.8540$** at $\kappa_0 = 3.0$ (and **`+0.0748`** at $\kappa_0 = 2.5$).

2. **Case 1B-I ($c > 0.15 n \implies c \ge 13$, $j \ge c + 2$):**
   Here $\kappa \ge 2\pi\sin(0.30\pi) = 5.0832$, $\frac{R_z}{R_y} \le e^{-2\pi/\sqrt{3}} = 0.02660$ (Lemma 4.1(3)), $R_y \le e^{-2\pi\cot(0.425\pi)} = 0.22125$ for even $c$ ($R_y \le \sqrt{0.22125} = 0.4704$ for odd $c$, so $R_y^2 \le 0.22125$ for both parities!), and $d = p + 2q \ge 5$ (odd) or $\ge 6$ (even) yields $R_y R_z \le e^{-\frac{\pi}{2\sqrt{3}} d}$ (when $a_1 = n-j-c$, since $2\pi\cot(0.425\pi) > \frac{\pi}{2\sqrt{3}}$; and $< 10^{-5}$ when $a_1 = j+c$). Since $\alpha_1 \in (0.30\pi, 2\pi/3)$, $|g(\theta_c)| \le \frac{\cos(0.30\pi)}{\sin(0.15\pi)}$ and $(q + 0.5 + 2|g(\theta_c)|) L_{\text{dy}}$ is controlled at $q = 2$ and $q \ge 3$ via $e^{-q\pi/\sqrt{3}}$, giving **Margin $= +0.8435$**.

3. **Case 1B-IIa ($0.15 n < c \le 0.25 n \implies c \ge 13$, $j = c + 1$):**
   Here $\kappa \ge 5.0832$, $\frac{R_z}{R_y} \le e^{-\pi\cot(0.25\pi + \frac{\pi}{168})} = 0.04579$, $a_1 - j > 0.15 n \implies R_y \le e^{-0.15 \cdot 84 \cdot \pi\cot(0.425\pi)} < 7.5\times 10^{-5}$, $\frac{r_2}{r_3} \le \cos(\frac{\pi}{84}) + \cot(0.15\pi)\sin(\frac{\pi}{84}) = 1.073$, giving **Margin $= +1.3934$**.

4. **Case 1B-IIb ($0.25 n < c \le 0.31 n \implies c \ge 22$, $j = c + 1$):**
   Here $\alpha_1 = \frac{(2c+1)\pi}{n} \in (\frac{\pi}{2}, 0.62\pi + \frac{\pi}{84}) \subset (\frac{\pi}{2}, \frac{2\pi}{3})$, so $\kappa > \pi\sqrt{3} = 5.4414$ and **$\cos\alpha_1 < 0 \implies g(\theta_c) > 0 \implies 0.5 < L_{\text{sum}} - c < 1.5$**! Also $\frac{R_z}{R_y} \le e^{-\pi\cot(0.31\pi+\frac{\pi}{168})} = 0.1160$ and $a_1 - j \ge 0.07n - 2$ with midpoint $< 0.375n$ gives $R_y \le e^{-(0.07\cdot 84 - 2)\pi(\sqrt{2}-1)} = 0.00642$, giving **Margin $= +0.6359$**.

5. **Case 1B-III even ($0.31 n < c \le (n-3)/3 \implies c \ge 27$, $j = c + 1$, even $c$):**
   Here $\kappa > \pi\sqrt{3}$, $\alpha_1 \in (0.62\pi, \frac{2\pi}{3}) \implies g(\theta_c) \in (-\cos(0.62\pi), \frac{1}{2}\frac{r_2}{r_3}) \implies 0.5 < L_{\text{sum}} - c \le 1.5 + 2\cos(0.62\pi) = 0.7638$. By parity $d = n - 3c \ge 4$, $a_1 - j = d - 2 \ge 2$, midpoint $\frac{n-c}{2} < 0.345 n$, giving $R_y \le e^{-2\pi\cot(0.345\pi)} = 0.03542$, $R_z < 0.004691$ (Lemma 4.1(4)), $d R_y R_z \le 4 e^{-5\pi\cot(0.345\pi)} < 0.00095$, yielding **Margin $= +0.5017$**.

6. **Case 1B-III odd ($0.31 n < c \le (n-3)/3 \implies c \ge 27$, $j = c + 1$, odd $c$):**
   Same $\kappa, L_{\text{sum}} - c \le 0.7638$, with $t = -R_y < 0 \implies |1 - t| = 1 + R_y \ge 1$, $|1 - t| - u_y R_y = 1 + (1 - u_y) R_y > 1$ (justifying $\text{omt}_{\min} = 1$ in Lemma 4.2), and $\frac{1+t^2}{1-t} \le 1$, $\frac{|1-2t^2|}{1-t} \le 1$ in $B_{\text{avg}}$. Here $d \ge 3 \implies R_y \le e^{-\pi\cot(0.345\pi)} = 0.1882$, $R_z \le e^{-2\pi/\sqrt{3}} = 0.02660$, $d R_y R_z \le 3 e^{-3\pi\cot(0.345\pi)} < 0.0200$, yielding **Margin $= +0.4699$**.

| Subcase ($\kappa_0 = 3.0$, $n \ge 84$) | $S_{\text{off}}^\sim$ | $t_1$ | $t_2$ | $M_1$ | $R_{\text{bdy}}$ | $B_{\min}$ | **Closed-Form Margin** |
|---|---:|---:|---:|---:|---:|---:|---:|
| **1A** ($c \le 0.15n$, $\kappa \ge 3.0$) | 0.4125 | 0.0011 | 0.0255 | 0.5513 | 0.0416 | 1.9796 | **`+0.8540`** |
| **1B-I** ($c > 0.15n$, $j \ge c+2$) | 0.0523 | 0.0545 | 0.2292 | 0.3993 | 0.5785 | 1.8732 | **`+0.8435`** |
| **1B-IIa** ($0.15n < c \le 0.25n$, $j=c+1$) | 0.0454 | 0.0000 | 0.1467 | 0.2118 | 0.3460 | 1.8903 | **`+1.3934`** |
| **1B-IIb** ($0.25n < c \le 0.31n$, $j=c+1$) | 0.0410 | 0.0001 | 0.2617 | 0.3537 | 1.1844 | 1.7623 | **`+0.6359`** |
| **1B-III even** ($c > 0.31n$, $j=c+1$, even) | 0.0495 | 0.0010 | 0.2107 | 0.2984 | 1.8726 | 1.6574 | **`+0.5017`** |
| **1B-III odd** ($c > 0.31n$, $j=c+1$, odd) | 0.0492 | 0.0206 | 0.2023 | 0.3126 | 1.8382 | 1.6699 | **`+0.4699`** |

*(At $\kappa_0 = 2.5$, Case 1A becomes $S_{\text{off}}^\sim = 0.6306, M_1 = 0.9329$, Margin $= +0.0748 > 0$, while Cases 1B-I through 1B-III are unchanged.)* $\square$

---

## 5. Complete Regime 2 Proofs: Piece G2 ($D = 0$ Major Arcs) & Piece G1 (Dirichlet Minor Arcs)

Throughout Section 5, fix even $n = 2m \ge 132$, strictly interior $2 \le c \le \lfloor (n-3)/3 \rfloor$, $j \in \mathcal{J}_{n,c}$ ($j \ge c + 1$), and $\kappa(n, c, j) = 2\pi \frac{r_1 r_2}{r_3} < 2.5$ with $c\kappa^2 < 42$. Set
$$d := \gcd(c, j) \ge 1, \qquad q := \frac{c}{d} \ge 1, \qquad p_0 := \frac{j}{d} \quad (\gcd(p_0, q) = 1), \qquad X := q^2\kappa, \qquad T := \frac{1}{2} c q \kappa = \frac{1}{2} d X.$$
By `@analyst`'s Universal Cubic Perturbation Reduction (R2.1, certified for all $\lambda \le 3$ inside $n \ge 132, \kappa < 2.5, c\kappa^2 < 42$), positivity $(-1)^j Q_c^{(2m)}(1, u_j) > 0$ reduces rigorously to exhibiting a contour radius $\lambda \in (0, 3]$ such that the exact unperturbed circle integral $I_{\text{pure}}(\lambda) := \frac{1}{2\pi}\int_{-\pi}^{\pi} \operatorname{Re}\bigl(B_0(\lambda e^{i\phi}) F_3^{(0)}(\lambda e^{i\phi})\bigr)\,d\phi$ satisfies
$$\frac{I_{\text{pure}}(\lambda)}{L^1(\lambda)} \;\ge\; 0.01 \;\gg\; 10^{-5}, \qquad L^1(\lambda) := \frac{1}{2\pi}\int_{-\pi}^{\pi}\bigl|B_0(\lambda e^{i\phi}) F_3^{(0)}(\lambda e^{i\phi})\bigr|\,d\phi.$$

Indexing the $n$ roots by $x_k = k + \text{shift} \in (-n/2, n/2]$ ($\text{shift} = 0$ for even $c$, $\frac{1}{2}$ for odd $c$, so $x_k - c/2 \in \mathbb{Z}$ for **both** parities of $c$), define $\eta(x) := \frac{2\pi x}{n}$, $u_k := (-1)^j R_z^{-1} (1 - s e^{i\eta(x_k)})^n$, and $a_k := u_k^{-1} = b(x_k)\, e^{2\pi i p_0 (x_k - c/2)/q}$, where
$$b(x) := \Bigl(\frac{r_3 / r_1}{1 - s e^{i\eta(x)}}\Bigr)^n e^{-2\pi i \frac{j}{c}(x - c/2)} (-1)^j = \exp\bigl(\log|b(x)| + i\,\delta\psi(x)\bigr), \qquad b(\pm c/2) = 1.$$

---

### 5.1 Piece G2: Complete 50-Digit Certified Proof of $D = 0$ Major Arcs

**Domain (Box G2):** $n \ge 132$, $\kappa < 2.5$, $c\kappa^2 < 42$, $d = \gcd(c, j) \ge 2$, $q = c/d \ge 1$, $X = q^2\kappa \le K_{\text{G2}, q=1} := 2.5$ (for $q = 1$, covering ALL of Regime 2 for $q = 1$ with zero crevice) and $q^2\kappa \le K_{\text{G2}} := 2.0$ (for $q \ge 2$), and $c q \kappa = d q^2\kappa \ge 0.85$ (so $T = \frac{1}{2} d q^2\kappa \ge 0.425$, strictly overlapping R2.3's $c q \kappa \le 0.90$ and Piece G1's $q \ge 2, q^2\kappa \ge 2.0$).

#### Theorem 5.1 (Exact Universal $D = 0$ Contour Identity — All Parities of $(d, q, p_0)$)
Set $z_{\text{var}} := Z^{-1} = (-1)^j w$ along $|z_{\text{var}}| = \lambda_*$. Partition the circle $|z_{\text{var}}| = \lambda_*$ into the $q$ equal sectors of width $2\pi/q$ centered at the $q$ roots of $z_{\text{var}}^q = -\lambda_*^q$:
$$z_{\text{var}} = \lambda_* \exp\!\Bigl(i\,\frac{(2m + 1)\pi + \psi}{q}\Bigr) \qquad (m = 0, 1, \dots, q - 1,\quad \psi \in [-\pi, \pi]).$$
Then for **every** even $n$, every $c = d q$, and every $j = d p_0$ with $\gcd(p_0, q) = 1$ (regardless of the parities of $d, q, p_0$):
$$B_0(Z)\, F_3^{(0)}(Z) \;=\; \Bigl(\prod_{\text{small}} |u_k|\Bigr)\,\lambda_*^{-c}\, e^{-i d\psi}\,\prod_{k=0}^{n-1}\bigl(1 - a_k z_{\text{var}}\bigr).$$
Moreover, within each sector $m$, the $q$ residue classes of $x_k - c/2 \pmod q$ permute bijectively onto the $q$ exact rays $\theta_t(\psi) \equiv \frac{(2t + 1)\pi + \psi}{q} \pmod{2\pi}$ ($t = 0, \dots, q-1$), whose product over all $q$ rays is identically $\prod_{t=0}^{q-1}(1 - b\lambda_* e^{i\theta_t(\psi)}) = 1 + (b\lambda_*)^q e^{i\psi}$.

*Proof.* Rewriting every factor of $B_0(Z) F_3^{(0)}(Z) = (-Z^{-1})(1 - Z)^2 \prod_{\text{small}}(1 - u_k Z)\prod_{\text{large}}(1 - u_k^{-1} Z^{-1})$ in terms of $z_{\text{var}} = Z^{-1}$ via $1 - u_k Z = -u_k z_{\text{var}}^{-1}(1 - a_k z_{\text{var}})$ and $(1 - Z)^2 = z_{\text{var}}^{-2}(1 - z_{\text{var}})^2$ produces the exact prefactor $(-1)^c(\prod_{\text{small}} u_k) z_{\text{var}}^{-c}$.
- If $q$ is **even**: $\gcd(p_0, q) = 1$ forces $p_0$ odd and $c = d q$ even, so $j = d p_0 \equiv d \pmod 2$. The unique unpaired small root $x = 0$ has $u_0 = (-1)^j |u_0|$, giving $(-1)^c \prod_{\text{small}} u_k = (-1)^j \prod |u_k| = (-1)^d \prod |u_k|$. Meanwhile $z_{\text{var}}^{-c} = \lambda_*^{-c} e^{-i d((2m+1)\pi + \psi)} = (-1)^d \lambda_*^{-c} e^{-i d\psi}$. The $(-1)^d$ factors cancel identically.
- If $q$ is **odd**: for even $c$ ($d$ even), $j$ is even and $(-1)^c (-1)^j = +1$, while $z_{\text{var}}^{-c} = (-1)^d \lambda_*^{-c} e^{-i d\psi} = +\lambda_*^{-c} e^{-i d\psi}$; for odd $c$ ($d$ odd), all $c - 1$ small roots pair up in conjugates so $(-1)^c \prod_{\text{small}} u_k = -\prod |u_k|$, while $z_{\text{var}}^{-c} = -\lambda_*^{-c} e^{-i d\psi}$. Again the signs cancel identically.

Since $x_k - c/2 \in \mathbb{Z}$ and $\gcd(p_0, q) = 1$, $p_0(x_k - c/2) \pmod q$ cycles through $\{0, \dots, q-1\}$, yielding the stated ray identity. $\square$

#### Lemma G2.1 (Universal Regime 2 Modulus, $Q_{\min} > 0.20$, Curvature, and Phase Laws)
For all $(n, c, j)$ in Regime 2 ($n \ge 132$, $\kappa < 2.5$, $j \ge c + 1$):
1. **Angles, Explicit $Q_{\min} > 0.20$ Cutoffs, & Dead Tail:** $\sin\alpha_1 = \frac{\kappa}{2\pi}\frac{r_3}{r_2} < \frac{2.5}{2\pi} \implies \alpha_1 < 0.13032\pi$, $\theta_c < \frac{\alpha_1}{2} < 0.06516\pi$, and $u_{Q_0} := \frac{\kappa}{\pi}\tan\frac{\theta_c}{2} < 0.0818$. Define $Q(x) := \frac{|1-se^{i\eta(x)}|^2}{(r_3/r_1)^2} - 1$. Using $Q(\rho n) \ge 8\cos\theta_c\sin(\rho\pi + \frac{\theta_c}{2})\sin(\rho\pi - \frac{\theta_c}{2})$ (and the sharper factor $4(3-4\sin^2\theta_c)(2\cos\theta_c)$ when $q = 1, j \ge 2c$), 50-digit `mpmath.iv` explicitly certifies $Q(\rho n) > Q_{\text{cut}} := 0.20$ at:
   - **All Regime 2 ($\rho = 0.100$):** $Q(0.100 n) \ge 0.2644 > 0.20$, $-\log|b(x)| \ge 15.484$, and $\sum_{Q(x) \ge 0.20} |b(x)|^q \le \text{tail\_sum\_Q20} := 132 \cdot 1.20^{-66} < 0.000785$.
   - **Box G2, $q = 1$ ($\rho = 0.050$, $j \ge 2c$, $\kappa < 2.5$):** $Q(0.050 n) \ge 0.5548 > 0.20$.
   - **Box G2, $q \ge 2$ ($\rho = 0.060$, $\kappa \le K_{\text{G2}}/4 = 0.50$):** $Q(0.060 n) \ge 0.2775 > 0.20$.
2. **Two-Sided Modulus Sandwich & Exact Global UHP Curvature Bound:** Writing $v_G(x) := \frac{\kappa}{c}|\frac{c^2}{4} - x^2|$, we have $\log|b(x)| = c_{\text{in}}(x)\, v_G(x)$ on $[0, c/2]$ and $-\log|b(x)| = c_{\text{out}}(x)\, v_G(x)$ on $[c/2, \rho n] \cap \{Q(x) < 0.20\}$. Moreover, for $z = x + i y$ ($y \ge 0$), $\xi = \frac{2\pi x}{n}$, $\zeta = \frac{2\pi y}{n} \ge 0$, and $\beta(z) := \frac{b'(z)}{b(z)}$, since $r \mapsto r/(1-r)^2$ is strictly increasing on $[0, 1)$,
   $$|\beta'(x + i y)| = \frac{4\pi^2}{n}\frac{s e^{-\zeta}}{|1 - s e^{-\zeta} e^{i\xi}|^2} \le \frac{4\pi^2 s}{n(1 - s)^2} \le C_{u''}\frac{2\kappa}{c} \qquad \Bigl(C_{u''} := \frac{1}{\operatorname{sinc}(\theta_c/2)\cos^2\alpha_1}\Bigr)$$
   **globally for all $x \in \mathbb{R}, y \ge 0$** (in particular $|u''(x)| \le C_{u''}\frac{2 q \kappa}{c}$ on $[0, n/2]$). Certified 50-digit `iv` constants:
   - **All Regime 2 ($\kappa < 2.5$):** $c_{\text{in}} \in [0.99130, 1.05056]$, $c_{\text{out}} \in [0.85699, 1.00702]$.
   - **Box G2, $q = 1$ ($\alpha_1 < 0.065\pi$, $\theta_c < 0.02167\pi$):** $c_{\text{in}} \in [0.999035, 1.01458]$, $c_{\text{out}} \in [0.903784, 1.00077]$, $C_{u'', q1} \le 1.04309$.
   - **Box G2, $q \ge 2$ ($\alpha_1 < 0.02537\pi$, $\theta_c < 0.01269\pi$):** $c_{\text{in}} \in [0.999669, 1.00185]$, $c_{\text{out}} \in [0.900742, 1.00026]$, $C_{u'', q2} \le 1.00645$.
3. **Exact Phase Curvature & Green's Function Bounds:** Since $\phi(\eta) := \arg(1 - s e^{i\eta})$ has $\phi''(\eta) = \frac{s(1-s^2)\sin\eta}{|1-se^{i\eta}|^4} > 0$ on $(0, \pi)$ and $1 - s^2 = \frac{r_3\sin(\alpha_1+\psi_j)}{r_1^2}$, the deviation $\delta\psi(x) = -n\phi(\eta(x)) - \frac{2\pi j}{c} x$ vanishes at $x = 0, c/2$ and satisfies
   $$0 \le -\delta\psi''(x) \le A_{\text{in}}\frac{\kappa^2}{c^2}\, x \quad (x \in [0, c/2]), \qquad 0 \le -\delta\psi''(x) \le A_{\text{out}}\frac{\kappa^2}{c^2}\, x \quad (x \in [c/2, n/2]),$$
   whence $0 \le \delta\psi(x) \le \frac{A_{\text{in}}}{6}\frac{\kappa^2}{c^2}\, x(\frac{c^2}{4} - x^2) \le \rho_{\text{in}}\,\kappa\log|b(x)|$ on $[0, c/2]$, $|\delta\psi'(c/2)| \le C_\psi\,\kappa^2$ ($C_\psi = A_{\text{in}}/12$), and $|\delta\psi(c/2 + m)| \le C_\psi\kappa^2 m + A_{\text{out}}\frac{\kappa^2}{c^2}(\frac{c m^2}{4} + \frac{m^3}{6})$ for $m > 0$, with certified constants:
   - Box G2, $q = 1$: $C_{\psi, q1} \le 0.0473395$.
   - Box G2, $q \ge 2$: $C_{\psi, q2} \le 0.080404$, $A_{\text{out}, q2} \le 0.9600$.

#### Lemma G2.2 (Saddle $\Lambda_* = \lambda_*^q$ & Curvature $H_q^+$ via Exact 2nd-Order Euler–Maclaurin)
Define the pure positive-root polynomial $G_q^+(W) := \prod_{k=0}^{n-1}(1 + |b(x_k)|^q W)$ and let $\Lambda_* = \lambda_*^q > 0$ solve $M_q^+(\Lambda_*) := \frac{1}{q}\sum_{k=0}^{n-1}\sigma(u_k) = d$ with curvature $H_q^+(\Lambda) := \frac{1}{q}\sum_{k=0}^{n-1} h(u_k)$, where $u_k := q\log(\lambda|b(x_k)|)$, $\sigma(u) = (1+e^{-u})^{-1}$, $h(u) = e^u/(1+e^u)^2$.
For both $g \in \{\sigma, h\}$, the 2nd-order Euler–Maclaurin remainder ($\max|B_2(\{x\})/2| = 1/12$) combined with Lemma G2.1's global bound $|u''(x)| \le C_{u''}\frac{2 q\kappa}{c}$ yields the **zero-hardcoded** 50-digit `iv` discretization error
$$\Bigl|\frac{1}{q}\sum_k g(u_k) - \frac{2}{q}\int_0^{\rho n} g(u(x))\,dx\Bigr| \le \frac{C_{u''}}{6}\kappa\Bigl(\|g''\|_{L^1(\mathbb{R})} + \sqrt{\frac{2}{c_{\text{out,lo}} T}}\, J_{\sqrt{}}(g; \log\Lambda) + \int_0^{2\rho n/c} |g'(u(t))|\,dt\Bigr) + \Lambda\cdot\text{tail\_sum\_Q20},$$
where $\|\sigma''\|_{L^1} = \frac{1}{2}$, $\|h''\|_{L^1} = \frac{2}{3\sqrt{3}}$, and $J_{\sqrt{}}(g; \ell) := \int_0^\infty \sqrt{s}\,|g''(\ell - s)|\,ds$ is explicitly certified in 50-digit `iv` over $\ell \in [-0.75, 0.90]$ as $J_{\sqrt{}}(\sigma) \le 0.49969$, $J_{\sqrt{}}(h) \le 0.41158$. Evaluating direct discrete sums for $q = 1, c \in \{2,\dots,10\}$ (which covers every $q = 1$ point with $\kappa \ge \sqrt{42/11} \approx 1.9540$) and 30 narrow $T$-bins on $[0.425, 14.0]$ plus an exact `iv`-certified FTC tail for $T \ge 14.0$ (with explicit $\Lambda\cdot\text{tail\_sum\_Q20}$ everywhere):
- **Subcase 2a ($q = 1, c = d \in \{2, \dots, 10\}$):** $H_1^+(\Lambda_*) \ge 0.68924$ ($c=2: 0.68924$, $c=3: 0.77078$, $c=4: 0.80396$, $c=5: 0.81388$, $c=6: 0.81356$, $c=7: 0.82480$, $c=8: 0.87341$, $c=9: 0.91971$, $c=10: 0.96386$).
- **Subcase 2b ($q = 1, c = d \ge 11 \implies \kappa < \sqrt{42/11} < 1.9541 < 2.0$):** $\Lambda_* \in [0.5717, 2.162]$ (tail $[0.645, 1.755]$), $H_{\text{cont,end}} \ge 0.84276$, $\max\text{disc}_H \le 0.24292 \implies H_1^+(\Lambda_*) \ge 0.62020$.
- **Subcase 2c ($q \ge 2, d \ge 2$):** $\Lambda_* \in [0.5519, 1.803]$ (tail $[0.7735, 1.492]$), $H_{\text{cont,end}} \ge 0.85844$, $\max\text{disc}_H \le 0.095147 \implies H_q^+(\Lambda_*) \ge 0.59428$.

To accommodate the Pigeonhole-Annulus shift $\tilde{\lambda} \in [\lambda_*(1 - \varepsilon_{\text{rad}}/n), \lambda_*(1 + \varepsilon_{\text{rad}}/n)]$ ($\varepsilon_{\text{rad}} = 10^{-4}$, Lemma G2.7), we pre-widen every $\Lambda$-interval by $e^{\pm\delta_{\text{rad}}}$ ($\delta_{\text{rad}} := 0.35\,\varepsilon_{\text{rad}} \ge |\log(\tilde{\Lambda}/\Lambda_*)|$ since $q \le n/3$), so Lemmas G2.3–G2.5 hold **literally at $\tilde{\Lambda}$**, while $|h'/h| = |\tanh(u/2)| < 1$ gives $H_q^+(\tilde{\Lambda}) \ge e^{-\delta_{\text{rad}}} H_q^+(\Lambda_*)$ ($H_1^+(\tilde{\Lambda}) \ge 0.62018$, $H_q^+(\tilde{\Lambda}) \ge 0.59425$, $H_{\text{cont,end}}(\tilde{\Lambda}) \ge 0.85841$).

#### Lemma G2.3 (Exact $P^2$-Monotone Conjugate-Pair Phase Factorization $E_{\text{phase}}$ at $\tilde{\Lambda}$)
Set $G_q(W) := \prod_{k=0}^{n-1}(1 + b(x_k)^q W)$. Pairing conjugate roots $(x_k, -x_k)$ with $b(-x_k) = \overline{b(x_k)}$ gives
$$\frac{G_q(\tilde{\Lambda} e^{i\psi})}{G_q^+(\tilde{\Lambda} e^{i\psi})} = \prod_{\text{pairs } k} \left(1 - \frac{4\sin^2(q\,\delta\psi_k/2)}{\bigl(e^{u_k/2 + i\psi/2} + e^{-u_k/2 - i\psi/2}\bigr)^2}\right) =: E_{\text{phase}}(\tilde{\Lambda} e^{i\psi}), \qquad |E_{\text{pair}, k} - 1| = \frac{\sin^2(q\,\delta\psi_k/2)}{\sinh^2(u_k/2) + \cos^2(\psi/2)}.$$
Using $\sin^2(P) \le P^2$ throughout (preserving exact linear $X = q^2\kappa$ monotonicity on every bin and outer tail $w \ge 20$ as $d \to \infty, X \to 0$), exact $T \ge 22$ inner split $[0, 8.8] \cup [8.8, 0.40 T] \cup \{t \in [0, \sqrt{0.20}]\}$ (strictly decreasing in $T \ge 22$), and explicit `iv` certification of $\Delta_{\max} < 1$:
- **$q \ge 2$ ($|\psi| \le 2\pi/3$):** $\Delta_{\max} \le 0.011303$ and $\bigl|\frac{1}{q}\log E_{\text{phase}}(\tilde{\Lambda} e^{i\psi})\bigr| \le 0.0210522$.
- **$q = 1$ ($|\psi| \le 1.75$):** $\Delta_{\max} \le 0.039652$ and $|\log E_{\text{phase}}(\tilde{\Lambda} e^{i\psi})| \le 0.121256$. Moreover, by AM–GM on conjugate pairs, $|(1 + u e^{i(\psi+d)})(1 + u e^{i(\psi-d)})|/(1+u)^2 \le 1 - 2 h(u)(1 - \cos\psi\cos d)$. For $|\psi| \in [1.75, 2\pi/3]$ this gives $\varepsilon_{\text{mid}} \le \frac{1}{2} S_{\text{AMGM}} \le 0.0443721$; for $|\psi| \in [2\pi/3, \pi]$ ($\cos\psi \in [-1, -1/2]$), the exact identity $1 - \cos\psi\cos d = 1 - \cos\psi + 2\cos\psi\sin^2(d/2) \ge \frac{3}{2} - 2\sin^2(d/2)$ holds for **all** $d \in \mathbb{R}$, giving $\log(\dots) \le -3 h + 4 h \sin^2(d/2) \le -3 h + \sin^2(d/2)$, i.e., decay rate $c_{\text{far}} = 1.50$ with $\varepsilon_{\text{far}} \le S_{\text{AMGM}} \le 0.0887442$.

#### Lemma G2.4 (100% Pole-Certified $\pi$-Line Meromorphic Residue & Cauchy–Riemann Level-Curve Bound on $\mathcal{R}_{\text{non-}q}$, $q \ge 2$)
Write $\log\prod_{k=0}^{n-1}(1 - a_k z_{\text{var}}) = \frac{1}{q}\log G_q(\tilde{\Lambda} e^{i\psi}) + \mathcal{R}_{\text{non-}q}(z_{\text{var}})$ and decompose $\mathcal{R}_{\text{non-}q}$ over the $q$ rays $\theta_t(\psi) = \frac{(2t+1)\pi + \psi}{q}$ via step-$q$ Poisson summation. Integrating each Fourier mode $\ell \ne 0$ by parts once eliminates all branch cuts in favor of the **meromorphic** integrand $\frac{\beta(z) u(z)}{1 - u(z)} e^{2\pi i \ell z / q}$ ($u(z) = \tilde{\lambda} b(z) e^{i\theta_t}$), whose residue at any simple zero $z_0 = x_0 + i y_0$ of $1 - u(z)$ is identically $-e^{2\pi i \ell z_0 / q}$, summing over $\ell \ge 1$ (plus the symmetric LHP pole $-z_0$) to at most $-2\log(1 - e^{-2\pi y_0/q})$.
1. **Exact 2-Line Cauchy–Riemann Level-Curve Theorem:** Writing $\log b(x + i y) = U(x, y) + i V(x, y)$ ($U$ even in $x$, $V$ odd in $x$), we have $U_x(x, y) = \frac{-2\pi s e^{-\zeta}\sin\xi}{|1-se^{-\zeta}e^{i\xi}|^2} < 0$ strictly on $(0, n/2) \times [0, \infty)$. Hence along the unique level curve $U(x_*(y), y) \equiv -\log\tilde{\lambda}$ ($x_* > 0$),
   $$\frac{d}{dy}\,V(x_*(y), y) = \frac{U_x^2 + U_y^2}{U_x} \;<\; 0 \qquad \text{strictly for all } y \ge 0.$$
   Moreover, when $U(0, 0) \le -\log\tilde{\lambda}$ (e.g. $T = 0.425, \tilde{\Lambda} < 1$), the level curve starts on the imaginary axis $(0, y_{\text{start}})$ where $V(0, y_{\text{start}}) \equiv 0$; when $U(0, 0) > -\log\tilde{\lambda}$, it starts at $(x_*(0), 0)$ (inner when $\tilde{\Lambda} < 1$, outer when $\tilde{\Lambda} > 1$) and 50-digit `iv` certifies $q|\delta\psi(x_*(0))| < 0.15 < \pi/3$. Consequently, for every ray $\theta_t \in (0, \pi]$ in $(-L_0, L_0) \times (0, H]$ ($L_0 := 0.25 n$):
   - **Case 1 ($m = 0$ rays on $|\psi| \le 2\pi/3$, where $\theta_t \le \frac{5\pi}{6}$, and all higher rays with $\theta_t \le 0.72\pi$):** Shifting to the $\pi$-line $H = H_\pi$ where $-V(x_*(H_\pi), H_\pi) = \pi$ encloses **exact zero count $= 1$** ($z_0 = x_0 + i y_0$ on $x > 0$; zero on $x < 0$ since $\theta_t + \pi \le \frac{11\pi}{6} < 2\pi$), with certified crossing phase distance $\delta_H = \pi - \theta_t \ge 0.28\pi$ ($|\psi| \le 1.35$) and $\ge \pi/6$ ($|\psi| \le 2\pi/3$).
   - **Case 2 (higher rays $m \ge 1$ with $\theta_t \in (0.72\pi, \pi]$):** Shifting to $H = H_{3\theta_t/4}$ where $-V(x_*(H), H) = \frac{3}{4}\theta_t$ encloses **exact zero count $= 0$**, with crossing phase distance $\delta_H \ge 0.18\pi$.
2. **Certified Lower Bound $E_{p,\text{lb}}(A) \le 2\pi y_p(A)/q$ & Explicit $\sqrt{A}$ Scaling:** On $x \in [c/2, n/2]$, dividing $-\delta\psi''(x)$ by $-U'(x, 0)$ yields $0 \le -\delta\psi''(x) \le C_{\text{ratio}}\frac{\kappa}{c}(-U'(x, 0))$, hence $0 \le U_y(x, 0) = -\delta\psi'(x) \le C_\psi\kappa^2 + C_{\text{ratio}}\frac{\kappa}{c}(-U(x, 0))$ and the exact global upper bound
   $$U(x, y) \;\le\; \Bigl(1 - C_{\text{ratio}}\frac{\kappa}{c}\, y\Bigr) U(x, 0) + C_\psi\kappa^2 y + C_{u''}\frac{\kappa}{c}\, y^2 \qquad \forall x \in [c/2, n/2],\; y \ge 0.$$
   Splitting $q = 2$ (where $\gcd(p_0, 2) = 1, p_0 \ge 3 \implies r_3/r_2 \le 0.668$) vs $q \ge 3$ ($\kappa \le 2/9$), this reduces $E_{p,\text{lb}}(A)$ on all 30 $T$-bins ($T \in [0.425, 14.0]$) to a 1-variable monotone bisection (plus exact $z = c/2$ Taylor expansion for $T \ge 14.0$, and explicit $w = \hat{y}/\sqrt{A}$ rescaling certifying $E_p(A) \ge c_{\sqrt{}}\sqrt{A}$ for all $A \ge A_0, q \ge \max(q_{\min}, A/\pi)$), giving $E_H(\pi q \ge 2\pi) \ge 8.30$ and $E_{H, C2}(\pi) \ge 5.28$.
3. **Zero-Hardcoded Top-Line, Side-Line, & Higher-Ray Tail Bounds:** Using the exact 1-line inequalities $U(x, y) - U(x, 0) \ge U(0, y) - U(0, 0) \ge -C_\psi\kappa^2 y + c_{y,\text{lo}}\frac{\kappa}{c} y^2$ (absorbing the linear term via AM–GM into both $c_{y,\text{eff}}$ and $B_0$) and $\frac{U_y(x, H)}{-U_x(x, H)} < \frac{(e^{\zeta_H}-1)(1 - se^{-\zeta_H})}{(1-s)\sin\xi} + \frac{1+s}{1-s}\tan\frac{\xi}{2}$ (with exact $c_{\text{sl,H,top}}$ combining $-U_x(t, H)/(-U_x(t, 0))$ and concavity of $\log(1+Q)$, plus $S_{\text{neg},U_y}$ for any $U_y < 0$), splitting $\frac{\beta u}{1-u} = -\beta - \frac{\beta}{u-1}$ on $|x| \le x_*(H)$ with 1st-IBP on $-\beta$, $Q \ge 1.00$ dead-zone ($2^{-\text{fac}_U n/2}$) + $L_0 = 0.25 n$ side-lines weighted by $q \le n/6$, and explicitly certified $\sqrt{q}/\sqrt{m}$ convergent higher-ray tails for all $q \ge 5, m \ge 2$, 50-digit `iv` certifies `line_C1 <= 0.015689` and
   $$|\mathcal{R}_{\text{non-}q}(z_{\text{var}})| \le 0.115416 \quad (|\psi| \le 1.25), \qquad |\mathcal{R}_{\text{non-}q}(z_{\text{var}})| \le 0.333745 \quad (|\psi| \le 2\pi/3).$$

#### Lemma G2.5 (Chebyshev-Mediant Pointwise-Split Far-Tail Decay on $|\psi| \in [2\pi/3, \pi]$, $q \ge 2$)
On $|\psi| = \pi - \beta$ ($\beta \in [0, \pi/3]$), replacing only the single near-zero ray $\theta_{t_0} = \beta/q$ by $\pi/(3q)$ restores the bound $|\mathcal{R}_{\text{non-}q}| \le 0.333745$ while multiplying each root $u_k$ by $e^{f(u_k; q, \beta)}$. By the Chebyshev mediant inequality $\frac{\sinh(qt)}{\sinh t} \ge q \ge \frac{\sin(q\theta)}{\sin\theta}$,
$$f(u; q, \beta) \le f_{\max}(u) := \frac{1}{2}\log\!\left(\Bigl(\frac{\tanh(u/2)}{u/2}\Bigr)^2\Bigl(\frac{u^2}{4} + \frac{\pi^2}{36}\Bigr)\right).$$
Splitting at $u_0 = 1.80$ with $c_0 = 1.52$ (using exact shifted threshold $u_0 + \log L_{\text{lo}}$ and wiring $H_{\text{cont,end}}(\tilde{\Lambda}) \ge 0.85841$ directly from Lemma G2.2c): (i) $f_{\max}(u) \le -1.52\,h(u)$ pointwise on $|u| \le 1.80$ (360 `iv` subintervals); (ii) on $|u| > 1.80$, $\max(\text{tail}/H_{\text{cont}}) \le 0.66073$ (`tail14 = 0.5903`), certifying
$$\frac{1}{q}\sum_{k=0}^{n-1} f(u_k; q, \beta) \;\le\; -c_{\text{far}}\, H_q^+(\tilde{\Lambda}) + \varepsilon_{\text{far,disc}}, \qquad c_{\text{far}} \ge 0.85927,\quad \varepsilon_{\text{far,disc}} \le 0.11787.$$

#### Lemma G2.7 (Exact Multiplicative Pigeonhole-Annulus Cubic Perturbation Lemma, ALL $n \ge 132$)
Rewriting Theorem 3.2's boundary factor $B(w)$ in terms of $z_{\text{var}} = Z^{-1} = (-1)^j w$ combines all $n$ roots ($c-1$ small, $2$ boundary, $n-c-1$ large) into a **purely multiplicative** ratio on $|w| = |z_{\text{var}}| = \tilde{\lambda} \le 2.40$ (zero additive $\Delta_B$ term and zero $\max|F_3|/\min|F_3|$ oscillation ratio):
$$\frac{F_{\text{full}}(w)}{F_{\text{pure}}(w)} = F_1(w)\, F_2(w)\, F_{\text{pref}}(w) \prod_{k=0}^{n-1}\frac{1 - \tilde{a}_k z_{\text{var}}}{1 - a_k z_{\text{var}}}, \qquad \tilde{a}_k := a_k\frac{1 + \delta_z}{(1 + a_k \varepsilon_s)^n},$$
where $F_{\text{pref}}(w) := \frac{(1+\varepsilon_+)(1+\varepsilon_-)}{1+\delta_z}\prod_{\text{small}}\frac{(1+a_k\varepsilon_s)^n}{1+\delta_z}$ and $|a_k| \le \frac{s}{1-s} < \frac{(R_y/R_z)^{1/n}}{\cos\alpha_1}$ for every $k = 0, \dots, n-1$.
Inside Box G2, $d = \gcd(c, j) \ge 2$ and $j > c$ force $j - c = d(p_0 - q) \ge 2$, while $c\kappa^2 < 42$ gives $\sin^2\alpha_1 \le \min(\sqrt{42}\sqrt{c}/n,\, \frac{42}{4\pi^2 c}) \implies \sin\alpha_1 \le (441/\pi^2)^{1/6} n^{-1/3}$, whence Lemma 4.1(2) yields
$$\frac{R_z}{R_y} \;\le\; \exp\!\Bigl(-2\pi\cot\frac{\alpha_1}{2}\Bigr) \;\le\; \exp\bigl(-c_{\text{se,G2}}\, n^{1/3}\bigr), \qquad c_{\text{se,G2}} \ge 6.5902.$$
By the Pigeonhole Principle, the interval $[\lambda_*(1 - \varepsilon_{\text{rad}}/n),\, \lambda_*(1 + \varepsilon_{\text{rad}}/n)]$ ($\varepsilon_{\text{rad}} := 10^{-4}$) contains a radius $\tilde{\lambda}$ separated from all $n$ root moduli $\{|a_k|^{-1}\}_{k=0}^{n-1}$ by at least $\frac{\varepsilon_{\text{rad}}\lambda_*}{n(n+1)}$. Hence $|a_k|$ cancels identically in $\bigl|\frac{1-\tilde{a}_k z_{\text{var}}}{1-a_k z_{\text{var}}} - 1\bigr| \le \frac{\tilde{\lambda}}{\lambda_*}\frac{n(n+1)}{\varepsilon_{\text{rad}}}\eta_{\text{root}}$, certifying in 50-digit `iv` across all $n \ge 132$:
$$\max_{|w| = \tilde{\lambda}} \Bigl|\frac{F_{\text{full}}(w)}{F_{\text{pure}}(w)} - 1\Bigr| \;\le\; E_{\text{cubic}} \;:=\; \begin{cases} 1.8182 \times 10^{-12} & (q = 1), \\ 1.6158 \times 10^{-31} & (q \ge 2). \end{cases}$$
Furthermore, since $H_q^+(\Lambda) \le \frac{n}{4q}$ everywhere and $|\log(\tilde{\Lambda}/\Lambda_*)| \le -q\log(1 - \varepsilon_{\text{rad}}/n)$, we have $|M_q^+(\tilde{\Lambda})| \le -\frac{N_0}{4}\log(1 - \varepsilon_{\text{rad}}/N_0)$ universally, so the linear saddle-shift phase drift on $[-\pi, \pi]$ is at most $\varepsilon_{\text{shift}} := -\frac{\pi N_0}{4}\log(1 - \varepsilon_{\text{rad}}/N_0) \le 7.854 \times 10^{-5}$.

#### Theorem 5.2 (50-Digit Certified FULL Exact Integral Positivity on ALL of Piece G2)
For every positive-root factor $1 + e^{u + i\psi}$, exact calculus yields (i) $-4\gamma(\psi)\, h(u)\sin^2(\psi/2) \le \log\frac{|1+e^{u+i\psi}|^2}{(1+e^u)^2} \le -4 h(u)\sin^2(\psi/2)$ with $\gamma(\psi) := -\log(\cos^2(\psi/2))/\sin^2(\psi/2)$; and (ii) the sharp phase bound
$$\Bigl|\frac{1}{q}\arg G_q^+(\tilde{\Lambda} e^{i\psi}) - d\psi\Bigr| \;\le\; H_q^+(\tilde{\Lambda})\, P_{\text{exact}}(\psi) + \varepsilon_{\text{shift}}, \qquad P_{\text{exact}}(\psi) := \begin{cases} \psi - \sin\psi & (0 \le \psi \le \pi/2), \\ (\frac{\pi}{2} - 1) - \log(2\cos^2(\psi/2)) & (\pi/2 < \psi < \pi). \end{cases}$$
Combining Lemmas G2.1–G2.7 in 50-digit outward-rounded `mpmath.iv` across all $H_q^+(\tilde{\Lambda}) \ge H_{\min}$ (subdivided up to $H = 60$ plus analytic Gaussian scaling for $H \ge 60$, using local per-bin $\mathcal{R}_{\text{non-}q}(\psi)$ on both $I_{\text{pos}}$ and $I_{\text{mid}}$ and subtracting $E_{\text{cubic}}$ directly in [`certify_G2_D0.py`](algebra/certify_G2_D0.py)) certifies for the **full exact cubic-root integral**:
- **Box G2, $q \ge 2$ ($H_q^+(\tilde{\Lambda}) \ge 0.59425$):** $\dfrac{I_{\text{full}}(\tilde{\lambda})}{L^1(\tilde{\lambda})} \;\ge\; \mathbf{+0.00176348} \;>\; 0$.
- **Box G2, $q = 1$ ($H_1^+(\tilde{\Lambda}) \ge 0.62018$, covering ALL $\kappa < 2.5$):** $\dfrac{I_{\text{full}}(\tilde{\lambda})}{L^1(\tilde{\lambda})} \;\ge\; \mathbf{+0.216651} \;>\; 0$.

All 28 exact finite-$n$ stress corners covering all four parity combinations of $(d, q)$ (including the $q = 1, \kappa \in (2.0, 2.5)$ band) pass Part B of [`certify_G2_D0.py`](algebra/certify_G2_D0.py) (`min I_pure/L1 = +0.97916`). $\square$

---

## 6. Computational Verification Scripts

All scripts live in `bardakov-kozlovskaya/algebra/`:
- [`certify_G2_D0.py`](algebra/certify_G2_D0.py): Complete 100% pure 50-digit outward-rounded `mpmath.iv` certification of **Piece G2 ($D = 0$ Major Arcs)** — Lemmas G2.1–G2.7 and Theorem 5.2 (zero hardcoded constants, explicit $Q_{\min} > 0.20$, exact 2nd-order Euler–Maclaurin, $P^2$-monotone $E_{\text{phase}}$, $\pi$-line meromorphic residue + Cauchy–Riemann level-curve $\mathcal{R}_{\text{non-}q}$, wired Chebyshev-mediant far-tail $c_{\text{far}}$, and exact multiplicative Pigeonhole-Annulus $E_{\text{cubic}}$) plus Part B finite-$n$ cross-check across all 4 parity combinations of $(d, q)$.
- [`verify_regime1_complete.py`](algebra/verify_regime1_complete.py): Computes all 6 closed-form Regime 1 margins for $\kappa_0 = 3.0$ and $\kappa_0 = 2.5$ (zero loops) AND cross-checks all 17 intermediate lemma inequalities across all 12,937 ($\kappa_0 = 3.0$) and 13,006 ($\kappa_0 = 2.5$) strictly interior points for $n \in \{84, 86, 88, 90, 100, 160, 320\}$ (0 failures).
- [`verify_exact_identities_mp.py`](algebra/verify_exact_identities_mp.py): 80-digit `mpmath` verification of Theorem 3.2 against exact integer coefficients of $Q_c^{(2m)}$ and Lemma 4.5's exact Cauchy integral identity for $B_{\text{avg}}$ across all 6 subcases at $n = 84$.
- [`verify_c1_and_probe_c2.py`](algebra/verify_c1_and_probe_c2.py): Verification of the $c = 1$ Rolle/Maclaurin proof and exact $\gamma_r$ identities.
- [`probe_depressed_cubic.py`](algebra/probe_depressed_cubic.py): Verification of Theorems 1.1–1.2 (Depressed-Cubic and Half-Size Lucas–Chebyshev resultants).
