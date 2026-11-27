# Compendium of Core Equations & Analytic Identities

> **Authors:** **[Atharsinh Jadav](https://athrix.me/)** ([X `@athrix_codes`](https://x.com/athrix_codes) · [GitHub `Atharvsinh-codez`](https://github.com/Atharvsinh-codez)), with the **Gemini 4 Neon Extra High multi-agent research swarm**  
> **Scope:** Every foundational algebraic, combinatorial, complex-analytic, and operator identity established or utilized across **Bardakov–Kozlovskaya Conjecture 3.3(2)** (Part I, solved) and **Gessel–Lin–Zeng Conjecture 15** (Part II, open).

---

## Part I. Bardakov–Kozlovskaya Conjecture 3.3(2)

### 1. Defining Polynomials, Symmetries & Classical Representations

**Buchstaber–Novikov $n$-valued group polynomial** ($n \ge 1$, $\zeta = e^{2\pi i/n}$, $\alpha^n = x, \beta^n = y, \gamma^n = z$):
$$p_n(x, y, z) = \prod_{k=0}^{n-1}\Bigl(z - (-1)^n(\alpha + \zeta^k \beta)^n\Bigr) = \prod_{r=0}^{n-1}\prod_{s=0}^{n-1}\bigl(\gamma + \zeta^r\alpha + \zeta^s\beta\bigr).$$

**Elementary symmetric expansion & exact indexing bijection** ($e_1 = x+y+z$, $e_2 = xy+yz+zx$, $e_3 = xyz$; $a = k_1 - k_2, b = k_2 - k_3, c = k_3$):
$$P_n(e_1, e_2, e_3) = \sum_{\substack{k_1 \ge k_2 \ge k_3 \ge 0 \\ k_1 + k_2 + k_3 = n}} A_{k_1,k_2,k_3}\, e_1^{k_1-k_2} e_2^{k_2-k_3} e_3^{k_3} = \sum_{\substack{a,b,c \ge 0 \\ a + 2b + 3c = n}} c_{a,b,c}^{(n)}\, e_1^a e_2^b e_3^c.$$

**Univariate $c$-slice polynomial** ($n = 2m$ even, $D_c := \lfloor (2m - 3c)/2 \rfloor$):
$$P_{2m}(e_1, e_2, e_3) = \sum_{c=0}^{\lfloor 2m/3 \rfloor} e_1^{2m-3c}\, Q_c^{(2m)}\!\Bigl(1, \frac{e_2}{e_1^2}\Bigr)\, e_3^c, \qquad Q_c^{(2m)}(1, u) := \sum_{b=0}^{D_c} c_{2m-3c-2b,\, b,\, c}^{(2m)}\, u^b.$$

**Universal Sign Law (The Jadav Sign-Law & Real-Rootedness Theorem — ALL even $n = 2m \ge 2$):**
$$\operatorname{sign}\!\bigl(c_{a,b,c}^{(2m)}\bigr) = (-1)^{b(c+1)+c} = (-1)^{b(a+1)+c} = (-1)^{k_1 k_2 + k_1 k_3 + k_2 k_3} = \begin{cases} (-1)^b & \text{if } c \text{ is even,} \\ -1 & \text{if } c \text{ is odd.} \end{cases}$$

**Resultant, Wendt $y$-circulant determinant, GRS discriminant & Newton power sums** ([`SOURCES_AND_STATE.md`](bardakov-kozlovskaya/SOURCES_AND_STATE.md), [`verify/engine.py`](bardakov-kozlovskaya/verify/engine.py)):
$$p_n(x, y, t^n) = \operatorname{Res}_u\bigl(u^n - x,\; (t + u)^n - (-1)^n y\bigr),$$
$$p_n(w^n; x, y) = \det \mathrm{Circ}_y\!\left(w^n + (-1)^{n+1}x + y,\; \binom{n}{1}w,\; \binom{n}{2}w^2,\; \dots,\; \binom{n}{n-1}w^{n-1}\right),$$
$$(-1)^n (n-1)^{2n-2} (xyz)^{n-2}\, p_n(x, y, z) = \operatorname{disc}_T\!\Bigl(\bigl((-1)^n x T^{n-1} + (-1)^n y\bigr)(1 + T)^{n-1} - T^{n-1} z\Bigr),$$
$$S_r(x, y) := \sum_{k=0}^{n-1} \bigl((-1)^n(\alpha + \zeta^k \beta)^n\bigr)^r = n(-1)^{nr}\sum_{a=0}^r \binom{nr}{na} x^{r-a} y^a.$$

---

### 2. Two Exact Structural Resultants & All Four Boundary Families

**Depressed-Cubic Resultant (Theorem 3.1):** With Waring polynomials
$$\Pi_{2m}(u, \eta) := \sum_{r=0}^{\lfloor m/3 \rfloor} \frac{2m}{m - r}\binom{m - r}{2r} u^{m - 3r}\eta^r, \qquad \Sigma_{2m}(u, \eta) := \sum_{j=0}^{\lfloor 2m/3 \rfloor} (-1)^j \frac{m}{m - j}\binom{2m - 2j}{j} u^{2m - 3j}\eta^j,$$
$$P_{2m}(e_1, e_2, \eta^m) = \operatorname{Res}_u\bigl(\Pi_{2m}(u, \eta) - e_1,\; \Sigma_{2m}(u, \eta) - e_2\bigr) \quad (3 \nmid m), \qquad P_{6K}(e_1, e_2, e_3) = \prod_{\omega^3 = 1} R_K(e_1, e_2, \omega e_3^{1/3}).$$

**Half-Size Lucas–Chebyshev Resultant (Theorem 3.2):** With $L_{2m}(u) := \sum_{r=0}^{m} (-1)^r \frac{2m}{2m - r}\binom{2m - r}{r} u^r$ satisfying $L_{2m}(t(1-t)) = t^{2m} + (1-t)^{2m}$:
$$P_{2m}(1 + q - p,\, q,\, p) = \operatorname{Res}_u\!\bigl(L_{2m}(u) + u^{2m} - q,\; u^{2m} - p\bigr) = \prod_{k=0}^{2m-1}\Bigl(L_{2m}(\zeta^k p^{1/(2m)}) + p - q\Bigr).$$

**Boundary $c = 0$ (Theorem 3.3):**
$$Q_0^{(2m)}(e_1, e_2) = (e_1^2 - 4e_2)^m \quad \Longrightarrow \quad c_{2m-2b,\, b,\, 0}^{(2m)} = (-4)^b \binom{m}{b}.$$

**Boundaries $a = 0$ and $b = 0$ (Theorem 3.4 — trigonometric evaluation & exact derivatives):**
$$L_{2m}\!\Bigl(\frac{1}{4\cos^2\theta}\Bigr) = \frac{2\cos(2m\theta)}{(2\cos\theta)^{2m}},$$
$$\frac{d}{d\theta}\!\left[\frac{2\cos(2m\theta)}{(2\cos\theta)^{2m}}\right] = \frac{-8m\sin((2m - 1)\theta)}{(2\cos\theta)^{2m+1}}, \qquad \frac{d}{d\theta}\!\Bigl[2(2\cos\theta)^{2m}\cos(2m\theta) + 1\Bigr] = -8m(2\cos\theta)^{2m-1}\sin((2m + 1)\theta).$$
At critical points $\theta_r = \frac{r\pi}{2m-1}$ ($a=0$) and $\tilde{\theta}_r = \frac{r\pi}{2m+1}$ ($b=0$):
$$h_a(\theta_r) = 1 + (-1)^r (2\cos\theta_r)^{1 - 2m}, \qquad h_b(\tilde{\theta}_r) = 1 + (-1)^r (2\cos\tilde{\theta}_r)^{2m + 1}.$$

**Boundary $c = 1$ (Theorem 3.5 — Rolle recurrence, palindromic factorization & alternating trig sum):**
$$(r - 1)\, F_r(u) = n u(1 - u)\, F'_{r-1}(u) + (r - 1)\bigl(1 + (n - 1)u\bigr)\, F_{r-1}(u), \qquad F_1(u) = 1,\;\; A_{2m-1}^{(2m)}(u) = F_{2m}(u),$$
$$\frac{A_{2m-1}^{(2m)}(u)}{1 + u} = \prod_{j=1}^{m-1}\bigl((1 + u)^2 + c_j u\bigr) = \sum_{r=0}^{m-1} \gamma_r (1 + u)^{2m-2-2r} u^r \quad \Longrightarrow \quad c_{2m-3-2b,\, b,\, 1}^{(2m)} = 2m\!\left((-4)^{b+1}\binom{m-1}{b+1} - \gamma_{b+1}\right),$$
$$\gamma_{m-1} = \frac{1}{2m}\sum_{k=0}^{m-1} (-1)^k \frac{\cos\theta_k}{\sin^{2m+1}\theta_k} > 4^{m-1} \qquad \left(\theta_k := \frac{(2k + 1)\pi}{4m}\right).$$
*(Special case $c = 1, b = 0$: $c_{n-3, 0, 1}^{(n)} = -n\bigl(\frac{1}{2}\binom{2n}{n} - 2 - (-1)^n\bigr)$.)*

---

### 3. Universal Sample Points & Exact 4-Factor 1D Circle Integral

**Universal trigonometric sample points ($j \in \mathcal{J}_{n,c} = \{c+1, \dots, m - \lceil c/2\rceil - \mathbf{1}_{c\text{ even}}\}$, Theorem 4.1):**
$$\theta_c = \frac{c\pi}{n}, \quad \psi_j = \frac{j\pi}{n}, \quad \alpha_1 = \frac{(j+c)\pi}{n}, \quad s_j = \frac{\sin\psi_j}{\sin\alpha_1}, \quad R_y = s_j^n, \quad R_z = \Bigl(\frac{\sin\theta_c}{\sin\alpha_1}\Bigr)^n, \quad t_j = (-1)^c R_y, \quad u_j = \frac{t_j}{(1 + t_j)^2},$$
$$\kappa(n, c, j) := 2\pi \frac{\sin\psi_j \sin\alpha_1}{\sin\theta_c}, \qquad d := \gcd(c, j), \qquad q := \frac{c}{d}, \qquad X := q^2 \kappa.$$

**Exact 1D Cauchy circle integral & 4-factor decomposition (Theorem 4.2):** On the cubic $T(T-1)(T-t) = t R_z w$ ($|w| = \lambda$, $Z = (-1)^j w^{-1}$):
$$(-1)^j (1 + t)^{n - 3c} Q_c^{(2m)}(1, u_j) = \frac{C_{m,c,j}}{2\pi}\int_{-\pi}^{\pi} F_1(\lambda e^{i\phi})\, F_2(\lambda e^{i\phi})\, F_3(\lambda e^{i\phi})\, B(\lambda e^{i\phi})\, d\phi,$$
$$B^{\text{pure}}(e^{i\phi}) = 2 - Z - Z^{-1} = 4\sin^2\!\Bigl(\frac{\phi + j\pi}{2}\Bigr) \ge 0, \qquad B_{\text{avg}} = \frac{1}{2\pi}\int_{-\pi}^\pi B(e^{i\phi})\,d\phi = 2 - (-1)^j\bigl(2 g(\theta_c)\delta'_y(0) - \delta'_z(0)\bigr).$$

**Jensen midpoint root-separation & off-boundary modulus identity (Lemmas 4.1 & 4.4):**
$$\Bigl(\frac{\sin u}{\sin v}\Bigr)^n \le \exp\!\Bigl(-n(v - u)\cot\frac{u + v}{2}\Bigr), \qquad \frac{|w_k|^2}{R_z^{2/n}} = 1 + \frac{2\kappa}{\pi\sin\theta_c}\sin\!\Bigl(\frac{|\eta_k|+\theta_c}{2}\Bigr)\sin\!\Bigl(\frac{|\eta_k|-\theta_c}{2}\Bigr).$$

---

### 4. Regime 2 Analytic Engines (Pieces G2, G1 & G3)

**Piece G2 ($D = 0$ Major Arcs) — All-parity contour identity, Cauchy–Riemann level-curve law & Pigeonhole-Annulus ratio (Prop. 5.4, Lemmas G2.4 & G2.7):**
$$F_{\text{pure}}(w) = \Bigl(\prod_{\text{small}} |u_k|\Bigr)\,\lambda^{-c}\, e^{-i d\psi}\,\prod_{k=0}^{n-1}(1 - a_k z_{\text{var}}), \qquad z_{\text{var}} = (-1)^j w = \lambda\exp\!\Bigl(i\frac{(2m+1)\pi + \psi}{q}\Bigr),$$
$$\frac{d}{dy}\, V(x_*(y), y) = \frac{U_x^2 + U_y^2}{U_x} < 0 \qquad \bigl(\log b(x+iy) = U + i V,\; U(x_*(y), y) \equiv -\log\tilde{\lambda}\bigr),$$
$$\frac{F_{\text{full}}(w)}{F_{\text{pure}}(w)} = F_1(w) F_2(w) F_{\text{pref}}(w) \prod_{k=0}^{n-1}\frac{1 - \tilde{a}_k z_{\text{var}}}{1 - a_k z_{\text{var}}}, \qquad \tilde{a}_k = a_k\frac{1 + \delta_z}{(1 + a_k\varepsilon_s)^n}.$$

**Piece G1 (Dirichlet Minor Arcs) — Hubbard–Stratonovich duality & 1-line $D_r^2$-cancellation identity (Lemmas G1.1–G1.2):**
$$2\pi^2 \nu^2 D(t\pi\nu, \nu) = \mathbb{E}_{Z \sim \mathcal{N}(0, 1)}\!\left[Z^2\, F\!\left(t, \sqrt{\frac{t}{\pi\nu}}\, Z\right)\right], \qquad \frac{1}{4\pi^2}\frac{c}{r D_r^2} H_D(t_r) = \frac{t_r^2 H_D(t_r)}{c\kappa^2 r^3} \le \frac{U_{2H}}{c\kappa^2 r^3}.$$

**Piece G3 ($D \ne 0$ Major Arcs, $r_0 j = c p_0 + D$) — Phase-unwrapped $\zeta(u)$, $(1-w)$-free IBP kernels, Harmonicity & Shifted-Line Poisson Identity (Sub-Lemmas G3.1–G3.3):**
$$\zeta(u) := r_0 n\log\!\left(\frac{r_3/r_1}{1 - s_j e^{2iu\theta_c}}\right) + i\pi r_0 j - 2\pi i c p_0\Bigl(u + \frac{1}{2}\Bigr), \qquad \zeta'(u) = \frac{A(u)}{1 - w(u)}, \qquad \zeta''(u) = \frac{-B w(u)}{(1 - w(u))^2},$$
$$W_{1,\text{fin}} = \frac{Bw}{A^2}, \qquad W_{2,\text{fin}} = \frac{3 B^2 w^2 + A B (2i\theta_c) w(1+w)}{A^4},$$
$$W_{3,\text{fin}} = \frac{15 B^3 w^3 + 10 A B^2 (2i\theta_c) w^2(1+w) - A^2 B (4\theta_c^2) w(1+4w+w^2)}{A^6},$$
$$\Psi_{\text{fin}}(\theta; u_0) := -\log\Lambda - i\theta + \int_{-L_n}^{L_n} \log(1 + \Lambda e^{\zeta(u)+i\theta})\,du, \qquad \frac{d}{du_0} g_{1,\text{fin}}(u_0) = \bigl(-\operatorname{Re}\Psi_{\text{fin}}''(0; u_0)\bigr)\bigl(-\operatorname{Re}\zeta'(u_0)\bigr),$$
$$\log F_{\text{pure}}(\lambda e^{i(\psi_\nu + \theta/r_0)}) - \tilde{c}\,\Psi_{\text{fin}}(\theta; u_0) = \sum_{p \in \text{strip}} \log\!\bigl(1 - \omega_p e^{2\pi i\sigma_p \tilde{c}(u_p+1/2)}\bigr) + \sum_{k\ne 0} \frac{e^{2\pi i\nu\mu(k)/r_0}}{2\pi i k}\int_{-L_n+i y_k}^{L_n+i y_k} e^{2\pi i k\tilde{c}(u+1/2)} \zeta'(u)\frac{y_*(u,\theta)^{\mu(k)}}{1 + \Lambda e^{\zeta(u)+i\theta}}\,du.$$

---

## Part II. Gessel–Lin–Zeng Conjecture 15

### 1. Recurrences, Operator Calculus & Pure-Power Sturm–Liouville Form

**Original insertion recurrence & Möbius-transformed $i$-independent recurrence** ($S_L f(x) = (1+Lx)f(x) + x(1-x)f'(x)$, $x = t/(1+t)$):
$$A_{k,i}(x) = S_{3k-2-i}\bigl(S_{3k-3-i} A_{k-1,i}(x) + A_{k-1,i-1}(x)\bigr), \qquad A_{1,0}(x) = 1 + x,\;\; A_{1,1}(x) = 1,$$
$$\tilde{h}_{k,i}(t) := (1+t)^{2k-i-1} A_{k,i}\!\Bigl(\frac{t}{1+t}\Bigr), \qquad \tilde{h}_{k,i} = D_{k+1}\bigl(D_{k+1}\tilde{h}_{k-1,i} + \tilde{h}_{k-1,i-1}\bigr), \qquad D_j := t(1+t)\partial_t + 1 + jt.$$

**Conjugation, commutation, polar intertwining & equal-degree lift** ($P_c f(x) = c f(x) - x f'(x)$, $E_{k,i} = D_{k+1}^i \tilde{h}_{k-1,i}$):
$$D_j f(t) = (1+t)^{2-j}\partial_t\bigl[t(1+t)^{j-1} f(t)\bigr], \qquad (1+t) D_{j+1}^r = D_j^r (1+t), \qquad [D_j, D_{j'}] = (j' - j) t(1+t),$$
$$P_{L+1} S_{L+1} = S_L P_{L+1}, \qquad \sum_{i=0}^k D_{k+1}^i \tilde{h}_{k,i}(t)\, y^i = D_{k+1}^2\!\Bigl((1+y)\sum_{i=0}^{k-1} E_{k,i}(t)\, y^i\Bigr).$$

**Pure-power Sturm–Liouville generating function ($x = -e^s$, $\psi_{k,i}(s) = (-t)(1+t)^{k+1}\tilde{h}_{k,i}(t)$):**
$$\Psi_k(s, y) := \sum_{i=0}^k \psi_{k,i}(s)\, y^i = \mathcal{T}_y^k\!\left(\frac{e^s}{1 + e^s}\right), \qquad \mathcal{T}_y := \frac{1}{1 + e^s}\partial_s(\partial_s + y).$$

---

### 2. Real-Stability, Plücker Quad, Wronskian & Turán Identities

**Borcea–Brändén symbol partial fraction for $D_m$ on degree $\le n$ ($D_m(t+w)^n = 0$):**
$$w = -\left(1 + \frac{n}{m}\right)t - \frac{n(m-1)}{m^2} + \frac{n(m-1)}{m^2(1 + mt)} \qquad \bigl(\operatorname{Im} w < 0\text{ on }\mathbb{H}\iff m \ge 1\bigr).$$

**Two-Regime Turán Dichotomy** ($\Phi(t) = 1 + mt$, $u(t) = t(1+t)$, $A = f^2 - gh$, $C_0 = f'^2 - g'h'$, $J_{\text{raw}} = \frac{4W_{gf}W_{fh}}{W_{gh}^2}$):
$$(D_m f)^2 - (D_m g)(D_m h) = A\Phi^2 + B\Phi + C, \qquad B = u A', \qquad C = u^2 C_0, \qquad B^2 - 4AC = u^2 W_{gh}^2\bigl(1 - J_{\text{raw}}(t)\bigr).$$

**Sharpest Quantitative Plücker Quad Inequality & Universal Cubic Wronskian Log-Concavity:**
$$J(A + B,\; B + C,\; C + E) \;\ge\; \frac{4(1 + \varepsilon)^2}{(2 + \varepsilon)^2} \qquad \Bigl(J(X, M, Y) := \frac{4 W(X, M) W(M, Y)}{W(X, Y)^2}\Bigr),$$
$$W(p_0, p_1)\, W(p_2, p_3) \;\le\; W(p_1, p_2)^2 \qquad \forall\text{ real-stable cubic } P(t, y) = \sum_{j=0}^3 p_j(t)\, y^j.$$

**Constant-$b$ $t = 0$ Wronskian asymptote & Exact GLZ Disproof Witnesses:**
$$J_i(0) = \frac{8}{9}\frac{(i+1)(k-i+1)}{i(k-i)}, \qquad \Bigl(B_{10,x}B_{10,y} - B_{10}B_{10,xy}\Bigr)\!\Bigl(0, -\frac{3}{4}\Bigr) = -\frac{2649516238079970755}{4294967296} < 0,$$
$$\bigl[R'(x)^2 - R(x)R''(x)\bigr]\!\Big|_{R = 932 A_{31,\{28..31\}} + A_{31,\{1..4\}},\, x=-23} \approx -9.065320\times 10^{279} < 0.$$

---

*See also: [`THEOREMS.md`](THEOREMS.md) · [`JOURNEY.md`](JOURNEY.md) · [`bardakov-kozlovskaya/PROOF.md`](bardakov-kozlovskaya/PROOF.md) · [`gessel-lin-zeng/AUDIT.md`](gessel-lin-zeng/AUDIT.md).*
