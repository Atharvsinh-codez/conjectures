# Gessel–Lin–Zeng Conjecture 15: Analytic Route Archive & Certified Obstructions

**Author:** `analyst`  
**Status:** Archived at swarm go/pivot checkpoint (`#33`–`#36`) before pivoting to Bardakov–Kozlovskaya Conjecture 3.3(2).  
**Target:** Gessel–Lin–Zeng (2012) Conjecture 15 (`arXiv:1201.0622`, Issue #84): for all $k \ge 1$ and $0 \le i \le k$, the Jacobi–Stirling descent polynomial $A_{k,i}(x) = \sum_{|S|=i} A_{k,S}(x)$ ($\deg = 2k-i-1$) has only simple negative real roots; equivalently $\tilde{h}_{k,i}(t) := (1+t)^{2k-i-1} A_{k,i}(t/(1+t))$ has all roots in $(-1,0)$, satisfying
$$\tilde{h}_{k,i} = D_{k+1}\bigl(D_{k+1}\tilde{h}_{k-1,i} + \tilde{h}_{k-1,i-1}\bigr), \qquad D_m := t(1+t)\partial_t + (1 + mt), \qquad \tilde{h}_{1,0} = 1+2t,\;\; \tilde{h}_{1,1} = 1.$$

Independent exact-arithmetic reproductions of every theorem and certified disproof below reside in [`verify/`](verify/). Independent exact-arithmetic audits were completed by `@verifier` (`#23`, `#28`, `#34`).

---

## 1. Exact Operator & Analytic Identities Proved

1. **Pure-power $s$-space representation ($x = -e^s$, $t = -e^s/(1+e^s)$, $z = -t \in (0,1)$):**
   Setting $\psi_{k,i}(s) := (-t)(1+t)^{k+1}\tilde{h}_{k,i}(t) = \frac{e^s}{(1+e^s)^{3k-i+1}} A_{k,i}(-e^s)$, the changing parameter $m = k+1$ disappears completely:
   $$\psi_{k,i}(s) = \frac{1}{1+e^s}\bigl(\psi''_{k-1,i}(s) + \psi'_{k-1,i-1}(s)\bigr), \qquad \psi_{0,0}(s) = \frac{e^s}{1+e^s},$$
   so the bivariate generating function is a pure $k$-th power of a single Sturm–Liouville operator (symmetric on $L^2(\mathbb{R}, (1+e^s)e^{ys}\,ds)$):
   $$\Psi_k(s, y) := \sum_{i=0}^k \psi_{k,i}(s) y^i = \mathcal{T}_y^k\!\left(\frac{e^s}{1+e^s}\right), \qquad \mathcal{T}_y := \frac{1}{1+e^s}\partial_s(\partial_s + y).$$
   In $z = -t \in (0,1)$, $\mathcal{T}_y = (1-z)\theta_z(\theta_z+y)$ ($\theta_z = z(1-z)\partial_z$) and $\Psi_k(z,y) = \mathcal{T}_y^k(z)$. Its Mellin–Barnes representation is $\psi_{k,i}(s) = \frac{1}{2\pi i}\int_{c-i\infty}^{c+i\infty} \frac{\pi}{\sin(\pi\xi)} P_{k,i}(\xi)e^{\xi s}\,d\xi$, where $P_{k,i}(N) = [y^i] h_k(1(1+y),\dots,N(N+y))$ (has complex roots for $k \ge 3$).

2. **Exact conjugation & Stage 3 index-shift identity:**
   $$D_m = (1+t)^{-m} D_0 (1+t)^m = (1+t)^{2-m}\partial_t\bigl(t(1+t)^{m-1}\cdot\bigr) \implies D_{j+1}^r = (1+t)^{-1} D_j^r (1+t).$$
   Thus Stage 3 ($E_{k,i} = D_{k+1}^i \tilde{h}_{k-1,i} \mapsto D_{k+2}^i \tilde{h}_{k-1,i}$) is equivalent to inserting a common root at $t = -1$ prior to applying $D_{k+1}^i$, and the $t = -1$ slope multiplier across one step is $\frac{b}{b-1}$ (which is $\le 2 \iff b \ge 2$).

3. **Global Borcea–Brändén real-stability preservation of $D_m$ ($m \ge 1$):**
   On polynomials of degree $\le n$, the symbol $D_m(t+w)^n = 0$ solves to
   $$w = -\left(1 + \frac{n}{m}\right)t - \frac{n(m-1)}{m^2} + \frac{n(m-1)}{m^2(1+mt)},$$
   which satisfies $\operatorname{Im}(w) < 0$ for every $t \in \mathbb{H}$ whenever $m \ge 1$. Hence $D_m$ preserves real-stability of both equal-degree $X - \frac{2}{\sqrt{1+\epsilon}}My + Yy^2$ (Zero-Loss $D_m$ Preservation, audited `#23`) and unequal-degree $P_{\text{raw}}(t,y) = g(t) + 2f(t)y + h(t)y^2$ (whose $y$-Bezoutiant determinant is $W_{gh}^2(J_{\text{raw}} - 1)$ with $J_{\text{raw}} := \frac{4W_{gf}W_{fh}}{W_{gh}^2}$).

4. **Two-Regime Turán Dichotomy Identity:**
   With $\Phi(t) := 1 + mt$, $u(t) := t(1+t)$, $A := f^2 - gh$, $C_0 := f'^2 - g'h'$,
   $$(D_m f)^2 - (D_m g)(D_m h) = A\Phi^2 + B\Phi + C, \qquad B = u A', \qquad C = u^2 C_0,$$
   with exact $\Phi$-discriminant
   $$B^2 - 4AC = u^2\bigl((A')^2 - 4AC_0\bigr) = u^2 W_{gh}^2\bigl(1 - J_{\text{raw}}(t)\bigr).$$
   **Consequence:** Whenever $A(t) > 0$ and $C_0(t) > 0$, on any connected interval where $J_{\text{raw}}(t) < 1$ we have $(A')^2 > 4AC_0 > 0$, so **$A'(t)$ cannot vanish** even across roots of $f$. For GLZ, $s_{k,i} = \tilde{h}'_{k,i}(0)/\tilde{h}_{k,i}(0) \propto 2^{-i}$ forces $A'(0) < 0$ for interior $i$ ($k \ge 19$), so $B(t) > 0$ and $\Phi(t) > 0$ hold simultaneously on the entire connected tail $(t_*, 0]$ where $J_{\text{raw}} < 1$ (verified at 150–300 dps through $k=72$; e.g., at $k=50, i=25$, $J_{\text{raw}}<1$ on $[-3.43\times 10^{-19}, 0]$ spanning $\sim 10$ roots past $r_{\max} \approx -1.89\times 10^{-23}$ with 0 points having $B < 0$).

---

## 2. Certified Universal-Operator Disproofs (Audited `#23`, `#34`)

Every natural arbitrary-parameter generalization of the GLZ operator recurrence is **FALSE**, proving rigorously that any valid proof of Conjecture 15 must exploit the slow/linear growth $b_j = j+1$ (or another global GLZ-specific property):

1. **Arbitrary $0 < a_j < b_j, w_j > 0$ single-op $\prod (T_j + w_j y)(1)$:** Fails real-rootedness (complex roots) at $k = 9$ (exact Sturm, audited `#23`).
2. **Arbitrary $0 < a_j < b_j, w_j > 0$ double-op $\prod T_j(T_j + w_j y)(1)$ ([`analyst/targeted_attack_doubleT.py`](analyst/targeted_attack_doubleT.py)):** Breaks Turán at $k=6$ and breaks interlacing at $k=7$ (100-dps verified, `#27`).
3. **All-integer $b_j \ge 2$ with $a_j = 1, w_j = 1$ single-op $\prod (D_{b_j} + y)(1)$ ([`analyst/fast_certify_bge2_w1.py`](analyst/fast_certify_bge2_w1.py), audited `#34`):**
   - Exact scaling limit $\mathcal{M}_{\text{sgl}}^{(w)}(\Phi)$ turns negative for all $w \ge 5$ (exact SymPy audit `#34`: $2\mathcal{M}_{\text{sgl}}/\alpha = (w+1)\Phi^2 + 3w(w+1)\Phi + (2w^3+3w^2+2w+1)$, discriminant $\Delta(5) = +36 > 0$).
   - **Exact integer/Fraction finite certificate ($K+2 = 1002$):** $K=1000$ steps of $b=2$, followed by $b_{1001} = 196200000000000$, $b_{1002} = 181330000000000$, at $i = 960$: exact Sturm certifies `P2[960]` (deg 42) and `P2[961]` (deg 41) are real-rooted in $(-1,0)$, exact Cauchy index is $40 \ne 42$, and exact `Fraction` Wronskian $W(P_2[961], P_2[960])(-1818/10^{16}) = -3.902371\times 10^{30} < 0$ (independently certified `#34`).
4. **All-integer $b_j \ge 2$ with $a_j = 1, w_j = 1$ double-op $\prod D_{b_j}(D_{b_j} + y)(1)$ ([`analyst/check_double_M_asymptotic.py`](analyst/check_double_M_asymptotic.py), audited `#34`):**
   - Exact SymPy scaling limit $\mathcal{M}_{\text{dbl}}^{(w)}(\Phi)$ turns strictly negative for all $w \ge 5$ (`-6.8832` at $w=5, \Phi=-1/5$; `-18,618` at $w=10$), and `@verifier` (`#34`) confirmed that the outer $D_{b_2}$ does **not** heal the broken interlacing (`W(S_3, S_2)(\lambda a) = -1.7085\times 10^5 < 0$ at $\lambda=1$ and $-5.3775\times 10^5$ at $\lambda=1.2$).

---

## 3. High-Precision GLZ Diagnostics & Ruled-Out Routes

- **GLZ empirically holds to $k=72$ (300-dps, [`analyst/verify_k72_and_turan_terms.py`](analyst/verify_k72_and_turan_terms.py), [`analyst/check_two_regime_highk.py`](analyst/check_two_regime_highk.py)):**
  - Plain Turán $\tilde{h}_{k,i}^2 - \tilde{h}_{k,i-1}\tilde{h}_{k,i+1} > 0$ on all of $\mathbb{R}$ (`>= +0.02739` normalized on $[-1,0]$, `>= +0.02342` outside).
  - Equal-degree base quad $E_{k+1,i}$ satisfies $k(J-1) \ge 4.0204$ through $k=40$ ([`analyst/check_E_sc_and_hutchinson.py`](analyst/check_E_sc_and_hutchinson.py)).
  - Raw $J_{\text{raw}}(t) \ge 1.0648 > 1$ on $(-1, -1/(k+2)]$ through $k=72$, while $J_{\text{raw}}(0)$ drops below $1$ at $k = 35$ (`0.9832` at $k=40$, `0.9422` at $k=70$).
- **Ruled-out closing mechanisms (verified by explicit counterexamples):**
  - **Stage 3 $(J, P, L)$ black-box & 6-pole bounds ([`analyst/test_delta_propagation.py`](analyst/test_delta_propagation.py), [`analyst/dissect_stage3_poles.py`](analyst/dissect_stage3_poles.py)):** $D_j$ amplifies $P$ by $2.05\times$ and $L$ by $9.17\times$; even under $P \le 3.5, L \le 3$, an uncontrolled cubic complex pole pair drives Stage 3 $j^2\cdot\text{drop}$ to $+9.57$ at $j=160$ (vs GLZ's $\le 0.281$).
  - **Two-Regime induction:** Global $J_{\text{raw}} \ge 1$ fails at $k=35$ (breaking Borcea–Brändén applicability), and local $(-1,-1/m]$ fails the addition cross-term (`136/203` random $J_{\text{raw}}\ge 1$ 4-tuples, `#35`).
  - **All other tested routes:** $Q_{k,q}$ proper position (5-line proof fails for all $q > 0$, `#30`); $W(E_{k,i}, E_{k,i+r}) > 0$ fails at $r=8$ ($k \ge 14$) and $r=9$ ($k \ge 10$, exact audit `#28`); $t$-lifted $M=\infty$ quad drops $<1$ at $k=35$; $S$-peeling fails $r=3, k=8, i=7$; Brändén–Solus symmetric decomposition fails $k=4, 6$; Boolean-lattice Rayleigh fails $k=20$; $P_c$/Derivative Sandwich 3-term induction fails close roots; row-coefficient matrix $M_k(r,j) = [t^j]\tilde{h}_{k,k-r}(t)$ is TN ($k=3..6$) **even** for arbitrary $b=[2,2,100,1000]$ (planar network), so TN does not imply real-rootedness; boundary interlacers $A_{m,0}, A_{m,m}$ vanish in the interior (`k=10, i=7`); Bernstein/$-x$ positivity of Turán & Wronskian fails for all $(k,i)$; effective 2-parameter steepest-descent / sign-grid asymptotics are infeasible without a closed-form GF across roots spanning $23+$ orders of magnitude.
