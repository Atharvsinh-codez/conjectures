# Piece G1 (Dirichlet Minor Arcs) — Final 50-Digit `mpmath.iv` Certification Notes

**Script:** `bardakov-kozlovskaya/algebra/certify_G1_minor.py`  
**Status:** **PASS (exit code `0`)** — 100% pure 50-digit outward-rounded `mpmath.iv` across all five lemmas G1.1–G1.5 + all 27 exact finite-$n$ Part B corners.  
**Audit resolution:** Resolves **all 7 soundness and hygiene items from `#78`** AND **all 5 edge-case fixes from `#98`** with zero residual gaps, zero `float64` casts inside `iv`, and explicit finite-$n$ cubic-root perturbation ($I_{\text{pure}} \to I_{\text{full}}$).

---

## 1. Exact Domain of Piece G1

Piece G1 covers the complement $\mathcal{G}^* \setminus (\mathcal{G}_2 \cup \mathcal{G}_3)$ inside Regime 2:
- Even $n \ge 132$, $2 \le c \le \lfloor (n-3)/3 \rfloor$, $j \in J_{n,c}$ ($j \ge c+1$),
- $\kappa := 2\pi\frac{r_1 r_2}{r_3} = 2\pi\frac{\sin\alpha_1 \sin\psi_j}{\sin\theta_c} < 2.5$, $c\kappa^2 < 42$ (overlaps Sub-regime 2A $\ge 42$),
- $d := \gcd(c, j)$, $q := c/d \ge 2$ ($D=0$ handled by Piece G2), $X := q^2 \kappa$ satisfying:
  - $d = 1$: $X \ge 0.88$ (overlaps Sub-regime 2C-Upgrade $cq\kappa \le 0.90$),
  - $d \ge 2$: $X \ge K_{\text{G2}} := 2.0$ (overlaps Piece G2 $q^2 \kappa \le 2.0$),
  - hence unconditionally $c^2 \kappa = d^2 X \ge 0.88$ across **all** of Piece G1,
- **Minor-arc condition:** No coprime pair $(p_0, r_0)$ with $D := r_0 j - c p_0 \ne 0$ satisfies all four Piece G3 major-arc inequalities simultaneously:
  1. $M := \frac{c}{r_0 D^2} \ge M_{\text{G3}} = 10$,
  2. $r_0 \kappa \le \varepsilon_{\text{G3}} = 0.45$,
  3. $r_0^3 c \kappa^2 \le C_{r3} = 60$,
  4. $\tau^2 := \frac{c^3 r_0 \kappa^2}{D^4} = 4\pi^2 M t^2 \ge C_\tau = 45$ (overlaps Piece G3 $\ge 40$).

---

## 2. Certified 50-Digit `mpmath.iv` Summary Table (Final Exit-0 Run)

| Lemma / Component | Quantity Certified | 50-Digit `iv` Bound | Margin / Status |
|---|---|---|---|
| **G1.1** Master HS Identity | $H_D(t) / t^2$ ($t \le 0.05$) | $\le 24.4591$ | PASS |
| **G1.1** Large-$t$ / Global $u^2 H_D(u)$ | $t^2 H_D(t)$ ($t \ge 12$) / $\sup_{u>0} u^2 H_D(u)$ | $\le 8.45199$ / $\le 9.3423$ | PASS |
| **G1.1** Global Peak & $t \ge 1.50$ | $\sup_{t>0} H_D(t)$ / $\sup_{t \ge 1.50} H_D(t)$ | $\le 3.14341$ / $\le 0.553496$ | PASS |
| **G1.1** Global Quadratic Envelope | $C_{\text{quad}} = \sup_{t>0} H_D(t)/t^2$ | $\le 26.8497$ | PASS |
| **G1.1** $D_0(T)$ Tails | $C_{D_0}(T \ge 30)$ / $C_{D_0}(T \ge 13)$ | $\le 4.34617$ / $\le 7.72465$ | PASS |
| **G1.2** $k \ne 0$ Shifts | $S_{k \ne 0}$ ($c \ge 30$ / $c \le 29$) | $\le 2.8563/c$ / $\le 0.21727$ | PASS |
| **G1.2** Subcase B1 ($d=1, X \ge 0.88$) | $S^{\text{Gauss}} + S_{k \ne 0}$ | $\le 1.25848$ (at $X \in [10.22, 10.84]$) | PASS |
| **G1.2** Subcase B1 ($d \ge 2, X \ge 2.0$) | $S^{\text{Gauss}} + S_{k \ne 0}$ | $\le 1.13689$ (at $d=3, X \ge 120$) | PASS |
| **G1.2** Subcase B2a ($X > 6\pi$) | $S^{\text{Gauss}} + S_{k \ne 0}$ | $\le 1.25112$ (at $d=1, t_1 \in [1.15, 1.20], \tilde{D}_1=2, a_1=1$) | PASS |
| **G1.2** Subcase B2b ($X > 6\pi$) | $S^{\text{Gauss}} + S_{k \ne 0}$ | $\le 1.23337$ (at $d=1, t_1 \in [0.65, 0.70], \tilde{D}_1=4, a_1 \ge 31$) | PASS |
| **G1.2** Case A ($t_0 > 1.50, X > 6\pi$) | $S^{\text{Gauss}} + S_{k \ne 0}$ | $\le 0.915136$ (at $s=1, a_1 \ge 10$) | PASS |
| **G1.3** $m$-Space Finite-$n$ Difference | Discrete $c \in \{2..19\}$ ($\kappa \ge 0.88/c^2$) & $\{20..39\}$ ($T_1 \ge 7$) | $\le 0.277972$ (at $c=5, \kappa \in [2.427, 2.50]$) | PASS |
| **G1.3** $m$-Space Zero-Loop | $c \ge 20$ ($T_1 \in [0.05, 7]$) & $c \ge 40$ ($T_1 \in [7, 50]$) | $\le 0.276524$ (at $T_1 \in [0.05, 0.054], c=20$) | PASS |
| **G1.3** $m$-Space Analytic Tails | $T_1 \le 0.05$ ($c \ge 20$) / $T_1 \ge 50$ ($c \ge 60$) / dead | $\le 0.232621$ / $\le 0.212014$ / $\le 0.00156907$ | PASS |
| **G1.3** Global $S_{m\text{-space}}$ | $\sup_{\text{Piece G1}} |\Phi^{\text{ex}} - \Phi^G|$ | $\le 0.277972$ | PASS |
| **G1.4** Grand Pure Phase Bound | $\sup_{\text{Piece G1}} |\Phi^{\text{ex}}(\phi; 1)|$ | **$\le 1.536448 < \frac{\pi}{2} = 1.57080$** | PASS |
| **G1.4** Pure Integral Lower Bound | $I_{\text{pure}} / L_1(1) \ge \cos(\Phi_{\max})$ | **$\ge +0.03434128$** | PASS |
| **G1.5** Cubic-Root Perturbation | $\Delta_{123}$ / $C_{\text{osc}} \ll C_{\text{exp,pert}}$ / $\text{bdy\_ratio}$ | $\le 5.1501 \times 10^{-6}$ / $16.913 \ll 37.5613$ / $\le 0.00077761$ | PASS |
| **G1.5** **FINAL FULL INTEGRAL** | **$I_{\text{full}} / L_1(1) \ge \cos(\Phi_{\max}) - \Delta_{123} - (1+\Delta_{123})\,\text{bdy\_ratio}$** | **$\ge +0.03355852 > 0$** | **PASS** |
| **Part B** Finite-$n$ Cross-Check | All 27 exact Minor-Arc corners | $\max|\Phi^{\text{ex}}| = 1.1542 < 1.53645$, $\min I_{\text{pure}}/L_1 = +0.88956$ | **27/27 PASS** |

---

## 3. Key Mathematical Upgrades Introduced for `#78`

1. **Symmetry Normalization ($\alpha \in (0, 1/2] \implies a_1 \ge 2$):**
   Since $D_r = c\,\mathrm{dist}(r j/c, \mathbb{Z})$ is identically invariant under $j/c \mapsto (c-j)/c$, replacing $\{j/c\}$ by $\alpha := \min(\{j/c\}, 1-\{j/c\}) \in (0, 1/2]$ guarantees that the standard continued fraction $\alpha = [0; a_1, \dots, a_s]$ (with $a_s \ge 2$) satisfies $a_1 = \lfloor 1/\alpha \rfloor \ge 2$. Consequently:
   - $s = 1 \implies q = a_s$, $c = d a_s$, $\mathrm{REST}(0) = 0$, no earlier levels, and $\mathrm{SEMI}(0)$ sums only $u = 1, \dots, \lfloor (a_s-1)/2 \rfloor$ (empty for $a_s = 2$).
   - $s \ge 2 \implies q_{s-1} \ge q_1 = a_1 \ge 2$ and $q = a_s q_{s-1} + q_{s-2} \ge 2 a_s + 1$.
   - For any level $k$ with $a_{k+1} = 1$: necessarily $k \ge 1$ ($q_k \ge 2$) and $k \le s - 2$ ($D_k \ge 2d$); if $a_k = 1$ then $k \ge 2$, $q_k \ge 3$, and $\theta_k = q_{k-1}/q_k > 1/2$ strictly.

2. **1-Line $D_r^2$-Cancellation Theorem for High Levels $k \ge k_2 + 1$ (Eliminates Subsequent Branch-4 Exclusion):**
   For *any* integer $r \ge 1$ with $D_r > 0$, substituting $t_r = \frac{r c \kappa}{2\pi D_r}$ cancels $D_r^2$ identically:
   $$\frac{1}{4\pi^2} \frac{c}{r D_r^2} H_D(t_r) = \frac{t_r^2 H_D(t_r)}{c \kappa^2 r^3} \le \frac{U_{2H}}{c \kappa^2 r^3}.$$
   Using $\frac{1}{c \kappa^2} = \frac{M_{k_2} q_{k_2}^3}{4\pi^2 t_{k_2}^2}$ and summing over $r \in \mathrm{CONV}(k) \cup \mathrm{SEMI}(k)$ gives, across **all 4 minor-arc branches simultaneously**:
   $$\sum_{k \ge k_2 + 1} \big(\mathrm{CONV}(k) + \mathrm{SEMI}(k)\big) \le 2 C_0 \frac{U_{2H}}{t_{k_2}^2} \frac{S_{\text{fib3}}(r_q)}{D_{k_2}},$$
   with certified 50-digit bounds $S_{\text{fib3}}(r_q \ge 2) \le 0.51776$ and $S_{\text{fib3}}(r_q \ge 1) \le 2.34517$. Applying the same identity directly from level $k_1$ in Subcases B2a/B2b yields the unconditional consecutive-level cap
   $$\mathrm{CONV}(k_2) + \mathrm{SEMI}(k_2) \le 2 C_0 U_{2H} \frac{M_1}{t_1^2 (q_{k_2}/q_{k_1})^3} \le 2 C_0 U_{2H} \frac{a_1 + 2}{D_1 t_1^2 (a_1 + \theta_{k_1})^3},$$
   summed over $k \ge k_2$ in B2b with $S_{\text{fib3}, 0} = \sum_{m \ge 0} F_{m+1}^{-3} \le 2.17259$.

3. **Exact Fibonacci-Duality Unconditional `REST(k)` Bound:**
   For $k_2 \le s - 2$ with depth $L := (s - 1) - k_2 \ge 1$, backward Fibonacci monotonicity gives $D_{k_2 + j} \ge F_{L - j + 2} d$ ($0 \le j \le L$) while forward Fibonacci growth gives $t_{k_2 + j} \ge F_{j+1}^2 t_{k_2}$ (or $F_{j+2} F_{j+1} t_{k_2}$ when $a_{k_2+1} \ge 2$) and $t_{s-1} \ge F_{L+1} F_{L+2} t_{k_2}$ with $a_s \ge 2$. Taking the maximum over exact depths $L \in \{1..5\}$ and analytic $L \ge 6$ (with the $O(\log t / t)$ tail `single_hi_t_rest_iv` for $j \ge 6$) bounds $\sum_{k \ge k_2} \mathrm{REST}(k)$ **unconditionally** (zero `high_kap` assumption).

4. **Zero-Cancellation Centered-Form $D_0(T)$ & Paired $\theta_+ + \theta_- = 1$ `S_SetI`:**
   - `D0_ub_exact_iv` evaluates $D_0(T_m)$ via midpoint quadrature with exact $\frac{h^3}{24}\sup|f''|$ error plus $\frac{T_b - T_a}{2}\max(g_1^{\text{ub}}(T_a), g_2^{\text{ub}}(T_a))$ with exact substitution tail $w = y + y^2/2$, eliminating $O(h)$ cancellation bloat.
   - `S_SetI_ub_iv` couples the opposite-sign residues $r_- = q - r_+$ ($\theta_+ + \theta_- = 1$, $\theta \in [0, 1/2]$) for $\nu \le 6$ (taking $\max$ with the $q \mid \nu$ case $\theta = 0$ for $\nu \ge 2$; $q \nmid 1$ since $q \ge 2$).

5. **Exact Subcase B1 Identities & Subcase B2 Coupling:**
   - Exact $M_{s-1} t_{s-1} = \frac{X}{2\pi d}$ (for $X \ge 120$, Branch 3 immediately forces $t_{s-1} \le \frac{C_\tau d}{240\pi}$ and $a_s \ge \lfloor \frac{120}{2\pi t_{b3}} \rfloor$).
   - Exact $t_{r_u} = \frac{X}{2\pi u} - t_{s-1}$ at $k = s - 1$, and exact $\theta_{\text{hi}} \implies a_{s-1} \ge \lfloor 1/\theta_{\text{hi}} \rfloor$ in `B1_earlier_iv`.
   - Exact 5-step $V_1 := C_0 M_1 H_D(t_1)$ ladder in B2, coupled 3-branch $(k_1 = 0 \text{ vs } a_{k_1} \ge 2 \text{ vs } a_{k_1} = 1)$ split with exact integer Branch-1 feasibility checks, exact $D_{r_u} = (u + 1/\tilde{D}_1) D_1$ in B2a `SEMI`, and coupled $a_{k_1 - 1} = 1$ ($\theta_{k_1} < 2/3$) vs $a_{k_1 - 1} \ge 2$ ($D_{k_1 - 2} \ge (3\tilde{D}_1 + 2)d$) in `B2_earlier_split_iv`.

6. **Universal Unit-Circle Cubic-Root Perturbation (`Lemma G1.5`, $I_{\text{pure}} \to I_{\text{full}}$):**
   Valid across **all** of Regime 2 ($n \ge 132$, $\kappa = \kappa_{\text{trig}} = 2\pi \frac{\sin\alpha_1 \sin\psi_j}{\sin\theta_c} < 2.5$) with **zero** reliance on $c\kappa^2 \ge 42$:
   - Uses exact $\gamma_{\text{in}} = \gamma_{\text{ext}}$, $n\kappa_{\text{trig}} \ge 15\pi^2 S_\alpha^2 \implies T_{\text{lg2}} \ge 8.88518$, and the 1-line Riemann-sum inequality $\sum_{m \ge 1} \log\coth(am/2) < \frac{\pi^2}{4a}$ for strictly decreasing positive integrands.
   - Certifies $\Delta_{123} \le 5.1501 \times 10^{-6}$, $C_{\text{osc,tot}} \le 16.913 \ll C_{\text{exp,pert}} \ge 37.5613$, and $\mathrm{bdy\_ratio} \le 0.00077761$, yielding $I_{\text{full}} / L_1(1) \ge +0.03355852 > 0$.

---

## 4. Complete 7-Item `#78` Audit Resolution Checklist

| `#78` Item | Issue Reported | Exact Resolution in `algebra/certify_G1_minor.py` |
|---|---|---|
| **1** | Missing explicit cubic-root perturbation `I_pure -> I_full` across Piece G1 | Added **Lemma G1.5** (L2825–L2915) in pure 50-digit `iv` using universal Regime-2 oscillation + Rouché bounds (exact $\kappa = \kappa_{\text{trig}}$, no $c\kappa^2 \ge 42$): certifies `I_full / L1(1) >= +0.03355852 > 0`. |
| **2** | `S_q(d, X)` non-monotonicity for $d \ge 2, X \in [2.0, 2.6]$ and `X in (30, 100)` / $d \ge 2$ gap | Added `S_q_all_d_ge_iv(d_min, Xa, Xb)` (L542–L556) checking every discrete $d$ with $d X_a/2 < 13$ + monotone `C_D0_13`/`C_D0_30` tail, and `S_q_largeX_iv(d_min, X_lo)` (L558–L577) covering all $X \ge X_{\text{lo}}$ and all $d \ge d_{\min}$ with zero gaps. |
| **3** | `M_lo = a/d_min` falsely dropping Branch 1 & Branch 3 when $d$ or $D_k$ is unbounded | `B1_eval_from_precomp(..., exact_d=False)`, B1 `X >= 120` ($d \ge 13$), B2 Part (i) ($d \ge 5$), B2 Part (ii) ($\tilde{D}_1 \ge 12$), and `hl_cs_ge2` (`D_eff_b1 = max(Dk, ceil(a/10))`) **never** drop Branch 1 or Branch 3 when $d$ or $D_k$ is unbounded. |
| **4** | CF inequalities: `(a+1)/D_k -> (a+2)/D_k` for $k < s-1$, and `t2_b` (`t_2/t_1 >= a_1`, not `1.5*a_1`) | Fixed `M_hi` to exact `(a+2)/D_k` in `Phi_tab_iv`, `Phi_amin_tab`, `hl_cs_ge2`, and B2 (with exact `(a + theta_hi)/D_1 + 1/(Dt*D_1)` in B2a where $D_{k_2}=d$). Fixed `t2_b` to exact `(Dt/(Dt-1)) * qr * ta` where `qr = a_1` (or `a_1 + 0.5` only when $a_{k_1} = 1$, where $\theta_{k_1} > 1/2$ strictly). |
| **5** | Missing unconditional `SEMI+REST` at subsequent Branch-4 levels | Replaced `hl_noB4` + `B4_later` completely by the 1-line $D_r^2$-cancellation theorem (covers `CONV + SEMI` across all 4 branches simultaneously for all $k \ge k_2 + 1$) + Fibonacci-duality unconditional `REST(k)` (L1164–L1425). |
| **6** | Lemma G1.3 `c >= 16` unclamped scaling in `drp` & outer tail | Split discrete $c \in \{2..19\}$ ($\kappa \ge 0.88/c^2$) + $\{20..39\}$ ($T_1 \ge 7$) vs zero-loop $c \ge 20$ ($T_1 \in [0.05, 7]$) & $c \ge 40$ ($T_1 \in [7, 50]$); uses unclamped `cout_err_uncl = max(cout_lo_err, cout_h - 1)` in `drp` & outer tail, with `assert cin_l.a > 0.992 and cin_h.b < 1.050` confirming strict monotonicity in $c$. |
| **7** | Hardcoded constants (`1.72`, `0.90`, `0.99`, L1681, L1796 peak, edge intervals, `assert tail_hit`, `2.6`/`2.5`) | Certified `M_3half_iv <= 0.883016`, `M_2_iv <= 1.11625`, `C_semi_sum_max < 2.6`, `C_rest_sum_max < 2.5`, `C_rest_m2_max = 4` in `iv`; replaced `0.99` with exact `cout_l_12`; fixed `>=`/`<=` edge intervals + self-test in `H_ub_iv`/`tH_ub_iv`; added `assert tail_hit`, separate inner/outer peaks, explicit discrete peak `g(15)` in $T_1 \ge 50$, and outward-rounded `make_geom_grid_iv`. |
