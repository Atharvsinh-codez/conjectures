# Lean 4 Core Formalization (`ConjecturesCore`)

**Zero external dependencies (`leanprover/lean4:v4.36.0` pure Core) · Zero `sorry` · Zero custom axioms · Zero `native_decide`**

This directory contains a standalone, machine-checked **Lean 4 Core** formalization of the algebraic, combinatorial, partition-completeness, and exact-integer verification backbone behind our resolution of **Bardakov–Kozlovskaya Conjecture 3.3(2)** (and the exact Rayleigh disproof certificate for **Gessel–Lin–Zeng Conjecture 15**).

Because `ConjecturesCore` uses **pure Lean 4 Core with zero Mathlib download**, anyone with Lean 4 (`elan`) installed can clone this repository and kernel-verify every theorem in seconds.

## Quick Verification (< 10 seconds)

```bash
cd lean
lake build
```

Running `lake build` compiles all five modules and executes [`ConjecturesCore/AxiomAudit.lean`](ConjecturesCore/AxiomAudit.lean), which runs `#print axioms` on every exported theorem directly in the Lean 4 kernel to prove that no theorem depends on `sorryAx`, `Lean.ofReduceBool` (`native_decide`), or any custom axiom (only Lean's standard core axioms `propext`, `Classical.choice`, `Quot.sound` where classical case splits are used, and zero axioms at all for the pure constructive/computational theorems).

## Module & Theorem Inventory

| Module | Key Kernel-Verified Theorems | Mathematical Role in the Proof |
| :--- | :--- | :--- |
| [`ConjecturesCore/SignLaw.lean`](ConjecturesCore/SignLaw.lean) | `parity_a_eq_c`, `mod2_bc_eq_ba`, `mod2_partition_eq_bc`, `sign_law_all_forms_equivalent` | Proves for **all** even $n = 2m$ and admissible $(a,b,c)$ that $a \equiv c \pmod 2$ and all three forms of the Jadav Sign Law (`(-1)^(b(c+1)+c)`, `(-1)^(b(a+1)+c)`, `(-1)^(b+c+bc)`) are identically equivalent to the piecewise parity law. |
| [`ConjecturesCore/BoundaryAndSample.lean`](ConjecturesCore/BoundaryAndSample.lean) | `intSign_c0Coeff`, `boundary_c0_poly_identity`, `high_c_is_boundary`, `sample_point_count_exact`, `lemma_jensen_index_bounds`, `g1_Dr2_exact_elimination` | Proves for **all** even $n = 2m \ge 2$: the closed-form $c = 0$ boundary polynomial identity $Q_0^{(2m)}(1, u) = (1 - 4u)^m$ and its exact sign $(-1)^b$; the high-$c$ boundary lemma ($c > \lfloor(2m-3)/3\rfloor \implies a = 0 \lor b = 0$); exact interior sample-point count $|J_{n,c}| = D_c + 1$; Lemma 4.3(1) index separation bounds; and the Piece G1 exact $D_r^2$-cancellation algebraic identity. |
| [`ConjecturesCore/ZeroCrevicePartition.lean`](ConjecturesCore/ZeroCrevicePartition.lean) | `strict_seam_overlaps`, `inner_G3_implies_outer_G3`, `interior_zero_crevice_complete`, `global_zero_crevice_partition_complete` | Proves **Proposition 1.3 (Global Zero-Crevice Partition Completeness)** across **all infinitely many** even $n = 2m \ge 2$, admissible triples $(a,b,c)$, interior indices $j$, and positive rational parameters $\kappa > 0$: strict overlaps $0.85 < 0.88 < 0.90$ and $40 < 45$, and zero uncovered crevices between Finite Base, the 4 Boundaries, Regime 1, Sub-regimes 2A/2C, and Pieces G2/G3/G1. |
| [`ConjecturesCore/ExactEngine.lean`](ConjecturesCore/ExactEngine.lean) | `fastChoose_matches_choose_up_to_12`, `kernel_verified_even_n2_to_n12`, `kernel_verified_small_expansions`, `kernel_verified_glz_B10_rayleigh_disproof` | Pure structurally recursive exact integer Newton power-sum + `(s1, s2, z) -> (e1, e2, e3)` peeling engine inside Lean 4 Core (`decide` / `rfl`, zero `native_decide`): kernel-verifies every admissible coefficient of $p_2, p_4, p_6, p_8, p_{10}, p_{12}$, the explicit Eq. (1.3) expansions of $p_2, p_4, p_6$, and the exact integer GLZ $B_{10}$ Rayleigh disproof numerator/denominator. |
| [`ConjecturesCore/AxiomAudit.lean`](ConjecturesCore/AxiomAudit.lean) | `#print axioms` audit across all 14 main theorems | Machine-checked guarantee of zero `sorryAx`, zero custom axioms, and zero `Lean.ofReduceBool`. |

## Honest Scope & Verification Boundary

To keep the formalization 100% transparent and free of exaggerated claims, here is the exact boundary between what is verified by the **Lean 4 kernel** vs. what is verified by **closed-form analytic proofs + 50-digit outward-rounded interval arithmetic (`mpmath.iv`)**:

- **Kernel-verified in Lean 4 Core (`lean/`):**
  - Universal sign-law algebraic equivalence across all even $n = 2m \ge 2$.
  - All-$n$ $c = 0$ boundary theorem, high-$c$ boundary classification, Theorem 4.1 sample-point counts, Lemma 4.3(1) index separation bounds, and the Piece G1 exact $D_r^2$-cancellation identity.
  - **Proposition 1.3 Global Zero-Crevice Partition Completeness** over the entire infinite parameter space $(m, a, b, c, j, \kappa)$ with exact rational arithmetic (`PosRat`).
  - Independent exact integer Newton + peeling engine certified inside the Lean 4 kernel for $n = 2, 4, 6, 8, 10, 12$, explicit Eq. (1.3) expansions, and the GLZ $B_{10}$ integer Rayleigh certificate.
- **Verified in [`bardakov-kozlovskaya/PROOF.md`](../bardakov-kozlovskaya/PROOF.md), [`bardakov-kozlovskaya/paper/main.tex`](../bardakov-kozlovskaya/paper/main.tex), and 50-digit `mpmath.iv` scripts (`bardakov-kozlovskaya/`):**
  - The full exact finite base $n = 2, 4, \dots, 130$ ($33{,}417$ exact integer coefficients via 3 independent Python engines).
  - The transcendental complex-analysis and trigonometric-integral bounds for $a = 0, b = 0, c = 1$, Regime 1 ($\kappa \ge 2.5$), and Regime 2 (Sub-regimes 2A/2C and Pieces G2/G3/G1), certified via closed-form analysis and 50-digit outward-rounded `mpmath.iv` interval arithmetic with locked SHA-256 checksums (see [`bardakov-kozlovskaya/AUDIT.md`](../bardakov-kozlovskaya/AUDIT.md)).
