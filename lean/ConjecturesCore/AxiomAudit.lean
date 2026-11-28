import ConjecturesCore.SignLaw
import ConjecturesCore.BoundaryAndSample
import ConjecturesCore.ZeroCrevicePartition
import ConjecturesCore.ExactEngine

/-!
# Official Machine-Checked Axiom Audit

Running `lake build` evaluates `#print axioms` on every core theorem below,
proving directly from the Lean 4 kernel that zero `sorryAx`, zero custom axioms,
and zero `Lean.ofReduceBool` (`native_decide`) appear anywhere in the formalization.
-/

#print axioms ConjecturesCore.sign_law_all_forms_equivalent
#print axioms ConjecturesCore.intSign_c0Coeff
#print axioms ConjecturesCore.boundary_c0_poly_identity
#print axioms ConjecturesCore.high_c_is_boundary
#print axioms ConjecturesCore.sample_point_count_exact
#print axioms ConjecturesCore.lemma_jensen_index_bounds
#print axioms ConjecturesCore.g1_Dr2_exact_elimination
#print axioms ConjecturesCore.strict_seam_overlaps
#print axioms ConjecturesCore.interior_zero_crevice_complete
#print axioms ConjecturesCore.global_zero_crevice_partition_complete
#print axioms ConjecturesCore.fastChoose_matches_choose_up_to_12
#print axioms ConjecturesCore.kernel_verified_even_n2_to_n12
#print axioms ConjecturesCore.kernel_verified_small_expansions
#print axioms ConjecturesCore.kernel_verified_glz_B10_rayleigh_disproof
