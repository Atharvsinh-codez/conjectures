import ConjecturesCore.SignLaw
import ConjecturesCore.BoundaryAndSample

/-!
# Proposition 1.3: Zero-Crevice Partition Completeness Theorem

Formal Lean 4 Core proof (zero `sorry`, zero custom axioms) that EVERY even integer
`n = 2m >= 2`, EVERY admissible triple `(a, b, c)` with `a + 2b + 3c = 2m`, EVERY
strictly interior sample index `j ∈ J_{n,c}`, and EVERY positive rational separation
parameter `κ > 0` is covered by at least one of the certified proof regimes of
Table 1, with strict seam overlaps (`0.85 < 0.88 < 0.90` and `40 < 45`) and zero crevices.
-/

namespace ConjecturesCore

/-- Exact strictly positive rational number `num / den` (`num > 0`, `den > 0`). -/
structure PosRat where
  num : Nat
  den : Nat
  num_pos : 0 < num
  den_pos : 0 < den

namespace PosRat

def le (x y : PosRat) : Prop := x.num * y.den ≤ y.num * x.den
def lt (x y : PosRat) : Prop := x.num * y.den < y.num * x.den

def ofPosNat (k : Nat) (hk : 0 < k) : PosRat := ⟨k, 1, hk, by omega⟩

def mul (x y : PosRat) : PosRat :=
  ⟨x.num * y.num, x.den * y.den, Nat.mul_pos x.num_pos y.num_pos, Nat.mul_pos x.den_pos y.den_pos⟩

def sq (x : PosRat) : PosRat := mul x x

def mulNat (k : Nat) (hk : 0 < k) (x : PosRat) : PosRat := mul (ofPosNat k hk) x

theorem lt_or_ge (x y : PosRat) : lt x y ∨ le y x :=
  Nat.lt_or_ge (x.num * y.den) (y.num * x.den)

theorem le_or_gt (x y : PosRat) : le x y ∨ lt y x :=
  (Nat.lt_or_ge (y.num * x.den) (x.num * y.den)).symm

theorem lt_trans {x y z : PosRat} (h1 : lt x y) (h2 : lt y z) : lt x z := by
  unfold lt at h1 h2 ⊢
  have hA : x.num * y.den * z.den < y.num * x.den * z.den :=
    Nat.mul_lt_mul_of_pos_right h1 z.den_pos
  have hB : y.num * z.den * x.den < z.num * y.den * x.den :=
    Nat.mul_lt_mul_of_pos_right h2 x.den_pos
  have hEq1 : y.num * x.den * z.den = y.num * z.den * x.den := by ac_rfl
  have hEq2 : x.num * y.den * z.den = (x.num * z.den) * y.den := by ac_rfl
  have hEq3 : z.num * y.den * x.den = (z.num * x.den) * y.den := by ac_rfl
  have hC : (x.num * z.den) * y.den < (z.num * x.den) * y.den := by
    rw [← hEq2, ← hEq3]
    exact Nat.lt_trans hA (hEq1 ▸ hB)
  exact Nat.lt_of_mul_lt_mul_right hC

theorem lt_imp_le {x y : PosRat} (h : lt x y) : le x y :=
  Nat.le_of_lt h

end PosRat

/-- Exact rational thresholds from Table 1 & Proposition 1.3. -/
def K_REG1 : PosRat := ⟨5, 2, by omega, by omega⟩        -- 2.5
def C_2A   : PosRat := ⟨42, 1, by omega, by omega⟩       -- 42
def C_2C   : PosRat := ⟨9, 10, by omega, by omega⟩       -- 0.90
def C_G2LO : PosRat := ⟨17, 20, by omega, by omega⟩      -- 0.85
def C_G1D1 : PosRat := ⟨22, 25, by omega, by omega⟩      -- 0.88
def K_G2Q2 : PosRat := ⟨2, 1, by omega, by omega⟩        -- 2.0
def R_G3K  : PosRat := ⟨9, 20, by omega, by omega⟩       -- 0.45
def R_G3CK : PosRat := ⟨60, 1, by omega, by omega⟩       -- 60
def C_TAU_OUTER : PosRat := ⟨40, 1, by omega, by omega⟩  -- 40 (Piece G3 outer threshold)
def C_TAU_INNER : PosRat := ⟨45, 1, by omega, by omega⟩  -- 45 (Piece G1/G3 overlap seam)

/-- Strict seam overlaps certified in Lean 4: `0.85 < 0.88 < 0.90` and `40 < 45`. -/
theorem strict_seam_overlaps :
    C_G2LO.lt C_G1D1 ∧ C_G1D1.lt C_2C ∧ C_TAU_OUTER.lt C_TAU_INNER := by
  refine ⟨?_, ?_, ?_⟩ <;> simp [PosRat.lt, C_G2LO, C_G1D1, C_2C, C_TAU_OUTER, C_TAU_INNER]

/-- Inner Piece G3 existence predicate with threshold `C_tau` (`45` for inner, `40` for outer). -/
def HasG3Approx (c j : Nat) (hc : 0 < c) (κ Ctau : PosRat) : Prop :=
  ∃ (p0 r0 : Nat) (hr0 : 0 < r0),
    Nat.gcd p0 r0 = 1 ∧
    ∃ (hDne : (r0 : Int) * (j : Int) - (c : Int) * (p0 : Int) ≠ 0),
    let Dabs := Int.natAbs ((r0 : Int) * (j : Int) - (c : Int) * (p0 : Int))
    10 * r0 * Dabs * Dabs ≤ c ∧
    (PosRat.mulNat r0 hr0 κ).le R_G3K ∧
    (PosRat.mulNat (r0 * r0 * r0 * c)
      (Nat.mul_pos (Nat.mul_pos (Nat.mul_pos hr0 hr0) hr0) hc) κ.sq).le R_G3CK ∧
    let hD : 0 < Dabs := Int.natAbs_pos.mpr hDne
    let hD4 : 0 < Dabs * Dabs * Dabs * Dabs := Nat.mul_pos (Nat.mul_pos (Nat.mul_pos hD hD) hD) hD
    (PosRat.mulNat (Dabs * Dabs * Dabs * Dabs) hD4 Ctau).le
      (PosRat.mulNat (c * c * c * r0) (Nat.mul_pos (Nat.mul_pos (Nat.mul_pos hc hc) hc) hr0) κ.sq)

/-- Inner Piece G3 (`C_tau = 45`) strictly implies outer Piece G3 (`C_tau = 40`). -/
theorem inner_G3_implies_outer_G3 {c j : Nat} {hc : 0 < c} {κ : PosRat}
    (h : HasG3Approx c j hc κ C_TAU_INNER) :
    HasG3Approx c j hc κ C_TAU_OUTER := by
  rcases h with ⟨p0, r0, hr0, hgcd, hDne, h1, h2, h3, h4⟩
  refine ⟨p0, r0, hr0, hgcd, hDne, h1, h2, h3, ?_⟩
  simp only [PosRat.le, PosRat.mulNat, PosRat.mul, PosRat.ofPosNat, C_TAU_OUTER, C_TAU_INNER] at h4 ⊢
  exact Nat.le_trans (Nat.mul_le_mul_right _ (Nat.mul_le_mul_left _ (by omega))) h4

/-- Exhaustive certified regimes for an interior point `(n, c, j, κ)` (`n >= 132`, `c >= 2`). -/
inductive InteriorRegime (c j : Nat) (hc : 2 ≤ c) (κ : PosRat) : Prop where
  | regime1 : K_REG1.le κ → InteriorRegime c j hc κ
  | sub2A   : κ.lt K_REG1 → C_2A.le (PosRat.mulNat c (by omega) κ.sq) → InteriorRegime c j hc κ
  | sub2C   :
      ∀ (hq : 0 < c / Nat.gcd c j),
      κ.lt K_REG1 →
      (PosRat.mulNat c (by omega) κ.sq).lt C_2A →
      (PosRat.mulNat (c * (c / Nat.gcd c j)) (Nat.mul_pos (by omega) hq) κ).le C_2C →
      InteriorRegime c j hc κ
  | pieceG2 :
      ∀ (hq : 0 < c / Nat.gcd c j),
      κ.lt K_REG1 →
      (PosRat.mulNat c (by omega) κ.sq).lt C_2A →
      2 ≤ Nat.gcd c j →
      C_G2LO.lt (PosRat.mulNat (c * (c / Nat.gcd c j)) (Nat.mul_pos (by omega) hq) κ) →
      (c / Nat.gcd c j = 1 ∨
        (PosRat.mulNat ((c / Nat.gcd c j) * (c / Nat.gcd c j)) (Nat.mul_pos hq hq) κ).le K_G2Q2) →
      InteriorRegime c j hc κ
  | pieceG3 :
      κ.lt K_REG1 →
      (PosRat.mulNat c (by omega) κ.sq).lt C_2A →
      HasG3Approx c j (by omega) κ C_TAU_OUTER →
      InteriorRegime c j hc κ
  | pieceG1 :
      ∀ (hq : 0 < c / Nat.gcd c j),
      κ.lt K_REG1 →
      (PosRat.mulNat c (by omega) κ.sq).lt C_2A →
      2 ≤ c / Nat.gcd c j →
      ((Nat.gcd c j = 1 ∧
          C_G1D1.le (PosRat.mulNat ((c / Nat.gcd c j) * (c / Nat.gcd c j)) (Nat.mul_pos hq hq) κ)) ∨
       (2 ≤ Nat.gcd c j ∧
          K_G2Q2.le (PosRat.mulNat ((c / Nat.gcd c j) * (c / Nat.gcd c j)) (Nat.mul_pos hq hq) κ))) →
      ¬ HasG3Approx c j (by omega) κ C_TAU_INNER →
      InteriorRegime c j hc κ

/-- Helper: `c >= 2` implies `d = gcd c j >= 1` and `q = c / d >= 1` with `d * q = c`. -/
theorem gcd_div_pos (c j : Nat) (hc : 2 ≤ c) :
    1 ≤ Nat.gcd c j ∧ 1 ≤ c / Nat.gcd c j ∧ Nat.gcd c j * (c / Nat.gcd c j) = c := by
  have hd_pos : 0 < Nat.gcd c j := Nat.gcd_pos_of_pos_left j (by omega)
  have hd_dvd : Nat.gcd c j ∣ c := Nat.gcd_dvd_left c j
  have hmul : Nat.gcd c j * (c / Nat.gcd c j) = c := Nat.mul_div_cancel' hd_dvd
  have hq_pos : 0 < c / Nat.gcd c j := by
    match hq : c / Nat.gcd c j with
    | 0 =>
      rw [hq, Nat.mul_zero] at hmul
      omega
    | _ + 1 => omega
  exact ⟨hd_pos, hq_pos, hmul⟩

/-- PROPOSITION 1.3 (Interior Zero-Crevice Partition Completeness — ALL `c >= 2`, `j`, `κ > 0`):
    Every interior sample point falls into at least one of `{Regime1, 2A, 2C, G2, G3, G1}`. -/
theorem interior_zero_crevice_complete (c j : Nat) (hc : 2 ≤ c) (κ : PosRat) :
    InteriorRegime c j hc κ := by
  rcases PosRat.lt_or_ge κ K_REG1 with hκ_lt | hκ_ge
  · rcases PosRat.lt_or_ge (PosRat.mulNat c (by omega) κ.sq) C_2A with h2A_lt | h2A_ge
    · have ⟨hd_pos, hq_pos, hdq⟩ := gcd_div_pos c j hc
      rcases PosRat.le_or_gt
          (PosRat.mulNat (c * (c / Nat.gcd c j)) (Nat.mul_pos (by omega) hq_pos) κ) C_2C with h2C_le | h2C_gt
      · exact InteriorRegime.sub2C hq_pos hκ_lt h2A_lt h2C_le
      · have h_gt_88 : C_G1D1.lt (PosRat.mulNat (c * (c / Nat.gcd c j)) (Nat.mul_pos (by omega) hq_pos) κ) :=
          PosRat.lt_trans strict_seam_overlaps.2.1 h2C_gt
        have h_gt_85 : C_G2LO.lt (PosRat.mulNat (c * (c / Nat.gcd c j)) (Nat.mul_pos (by omega) hq_pos) κ) :=
          PosRat.lt_trans strict_seam_overlaps.1 h_gt_88
        by_cases hq1 : c / Nat.gcd c j = 1
        · have hd2 : 2 ≤ Nat.gcd c j := by
            rw [hq1, Nat.mul_one] at hdq
            omega
          exact InteriorRegime.pieceG2 hq_pos hκ_lt h2A_lt hd2 h_gt_85 (Or.inl hq1)
        · have hq2 : 2 ≤ c / Nat.gcd c j := by omega
          by_cases hd1 : Nat.gcd c j = 1
          · have hqc : c / Nat.gcd c j = c := by
              rw [hd1]
              omega
            have hX_ge : C_G1D1.le
                (PosRat.mulNat ((c / Nat.gcd c j) * (c / Nat.gcd c j)) (Nat.mul_pos hq_pos hq_pos) κ) := by
              have hlt : C_G1D1.lt
                  (PosRat.mulNat ((c / Nat.gcd c j) * (c / Nat.gcd c j)) (Nat.mul_pos hq_pos hq_pos) κ) := by
                simpa [hqc] using h_gt_88
              exact PosRat.lt_imp_le hlt
            by_cases hG3 : HasG3Approx c j (by omega) κ C_TAU_INNER
            · exact InteriorRegime.pieceG3 hκ_lt h2A_lt (inner_G3_implies_outer_G3 hG3)
            · exact InteriorRegime.pieceG1 hq_pos hκ_lt h2A_lt hq2 (Or.inl ⟨hd1, hX_ge⟩) hG3
          · have hd2 : 2 ≤ Nat.gcd c j := by omega
            rcases PosRat.le_or_gt
                (PosRat.mulNat ((c / Nat.gcd c j) * (c / Nat.gcd c j)) (Nat.mul_pos hq_pos hq_pos) κ) K_G2Q2 with hX_le | hX_gt
            · exact InteriorRegime.pieceG2 hq_pos hκ_lt h2A_lt hd2 h_gt_85 (Or.inr hX_le)
            · by_cases hG3 : HasG3Approx c j (by omega) κ C_TAU_INNER
              · exact InteriorRegime.pieceG3 hκ_lt h2A_lt (inner_G3_implies_outer_G3 hG3)
              · exact InteriorRegime.pieceG1 hq_pos hκ_lt h2A_lt hq2 (Or.inr ⟨hd2, PosRat.lt_imp_le hX_gt⟩) hG3
    · exact InteriorRegime.sub2A hκ_lt h2A_ge
  · exact InteriorRegime.regime1 hκ_ge

/-- Global master partition covering ALL even integers `n = 2m >= 2` and ALL admissible `(a, b, c)`:
    either Finite Base (`n <= 130`), or one of the 4 Boundaries (`c=0, c=1, a=0, b=0`),
    or an Interior Slice (`n >= 132`, `2 <= c <= (n-3)/3`) where every `(j, κ)` is covered
    by `interior_zero_crevice_complete`. -/
inductive GlobalCoverage (m a b c : Nat) : Prop where
  | finiteBase : 2 * m ≤ 130 → GlobalCoverage m a b c
  | boundary   : (c = 0 ∨ c = 1 ∨ a = 0 ∨ b = 0) → GlobalCoverage m a b c
  | interior   :
      132 ≤ 2 * m →
      (hc : 2 ≤ c) →
      c ≤ (2 * m - 3) / 3 →
      (∀ (j : Nat) (κ : PosRat), InteriorRegime c j hc κ) →
      GlobalCoverage m a b c

/-- PROPOSITION 1.3 (Global Zero-Crevice Partition Completeness Theorem — ALL even `n = 2m >= 2`):
    Zero crevices across the entire infinite parameter space. -/
theorem global_zero_crevice_partition_complete {m a b c : Nat}
    (hm : 1 ≤ m) (hadm : Admissible (2 * m) a b c) :
    GlobalCoverage m a b c := by
  by_cases hn : 2 * m ≤ 130
  · exact GlobalCoverage.finiteBase hn
  · have hn132 : 132 ≤ 2 * m := by omega
    by_cases hbdy : c = 0 ∨ c = 1 ∨ a = 0 ∨ b = 0
    · exact GlobalCoverage.boundary hbdy
    · have hc2 : 2 ≤ c := by omega
      by_cases hchi : c ≤ (2 * m - 3) / 3
      · exact GlobalCoverage.interior hn132 hc2 hchi
          (fun j κ => interior_zero_crevice_complete c j hc2 κ)
      · have hab : a = 0 ∨ b = 0 := high_c_is_boundary hm hadm (by omega)
        rcases hab with ha0 | hb0
        · exact GlobalCoverage.boundary (Or.inr (Or.inr (Or.inl ha0)))
        · exact GlobalCoverage.boundary (Or.inr (Or.inr (Or.inr hb0)))

end ConjecturesCore
