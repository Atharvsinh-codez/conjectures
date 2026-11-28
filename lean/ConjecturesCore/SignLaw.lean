/-!
# Universal Sign-Law Equivalence (Bardakov-Kozlovskaya Conjecture 3.3(2))

Formal Lean 4 Core verification (zero `sorry`, zero custom axioms) of the exact
algebraic equivalence of all four formulations of the Jadav Sign Law across every
even integer `n = 2m >= 2` and every admissible triple `(a, b, c)` with
`a + 2b + 3c = 2m` and partition indexing `k1 = a + b + c`, `k2 = b + c`, `k3 = c`.
-/

namespace ConjecturesCore

/-- Exact sign function on integers: `+1`, `-1`, or `0`. -/
def intSign (x : Int) : Int :=
  if x > 0 then 1 else if x < 0 then -1 else 0

/-- Exact power `(-1)^k` as an integer. -/
def signPow (k : Nat) : Int :=
  if k % 2 == 0 then 1 else -1

/-- Admissibility condition `a + 2b + 3c = n` for monomial `e1^a e2^b e3^c`. -/
def Admissible (n a b c : Nat) : Prop :=
  a + 2 * b + 3 * c = n

/-- For every even degree `n = 2m`, admissibility forces `a ≡ c (mod 2)`. -/
theorem parity_a_eq_c {m a b c : Nat} (h : Admissible (2 * m) a b c) :
    a % 2 = c % 2 := by
  unfold Admissible at h
  omega

/-- `signPow k` is never zero. -/
theorem signPow_ne_zero (k : Nat) : signPow k ≠ 0 := by
  unfold signPow
  split <;> decide

/-- `signPow` depends only on `k % 2`. -/
theorem signPow_eq_of_mod_eq {k1 k2 : Nat} (h : k1 % 2 = k2 % 2) :
    signPow k1 = signPow k2 := by
  unfold signPow
  rw [h]

/-- Exact parity identity between `b*(c+1) + c` and `b*(a+1) + c` whenever `a ≡ c (mod 2)`. -/
theorem mod2_bc_eq_ba {a b c : Nat} (hac : a % 2 = c % 2) :
    (b * (c + 1) + c) % 2 = (b * (a + 1) + c) % 2 := by
  have h1 : (b * (c + 1)) % 2 = ((b % 2) * ((c + 1) % 2)) % 2 := Nat.mul_mod b (c + 1) 2
  have h2 : (b * (a + 1)) % 2 = ((b % 2) * ((a + 1) % 2)) % 2 := Nat.mul_mod b (a + 1) 2
  have h3 : (c + 1) % 2 = (a + 1) % 2 := by omega
  have h4 : (b * (c + 1)) % 2 = (b * (a + 1)) % 2 := by rw [h1, h2, h3]
  omega

/-- Exact parity identity linking `(a, b, c)` indexing to partition indexing
    `k1 = a + b + c`, `k2 = b + c`, `k3 = c`:
    `k1*k2 + k1*k3 + k2*k3 ≡ b*(c + 1) + c (mod 2)` whenever `a ≡ c (mod 2)`. -/
theorem mod2_partition_eq_bc {a b c : Nat} (hac : a % 2 = c % 2) :
    let k1 := a + b + c
    let k2 := b + c
    let k3 := c
    (k1 * k2 + k1 * k3 + k2 * k3) % 2 = (b * (c + 1) + c) % 2 := by
  dsimp only
  have h1 : ((a + b + c) * (b + c)) % 2 = (((a + b + c) % 2) * ((b + c) % 2)) % 2 :=
    Nat.mul_mod (a + b + c) (b + c) 2
  have h2 : ((a + b + c) * c) % 2 = (((a + b + c) % 2) * (c % 2)) % 2 :=
    Nat.mul_mod (a + b + c) c 2
  have h3 : ((b + c) * c) % 2 = (((b + c) % 2) * (c % 2)) % 2 :=
    Nat.mul_mod (b + c) c 2
  have h4 : (b * (c + 1)) % 2 = ((b % 2) * ((c + 1) % 2)) % 2 :=
    Nat.mul_mod b (c + 1) 2
  have hmodA : a % 2 = 0 ∨ a % 2 = 1 := by omega
  have hmodB : b % 2 = 0 ∨ b % 2 = 1 := by omega
  have hmodC : c % 2 = 0 ∨ c % 2 = 1 := by omega
  rcases hmodA with hA | hA <;> rcases hmodB with hB | hB <;> rcases hmodC with hC | hC
  all_goals
    first
    | omega
    | (rw [Nat.add_mod, Nat.add_mod ((a + b + c) * (b + c)), h1, h2, h3, Nat.add_mod (b * (c + 1)), h4]
       simp [Nat.add_mod, hA, hB, hC])

/-- Master equivalence of all four forms of the Universal Sign Law (Theorem 1.1(1)). -/
theorem sign_law_all_forms_equivalent {m a b c : Nat} (h : Admissible (2 * m) a b c) :
    let k1 := a + b + c
    let k2 := b + c
    let k3 := c
    signPow (b * (c + 1) + c) = signPow (b * (a + 1) + c) ∧
    signPow (b * (c + 1) + c) = signPow (k1 * k2 + k1 * k3 + k2 * k3) ∧
    signPow (b * (c + 1) + c) = (if c % 2 == 0 then signPow b else -1) := by
  have hac : a % 2 = c % 2 := parity_a_eq_c h
  refine ⟨?_, ?_, ?_⟩
  · exact signPow_eq_of_mod_eq (mod2_bc_eq_ba hac)
  · exact signPow_eq_of_mod_eq (mod2_partition_eq_bc hac).symm
  · have hmodC : c % 2 = 0 ∨ c % 2 = 1 := by omega
    rcases hmodC with hC | hC
    · have h1 : (b * (c + 1) + c) % 2 = b % 2 := by
        rw [Nat.add_mod, Nat.mul_mod]
        have : (c + 1) % 2 = 1 := by omega
        rw [this, hC]
        omega
      simp only [hC, beq_self_eq_true, ↓reduceIte]
      exact signPow_eq_of_mod_eq h1
    · have h1 : (b * (c + 1) + c) % 2 = 1 := by
        rw [Nat.add_mod, Nat.mul_mod]
        have : (c + 1) % 2 = 0 := by omega
        rw [this, hC]
        omega
      simp only [hC, show (1 == 0) = false from rfl, Bool.false_eq_true, ↓reduceIte]
      unfold signPow
      rw [h1]
      rfl

end ConjecturesCore
