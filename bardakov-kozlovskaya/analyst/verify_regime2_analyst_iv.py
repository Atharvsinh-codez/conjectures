"""
MASTER 50-DIGIT OUTWARD-ROUNDED INTERVAL (mpmath.iv) VERIFICATION RUNNER
FOR ALL @analyst REGIME 2 COMPONENTS OF BARDAKOV-KOZLOVSKAYA CONJECTURE 3.3(2):
  1. Lemma R2.1 (Universal Regime 2 Cubic-Root Perturbation) + Lemma R2.2 (Sub-regime 2A, c*kap^2 >= 42)
  2. Lemma R2.3 (Sub-regime 2C-Upgrade, c*q*kap <= 0.90, ALL q >= 1, g >= 1)
  3. Lemma R2.5 Piece G3 Cases G3-A1, G3-A2, G3-B, G3-C (zero-hardcoded 50-digit iv)
  4. Lemma R2.5 Piece G3 Case G3-D (|D| = 1, tau in [1.05, 40.0], 50-digit iv)
  5. 32-Corner Finite-n Empirical Cross-Check across ALL 5 Piece G3 sub-cases (up to n ~ 10^8)

Exits with code 0 if and only if EVERY 50-digit interval assertion and empirical check passes.
"""
import subprocess
import sys
import time
from pathlib import Path

BASE = Path(__file__).resolve().parent
SCRIPTS = [
    ("R2.1 + R2.2 (Sub-regime 2A, c*kap^2 >= 42)",      "certify_R2_1_and_2A_iv.py"),
    ("R2.3 (Sub-regime 2C-Upgrade, c*q*kap <= 0.90)",   "certify_R2_3_iv.py"),
    ("R2.5 Piece G3 (Cases A1, A2, B, C)",              "certify_G3_ABC_full_iv.py"),
    ("R2.5 Piece G3 (Case D, |D|=1 tau in [1.05, 40])", "certify_G3_D_full_iv.py"),
    ("32-Corner Finite-n Empirical Cross-Check (G3)",   "verify_G3_empirical_corners.py"),
]

t_all = time.time()
print("==========================================================================")
print("MASTER 50-DIGIT mpmath.iv VERIFICATION: ALL @analyst REGIME 2 PROOFS")
print("==========================================================================")
for name, fname in SCRIPTS:
    t0 = time.time()
    print(f"\n>>> Running: {name} ({fname}) ...")
    res = subprocess.run([sys.executable, "-u", str(BASE / fname)], check=False)
    dt = time.time() - t0
    if res.returncode != 0:
        print(f"[FAIL] {name} exited with code {res.returncode} after {dt:.1f}s!")
        sys.exit(res.returncode)
    print(f"[PASS] {name} completed cleanly in {dt:.1f}s.")

print("\n==========================================================================")
print(f"ALL {len(SCRIPTS)} @analyst REGIME 2 VERIFICATION SUITES PASSED (EXIT CODE 0)!")
print(f"Total elapsed: {time.time() - t_all:.1f}s")
print("==========================================================================")
