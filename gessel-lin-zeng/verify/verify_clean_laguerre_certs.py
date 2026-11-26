#!/usr/bin/env python3
"""Zero-dependency exact integer/Fraction verification of the clean Ma-Wang
Theorem 2 sharpness certificates (Laguerre inequality violations R'(x)^2 - R(x)R''(x) < 0)
across all remaining boundary families i=4, i=5, and i=k-3.
"""
from fractions import Fraction


def S_L(p: list[int], L: int) -> list[int]:
    """Apply S_L(p)(x) = (1 + L*x)*p(x) + x*(1 - x)*p'(x) in exact integer arithmetic."""
    d = len(p) - 1
    q = [0] * (d + 2)
    for j, c in enumerate(p):
        q[j] += (j + 1) * c
        q[j + 1] += (L - j) * c
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q


def A_kS(k: int, deleted_bars) -> list[int]:
    """Exact descent polynomial A_{k,S}(x) via Ma-Wang eqs. (8)-(10)."""
    S = set(deleted_bars)
    p = [1]
    L = 0
    for j in range(1, k + 1):
        if j not in S:
            p = S_L(p, L)
            L += 1
        p = S_L(p, L)
        L += 2
    return p


def eval_laguerre(p: list[int], x) -> tuple:
    v = sum(c * (x ** j) for j, c in enumerate(p))
    d1 = sum(j * c * (x ** (j - 1)) for j, c in enumerate(p) if j >= 1)
    d2 = sum(j * (j - 1) * c * (x ** (j - 2)) for j, c in enumerate(p) if j >= 2)
    return v, d1, d2, d1 * d1 - v * d2


def main() -> None:
    cases = [
        ("i=4, k=31 (first fail)", 31, range(28, 32), range(1, 5), 932, -23),
        ("i=5, k=16 (first fail)", 16, range(12, 17), range(1, 6), 269, -11),
        ("i=k-3, k=28 (first fail)", 28, range(4, 29), range(1, 26), 272, Fraction(-73, 10)),
        ("i=k-3, k=29 (integer pt)", 29, range(4, 30), range(1, 27), 349, -8),
    ]
    for label, k, s_late, s_early, lam, x0 in cases:
        f = A_kS(k, s_late)
        g = A_kS(k, s_early)
        assert len(f) == len(g), (label, len(f), len(g))
        R = [lam * a + b for a, b in zip(f, g)]
        _, _, _, lag = eval_laguerre(R, x0)
        assert lag < 0, f"Expected negative Laguerre form for {label}"
        print(f"PASS {label}: deg={len(R)-1}, lambda={lam}, x={x0}, "
              f"Laguerre ~= {float(lag):.6e} < 0 (exact)")


if __name__ == "__main__":
    main()
