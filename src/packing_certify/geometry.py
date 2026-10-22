"""The two AlphaEvolve packing problems as constrained optimisation problems.

Variables v: (x_i, y_i, r_i) for each circle, plus the rectangle width W for the perimeter-4 problem.
Objective: maximise sum r_i. Constraints g_k(v) >= 0:
  - disjoint:  dist(c_i, c_j) - r_i - r_j
  - left/bottom walls: x_i - r_i, y_i - r_i
  - right/top walls:   W - x_i - r_i, H - y_i - r_i, with W = H = 1 (square) or H = 2 - W (rectangle)

Every function takes the number type's own sqrt, so the same code runs on floats, mpmath and Arb balls.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from fractions import Fraction
from typing import Callable, Sequence

Constraint = tuple  # ("pair", i, j) or ("wall", i, side) with side in 0..3


@dataclass(frozen=True)
class Problem:
    kind: str  # "square" | "rectangle"
    n: int

    @property
    def nvar(self) -> int:
        return 3 * self.n + (1 if self.kind == "rectangle" else 0)

    def constraints(self) -> list[Constraint]:
        out: list[Constraint] = [("pair", i, j) for i, j in itertools.combinations(range(self.n), 2)]
        out += [("wall", i, s) for i in range(self.n) for s in range(4)]
        return out

    def box(self, v: Sequence):
        """(W, H) of the container in the same number type as v."""
        if self.kind == "square":
            one = v[0] * 0 + 1
            return one, one
        w = v[3 * self.n]
        return w, 2 - w

    def value(self, k: Constraint, v: Sequence, sqrt: Callable):
        if k[0] == "pair":
            _, i, j = k
            dx, dy = v[3 * i] - v[3 * j], v[3 * i + 1] - v[3 * j + 1]
            return sqrt(dx * dx + dy * dy) - v[3 * i + 2] - v[3 * j + 2]
        _, i, s = k
        x, y, r = v[3 * i], v[3 * i + 1], v[3 * i + 2]
        W, H = self.box(v)
        return (x - r, y - r, W - x - r, H - y - r)[s]

    def gradient(self, k: Constraint, v: Sequence, sqrt: Callable) -> dict[int, object]:
        """Non-zero partial derivatives of constraint k (sparse row of the Jacobian)."""
        if k[0] == "pair":
            _, i, j = k
            dx, dy = v[3 * i] - v[3 * j], v[3 * i + 1] - v[3 * j + 1]
            d = sqrt(dx * dx + dy * dy)
            return {3 * i: dx / d, 3 * i + 1: dy / d, 3 * j: -dx / d, 3 * j + 1: -dy / d, 3 * i + 2: -1, 3 * j + 2: -1}
        _, i, s = k
        row = {
            0: {3 * i: 1, 3 * i + 2: -1},
            1: {3 * i + 1: 1, 3 * i + 2: -1},
            2: {3 * i: -1, 3 * i + 2: -1},
            3: {3 * i + 1: -1, 3 * i + 2: -1},
        }[s]
        if self.kind == "rectangle" and s in (2, 3):
            row[3 * self.n] = 1 if s == 2 else -1
        return row

    def objective_gradient(self) -> list[int]:
        return [1 if c < 3 * self.n and c % 3 == 2 else 0 for c in range(self.nvar)]

    def total(self, v: Sequence):
        return sum(v[3 * i + 2] for i in range(self.n))


def from_circles(kind: str, circles: Sequence[Sequence[Fraction]]) -> tuple[Problem, list[Fraction]]:
    """Variables for a published construction. Rectangle constructions are translated so the
    bounding box starts at the origin; W is then its width (the height is checked separately)."""
    p = Problem(kind, len(circles))
    if kind == "square":
        return p, [c for circle in circles for c in circle]
    x0 = min(x - r for x, y, r in circles)
    y0 = min(y - r for x, y, r in circles)
    moved = [(x - x0, y - y0, r) for x, y, r in circles]
    width = max(x + r for x, y, r in moved)
    return p, [c for circle in moved for c in circle] + [width]


def circles_of(p: Problem, v: Sequence) -> list[tuple]:
    return [(v[3 * i], v[3 * i + 1], v[3 * i + 2]) for i in range(p.n)]
