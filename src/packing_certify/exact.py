"""Exact verification in rational arithmetic: no floating point, no tolerance, no square roots.

Each published number is read as the exact decimal it is written as. Two circles are disjoint iff
d^2 >= (r_i + r_j)^2 (both sides non-negative), and every wall constraint is linear, so the whole
check is exact over the rationals.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import mpmath as mp

from .geometry import Problem, from_circles

DATA = Path(__file__).resolve().parents[2] / "data" / "constructions.json"


def load(name: str, path: Path = DATA) -> tuple[Problem, list[Fraction], dict]:
    meta = json.loads(path.read_text())[name]
    circles = [tuple(Fraction(Decimal(s)) for s in c) for c in meta["circles"]]
    p, v = from_circles(meta["problem"], circles)
    return p, v, meta


def names(path: Path = DATA) -> list[str]:
    return list(json.loads(path.read_text()))


@dataclass
class ExactReport:
    valid: bool
    total: Fraction
    violations: list[str] = field(default_factory=list)
    min_pair_gap: float = 0.0  # in length units, computed at 60 digits from the exact values
    min_wall_gap: float = 0.0


def _gap(d2: Fraction, s: Fraction) -> mp.mpf:
    with mp.workdps(60):
        return mp.sqrt(mp.mpf(d2.numerator) / d2.denominator) - mp.mpf(s.numerator) / s.denominator


def verify(p: Problem, v: list[Fraction]) -> ExactReport:
    bad: list[str] = []
    pair_gaps, wall_gaps = [], []
    for k in p.constraints():
        if k[0] == "pair":
            _, i, j = k
            dx, dy = v[3 * i] - v[3 * j], v[3 * i + 1] - v[3 * j + 1]
            s = v[3 * i + 2] + v[3 * j + 2]
            d2 = dx * dx + dy * dy
            if d2 < s * s:
                bad.append(f"circles {i} and {j} overlap")
            pair_gaps.append(_gap(d2, s))
        else:
            g = p.value(k, v, sqrt=None)  # linear, no sqrt needed
            if g < 0:
                bad.append(f"circle {k[1]} crosses wall {k[2]}")
            wall_gaps.append(g)
    for i in range(p.n):
        if v[3 * i + 2] <= 0:
            bad.append(f"circle {i} has non-positive radius")
    return ExactReport(
        valid=not bad,
        total=p.total(v),
        violations=bad,
        min_pair_gap=float(min(pair_gaps)),
        min_wall_gap=float(min(wall_gaps)),
    )
