"""From a published construction to a high-precision local optimum.

1. converge: SLSQP in float64 from the published point (published points are often not fully
   converged; AlphaEvolve's are written with 8-decimal centres).
2. active_set: the constraints that are touching, chosen at the largest gap in the sorted slacks so
   the choice is unambiguous (the report shows the separation).
3. polish: Newton's method at 60 digits on "every active constraint holds with equality". When the
   active constraints are as many as the variables this is a square system with a unique solution.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import mpmath as mp
import numpy as np
from scipy.optimize import minimize

from .geometry import Problem


def _f(v):
    return [float(x) for x in v]


def converge(p: Problem, v0) -> tuple[np.ndarray, int]:
    cons = p.constraints()
    x0 = np.array(_f(v0))

    def g(v):
        return np.array([p.value(k, v, math.sqrt) for k in cons])

    def jac(v):
        J = np.zeros((len(cons), p.nvar))
        for r, k in enumerate(cons):
            for c, val in p.gradient(k, v, math.sqrt).items():
                J[r, c] = val
        return J

    grad = -np.array(p.objective_gradient(), dtype=float)
    res = minimize(
        lambda v: -float(p.total(v)), x0, jac=lambda v: grad, method="SLSQP",
        constraints=[{"type": "ineq", "fun": g, "jac": jac}],
        options={"ftol": 1e-16, "maxiter": 2000},
    )
    return res.x, int(res.nit)


@dataclass
class ActiveSet:
    keys: list
    largest_active_slack: float
    smallest_inactive_slack: float

    @property
    def separation_orders(self) -> float:
        return math.log10(self.smallest_inactive_slack / max(self.largest_active_slack, 1e-14))


def active_set(p: Problem, v, ceiling: float = 1e-5) -> ActiveSet:
    """Split constraints at the largest jump (in orders of magnitude) among slacks below `ceiling`."""
    slacks = sorted((abs(p.value(k, v, math.sqrt)), k) for k in p.constraints())
    # Slacks below float64 resolution are all "touching": floor them so 0.0 vs 1e-15 is not a jump.
    logs = [math.log10(max(s, 1e-14)) for s, _ in slacks]
    cut, best = None, -1.0
    for i in range(len(slacks) - 1):
        if slacks[i][0] > ceiling:
            break
        jump = logs[i + 1] - logs[i]
        if jump > best:
            best, cut = jump, i + 1
    keys = [k for _, k in slacks[:cut]]
    return ActiveSet(keys, slacks[cut - 1][0], slacks[cut][0])


def _jacobian(p: Problem, v, keys) -> mp.matrix:
    J = mp.matrix(len(keys), p.nvar)
    for r, k in enumerate(keys):
        for c, val in p.gradient(k, v, mp.sqrt).items():
            J[r, c] = val
    return J


@dataclass
class Polished:
    v: list  # mpmath values
    residual: mp.mpf
    total: mp.mpf
    multipliers: list
    dps: int


def polish(p: Problem, v_float, keys, dps: int = 60, steps: int = 40) -> Polished:
    if len(keys) != p.nvar:
        raise ValueError(f"{len(keys)} active constraints for {p.nvar} variables: not a square system")
    with mp.workdps(dps):
        v = [mp.mpf(x) for x in _f(v_float)]
        tol = mp.mpf(10) ** (-(dps - 8))
        for _ in range(steps):
            g = mp.matrix([p.value(k, v, mp.sqrt) for k in keys])
            res = max(abs(x) for x in g)
            if res < tol:
                break
            step = mp.lu_solve(_jacobian(p, v, keys), -g)
            v = [v[c] + step[c] for c in range(p.nvar)]
        g = mp.matrix([p.value(k, v, mp.sqrt) for k in keys])
        lam = mp.lu_solve(_jacobian(p, v, keys).T, -mp.matrix(p.objective_gradient()))
        return Polished(v, max(abs(x) for x in g), p.total(v), [lam[i] for i in range(len(keys))], dps)
