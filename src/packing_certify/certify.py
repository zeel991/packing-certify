"""Rigorous local-optimality certificate in ball (interval) arithmetic, using Arb via python-flint.

Given a polished point x~ and an active set with as many constraints as variables, it proves:

1. Existence and uniqueness (Krawczyk). Let B = x~ + [-eps, eps]^m, G the active constraints and
   Y an approximate inverse of G'(x~). If K(B) = x~ - Y G(x~) + (I - Y G'(B)) (B - x~) lies strictly
   inside B, then G has exactly one zero x* in B, and every matrix in G'(B) is nonsingular.
2. Strict complementarity. The multipliers solving G'(x*)^T lambda = -grad f are enclosed by solving
   with the interval matrix G'(B); every lower bound must be > 0.
3. Feasibility. Every other constraint has a positive lower bound over all of B.

Then x* is feasible, the active gradients are linearly independent, and the multipliers are strictly
positive with m active constraints in m variables. So the only direction that keeps all active
constraints satisfied to first order is zero, and x* is a strict local maximum of sum(r). Its value
is enclosed by sum(r) over B.
"""

from __future__ import annotations

from dataclasses import dataclass

import mpmath as mp
from flint import arb, arb_mat, ctx

from .geometry import Problem

SQRT = lambda t: t.sqrt()  # noqa: E731


@dataclass
class Certificate:
    certified: bool
    reason: str
    eps: float
    value_lo: str
    value_hi: str
    krawczyk_margin: float  # max |K_i - x~_i| / eps; must be < 1
    min_multiplier_lo: float
    min_inactive_lo: float


def _arb(x: mp.mpf) -> arb:
    return arb(mp.nstr(x, 90, min_fixed=-mp.inf, max_fixed=mp.inf)).mid()


def certify(p: Problem, v: list, keys: list, eps: str = "1e-40", prec_bits: int = 320) -> Certificate:
    m = p.nvar
    if len(keys) != m:
        return Certificate(False, f"{len(keys)} active constraints for {m} variables", float(eps), "", "", 0, 0, 0)
    old = ctx.prec
    ctx.prec = prec_bits
    try:
        e = arb(eps)
        xt = [_arb(x) for x in v]
        B = [arb(x, e) for x in xt]

        G0 = arb_mat([[p.value(k, xt, SQRT)] for k in keys])

        def jac(at) -> arb_mat:
            rows = []
            for k in keys:
                row = [arb(0)] * m
                for c, val in p.gradient(k, at, SQRT).items():
                    row[c] = arb(val) if isinstance(val, int) else val
                rows.append(row)
            return arb_mat(rows)

        Jt = jac(xt)
        Y = Jt.inv()
        Y = arb_mat([[Y[i, j].mid() for j in range(m)] for i in range(m)])  # any approximate inverse works
        JB = jac(B)
        I = arb_mat([[1 if i == j else 0 for j in range(m)] for i in range(m)])
        delta = arb_mat([[arb(0, e)] for _ in range(m)])
        K = (-(Y * G0)) + (I - Y * JB) * delta  # this is K(B) - x~
        margin = max(float(K[i, 0].abs_upper() / e) for i in range(m))
        if not margin < 1:
            return Certificate(False, "Krawczyk test failed (K(B) not inside B)", float(eps), "", "", margin, 0, 0)

        grad_f = arb_mat([[-arb(gi)] for gi in p.objective_gradient()])
        lam = JB.transpose().solve(grad_f)
        lam_lo = min(float(lam[i, 0].lower()) for i in range(m))

        active = set(keys)
        inactive_lo = min(float(p.value(k, B, SQRT).lower()) for k in p.constraints() if k not in active)

        total = p.total(B)
        ok = lam_lo > 0 and inactive_lo > 0
        reason = "strict local maximum" if ok else ("a multiplier is not provably positive" if lam_lo <= 0 else "an inactive constraint is not provably positive")
        return Certificate(
            ok, reason, float(eps),
            total.lower().str(40, radius=False), total.upper().str(40, radius=False),
            margin, lam_lo, inactive_lo,
        )
    finally:
        ctx.prec = old
