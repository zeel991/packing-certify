import math
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache

import mpmath as mp
import numpy as np
import pytest

from packing_certify import certify, exact, exploit, polish
from packing_certify.geometry import Problem

# --- exact verification -------------------------------------------------------------------------


def test_exact_accepts_touching_and_rejects_a_tiny_overlap():
    p = Problem("square", 2)
    q = Fraction(1, 4)
    touching = [q, q, q, 3 * q, q, q]  # centres 1/2 apart, radii 1/4: they touch exactly
    assert exact.verify(p, touching).valid
    overlap = touching[:5] + [q + Fraction(1, 10**30)]
    rep = exact.verify(p, overlap)
    assert not rep.valid and "circles 0 and 1 overlap" in rep.violations


def test_published_constructions_are_exactly_valid():
    for name in exact.names():
        p, v, _ = exact.load(name)
        assert exact.verify(p, v).valid, name


# --- geometry -----------------------------------------------------------------------------------


@pytest.mark.parametrize("kind", ["square", "rectangle"])
def test_gradients_match_finite_differences(kind):
    rng = np.random.default_rng(0)
    p = Problem(kind, 4)
    v = list(rng.uniform(0.2, 0.8, p.nvar))
    h = 1e-7
    for k in p.constraints():
        g = p.gradient(k, v, math.sqrt)
        for c in range(p.nvar):
            w = v[:]
            w[c] += h
            fd = (p.value(k, w, math.sqrt) - p.value(k, v, math.sqrt)) / h
            assert abs(fd - g.get(c, 0)) < 1e-5, (k, c)


# --- certificate --------------------------------------------------------------------------------


@lru_cache
def _aemon():
    p, v, _ = exact.load("square26_aemon")
    vf, _ = polish.converge(p, v)
    act = polish.active_set(p, vf)
    return p, polish.polish(p, vf, act.keys), act


def test_aemon_optimum_is_certified():
    p, pol, act = _aemon()
    c = certify.certify(p, pol.v, act.keys)
    assert c.certified and c.krawczyk_margin < 1e-6
    assert c.value_lo.startswith("2.6359830849176077831865694854")
    assert c.min_multiplier_lo > 0.01


def test_certificate_rejects_a_wrong_objective(monkeypatch):
    p, pol, act = _aemon()
    monkeypatch.setattr(Problem, "objective_gradient", lambda self: [-1 if c % 3 == 2 else 0 for c in range(self.nvar)])
    c = certify.certify(p, pol.v, act.keys)
    assert not c.certified and "multiplier" in c.reason


def test_certificate_rejects_a_wrong_contact_set():
    p, pol, act = _aemon()
    wrong = [k for k in p.constraints() if k not in set(act.keys)][0]
    c = certify.certify(p, pol.v, act.keys[:-1] + [wrong])
    assert not c.certified


def test_certificate_rejects_a_point_outside_its_box():
    p, pol, act = _aemon()
    with mp.workdps(pol.dps):
        moved = [x + mp.mpf("1e-20") for x in pol.v]  # far outside the 1e-40 box
    c = certify.certify(p, moved, act.keys)
    assert not c.certified and "Krawczyk" in c.reason


def test_alphaevolve_26_polishes_to_a_different_lower_optimum():
    p, v, _ = exact.load("square26_alphaevolve")
    vf, _ = polish.converge(p, v)
    act = polish.active_set(p, vf)
    c = certify.certify(p, polish.polish(p, vf, act.keys).v, act.keys)
    assert c.certified and Decimal(c.value_hi) < Decimal(_aemon_value())


def _aemon_value() -> str:
    p, pol, act = _aemon()
    return certify.certify(p, pol.v, act.keys).value_lo


# --- verifiers and exploits ---------------------------------------------------------------------


def test_reproduced_deepmind_verifier_accepts_the_published_data():
    for name in exact.names():
        p, v, meta = exact.load(name)
        circles = np.array([[float(x) for x in c] for c in meta["circles"]])
        check = exploit.deepmind_square if p.kind == "square" else exploit.deepmind_rectangle
        assert check(circles), name


def test_tolerance_evaluator_accepts_an_invalid_packing_above_the_certified_maximum():
    p, pol, act = _aemon()
    e = exploit.tolerance(p, pol.v)
    cert = certify.certify(p, pol.v, act.keys)
    assert e.accepted and not e.exact_valid
    assert float(e.exact_sum) - float(cert.value_hi) > 1e-5


def test_strict_float_verifier_accepts_an_overlapping_packing():
    p, pol, _ = _aemon()
    e = exploit.float_rounding(p, pol.v)
    assert e.accepted and not e.exact_valid and e.worst_violation < 1e-15
