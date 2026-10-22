# packing-certify: proof-checking the public AlphaEvolve circle-packing records

LLM-driven search systems (AlphaEvolve, and the open-source and commercial systems that followed) report new records on circle-packing problems. A record is only as trustworthy as the check behind it, and the published checks use floating point, sometimes with a 1e-6 tolerance.

This tool takes any published construction and, without floating point where it matters:

1. **Verifies it exactly.** Every published decimal is read as the exact rational it is written as. Disjointness becomes `d² ≥ (rᵢ + rⱼ)²` over the rationals, with no square roots and no tolerance.
2. **Finds the local optimum it sits next to and polishes it to 60 digits.** A local optimiser (SLSQP) reaches the optimum from the published point. The touching constraints are identified where the gap between touching and free constraints is unambiguous, and Newton's method then solves "touching constraints hold with equality" at 60 digits.
3. **Proves that optimum is a strict local maximum, rigorously.** In Arb ball arithmetic:
   - a Krawczyk test proves a unique exact optimum lies in a box of radius 1e-40;
   - every Lagrange multiplier has a positive lower bound;
   - every non-touching constraint stays positive over the whole box.

   The value comes out enclosed in an interval about 5e-39 wide.
4. **Attacks the published verifiers.** From each certified optimum it grows radii for as long as a verifier still accepts, then checks the result exactly.

## Results (`reports/REPORT.md`)

All four public constructions: AlphaEvolve's B.12 (n=26, n=32) and B.13 (n=21, perimeter-4 rectangle), plus Aemon's n=26 record.

| Construction | Published sum | Valid, checked exactly | Certified strict local maximum next to it |
|---|---|---|---|
| AlphaEvolve, n=26 | 2.635862756413698 | yes | 2.635977394754392752218631205141… |
| Aemon, n=26 | 2.635983082278162 | yes | 2.635983084917607783186569485443… |
| AlphaEvolve, n=32 | 2.937944526205518 | yes | 2.939572771206323634996322323460… |
| AlphaEvolve, rectangle n=21 | 2.365832133416763 | yes | 2.365832375910832741572079722291… |

![Certified local optima with their contact graphs](reports/certified-optima.png)

What the certificates show:

- **All four published constructions are valid in exact arithmetic.** In each, every circle is pinned: touching constraints equal the variables (78/78, 78/78, 96/96, 64/64), and the next constraint is at least 10 orders of magnitude further away.
- **Aemon's n=26 packing is a different and strictly better local optimum than AlphaEvolve's.** It's better even after AlphaEvolve's configuration is fully polished, by 5.69e-6. Most of the headline gap between the two published numbers (1.203e-4) is polishing, though: 1.146e-4 of it comes from locally optimising AlphaEvolve's own published configuration. Aemon's published value sits 2.6e-9 below its exact optimum, which looks like a deliberate safety margin.
- **AlphaEvolve's published n=32 construction is 1.63e-3 below the local optimum next to it.** Polishing alone reaches 2.939572771206…, the size of the later "2.939+" listed for n=32 on [Erich Friedman's packing pages](https://erich-friedman.github.io/packing/cirRsqu/). Whether that listed result has the same configuration isn't checked here.

What the published verifiers accept:

| Verifier | Accepts an invalid packing? | Score it can fake |
|---|---|---|
| DeepMind's notebook (strict float64 comparison) | Yes: 11–22 overlapping pairs of about 4e-17 | None measurable: totals stay below the certified maximum |
| OpenEvolve example / CodeEvolve problem (1e-6 tolerance) | Yes: every contact overlaps by up to 1e-6 | **+1.3e-5 (n=26), +1.6e-5 (n=32)**. From Aemon's optimum it scores 2.6359961, above every published n=26 value. |

The tolerance verifiers are reproduced line for line from pinned commits; `exploit.py` cites the exact files and lines. **This is not a claim that any published result used the tolerance.** It shows that a search optimising against such an evaluator is rewarded for overlaps worth about 1e-5, the same size as recent record margins.

## What has been checked

- **12 tests pass**, covering:
  - exact checks on hand-built touching and overlapping pairs (1e-30 overlap rejected);
  - analytic gradients against finite differences for both problem types;
  - the Aemon certificate value;
  - the reproduced DeepMind verifier accepting all published data, as the notebook does;
  - both exploits being accepted yet exactly invalid.
- **Three tests break the certificate on purpose,** and each must come back "not certified":
  - a reversed objective (caught by the multiplier check);
  - a wrong contact set;
  - a polished point moved outside its box (caught by the Krawczyk test).

## Limits

- **Certificates are local.** Each proves a strict local maximum next to a published point, not global optimality.
- **The optimisation path matters.** The "polished" optimum for an under-converged construction is the one SLSQP reaches from the published point; a different local method could stop at a different local optimum.
- **Only fully pinned configurations are certified,** meaning as many touching constraints as variables. All four cases here are like that. A configuration with fewer contacts needs a second-order check, which isn't implemented.
- **Arb's interval arithmetic is rigorous, but this code is not formally verified.** The negative tests are there to show the certificate can fail.
- **The rectangle problem has no tolerance exploit,** because no tolerance-based rectangle verifier was reproduced.

## Run it

```bash
uv sync
uv run packing-certify report                 # all constructions: reports/REPORT.md, results.json, certified-optima.png (about 7 s)
uv run packing-certify check square26_aemon   # one construction
uv run pytest
```

Data: AlphaEvolve constructions from [google-deepmind/alphaevolve_results](https://github.com/google-deepmind/alphaevolve_results) at commit 4226acb (CC BY 4.0, Google LLC). Aemon's construction is as published on [aemon.ai](https://aemon.ai/) for public verification, retrieved 2026-10-22. Both are stored as their original decimal strings in [data/constructions.json](data/constructions.json).
