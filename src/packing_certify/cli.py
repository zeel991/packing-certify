"""packing-certify: exact verification, high-precision polishing, rigorous local-optimality certificates
and verifier exploits for the public AlphaEvolve circle-packing constructions."""

from __future__ import annotations

import argparse
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="packing-certify", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("report", help="run every stage on every construction (writes reports/)")
    r.add_argument("--out", type=Path, default=Path("reports"))
    c = sub.add_parser("check", help="exact check, polish and certificate for one construction")
    c.add_argument("name", help="e.g. square26_aemon")
    args = p.parse_args(argv)

    if args.cmd == "report":
        from . import report

        report.build(args.out)
        print(f"Wrote {args.out / 'REPORT.md'}, results.json and certified-optima.png")
        return 0

    from . import certify, exact, polish

    prob, v, meta = exact.load(args.name)
    rep = exact.verify(prob, v)
    print(f"published: exactly valid={rep.valid}, sum={float(rep.total):.15f}")
    vf, _ = polish.converge(prob, v)
    act = polish.active_set(prob, vf)
    pol = polish.polish(prob, vf, act.keys)
    cert = certify.certify(prob, pol.v, act.keys)
    print(f"active constraints {len(act.keys)}/{prob.nvar}, separated by {act.separation_orders:.1f} orders")
    print(f"certificate: {cert.reason}; value in [{cert.value_lo}, {cert.value_hi}]")
    return 0 if cert.certified else 1


if __name__ == "__main__":
    raise SystemExit(main())
