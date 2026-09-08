#!/usr/bin/env python3
"""M0 multireference self-gate: the <S**2>-based pre-flight check.

Parses unrestricted-DFT outputs (Gaussian-format spin-annihilation lines) and
applies the M0 pass/flag rule of the spin-aware AutoMeKin extension:

  PASS   : annihilation LOWERS <S**2> toward S(S+1)  AND
           |<S**2>_annihilated - S(S+1)| <= PASS_TOL          (per functional)
  FLAG   : annihilation RAISES <S**2>  OR
           |<S**2>_annihilated - S(S+1)| >= FLAG_TOL          (per functional)
  REVIEW : anything between (reported for the operator, never silently classified)

Aggregate verdict over the functional cross-check set:
  FLAG requires the per-functional flag to persist across >= MIN_FLAG functionals
  (a single functional's pathology cannot trigger a refusal); PASS requires every
  functional to pass; everything else is REVIEW.

Thresholds calibrated on the populations observed in this work: every passing
state sits within 0.4 of ideal after annihilation (most within 0.005); every
refused state sits >= 2.5 above ideal with annihilation raising.

Usage:
  python m0_gate.py --mult 3 triplet_ub3lyp.log triplet_opbe.log triplet_tpssh.log
  python m0_gate.py --mult 3 --json gate_verdict.json logs/*.log
"""
import argparse
import json
import re
import sys

PASS_TOL = 0.5      # |<S**2>_ann - S(S+1)| for a per-functional PASS
FLAG_TOL = 1.0      # |<S**2>_ann - S(S+1)| at/above which a functional FLAGs
MIN_FLAG = 2        # functionals that must flag for an aggregate FLAG

ANNIH_RE = re.compile(
    r"S\*\*2 before annihilation\s+([-\d.]+),?\s+after\s+([-\d.]+)")


def parse_annihilation(path):
    """Return (before, after) from the LAST annihilation line of a log."""
    hits = ANNIH_RE.findall(open(path, errors="replace").read())
    if not hits:
        raise ValueError(f"no 'S**2 before annihilation' line in {path}")
    before, after = map(float, hits[-1])
    return before, after


def classify(before, after, ideal):
    """Per-functional verdict: PASS / FLAG / REVIEW."""
    dev_after = abs(after - ideal)
    lowers = abs(after - ideal) < abs(before - ideal)
    if lowers and dev_after <= PASS_TOL:
        return "PASS"
    if (not lowers) or dev_after >= FLAG_TOL:
        return "FLAG"
    return "REVIEW"


def gate(paths, mult):
    spin = (mult - 1) / 2.0
    ideal = spin * (spin + 1.0)
    rows = []
    for p in paths:
        before, after = parse_annihilation(p)
        rows.append({
            "log": p, "s2_before": before, "s2_after": after,
            "ideal": ideal, "dev_after": round(abs(after - ideal), 4),
            "annihilation": "lowers" if abs(after - ideal) < abs(before - ideal)
                            else "raises",
            "verdict": classify(before, after, ideal),
        })
    n_flag = sum(r["verdict"] == "FLAG" for r in rows)
    if n_flag >= MIN_FLAG:
        overall = "FLAG (multireference — single-determinant machinery withheld)"
    elif all(r["verdict"] == "PASS" for r in rows):
        overall = "PASS (single-reference-clean)"
    else:
        overall = "REVIEW (mixed signals — operator decision required)"
    return {"multiplicity": mult, "ideal_s2": ideal,
            "thresholds": {"pass_tol": PASS_TOL, "flag_tol": FLAG_TOL,
                           "min_flag_functionals": MIN_FLAG},
            "functionals": rows, "overall": overall}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("logs", nargs="+",
                    help="unrestricted-DFT logs, one per functional")
    ap.add_argument("--mult", type=int, required=True,
                    help="spin multiplicity 2S+1 of the state under test")
    ap.add_argument("--json", metavar="FILE",
                    help="also write the verdict record as JSON")
    args = ap.parse_args(argv)
    result = gate(args.logs, args.mult)
    for r in result["functionals"]:
        print(f"{r['log']}: <S**2> {r['s2_before']:.4f} -> {r['s2_after']:.4f} "
              f"(ideal {r['ideal']:.2f}, dev {r['dev_after']:.4f}, "
              f"annihilation {r['annihilation']}) : {r['verdict']}")
    print(f"M0 verdict: {result['overall']}")
    if args.json:
        json.dump(result, open(args.json, "w"), indent=2)
    return 0 if result["overall"].startswith("PASS") else 1


if __name__ == "__main__":
    sys.exit(main())
