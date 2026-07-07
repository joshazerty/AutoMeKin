#!/usr/bin/env python3
"""
spin_enumerate.py — Stage M1 of the MECP-aware-node feature (Project 18, Phase 2 Weeks 2+).

PURPOSE
-------
Given the spin-resolved networks AutoMeKin now builds (Week-0 engine + §3 TS/IRC), find the
nodes where a *different* spin state lies close in energy — the candidate spin-crossing
points where an MECP should be localised (Stage M2, locate_mecp.py).

For each node (minimum or TS geometry) carrying a native multiplicity `mult`, this evaluates
the neighbouring spin states (`mult ± 2`, i.e. ΔS = 1) at the SAME geometry (a vertical
spin-gap) using the spin-capable GFN2-xTB engine `XTBamk` from Week 0, and flags the node
when the smallest vertical gap is within a tunable window.

TWO-TIER ACCURACY (ESC-009)
---------------------------
GFN2's Fe spin-gap *energetics* are unreliable, so GFN2 is used only to **flag candidates
cheaply** (a generous window catches anything plausibly crossing). The *decision* and the
crossing energy come from DFT in Stage M2 (easyMECP runs DFT). So set the M1 window WIDE
(default 15 kcal/mol, not 5): better to over-flag cheaply at GFN2 and let DFT/easyMECP reject,
than to miss a crossing because GFN2 mis-ordered the gap. The 5 kcal/mol "close-lying" decision
is applied to the DFT numbers downstream.

INPUT
-----
A set of node geometries (`.xyz` / `.rxyz`), e.g. the per-structure files AutoMeKin writes for
the minima/TSs of a run, or any directory of geometries. Native multiplicity is the run's
`mult` (pass --mult). For the Tp-Fe(II)-ethyl validation you can point it straight at the
`results/V3_betaH/gfn2/{triplet,quintet}/*.xyz` geometries.

OUTPUT
------
`mecp_candidates.json`: one record per flagged node:
  {node, xyz, charge, base_mult, alt_mult, dE_vert_kcal, note}
which Stage M2 (locate_mecp.py) consumes to seed an MECP optimisation between base_mult and
alt_mult at that geometry.

Requires: amk_xtb.py (Week-0 engine) on PYTHONPATH and the `xtb` binary, for the real energy
evaluation. The flagging/IO logic is import-light and unit-tested without xtb.
"""
import os
import sys
import json
import glob
import argparse

HARTREE_KCAL = 627.5094740631


# --------------------------------------------------------------------------- #
# geometry IO (xyz / rxyz; rxyz = xyz with an energy/comment line, same layout)
# --------------------------------------------------------------------------- #
def read_xyz(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    syms, coords = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        syms.append(p[0])
        coords.append((float(p[1]), float(p[2]), float(p[3])))
    return syms, coords


# --------------------------------------------------------------------------- #
# energy evaluation (real one uses XTBamk; injectable for testing)
# --------------------------------------------------------------------------- #
def vertical_energy_xtb(syms, coords, charge, mult, method="gfn2"):
    """GFN2 single-point energy (Eh) at a fixed geometry for the given multiplicity.
    Imported lazily so the flagging logic is testable without xtb/ase installed."""
    from amk_xtb import xtb_sp_grad          # Week-0 engine, functional API
    e_eh, _ = xtb_sp_grad(syms, coords, charge=charge, mult=mult, method=method)
    return e_eh


def alt_multiplicities(base_mult):
    """ΔS = 1 neighbours (mult ± 2), keeping mult >= 1."""
    return [m for m in (base_mult - 2, base_mult + 2) if m >= 1]


def evaluate_node(path, charge, base_mult, method="gfn2", energy_fn=vertical_energy_xtb):
    """Return (base_E_Eh, {alt_mult: alt_E_Eh}) for one node geometry."""
    syms, coords = read_xyz(path)
    base_e = energy_fn(syms, coords, charge, base_mult, method)
    alts = {}
    for m in alt_multiplicities(base_mult):
        try:
            alts[m] = energy_fn(syms, coords, charge, m, method)
        except Exception as e:                # a single failed spin state is non-fatal
            alts[m] = None
            print(f"  [warn] {os.path.basename(path)} mult {m}: {e}", file=sys.stderr)
    return base_e, alts


def flag_candidates(nodes, window_kcal=15.0):
    """nodes: list of dicts {node,xyz,charge,base_mult,base_E,alts:{m:E}}.
    Return the flagged candidate records (smallest |gap| within window)."""
    out = []
    for nd in nodes:
        best = None
        for m, e in nd["alts"].items():
            if e is None:
                continue
            d = (e - nd["base_E"]) * HARTREE_KCAL
            if best is None or abs(d) < abs(best[1]):
                best = (m, d)
        if best is None:
            continue
        alt_mult, dE = best
        if abs(dE) <= window_kcal:
            out.append({
                "node": nd["node"], "xyz": nd["xyz"], "charge": nd["charge"],
                "base_mult": nd["base_mult"], "alt_mult": alt_mult,
                "dE_vert_kcal": round(dE, 2),
                "note": ("alt spin BELOW base (likely crossing)" if dE < 0
                         else "alt spin within window above base"),
            })
    # most-crossing first (most negative gap = alt spin most stabilised)
    out.sort(key=lambda r: r["dE_vert_kcal"])
    return out


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def run(paths, charge, base_mult, window, method, out_json, energy_fn=vertical_energy_xtb):
    nodes = []
    for p in paths:
        base_e, alts = evaluate_node(p, charge, base_mult, method, energy_fn)
        nodes.append({"node": os.path.splitext(os.path.basename(p))[0], "xyz": os.path.abspath(p),
                      "charge": charge, "base_mult": base_mult, "base_E": base_e, "alts": alts})
        gaps = {m: (None if e is None else round((e - base_e) * HARTREE_KCAL, 2))
                for m, e in alts.items()}
        print(f"{os.path.basename(p):30s} base mult {base_mult}: vert gaps (kcal) {gaps}")
    cands = flag_candidates(nodes, window)
    with open(out_json, "w") as fh:
        json.dump({"window_kcal": window, "method": method, "level": "gfn2-flag",
                   "candidates": cands}, fh, indent=2)
    print(f"\n{len(cands)} candidate(s) within {window} kcal/mol -> {out_json}")
    for c in cands:
        print(f"  {c['node']}: base mult {c['base_mult']} <-> alt {c['alt_mult']}  "
              f"dE_vert={c['dE_vert_kcal']} kcal  ({c['note']})")
    return cands


def main():
    ap = argparse.ArgumentParser(description="M1: spin-state enumeration + close-lying detection")
    ap.add_argument("geoms", nargs="*", help="node geometry files or globs (.xyz/.rxyz)")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, required=False, help="native multiplicity of this network")
    ap.add_argument("--window", type=float, default=15.0,
                    help="GFN2 flag window (kcal/mol); keep WIDE — DFT decides downstream (ESC-009)")
    ap.add_argument("--method", default="gfn2")
    ap.add_argument("--out", default="mecp_candidates.json")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return _selftest()
    if args.mult is None:
        ap.error("--mult is required (native multiplicity of the network)")
    paths = []
    for g in args.geoms:
        paths.extend(sorted(glob.glob(g)) if any(c in g for c in "*?[") else [g])
    if not paths:
        ap.error("no geometry files found")
    run(paths, args.charge, args.mult, args.window, args.method, args.out)


def _selftest():
    """Exercise flag logic with synthetic energies (no xtb needed)."""
    # fake: a quintet node whose triplet lies 3 kcal below (crossing), and one whose
    # triplet lies 40 kcal above (no crossing).
    def fake_energy(syms, coords, charge, mult, method):
        table = {5: -100.000000, 3: -100.000000 + 3.0 / HARTREE_KCAL}   # triplet 3 kcal BELOW quintet
        return table[mult]
    nodes = [
        {"node": "react_quintet", "xyz": "a.xyz", "charge": 0, "base_mult": 5,
         "base_E": -100.0, "alts": {3: -100.0 + 3.0 / HARTREE_KCAL, 7: -100.0 + 80.0 / HARTREE_KCAL}},
        {"node": "far_node", "xyz": "b.xyz", "charge": 0, "base_mult": 5,
         "base_E": -50.0, "alts": {3: -50.0 + 40.0 / HARTREE_KCAL}},
    ]
    cands = flag_candidates(nodes, window_kcal=15.0)
    assert len(cands) == 1 and cands[0]["node"] == "react_quintet", cands
    assert cands[0]["alt_mult"] == 3 and abs(cands[0]["dE_vert_kcal"] - 3.0) < 0.1, cands
    assert alt_multiplicities(5) == [3, 7] and alt_multiplicities(1) == [3], alt_multiplicities(1)
    print("M1 selftest PASS:", cands[0])
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
