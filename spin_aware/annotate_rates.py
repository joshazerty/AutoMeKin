#!/usr/bin/env python3
"""
annotate_rates.py — Stage M4: attach NA-TST spin-crossing rates to the M3 graph.

Reads the M3 spin-aware graph (`*_graph.json`) and a per-MECP physical-inputs file
(`mecp_rate_inputs.json`), computes the Landau–Zener / NA-TST rate for each accepted MECP
(via na_tst.py), and writes the rate onto that MECP node and BOTH its
`nonadiabatic_crossing` edges. Output is an M4-annotated graph (same schema + a `rate` block),
plus a rates summary and an optional SOC sensitivity sweep.

rate inputs schema (mecp_rate_inputs.json):
{
  "MECP_5_3_0": {
     "dE_kcal": 19.56,            # MECP vs reactant (elec; set add_zpe + freqs to ZPE-correct)
     "H_SO_cm": 200,              # spin-orbit coupling at the MECP (ORCA SOMF; see INTEGRATION_M4)
     "dF_eh_bohr": 0.05,          # |grad_A - grad_B| at the MECP (from the two surface gradients)
     "mu_amu": 10.0,              # effective hop mass (sensitivity-swept; default 1.0)
     "freqs_react": [...],        # cm^-1, 3N-6
     "freqs_mecp": [...],         # cm^-1, 3N-7 (seam)
     "add_zpe": false
  }
}

Only MECP nodes with `accepted=true` (ESC-012 gap test, set by M3) AND an inputs entry are rated.
"""
import os
import sys
import json
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import na_tst


def annotate(graph, inputs, T=298.15, tsweep=None, passage="double", soc_sweep=None):
    temps = tsweep or [T]
    by_id = {n["id"]: n for n in graph["nodes"]}
    rated = []
    for nid, n in by_id.items():
        if n.get("role") != "mecp":
            continue
        if not n.get("accepted"):
            n["rate"] = {"status": "skipped", "reason": "MECP not accepted (gap>tol, ESC-012)"}
            continue
        d = inputs.get(nid)
        if not d:
            n["rate"] = {"status": "skipped", "reason": "no rate inputs provided"}
            continue
        kw = dict(dE_kcal=d["dE_kcal"], H_SO_cm=d["H_SO_cm"], dF_eh_bohr=d["dF_eh_bohr"],
                  mu_amu=d.get("mu_amu", 1.0), freqs_react=d.get("freqs_react"),
                  freqs_mecp=d.get("freqs_mecp"), passage=passage,
                  add_zpe=d.get("add_zpe", False))
        sweep = na_tst.k_sweep(temps, **kw)
        rate = {"status": "ok", "model": "NA-TST / Landau-Zener (%s passage)" % passage,
                "inputs": {k: d.get(k) for k in ("dE_kcal", "H_SO_cm", "dF_eh_bohr", "mu_amu")},
                "k_of_T": [{"T": r["T"], "k_s^-1": r["k_s^-1"], "kappa": r["kappa"]} for r in sweep]}
        if soc_sweep:
            rate["soc_sweep_T"] = temps[0]
            rate["soc_sweep"] = []
            for hso in soc_sweep:
                kw2 = dict(kw); kw2["H_SO_cm"] = hso
                rr = na_tst.k_natst(temps[0], **kw2)
                rate["soc_sweep"].append({"H_SO_cm": hso, "kappa": rr["kappa"], "k_s^-1": rr["k_s^-1"]})
        n["rate"] = rate
        # stamp both crossing edges
        for e in graph["edges"]:
            if e["kind"] == "nonadiabatic_crossing" and e["source"] == nid:
                e["rate_k_s^-1"] = sweep[0]["k_s^-1"]
                e["rate_kappa"] = sweep[0]["kappa"]
                e["rate_T"] = sweep[0]["T"]
        rated.append((nid, sweep))
    return rated


def summary_md(rated, passage):
    out = ["# M4 NA-TST rate annotation summary", "",
           f"Model: Landau–Zener / NA-TST ({passage} passage). k in s⁻¹.", ""]
    for nid, sweep in rated:
        out.append(f"## {nid}")
        for r in sweep:
            out.append(f"- T={r['T']:.2f} K: κ={r['kappa']:.3f}, "
                       f"k={r['k_s^-1']:.3e} s⁻¹, prefactor={r['prefactor_s^-1']:.3e} s⁻¹")
        out.append("")
    if not rated:
        out.append("_No MECP nodes rated (none accepted with inputs)._")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description="M4: annotate M3 graph with NA-TST spin-crossing rates")
    ap.add_argument("--graph", help="M3 *_graph.json")
    ap.add_argument("--inputs", help="mecp_rate_inputs.json")
    ap.add_argument("--T", type=float, default=298.15)
    ap.add_argument("--tsweep", help="comma T list")
    ap.add_argument("--passage", default="double", choices=["single", "double"])
    ap.add_argument("--soc-sweep", help="comma H_SO list (cm^-1)")
    ap.add_argument("--out", default="graph_rated")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return _selftest()
    if not (a.graph and a.inputs):
        ap.error("--graph and --inputs are required (or use --selftest)")
    graph = json.load(open(a.graph))
    inputs = json.load(open(a.inputs))
    temps = [float(x) for x in a.tsweep.split(",")] if a.tsweep else [a.T]
    soc = [float(x) for x in a.soc_sweep.split(",")] if a.soc_sweep else None
    rated = annotate(graph, inputs, a.T, temps, a.passage, soc)
    json.dump(graph, open(a.out + ".json", "w"), indent=2, default=str)
    open(a.out + "_summary.md", "w").write(summary_md(rated, a.passage))
    print(summary_md(rated, a.passage))
    print(f"-> {a.out}.json / _summary.md  ({len(rated)} MECP node(s) rated)")


def _selftest():
    graph = {"nodes": [
        {"id": "MECP_5_3_0", "role": "mecp", "accepted": True, "mult": "5/3", "rel_kcal": 19.56},
        {"id": "q_react", "role": "reactant", "mult": 5, "rel_kcal": 0.0},
        {"id": "MECP_bad", "role": "mecp", "accepted": False, "mult": "5/3"}],
        "edges": [
        {"source": "MECP_5_3_0", "target": "q_react", "kind": "nonadiabatic_crossing", "surface_mult": 5},
        {"source": "MECP_5_3_0", "target": "t_ts", "kind": "nonadiabatic_crossing", "surface_mult": 3}]}
    inputs = {"MECP_5_3_0": {"dE_kcal": 19.56, "H_SO_cm": 200, "dF_eh_bohr": 0.05, "mu_amu": 10.0,
                             "freqs_react": [200, 400, 600, 900, 1200],
                             "freqs_mecp": [210, 410, 620, 950]}}
    rated = annotate(graph, inputs, T=298.15, soc_sweep=[100, 400])
    assert len(rated) == 1, rated
    mecp = [n for n in graph["nodes"] if n["id"] == "MECP_5_3_0"][0]
    assert mecp["rate"]["status"] == "ok" and mecp["rate"]["k_of_T"][0]["k_s^-1"] > 0
    bad = [n for n in graph["nodes"] if n["id"] == "MECP_bad"][0]
    assert bad["rate"]["status"] == "skipped"           # not accepted -> skipped
    xed = [e for e in graph["edges"] if e["kind"] == "nonadiabatic_crossing"]
    assert all("rate_k_s^-1" in e for e in xed)          # both edges stamped
    assert len(mecp["rate"]["soc_sweep"]) == 2
    print("annotate_rates selftest PASS")
    print(f"  MECP k(298)={mecp['rate']['k_of_T'][0]['k_s^-1']:.3e} s^-1, "
          f"both crossing edges stamped; rejected MECP skipped; SOC sweep ok")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
