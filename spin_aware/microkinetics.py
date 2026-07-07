#!/usr/bin/env python3
"""
microkinetics.py — Stage M4 post-analysis: kinetic role of the spin crossing.

Takes the M4 rate-annotated graph (`annotate_rates.py` output) and answers the question the
whole project is built to answer: *does the spin crossing govern the rate, or does it merely
open access to a lower chemical TS?*

It does this with an energy-span (Kozuch–Shaik) analysis over the two-state path, where the
NA-TST surface-hopping coefficient κ enters as an effective free-energy penalty on the MECP:

    ΔG‡_cross,eff = E_MECP − R T ln κ          (so k_cross = (kT/h) e^{−ΔG‡_cross,eff/RT})

The turnover-determining transition state (TDTS) is the highest effective barrier along the
path; the rate-determining question is simply whether that is the MECP or a chemical TS. Also
emits Eyring rates per step and a mikimo/MicroKatc-style step list for the production engine.

Energies are taken from each node's `rel_kcal` (DFT, vs the reference set in M3 — ESC-009).
"""
import os
import sys
import json
import math
import argparse

KCAL = 4184.0 / 6.02214076e23
kB = 1.380649e-23
h = 6.62607015e-34
R_kcal = 8.314462618 / 4184.0          # kcal/(mol K)


def eyring_k(dG_kcal, T):
    return (kB * T / h) * math.exp(-dG_kcal / (R_kcal * T))


def effective_cross_barrier(e_mecp_kcal, kappa, T):
    """Fold κ into an apparent ΔG‡: E_MECP − RT ln κ (κ<1 raises the effective barrier)."""
    return e_mecp_kcal - R_kcal * T * math.log(max(kappa, 1e-300))


def rate_at_T(node, T):
    """Pick the MECP node's k_of_T entry closest to T (a T-sweep makes index 0 ≠ T)."""
    rows = node.get("rate", {}).get("k_of_T", [])
    if not rows:
        return None
    return min(rows, key=lambda r: abs(r["T"] - T))


def analyse(graph, T=298.15):
    """Route-aware two-state energy-span analysis.

    Compares the two competing routes from the reactant ground-state surface to product:
      * same-spin route : the chemical TS on the reactant's own spin surface;
      * spin-cross route: hop at the MECP onto the other surface, over that surface's TS;
        effective barrier = max(MECP_eff, other-surface TS).
    The operative route is the lower-barrier one; the TDTS is its highest point. The spin
    crossing is 'rate-determining' only if MECP_eff is the highest point on the operative route.
    """
    nodes = {n["id"]: n for n in graph["nodes"]}
    reac = next((n for n in graph["nodes"] if n.get("role") == "reactant"), None)
    reac_mult = reac.get("mult") if reac else None

    # all TS-like points for the ranking table
    points = []
    for n in graph["nodes"]:
        rel = n.get("rel_kcal")
        if rel is None:
            continue
        if n.get("role") == "ts":
            points.append({"id": n["id"], "type": "chem_TS", "mult": n.get("mult"),
                           "barrier_kcal": rel})
        elif n.get("role") == "mecp" and n.get("rate", {}).get("status") == "ok":
            row = rate_at_T(n, T)
            kap = row["kappa"]
            points.append({"id": n["id"], "type": "MECP(eff)", "mult": n.get("mult"),
                           "barrier_kcal": round(effective_cross_barrier(rel, kap, T), 2),
                           "kappa": kap, "E_mecp_kcal": rel, "k_cross_s^-1": row["k_s^-1"]})
    if not points:
        return {"error": "no energy-bearing TS-like points (need DFT rel_kcal on TS/MECP nodes)"}

    # --- build the competing routes (needs reactant surface + MECP crossing edges to TSs) ---
    routes = []
    same_spin = [p for p in points if p["type"] == "chem_TS" and p["mult"] == reac_mult]
    for p in same_spin:
        routes.append({"route": "same-spin", "barrier_kcal": p["barrier_kcal"],
                       "rds": p["id"], "rds_type": "chem_TS"})
    for n in graph["nodes"]:
        if n.get("role") != "mecp" or n.get("rate", {}).get("status") != "ok":
            continue
        mecp_eff = effective_cross_barrier(n["rel_kcal"], rate_at_T(n, T)["kappa"], T)
        # other-surface TS bridged by this MECP's crossing edges
        for e in graph["edges"]:
            if e["kind"] == "nonadiabatic_crossing" and e["source"] == n["id"]:
                tgt = nodes.get(e["target"], {})
                if tgt.get("role") == "ts" and tgt.get("mult") != reac_mult and tgt.get("rel_kcal") is not None:
                    bar = max(mecp_eff, tgt["rel_kcal"])
                    rds, rtype = ((n["id"], "MECP(eff)") if mecp_eff >= tgt["rel_kcal"]
                                  else (tgt["id"], "chem_TS"))
                    routes.append({"route": f"spin-cross via {n['id']}",
                                   "barrier_kcal": round(bar, 2), "rds": rds, "rds_type": rtype,
                                   "mecp_eff_kcal": round(mecp_eff, 2),
                                   "cross_TS": tgt["id"], "cross_TS_kcal": tgt["rel_kcal"]})

    res = {"T": T, "reactant_mult": reac_mult,
           "points": sorted(points, key=lambda p: p["barrier_kcal"], reverse=True),
           "routes": sorted(routes, key=lambda r: r["barrier_kcal"])}
    if routes:
        op = res["routes"][0]                          # operative = lowest-barrier route
        res["operative_route"] = op
        res["energy_span_kcal"] = op["barrier_kcal"]
        res["TOF_s^-1"] = eyring_k(op["barrier_kcal"], T)
        res["TDTS"] = {"id": op["rds"], "type": op["rds_type"]}
        spin_rds = op["rds_type"] == "MECP(eff)"
        enables = (op["route"].startswith("spin-cross") and len(res["routes"]) > 1 and
                   res["routes"][1]["barrier_kcal"] > op["barrier_kcal"])
        if spin_rds:
            res["verdict"] = (f"spin-crossing IS rate-determining (MECP effective barrier "
                              f"{op['barrier_kcal']} kcal is the highest point on the operative route)")
        else:
            res["verdict"] = (f"spin-crossing is NOT rate-determining — the operative route is the "
                              f"spin-cross route (barrier {op['barrier_kcal']} kcal, RDS {op['rds']}); "
                              f"the MECP is below it" +
                              (f", and the crossing opens this route below the same-spin route "
                               f"({res['routes'][1]['barrier_kcal']} kcal)" if enables else ""))
    else:
        tdts = max(points, key=lambda p: p["barrier_kcal"])
        res.update({"TDTS": tdts, "energy_span_kcal": tdts["barrier_kcal"],
                    "TOF_s^-1": eyring_k(tdts["barrier_kcal"], T),
                    "verdict": "routes undetermined (no MECP→TS crossing edges); showing global ranking"})
    return res


def mikimo_export(graph, T=298.15):
    """A minimal step list (ΔG‡ per step vs reactant) for mikimo / MicroKatc / Cantera-style input.
    Production engines want a full mechanism; this is the spin-aware skeleton to extend."""
    lines = ["# step  type  from->to  dG_act_kcal  k_s^-1  (T=%.2f K)" % T]
    for n in graph["nodes"]:
        rel = n.get("rel_kcal")
        if n.get("role") == "ts" and rel is not None:
            lines.append(f"chemTS  {n['id']}  dG‡={rel:.2f}  k={eyring_k(rel,T):.3e}")
        elif n.get("role") == "mecp" and n.get("rate", {}).get("status") == "ok" and rel is not None:
            row = rate_at_T(n, T)
            eff = effective_cross_barrier(rel, row["kappa"], T)
            lines.append(f"MECP    {n['id']}  dG‡eff={eff:.2f} (E_MECP={rel:.2f}, κ={row['kappa']:.3f})  "
                         f"k={row['k_s^-1']:.3e}")
    return "\n".join(lines) + "\n"


def report_md(a):
    if "error" in a:
        return "# M4 microkinetics\n\n" + a["error"] + "\n"
    out = ["# M4 microkinetic post-analysis (route-aware energy-span)", "",
           f"T = {a['T']:.2f} K. Energies vs reactant GS (DFT, ESC-009). "
           f"Reactant surface: mult {a.get('reactant_mult')}.", "",
           "## TS-like points (ranked)", "",
           "| rank | point | type | mult | eff. barrier (kcal/mol) |",
           "|---|---|---|---|---|"]
    for i, p in enumerate(a["points"], 1):
        out.append(f"| {i} | {p['id']} | {p['type']} | {p.get('mult')} | {p['barrier_kcal']} |")
    if a.get("routes"):
        out += ["", "## Competing routes (reactant → product)", "",
                "| route | barrier (kcal/mol) | rate-determining point |", "|---|---|---|"]
        for r in a["routes"]:
            mark = "  ← operative" if r is a["routes"][0] else ""
            out.append(f"| {r['route']} | {r['barrier_kcal']} | {r['rds']} ({r['rds_type']}){mark} |")
    out += ["", f"**TDTS:** {a['TDTS']['id']} ({a['TDTS']['type']}), "
            f"energy span **{a['energy_span_kcal']:.2f} kcal/mol** → "
            f"TOF ≈ **{a['TOF_s^-1']:.3e} s⁻¹**.", "",
            f"**Verdict:** {a['verdict']}." if a.get("verdict") else ""]
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description="M4 microkinetics: spin-crossing kinetic role (energy span)")
    ap.add_argument("--graph", help="rate-annotated graph json (annotate_rates output)")
    ap.add_argument("--T", type=float, default=298.15)
    ap.add_argument("--out", default="microkinetics")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return _selftest()
    if not a.graph:
        ap.error("--graph is required (or use --selftest)")
    graph = json.load(open(a.graph))
    res = analyse(graph, a.T)
    json.dump(res, open(a.out + ".json", "w"), indent=2)
    open(a.out + "_report.md", "w").write(report_md(res))
    open(a.out + "_mikimo.txt", "w").write(mikimo_export(graph, a.T))
    print(report_md(res))


def _selftest():
    # Tp-Fe-like two-state graph with MECP crossing edges to the β-H TS on each surface.
    graph = {"nodes": [
        {"id": "q_react", "role": "reactant", "mult": 5, "rel_kcal": 0.0},
        {"id": "t_betaH_TS", "role": "ts", "mult": 3, "rel_kcal": 26.93},
        {"id": "q_betaH_TS", "role": "ts", "mult": 5, "rel_kcal": 40.41},
        {"id": "MECP_5_3_0", "role": "mecp", "mult": "5/3", "rel_kcal": 19.56,
         "rate": {"status": "ok", "k_of_T": [{"T": 298.15, "k_s^-1": 1.5e-2, "kappa": 0.573}]}}],
        "edges": [
        {"source": "MECP_5_3_0", "target": "q_betaH_TS", "kind": "nonadiabatic_crossing", "surface_mult": 5},
        {"source": "MECP_5_3_0", "target": "t_betaH_TS", "kind": "nonadiabatic_crossing", "surface_mult": 3}]}
    a = analyse(graph, 298.15)
    # Two routes: same-spin (quintet TS 40.41) vs spin-cross (max(MECP_eff~19.9, triplet TS 26.93)=26.93).
    # Operative = spin-cross route (26.93 < 40.41); RDS = the triplet chemical TS, NOT the MECP.
    assert a["operative_route"]["route"].startswith("spin-cross"), a["operative_route"]
    assert abs(a["energy_span_kcal"] - 26.93) < 0.01, a["energy_span_kcal"]
    assert a["TDTS"]["id"] == "t_betaH_TS" and a["TDTS"]["type"] == "chem_TS", a["TDTS"]
    assert "NOT rate-determining" in a["verdict"], a["verdict"]
    assert "opens this route below the same-spin route" in a["verdict"], a["verdict"]
    assert a["TOF_s^-1"] > 0
    eff = effective_cross_barrier(19.56, 0.573, 298.15)
    assert eff < 26.93, eff
    mk = mikimo_export(graph, 298.15)
    assert "MECP" in mk and "κ=0.573" in mk
    print("microkinetics selftest PASS")
    print(f"  routes: same-spin 40.41 vs spin-cross 26.93 -> operative = spin-cross")
    print(f"  TDTS={a['TDTS']['id']} (triplet chemical TS), span={a['energy_span_kcal']:.2f} kcal, "
          f"TOF={a['TOF_s^-1']:.2e} s^-1; MECP_eff={eff:.2f} < 26.93 -> crossing not rate-limiting")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
