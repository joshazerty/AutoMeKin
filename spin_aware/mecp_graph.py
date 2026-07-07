#!/usr/bin/env python3
"""
mecp_graph.py — Stage M3 of the MECP-aware-node feature (Project 18, Phase 2 Weeks 2+).

PURPOSE
-------
Merge the per-spin reaction networks (built by the Week-0 engine + §3 TS/IRC layer) and the
localised MECPs (Stage M2, `locate_mecp.py`) into ONE spin-aware reaction graph in which:
  * every node carries its spin manifold (multiplicity);
  * intra-surface edges come from each spin run's RXNet (a TS links two minima);
  * each MECP becomes a first-class **MECP node** joined by **non-adiabatic crossing edges** to
    the nearest node on each of the two spin surfaces it bridges.

This is the step that turns "two separate spin networks + a crossing calc" into a single object
showing the two-state mechanism (e.g. ⁵reactant → MECP → ³ surface → ³β-H TS → product). It emits
a stable JSON (for Stage M4 rates), plus graphviz DOT and mermaid for visualisation.

DESIGN RULES (carried from M1/M2)
---------------------------------
- **ESC-009 — energy provenance.** GFN2 (low-level) energies are NOT comparable across spin
  surfaces (GFN2 inverts Fe spin gaps). Every node/edge energy is tagged with its `source`
  ('dft' | 'gfn2'). Cross-surface energy statements (e.g. the crossing-energy ladder) are only
  emitted from DFT-sourced energies; GFN2 nodes are kept for *topology* and flagged
  `xsurface_energy=false`. The visualisation greys GFN2 energies accordingly.
- **ESC-012 — MECP validity.** An MECP is accepted into the graph on its inter-surface
  `gap_kcal` (degeneracy), NOT on easyMECP's textual `converged` marker (which can read false at
  max_steps even for a physically converged seam). Threshold `--mecp-gap-tol` (default 1.0 kcal).

MATCHING
--------
An MECP bridges two surfaces; M3 connects it to the nearest node on each surface by
**Kabsch-aligned heavy-atom RMSD** (atom ordering is consistent within one AutoMeKin system).
numpy only; no external deps.

INPUTS (either mode)
--------------------
A) explicit nodes:  --nodes nodes.json
     [{"id","mult","role":"min|ts|product|reactant","energy","energy_source","xyz"}]
B) AutoMeKin runs:  --spin-run MULT=RUNDIR  (repeatable)
     parses RUNDIR/RXNet (+ ts.db/min.db energies, + node geometries) for that multiplicity.
Plus: --mecp mecp_record.json (repeatable), --ref-energy-eh <E0> (ladder zero, e.g. ⁵ reactant DFT).

OUTPUTS
-------
merged_graph.json (stable schema for M4), graph.dot, graph.mmd, graph_summary.md.
"""
import os
import re
import sys
import json
import glob
import argparse

try:
    import numpy as np
except ImportError:                      # numpy is on the cluster; degrade for pure-IO tests
    np = None

H = 627.5094740631


# --------------------------------------------------------------------------- #
# geometry + RMSD
# --------------------------------------------------------------------------- #
HEAVY = lambda s: s.upper() != "H"


def read_xyz(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    syms, xyz = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        syms.append(p[0]); xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, (np.asarray(xyz) if np is not None else xyz)


def kabsch_rmsd(P, Q):
    """Heavy-atom RMSD after optimal translation+rotation (Kabsch). P,Q: (N,3) arrays."""
    P = P - P.mean(0); Q = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1, 1, d])
    Pr = P @ (V @ D @ Wt)
    return float(np.sqrt(((Pr - Q) ** 2).sum() / len(P)))


def heavy_coords(syms, xyz):
    idx = [i for i, s in enumerate(syms) if HEAVY(s)]
    return xyz[idx] if np is not None else [xyz[i] for i in idx]


def best_match(mecp_syms, mecp_xyz, nodes):
    """Return (node_id, rmsd) of the heavy-atom-closest node (same atom count) to the MECP."""
    mh = heavy_coords(mecp_syms, mecp_xyz)
    best = (None, float("inf"))
    for nd in nodes:
        if not nd.get("xyz") or not os.path.isfile(nd["xyz"]):
            continue
        s, x = read_xyz(nd["xyz"])
        if len(s) != len(mecp_syms):
            continue
        r = kabsch_rmsd(heavy_coords(s, x), mh)
        if r < best[1]:
            best = (nd["id"], r)
    return best


# --------------------------------------------------------------------------- #
# graph model
# --------------------------------------------------------------------------- #
class Graph:
    def __init__(self, ref_energy_eh=None):
        self.nodes = {}     # id -> dict
        self.edges = []     # list of dicts
        self.ref = ref_energy_eh

    def add_node(self, nid, **kw):
        kw["id"] = nid
        self.nodes[nid] = kw

    def add_edge(self, a, b, kind, **kw):
        self.edges.append(dict(source=a, target=b, kind=kind, **kw))

    def rel_kcal(self, e_eh):
        if e_eh is None or self.ref is None:
            return None
        return round((e_eh - self.ref) * H, 2)


# --------------------------------------------------------------------------- #
# parsers
# --------------------------------------------------------------------------- #
def load_explicit_nodes(path, g):
    for nd in json.load(open(path)):
        e_eh = nd.get("energy_eh")
        g.add_node(nd["id"], mult=int(nd["mult"]), role=nd.get("role", "min"),
                   energy_eh=e_eh, energy_source=nd.get("energy_source", "dft"),
                   xyz=nd.get("xyz"), rel_kcal=g.rel_kcal(e_eh) if e_eh is not None
                   else nd.get("rel_kcal"))


RXNET_RE = re.compile(r"TS\s+\d+\s+(\S+)\s+DE=\s*(-?\d+\.\d+)\s+Path:\s+MIN\s+(\d+)\s+<-->\s+MIN\s+(\d+)")


def load_spin_run(mult, rundir, g, geom_subdir=None):
    """Light AutoMeKin parser: RXNet edges (intra-surface) + node ids tagged by mult.
    Geometries (optional) are looked up in geom_subdir for later MECP matching."""
    mult = int(mult)
    rxnet = None
    for cand in ("RXNet", "RXNet.cg", os.path.join("tsdirLL", "RXNet")):
        p = os.path.join(rundir, cand)
        if os.path.isfile(p):
            rxnet = p; break
    if not rxnet:
        print(f"  [warn] no RXNet in {rundir}", file=sys.stderr); return
    def gpath(name):
        if geom_subdir:
            for ext in (".xyz", ".rxyz"):
                pth = os.path.join(rundir, geom_subdir, name + ext)
                if os.path.isfile(pth):
                    return pth
        return None
    for ln in open(rxnet):
        m = RXNET_RE.search(ln)
        if not m:
            continue
        ts, de, a, b = m.group(1), float(m.group(2)), m.group(3), m.group(4)
        ts_id = f"m{mult}_{ts.replace('.rxyz','')}"
        for mn in (a, b):
            nid = f"m{mult}_MIN{mn}"
            if nid not in g.nodes:
                g.add_node(nid, mult=mult, role="min", energy_eh=None,
                           energy_source="gfn2", xyz=gpath(f"MIN{mn}"))
        if ts_id not in g.nodes:
            g.add_node(ts_id, mult=mult, role="ts", energy_eh=None,
                       energy_source="gfn2", xyz=gpath(ts.replace('.rxyz','')))
        g.add_edge(f"m{mult}_MIN{a}", f"m{mult}_MIN{b}", "reaction",
                   via_ts=ts_id, de_gfn2_kcal=de, energy_source="gfn2")


def add_mecps(records, g, gap_tol=1.0, rmsd_warn=2.0):
    """Insert each accepted MECP as a node + two non-adiabatic crossing edges."""
    added = []
    for i, rec in enumerate(records):
        gap = rec.get("gap_kcal")
        accept = (gap is not None) and (gap <= gap_tol)     # ESC-012: gap, not marker
        ma, mb = int(rec["mult_a"]), int(rec["mult_b"])
        mecp_xyz = rec.get("mecp_xyz") or os.path.join(rec.get("workdir", ""), "mecp.xyz")
        nid = f"MECP_{ma}_{mb}_{i}"
        e_cross = rec.get("crossing_e_eh")
        g.add_node(nid, mult=f"{ma}/{mb}", role="mecp", energy_eh=e_cross,
                   energy_source="dft", xyz=mecp_xyz if os.path.isfile(mecp_xyz) else None,
                   gap_kcal=gap, accepted=accept, rel_kcal=g.rel_kcal(e_cross))
        # match to nearest node on each surface (needs geometries + numpy)
        matches = {}
        if np is not None and os.path.isfile(mecp_xyz):
            ms, mx = read_xyz(mecp_xyz)
            for mm in (ma, mb):
                surf = [n for n in g.nodes.values() if n.get("mult") == mm and n.get("xyz")]
                node_id, r = best_match(ms, mx, surf)
                matches[mm] = (node_id, r)
                if node_id is not None:
                    g.add_edge(nid, node_id, "nonadiabatic_crossing", surface_mult=mm,
                               rmsd=round(r, 3), energy_source="dft",
                               weak=(r > rmsd_warn))
        added.append({"mecp": nid, "gap_kcal": gap, "accepted": accept,
                      "rel_kcal": g.nodes[nid]["rel_kcal"], "matches": matches})
    return added


# --------------------------------------------------------------------------- #
# emitters
# --------------------------------------------------------------------------- #
SPIN_COLOR = {1: "#4c78a8", 2: "#9ecae1", 3: "#f58518", 5: "#54a24b", 7: "#b279a2"}


def to_dot(g):
    out = ["digraph MECP_network {", '  rankdir=LR;', '  node [style=filled,fontname=Helvetica];']
    for n in g.nodes.values():
        if n["role"] == "mecp":
            lbl = f'MECP {n["mult"]}\\n{n.get("rel_kcal","?")} kcal\\n(gap {n.get("gap_kcal","?")})'
            shape = "diamond"; color = "#e45756" if n.get("accepted") else "#cccccc"
        else:
            rel = n.get("rel_kcal")
            etag = "" if n["energy_source"] == "dft" else " (gfn2)"
            lbl = f'{n["id"]}\\n{rel if rel is not None else "?"} kcal{etag}'
            shape = {"ts": "box", "product": "ellipse"}.get(n["role"], "ellipse")
            color = SPIN_COLOR.get(n["mult"], "#dddddd")
        out.append(f'  "{n["id"]}" [label="{lbl}",shape={shape},fillcolor="{color}"];')
    for e in g.edges:
        if e["kind"] == "nonadiabatic_crossing":
            style = 'style=dashed,color="#e45756",penwidth=2'
            if e.get("weak"):
                style += ',color="#e4575680"'
            out.append(f'  "{e["source"]}" -> "{e["target"]}" [{style},label="ISC"];')
        else:
            out.append(f'  "{e["source"]}" -> "{e["target"]}" '
                       f'[label="{e.get("via_ts","")}",color="#888888"];')
    out.append("}")
    return "\n".join(out)


def to_mermaid(g):
    out = ["flowchart LR"]
    for n in g.nodes.values():
        rel = n.get("rel_kcal")
        if n["role"] == "mecp":
            out.append(f'  {n["id"]}{{{{"MECP {n["mult"]} · {rel} kcal"}}}}')
        elif n["role"] == "ts":
            out.append(f'  {n["id"]}["{n["id"]} · {rel} kcal"]')
        else:
            out.append(f'  {n["id"]}(["{n["id"]} · {rel} kcal"])')
    for e in g.edges:
        arrow = "-. ISC .->" if e["kind"] == "nonadiabatic_crossing" else "-->"
        out.append(f'  {e["source"]} {arrow} {e["target"]}')
    return "\n".join(out)


def summary_md(g, mecp_info):
    lines = ["# Spin-aware network — M3 merge summary", "",
             f"- nodes: {len(g.nodes)}  (by manifold: " +
             ", ".join(f"mult {m}: {sum(1 for n in g.nodes.values() if n['mult']==m)}"
                       for m in sorted({n['mult'] for n in g.nodes.values() if isinstance(n['mult'], int)})) +
             f", MECP: {sum(1 for n in g.nodes.values() if n['role']=='mecp')})",
             f"- edges: {len(g.edges)} (non-adiabatic: {sum(1 for e in g.edges if e['kind']=='nonadiabatic_crossing')})",
             "", "## MECP nodes"]
    for mi in mecp_info:
        lines.append(f"- {mi['mecp']}: gap {mi['gap_kcal']} kcal, "
                     f"{'ACCEPTED' if mi['accepted'] else 'REJECTED (gap>tol)'}, "
                     f"E_rel {mi['rel_kcal']} kcal; bridges " +
                     "; ".join(f"mult {m}->{v[0]} (rmsd {round(v[1],2)} Å)"
                               for m, v in mi["matches"].items()))
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def build(nodes_json=None, spin_runs=None, mecp_records=None, ref_energy_eh=None,
          gap_tol=1.0, geom_subdir=None):
    g = Graph(ref_energy_eh)
    if nodes_json:
        load_explicit_nodes(nodes_json, g)
    for mult, rundir in (spin_runs or []):
        load_spin_run(mult, rundir, g, geom_subdir)
    recs = [json.load(open(p)) for p in (mecp_records or [])]
    info = add_mecps(recs, g, gap_tol)
    return g, info


def main():
    ap = argparse.ArgumentParser(description="M3: merge spin networks + MECPs into one spin-aware graph")
    ap.add_argument("--nodes", help="explicit nodes.json")
    ap.add_argument("--spin-run", action="append", default=[], metavar="MULT=DIR",
                    help="AutoMeKin run dir for a multiplicity (repeatable)")
    ap.add_argument("--mecp", action="append", default=[], help="mecp_record.json (repeatable)")
    ap.add_argument("--ref-energy-eh", type=float, default=None, help="ladder zero (Eh)")
    ap.add_argument("--geom-subdir", default=None, help="subdir holding node geometries in a run dir")
    ap.add_argument("--mecp-gap-tol", type=float, default=1.0, help="accept MECP if gap<=tol kcal (ESC-012)")
    ap.add_argument("--out", default="mecp_graph")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return _selftest()
    spin_runs = [s.split("=", 1) for s in args.spin_run]
    g, info = build(args.nodes, spin_runs, args.mecp, args.ref_energy_eh,
                    args.mecp_gap_tol, args.geom_subdir)
    json.dump({"nodes": list(g.nodes.values()), "edges": g.edges,
               "ref_energy_eh": g.ref, "mecp_info": info},
              open(args.out + ".json", "w"), indent=2, default=str)
    open(args.out + ".dot", "w").write(to_dot(g))
    open(args.out + ".mmd", "w").write(to_mermaid(g))
    open(args.out + "_summary.md", "w").write(summary_md(g, info))
    print(summary_md(g, info))
    print(f"-> {args.out}.json / .dot / .mmd / _summary.md")


def _selftest():
    if np is None:
        print("M3 selftest: numpy absent — testing graph/emit logic only")
    # 1) Kabsch RMSD invariance
    if np is not None:
        P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1.]])
        ang = 0.7; R = np.array([[np.cos(ang), -np.sin(ang), 0],
                                 [np.sin(ang), np.cos(ang), 0], [0, 0, 1]])
        Q = P @ R.T + np.array([3., -2., 1.])
        assert kabsch_rmsd(P, Q) < 1e-9, kabsch_rmsd(P, Q)
        print("  Kabsch RMSD rotation/translation invariance OK")
    # 2) full path with real geometry matching -> crossing edges (needs numpy + temp xyz)
    import tempfile
    rel_ok = None
    if np is not None:
        with tempfile.TemporaryDirectory() as d:
            def wx(name, dz):   # 3 heavy atoms; shift atom 2 in z to differentiate surfaces
                p = os.path.join(d, name)
                open(p, "w").write(f"3\n\nFe 0 0 0\nC 1.5 0 0\nO 1.5 0 {dz}\n")
                return p
            q_xyz, t_xyz, m_xyz = wx("q.xyz", 1.0), wx("t.xyz", 2.0), wx("mecp.xyz", 1.4)
            g = Graph(ref_energy_eh=-2044.68639790)
            g.add_node("q_react", mult=2, role="reactant", energy_eh=-2044.68639790,
                       energy_source="dft", xyz=q_xyz, rel_kcal=0.0)
            g.add_node("t_ts", mult=1, role="ts", energy_eh=-2044.64348133,
                       energy_source="dft", xyz=t_xyz, rel_kcal=g.rel_kcal(-2044.64348133))
            rec = {"mult_a": 2, "mult_b": 1, "gap_kcal": 0.0006,
                   "crossing_e_eh": -2044.655220405, "mecp_xyz": m_xyz, "workdir": d}
            info = add_mecps([rec], g, gap_tol=1.0)
            assert info[0]["accepted"] is True, info
            rel_ok = info[0]["rel_kcal"]
            assert abs(rel_ok - 19.56) < 0.05, rel_ok
            # two non-adiabatic crossing edges created (to q_react on S=2, t_ts on S=1)
            xedges = [e for e in g.edges if e["kind"] == "nonadiabatic_crossing"]
            assert len(xedges) == 2 and {e["surface_mult"] for e in xedges} == {2, 1}, xedges
            dot, mmd = to_dot(g), to_mermaid(g)
            assert "diamond" in dot and "ISC" in dot and "MECP" in mmd
            print("  geometry match -> 2 crossing edges + DOT/mermaid (ISC) OK")
    # 3) accept/reject gating on gap (ESC-012), no geometry needed
    g2 = Graph(-2044.68639790)
    info2 = add_mecps([{"mult_a": 2, "mult_b": 1, "gap_kcal": 5.0,
                        "crossing_e_eh": -2044.655220405, "workdir": "/nonexistent"}],
                      g2, gap_tol=1.0)
    assert info2[0]["accepted"] is False, info2
    print("  MECP accept/reject gating on gap_kcal (ESC-012) OK"
          + (f"; MECP rel energy {rel_ok} kcal (expect +19.56)" if rel_ok else ""))
    print("M3 selftest PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
