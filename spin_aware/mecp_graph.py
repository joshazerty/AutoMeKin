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

CHANGES IN THIS COPY (not committed; review copy of spin_aware/mecp_graph.py @ ca2aa835)
---------------------------------------------------------------------------------------
1. Native AutoMeKin RXNet (the `TS #  DE(kcal/mol)  Reaction path information` table written by
   a run's FINAL_* directory) is read directly, including product-channel lines
   (`MIN 7 ----> PR1: CO + C4FeO4`, `PR4: ... <---> PR5: ...`). Product channels become nodes with
   role "product" (id m<mult>_PR<k>, geometry from <geom_subdir>/PR<k>.xyz when present). The
   previous one-line format is still read. If a run dir holds a native MINinfo, minima that occur in no
   TS line are added as isolated min nodes.
2. MECP records that are the same structure are merged into one MECP node: heavy-atom RMSD with atom
   permutation (same elements) < --dedupe-rmsd (default 0.30 A) and |dE| < --dedupe-de (default
   0.5 kcal/mol). Records without a readable geometry are never merged. --no-dedupe restores the old
   behaviour; merged records are listed in mecp_info[*]["merged"].
3. Optional --perm-match: nearest-node matching uses the permutation-aware RMSD instead of the
   fixed-atom-order Kabsch RMSD (default off, so existing graphs are unchanged).
4. Optional --partner-role MULT=ROLE (repeatable): on surface MULT the MECP is bridged to the nearest
   node of that role (e.g. 3=product for a ligand-recombination MECP bridging the singlet minimum and
   the triplet fragment channel).
"""
import os
import re
import sys
import json
import glob
import argparse
import itertools

try:
    import numpy as np
except ImportError:                      # numpy is on the cluster; degrade for pure-IO tests
    np = None
try:
    from scipy.optimize import linear_sum_assignment   # only for the permutation-aware RMSD
except ImportError:
    linear_sum_assignment = None

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


def _rand_rotations(n, seed=0):
    rng = np.random.default_rng(seed)
    out = [np.eye(3)]
    for _ in range(n - 1):
        q = rng.normal(size=4); q /= np.linalg.norm(q); a, b, c, d = q
        out.append(np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                             [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                             [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]]))
    return out


def _kabsch_rot(P, Q):
    U, S, Vt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(U @ Vt))
    return U @ np.diag([1, 1, d]) @ Vt


def rmsd_perm(syms1, X1, syms2, X2, heavy_only=True, nstart=30):
    """Heavy-atom RMSD minimised over rotation AND over the permutation of like atoms
    (alternating Hungarian assignment / Kabsch from `nstart` starting orientations).
    inf if the element lists differ. Needs numpy and scipy."""
    if linear_sum_assignment is None:
        raise RuntimeError("rmsd_perm needs scipy (linear_sum_assignment)")
    s1 = [s for s in syms1]; s2 = [s for s in syms2]
    X1 = np.asarray(X1, float); X2 = np.asarray(X2, float)
    if heavy_only:
        i1 = [i for i, s in enumerate(s1) if HEAVY(s)]; i2 = [i for i, s in enumerate(s2) if HEAVY(s)]
        s1 = [s1[i] for i in i1]; X1 = X1[i1]; s2 = [s2[i] for i in i2]; X2 = X2[i2]
    if sorted(s1) != sorted(s2):
        return float("inf")
    e1 = np.array(s1); e2 = np.array(s2)
    X1 = X1 - X1.mean(0); X2 = X2 - X2.mean(0)
    best = float("inf")
    for R0 in _rand_rotations(nstart):
        P = X1 @ R0; prev = None
        for _ in range(40):
            perm = np.empty(len(e1), int)
            for el in set(s1):
                a = np.where(e1 == el)[0]; b = np.where(e2 == el)[0]
                D = ((P[a][:, None, :] - X2[b][None, :, :]) ** 2).sum(-1)
                r, c = linear_sum_assignment(D); perm[a[r]] = b[c]
            Q = X2[perm]; P = P @ _kabsch_rot(P, Q)
            v = float(np.sqrt(((P - Q) ** 2).sum(1).mean()))
            if prev is not None and abs(prev - v) < 1e-9:
                break
            prev = v
        best = min(best, v)
    return best


def best_match(mecp_syms, mecp_xyz, nodes, perm=False):
    """Return (node_id, rmsd) of the heavy-atom-closest node (same atom count) to the MECP."""
    mh = heavy_coords(mecp_syms, mecp_xyz)
    best = (None, float("inf"))
    for nd in nodes:
        if not nd.get("xyz") or not os.path.isfile(nd["xyz"]):
            continue
        s, x = read_xyz(nd["xyz"])
        if len(s) != len(mecp_syms):
            continue
        if perm:
            r = rmsd_perm(s, x, mecp_syms, mecp_xyz)
        else:
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
# native AutoMeKin RXNet line:  "  <ts#>  <DE>   MIN  a <---> MIN  b"  /  "MIN  a ----> PRk:  CO + C4FeO4"
# /  "PR4:  CO + C4FeO4 <---> PR5:  CO + C4FeO4"
_END = r"(?:MIN\s+(\d+)|PR(\d+):\s*(.*?))"
NATIVE_RE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+" + _END + r"\s*(<--->|---->)\s*" + _END + r"\s*$")


def _native_end(m, i):
    """-> ('MIN', n, None) or ('PR', k, formula) from regex groups starting at index i."""
    if m.group(i) is not None:
        return ("MIN", m.group(i), None)
    return ("PR", m.group(i + 1), m.group(i + 2).strip())


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
    def ensure(kind, idx, formula):
        if kind == "MIN":
            nid = f"m{mult}_MIN{idx}"
            if nid not in g.nodes:
                g.add_node(nid, mult=mult, role="min", energy_eh=None,
                           energy_source="gfn2", xyz=gpath(f"MIN{idx}"))
        else:
            nid = f"m{mult}_PR{idx}"
            if nid not in g.nodes:
                g.add_node(nid, mult=mult, role="product", energy_eh=None,
                           energy_source="gfn2", xyz=gpath(f"PR{idx}"), label=formula)
        return nid

    for ln in open(rxnet):
        m = RXNET_RE.search(ln)
        if m:                                   # one-line format (unchanged behaviour)
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
            continue
        m = NATIVE_RE.match(ln)
        if not m:
            continue
        tsn, de = m.group(1), float(m.group(2))
        e1 = _native_end(m, 3); e2 = _native_end(m, 7)
        ts_id = f"m{mult}_ts{tsn}"
        n1 = ensure(*e1); n2 = ensure(*e2)
        if ts_id not in g.nodes:
            g.add_node(ts_id, mult=mult, role="ts", energy_eh=None,
                       energy_source="gfn2", xyz=gpath(f"ts{tsn}"))
        g.add_edge(n1, n2, "reaction", via_ts=ts_id, de_gfn2_kcal=de, energy_source="gfn2",
                   product_channel=(e1[0] == "PR" or e2[0] == "PR"))
    # native MINinfo (all minima, incl. those in no TS line)
    minfo = os.path.join(os.path.dirname(rxnet), "MINinfo")
    if os.path.isfile(minfo):
        for ln in open(minfo):
            mm = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s*$", ln)
            if mm:
                nid = ensure("MIN", mm.group(1), None)
                g.nodes[nid]["de_ll_kcal"] = float(mm.group(2))


def _mecp_geom(rec):
    p = rec.get("mecp_xyz") or os.path.join(rec.get("workdir", ""), "mecp.xyz")
    return p if os.path.isfile(p) else None


def dedupe_records(records, rmsd_tol=0.30, de_tol=0.5):
    """Merge records that are the same structure: heavy-atom RMSD with atom permutation < rmsd_tol
    and |dE| < de_tol kcal/mol (crossing_e_eh) and the same spin pair. Records without a readable
    geometry or energy are kept. Returns (kept_records, {kept_index_in_input: [merged input indices]})."""
    kept, merged = [], {}
    geoms = []
    for i, rec in enumerate(records):
        gp = _mecp_geom(rec)
        geoms.append(read_xyz(gp) if (gp and np is not None) else None)
    for i, rec in enumerate(records):
        dup = None
        if geoms[i] is not None and rec.get("crossing_e_eh") is not None:
            for k in kept:
                rk = records[k]
                if geoms[k] is None or rk.get("crossing_e_eh") is None:
                    continue
                if (int(rk["mult_a"]), int(rk["mult_b"])) != (int(rec["mult_a"]), int(rec["mult_b"])):
                    continue
                if abs(rk["crossing_e_eh"] - rec["crossing_e_eh"]) * H >= de_tol:
                    continue
                if rmsd_perm(geoms[k][0], geoms[k][1], geoms[i][0], geoms[i][1]) < rmsd_tol:
                    dup = k; break
        if dup is None:
            kept.append(i); merged[i] = []
        else:
            merged[dup].append(i)
    return [records[k] for k in kept], {k: merged[k] for k in kept}, kept


def add_mecps(records, g, gap_tol=1.0, rmsd_warn=2.0, dedupe=True, dedupe_rmsd=0.30,
              dedupe_de=0.5, perm_match=False, partner_roles=None):
    """Insert each accepted MECP as a node + two non-adiabatic crossing edges.
    dedupe: merge records that are the same structure (see dedupe_records).
    perm_match: use the permutation-aware RMSD for nearest-node matching.
    partner_roles: {mult: role} restricts the bridged node on that surface to that role."""
    added = []
    merged_map = {}
    if dedupe and len(records) > 1 and np is not None and linear_sum_assignment is not None:
        records, merged_map, kept_idx = dedupe_records(records, dedupe_rmsd, dedupe_de)
    else:
        kept_idx = list(range(len(records)))
    partner_roles = partner_roles or {}
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
                if mm in partner_roles:
                    surf = [n for n in surf if n.get("role") == partner_roles[mm]]
                node_id, r = best_match(ms, mx, surf, perm=perm_match)
                matches[mm] = (node_id, r)
                if node_id is not None:
                    g.add_edge(nid, node_id, "nonadiabatic_crossing", surface_mult=mm,
                               rmsd=round(r, 3), energy_source="dft",
                               weak=(r > rmsd_warn))
        added.append({"mecp": nid, "gap_kcal": gap, "accepted": accept,
                      "rel_kcal": g.nodes[nid]["rel_kcal"], "matches": matches,
                      "merged": list(merged_map.get(kept_idx[i], []))})
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
                               for m, v in mi["matches"].items())
                     + (f"; merged duplicate record(s) {mi['merged']}" if mi.get("merged") else ""))
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def build(nodes_json=None, spin_runs=None, mecp_records=None, ref_energy_eh=None,
          gap_tol=1.0, geom_subdir=None, dedupe=True, dedupe_rmsd=0.30, dedupe_de=0.5,
          perm_match=False, partner_roles=None):
    g = Graph(ref_energy_eh)
    if nodes_json:
        load_explicit_nodes(nodes_json, g)
    for mult, rundir in (spin_runs or []):
        load_spin_run(mult, rundir, g, geom_subdir)
    recs = [json.load(open(p)) for p in (mecp_records or [])]
    info = add_mecps(recs, g, gap_tol, dedupe=dedupe, dedupe_rmsd=dedupe_rmsd,
                     dedupe_de=dedupe_de, perm_match=perm_match, partner_roles=partner_roles)
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
    ap.add_argument("--no-dedupe", action="store_true", help="keep MECP records that are the same structure")
    ap.add_argument("--dedupe-rmsd", type=float, default=0.30, help="same-structure RMSD (A, permutation-aware)")
    ap.add_argument("--dedupe-de", type=float, default=0.5, help="same-structure |dE| (kcal/mol)")
    ap.add_argument("--perm-match", action="store_true", help="permutation-aware RMSD for node matching")
    ap.add_argument("--partner-role", action="append", default=[], metavar="MULT=ROLE",
                    help="bridge the MECP to the nearest node of ROLE (min|ts|product|reactant) on surface MULT")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return _selftest()
    spin_runs = [s.split("=", 1) for s in args.spin_run]
    roles = {int(a): b for a, b in (x.split("=", 1) for x in args.partner_role)}
    g, info = build(args.nodes, spin_runs, args.mecp, args.ref_energy_eh,
                    args.mecp_gap_tol, args.geom_subdir, dedupe=not args.no_dedupe,
                    dedupe_rmsd=args.dedupe_rmsd, dedupe_de=args.dedupe_de,
                    perm_match=args.perm_match, partner_roles=roles)
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
    # 4) native AutoMeKin RXNet incl. product channels + MINinfo; MECP dedupe (same structure) and
    #    partner-role bridging to a product node
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        open(os.path.join(d, "RXNet"), "w").write(
            " TS #   DE(kcal/mol)           Reaction path information\n"
            " ====   ============           =========================\n"
            "    1       53.5             MIN    1 <--->             MIN    4            \n"
            "    2       86.2             MIN    4 ---->    PR1:  CO + C4FeO4\n"
            "    3       87.4    PR2:  CO + C4FeO4 <--->    PR3:  CO + C4FeO4\n")
        open(os.path.join(d, "MINinfo"), "w").write("MIN #    DE(kcal/mol)\n    1       -1.125\n    2        0.000\n"
                                                    "    4       54.115\nConformational isomers are listed in the same line:\n1 2 \n")
        g = Graph()
        load_spin_run(3, d, g)
        roles = sorted(n["role"] for n in g.nodes.values())
        assert roles.count("min") == 3 and roles.count("product") == 3 and roles.count("ts") == 3, roles
        assert any(e.get("product_channel") for e in g.edges) and len(g.edges) == 3, g.edges
        assert g.nodes["m3_PR1"]["label"] == "CO + C4FeO4" and g.nodes["m3_MIN2"]["de_ll_kcal"] == 0.0
        print("  native RXNet incl. product channels + MINinfo OK")
        if np is not None and linear_sum_assignment is not None:
            def wx(name, pts):
                pth = os.path.join(d, name)
                open(pth, "w").write(f"{len(pts)}\n\n" + "".join(f"{e} {x} {y} {z}\n" for e, x, y, z in pts))
                return pth
            base = [("Fe", 0, 0, 0), ("C", 1.8, 0, 0), ("O", 2.9, 0, 0), ("C", -1.8, 0, 0), ("O", -2.9, 0, 0)]
            perm = [base[0], base[3], base[4], base[1], base[2]]            # same structure, atoms permuted
            far = [base[0], base[1], base[2], ("C", -3.5, 0.2, 0), ("O", -4.6, 0.2, 0)]
            m1, m2 = wx("m1.xyz", base), wx("m2.xyz", perm)
            recs = [{"mult_a": 1, "mult_b": 3, "gap_kcal": 0.01, "crossing_e_eh": -10.0, "mecp_xyz": m1},
                    {"mult_a": 1, "mult_b": 3, "gap_kcal": 0.02, "crossing_e_eh": -10.0 + 0.0003, "mecp_xyz": m2},
                    {"mult_a": 1, "mult_b": 3, "gap_kcal": 0.02, "crossing_e_eh": -10.0 + 0.002, "mecp_xyz": m2}]
            g = Graph(ref_energy_eh=-10.0)
            g.add_node("s_min", mult=1, role="min", energy_eh=-10.1, energy_source="dft", xyz=m1)
            g.add_node("t_min", mult=3, role="min", energy_eh=-10.05, energy_source="dft", xyz=m2)
            g.add_node("t_prod", mult=3, role="product", energy_eh=None, energy_source="gfn2", xyz=wx("p.xyz", far))
            info = add_mecps(recs, g, perm_match=True, partner_roles={3: "product"})
            mecps = [n for n in g.nodes.values() if n["role"] == "mecp"]
            # record 2 (same structure, |dE| 0.19 kcal) merged; record 3 (|dE| 1.3 kcal) kept
            assert len(mecps) == 2 and info[0]["merged"] == [1], (len(mecps), info)
            tgt = {e["surface_mult"]: e["target"] for e in g.edges if e["source"] == mecps[0]["id"]}
            assert tgt == {1: "s_min", 3: "t_prod"}, tgt
            print("  MECP dedupe (permuted duplicate merged, distinct energy kept) + partner-role bridge OK")
    print("M3 selftest PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
