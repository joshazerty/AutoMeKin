#!/usr/bin/env python3
"""
dissociation_scan.py - M1 (ligand-loss channels): spin crossings with no saddle point.

WHY
---
spin_enumerate.py (M1) checks vertical spin gaps only at stationary points
(minima and TSs). A spin crossing that lies on a barrierless dissociation
channel - the textbook case is 1[Fe(CO)5] -> 3[Fe(CO)4] + CO - has no saddle
point on that channel, so it is never presented to M1 and cannot be flagged.

WHAT
----
For each minimum and each metal-ligand bond (ligand = connected fragment left
when the metal is removed; donor = ligand atom bonded to the metal), run a
GFN2-xTB relaxed scan of the M-donor distance on the base multiplicity and on
each Delta S = 1 neighbour (mult +/- 2). Flag the channel when
  (a) the two relaxed curves cross inside the scan (sign change of
      E_alt - E_base), or
  (b) the smallest |E_alt - E_base| along the scan is within the M1 window.
The candidate record has the same schema as spin_enumerate.py output, with
the seed geometry = the scan frame closest to the crossing, so locate_mecp.py
(M2) consumes it unchanged. Symmetry-equivalent ligands (same scan profile to
0.1 kcal/mol) are reported once.

GFN2 decides nothing here (as for M1): it only proposes where to look. The
crossing energy and geometry come from M2 at DFT.

PIPELINE
--------
  1. dissociation_scan.py      GFN2 relaxed M-L scans per spin -> mecp_candidates_dissoc.json
  2. dissociation_seed_dft.py  DFT single points on the scan frames; seed at the DFT
                               crossing -> mecp_candidates_dissoc_dft.json
  3. locate_mecp.py --from-candidates mecp_candidates_dissoc_dft.json   (M2, unchanged)
Step 2 is required whenever GFN2 may misorder the spin states of the unsaturated
fragment (Fe(CO)4: GFN2 puts the singlet below the triplet, so the GFN2 curves never
cross). The GFN2 candidates of step 1 are not meant to be passed to M2 directly.

Each candidate record carries `scan_files` ({mult: path}) so step 2 can be run in batch
mode on the step-1 JSON.

Validation: 1Fe(CO)5 -> 3Fe(CO)4 + CO (Harvey & Aschi, Faraday Discuss. 2003, 124, 129);
see examples/feco5_co_loss/.
Requires: the `xtb` binary (6.x, GFN2) on PATH or via --xtb.
"""
import os, sys, json, argparse, subprocess, tempfile, shutil
import numpy as np

HARTREE_KCAL = 627.5094740631
TM = {"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Y","Zr","Nb","Mo","Tc",
      "Ru","Rh","Pd","Ag","Cd","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"}
COV = {"H":0.31,"B":0.84,"C":0.76,"N":0.71,"O":0.66,"F":0.57,"P":1.07,"S":1.05,
       "Cl":1.02,"Br":1.20,"I":1.39,"Si":1.11,
       "Fe":1.32,"Co":1.26,"Ni":1.24,"Mn":1.39,"Cr":1.39,"Ru":1.46,"Rh":1.42,
       "Ir":1.41,"Pt":1.36,"Pd":1.39,"Cu":1.32,"Mo":1.54,"W":1.62,"Os":1.44,
       "Sc":1.70,"Ti":1.60,"V":1.53,"Zn":1.22,"Ag":1.45,"Au":1.36,"Re":1.51,"As":1.19,"Se":1.20}


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0].split()[0])
    s, x = [], []
    for ln in L[2:2+n]:
        p = ln.split(); s.append(p[0]); x.append([float(v) for v in p[1:4]])
    return s, np.array(x)


def write_xyz(path, syms, X, comment=""):
    with open(path, "w") as fh:
        fh.write(f"{len(syms)}\n{comment}\n")
        for a, r in zip(syms, X):
            fh.write(f"{a:3s} {r[0]: .8f} {r[1]: .8f} {r[2]: .8f}\n")


def bonds(syms, X, scale=1.25):
    n = len(syms); B = set()
    for i in range(n):
        for j in range(i+1, n):
            if np.linalg.norm(X[i]-X[j]) < scale*(COV.get(syms[i],1.5)+COV.get(syms[j],1.5)):
                B.add((i, j))
    return B


def ligands(syms, X):
    """[(metal, donor, ligand_atom_indices)] for monodentate metal-ligand bonds."""
    B = bonds(syms, X)
    metals = [i for i, a in enumerate(syms) if a in TM]
    adj = {i: set() for i in range(len(syms))}
    for i, j in B:
        if i in metals or j in metals:
            continue
        adj[i].add(j); adj[j].add(i)
    seen, frags = set(), []
    for i in range(len(syms)):
        if i in metals or i in seen:
            continue
        st, comp = [i], set()
        while st:
            k = st.pop()
            if k in comp: continue
            comp.add(k); st.extend(adj[k] - comp)
        seen |= comp; frags.append(sorted(comp))
    out = []
    for m in metals:
        for fr in frags:
            donors = [d for d in fr if (min(m, d), max(m, d)) in B]
            if len(donors) == 1:              # monodentate only (kappa-1 loss)
                out.append((m, donors[0], fr))
    return out


def xtb_opt(syms, X, charge, mult, constrain=None, xtb="xtb", level="tight"):
    """GFN2 optimisation, optionally with a fixed distance (i, j, r_A).
    Returns (E_Eh, X_opt)."""
    d = tempfile.mkdtemp(prefix="dscan_")
    try:
        write_xyz(os.path.join(d, "in.xyz"), syms, X)
        args = [xtb, "in.xyz", "--gfn", "2", "--chrg", str(charge),
                "--uhf", str(mult-1), "--opt", level]
        if constrain:
            i, j, r = constrain
            open(os.path.join(d, "xc.inp"), "w").write(
                f"$constrain\n force constant=5.0\n distance: {i+1}, {j+1}, {r:.4f}\n$end\n")
            args += ["--input", "xc.inp"]
        env = dict(os.environ, OMP_NUM_THREADS="1", OMP_STACKSIZE="1G")
        p = subprocess.run(args, cwd=d, capture_output=True, text=True, env=env)
        E = None
        for ln in p.stdout.splitlines():
            if "TOTAL ENERGY" in ln:
                E = float(ln.split()[3])
        if E is None or not os.path.isfile(os.path.join(d, "xtbopt.xyz")):
            raise RuntimeError(f"xtb failed (mult {mult}): {p.stdout[-400:]}")
        _, Xo = read_xyz(os.path.join(d, "xtbopt.xyz"))
        return E, Xo
    finally:
        shutil.rmtree(d, ignore_errors=True)


def scan_channel(syms, X0, charge, m, dnr, mults, dr=0.1, stretch=2.5, xtb="xtb"):
    """Relaxed M-donor scan on each multiplicity, each surface propagated from
    its own previous frame. Returns r grid, {mult: [E]}, {mult: [X]}."""
    r0 = float(np.linalg.norm(X0[m]-X0[dnr]))
    grid = np.round(np.arange(r0, r0+stretch+1e-9, dr), 4)
    E, G = {}, {}
    for mu in mults:
        Xc = X0.copy(); E[mu], G[mu] = [], []
        for r in grid:
            e, Xc = xtb_opt(syms, Xc, charge, mu, constrain=(m, dnr, r), xtb=xtb)
            E[mu].append(e); G[mu].append(Xc.copy())
    return grid, E, G


def analyse(grid, Eb, Ea):
    """Gap profile (kcal), crossing radius (linear interpolation) and seed index."""
    g = (np.array(Ea)-np.array(Eb))*HARTREE_KCAL
    cross = None
    for k in range(len(g)-1):
        if g[k] == 0 or g[k]*g[k+1] < 0:
            t = g[k]/(g[k]-g[k+1]); cross = float(grid[k]+t*(grid[k+1]-grid[k])); break
    seed = int(np.argmin(np.abs(g))) if cross is None else \
        int(np.argmin(np.abs(grid-cross)))
    return g, cross, seed


def run(xyz, charge, base_mult, window=30.0, dr=0.1, stretch=2.5, outdir="dscan",
        out_json="mecp_candidates_dissoc.json", xtb="xtb", node=None):
    os.makedirs(outdir, exist_ok=True)
    syms, X = read_xyz(xyz)
    node = node or os.path.splitext(os.path.basename(xyz))[0]
    alts = [mu for mu in (base_mult-2, base_mult+2) if mu >= 1]
    # relax the minimum on the base surface first (scan starts at equilibrium)
    _, X = xtb_opt(syms, X, charge, base_mult, xtb=xtb)
    chans = ligands(syms, X)
    print(f"{node}: {len(chans)} monodentate M-L channel(s); base mult {base_mult}, alt {alts}")
    profiles, cands, report = [], [], []
    for (m, dnr, fr) in chans:
        lig = "".join(sorted(syms[i] for i in fr))
        grid, E, G = scan_channel(syms, X, charge, m, dnr, [base_mult]+alts, dr, stretch, xtb)
        # symmetry-equivalence: same equilibrium M-donor distance (0.02 A) and same
        # base-surface asymptote (0.5 kcal/mol)
        prof_key = (float(grid[0]), float((E[base_mult][-1]-E[base_mult][0])*HARTREE_KCAL))
        if any(abs(prof_key[0]-p[0]) < 0.02 and abs(prof_key[1]-p[1]) < 0.5 for p in profiles):
            report.append({"channel": f"{syms[m]}{m+1}-{syms[dnr]}{dnr+1} ({lig})",
                           "equivalent_to_earlier": True}); continue
        profiles.append(prof_key)
        scan_tag = f"{node}_L{syms[dnr]}{dnr+1}"
        scan_files = {}
        for mu in [base_mult]+alts:     # full scan trajectories (for DFT seed refinement)
            scan_files[mu] = os.path.abspath(os.path.join(outdir, f"{scan_tag}_mult{mu}_scan.xyz"))
            with open(scan_files[mu], "w") as fh:
                for r, Xk, e in zip(grid, G[mu], E[mu]):
                    fh.write(f"{len(syms)}\nr={r:.4f} E_gfn2={e:.8f} mult={mu} metal={m+1} donor={dnr+1}\n")
                    for at, p in zip(syms, Xk):
                        fh.write(f"{at:3s} {p[0]: .8f} {p[1]: .8f} {p[2]: .8f}\n")
        for a in alts:
            g, cross, k = analyse(grid, E[base_mult], E[a])
            tag = f"{node}_L{syms[dnr]}{dnr+1}_m{base_mult}to{a}"
            rec = {"channel": f"{syms[m]}{m+1}-{syms[dnr]}{dnr+1} ({lig})", "alt_mult": a,
                   "r_grid_A": grid.tolist(), "gap_kcal": np.round(g, 2).tolist(),
                   "E_base_rel_kcal": np.round((np.array(E[base_mult])-E[base_mult][0])*HARTREE_KCAL, 2).tolist(),
                   "E_alt_rel_kcal": np.round((np.array(E[a])-E[base_mult][0])*HARTREE_KCAL, 2).tolist(),
                   "crossing_r_A": cross, "min_abs_gap_kcal": float(np.min(np.abs(g)))}
            report.append(rec)
            if cross is not None or np.min(np.abs(g)) <= window:
                seed_xyz = os.path.abspath(os.path.join(outdir, tag + "_seed.xyz"))
                write_xyz(seed_xyz, syms, G[base_mult][k],
                          f"dissociation-channel seed r={grid[k]:.3f} A gap={g[k]:.2f} kcal (GFN2)")
                cands.append({"node": tag, "xyz": seed_xyz, "charge": charge,
                              "base_mult": base_mult, "alt_mult": a,
                              "dE_vert_kcal": round(float(g[k]), 2),
                              "kind": "dissociation_channel",
                              "channel": rec["channel"], "r_seed_A": float(grid[k]),
                              "crossing_r_A": cross,
                              "scan_tag": scan_tag,
                              "scan_files": {str(base_mult): scan_files[base_mult],
                                             str(a): scan_files[a]},
                              "note": ("relaxed curves cross on ligand-loss channel" if cross is not None
                                       else "relaxed curves within window on ligand-loss channel")})
            # asymptote fragment geometry (last frame, alt surface) for ordinary M1
            write_xyz(os.path.join(outdir, tag + "_asymptote_alt.xyz"), syms, G[a][-1],
                      f"last scan frame, mult {a}")
        write_xyz(os.path.join(outdir, f"{node}_L{syms[dnr]}{dnr+1}_asymptote_base.xyz"),
                  syms, G[base_mult][-1], f"last scan frame, mult {base_mult}")
    json.dump({"window_kcal": window, "method": "gfn2", "level": "gfn2-flag",
               "stage": "M1-dissociation", "candidates": cands}, open(out_json, "w"), indent=2)
    json.dump(report, open(os.path.join(outdir, "scan_report.json"), "w"), indent=1)
    print(f"{len(cands)} dissociation-channel candidate(s) -> {out_json}")
    for c in cands:
        print(f"  {c['node']}: {c['channel']} {c['base_mult']}<->{c['alt_mult']} "
              f"r_seed={c['r_seed_A']:.2f} A gap={c['dE_vert_kcal']} kcal cross_r={c['crossing_r_A']}")
    return cands, report


def _selftest():
    g = np.array([2.0, 2.1, 2.2]); c = analyse(g, [0, 0, 0], [5/HARTREE_KCAL, 1/HARTREE_KCAL, -3/HARTREE_KCAL])
    assert c[1] is not None and 2.1 < c[1] < 2.2 and c[2] in (1, 2), c
    s = ["Fe", "C", "O", "C", "O"]
    X = np.array([[0, 0, 0], [1.8, 0, 0], [2.95, 0, 0], [-1.8, 0, 0], [-2.95, 0, 0]], float)
    L = ligands(s, X); assert len(L) == 2 and all(len(f) == 2 for _, _, f in L), L
    print("dissociation_scan selftest PASS"); return 0


def main():
    ap = argparse.ArgumentParser(description="M1 extension: spin crossings on ligand-loss channels")
    ap.add_argument("xyz", nargs="?"); ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int); ap.add_argument("--window", type=float, default=30.0)
    ap.add_argument("--dr", type=float, default=0.1); ap.add_argument("--stretch", type=float, default=2.5)
    ap.add_argument("--outdir", default="dscan"); ap.add_argument("--out", default="mecp_candidates_dissoc.json")
    ap.add_argument("--xtb", default="xtb"); ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest: return _selftest()
    if a.xyz is None or a.mult is None:
        ap.error("provide a minimum geometry (xyz) and --mult (its multiplicity)")
    run(a.xyz, a.charge, a.mult, a.window, a.dr, a.stretch, a.outdir, a.out, a.xtb)


if __name__ == "__main__":
    sys.exit(main() or 0)
