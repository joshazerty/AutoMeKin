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
each Delta S = 1 neighbour (mult +/- 2). Which bonds are scanned:
  kappa-1 (one donor atom), even-electron fragment (CO, alkene bound through one
      atom, phosphine, N-heterocycle ...): scanned, variable = M-donor distance.
  eta-2 (two donor atoms bonded to each other, even-electron fragment, e.g. an
      alkene bound side-on): scanned, variable = distance from the metal to the
      midpoint of the two donors. Both M-donor distances are constrained together,
      d_i(r) = sqrt(r**2 + d_i0**2 - r0**2), with d_i0 and r0 taken from the
      starting geometry, so the ligand keeps its side-on orientation.
  X-type (odd-electron neutral fragment: H, alkyl, aryl, allyl ...): skipped unless
      --include-x. Homolysis of M-X is a bond stretch, not ligand loss.
  everything else (two or more donors that are not an eta-2 pair, e.g. chelates,
      kappa-3 Tp, the bis(carbene)borate): skipped, and written to the output as a
      skipped bond with the reason; nothing is dropped silently.
Flag the channel when
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
Batch use: several minima in one call, with an energy window (--emax, default 20
kcal/mol) measured from the lowest minimum supplied. Energies come from --energies
(file with `name  energy_kcal` per line, e.g. the network low-level energies) or,
if absent, from the GFN2 energy of each relaxed minimum. Minima outside the window
are listed in the summary and not scanned.
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


ATOMIC_NUMBER = {"H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Si": 14, "P": 15, "S": 16,
                 "Cl": 17, "As": 33, "Se": 34, "Br": 35, "I": 53}


def is_x_type(syms, fragment):
    """True if the neutral free fragment has an odd electron count (H, alkyl, aryl, allyl ...).
    Cleaving M-X then leaves a radical: that is a bond homolysis, not a ligand-loss channel.
    Even-electron fragments (CO, alkene, phosphine, N-heterocycle ...) are L-type."""
    try:
        return sum(ATOMIC_NUMBER[syms[i]] for i in fragment) % 2 == 1
    except KeyError as e:
        raise KeyError(f"no atomic number for {e}; add it to ATOMIC_NUMBER") from None


def metal_bonds(syms, X, include_x=False):
    """Classify every ligand fragment bonded to a metal (one record per metal and fragment).
    Record: metal, donors (sorted atom indices bonded to the metal), fragment (atom indices),
    formula, kind in {kappa1, eta2, polydentate}, x_type, scan (bool), reason (why not scanned)."""
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
            donors = sorted(d for d in fr if (min(m, d), max(m, d)) in B)
            if not donors:
                continue
            xt = is_x_type(syms, fr)
            formula = "".join(sorted(syms[i] for i in fr))
            if len(donors) == 1:
                kind = "kappa1"
            elif len(donors) == 2 and donors[1] in adj[donors[0]]:
                kind = "eta2"
            else:
                kind = "polydentate"
            reason = ""
            if kind == "polydentate":
                reason = (f"{len(donors)} donor atoms bonded to the metal and not an adjacent pair "
                          "(chelate / kappa-n ligand): not scanned")
            elif xt and not include_x:
                reason = ("odd-electron fragment (X-type): M-X cleavage is a homolysis, "
                          "not ligand loss (use --include-x to scan it)")
            out.append({"metal": m, "donors": donors, "fragment": fr, "formula": formula,
                        "kind": kind, "x_type": xt, "scan": reason == "", "reason": reason})
    return out


def ligands(syms, X, include_x=False):
    """[(metal, donor, ligand_atom_indices)] for monodentate (kappa-1) metal-ligand bonds.
    X-type (odd-electron) fragments are left out unless include_x. eta-2 ligands are
    returned by metal_bonds(), not here."""
    return [(r["metal"], r["donors"][0], r["fragment"]) for r in metal_bonds(syms, X, include_x)
            if r["kind"] == "kappa1" and r["scan"]]


def relaxed_away(recs_in, recs_out):
    """Fragments bonded to the metal in the input geometry but no longer bonded after the GFN2
    relaxation (dissociated or weakly bound): same metal and same atom set, no record in recs_out."""
    have = {(r["metal"], tuple(r["fragment"])) for r in recs_out}
    return [r for r in recs_in if (r["metal"], tuple(r["fragment"])) not in have]


def channel_label(syms, rec):
    m = rec["metal"]
    if rec["kind"] == "eta2":
        d = "+".join(f"{syms[i]}{i+1}" for i in rec["donors"])
        return f"{syms[m]}{m+1}-eta2({d}) ({rec['formula']})"
    d = rec["donors"][0]
    return f"{syms[m]}{m+1}-{syms[d]}{d+1} ({rec['formula']})"


def channel_tag(syms, rec):
    return "L" + "".join(f"{syms[i]}{i+1}" for i in rec["donors"])


def scan_variable(X, rec):
    """Scan coordinate r0 (A) of a channel: M-donor distance (kappa1) or M-midpoint distance (eta2)."""
    m = rec["metal"]
    if rec["kind"] == "eta2":
        return float(np.linalg.norm(X[m] - X[rec["donors"]].mean(axis=0)))
    return float(np.linalg.norm(X[m] - X[rec["donors"][0]]))


def constraints_at(X0, rec, r):
    """Distance constraints [(i, j, r_A), ...] that put the scan variable at r.
    kappa1: M-donor = r. eta2: both M-donor distances move together,
    d_i = sqrt(r**2 + d_i0**2 - r0**2) (r0 = M-midpoint distance of X0), so the
    midpoint distance follows r while the side-on geometry is preserved."""
    m = rec["metal"]
    if rec["kind"] != "eta2":
        return [(m, rec["donors"][0], r)]
    r0 = scan_variable(X0, rec); out = []
    for d in rec["donors"]:
        d0 = float(np.linalg.norm(X0[m] - X0[d]))
        out.append((m, d, float(np.sqrt(max(r**2 + d0**2 - r0**2, 0.25)))))
    return out


def xtb_opt(syms, X, charge, mult, constrain=None, xtb="xtb", level="tight"):
    """GFN2 optimisation, optionally with fixed distance(s): one (i, j, r_A) or a list of them.
    Returns (E_Eh, X_opt)."""
    d = tempfile.mkdtemp(prefix="dscan_")
    try:
        write_xyz(os.path.join(d, "in.xyz"), syms, X)
        args = [xtb, "in.xyz", "--gfn", "2", "--chrg", str(charge),
                "--uhf", str(mult-1), "--opt", level]
        if constrain:
            cl = [constrain] if isinstance(constrain[0], (int, np.integer)) else list(constrain)
            body = "".join(f" distance: {i+1}, {j+1}, {r:.4f}\n" for i, j, r in cl)
            open(os.path.join(d, "xc.inp"), "w").write(
                f"$constrain\n force constant=5.0\n{body}$end\n")
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


def scan_channel(syms, X0, charge, rec, mults, dr=0.1, stretch=2.5, xtb="xtb"):
    """Relaxed scan of one channel on each multiplicity, each surface propagated from
    its own previous frame. Returns r grid, {mult: [E]}, {mult: [X]}."""
    r0 = scan_variable(X0, rec)
    grid = np.round(np.arange(r0, r0+stretch+1e-9, dr), 4)
    E, G = {}, {}
    for mu in mults:
        Xc = X0.copy(); E[mu], G[mu] = [], []
        for r in grid:
            e, Xc = xtb_opt(syms, Xc, charge, mu, constrain=constraints_at(X0, rec, float(r)), xtb=xtb)
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
        out_json="mecp_candidates_dissoc.json", xtb="xtb", node=None, include_x=False):
    os.makedirs(outdir, exist_ok=True)
    syms, X = read_xyz(xyz)
    node = node or os.path.splitext(os.path.basename(xyz))[0]
    alts = [mu for mu in (base_mult-2, base_mult+2) if mu >= 1]
    recs_in = metal_bonds(syms, X, include_x)
    # relax the minimum on the base surface first (scan starts at equilibrium)
    _, X = xtb_opt(syms, X, charge, base_mult, xtb=xtb)
    recs = metal_bonds(syms, X, include_x)
    chans = [r for r in recs if r["scan"]]
    skipped = [r for r in recs if not r["scan"]]
    nk = sum(r["kind"] == "kappa1" for r in chans); ne = sum(r["kind"] == "eta2" for r in chans)
    print(f"{node}: {len(chans)} channel(s) to scan ({nk} kappa-1, {ne} eta-2), {len(skipped)} metal bond(s) "
          f"skipped; base mult {base_mult}, alt {alts}")
    report, cands = [], []
    for r in skipped:
        lab = channel_label(syms, r)
        print(f"  [skip] {lab}: {r['reason']}")
        report.append({"channel": lab, "skipped": True, "kind": r["kind"], "x_type": r["x_type"],
                       "donors": [f"{syms[i]}{i+1}" for i in r["donors"]], "formula": r["formula"],
                       "reason": r["reason"]})
    for r in relaxed_away(recs_in, recs):
        lab = channel_label(syms, r)
        why = ("bonded to the metal in the input geometry but not after the GFN2 relaxation "
               "(dissociated or weakly bound; M-donor distance above 1.25 x covalent radii): not scanned")
        print(f"  [skip] {lab}: {why}")
        report.append({"channel": lab, "skipped": True, "kind": r["kind"], "x_type": r["x_type"],
                       "donors": [f"{syms[i]}{i+1}" for i in r["donors"]], "formula": r["formula"],
                       "reason": why})
    profiles = []
    for rec in chans:
        m = rec["metal"]; lab = channel_label(syms, rec)
        grid, E, G = scan_channel(syms, X, charge, rec, [base_mult]+alts, dr, stretch, xtb)
        # symmetry-equivalence: same starting scan variable (0.02 A) and same base-surface asymptote
        # (0.5 kcal/mol), same kind of channel
        prof_key = (rec["kind"], float(grid[0]), float((E[base_mult][-1]-E[base_mult][0])*HARTREE_KCAL))
        if any(p[0] == prof_key[0] and abs(prof_key[1]-p[1]) < 0.02 and abs(prof_key[2]-p[2]) < 0.5 for p in profiles):
            report.append({"channel": lab, "equivalent_to_earlier": True, "kind": rec["kind"]}); continue
        profiles.append(prof_key)
        ctag = channel_tag(syms, rec)
        scan_tag = f"{node}_{ctag}"
        dn = "+".join(str(i+1) for i in rec["donors"])
        scan_files = {}
        for mu in [base_mult]+alts:     # full scan trajectories (for DFT seed refinement)
            scan_files[mu] = os.path.abspath(os.path.join(outdir, f"{scan_tag}_mult{mu}_scan.xyz"))
            with open(scan_files[mu], "w") as fh:
                for r, Xk, e in zip(grid, G[mu], E[mu]):
                    fh.write(f"{len(syms)}\nr={r:.4f} E_gfn2={e:.8f} mult={mu} metal={m+1} donor={dn} kind={rec['kind']}\n")
                    for at, p in zip(syms, Xk):
                        fh.write(f"{at:3s} {p[0]: .8f} {p[1]: .8f} {p[2]: .8f}\n")
        for a in alts:
            g, cross, k = analyse(grid, E[base_mult], E[a])
            tag = f"{scan_tag}_m{base_mult}to{a}"
            row = {"channel": lab, "kind": rec["kind"], "x_type": rec["x_type"], "alt_mult": a,
                   "r_grid_A": grid.tolist(), "gap_kcal": np.round(g, 2).tolist(),
                   "E_base_rel_kcal": np.round((np.array(E[base_mult])-E[base_mult][0])*HARTREE_KCAL, 2).tolist(),
                   "E_alt_rel_kcal": np.round((np.array(E[a])-E[base_mult][0])*HARTREE_KCAL, 2).tolist(),
                   "crossing_r_A": cross, "min_abs_gap_kcal": float(np.min(np.abs(g)))}
            report.append(row)
            if cross is not None or np.min(np.abs(g)) <= window:
                seed_xyz = os.path.abspath(os.path.join(outdir, tag + "_seed.xyz"))
                write_xyz(seed_xyz, syms, G[base_mult][k],
                          f"dissociation-channel seed r={grid[k]:.3f} A gap={g[k]:.2f} kcal (GFN2)")
                cands.append({"node": tag, "xyz": seed_xyz, "charge": charge,
                              "base_mult": base_mult, "alt_mult": a,
                              "dE_vert_kcal": round(float(g[k]), 2),
                              "kind": "dissociation_channel",
                              "hapticity": rec["kind"], "donors": [f"{syms[i]}{i+1}" for i in rec["donors"]],
                              "fragment_formula": rec["formula"],
                              "channel": lab, "r_seed_A": float(grid[k]),
                              "crossing_r_A": cross,
                              "scan_tag": scan_tag,
                              "scan_files": {str(base_mult): scan_files[base_mult],
                                             str(a): scan_files[a]},
                              "note": ("relaxed curves cross on ligand-loss channel" if cross is not None
                                       else "relaxed curves within window on ligand-loss channel")})
            # asymptote fragment geometry (last frame, alt surface) for ordinary M1
            write_xyz(os.path.join(outdir, tag + "_asymptote_alt.xyz"), syms, G[a][-1],
                      f"last scan frame, mult {a}")
        write_xyz(os.path.join(outdir, f"{scan_tag}_asymptote_base.xyz"),
                  syms, G[base_mult][-1], f"last scan frame, mult {base_mult}")
    skipped_json = [x for x in report if x.get("skipped")]
    json.dump({"window_kcal": window, "method": "gfn2", "level": "gfn2-flag",
               "stage": "M1-dissociation", "include_x": include_x,
               "skipped_bonds": skipped_json, "candidates": cands}, open(out_json, "w"), indent=2)
    json.dump(report, open(os.path.join(outdir, "scan_report.json"), "w"), indent=1)
    print(f"{len(cands)} dissociation-channel candidate(s) -> {out_json}")
    for c in cands:
        print(f"  {c['node']}: {c['channel']} {c['base_mult']}<->{c['alt_mult']} "
              f"r_seed={c['r_seed_A']:.2f} A gap={c['dE_vert_kcal']} kcal cross_r={c['crossing_r_A']}")
    return cands, report


def select_window(energies, emax):
    """{name: energy_kcal} -> (in_window names, {name: energy relative to the lowest}). The window is
    measured from the lowest minimum supplied; a minimum at exactly emax is inside."""
    lo = min(energies.values())
    rel = {k: v - lo for k, v in energies.items()}
    return [k for k, v in rel.items() if v <= emax + 1e-9], rel


def _batch_worker(args):
    (xyz, node, charge, mult, window, dr, stretch, d, xtb, include_x) = args
    import contextlib, io
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            cands, rep = run(xyz, charge, mult, window, dr, stretch, d, os.path.join(d, "mecp_candidates_dissoc.json"),
                             xtb, node, include_x)
        open(os.path.join(d, "log.txt"), "w").write(buf.getvalue())
        return node, "ok", cands, [x for x in rep if x.get("skipped")], \
            sum(1 for x in rep if x.get("alt_mult") is not None and not x.get("skipped")), ""
    except Exception as ex:                      # xtb SCF failure etc.: report, do not stop the batch
        open(os.path.join(d, "log.txt"), "w").write(buf.getvalue() + "\nFAILED: " + repr(ex)[:500])
        return node, "failed", [], [], 0, repr(ex)[:300]


def run_batch(xyzs, charge, base_mult, energies=None, emax=20.0, window=30.0, dr=0.1, stretch=2.5,
              outdir="dscan", out_json="mecp_candidates_dissoc.json", xtb="xtb", include_x=False, jobs=1):
    """Scan every minimum within emax (kcal/mol) of the lowest one supplied. energies = {name: kcal};
    when None the GFN2 energy of each relaxed minimum is used."""
    os.makedirs(outdir, exist_ok=True)
    names = {os.path.splitext(os.path.basename(x))[0]: x for x in xyzs}
    if energies is None:
        energies = {}
        for n, x in names.items():
            syms, X = read_xyz(x)
            e, _ = xtb_opt(syms, X, charge, base_mult, xtb=xtb)
            energies[n] = e*HARTREE_KCAL
    else:
        missing = [n for n in names if n not in energies]
        if missing:
            raise KeyError(f"no energy for {missing}")
        energies = {n: energies[n] for n in names}
    inside, rel = select_window(energies, emax)
    print(f"{len(inside)} of {len(names)} minima within {emax} kcal/mol of the lowest ({min(energies, key=energies.get)})")
    work = [(names[n], n, charge, base_mult, window, dr, stretch, os.path.join(outdir, n), xtb, include_x)
            for n in names if n in inside]
    for w in work: os.makedirs(w[7], exist_ok=True)
    if jobs > 1:
        from multiprocessing import Pool
        with Pool(jobs) as p: res = p.map(_batch_worker, work)
    else:
        res = [_batch_worker(w) for w in work]
    summ, allc = [], []
    for n in names:
        row = {"node": n, "energy_kcal": round(energies[n], 3), "rel_kcal": round(rel[n], 3),
               "in_window": n in inside}
        r = next((x for x in res if x[0] == n), None)
        if r:
            row.update({"status": r[1], "n_candidates": len(r[2]), "n_skipped_bonds": len(r[3]),
                        "skipped_bonds": r[3], "error": r[5]}); allc += r[2]
        else:
            row["status"] = "outside energy window (not scanned)"
        summ.append(row)
    json.dump({"emax_kcal": emax, "window_kcal": window, "include_x": include_x,
               "minima": summ}, open(os.path.join(outdir, "scan_batch_summary.json"), "w"), indent=1)
    json.dump({"window_kcal": window, "method": "gfn2", "level": "gfn2-flag", "stage": "M1-dissociation",
               "include_x": include_x, "candidates": allc}, open(out_json, "w"), indent=2)
    print(f"{len(allc)} candidate(s) from {len(work)} minima -> {out_json}")
    return allc, summ


def read_energies(path):
    """`name  energy_kcal` per line (blank lines and # comments ignored)."""
    out = {}
    for ln in open(path):
        p = ln.split("#")[0].split()
        if len(p) >= 2: out[p[0]] = float(p[1])
    return out


def _selftest():
    g = np.array([2.0, 2.1, 2.2]); c = analyse(g, [0, 0, 0], [5/HARTREE_KCAL, 1/HARTREE_KCAL, -3/HARTREE_KCAL])
    assert c[1] is not None and 2.1 < c[1] < 2.2 and c[2] in (1, 2), c
    s = ["Fe", "C", "O", "C", "O"]
    X = np.array([[0, 0, 0], [1.8, 0, 0], [2.95, 0, 0], [-1.8, 0, 0], [-2.95, 0, 0]], float)
    L = ligands(s, X); assert len(L) == 2 and all(len(f) == 2 for _, _, f in L), L
    # --- X-type: odd-electron fragments (H, ethyl) skipped by default, kept with include_x
    assert is_x_type(["H"], [0]) and not is_x_type(["C", "O"], [0, 1])
    assert is_x_type(["C", "C", "H", "H", "H", "H", "H"], list(range(7)))          # ethyl
    assert not is_x_type(["C", "C", "H", "H", "H", "H"], list(range(6)))           # ethylene
    # Fe(CO)(H)(side-on C2H4) + a chelate (two non-adjacent N donors on one fragment)
    s2 = ["Fe", "C", "O", "H", "C", "C", "H", "H", "H", "H", "N", "C", "C", "C", "N"]   # CO, H, C2H4, N-C-C-C-N chain
    X2 = np.array([[0, 0, 0], [1.8, 0, 0], [2.95, 0, 0], [0, 1.55, 0],
                   [-0.1, -0.1, 2.0], [0.1, 1.2, 2.05],        # alkene C-C along y, 2.0 A above Fe
                   [-0.5, -0.7, 2.6], [-0.9, 0.0, 2.6], [0.7, 1.8, 2.6], [0.9, 1.0, 2.6],
                   [-1.5, -1.4, -0.3], [-2.5, -1.2, -0.3], [-3.4, 0.0, -0.3], [-2.5, 1.2, -0.3],
                   [-1.5, 1.4, -0.3]], float)
    recs2 = metal_bonds(s2, X2)
    kinds = {r["kind"] for r in recs2}
    assert kinds == {"kappa1", "eta2", "polydentate"}, kinds
    by_formula = {r["formula"]: r for r in recs2}
    assert by_formula["H"]["x_type"] and not by_formula["H"]["scan"] and "X-type" in by_formula["H"]["reason"]
    assert by_formula["CCHHHH"]["kind"] == "eta2" and by_formula["CCHHHH"]["scan"]
    pol = [r for r in recs2 if r["kind"] == "polydentate"][0]
    assert not pol["scan"] and "chelate" in pol["reason"] and len(pol["donors"]) == 2
    assert [r["formula"] for r in recs2 if r["scan"]] == ["CO", "CCHHHH"] or \
        sorted(r["formula"] for r in recs2 if r["scan"]) == ["CCHHHH", "CO"]
    assert len(ligands(s2, X2)) == 1 and len(ligands(s2, X2, include_x=True)) == 2
    # --- eta-2 constraints: both distances follow the midpoint distance, geometry at r0 is reproduced
    e2 = by_formula["CCHHHH"]; r0 = scan_variable(X2, e2)
    c0 = constraints_at(X2, e2, r0)
    for (i, j, d), dd in zip(c0, e2["donors"]):
        assert abs(d - np.linalg.norm(X2[0] - X2[dd])) < 1e-9
    c1 = constraints_at(X2, e2, r0 + 1.0)
    assert all(b[2] > a[2] for a, b in zip(c0, c1)) and len(c1) == 2
    # a scan of kappa-1 uses one constraint
    assert len(constraints_at(X, {"metal": 0, "donors": [1], "kind": "kappa1"}, 2.0)) == 1
    # --- ligand bonded in the input but gone after relaxation is reported
    assert relaxed_away(recs2, [r for r in recs2 if r["kind"] != "eta2"]) == [by_formula["CCHHHH"]]
    # --- energy window
    ins, rel = select_window({"a": -10.0, "b": 5.0, "c": 10.0, "d": 25.0}, 20.0)
    assert sorted(ins) == ["a", "b", "c"] and abs(rel["d"] - 35.0) < 1e-9, (ins, rel)
    ins, _ = select_window({"a": 0.0, "b": 20.0, "c": 20.01}, 20.0); assert sorted(ins) == ["a", "b"]
    print("dissociation_scan selftest PASS"); return 0


def main():
    ap = argparse.ArgumentParser(description="M1 extension: spin crossings on ligand-loss channels")
    ap.add_argument("xyz", nargs="*", help="one minimum, or several (batch with an energy window)")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int); ap.add_argument("--window", type=float, default=30.0)
    ap.add_argument("--dr", type=float, default=0.1); ap.add_argument("--stretch", type=float, default=2.5)
    ap.add_argument("--outdir", default="dscan"); ap.add_argument("--out", default="mecp_candidates_dissoc.json")
    ap.add_argument("--xtb", default="xtb"); ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--include-x", action="store_true",
                    help="also scan odd-electron (X-type) fragments: homolytic M-X cleavage")
    ap.add_argument("--emax", type=float, default=20.0,
                    help="batch: scan only minima within this energy (kcal/mol) of the lowest supplied")
    ap.add_argument("--energies", help="batch: file with `name energy_kcal` per line (default: GFN2 energies)")
    ap.add_argument("--jobs", type=int, default=1, help="batch: parallel minima")
    a = ap.parse_args()
    if a.selftest: return _selftest()
    if not a.xyz or a.mult is None:
        ap.error("provide at least one minimum geometry (xyz) and --mult (its multiplicity)")
    if len(a.xyz) == 1 and not a.energies:
        run(a.xyz[0], a.charge, a.mult, a.window, a.dr, a.stretch, a.outdir, a.out, a.xtb, include_x=a.include_x)
    else:
        run_batch(a.xyz, a.charge, a.mult, read_energies(a.energies) if a.energies else None, a.emax,
                  a.window, a.dr, a.stretch, a.outdir, a.out, a.xtb, a.include_x, a.jobs)


if __name__ == "__main__":
    sys.exit(main() or 0)
