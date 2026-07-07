#!/usr/bin/env python3
"""
locate_mecp.py — Stage M2 of the MECP-aware-node feature (Project 18, Phase 2 Weeks 2+).

PURPOSE
-------
Take a spin-crossing candidate flagged by Stage M1 (`spin_enumerate.py`) and **localise the
minimum-energy crossing point (MECP)** between the two spin surfaces at DFT, by driving
**easyMECP** (Jaime Rodríguez-Guerra's self-contained wrapper around J. N. Harvey's MECP
Fortran, Gaussian-backed). DFT is the level that decides — GFN2 only flagged the candidate
(ESC-009); the MECP energy/geometry reported here is the trustworthy result.

WHY easyMECP
------------
ORCA 5 on `dirac` has no native MECP and ORCA 6 had a shared-lib issue (Phase1_Decision);
Gaussian 16 is the validated high-level. easyMECP is a single file, needs only g16 + gfortran,
and uses the original Harvey MECP algorithm — the Harvey-group lineage this project sits in.
(MECPro is the Python-native fallback if easyMECP's gfortran build is unavailable.)

WHAT IT DOES
------------
1. Generates easyMECP's single-file Gaussian input (`system.gjf`) from a geometry + the two
   multiplicities, using the `{A,B}` brace syntax for the divergent values (chk, title, mult).
2. Runs `python easymecp.py -f system.gjf --gaussian_exe g16` in a per-candidate workdir.
3. Collects the converged MECP geometry (easyMECP's `geom`) and the two surface energies at
   the MECP (from the final Gaussian logs), reporting their gap (|E_A - E_B|, the MECP
   convergence indicator) and the crossing energy.
4. Emits `mecp_record.json` (+ `mecp.xyz`). Batch mode consumes M1's `mecp_candidates.json`.

VALIDATION TARGET
-----------------
The Tp-Fe(II)-ethyl quintet(S=2)/triplet(S=1) crossing on the β-H path (V3). Seed from the
flagged geometry; easyMECP localises the ⁵/³ MECP; compare to Harvey 2013's spin-crossover.
Re-discovering a published MECP automatically is the project's proof-of-concept (aim #2).

Requires: easymecp.py (or `easymecp` on PATH), g16, gfortran. Input generation + log parsing
are unit-tested without Gaussian.
"""
import os
import re
import sys
import json
import glob
import argparse
import subprocess

PERIODIC = {  # minimal Z->symbol for writing geometries (extend as needed)
    1: "H", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S", 17: "Cl",
    26: "Fe", 27: "Co", 28: "Ni",
}
HARTREE_KCAL = 627.5094740631


# --------------------------------------------------------------------------- #
# input generation (easyMECP single-file {A,B} format)
# --------------------------------------------------------------------------- #
def read_xyz(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    rows = []
    for ln in lines[2:2 + n]:
        p = ln.split()
        rows.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    return rows


def generate_gjf(geom_rows, charge, mult_a, mult_b, method="UB3LYP", basis="def2SVP",
                 mem="16GB", nproc=8, max_steps=80, extra_route="", out="system.gjf"):
    """Write the easyMECP single-file input. Divergent values use {A,B} braces."""
    route = f"#n {method}/{basis} force guess(read)"
    if extra_route:
        route += " " + extra_route.strip()
    lines = [
        f"! easymecp: max_steps={max_steps}",
        f"%mem={mem}",
        f"%nproc={nproc}",
        "%chk={A,B}.chk",
        route,
        "",
        "{StateA,StateB} MECP",
        "",
        f"{charge} {{{mult_a},{mult_b}}}",
    ]
    for sym, x, y, z in geom_rows:
        lines.append(f" {sym:3s} {x: .8f} {y: .8f} {z: .8f}")
    lines.append("")            # trailing blank required by Gaussian
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return out


# --------------------------------------------------------------------------- #
# result collection
# --------------------------------------------------------------------------- #
def parse_scf_done(logpath):
    """Last 'SCF Done: E(...) = <Eh>' from a Gaussian log."""
    e = None
    with open(logpath, errors="ignore") as fh:
        for ln in fh:
            m = re.search(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)", ln)
            if m:
                e = float(m.group(1))
    return e


def read_easymecp_geom(geom_path):
    """easyMECP writes the latest geometry to `geom` (element symbol OR Z, then x y z)."""
    rows = []
    with open(geom_path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 4:
                continue
            tok = p[0]
            sym = PERIODIC.get(int(tok), tok) if tok.isdigit() else tok
            try:
                rows.append((sym, float(p[1]), float(p[2]), float(p[3])))
            except ValueError:
                continue
    return rows


def write_xyz(rows, path, comment=""):
    with open(path, "w") as fh:
        fh.write(f"{len(rows)}\n{comment}\n")
        for sym, x, y, z in rows:
            fh.write(f"{sym:3s} {x: .8f} {y: .8f} {z: .8f}\n")


def collect_result(workdir):
    """Gather converged MECP geometry + the two surface energies at the MECP.

    Energies: prefer the two final-step Gaussian logs (A.log/B.log if present, else the two
    most-recently-modified *.log in workdir). Gap |E_A - E_B| is the convergence indicator
    (small at a real MECP). Crossing energy = mean of the two."""
    res = {"workdir": os.path.abspath(workdir), "converged": None,
           "e_surface_a": None, "e_surface_b": None, "gap_kcal": None, "crossing_e_eh": None,
           "n_geom_atoms": None, "report_tail": None}
    # geometry
    geom = os.path.join(workdir, "geom")
    if os.path.isfile(geom):
        rows = read_easymecp_geom(geom)
        if rows:
            write_xyz(rows, os.path.join(workdir, "mecp.xyz"), "MECP (easyMECP converged)")
            res["n_geom_atoms"] = len(rows)
    # energies: A.log / B.log, else two newest logs
    logs = []
    for tag in ("A", "B"):
        for cand in (f"{tag}.log", os.path.join("JOBS", f"{tag}.log")):
            p = os.path.join(workdir, cand)
            if os.path.isfile(p):
                logs.append(p); break
    if len(logs) < 2:
        logs = sorted(glob.glob(os.path.join(workdir, "*.log")),
                      key=os.path.getmtime)[-2:]
    energies = [parse_scf_done(p) for p in logs] if len(logs) >= 2 else []
    energies = [e for e in energies if e is not None]
    if len(energies) == 2:
        res["e_surface_a"], res["e_surface_b"] = energies
        res["gap_kcal"] = round(abs(energies[0] - energies[1]) * HARTREE_KCAL, 3)
        res["crossing_e_eh"] = round(sum(energies) / 2, 6)
    # report tail (for the agent to confirm convergence wording on first run)
    rep = os.path.join(workdir, "ReportFile")
    if os.path.isfile(rep):
        res["report_tail"] = "".join(open(rep, errors="ignore").readlines()[-15:])
        low = res["report_tail"].lower()
        res["converged"] = ("converged" in low) or ("mecp" in low and "found" in low)
    return res


# --------------------------------------------------------------------------- #
# run flow
# --------------------------------------------------------------------------- #
def locate_one(geom_xyz, charge, mult_a, mult_b, workdir, easymecp="easymecp.py",
               gaussian_exe="g16", method="UB3LYP", basis="def2SVP", max_steps=80,
               mem="16GB", nproc=8, extra_route="", dry_run=False):
    os.makedirs(workdir, exist_ok=True)
    rows = read_xyz(geom_xyz)
    gjf = generate_gjf(rows, charge, mult_a, mult_b, method, basis, mem, nproc, max_steps,
                       extra_route, out=os.path.join(workdir, "system.gjf"))
    cmd = ([sys.executable, easymecp] if easymecp.endswith(".py") else [easymecp]) + \
          ["-f", os.path.basename(gjf), "--gaussian_exe", gaussian_exe]
    record = {"geom_in": os.path.abspath(geom_xyz), "charge": charge,
              "mult_a": mult_a, "mult_b": mult_b, "method": method, "basis": basis,
              "cmd": " ".join(cmd), "workdir": os.path.abspath(workdir)}
    if dry_run:
        record["dry_run"] = True
        return record
    with open(os.path.join(workdir, "easymecp.out"), "w") as out:
        p = subprocess.run(cmd, cwd=workdir, stdout=out, stderr=subprocess.STDOUT)
    record["returncode"] = p.returncode
    record.update(collect_result(workdir))
    with open(os.path.join(workdir, "mecp_record.json"), "w") as fh:
        json.dump(record, fh, indent=2)
    return record


def main():
    ap = argparse.ArgumentParser(description="M2: localise an MECP between two spin states (easyMECP/DFT)")
    sub = ap.add_argument
    sub("--geom", help="seed geometry (.xyz)")
    sub("--charge", type=int, default=0)
    sub("--mult-a", type=int, help="multiplicity A (e.g. 5 quintet)")
    sub("--mult-b", type=int, help="multiplicity B (e.g. 3 triplet)")
    sub("--from-candidates", help="mecp_candidates.json from spin_enumerate.py (batch mode)")
    sub("--workdir", default="mecp_run")
    sub("--easymecp", default="easymecp.py")
    sub("--gaussian-exe", default="g16")
    sub("--method", default="UB3LYP")
    sub("--basis", default="def2SVP")
    sub("--max-steps", type=int, default=80)
    sub("--dry-run", action="store_true", help="generate inputs only; do not run Gaussian")
    sub("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return _selftest()

    runs = []
    if args.from_candidates:
        cands = json.load(open(args.from_candidates))["candidates"]
        for i, c in enumerate(cands):
            runs.append(dict(geom_xyz=c["xyz"], charge=c["charge"],
                             mult_a=c["base_mult"], mult_b=c["alt_mult"],
                             workdir=os.path.join(args.workdir, f"cand{i}_{c['node']}")))
    else:
        if not (args.geom and args.mult_a and args.mult_b):
            ap.error("provide --geom --mult-a --mult-b, or --from-candidates")
        runs.append(dict(geom_xyz=args.geom, charge=args.charge,
                         mult_a=args.mult_a, mult_b=args.mult_b, workdir=args.workdir))
    for r in runs:
        rec = locate_one(easymecp=args.easymecp, gaussian_exe=args.gaussian_exe,
                         method=args.method, basis=args.basis, max_steps=args.max_steps,
                         dry_run=args.dry_run, **r)
        print(json.dumps({k: rec.get(k) for k in
                          ("workdir", "returncode", "converged", "gap_kcal",
                           "crossing_e_eh", "n_geom_atoms")}, indent=2))


def _selftest():
    import tempfile
    # 1) input generation
    rows = [("Fe", 0.0, 0.0, 0.0), ("C", 0.0, 0.0, 1.8), ("O", 0.0, 0.0, 2.95)]
    with tempfile.TemporaryDirectory() as d:
        gjf = generate_gjf(rows, 0, 5, 3, out=os.path.join(d, "system.gjf"))
        txt = open(gjf).read()
        assert "%chk={A,B}.chk" in txt and "0 {5,3}" in txt and "{StateA,StateB} MECP" in txt
        assert "UB3LYP/def2SVP" in txt and txt.rstrip().endswith("2.95000000")
        # 2) SCF parser
        log = os.path.join(d, "A.log")
        open(log, "w").write("blah\n SCF Done:  E(UB3LYP) =  -2044.686398 A.U. after 12 cycles\n end\n")
        assert abs(parse_scf_done(log) - (-2044.686398)) < 1e-6
        # 3) easyMECP geom reader (atomic numbers -> symbols)
        g = os.path.join(d, "geom")
        open(g, "w").write("26 0.0 0.0 0.0\n6 0.0 0.0 1.8\n8 0.0 0.0 2.95\n")
        gr = read_easymecp_geom(g)
        assert gr[0][0] == "Fe" and gr[2][0] == "O" and len(gr) == 3
    print("M2 selftest PASS: gjf generation, SCF parse, geom(Z->sym) read all OK")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
