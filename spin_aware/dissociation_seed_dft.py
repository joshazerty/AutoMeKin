#!/usr/bin/env python3
"""
dissociation_seed_dft.py - M1 (ligand-loss channels), step 2: DFT refinement of the seed.

GFN2 can misorder the spin states of the unsaturated fragment (for Fe(CO)4 it
puts the singlet below the triplet, the opposite of CCSD(T)), so the GFN2
relaxed curves need not cross even when the DFT curves do. This step takes the
GFN2 relaxed-scan frames written by dissociation_scan.py, runs DFT single
points on each surface's own frames (Gaussian), and locates the crossing of
the two DFT curves E_base(r) and E_alt(r). The seed passed to M2 is the
base-surface frame nearest that crossing; if the DFT curves do not cross inside
the scan, the frame with the smallest |gap| is used. Output schema = M1
candidates, so `locate_mecp.py --from-candidates` consumes it unchanged.

Usage
-----
Batch (all channels flagged by dissociation_scan.py):
  python dissociation_seed_dft.py --from-candidates mecp_candidates_dissoc.json \
      --route "UB3LYP/def2SVP" --workdir dft_seed
Single channel (scan files <scan-dir>/<tag>_mult<m>_scan.xyz):
  python dissociation_seed_dft.py --scan-dir dscan --tag feco5_s_LC2 \
      --base-mult 1 --alt-mult 3 --route "UB3LYP/def2SVP" --workdir dft_seed

Single points already present (Normal termination + SCF Done) are reused, so a
re-run only computes missing frames. --reuse-only never launches Gaussian and
fails if a frame is missing (regression checks against archived runs).

Requires: Gaussian 16 (g16) for new single points.
"""
import os, re, sys, json, argparse, subprocess
import numpy as np

HARTREE_KCAL = 627.5094740631


def read_frames(path):
    """Multi-frame xyz written by dissociation_scan.py -> [(r_A, rows)]."""
    L = open(path).read().splitlines(); out = []; i = 0
    while i < len(L):
        if not L[i].strip():
            i += 1; continue
        n = int(L[i]); com = L[i+1]
        r = float(re.search(r"r=([\d.]+)", com).group(1))
        rows = [ln.split() for ln in L[i+2:i+2+n]]
        out.append((r, [(p[0], float(p[1]), float(p[2]), float(p[3])) for p in rows])); i += n+2
    return out


def gjf(rows, charge, mult, route, nproc, mem, title, chk=None, guess_read=False):
    """Gaussian input. chk: checkpoint name written by this single point; guess_read: start the SCF from
    the orbitals in that checkpoint (left by the previous frame of the same surface)."""
    s = [f"%nprocshared={nproc}", f"%mem={mem}"] + ([f"%chk={chk}"] if chk else []) + \
        [f"#p {route} SCF=(XQC,MaxCycle=300) Int=UltraFine" + (" Guess=Read" if guess_read else ""),
         "", title, "", f"{charge} {mult}"]
    s += [f" {a:3s} {x: .8f} {y: .8f} {z: .8f}" for a, x, y, z in rows]
    return "\n".join(s) + "\n\n"


def scf(log):
    """Last SCF energy of a normally terminated Gaussian log, else None."""
    e = None; ok = False
    for ln in open(log, errors="ignore"):
        m = re.search(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)", ln)
        if m: e = float(m.group(1))
        if "Normal termination" in ln: ok = True
    return e if ok else None


def crossing_and_seed(R, g):
    """First sign change of the DFT gap g(r) (linear interpolation) and the index of
    the frame nearest it; without a sign change, the frame with min |g|."""
    cross, k = None, int(np.nanargmin(np.abs(g)))
    for i in range(len(g)-1):
        if np.isfinite(g[i]) and np.isfinite(g[i+1]) and g[i]*g[i+1] < 0:
            t = g[i]/(g[i]-g[i+1]); cross = float(R[i]+t*(R[i+1]-R[i]))
            k = i if abs(g[i]) <= abs(g[i+1]) else i+1; break
    return cross, k


def run(scan_dir, tag, charge, base_mult, alt_mult, route, workdir, g16="g16",
        nproc=8, mem="8GB", every=1, out_json="mecp_candidates_dissoc_dft.json", dry=False,
        scan_files=None, reuse_only=False, write_json=True, m2_window=30.0, guess_read=True):
    """Refine one channel. scan_files ({mult: path}) overrides scan_dir/tag lookup."""
    os.makedirs(workdir, exist_ok=True)
    if scan_files is None:
        scan_files = {mu: os.path.join(scan_dir, f"{tag}_mult{mu}_scan.xyz")
                      for mu in (base_mult, alt_mult)}
    fr = {mu: read_frames(scan_files[mu])[::every] for mu in (base_mult, alt_mult)}
    E = {}
    for mu, frames in fr.items():
        E[mu] = []
        for k, (r, rows) in enumerate(frames):
            stem = os.path.join(workdir, f"sp_m{mu}_{k:02d}")
            if not (os.path.isfile(stem+".log") and scf(stem+".log") is not None):
                if reuse_only:
                    raise FileNotFoundError(f"--reuse-only: no completed log {stem}.log")
                chk = f"sp_m{mu}.chk"          # one checkpoint per surface, overwritten frame by frame
                title = f"{tag} r={r:.3f} mult {mu}"
                use_guess = guess_read and k > 0 and os.path.isfile(os.path.join(workdir, chk))
                open(stem+".gjf", "w").write(gjf(rows, charge, mu, route, nproc, mem, title,
                                                 chk=chk if guess_read else None, guess_read=use_guess))
                if not dry:
                    with open(stem+".gjf") as fi, open(stem+".log", "w") as fo:
                        subprocess.run([g16], stdin=fi, stdout=fo, stderr=subprocess.STDOUT, cwd=workdir)
                    if use_guess and scf(stem+".log") is None:   # guess from the neighbouring frame failed: plain guess
                        open(stem+".gjf", "w").write(gjf(rows, charge, mu, route, nproc, mem, title,
                                                         chk=chk, guess_read=False))
                        with open(stem+".gjf") as fi, open(stem+".log", "w") as fo:
                            subprocess.run([g16], stdin=fi, stdout=fo, stderr=subprocess.STDOUT, cwd=workdir)
            E[mu].append(None if dry else scf(stem+".log"))
    if dry:
        return None
    R = np.array([r for r, _ in fr[base_mult]])
    eb = np.array(E[base_mult], float); ea = np.array(E[alt_mult], float)
    g = (ea-eb)*HARTREE_KCAL
    cross, k = crossing_and_seed(R, g)
    seed = os.path.abspath(os.path.join(workdir, f"{tag}_dftseed.xyz"))
    rows = fr[base_mult][k][1]
    with open(seed, "w") as fh:
        fh.write(f"{len(rows)}\nDFT-refined dissociation seed r={R[k]:.3f} A gap={g[k]:.2f} kcal ({route})\n")
        for a, x, y, z in rows: fh.write(f"{a:3s} {x: .8f} {y: .8f} {z: .8f}\n")
    prof = {"r_A": R.tolist(), "gap_dft_kcal": np.round(g, 2).tolist(),
            "E_base_rel_kcal": np.round((eb-eb[0])*HARTREE_KCAL, 2).tolist(),
            "E_alt_rel_kcal": np.round((ea-eb[0])*HARTREE_KCAL, 2).tolist(),
            "crossing_r_A": cross, "seed_index": k, "route": route}
    json.dump(prof, open(os.path.join(workdir, "dft_profile.json"), "w"), indent=1)
    min_gap = float(np.nanmin(np.abs(g)))
    cand = {"node": f"{tag}_m{base_mult}to{alt_mult}_dft", "xyz": seed, "charge": charge,
            "base_mult": base_mult, "alt_mult": alt_mult, "dE_vert_kcal": round(float(g[k]), 2),
            "kind": "dissociation_channel", "seed_level": route, "r_seed_A": float(R[k]),
            "crossing_r_A": cross,
            "min_abs_gap_dft_kcal": round(min_gap, 2),
            "m2_eligible": bool(cross is not None or min_gap <= m2_window),
            "m2_window_kcal": m2_window,
            "note": "DFT curves cross on ligand-loss channel" if cross is not None
                    else "no DFT crossing inside scan; seed at min |gap|"}
    if write_json:
        json.dump({"stage": "M1-dissociation-dft", "candidates": [cand]}, open(out_json, "w"), indent=2)
    print(json.dumps(cand, indent=1)); print("gap profile (kcal):", np.round(g, 1).tolist())
    return cand


def run_batch(cand_json, route, workdir, g16="g16", nproc=8, mem="8GB", every=1,
              out_json="mecp_candidates_dissoc_dft.json", dry=False, reuse_only=False,
              m2_window=30.0, eligible_out=None, guess_read=True):
    """Refine every candidate of a dissociation_scan.py JSON (one subdir per channel)."""
    cands = json.load(open(cand_json))["candidates"]
    out = []
    for c in cands:
        if c.get("kind") != "dissociation_channel" or "scan_files" not in c:
            print(f"  [skip] {c.get('node')}: not a dissociation_scan.py candidate", file=sys.stderr)
            continue
        b, a = int(c["base_mult"]), int(c["alt_mult"])
        sf = {int(m): p for m, p in c["scan_files"].items()}
        tag = c.get("scan_tag") or c["node"].rsplit("_m", 1)[0]
        res = run(None, tag, c["charge"], b, a, route,
                  os.path.join(workdir, f"{tag}_m{b}to{a}"), g16, nproc, mem, every,
                  out_json, dry, scan_files=sf, reuse_only=reuse_only, write_json=False,
                  m2_window=m2_window, guess_read=guess_read)
        if res is not None:
            res["channel"] = c.get("channel")
            out.append(res)
    if not dry:
        json.dump({"stage": "M1-dissociation-dft", "route": route, "candidates": out},
                  open(out_json, "w"), indent=2)
        print(f"{len(out)} DFT-refined dissociation seed(s) -> {out_json}")
        if eligible_out:       # only seeds that pass the M2 rule: DFT curves cross or min |gap| <= m2_window
            el = [c for c in out if c["m2_eligible"]]
            json.dump({"stage": "M1-dissociation-dft", "route": route, "candidates": el},
                      open(eligible_out, "w"), indent=2)
            print(f"{len(el)} of {len(out)} seed(s) eligible for M2 (DFT crossing or min |gap| <= {m2_window} kcal/mol) -> {eligible_out}")
    return out


def _selftest():
    """Synthetic scan + fake Gaussian logs: crossing found, seed picked, reuse works."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        rows = [("Fe", 0.0, 0.0, 0.0), ("C", 0.0, 0.0, 1.8), ("O", 0.0, 0.0, 2.95)]
        R = [1.8, 2.0, 2.2, 2.4]
        gaps = [6.0, 2.0, -1.0, -4.0]                       # kcal; crossing between 2.0 and 2.2
        sf = {}
        for mu in (1, 3):
            sf[mu] = os.path.join(d, f"t_mult{mu}_scan.xyz")
            with open(sf[mu], "w") as fh:
                for r in R:
                    fh.write(f"3\nr={r:.4f} E_gfn2=0 mult={mu}\n")
                    for a, x, y, z in rows: fh.write(f"{a} {x} {y} {z+r-1.8}\n")
        wd = os.path.join(d, "w"); os.makedirs(wd)
        for k, gp in enumerate(gaps):
            for mu, e in ((1, -100.0), (3, -100.0 + gp/HARTREE_KCAL)):
                open(os.path.join(wd, f"sp_m{mu}_{k:02d}.log"), "w").write(
                    f" SCF Done:  E(UB3LYP) =  {e:.10f}     A.U. after 10 cycles\n Normal termination of Gaussian 16\n")
        cj = os.path.join(d, "cands.json")
        json.dump({"candidates": [{"node": "t_m1to3", "xyz": "x", "charge": 0, "base_mult": 1,
                                   "alt_mult": 3, "kind": "dissociation_channel", "scan_tag": "t",
                                   "scan_files": {"1": sf[1], "3": sf[3]}}]}, open(cj, "w"))
        # batch mode reads <workdir>/<tag>_m1to3/: move the fake logs there
        bw = os.path.join(d, "batch", "t_m1to3"); os.makedirs(bw)
        for f in os.listdir(wd): os.replace(os.path.join(wd, f), os.path.join(bw, f))
        out = run_batch(cj, "UB3LYP/def2SVP", os.path.join(d, "batch"),
                        out_json=os.path.join(d, "o.json"), reuse_only=True)
        c = out[0]
        assert c["crossing_r_A"] is not None and 2.0 < c["crossing_r_A"] < 2.2, c
        assert abs(c["r_seed_A"] - 2.2) < 1e-9 and c["dE_vert_kcal"] == -1.0, c   # |-1| < |2|
        assert json.load(open(os.path.join(d, "o.json")))["candidates"][0]["node"] == "t_m1to3_dft"
        assert c["m2_eligible"] is True and c["min_abs_gap_dft_kcal"] == 1.0, c
        # M2 rule: no crossing and min |gap| above the window -> not eligible; written file keeps only eligible seeds
        far = [60.0, 55.0, 52.0, 53.0]
        wd2 = os.path.join(d, "batch2", "t_m1to3"); os.makedirs(wd2)
        for k, gp in enumerate(far):
            for mu, e in ((1, -100.0), (3, -100.0 + gp/HARTREE_KCAL)):
                open(os.path.join(wd2, f"sp_m{mu}_{k:02d}.log"), "w").write(
                    f" SCF Done:  E(UB3LYP) =  {e:.10f}     A.U. after 10 cycles\n Normal termination of Gaussian 16\n")
        eo = os.path.join(d, "elig.json")
        out2 = run_batch(cj, "UB3LYP/def2SVP", os.path.join(d, "batch2"), out_json=os.path.join(d, "o2.json"),
                         reuse_only=True, eligible_out=eo)
        t = gjf([("Fe", 0.0, 0.0, 0.0)], 0, 3, "UB3LYP/def2SVP", 4, "4GB", "x", chk="sp_m3.chk", guess_read=True)
        assert "%chk=sp_m3.chk" in t and "Guess=Read" in t and "XQC" in t
        t0 = gjf([("Fe", 0.0, 0.0, 0.0)], 0, 3, "UB3LYP/def2SVP", 4, "4GB", "x"); assert "chk" not in t0 and "Guess" not in t0
        assert out2[0]["m2_eligible"] is False and out2[0]["crossing_r_A"] is None, out2[0]
        assert json.load(open(eo))["candidates"] == []
        # no sign change -> min |gap| frame, crossing None
        cr, k = crossing_and_seed(np.array(R), np.array([9.0, 5.0, 3.0, 4.0]))
        assert cr is None and k == 2
        # reuse-only refuses to launch Gaussian
        try:
            run(None, "t", 0, 1, 3, "x", os.path.join(d, "empty"), scan_files=sf, reuse_only=True)
            raise AssertionError("reuse_only did not raise")
        except FileNotFoundError:
            pass
    print("dissociation_seed_dft selftest PASS")
    return 0


def main():
    ap = argparse.ArgumentParser(description="M1 ligand-loss channels: DFT refinement of the scan seed")
    ap.add_argument("--from-candidates", help="mecp_candidates_dissoc.json from dissociation_scan.py (batch)")
    ap.add_argument("--scan-dir"); ap.add_argument("--tag")
    ap.add_argument("--charge", type=int, default=0); ap.add_argument("--base-mult", type=int)
    ap.add_argument("--alt-mult", type=int); ap.add_argument("--route")
    ap.add_argument("--workdir", default="dft_seed"); ap.add_argument("--g16", default="g16")
    ap.add_argument("--nproc", type=int, default=8); ap.add_argument("--mem", default="8GB")
    ap.add_argument("--every", type=int, default=1, help="use every n-th scan frame")
    ap.add_argument("--out", default="mecp_candidates_dissoc_dft.json")
    ap.add_argument("--dry-run", action="store_true", help="write .gjf inputs only")
    ap.add_argument("--reuse-only", action="store_true",
                    help="never run Gaussian; fail if a single point is missing")
    ap.add_argument("--m2-window", type=float, default=30.0,
                    help="seed is eligible for M2 if the DFT curves cross or min |gap| <= this (kcal/mol)")
    ap.add_argument("--eligible-out", help="batch: also write the M2-eligible seeds to this JSON")
    ap.add_argument("--no-guess-read", action="store_true",
                    help="start every single point from the default guess (default: reuse the orbitals of the "
                         "previous frame of the same surface, with a fallback to the default guess)")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return _selftest()
    if not a.route:
        ap.error("--route is required (e.g. 'UB3LYP/def2SVP')")
    if a.from_candidates:
        run_batch(a.from_candidates, a.route, a.workdir, a.g16, a.nproc, a.mem, a.every,
                  a.out, a.dry_run, a.reuse_only, a.m2_window, a.eligible_out, not a.no_guess_read)
    else:
        if not (a.scan_dir and a.tag and a.base_mult and a.alt_mult):
            ap.error("provide --from-candidates, or --scan-dir --tag --base-mult --alt-mult")
        run(a.scan_dir, a.tag, a.charge, a.base_mult, a.alt_mult, a.route, a.workdir, a.g16,
            a.nproc, a.mem, a.every, a.out, a.dry_run, reuse_only=a.reuse_only, m2_window=a.m2_window,
            guess_read=not a.no_guess_read)


if __name__ == "__main__":
    sys.exit(main() or 0)
