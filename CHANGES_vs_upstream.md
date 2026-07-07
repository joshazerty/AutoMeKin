# Changes vs upstream AutoMeKin

This fork of [emartineznunez/AutoMeKin](https://github.com/emartineznunez/AutoMeKin)
adds a **spin-aware layer** as a new, self-contained `spin_aware/` package. The
upstream exploration core is unmodified; all additions are new files.

## New files
- `spin_aware/amk_xtb.py` — spin-capable GFN2-xTB engine (multiplicity → `--uhf`).
- `spin_aware/spin_enumerate.py` — M0 ⟨S²⟩ self-gate + M1 crossing detection.
- `spin_aware/locate_mecp.py` — M2 MECP localisation (easyMECP/Gaussian).
- `spin_aware/mecp_graph.py` — M3 MECP graph-node insertion.
- `spin_aware/na_tst.py` — M4 nonadiabatic TST rate core.
- `spin_aware/microkinetics.py` — M4 route-aware microkinetics.
- `spin_aware/annotate_rates.py` — M4 rate annotation onto the M3 graph.
- `environment.yml` — reproducible conda environment.
- `README_spin_aware.md` — pipeline overview and validation cases.

## Motivation
Stock AutoMeKin's MOPAC/PM7 low level applies multiplicity only at the Gaussian
high-level stage; every PM7 input runs closed-shell RHF, so a two-state search is
blind at the exploration stage. The `amk_xtb` engine plus the M0–M4 stages make the
spin degree of freedom explicit end-to-end, from sampling to nonadiabatic rate.
