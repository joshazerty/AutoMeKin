# Changelog

All notable changes to the spin-aware AutoMeKin layer. Versions follow the release tags
in <https://github.com/joshazerty/AutoMeKin>.

## v1.0.0 — first public release

Tag `v1.0.0` on branch `spin-aware-release`. The spin-aware layer was developed through
commit `afc8714e6f6841fa77e3a9abb195204867b44207`; the tagged commit adds `m0_gate.py`
and the release metadata (`NOTICE`, `CHANGELOG.md`, `.zenodo.json`, and `CITATION.cff`
updated to this version). The release accompanying the methods paper; archived with its
supporting data at Zenodo.

### Added — spin-aware layer (`spin_aware/`)

The upstream AutoMeKin exploration core is unmodified; every addition is a new file.
File-by-file detail is in `code/CHANGES_vs_upstream.md`.

- `amk_xtb.py` — spin-capable GFN2-xTB low-level engine. Stock AutoMeKin applies
  multiplicity only at the Gaussian high-level stage, so its MOPAC/PM7 exploration runs
  closed-shell RHF and a two-state search is blind while sampling. This calculator maps
  multiplicity onto the xtb `--uhf` argument, making the spin degree of freedom explicit
  from the sampling stage onward.
- `spin_enumerate.py` — M0 multireference diagnostic gate and M1 detection of close-lying
  spin states along the explored coordinates.
- `locate_mecp.py` — M2 automatic minimum-energy-crossing-point localisation at flagged
  crossings, driving easyMECP over Gaussian.
- `mecp_graph.py` — M3 insertion of MECP nodes and non-adiabatic edges into the reaction
  graph.
- `na_tst.py`, `microkinetics.py`, `annotate_rates.py` — M4 nonadiabatic
  transition-state-theory rate core, route-aware microkinetics, and annotation of the
  resulting rates onto the M3 graph.
- `m0_gate.py` — standalone implementation of the M0 pass/flag rule, so a verdict can be
  regenerated from unrestricted-DFT output without running the pipeline.
- `environment.yml` — conda environment for the layer.
- `README_spin_aware.md`, `validation/tpfe_bh_elim/README.md`,
  `validation/fe1_isomerisation/README.md` — pipeline overview and the two validation
  walkthroughs.
