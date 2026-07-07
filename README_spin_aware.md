# spin_aware — spin-aware extensions to AutoMeKin

This directory adds a **spin-aware layer** to [AutoMeKin](https://github.com/emartineznunez/AutoMeKin)
(MIT-licensed) for the automated discovery and rate-quantification of **two-state
(spin-crossing) reactivity** in transition-metal chemistry. It is the software
accompanying:

> J. Sims and J. N. Harvey, *Spin-Aware Automated Reaction-Network Exploration with
> Automatic Minimum-Energy-Crossing-Point Localisation and a Multireference Self-Gate*
> (submitted to *J. Chem. Theory Comput.*, 2026).

## The pipeline (M0–M4)

| Stage | Module | What it does |
|---|---|---|
| explorer | `amk_xtb.py` | Spin-capable GFN2-xTB engine — maps multiplicity to the xtb `--uhf` flag so each spin state is sampled explicitly (stock AutoMeKin's PM7 low level ignores multiplicity). |
| **M0** | `spin_enumerate.py` | ⟨S²⟩ **self-gate**: checks single-reference validity before committing DFT, and refuses multireference cases rather than returning a misleading single-determinant answer. |
| **M1** | `spin_enumerate.py` | Crossing detection: flags nodes where alternative-multiplicity states lie close in energy. |
| **M2** | `locate_mecp.py` | MECP localisation via the Harvey MECP algorithm (easyMECP + Gaussian). |
| **M3** | `mecp_graph.py` | Inserts the located MECP as a first-class node in the reaction-network graph (Kabsch-RMSD matching; DOT/Mermaid/Markdown export). |
| **M4** | `na_tst.py`, `microkinetics.py`, `annotate_rates.py` | Nonadiabatic (Landau–Zener/WKB) TST hop rate, route-aware microkinetics, and rate annotation onto the graph. |

Each module runs standalone (`python -m spin_aware.<module> --help`) and carries a
`_selftest`.

## Installation

```bash
conda env create -f ../environment.yml
conda activate automekin-spin-aware
```

External QM engines are **not** conda-installable and must be provided separately:
Gaussian 16 rev. A.03 (DFT + MECP single points), ORCA 6.1 (CASSCF/QDPT-SOMF
spin–orbit coupling), and optionally pyscf 2.13.0 (SA-CASSCF+NEVPT2 crossing map).

## Worked validation cases

See `validation/` for two end-to-end walkthroughs:
- **`validation/tpfe_bh_elim/`** — Tp-Fe(II)-ethyl β-H elimination (re-discovers the
  published quintet/triplet MECP; the ⁵/³ crossing undercuts the same-spin barrier).
- **`validation/fe1_isomerisation/`** — catalytic Fe(I) alkene isomerisation
  (two-state σ-base resistance; the ⁴/² crossing is co-located with oxidative addition).

## Relationship to upstream AutoMeKin

This is a fork of `emartineznunez/AutoMeKin` (MIT). The spin-aware layer is additive
— it lives entirely in `spin_aware/` and does not modify the upstream exploration
core. See `CHANGES_vs_upstream.md`. Licence: MIT (inherited).
