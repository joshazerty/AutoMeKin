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
| **M1** (ligand-loss channels) | `dissociation_scan.py`, `dissociation_seed_dft.py` | Crossings with no saddle point: GFN2-xTB relaxed scan of each monodentate metal–ligand bond of a minimum on the base multiplicity and on mult ± 2, then DFT single points on the scan frames; the frame nearest the DFT crossing is passed to M2. |
| **M2** | `locate_mecp.py` | MECP localisation via the Harvey MECP algorithm (easyMECP + Gaussian). |
| **M3** | `mecp_graph.py` | Inserts the located MECP as a first-class node in the reaction-network graph (Kabsch-RMSD matching; DOT/Mermaid/Markdown export). |
| **M4** | `na_tst.py`, `microkinetics.py`, `annotate_rates.py` | Nonadiabatic (Landau–Zener/WKB) TST hop rate, route-aware microkinetics, and rate annotation onto the graph. |

### Crossings on ligand-loss channels

`spin_enumerate.py` evaluates vertical gaps only at the minima and saddle points it is
given. A crossing on a barrierless dissociation channel (for example
¹Fe(CO)₅ → ³Fe(CO)₄ + CO) has no saddle point and is not seen by that check. For
each minimum, run:

```bash
# 1. GFN2 relaxed M–L scans per spin (r0 -> r0 + 2.5 Å, 0.1 Å steps; window 30 kcal/mol)
python -m spin_aware.dissociation_scan min.xyz --charge 0 --mult 1 --outdir dscan \
       --out mecp_candidates_dissoc.json
# 2. DFT single points on each surface's scan frames; seed at the DFT crossing
python -m spin_aware.dissociation_seed_dft --from-candidates mecp_candidates_dissoc.json \
       --route "UB3LYP/def2SVP" --workdir dft_seed --out mecp_candidates_dissoc_dft.json \
       --eligible-out mecp_candidates_dissoc_m2.json
# 3. M2 on the seeds whose DFT curves cross or whose smallest DFT gap is <= 30 kcal/mol
python -m spin_aware.locate_mecp --from-candidates mecp_candidates_dissoc_m2.json \
       --method UB3LYP --basis def2SVP --extra-route "Int=UltraFine SCF=(XQC,MaxCycle=300)" \
       --nproc 12 --mem 12GB --max-steps 300 --workdir m2
```

Step 2 is needed whenever GFN2 can misorder the spin states of the unsaturated
fragment. For Fe(CO)₄, GFN2 puts the singlet below the triplet, so the GFN2 curves never
cross while the DFT curves do. Pass the step-2 output to M2, not the step-1 GFN2 seeds.

Each M2 run writes `mecp_record.json` with `status` = `converged` (all five easyMECP
criteria met), `located` (the search stopped with final |ΔE| ≤ `--located-tol`, default
0.5 kcal/mol: a point on the crossing seam, whose energy is an upper bound to the MECP
energy) or `not_converged`, together with the final criteria, the step count and any
easyMECP error message.
Worked case: `examples/feco5_co_loss/`.

Each module runs standalone (`python -m spin_aware.<module> --help`) and carries a
`_selftest` (`for m in spin_enumerate dissociation_scan dissociation_seed_dft locate_mecp mecp_graph na_tst microkinetics annotate_rates; do python -m spin_aware.$m --selftest; done`).

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
