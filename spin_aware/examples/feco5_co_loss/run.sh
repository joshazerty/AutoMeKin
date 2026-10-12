#!/bin/bash
# 1Fe(CO)5 -> 3Fe(CO)4 + CO: spin-forbidden CO loss / recombination MECP.
# Needs xtb (step 1), g16 (step 2) and easyMECP + g16 (step 3). Run from this directory
# with the repository root on PYTHONPATH.
set -e
python -m spin_aware.dissociation_scan feco5_s.xyz --charge 0 --mult 1 --window 30 \
       --outdir dscan --out mecp_candidates_dissoc.json
python -m spin_aware.dissociation_seed_dft --from-candidates mecp_candidates_dissoc.json \
       --route "UB3LYP/def2SVP" --workdir dft_seed --nproc 12 --mem 12GB \
       --out mecp_candidates_dissoc_dft.json
python -m spin_aware.locate_mecp --from-candidates mecp_candidates_dissoc_dft.json \
       --method UB3LYP --basis def2SVP --extra-route "Int=UltraFine SCF=(XQC,MaxCycle=300)" \
       --nproc 12 --mem 12GB --max-steps 150 --easymecp easymecp --gaussian-exe g16 --workdir m2
