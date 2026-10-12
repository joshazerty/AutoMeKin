# Fe(CO)5 CO loss: MECP on a ligand-loss channel

Input: `feco5_s.xyz`, the GFN2-xTB singlet Fe(CO)5 minimum. `run.sh` runs the three steps.

Reference results (UB3LYP/def2-SVP; Gaussian 16 A.03, xtb 6.7.1, easyMECP 0.3.2):

| Step | Axial CO (Fe1–C2) | Equatorial CO (Fe1–C6) |
|---|---|---|
| 1. GFN2 scan, min. gap 1↔3 (no GFN2 crossing) | 9.01 kcal/mol at 3.77 Å | 23.94 kcal/mol at 4.30 Å |
| 2. DFT crossing on the scan / seed frame | 2.90 Å / 2.87 Å (gap 0.76) | none / 2.80 Å (gap 47.53, min. gap) |
| 3. M2 | 74 steps, gap 0.006 kcal/mol | 128 steps, gap 0.001 kcal/mol |

Both seeds converge to the same Cs MECP between ¹Fe(CO)5 and ³Fe(CO)4 + CO
(Fe–C 2.31 Å, Fe–C–O 140°), +1.38 kcal/mol above ³Fe(CO)4 + CO; the other two GFN2
channels are symmetry-equivalent to the axial/equatorial ones and are merged in step 1.
At B3PW91* (15 % exact exchange)/def2-TZVP the axial seed gives +1.57 kcal/mol, and the
equatorial seed converges instead to the ³/¹ crossing within Fe(CO)4 (+10.2 kcal/mol).
Literature: Cs MECP 0.43 kcal/mol above ³Fe(CO)4 + CO (J. N. Harvey and M. Aschi,
Faraday Discuss. 124 (2003): 129, 10.1039/b211871h). Energies are electronic, without ZPE.
