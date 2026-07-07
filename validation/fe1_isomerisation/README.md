# Validation case 2 — Fe(I) alkene isomerisation (two-state σ-base resistance)

Reproduces the two-state mechanism of a Fe(I) bis(carbene)borate alkene-isomerisation
catalyst (Lutz, Hickey, Gao, Chen & Smith, *J. Am. Chem. Soc.* **2020**, *142*, 15527;
DOI 10.1021/jacs.0c07300): substrate binding on the quartet (S=3/2) surface,
oxidative addition to an η¹-allyl only on the doublet (S=1/2) surface.

## Reference energetics (OPBE/def2-TZVP // def2-SVP, ΔG₂₉₈ vs ⁴ π-complex, kcal/mol)

| species | ΔG |
|---|---|
| ⁴ OA-TS | +13.7 |
| ² OA-TS | +17.2 |
| allyl / Fe–H (² product) | +11.2 |

The two OA-TSs are near-degenerate (ΔΔG‡ = +3.5 kcal/mol, quartet lower), robust to
dispersion (ΔΔG‡ = +3.4 with D3(BJ)). **Functional sensitivity:** UB3LYP hides the
crossing; OPBE/TPSSh reveal it — hence OPBE is adopted as the reference functional.

## ⁴/² crossing location (balanced multireference)
SA-CASSCF(7,5)+SC-NEVPT2/def2-SVP on the trusted DFT geometries gives the vertical
doublet−quartet gap +22.1/+50.0 → +15.4/+43.2 → −29.9/−13.9 kcal (CASSCF/NEVPT2) at
the ⁴reactant, ⁴OA-TS and ²OA-TS. The gap changes sign in the OA-TS region at both
levels: the ⁴/² crossing is **co-located with oxidative addition**. (The minimal
Fe-3d active space gives reliable spin gaps / crossing location, not absolute
barriers — those come from DFT.)

## Walkthrough
1. **M0** — `python -m spin_aware.spin_enumerate` on the π-complex confirms
   single-reference validity (⟨S²⟩ clean; PASS).
2. **M1–M3** — crossing detection along the OA coordinate, then graph insertion.
3. **M4** — NA-TST / microkinetics with the quartet and doublet OA-TS branches.

SOC (ORCA 6.1, CASSCF(7,5)/QDPT-SOMF): H_SO ≈ 90 cm⁻¹ (RMS).
