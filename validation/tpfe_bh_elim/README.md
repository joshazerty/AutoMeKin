# Validation case 1 — Tp-Fe(II)-ethyl β-hydride elimination

Re-discovers the published quintet(S=2)/triplet(S=1) spin crossing on the β-H
elimination path of a Tp-supported Fe(II)-ethyl complex (Bellows, Cundari & Holland,
*Organometallics* **2013**, *32*, 4741; DOI 10.1021/om400325x).

## Reference energetics (UB3LYP, relative to the ⁵ reactant, kcal/mol)

| basis | ³reactant | ³β-H TS | ⁵β-H TS | ⁵/³ MECP |
|---|---|---|---|---|
| def2-SVP        | 22.3 | 26.9 | 40.4 | 19.6 |
| def2-TZVP       | 21.9 | 27.1 | 40.3 | 17.1 |
| def2-TZVP+D3(BJ)| 18.3 | 25.4 | 36.7 | 13.4 |

The spin-cross route (via the ⁵/³ MECP) undercuts the same-spin quintet barrier by
~13.5 kcal/mol; the operative TS is the triplet β-H TS.

## Walkthrough
1. **M0/M1** — `python -m spin_aware.spin_enumerate --geom reactant.xyz --charge 0 --base-mult 5`
   flags the close-lying triplet along the β-H coordinate.
2. **M2** — `python -m spin_aware.locate_mecp --geom seed.xyz --charge 0 --mult-a 5 --mult-b 3 --method UB3LYP --basis def2TZVP`
   localises the MECP (converges to a surface gap ~1×10⁻³ kcal/mol).
3. **M3** — `python -m spin_aware.mecp_graph` inserts the MECP node.
4. **M4** — `python -m spin_aware.na_tst` / `microkinetics` returns κ ≈ 0.77,
   k(298 K) ≈ 1.5×10⁻³ s⁻¹ for the hop (~4 orders faster than turnover), and the
   spin-cross route as operative.

SOC (ORCA 6.1, CASSCF(6,5)/QDPT-SOMF): H_SO ≈ 279 cm⁻¹.
