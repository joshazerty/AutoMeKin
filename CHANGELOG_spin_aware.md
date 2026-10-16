## spin_aware 1.1.0 (tag v1.1.0, 10/16/2026; manuscript revision)

The package version now follows the repository tag. Includes everything listed under 0.2.0 below.

| Date | Changes |
| ---------- | ----------|
|10/16/2026| <code>locate_mecp.py</code>: default <code>--max-steps</code> raised from 80 to 300, the single step limit used for every MECP search of the revised paper. <code>mecp_record.json</code> now records <code>max_steps</code>, the step count, the final value/threshold/met flag of the five easyMECP criteria, <code>final_de_kcal</code>, any easyMECP error message (previously a step that failed with "MECP returned code 2" was reported as returncode 0) and a <code>status</code> label: <code>converged</code> (all five criteria met), <code>located</code> (final \|ΔE\| ≤ <code>--located-tol</code>, default 0.5 kcal/mol, the conservative MECP-energy error of Harvey and Aschi, Faraday Discuss. 2003, 124, 129) or <code>not_converged</code>. The convergence test now matches easyMECP's "has CONVERGED" message instead of any occurrence of "converged".|
|10/16/2026| <code>dissociation_scan.py</code>: X-type (homolysis) bonds are skipped by default (<code>--include-x</code> keeps them) and listed; η² ligands are scanned along the metal–midpoint distance; <code>--emax</code> restricts the scan to minima within a window of the lowest minimum (20 kcal/mol in the paper). <code>dissociation_seed_dft.py</code>: <code>Guess=Read</code> along each surface, and the M2 rule (DFT curves cross or smallest DFT gap ≤ 30 kcal/mol) written to <code>--eligible-out</code>. This is the code used for the 20 kcal/mol ligand-loss campaign.|
|10/16/2026| <code>mecp_graph.py</code>: reads native AutoMeKin RXNet files including PR product lines and optional MINinfo, de-duplicates MECP records (permutation-aware RMSD &lt; 0.3 Å and \|ΔE\| &lt; 0.5 kcal/mol), new options <code>--perm-match</code> and <code>--partner-role</code>. The Tp-Fe example graph is reproduced unchanged.|
|10/16/2026| <code>m0_gate.py</code>: CLI options <code>--pass-tol</code>, <code>--flag-tol</code> and <code>--min-flag</code> (defaults unchanged: 0.5, 1.0, 2), so the thresholds can be varied without editing the code. All verdicts reported in the paper are unchanged.|
|10/16/2026| <code>examples/feco5_co_loss/run.sh</code> and <code>README_spin_aware.md</code>: M2 runs on the <code>--eligible-out</code> seeds with 300 steps.|

## spin_aware 0.2.0 (merged into 1.1.0)

| Date | Changes |
| ---------- | ----------|
|10/12/2026| New M1 modules <code>dissociation_scan.py</code> (GFN2 relaxed metal–ligand scans per multiplicity; flags crossings on ligand-loss channels, which have no saddle point) and <code>dissociation_seed_dft.py</code> (DFT single points on the scan frames; seed at the DFT crossing; batch <code>--from-candidates</code>, <code>--reuse-only</code>). Output is M1 candidate schema, consumed unchanged by <code>locate_mecp.py</code>. Worked case <code>spin_aware/examples/feco5_co_loss/</code>.|
|10/12/2026| <code>locate_mecp.py</code>: CLI options <code>--nproc</code>, <code>--mem</code>, <code>--extra-route</code> (previously only in the Python API); result collection reads the final-step easyMECP 0.3.2 logs <code>JOBS/Job&lt;N&gt;_{A,B}.log</code> and falls back to the last frame of <code>JOBS/trajectory.xyz</code> when no <code>geom</code> file is written. This is the version used for the Tp-Fe def2-TZVP and Fe(CO)5 MECP runs.|
|10/12/2026| <code>na_tst.py</code>: the <code>trapz</code> fallback evaluated <code>np.trapz</code> eagerly and failed to import under NumPy 2; <code>na_tst</code> and <code>annotate_rates</code> now run on NumPy 1.x and 2.x (numerics unchanged).|

## Revision number 1142 (September 26, 2024) 

| Date | 	Changes |
| ----------	 | ----------|
|10/18/2024|	createMat.py, Heuristics.py, mopacamk.py and tag_prod.py have been updated to avoid futurewarning with networkx3.0|
|09/26/2024|	New keyword: <code>timeout</code>. This is used for jobs taking more than a certain time, and known to fail.|
|06/24/2024|	Fixed a bug in <code>IRC.sh</code> and added keyword <code>MaxBO</code> to set the maximum bond order for breaking a bond |
|06/05/2024|    <code>mopacamk</code> changed to make it compatible with ASE's newest version. FutureWarnings from locate_barrless to barrless.err|
| 05/23/2024|   Bugfixes in <code>amk_parallel.sh</code> and <code>sel_mol.sh</code>. locate_barrless modified to allow multiple searchs.| 
| 03/06/2024| 	Removed dependence with <code>zenity</code> and <code>gnuplot</code>|  
| 01/15/2024| 	Fixed some bugs with <code>llcalcs.sh</code> and improved ChemKnow algorithm (new keyword <code>CK_minima</code>)| 
| 11/07/2023| 	Fixed some issues with the HL calcs of fragments| 
| 05/12/2022| 	Fixed some issues with <code>qcore</code> calcs| 
| 04/26/2022| 	Bugfix in <code>utils.sh</code> and improved performace hl calcs| 
| 04/13/2022| 	Removed STOP in <code>diag.f</code>.| 
| 03/17/2022| 	Bugfix in <code>final.sh</code>.| 
| 03/11/2022| 	Interface with <code>Gaussian16</code>.| 
| 02/03/2022| 	Min and Max temperatures for the kinetics set to 100 K and 9999 K, respectively.| 
| 12/10/2021| 	Bugfixes in Python scripts that read inputfile (charge and long mopac inputs were not read correctly).| 
| 12/01/2021| 	Bugfix in <code>LocateTS.py</code>.| 
| 11/18/2021| 	Bugfix in ChemKnow. Improved torsional search| 
| 11/16/2021| 	ChemKnow and barrierless search improved. | 
| 10/17/2021| 	Updated tables with the geometry orientation for normal mode calcls in <code>FINALDIR</code> (low-level one). | 
| 09/17/2021| 	<code>select.sh</code> now exectuted from WRKDIR| 
| 07/24/2021| 	Bugfixes in <code>screening.sh</code>|
| 07/23/2021| 	Barrierless process are now included in the search (not in the kinetics)|
| 06/25/2021| 	Screening and rxnetwork build optimized. Adaptive selection of temp in MD optimized.|
| 06/18/2021| 	Bugfix in select_assoc script|
| 06/17/2021| 	Systems with charge can now be modeled with ChemKnow and BXDE|
| 04/06/2021| 	Bugfixes in <code>irc_analysis.sh</code> (mopac2021) and Heuristics (crossed bonds)|
| 05/28/2021| 	New keyword: <code>Use_LET</code> for mopac TS optimizations|
| 05/27/2021| 	New keywords: <code>recalc</code> and <code>Hookean</code>, and improved efficiency|
| 05/14/2021| 	ChemKnow improved and calcfc in <code>g09</code> ts opt calcs|
| 05/10/2021| 	New keyword <code>tight_ts</code> and bug in bxde|
| 05/02/2021| 	kinetics module improved and bugs corrected|
| 04/07/2021| 	In <code>tors.sh</code> do not consider rotation about bonds that belong to a ring|
| 03/16/2021| 	mopac calculator updated to ase3.21.1|
| 03/14/2021| 	Bugfixes|
| 02/15/2021| 	2021 release. Revision 1007--> qcore high-level calcs available|
| 10/20/2020| 	linked_paths in python and 1D rotors in rrkm|
| 10/14/2020| 	tutorial updated and bug fixes in <code>DVV.py</code>|
| 10/05/2020| 	Density matrix read in freq calc (MOPAC) and new xtb parameters for better convergence|
| 07/15/2020| 	Implemented interface with <code>Entos Qcore</code>|
| 04/21/2020| 	Implemented ExtForce and fixed some bugs|
| 01/28/2020| 	Simplified adjacency matrix from XYZ|
| 01/21/2020| 	2020 version|
| 12/02/2019| 	Bug in <code>FINAL.sh</code>|
| 11/22/2019| 	Maximum number of paths set to 50 (in <code>bbfs.f90</code>).|
| 11/16/2019| 	Check of input structure in <code>amk.sh</code>. Fragmented molecules are no longer valid.|
| 09/16/2019| 	pdfs are now also generated in FINAL_HL|
| 07/09/2019| 	if name of working dir is too long, name--->wrkdir|
| 06/30/2019| 	amk acronym replaces old tsscds acronym|
| 04/18/2019| 	MIT license|
| 04/17/2019| 	The label of the starting min in the kmc simulations is in <code>tsdirll/KMC/starting_minimum</code>|
| 04/15/2019| 	A bug in <code>get_energy_g09_MP2.sh</code> was corrected|
		
