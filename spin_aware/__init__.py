"""
spin_aware — spin-aware extensions to AutoMeKin for two-state reactivity.

Modules
-------
amk_xtb          : spin-capable GFN2-xTB low-level engine (multiplicity -> --uhf)
spin_enumerate   : M0/M1 — <S^2> gate + close-lying-spin-state (crossing) detection
locate_mecp      : M2 — minimum-energy crossing-point localisation (easyMECP/Gaussian wrapper)
mecp_graph       : M3 — insert MECP nodes into the reaction-network graph
na_tst           : M4 physics core — nonadiabatic (Landau-Zener/WKB) TST hop rate
microkinetics    : M4 post-analysis — kinetic role of the spin crossing
annotate_rates   : M4 — attach NA-TST rates to the M3 graph

See README_spin_aware.md for the pipeline overview and the two worked
validation cases (Tp-Fe(II)-ethyl beta-H elimination; Fe(I) alkene isomerisation).
"""
__all__ = ["amk_xtb","spin_enumerate","locate_mecp","mecp_graph",
           "na_tst","microkinetics","annotate_rates"]
__version__ = "0.1.0"
