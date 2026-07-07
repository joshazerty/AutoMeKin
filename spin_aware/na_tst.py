#!/usr/bin/env python3
"""
na_tst.py — Stage M4 physics core: nonadiabatic TST rate for a spin crossing (Project 18).

Computes the rate of a spin-forbidden step at a localised MECP (from M2/M3) using
nonadiabatic transition state theory with a Landau–Zener surface-hopping transmission
coefficient — the Harvey framework (Harvey, PCCP 2007, 9, 331; Lykhin et al., IJQC 2016).

WORKING EQUATION (canonical, unimolecular spin flip)
----------------------------------------------------
    k(T) = kappa(T) · (k_B T / h) · (Q_MECP / Q_react) · exp(-dE_MECP / k_B T)

  * dE_MECP : MECP energy relative to the reactant, INCLUDING ZPE (kcal/mol in, J inside).
  * Q_react : reactant vibrational partition function (3N-6 real modes).
  * Q_MECP  : MECP "seam" vibrational partition function (3N-7 real modes — the coordinate
              orthogonal to the seam, along the gradient-difference vector, is removed).
  * kappa(T): Boltzmann-averaged Landau–Zener hopping (transmission) coefficient, <= 1:
                kappa(T) = (1/k_B T) ∫_0^inf p_sh(E_perp) e^{-E_perp/k_B T} dE_perp
              with the single-passage LZ probability of STAYING diabatic
                P_LZ(E_perp) = exp(-2*pi*Gamma),   2*pi*Gamma = 2*pi*H_SO^2 / (hbar |dF| v),
                v = sqrt(2 E_perp / mu),
              and the surface-hop probability per attempt:
                single passage : p_sh = 1 - P_LZ
                double passage : p_sh = (1 - P_LZ)(1 + P_LZ) = 1 - P_LZ^2   (default; Harvey)

INPUTS (physical, with the units the CLI/JSON expects)
  H_SO   spin-orbit coupling matrix element   [cm^-1]
  dF     |grad_A - grad_B| at the MECP        [Eh/Bohr]   (difference of the two surface gradients)
  mu     effective mass for the hop coordinate [amu]      (default 1.0; sensitivity-swept)
  dE     MECP energy vs reactant (elec or +ZPE) [kcal/mol]
  freqs_react, freqs_mecp                      [cm^-1]    (real modes; MECP has one fewer)

Everything is converted to SI internally. numpy only.
"""
import json
import sys
import argparse
import numpy as np

# ---- SI constants ----
h = 6.62607015e-34          # J s
hbar = 1.054571817e-34      # J s
kB = 1.380649e-23           # J/K
c_cm = 2.99792458e10        # cm/s
amu = 1.66053906660e-27     # kg
Eh = 4.3597447222071e-18    # J
a0 = 5.29177210903e-11      # m
KCAL = 4184.0 / 6.02214076e23      # J per molecule per kcal/mol
CM1_J = h * c_cm                    # J per cm^-1
EHB_N = Eh / a0                     # N per (Eh/Bohr)


# --------------------------------------------------------------------------- #
# Landau–Zener
# --------------------------------------------------------------------------- #
def p_lz_stay(E_perp_J, H_SO_J, dF_N, mu_kg):
    """Single-passage Landau–Zener probability of REMAINING on the diabatic surface."""
    E = np.maximum(E_perp_J, 1e-20 * Eh)        # avoid v=0 singularity
    v = np.sqrt(2.0 * E / mu_kg)
    two_pi_gamma = 2.0 * np.pi * H_SO_J**2 / (hbar * dF_N * v)
    return np.exp(-two_pi_gamma)


def p_hop(E_perp_J, H_SO_J, dF_N, mu_kg, passage="double"):
    P = p_lz_stay(E_perp_J, H_SO_J, dF_N, mu_kg)
    if passage == "single":
        return 1.0 - P
    return (1.0 - P) * (1.0 + P)               # double passage (default)


def kappa(T, H_SO_cm, dF_eh_bohr, mu_amu, passage="double", n=4000):
    """Boltzmann-averaged surface-hopping transmission coefficient (dimensionless, <=1)."""
    H_SO_J = H_SO_cm * CM1_J
    dF_N = dF_eh_bohr * EHB_N
    mu_kg = mu_amu * amu
    kT = kB * T
    # integrate E_perp from ~0 to 60 kT
    E = np.linspace(1e-6 * kT, 60.0 * kT, n)
    integrand = p_hop(E, H_SO_J, dF_N, mu_kg, passage) * np.exp(-E / kT)
    trapz = getattr(np, "trapezoid", np.trapz)
    return float(trapz(integrand, E) / kT)


# --------------------------------------------------------------------------- #
# partition functions / ZPE
# --------------------------------------------------------------------------- #
def q_vib(freqs_cm, T):
    """Vibrational partition function (ground-state referenced, harmonic)."""
    q = 1.0
    for v in freqs_cm:
        if v <= 0:
            continue                            # skip imaginary/zero (seam coordinate already removed)
        x = CM1_J * v / (kB * T)
        q *= 1.0 / (1.0 - np.exp(-x))
    return q


def zpe_kcal(freqs_cm):
    return 0.5 * sum(v for v in freqs_cm if v > 0) * CM1_J / KCAL


# --------------------------------------------------------------------------- #
# the rate
# --------------------------------------------------------------------------- #
def k_natst(T, dE_kcal, H_SO_cm, dF_eh_bohr, mu_amu,
            freqs_react=None, freqs_mecp=None, passage="double",
            add_zpe=False):
    """NA-TST spin-crossing rate (s^-1). If add_zpe, dE_kcal is electronic and ZPE
    (from freqs) is added; otherwise dE_kcal is taken as already ZPE-corrected."""
    kp = kappa(T, H_SO_cm, dF_eh_bohr, mu_amu, passage)
    qr = q_vib(freqs_react or [], T)
    qm = q_vib(freqs_mecp or [], T)
    dE = dE_kcal
    if add_zpe and freqs_react and freqs_mecp:
        dE = dE_kcal + (zpe_kcal(freqs_mecp) - zpe_kcal(freqs_react))
    pref = (kB * T / h) * (qm / qr)
    k = kp * pref * np.exp(-dE * KCAL / (kB * T))
    return {"T": T, "k_s^-1": k, "kappa": kp, "prefactor_s^-1": pref,
            "Qmecp_over_Qreact": qm / qr, "dE_used_kcal": dE,
            "passage": passage}


def k_sweep(temps, **kw):
    return [k_natst(T, **kw) for T in temps]


# --------------------------------------------------------------------------- #
# CLI + selftest
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description="M4 NA-TST/Landau-Zener spin-crossing rate")
    ap.add_argument("--inputs", help="JSON with dE_kcal,H_SO_cm,dF_eh_bohr,mu_amu,freqs_react,freqs_mecp")
    ap.add_argument("--T", type=float, default=298.15)
    ap.add_argument("--tsweep", help="comma T list, e.g. 200,250,298.15,350")
    ap.add_argument("--passage", default="double", choices=["single", "double"])
    ap.add_argument("--soc-sweep", help="comma H_SO cm^-1 list for sensitivity")
    ap.add_argument("--out", default=None)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return _selftest()
    d = json.load(open(a.inputs))
    base = dict(dE_kcal=d["dE_kcal"], H_SO_cm=d["H_SO_cm"], dF_eh_bohr=d["dF_eh_bohr"],
                mu_amu=d.get("mu_amu", 1.0), freqs_react=d.get("freqs_react"),
                freqs_mecp=d.get("freqs_mecp"), passage=a.passage,
                add_zpe=d.get("add_zpe", False))
    temps = [float(x) for x in a.tsweep.split(",")] if a.tsweep else [a.T]
    res = {"base": k_sweep(temps, **base)}
    if a.soc_sweep:
        res["soc_sweep_at_T"] = temps[0]
        res["soc_sweep"] = []
        for hso in [float(x) for x in a.soc_sweep.split(",")]:
            b2 = dict(base); b2["H_SO_cm"] = hso
            r = k_natst(temps[0], **b2); r["H_SO_cm"] = hso
            res["soc_sweep"].append(r)
    out = json.dumps(res, indent=2)
    print(out)
    if a.out:
        open(a.out, "w").write(out)


def _selftest():
    # 1) units/sanity: a reasonable Fe crossing gives a finite, sub-TST rate
    r = k_natst(298.15, dE_kcal=19.56, H_SO_cm=200.0, dF_eh_bohr=0.05, mu_amu=10.0,
                freqs_react=[200, 400, 600, 900, 1200, 1500],
                freqs_mecp=[210, 410, 620, 950, 1250])
    assert 0 < r["kappa"] <= 1.0, r["kappa"]
    assert r["k_s^-1"] > 0
    # 2) weak-coupling scaling: at FIXED energy, p_hop ∝ H_SO^2 (LZ linearisation).
    #    (kappa itself does NOT scale cleanly because its Boltzmann average always
    #     includes the low-velocity region where p_hop saturates to 1 regardless of H_SO.)
    E_fixed = 5.0 * kB * 298.15            # well inside the weak regime for small H_SO
    p1 = p_hop(E_fixed, 1.0 * CM1_J, 0.05 * EHB_N, 10.0 * amu)
    p2 = p_hop(E_fixed, 2.0 * CM1_J, 0.05 * EHB_N, 10.0 * amu)
    ratio = p2 / p1
    assert 3.9 < ratio < 4.1, f"weak-coupling H_SO^2 scaling broken: {ratio}"
    # 3) adiabatic (strong-coupling) limit: huge H_SO -> kappa -> 1 (double passage saturates)
    k_strong = kappa(298.15, 1e5, 0.05, 10.0)
    assert k_strong > 0.99, k_strong
    # 4) monotonic in T for an activated process: higher T -> higher k
    ks = [k_natst(T, 19.56, 200.0, 0.05, 10.0, [200,400,600,900,1200,1500],
                  [210,410,620,950,1250])["k_s^-1"] for T in (250, 298.15, 350)]
    assert ks[0] < ks[1] < ks[2], ks
    # 5) double passage >= single passage hop probability
    assert kappa(298.15,200,0.05,10.0,"double") >= kappa(298.15,200,0.05,10.0,"single")
    print("na_tst selftest PASS")
    print(f"  example Fe crossing @298 K: kappa={r['kappa']:.3f}, "
          f"k={r['k_s^-1']:.3e} s^-1, Qmecp/Qreact={r['Qmecp_over_Qreact']:.3f}")
    print(f"  weak-coupling H_SO^2 scaling ratio (x2 SOC) = {ratio:.2f} (expect ~4)")
    print(f"  strong-coupling kappa -> {k_strong:.4f} (expect ~1)")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
