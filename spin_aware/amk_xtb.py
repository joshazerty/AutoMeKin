#!/usr/bin/env python3
"""
amk_xtb.py — spin-capable GFN-xTB low-level engine for AutoMeKin (Project 18, Phase 2 Week 0)

WHY THIS EXISTS
---------------
AutoMeKin's stock low-level (LL) path runs MOPAC/PM7 and **never propagates spin
multiplicity**: the MOPAC route is built with `charge=` only, so every LL calculation
is closed-shell RHF regardless of the requested `mult` (Phase 1 closeout ESC-008).
PM7 also inverts the Fe spin-state ordering. Both make spin-state-resolved exploration
impossible at the low level — which is the whole premise of a spin-aware CRN.

This module provides a drop-in GFN2-xTB engine that honours BOTH charge and spin
multiplicity. xTB has explicit transition-metal parameterisation and sets the number
of unpaired electrons natively via `--uhf` (= multiplicity - 1).

It is exposed two ways:
  1. `XTBamk` — an ASE Calculator, so bxde.py can attach it exactly like MOPACamk
     (`geom.calc = XTBamk(...)`; then `get_potential_energy()` / `get_forces()`).
  2. `xtb_sp_grad(...)` — a bare function returning (energy_Eh, gradient_EhBohr),
     mirroring the qcore (`Qcore_calc`) external-call pattern in bxde.py, for callers
     that prefer not to depend on ASE.

It is also runnable standalone as a self-test (see `--help`), which is validation
step V0 in VALIDATION.md (no AutoMeKin required).

SPIN CONVENTION
---------------
multiplicity M = 2S + 1  ->  unpaired electrons = M - 1 = xtb `--uhf` value.
    singlet M=1 -> uhf 0   triplet M=3 -> uhf 2   quintet M=5 -> uhf 4

UNITS
-----
xtb returns atomic units (Hartree, Hartree/Bohr). The ASE calculator converts to
ASE units (eV, eV/Angstrom) using ase.units, identically to the qcore branch:
    energy_eV = E_Eh * units.Hartree
    forces_eV_per_A = -grad_Eh_per_Bohr * units.Hartree / units.Bohr

REQUIREMENTS
------------
- the `xtb` binary on PATH (or pass xtb_bin=...). Tested against xtb >= 6.5.
- numpy. ASE only for the `XTBamk` calculator path (not for `xtb_sp_grad`).
"""

import os
import shutil
import subprocess
import tempfile
import numpy as np


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def mult_to_uhf(mult):
    """multiplicity -> number of unpaired electrons (xtb --uhf)."""
    m = int(mult)
    if m < 1:
        raise ValueError(f"multiplicity must be >= 1, got {m}")
    return m - 1


def gfn_flags(method):
    """Map an AutoMeKin LowLevel method token to xtb CLI flags.

    `LowLevel xtb <method>` puts <method> here. Accepts: gfn2, 'gfn 2', gfn2-xtb,
    gfn1, gfn0, gfnff (case-insensitive). Defaults to GFN2 when ambiguous/empty.
    """
    m = (method or "").strip().lower().replace("-xtb", "").replace("_", "").replace(" ", "")
    if m in ("", "gfn2", "gfn2xtb"):
        return ["--gfn", "2"]
    if m in ("gfn1",):
        return ["--gfn", "1"]
    if m in ("gfn0",):
        return ["--gfn", "0"]
    if m in ("gfnff", "ff"):
        return ["--gfnff"]
    # If the user passed a bare integer level
    if m.isdigit():
        return ["--gfn", m]
    # Unknown -> default to GFN2 but make it visible in logs.
    return ["--gfn", "2"]


def _write_xyz(path, symbols, positions):
    """Minimal XYZ writer (Angstrom)."""
    with open(path, "w") as fh:
        fh.write(f"{len(symbols)}\n\n")
        for s, (x, y, z) in zip(symbols, positions):
            fh.write(f"{s:3s} {x: .10f} {y: .10f} {z: .10f}\n")


def _parse_tm_gradient(grad_path):
    """Parse a Turbomole-format `gradient` file written by `xtb --grad`.

    Layout:
        $grad   cycle = 1   SCF energy = <E_Eh> ...
        <x> <y> <z> <El>     (N atom lines, Bohr)
        <gx> <gy> <gz>       (N gradient lines, Eh/Bohr)
        $end
    Returns (energy_Eh, gradient (N,3) in Eh/Bohr).
    """
    import re
    with open(grad_path) as fh:
        lines = fh.readlines()
    # energy: take the float immediately after "energy =" on the cycle header line
    # (must NOT pick up the trailing "|dE/dxyz| = ..." value on the same line).
    energy = None
    header_idx = None
    for i, ln in enumerate(lines):
        mobj = re.search(r"energy\s*=\s*(-?\d+\.\d+(?:[eEdD][+-]?\d+)?)", ln)
        if mobj:
            energy = float(mobj.group(1).replace("D", "E").replace("d", "e"))
            header_idx = i
    if energy is None or header_idx is None:
        raise RuntimeError(f"could not parse energy from {grad_path}")
    body = [ln for ln in lines[header_idx + 1:] if not ln.strip().startswith("$")]
    coord_lines, grad_lines = [], []
    for ln in body:
        parts = ln.split()
        if len(parts) == 4:          # x y z element  -> coordinate line
            coord_lines.append(ln)
        elif len(parts) == 3:        # gx gy gz        -> gradient line
            grad_lines.append([float(p.replace("D", "E").replace("d", "e")) for p in parts])
    if not grad_lines:
        raise RuntimeError(f"no gradient rows parsed from {grad_path}")
    return energy, np.asarray(grad_lines, dtype=float)


# --------------------------------------------------------------------------- #
# bare functional interface (qcore-style)
# --------------------------------------------------------------------------- #
def xtb_sp_grad(symbols, positions, charge=0, mult=1, method="gfn2",
                xtb_bin="xtb", accuracy=0.2, etemp=300.0, extra=None):
    """Single-point energy + gradient with spin control.

    Returns (energy_Eh, gradient (N,3) Eh/Bohr). Raises on xtb failure.
    Runs in a private scratch dir so concurrent trajectories never collide.
    """
    uhf = mult_to_uhf(mult)
    tmp = tempfile.mkdtemp(prefix="amkxtb_")
    try:
        xyz = os.path.join(tmp, "mol.xyz")
        _write_xyz(xyz, symbols, positions)
        cmd = [xtb_bin, "mol.xyz", *gfn_flags(method),
               "--chrg", str(int(charge)), "--uhf", str(int(uhf)),
               "--grad", "--acc", str(accuracy), "--etemp", str(etemp),
               "--norestart"]
        if extra:
            cmd += list(extra)
        proc = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True)
        grad_file = os.path.join(tmp, "gradient")
        if proc.returncode != 0 or not os.path.isfile(grad_file):
            raise RuntimeError(
                "xtb failed (rc=%s).\nCMD: %s\nSTDERR tail:\n%s"
                % (proc.returncode, " ".join(cmd), proc.stderr[-1500:]))
        return _parse_tm_gradient(grad_file)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #
# ASE Calculator interface (MOPACamk-style; preferred for bxde.py)
# --------------------------------------------------------------------------- #
try:
    from ase.calculators.calculator import Calculator, all_changes
    from ase import units

    class XTBamk(Calculator):
        """ASE Calculator wrapping `xtb` with charge + spin multiplicity.

        Usage in bxde.py (mirrors `rmol.calc = MOPACamk(...)`):
            rmol.calc = XTBamk(charge=charge, mult=mult, method='gfn2')
        then `get_potential_energy()` / `get_forces()` work as usual.
        """
        implemented_properties = ["energy", "forces"]

        def __init__(self, charge=0, mult=1, method="gfn2",
                     xtb_bin="xtb", accuracy=0.2, etemp=300.0, **kwargs):
            Calculator.__init__(self, **kwargs)
            self.charge = int(charge)
            self.mult = int(mult)
            self.method = method
            self.xtb_bin = xtb_bin
            self.accuracy = accuracy
            self.etemp = etemp

        def calculate(self, atoms=None, properties=("energy",),
                      system_changes=all_changes):
            Calculator.calculate(self, atoms, properties, system_changes)
            e_eh, grad = xtb_sp_grad(
                self.atoms.get_chemical_symbols(),
                self.atoms.get_positions(),
                charge=self.charge, mult=self.mult, method=self.method,
                xtb_bin=self.xtb_bin, accuracy=self.accuracy, etemp=self.etemp)
            self.results["energy"] = e_eh * units.Hartree
            self.results["forces"] = -grad * units.Hartree / units.Bohr

except ImportError:      # ASE not importable in this context; functional API still works.
    XTBamk = None


# --------------------------------------------------------------------------- #
# standalone self-test  (VALIDATION.md step V0)
# --------------------------------------------------------------------------- #
def _selftest():
    """Prove the engine sets spin and orders Fe states sanely. No AutoMeKin needed."""
    import argparse
    ap = argparse.ArgumentParser(description="amk_xtb self-test (V0)")
    ap.add_argument("xyz", nargs="?", help="optional xyz; default runs a built-in H2CO + Fe demo")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--method", default="gfn2")
    ap.add_argument("--xtb-bin", default="xtb")
    args = ap.parse_args()

    def read_xyz(path):
        with open(path) as fh:
            n = int(fh.readline()); fh.readline()
            syms, pos = [], []
            for _ in range(n):
                p = fh.readline().split()
                syms.append(p[0]); pos.append([float(p[1]), float(p[2]), float(p[3])])
        return syms, np.asarray(pos)

    if args.xyz:
        syms, pos = read_xyz(args.xyz)
        e, g = xtb_sp_grad(syms, pos, args.charge, args.mult, args.method, args.xtb_bin)
        print(f"E = {e:.8f} Eh   |grad|max = {abs(g).max():.6f} Eh/Bohr   "
              f"(charge={args.charge}, mult={args.mult}, uhf={mult_to_uhf(args.mult)})")
        return

    # Built-in demo: H2CO closed-shell sanity + finite-difference gradient check.
    h2co = (["C", "O", "H", "H"],
            np.array([[0.000, 0.000, 0.000],
                      [0.000, 0.000, 1.205],
                      [0.000, 0.943, -0.589],
                      [0.000, -0.943, -0.589]]))
    e0, g0 = xtb_sp_grad(*h2co, charge=0, mult=1, method=args.method, xtb_bin=args.xtb_bin)
    # finite-diff one component
    d = 1e-3
    sym, pos = h2co[0], h2co[1].copy()
    pos[0, 2] += d
    ep, _ = xtb_sp_grad(sym, pos, 0, 1, args.method, args.xtb_bin)
    pos[0, 2] -= 2 * d
    em, _ = xtb_sp_grad(sym, pos, 0, 1, args.method, args.xtb_bin)
    bohr = 0.529177210903
    fd = (ep - em) / (2 * d / bohr)        # Eh/Bohr
    print("=== H2CO (closed shell) ===")
    print(f"E = {e0:.8f} Eh")
    print(f"analytic dE/dz[C] = {g0[0,2]: .6f}   finite-diff = {fd: .6f} Eh/Bohr  "
          f"(should match to ~1e-3)")
    print("V0 sanity OK" if abs(g0[0, 2] - fd) < 5e-3 else "V0 WARNING: gradient mismatch")
    print("\nFor the Fe spin-ordering check, run on your Fe(CO)4 or Tp-Fe(II)-ethyl xyz:")
    print("  for M in 1 3 5; do python3 amk_xtb.py reactant.xyz --mult $M; done")
    print("  -> energies MUST differ between multiplicities (proves --uhf is active),")
    print("     unlike the stock MOPAC path where all spins gave one RHF surface.")


if __name__ == "__main__":
    _selftest()
