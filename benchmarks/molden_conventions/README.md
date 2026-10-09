# Checking a program's Molden conventions

Programs differ in how they write Molden files: Cartesian or spherical functions, the order and
sign of the components, the normalisation of Cartesian d/f/g functions, and the order of the
orbitals. `pyxscat.interfaces.molden.scf_energy` rebuilds the SCF energy from a Molden file with
exact integrals. It reproduces the program's SCF energy only if every one of these conventions is
read correctly, so it is a complete test.

Run the same closed-shell Hartree–Fock calculation (water, cc-pVTZ, which has d and f functions,
no symmetry) in each program, then:

```bash
python check.py water.scf.molden -76.057127420     # the program's SCF energy
```

| Program | Input | Molden file it writes | Expected agreement |
|---|---|---|---|
| OpenMolcas | `water_molcas.input` | `water.scf.molden` (spherical, `[5D7F]`) | 1e-7 Eh (exact integrals, no `ricd`) |
| BAGEL | `water_bagel.json` | `water_bagel.molden` | ~1e-4 Eh (BAGEL fits integrals) |
| Molpro | `water_molpro.inp` | `water_molpro.molden` (Cartesian) | 1e-7 Eh |

A convention error shows up as an energy error of millihartree or more, not 1e-5.
PySCF-written files already pass (`tests/test_molden.py::test_scf_energy_check`).
