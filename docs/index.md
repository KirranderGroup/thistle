# THISTLE

**Ab initio X-ray and electron scattering from electronic wavefunctions.**

THISTLE computes what gas-phase X-ray scattering and electron diffraction experiments measure,
directly from the wavefunction of a quantum chemistry calculation, without the independent atom
model.  It takes a Gaussian basis and the one- and two-particle density matrices of any method
(Hartree–Fock, CASSCF, CI, coupled cluster) and returns rotationally averaged, partially aligned
and fixed-orientation scattering for X-rays and electrons, with the independent atom model as a
reference.

The code has three layers:

* `thistle`, the simple interface: load an output file, get a `Molecule`, ask it for scattering,
  plot and save the result; input files, the `thistle` command and Slurm submission.
* `pyxscat`, the kernel-level API: basis sets, density matrices, every scattering kind, the
  adapters for PySCF, Molden, OpenMolcas and Molpro, CI and CCSD densities.
* The Fortran kernels (`src/fortran`), compiled into `libpyxscat` and called through a C
  interface: the analytic integrals, screened and parallel (OpenMP).

```{toctree}
:maxdepth: 2
:caption: Using THISTLE

QUICKSTART
FOR_TESTERS
THEORY
```

```{toctree}
:maxdepth: 2
:caption: Reference

api/index
kernels
```

```{toctree}
:maxdepth: 1
:caption: Development

STATUS
```
