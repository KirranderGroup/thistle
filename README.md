# THISTLE

**Ab initio X-ray and electron scattering from electronic wavefunctions.**

THISTLE computes what gas-phase X-ray scattering and electron diffraction experiments measure,
directly from the wavefunction of a quantum chemistry calculation, without the independent atom
model. It takes a Gaussian basis and the one- and two-particle density matrices of any method
(Hartree–Fock, CASSCF, CI, coupled cluster) and returns:

| | X-rays | Electrons (UED, GED) |
|---|---|---|
| Rotational average (gas phase) | total, elastic, inelastic; Hartree, exchange and correlation parts | dσ/dΩ: total, elastic, inelastic; nucleus–nucleus, electron–nucleus, electron–electron parts |
| Partially aligned molecules | j₂ moment S₂(q) | j₂ moment |
| Fixed orientation | detector images, any q vectors | detector images for a given beam energy |
| Reference | independent atom model (ab initio or tables), % difference | IAM, % difference, modified intensity sM(s) |

All integrals are analytic (McMurchie–Davidson with spherical Bessel derivatives, Hobson's
theorem for j₂, analytic Fourier transforms at fixed orientation), screened and parallel
(Fortran + OpenMP), and tested against independent routes to 10⁻⁸ or better. Wavefunctions come
from PySCF objects, Molden files, OpenMolcas and Molpro outputs, or your own CI vectors.

```bash
pip install .                              # compiles the Fortran kernels (needs gfortran)
thistle run lif.out --state 2 --iam        # Molpro/Molcas output or Molden file -> CSV + PNG
thistle run lif.out --electrons --iam      # electron diffraction, with sM(s)
thistle new lif.toml; thistle run lif.toml # the same with a commented input file
thistle submit lif.toml --now              # ... as a Slurm job
thistle plot results/lif_scattering.csv --iam results/lif_iam.csv   # plots from the result files
```

```python
import thistle as th
lif = th.load("lif.out", state=2)                    # or a Molden file, or a PySCF object
lif.scattering(qmax=4).plot("lif.png")               # X-rays, q in 1/Angstrom
th.compare(lif.scattering(), lif.iam()).plot()       # % difference from free atoms
e = lif.scattering(qmax=10, probe="electron")        # electron diffraction cross section
th.sm(e, lif.iam(qmax=10, probe="electron")).plot()  # modified intensity sM(s)
```

* **New users:** [docs/QUICKSTART.md](docs/QUICKSTART.md)
* **The physics, the algorithms and the validation:** [docs/THEORY.md](docs/THEORY.md)
* **Development log:** [docs/STATUS.md](docs/STATUS.md)

The rest of this file describes `pyxscat`, the layer underneath `thistle`, for full control.

## PyXSCAT: the kernels

The scattering kernels take a Gaussian basis set and density matrices and nothing else, so any
electronic structure package can be used through a thin adapter:

| Source | Adapter | Gives |
|---|---|---|
| PySCF (in memory) | `interfaces.pyscf` | SCF, CCSD, CISD, FCI, CASSCF: 1RDM and 2RDM |
| Any Molden file (Molpro, Molcas, BAGEL, ORCA, PySCF, ...) | `interfaces.molden` | basis, orbitals, occupations: 1RDM |
| OpenMolcas RASSCF/CASSCF output + its Molden file | `interfaces.molcas` | 1RDM and 2RDM of any printed root |
| Molpro MULTI output (`config,det`) or punch file + its Molden file | `interfaces.molpro` | 1RDM, 2RDM and transition 1RDMs of every state |
| Your own CI vector | `pyxscat.ci` | 1RDM and 2RDM from determinants or GUGA CSFs |

The BAGEL and NEXMD readers of the earlier code base are still to be ported.

```python
import numpy as np
from pyscf import gto, scf, cc
import pyxscat as xs

mol = gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="cc-pvdz")
mycc = cc.CCSD(scf.RHF(mol).run()).run()

q = np.linspace(0, 10, 101)                      # inverse bohr
wf = xs.interfaces.pyscf.from_cc(mycc)
total = xs.scattering(wf, q, "total")            # also "elastic", "inelastic", "two_electron"
parts = xs.decompose(wf, q)                      # self + hartree + exchange + cumulant
```

From the command line (batch and cluster jobs):

```bash
pyxscat --molden scf.molden --kind elastic --qmax 10 --npoints 201 -o elastic.dat
pyxscat --molden casscf.rasscf.molden --molcas-log casscf.log --root 1 --kind decompose -o cas.dat
pyxscat --molden lif.mld --molpro-out lif.out --root 2 --qvec q.txt --kind total -o lif_s2.dat
pyxscat --molden sf6.molden --molpro-punch sf6.pun --molpro-out sf6.out --kind decompose -o sf6.dat
```

## Install

```bash
pip install ".[all]"      # compiles python/pyxscat/libpyxscat.{so,dylib} (gfortran + OpenMP) and installs
```

For development, build in place and install in editable mode:

```bash
make                      # builds python/pyxscat/libpyxscat.{so,dylib}
pip install -e ".[test]"
python -m pytest          # about 4 min; PYXSCAT_SLOW=1 adds the large-basis tests
python -m build --sdist   # the source package (thistle-<version>.tar.gz) to hand out
```

Intel oneAPI: load its environment (`source /opt/intel/oneapi/setvars.sh` or `module load intel`),
then `FC=ifx pip install .` (or `make FC=ifx`); `EXTRA=-xHost` adds flags. The Intel build uses
`-fp-model precise` and links the Intel runtime statically; keep the oneAPI environment loaded at
run time for the OpenMP runtime (libiomp5).

On macOS, point `FC` at a gfortran that supports OpenMP (Homebrew `gcc` or xPack).

**Testers:** see [docs/FOR_TESTERS.md](docs/FOR_TESTERS.md).

## What is computed

All quantities are rotationally averaged (j0) and given in atomic units.

| Function | Quantity | Needs |
|---|---|---|
| `elastic(basis, q, D)` | ⟨\|F(q)\|²⟩ | 1RDM |
| `elastic(basis, q, Da, Db)` | ⟨F_a F_b*⟩ (transition and coherent mixed terms) | two 1-particle matrices |
| `two_electron(basis, q, Γ)` | Σ Γ ⟨ij\|j0(q r12)\|kr⟩ | 2RDM |
| `mean_field(basis, q, D)` | two-electron term of DD − ½DD (exact for closed-shell HF) | 1RDM only |
| `exchange(basis, q, D)` | −½ D D exchange part | 1RDM only |
| `scattering(wf, q, kind)` | total = N + two-electron, elastic, inelastic | |
| `decompose(wf, q)` | total = N + Hartree + exchange + cumulant | 1RDM and 2RDM |

Conventions: densities are spin-summed, in the AO basis of the `Basis` object. The 2RDM is in
chemists' notation with trace N(N−1), the same as PySCF's `make_rdm2`.

The 2RDM may also be a `FactoredDM2`, a sum of products of one-particle matrices
(Γ = Σ A_k⊗B_k in Coulomb and exchange pairings) that the kernels contract without forming the
four-index array. `pyxscat.ccsd.factored_rdm2(mycc)` writes the CCSD density in that form
(7 n_occ² + 4 terms, memory 2 K n² instead of n⁴, exact to 10⁻¹³); `interfaces.pyscf.from_cc`
uses it by default beyond 60 basis functions (`factored=True/False` to choose).

## Partially aligned molecules: the j2 moment

```python
s2 = xs.scattering(wf, q, "total", l=2)       # -sum Gamma <j2(q r12) P2(cos theta_12)> (i^2 = -1), z = alignment axis
iam2 = iam.for_wavefunction(wf, q, l=2)["total"]
# I(q, theta_q) = I0(q) + 5 <P2(cos beta)> P2(cos theta_q) S2(q) + ...
```

On a detector grid (the polar grid of `grid.f90`, a q-grid file or any list of vectors):

```python
from pyxscat import grids
qlab = grids.qcoord3(10.0, 60.0, 49, 73)           # e0 (keV), thsmax (deg), n1, n2: as grid.f90
qn, uq, inv = grids.unique_norms(qlab)             # one rotational average per ring
i0 = xs.scattering(wf, uq, "total"); s2 = xs.scattering(wf, uq, "total", l=2)
img = grids.aligned_ensemble(qlab, axis=(0, 1, 0), i0=i0[inv], s2=s2[inv], p2=0.5)
# or simply: thistle.load(...).aligned_ensemble(axis=(0, 1, 0), p2=0.5, grid="polar", energy_kev=10)
# the IAM on the same grid: model="iam"; thistle.compare(pattern, iam_pattern) gives the % map
```

## Fixed orientation (aligned molecules, detector images)

```python
from pyxscat import aligned
k0 = 2 * np.pi / (1.0 / 0.529177)                  # 1 Angstrom X-rays, bohr^-1
q, on_detector = aligned.detector_q(k0, qmax=4.0, npix=101)
R = aligned.rotation_matrix(0.0, np.pi / 2, 0.0)   # molecular orientation (zyz Euler angles)
img = aligned.decompose(wf, aligned.rotate_q(q[on_detector], R))   # total, elastic, cumulant, ...
```

The same quantities as the rotational average are available (`elastic`, `two_electron`,
`exchange`, `mean_field`, `scattering`, `decompose`), built from analytic Fourier transforms of
basis-function products.

## Electron scattering

```python
from pyxscat import electron
d = electron.decompose(wf, q)                     # nn, ne, ee_elastic, ee_total, elastic, inelastic, total
cs = electron.cross_section(d["total"], q)        # 4 [...] / q^4: dsigma/dOmega in bohr^2/sr
s2 = electron.decompose(wf, q, l=2)["total"]      # j2 moment about z (aligned molecules)
img = electron.aligned(wf, qvec)                  # fixed orientation: |Z(q) - F(q)|^2 + inelastic
ne = electron.electron_nucleus(wf.basis, q, wf.dm1, charges, coords)   # the new one-electron integrals
```

The bracket I_nn − 2 I_ne + I_ee is the quantity of the earlier code; `hf=True` gives total and
inelastic for a closed-shell determinant without a 2RDM. Command line: `pyxscat ... --probe electron`.

## Independent atom model

```python
from pyxscat import iam
ref = iam.for_wavefunction(wf, q, basis="aug-cc-pvtz")          # free-atom HF, same basis
pd = iam.percent_difference(xs.scattering(wf, q, "total"), ref["total"])
ref_al = iam.for_wavefunction(wf, qvec, aligned=True)            # fixed orientation
f_C = iam.atomic_form_factor("C", q, source="tables")             # Waasmaier-Kirfel (xraydb)
e_ref = iam.electron_for_wavefunction(wf, q)                     # electrons: Mott-Bethe amplitudes Z - f
```

## Numerical controls

| Option | Default | Meaning |
|---|---|---|
| `screen_tol` | 1e-12 | Shell-quartet and primitive-quartet screening (0 = none) |
| `cutoffcentre` | 1e-10 bohr | Quartet route: pair-centre distance below which the h = 0 closed form is used (the radial functions are stable at any distance, so this is not an approximation) |
| `method` (elastic) | `"pairs"` | `"pairs"`: density folded into each pair first (exact at any pair distance, faster); `"quartets"`: the one-loop shell-quartet kernel |
| `cutoffz`, `cutoffmd` | 1e-20 | Pair-density and Hermite-product cutoffs |

## Layout

```
src/fortran/   xs_hermite (Hermite coefficients, any l), xs_bessel_j0, xs_p0 (h = 0 limits),
               xs_kernel_j0 (screened shell-quartet kernel), xs_pairft (fixed-q pair
               transforms), xs_elastic (elastic by density-weighted pairs), xs_nuclear
               (electron-nucleus integrals), xs_radial (j_n(x)/x^n), xs_capi (C interface)
python/thistle the simple interface: load / Molecule / Result, input files, `thistle` command,
               plotting (a plotter for every calculation, reading back every file written)
python/pyxscat basis, core (API), aligned (fixed orientation), grids (detector grids, aligned ensembles), electron, iam, ci (CI vectors -> RDMs), ccsd (factored CCSD density), cli,
               interfaces/{pyscf, molden, molcas, molpro}, transforms
docs/          QUICKSTART (users), THEORY (physics, algorithms, validation), STATUS (development log)
tests/         unit, physics, consistency and regression tests
examples/      quickstart.py (PySCF), lif.toml and lif_electrons.toml (input files), lower-level scripts
```

## Authors

Andrés Moreno Carrascosa, Mats Simmermacher, Patrick Yuheng Wang, Eirik Moen Liane, Joseph Cooper, Nikola Zotev, Thomas Northey and Adam Kirrander.

## Licence

Apache License 2.0 (see `LICENSE`); copyright 2026 Andrés Moreno Carrascosa. Parts of
`pyxscat/ccsd.py` are derived from PySCF (Apache 2.0, see `NOTICE`). To cite THISTLE, see
`CITATION.cff`.
