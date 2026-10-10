# THISTLE quick start

*Authors: Andrés Moreno Carrascosa, Mats Simmermacher, Patrick Yuheng Wang, Eirik Moen Liane, Joseph Cooper, Nikola Zotev, Thomas Northey and Adam Kirrander.*

THISTLE computes X-ray and electron scattering (gas-phase X-ray scattering, ultrafast and gas
electron diffraction) from the wavefunction of a quantum chemistry calculation.
It can be used in three ways. Pick whichever suits you; they all give the same numbers.

| You want to... | Use |
|---|---|
| get a curve from an output file you already have, now | `thistle run my.out` |
| keep a record of the settings, or run on a cluster | an input file: `thistle new`, `thistle run`, `thistle submit` |
| script, loop over geometries, make your own figures | Python: `import thistle as th` |

## 1. Install (once)

You need Python 3.10 or newer and a Fortran compiler (gfortran).

```bash
# get gfortran:  macOS: brew install gcc   Ubuntu: sudo apt install gfortran   cluster: module load gcc
pip install thistle-0.2.0.dev0.tar.gz          # or, from the source folder: pip install .
# Intel compiler instead:  source /opt/intel/oneapi/setvars.sh;  FC=ifx pip install thistle-0.2.0.dev0.tar.gz
pip install pyscf scipy xraydb                 # optional: PySCF input and the independent atom model
thistle --help
```

`pip` compiles the Fortran kernels for your machine during the install. On a cluster, do this
once in your home folder (`python -m venv ~/thistle-env; source ~/thistle-env/bin/activate`
first), then reuse the environment in your jobs.

## 2. Straight from an output file

```bash
thistle run lif.out                     # Molpro: lif.mld (and lif.pun, if there) are found by themselves
thistle run lif.out --state 2 --iam     # second state, and the comparison with free atoms
thistle run casscf.rasscf.molden        # any Molden file: elastic scattering only
thistle run water_scf.molden --hf       # Molden file from closed-shell Hartree-Fock: total and inelastic too
thistle run c7h8.log --all-parts        # OpenMolcas: also Hartree, exchange and correlation parts
thistle run lif.out --electrons --iam   # electron diffraction: cross section, and sM(s) against the IAM
thistle run lif.out --probe both        # X-rays and electrons in one go
```

Each run writes a table (`lif_scattering.csv`, opens in Excel or any plotting program) and a
picture (`lif_scattering.png`). q is in inverse Angstrom unless you add `--unit 1/bohr`.

What THISTLE reads, and what it can compute from it:

| File | What it contains | Gives |
|---|---|---|
| Molden file | orbitals and occupations | elastic; with `--hf` (closed-shell Hartree–Fock only) also total and inelastic |
| Molpro MULTI output (`config,det`) or `.pun`, with its Molden file (the punch and Molden files alone are enough) | CI vectors | total, elastic, inelastic, all parts, every state |
| OpenMolcas RASSCF output, with its `.rasscf.molden` | CI vectors | the same |
| PySCF object (Python only) | SCF, CCSD, CISD, CASSCF | the same |
| thistle-dev input deck (`deck.h5`, orbital or dense-AO form) | basis, 1RDM, 2RDM | total, elastic, inelastic, all parts |

Why `--hf`: a Molden file holds orbitals and occupations, i.e. the one-electron density, which is
all the elastic scattering needs. Total and inelastic scattering need the pair density. For a
closed-shell Hartree–Fock determinant the pair density follows exactly from the orbitals, but
for natural orbitals of a correlated calculation (CASSCF, CI) it does not, and the file does not
say which kind it holds. So you say it: `--hf` (or `--scf`; `hf = true` in an input file). THISTLE
refuses it if the occupations are not 0 or 2.

If you are unsure that a file is read correctly:

```bash
thistle check lif.mld --scf-energy <E from output>  # orbitals orthonormal? energy reproduced? I(0) = N^2?
```

## Electrons

Everything works for electrons too: add `--electrons` on the command line, `probe = "electron"`
(or `"both"`) in an input file, or `probe="electron"` in Python. What changes:

* The result is the differential cross section dσ/dΩ in Å²/sr (bohr²/sr with `--unit 1/bohr`),
  first Born approximation. It diverges at q = 0, so the first row is empty (NaN).
* `--energy 3000` (keV) multiplies by the relativistic factor γ² of the beam electrons and sets
  the wavelength of detector images. Without it the curves are the non-relativistic ones.
* With `--iam` you also get `_electron_sM.csv/.png`: the modified scattering
  sM(s) = s (I − I_at)/I_at that UED experiments report, for the molecule and for the
  independent atom model. Their difference is the effect of chemical bonding (and correlation).
* `--all-parts` writes the nucleus–nucleus, electron–nucleus and electron–electron terms to the
  table. They cancel almost completely at small q, so they are not plotted.
* Elastic electron scattering needs only the orbitals; total and inelastic need the pair
  density, exactly as for X-rays (so `--hf` again for a Hartree–Fock Molden file).

## 3. With an input file

```bash
thistle new lif.toml      # writes a commented file: fill in the name of your output
thistle run lif.toml      # computes everything it asks for
```

The input file is plain text (see `examples/lif.toml`):

```toml
[input]
file  = "lif_sa4_casscf.out"
state = 1

[q]
unit   = "1/A"
max    = 4.0
points = 201

[compute]
probe      = "xray"                # "xray", "electron" or "both"
# energy_kev = 3000                # electron beam energy (relativistic factor, detector wavelength)
scattering = true
parts      = "standard"            # "all" adds Hartree, exchange, correlation
iam        = true                  # independent atom model and % difference
iam_basis  = { Li = "cc-pvtz", F = "aug-cc-pvtz" }
anisotropy = true                  # j2 moment, for partially aligned molecules

[detector]                         # optional: one fixed orientation
wavelength  = 1.0
orientation = [0, 90, 0]
# grid = "polar"                   # the rings of grid.f90 (theta_max, rings, azimuths, photon_kev),
                                   # or grid = "qgrid.dat" for your own q vectors (bohr^-1)

# [alignment]                      # optional: partially aligned molecules (laser alignment)
# axis = [0, 1, 0]                 # laboratory alignment axis
# p2   = 0.5                       # <P2(cos beta)> of the molecular z axis
# grid = "polar"

[output]
prefix = "results/lif"
```

`thistle run lif.toml` then writes, in `results/`:

| File | Contents |
|---|---|
| `lif_scattering.csv` / `.png` | total, elastic, inelastic (and parts) |
| `lif_iam.csv` | the same for free atoms |
| `lif_vs_iam.csv` / `.png` | 100 (I_molecule − I_IAM) / I_IAM |
| `lif_anisotropy.csv` / `.png` | the j2 moment S2(q) |
| `lif_detector.npz` / `.png` | the detector image (`.dat` with columns qx qy qz \|q\| ... for a polar grid or a q-grid file) |
| `lif_aligned.dat` / `.png` | the partially aligned ensemble on the detector grid |
| `lif_detector_iam.dat`, `lif_detector_vs_iam.dat` / `.png` (and `lif_aligned_...`) | with `iam = true`: the IAM on the same grid and the % difference map |
| `lif_electron_*.csv` / `.png` | the same for electrons (`probe = "electron"` or `"both"`), plus `lif_electron_sM` |

| `lif_plots.py` | a short script that redraws every plot from these files: edit it to change them |
| `lif_run.json` | the settings, the time of each step and the list of files |

The first line of every CSV records where the numbers came from (file, state, date, version).
At the end of a run THISTLE prints the time of each step and every file it wrote, with its folder:

```
Finished in 57.0 s.
  read the wavefunction                     0.8 s
  X-ray scattering, rotational average      5.3 s
  ...
  plots (15 figures)                        6.0 s

Files written in /home/me/lif/results/
  lif_scattering.csv     X-ray scattering: total, elastic, inelastic (201 q points)
  lif_scattering.png     plot
  lif_vs_iam.png         plot: molecule and IAM, and the % difference
  ...
  lif_plots.py           redraws every plot from these files: python lif_plots.py
```

### Plots

Every calculation has a plotter in `thistle.plotting`, and `thistle plot` uses them on files:

```bash
thistle plot results/lif_scattering.csv                           # any CSV / .dat / .npz THISTLE wrote
thistle plot results/lif_scattering.csv --iam results/lif_iam.csv # molecule vs IAM, % difference below
thistle plot results/lif_detector.dat --iam results/lif_detector_iam.dat   # three panels
thistle plot state1_scattering.csv state2_scattering.csv --curve elastic   # overlay of several runs
thistle plot results/lif_scattering.csv --curves total elastic --log -o fig.pdf
```

```python
import thistle.plotting as tp
res = tp.read("results/lif_scattering.csv")      # back as the object that wrote it
fig = tp.scattering(res, curves=["total", "elastic"], path="fig.png")
fig.axes[0].set_xlim(0, 3); fig.savefig("fig.pdf")   # a matplotlib Figure: change anything
tp.versus_iam(res, tp.read("results/lif_iam.csv"))
tp.sm(...), tp.anisotropy(...), tp.detector(...), tp.detector_versus_iam(...), tp.overlay([...], "total")
```

Each quantity keeps its colour in every plot (total blue, elastic orange, inelastic green,
Hartree yellow, exchange pink, correlation green; the IAM is a dashed grey line), from a palette
checked for colour-blind readers; the footer of each figure says which file, state and probe it shows.

## 4. On a cluster

Add (or edit) the `[cluster]` section of the input file:

```toml
[cluster]
time   = "04:00:00"
cpus   = 32
memory = "32G"
setup  = ["module load gcc", "source ~/thistle-env/bin/activate"]
```

```bash
thistle submit lif.toml          # writes lif.slurm, so you can look at it first
thistle submit lif.toml --now    # writes it and sends it to the queue (sbatch)
```

No Slurm knowledge is needed beyond these four lines. The job runs `thistle run lif.toml` in the
folder of the input file, on all the cores it asked for.

## 5. From Python

```python
import thistle as th

lif = th.load("lif.out", state=2)        # or th.load("scf.molden", hf=True), th.load(my_pyscf_ccsd)
print(lif)                               # what was read

res = lif.scattering(qmax=4)             # q in 1/Angstrom; parts="all" for the decomposition
res.plot("lif.png"); res.save("lif.csv")
res["total"], res.q                      # plain numpy arrays

th.compare(res, lif.iam()).plot()        # % difference from the independent atom model
lif.anisotropy().plot()                  # j2 moment about the molecular z axis
lif.detector(wavelength=1.0, orientation=(0, 90, 0)).plot("image.png")
pat = lif.detector(grid="polar", energy_kev=10, theta_max=60, rings=49, azimuths=72)  # grid.f90 grid
pat.save("pattern.dat"); pat["total"], pat.qvec          # values per point, lab q vectors (bohr^-1)
lif.detector(grid="qgrid.dat", orientation=(0, 90, 0))   # any q-grid file
al = lif.aligned_ensemble(axis=(0, 1, 0), p2=0.5, grid="polar", energy_kev=10)
al.plot("aligned.png")
ref = lif.aligned_ensemble(axis=(0, 1, 0), p2=0.5, grid="polar", energy_kev=10, model="iam")  # same grid
th.compare(al, ref, "total").plot("aligned_vs_iam.png")      # % difference map (also for detector())
th.write_qgrid("qgrid.dat", th.qcoord3(10.0, 60.0, 49, 73))   # exactly grid.f90's file

# electrons: the same calls with probe="electron"
e = lif.scattering(qmax=10, probe="electron", energy_kev=3000)     # dsigma/dOmega in A^2/sr
th.sm(e, lif.iam(qmax=10, probe="electron")).plot("sm.png")        # sM(s): molecule and IAM
lif.detector(probe="electron", energy_kev=3000, orientation=(0, 90, 0)).plot("ued.png")
```

From the command line: `thistle run lif.out --orientation 0 90 0 --grid polar --photon-kev 10`
(fixed orientation) or `thistle run lif.out --align 0 1 0 --p2 0.5` (aligned ensemble on the
polar grid); `--theta-max`, `--rings`, `--azimuths` set the grid, `--grid qgrid.dat` reads one.

### Coupled cluster on larger molecules

A CCSD two-particle density is an n⁴ array: 120 basis functions already need more than 6 GB,
300 would need hundreds of GB. THISTLE therefore keeps the CCSD density as a sum of
7 n_occ² + 4 products of one-particle matrices (memory ≈ 112 (n_occ·n_AO)² bytes: 1 GB at
n_occ·n_AO ≈ 3000, 10 GB at 9500, 100 GB at 30 000) and contracts it directly. This happens
automatically for `th.load(mycc)` beyond 60 basis functions; the numbers are identical to 10⁻¹³.

```python
from pyscf import gto, scf, cc
mol = gto.M(atom="...", basis="aug-cc-pvtz")
mycc = cc.CCSD(scf.RHF(mol).run()).run()          # lambda equations are solved when needed
mol_th = th.load(mycc)                             # factored density (60+ functions)
import pyxscat as xs
wf = xs.interfaces.pyscf.from_cc(mycc, factored=True)   # or False: PySCF's four-index array
```

The limit is then usually PySCF's own CCSD (it stores the virtual–virtual integrals: ~300 GB for
640 functions, ~900 GB for 840), not the scattering. Rough guide for a 10-heavy-atom molecule
(30 occupied orbitals): aug-cc-pVTZ needs ~40 GB of factors and ~2 h of scattering on 64 cores;
cc-pVQZ ~70 GB and ~5 h; the time grows as n_occ² n_AO⁴. Frozen-core CCSD is supported.

CASSCF, CASCI and FCI densities (from PySCF, Molpro or Molcas) are handled the same way beyond
60 basis functions: the core and active blocks are kept as products of one-particle matrices
(at most n_act(n_act+1)/2 + 6 of them), so the basis set no longer limits memory: benzene
CAS(6,6)/cc-pVTZ runs in 0.33 GB instead of needing a 65 GB array. `factored=False` in
`from_casscf`, `from_fci`, `molpro.from_multi` or `molcas.from_rasscf` keeps the array.

A full example with PySCF is in `examples/quickstart.py`. Everything underneath (cutoffs, raw
kernels, your own density matrices) is in the `pyxscat` package, documented in the README.
