# Testing THISTLE 0.2.0.dev0

*Authors: Andrés Moreno Carrascosa, Mats Simmermacher, Patrick Yuheng Wang, Eirik Moen Liane, Joseph Cooper, Nikola Zotev, Thomas Northey and Adam Kirrander.*

Thank you for trying the code. This is a development version: the numbers have been validated
extensively (see `docs/THEORY.md`, section 5). It has been built and checked with gfortran 9-14
and with Intel ifx 2026.1 (Linux x86_64, 64 threads); other compilers and platforms are untested. Problems with your
compiler, platform or input files are exactly what we want to hear about.

## 1. Install

Python 3.10 or newer and a Fortran compiler with OpenMP.

```bash
tar xzf thistle-0.2.0.dev0.tar.gz
cd thistle-0.2.0.dev0
python -m venv ~/thistle-env && source ~/thistle-env/bin/activate     # recommended

pip install ".[test]"                      # gfortran (default)
# or, with Intel oneAPI:
source /opt/intel/oneapi/setvars.sh        # or: module load intel
FC=ifx pip install ".[test]"               # FC=ifx EXTRA=-xHost ... optimises for this CPU
```

`pip` compiles the Fortran library during the install (look for `built ... with gfortran/ifx`
in `pip install -v` output). `[test]` also installs PySCF, SciPy, xraydb and pytest.

No Intel compiler on the cluster and no root access? Install ifx for your user only with conda
(Miniforge, if you have no conda yet):

```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3 && source $HOME/miniforge3/bin/activate
conda create -n thistle-ifx -c conda-forge python=3.12 ifx_linux-64
conda activate thistle-ifx
ifx --version
FC=ifx pip install -v ".[test]"
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH     # if importing pyxscat cannot find libiomp5
```

With ifx, keep the oneAPI environment loaded whenever you run THISTLE (the library needs the
Intel OpenMP runtime, libiomp5). If Python reports `OMP: Error #15 ... already initialized`,
two OpenMP runtimes were loaded; `export KMP_DUPLICATE_LIB_OK=TRUE` is a workaround, and please
report it.

## 2. Quick check (seconds)

```bash
thistle selftest                           # build info + every kernel against sum rules and an independent route
thistle selftest my_calculation.out        # the same sum rules for your own file (Python side and Fortran side)
```

The first lines say which compiler and flags built the library that Python actually loaded.
If anything says FAIL, the full output of these two commands is the most useful bug report.

To switch compiler, just reinstall (`FC=ifx pip install .` or `pip install .`): the build always
starts from scratch, so a library from the other compiler or an older version is never reused.

## 3. Run the test suite

From the unpacked folder:

```bash
python -m pytest -q                        # about 5 minutes, 217 tests, 5 skipped
PYXSCAT_SLOW=1 python -m pytest -q         # adds the large-basis tests (about 10 minutes)
OMP_NUM_THREADS=1 python -m pytest -q tests/test_elastic_pairs.py    # thread-count check
```

Everything should pass. Several tests compare two independent routes to 10⁻⁹–10⁻¹², so a
failure by a small margin with ifx (or on another CPU) is informative: please send the full
pytest output.

### Comparing builds or machines, and parallel scaling

```bash
python benchmarks/cluster_check.py reference ifx.npz      # every kernel and probe, ~20 s
# ... reinstall with the other compiler (or on another machine), then:
python benchmarks/cluster_check.py reference gnu.npz
python benchmarks/cluster_check.py compare ifx.npz gnu.npz  # should agree to < 1e-8
python benchmarks/cluster_check.py scaling                # timings for 1, 2, 4, ... threads
```

## 4. Try it on your own calculations

```bash
thistle run my_molpro.out --state 2 --iam        # Molpro MULTI output (config,det) + its Molden file
thistle run my_molpro.pun --state 2              # or just the punch file + Molden file (no output)
thistle run my_molcas.log --all-parts            # OpenMolcas RASSCF output + .rasscf.molden
thistle run scf.molden --hf                      # closed-shell Hartree-Fock Molden file
thistle run my_molpro.out --electrons --iam      # electron diffraction, with sM(s)
thistle check my.molden --scf-energy -76.0266    # is the file read correctly?
thistle run my_molpro.out --orientation 0 90 0 --grid polar --photon-kev 10   # one orientation, grid.f90 grid
thistle run my_molpro.out --align 0 1 0 --p2 0.5 --grid qgrid.dat             # aligned ensemble on your grid
```

or in Python (see `docs/QUICKSTART.md` and `examples/quickstart.py`):

```python
import thistle as th
mol = th.load("my_molpro.out", state=2)
mol.scattering(qmax=4).plot("xray.png")
```

Each run writes CSV tables and PNG plots, a `*_plots.py` script that redraws the plots from the
files, and a `*_run.json` log; it ends by listing every file with the time of each step. The
first line of each CSV records the input file, state, settings and version. `thistle plot FILE`
redraws any result file (`--iam` for the comparison with the IAM).

## 5. What to report

Please send:

* the output of `thistle selftest` (it prints the compiler and options the library was built
  with, the thread count and the versions), and how you installed;
* the full output of any failing test or command;
* for your own calculations: the input files (or a small case that shows the problem) and what
  you expected, e.g. a curve from the earlier code;
* timings that look wrong, with `OMP_NUM_THREADS` and the number of cores.

CCSD from PySCF beyond 60 basis functions uses a factored two-particle density (no n⁴ array;
memory ≈ 112 (n_occ n_AO)² bytes), the same numbers to 10⁻¹³; `from_cc(mycc, factored=False)`
forces the array for comparison.

Known limits: Molden files stop at g functions (h works through PySCF); Molpro MULTI needs
`config,det` and no point-group symmetry; Molcas needs the RASSCF CI printout. See the end of
`docs/THEORY.md`.
