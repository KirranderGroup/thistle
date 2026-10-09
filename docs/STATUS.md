# Status of the refactoring

## Milestone 1: package-independent kernel (4 October)

**One kernel, called in memory.** Rotationally averaged scattering is computed from a basis set
and density matrices, called from Python through a C interface. There are no intermediate files.

**Density providers.**
- Explicit 2RDM (any correlated method).
- Built on the fly from one-particle matrices: elastic, transition and coherent cross terms,
  exchange, and Hartree–Fock mean field.

**Hermite coefficients by recurrence.** Valid for any angular momentum; the overall limit (g
functions) comes from the Bessel and h = 0 tables.

**PySCF adapter.** Basis taken from the `Mole` object, spherical basis sets included. Adapters
for SCF, CCSD, CISD, FCI and CASSCF.

## Milestone 2: other packages (night of 4 October)

**Molden reader (`interfaces.molden`).** Pure numpy, no PySCF needed.
- Handles Angstrom and bohr geometries, `sp` shells, Fortran `D` exponents, and the
  `[5D]`/`[7F]`/`[9G]` flags with spherical d, f and g functions.
- Molden component order and signs are handled; general contractions written as separate
  shells are merged into shared primitive shells.
- Cartesian normalisation differs between programs. With the default (`cart_norm="auto"`) the
  reader picks the convention that makes the MOs orthonormal.
- Checked by round trip through PySCF (cc-pVTZ and cc-pVQZ spherical up to g, 6-31G*
  Cartesian, def2-SVP). Density on a grid agrees to 1e-11 and the scattering to 1e-9.
- Checked on the three Molden files in `Scattering_eirik` (OpenMolcas CASSCF, a C5H6 file, a
  PySCF UHF file): MOs orthonormal to the printed precision, electron count exact.

**CI vectors to RDMs (`pyxscat.ci`).**
- Spin-summed 1RDM and 2RDM from any determinant expansion; equal to PySCF FCI to 1e-13.
- GUGA CSF → determinant conversion (Molcas step strings `2ud0`), ported from the old reader
  with the reordering phase made explicit. Every CSF tested is a normalised S² eigenfunction
  for every M_S, and CSFs of one configuration are orthogonal.
- Core embedding of active-space RDMs; equal to PySCF CASCI densities to 1e-10.

**OpenMolcas adapter (`interfaces.molcas`).** Reads the RASSCF output (orbital spaces, spin,
CI printout of every root, natural occupations) and the `.rasscf.molden` file. Validated on the
stored CASSCF(8,8) run of C7H8 (1764 CSFs):
- The 1RDM is diagonal in the printed natural orbitals, with Molcas' occupations to 5e-6.
- The converted vector is PySCF's CASCI ground state in the same orbitals (overlap 0.999996;
  the deficit comes from the five printed decimals).

**Command line (`pyxscat`).** Molden alone (elastic, `hf_total`) or Molden plus Molcas output
(all kinds, including `decompose`).

**Tests.** 66 run by default (about 2 min). Five more on the large ANO basis run with
`PYXSCAT_SLOW=1`.

## Milestone 3: aligned molecules (night of 4 October)

**New approach.** The plane wave factorises, exp(iq·(r1−r2)) = exp(iq·r1) exp(−iq·r2), so
fixed-orientation scattering never needs the four-index quartet loop. One Fortran routine
(`xs_pairft.f90`) computes the analytic Fourier transforms of all basis-function products,
F_μν(q), at any list of q vectors. Everything else is contraction with density matrices:
- elastic: |Σ D F|²
- two-electron: Σ Γ F F*, a single matrix product
- exchange and mean field: from the 1RDM, at O(n³) per q

There is no Bessel expansion, so the close-centre problem does not arise here.

**Python (`pyxscat.aligned`).**
- `pair_transform`, `form_factor`, `elastic` (also transition/mixed terms), `two_electron`,
  `exchange`, `mean_field`, `scattering`, `decompose`.
- Geometry helpers: `ewald_q` (Ewald sphere, beam along +z), `detector_q` (q-space detector
  image), `rotation_matrix` (zyz Euler angles) and `rotate_q` (lab → molecular frame).
- Command line: `--qvec file [--orientation α β γ]`.

**Validation (14 tests).**
- Pair transforms equal PySCF's `ft_aopair` (complex conjugate, since PySCF uses exp(−iq·r)) to
  1e-12, up to g functions. Spherical basis sets are checked through the form factor.
- Averaging the aligned elastic, exchange and CCSD two-electron terms over a 131-point Lebedev
  grid reproduces the rotationally averaged kernel to 1e-8. These are two independent code
  paths.
- Exact q = 0 limits; the explicit HF 2RDM equals the on-the-fly mean field; the decomposition
  sums to the total; Friedel symmetry I(q) = I(−q).
- Rotating q is equivalent to rotating the molecule and re-running SCF.

**Speed.** A 101×101 detector image of the C7H8 CASSCF wavefunction (ANO-L, full 2RDM, all
terms of the decomposition) takes 7 s on 2 cores. Three rotationally averaged q points of the
same system take about 3 minutes.

**Not yet compared with Eirik's aligned code.** I have not run it, so I can't yet say which bug
it had. Reading it, the structure is the same (Hermite expansion at the pair centre, phase
exp(iq·P)). It loops over all primitive quartets with no symmetry, with the 2RDM entering
through the old Z1/Z2 factors scaled by 1/8, and the powers (iq)^t are built with a
single-precision `(0,1.00)` constant. The q grids (`QCOORD3`, mode 4) are hard-coded. An input
and output where it went wrong would let me pin it down.

## Molpro (5 October)

**Molpro MULTI adapter (`interfaces.molpro`).** Reads the CI vectors of all states from a MULTI
run with `config,det` and `gprint,civector`, with orbitals from the Molden file written after
MULTI. Also gives transition 1RDMs between states (`transition_dm1`, `ci.make_transition_rdm1`)
for inelastic and coherent mixed terms. Command line: `--molpro-out file --root k`.

**Determinant phase.** Molpro orders the creation operators orbital by orbital (1α 1β 2α 2β …),
while pyxscat uses all-α-first. Without the phase correction (`ci.interleaved_phase`) the state
energies are off by 2–6 mEh; with it they agree to 1e-8.

**Validation (`tests/test_molpro.py`).** Test case: LiF from `Dropbox/PyXSCAT/src/EQ-Minimal`, a
state-averaged CASSCF(6,4) over four singlets.
- The energies of all four states, rebuilt from our RDMs with PySCF integrals in the Molden
  orbitals, match Molpro to 2e-7.
- The state-averaged natural occupations are 1.75/1.75/1.75/0.75.
- The ground-state dipole and the transition dipoles ⟨1|x|3⟩, ⟨1|y|2⟩ and ⟨1|z|4⟩ match Molpro to
  1e-7, signs included.
- Fixed-orientation scattering has exact q = 0 limits, and the transition form factor vanishes
  at q = 0.

**Still to do for Molpro.**
- The CSF printout (`/\` notation, the default without `config,det`) is not supported; Molpro's
  CSF coupling convention would need checking.
- Point-group symmetry is not supported.
- MRCI is not covered. `MRCI_CI.dat` holds reference CSFs only; MRCI densities would have to come
  from Molpro's density matrices instead.

## Independent atom model and more Molpro cases (5 October)

**IAM (`pyxscat.iam`).** Free spherical atoms at the nuclear positions, rotationally averaged
(Debye formula) or at fixed orientation (`aligned=True`). It returns elastic, total, inelastic,
atomic and interference parts.
- **Ab initio atoms:** free-atom UHF (or CCSD for closed shells) in any PySCF basis, run through
  the same kernels as the molecule. Elastic and total.
  - The form factor comes from the spherically averaged density, by a Lebedev average of the
    fixed-q transform; the atomic total from the rotationally averaged two-electron kernel.
- **Tables:** Waasmaier–Kirfel form factors through `xraydb`. Elastic only.
- **Validation (11 tests).**
  - HF/aug-cc-pVQZ form factors of H, C, O, F, Ne and S agree with Waasmaier–Kirfel within 1 %
    of Z, and typically 0.1 %.
  - A closed-shell atom reproduces its own IAM to 1e-8.
  - The aligned IAM averaged over a Lebedev grid equals the Debye IAM to 1e-9.
  - Water CCSD/aug-cc-pVTZ tends to the IAM at high q (within 0.5 % at q = 14 bohr⁻¹) and differs
    by several per cent at q ≈ 1.

**First comparisons with the IAM** (percent difference, total scattering, q = 1 bohr⁻¹):
- LiF ground state (Molpro SA-CASSCF): +16 %. It is ionic (Li⁺F⁻), which a neutral-atom model
  cannot describe.
- LiF state 4: −2 %, with the dipole reversed.
- H₂O CCSD: +6 %.

**More Molpro.**
- The older printout format (`CI vector … TOTAL ENERGIES`) is read, as in the SF6, SF6⁺, SF6²⁺
  and He runs in Dropbox.
- These runs used Molpro's default print threshold, so only 3–9 determinants of each vector were
  printed. The printed weight is 0.974 for SF6 CAS(6,6), 0.999 for SF6⁺ and 0.996 for SF6²⁺.
- Densities built from such printouts give energies 3–52 mEh above Molpro's. **Any earlier
  scattering computed from these printed vectors used truncated wavefunctions.** The reader now
  warns, and `gthresh,printci=0` avoids the problem.
- The printed coefficients themselves are exact: their overlap with the PySCF CASCI eigenvector
  equals the printed weight, for singlet, doublet and triplet. So the phase convention holds for
  open shells and for the older format too.
- `check_molden_consistency` compares the state-averaged 1RDM with the Molden occupations. It
  flags the MateoBottaro He run, where the CI vector and the Molden orbitals are in a different
  order (two near-degenerate active orbitals appear swapped).

## Punch files, orbital order, Molden conventions (5 October, night)

**Molpro punch files** (`molpro.read_punch`, `from_multi(..., punch=...)`, command line
`--molpro-punch`).
- Older Molpro versions write the MULTI CI vector to the `.pun` file with seven decimals and a
  low threshold.
- For the SF6, SF6⁺ and SF6²⁺ runs the punch files hold 120–299 determinants (weight 1.0000000).
  The state energies rebuilt from them match Molpro to 3e-9 Eh, compared with 3–52 mEh off from
  the output printout.
- The closed-shell count is not in the punch file. It comes from the output, from `nclosed`, or
  from the electron count of the Molden file.

**Molpro orbital order.** Molpro may write natural orbitals to the Molden file out of its
internal order: the He run lists `Sym= 1.1 2.1 6.1 3.1 4.1 5.1 …`. The CI vector refers to the
internal order, so the Molden reader now restores it from the `n.1` labels (`reorder="auto"`).
With this fix the two He state energies match Molpro to 5e-9 Eh; read in file order they were 33
and 43 mEh off. This was the cause of the "orbitals out of step" flag.

**Molden conventions check (`molden.scf_energy`).**
- Rebuilds the SCF energy from a Molden file with exact integrals. It equals the program's SCF
  energy only if the basis, normalisation, spherical → Cartesian conversion, component order and
  signs are all read correctly.
- Passes to 1e-12 on PySCF files (spherical up to g, Cartesian).
- `benchmarks/molden_conventions/` has water/cc-pVTZ HF inputs for OpenMolcas, BAGEL and Molpro,
  plus `check.py`, to confirm each program's conventions.
- None of the Molden files in Dropbox or the code folders uses spherical d or f functions: Molpro
  writes Cartesian functions, and the one Molcas file has no d functions. So the spherical path is
  so far verified only on PySCF output.

## j2: second Legendre moment (5 October, night)

**What the old code computes.** Eirik's j2 modes evaluate
Σ Γ ⟨j₂(q r₁₂) Y₂₀(r̂₁₂)⟩, the l = 2 Legendre component of the scattering about the molecular z
axis. They use their own expansion of the derivatives of j₂(qh)(3h_z² − h²)/h² (`rrdj2a/b`,
`BesselDeriv2j`, `bessels2rr`) and a j2 table of coincident-centre limits. Only the elastic mode
(8) was active; the total mode (7) is commented out in the driver.

**New implementation.** The j2 kernel follows from the j0 one by Hobson's theorem: the operator
3∂²_z − ∇² acting on j₀(qr) gives (3z² − r²) times a radial factor, so
j₂(qr) P₂(cos θ) = [3 ∂²_z j₀(qr) + q² j₀(qr)] / 2q².
- So the existing j0 quartet kernel only carries two extra z-derivatives: Hermite index v + 2
  for distant pairs, and the h = 0 limits for coincident ones.
- The h = 0 limits are now a closed formula for any order, replacing the hard-coded table
  (limited to L + M + N ≤ 16).
- All density providers work: total, elastic, exchange, mean field, decomposition, transition
  terms.
- Python: `l=2` in `scattering`, `decompose`, `elastic`, `two_electron`, `exchange`,
  `mean_field`, and in `iam.scattering`. Command line: `--legendre 2`.

**Normalisation.** pyxscat uses Legendre normalisation, S₂(q) = Σ Γ ⟨j₂ P₂⟩. For molecules
aligned along a lab axis,
I(q, θ_q) = I₀(q) − 5 ⟨P₂(cos β)⟩ P₂(cos θ_q) S₂(q) + …. Eirik's (Y₂₀) values are
√(5/4π) S₂.

**Validation (9 tests, `tests/test_j2.py`).**
- **Projection:** the l = 2 kernel equals −⟨P₂(cos θ_q) I(q n̂)⟩ computed from the
  fixed-orientation code on a Lebedev grid. This holds for elastic, exchange and CCSD
  two-electron terms to 5e-8 (water, 6-31G*, spherical), across two independent code paths.
- **Eirik's code:** his elastic j2 (`elastic_integration_j2`), driven directly from Fortran,
  equals √(5/4π) × ours to 1e-6 relative.
- **IAM j2:** S₂^IAM = Σ_{A≠B} f_A f_B j₂(qR_AB) P₂(cos θ_AB); atomic and inelastic parts vanish
  for spherical atoms. Two Ne atoms 14 bohr apart on an oblique axis reproduce it from the
  molecular code to 1e-6, both elastic and total. The analytic IAM j2 equals the P₂ projection of
  the aligned IAM to 1e-10.
- **Limits:** a spherical atom gives S₂ = 0; S₂ ∝ q² at small q; the decomposition sums.

**Small q.** The formula divides by q². Below q ≈ 0.3 bohr⁻¹ the absolute error is about 1e-8,
growing roughly as 1/q²; values below q = 1e-5 are set to zero.

## Simple interface and packaging (5 October, night)

`thistle` is a thin layer over `pyxscat` for people who do not want to touch density matrices
or q units:

* `th.load(x)` takes a PySCF object, a Molden file, an OpenMolcas output or a Molpro output or
  punch file, recognises it by content, and finds the companion Molden (and `.pun`) file by
  itself. Several Molden candidates in one folder is an error, not a guess.
* `Molecule.scattering / iam / anisotropy / detector` return `Result` / `DetectorImage`
  objects with `.plot()` and `.save()` (CSV whose first line records the provenance). q in
  1/Å by default, `unit="1/bohr"` available.
* Orbitals-only input (Molden) gives elastic only. Total and inelastic only with the explicit
  option `--hf` / `--scf` (`hf = true`, `th.load(..., hf=True)`), which is refused unless the
  natural occupations (Löwdin) are 0 or 2. PySCF RHF objects are recognised by themselves.
* `thistle new | run | check | submit`: commented TOML template; `run` also takes a calculation
  file directly with `--state --qmax --iam --all-parts ...`; `check` verifies orthonormality,
  SCF energy and I(0) = N²; `submit` writes a Slurm script (cd, setup lines, OMP threads) and
  sends it with `--now`.
* `pip install .` (or the sdist `thistle-0.2.0.dev0.tar.gz`) compiles the Fortran through
  `make` from `setup.py`; the distribution is now called `thistle` and installs both the
  `thistle` and `pyxscat` commands. Tested from the sdist in a clean virtual environment.
* `tests/test_thistle.py` (17 tests). Full suite: 133 passed, 5 skipped.

To decide: default q unit (1/Å) and the distribution name (`thistle`; the import name of the
kernels stays `pyxscat`).

## Comparison with thistle-dev (7 October)

Patrick Wang's `KirranderGroup/thistle-dev` (private, beta "pythistle 1.0.0b1") builds a
distributable from the same legacy kernels (Eirik's `Integrals.f90` MO path + Patrick's dense-AO
`TSj0groupsfast`), with a Julia engine and an HDF5 input "deck". On its five test decks our
kernels reproduce its reference curves to 4e-13 when run with its `cutoffcentre = 0.01`; our exact
default differs from them by 1.6e-5 (the coincident-centre approximation). Its open j₂ question
(Eirik's vs its Julia engine's) is settled: Eirik's is right. Its v1 scope is X-ray total and
elastic; ours is a superset except the Julia engine, the HDF5 deck, and the generic PySCF
`make_rdm12` route (DMRG-SCF, SHCI-SCF). Details and scripts: `benchmarks/thistle_dev_crosscheck/`.
To ease the merge, `th.load("deck.h5")` now reads thistle-dev decks natively (orbital and dense-AO
forms; `pyxscat.interfaces.deck`; `tests/test_deck.py`). Our code is on the branch
`amc/thistle-core` of thistle-dev, in its own directory `thistle_core/`.

## CCSD without the four-index density (6 October)

Checked first as asked: the CCSD route took PySCF's full 2RDM, transformed it to the AO basis
and read it in the kernel. Measured: water/cc-pVTZ 0.5 GB peak; benzene/cc-pVDZ (120 Cartesian
functions) killed at 6 GB (the AO array alone is 1.35 GB, PySCF's `make_rdm2` and the transform
hold several copies). So anything beyond ~150 functions was out of reach.

* `pyxscat.ccsd.factored_rdm2`: the CCSD density as `FactoredDM2`, Σ_k A_k⊗B_k in Coulomb and
  exchange pairings, 7 n_occ² + 4 terms, from t₁, t₂, λ₁, λ₂ and the small PySCF intermediates;
  the n_vir⁴ and n_occ n_vir³ blocks are never formed. Every block of the MO density agrees
  with `make_rdm2` to 10⁻¹⁵; frozen core handled.
* Kernel mode 3 contracts the factored form: Coulomb terms in the canonical loop; exchange terms
  in a new pass ordered by the (μ,λ) shell pair with batches of shells, eight matrix products
  over k per batch (dgemm from SciPy through a function pointer, matmul fallback). Validated
  against mode 1 with random non-symmetric factors to 10⁻¹³ (l = 0, 2); all routes (rotational,
  j₂, fixed orientation, electrons) agree with the dense CCSD route to 10⁻¹³.
* The old `Total_scattering_optimized_ccsd.f90` had a factored cumulant too, but its
  `build_tau_ij` stops at the "dominant terms" (its own TODO) and was never validated; this
  one is complete.
* Timings (2 cores): water/cc-pVTZ 4.3 s (dense 2.2 s); benzene/cc-pVDZ 120 s in 2.4 GB (the
  integrals take 16 s; the exchange contraction dominates and is compute-bound in dgemm, so it
  scales with cores). Estimate for 500 functions / 40 occupied: ~45 GB of factors, tens of
  minutes on a 64-core node.
* Default: factored beyond 60 basis functions (`from_cc(mycc, factored=...)`).
* Memory rule: factors ≈ 112 (n_occ n_AO)² bytes. Estimates for dimethylpyrazine (29 occupied)
  on the 2 TB node: aug-cc-pVTZ (640 Cartesian functions) 39 GB of factors, ~1.7 h of
  scattering on 64 cores, PySCF CCSD ~300 GB; cc-pVQZ (840) 66 GB, ~5 h, CCSD ~870 GB (the
  practical maximum); aug-cc-pVQZ (1280) would need density-fitted CCSD (conventional vvvv
  storage 2.4–5 TB). Suggested headline CCSD calculation for the paper: aug-cc-pVTZ.
* Other programs' CCSD amplitudes (Molpro, ...) would need a reader to the same conventions
  (closed-shell t₁, t₂, λ₁, λ₂ in canonical orbitals); not written for lack of files.
  `tests/test_ccsd_factored.py` (9 tests); full suite 217 passed, 5 skipped.

## Plotters, run report (6 October)

* `thistle.plotting`: a plotter per calculation (`scattering`, `versus_iam` with the % difference
  below, `sm`, `anisotropy`, `detector` for images and patterns, `detector_versus_iam` in three
  panels, `overlay` of several runs, `difference`), `read` for every file THISTLE writes (CSV,
  .dat, .npz, all round-trip to 10⁻¹¹) and `plot_file`; `thistle plot FILE [--iam FILE]`. Fixed
  colour per quantity (validated colour-blind-safe palette), provenance footer, identical curves
  drawn once ("elastic = hartree"), polar patterns triangulated along their rings (a Delaunay
  triangulation of co-circular points changed with the last digit), robust colour range for
  electron patterns (1/q⁴ at the centre).
* `thistle run` prints each step with its start time and duration, and at the end the time of
  every step and every file written (with its folder and a line on its contents); it also writes
  `<prefix>_plots.py` (redraws every plot from the files, to edit) and `<prefix>_run.json`.
* Electrons on the polar grid: the default largest angle reaches the run's q max (60° made no
  sense for fast electrons: q up to 1000 Å⁻¹).
* `tests/test_plotting.py` (4 tests); full suite 208 passed, 5 skipped.

## Molpro without the output file (6 October)

A punch file and the Molden file are enough: `thistle run calc.pun` (or `th.load("calc.pun")`)
gives the same RDMs as the output to 10⁻¹² (tested). New: when the named output is missing,
the punch file of the same name is used (with a warning); a Molden file with `punch=` /
`--punch` / `[input] punch = ...` is read as the MULTI calculation; a Molden file alone with a
punch file next to it warns that only the elastic scattering of the (state-averaged) density
follows from it.

## 2RDMs from long determinant expansions (6 October)

The builder used so far (`ci.make_rdm12`) was the string-driven one: the CI vector embedded in
the complete active space and Γ = ⟨E_qp Ψ|E_rs Ψ⟩ − δγ in one matrix product. It is exact and
fast for CAS vectors (CAS(12,12), 8.5·10⁵ determinants, 1.6 s), but its cost and memory follow
N_FCI, not the number of determinants given: a truncated or selected expansion in a large
orbital space (MCCI, a printed threshold) was impossible (CAS(16,16): MemoryError for 1000
determinants). Now:

* `method="bitwise"`: the determinant-pair algorithm of `twordm_jeremy.f90`, in Fortran with
  OpenMP (`src/fortran/xs_rdm.f90`): popcnt(ieor) screening, Slater–Condon phases from the
  bit strings, all spin cases, the Hermitian partner added at the end; up to 64 orbitals.
  Built with the hardware popcnt on x86-64 (`-mpopcnt`, ifx `-march=x86-64-v2`, this file
  only): without it the screening is 7× slower.
* `method="sparse"`: the same intermediates restricted to the determinants reached, linear in
  N_det, in memory-bounded passes (NumPy/SciPy).
* `method="dense"`: as before, now in blocks of β strings (memory ~ the CI vector).
* `"auto"` (default everywhere) picks the cheapest by a measured cost model.

Timings (2 cores here): 27 orbitals, 12 electrons: CISD-like 22 429 determinants, bitwise
0.4 s, sparse 1.4 s; 10⁵ random selected determinants, bitwise 2.9 s, sparse 9.0 s;
4·10⁵: bitwise 264 s (pair loop), sparse 55 s. With 64 cores the pair loop is ≈ 30× faster.
Validated against PySCF FCI and selected CI (27 orbitals) to 10⁻¹⁵, the three against each
other up to 64 orbitals; `thistle selftest` now checks the Fortran builder too.
`tests/test_rdm_builders.py` (22 tests); full suite 203 passed, 5 skipped.

## Detector grids of grid.f90; partially aligned ensembles (6 October)

* `pyxscat.grids`: the polar grid of `grid.f90`/`qgrid3d.f90` (`qcoord3`), point for point
  (`write_qgrid` gives a byte-identical `qgrid.dat`), q-grid files, photon and electron k₀. The
  original sources are kept in `benchmarks/legacy_grid/`. Note: `grid.f90`'s n2 = 73 means
  72 points per ring (plus q = 0).
* `Molecule.detector(grid="polar" | file | array, ...)` for one fixed orientation, returning a
  `DetectorPattern` (text output with columns qx qy qz |q| curves, as qgrid.dat).
* `Molecule.aligned_ensemble(axis, p2, grid=...)`: I = I₀ + 5⟨P₂⟩P₂(q̂·â)S₂ on any grid (sign as of 8 October), X-rays
  and electrons, I₀ and S₂ computed once per ring.
* Input files: `[detector] grid/theta_max/rings/azimuths/photon_kev`, new `[alignment]
  axis/p2/grid`; `thistle run ... --orientation A B G | --align X Y Z --p2 P --grid polar ...`;
  `pyxscat ... --polar-grid E0 THSMAX N1 N2`.
* Validation: the ensemble formula equals an explicit average of fixed-orientation patterns over
  the distribution 1 + 5⟨P₂⟩P₂(cos β) to 10⁻¹²–10⁻¹⁰ (X-ray elastic and HF total, electron
  elastic).
* The IAM on the same grids (`model="iam"` in `detector` and `aligned_ensemble`; new
  fixed-orientation electron IAM), `th.compare` for patterns and images (% difference maps), and
  with `iam = true` the input files write `_detector_iam`, `_detector_vs_iam`, `_aligned_iam`,
  `_aligned_vs_iam`. First result: LiF SA-CASSCF at 10 keV on the 49 × 72 grid differs from the
  IAM by up to 25 % (fixed orientation, along the bond near q = 1.6 Å⁻¹, where the IAM has its
  interference minimum) and 21 % for ⟨P₂⟩ = 0.5. `tests/test_grids.py` (14 tests); full suite 181 passed, 5 skipped.
* Authors: Andrés Moreno Carrascosa, Mats Simmermacher, Nikola Zotev, Joseph Cooper, Patrick
  Yuheng Wang, Thomas Northey (pyproject, README, docs, source headers). Licence still open.

## Tester report: total scattering wrong on the cluster (5 October, evening)

Report: LiF (Molpro test case), total(0) = 41.8 instead of 144 (elastic correct), "with ifx and
with gfortran". The same tarball gives 144 here (gfortran 13, 1-32 threads, latest NumPy/SciPy/
PySCF in a fresh environment). Elastic goes through the new pair kernel, total through the
quartet kernel, so the quartet kernel was wrong in that build.

Likely mechanism for "the compiler makes no difference": the Makefile rebuilt by file times
only. Unpacking a tarball over an earlier build gives sources that look older than the existing
objects, and changing FC did not trigger a rebuild, so the second "gfortran" install probably
reused the ifx library (or objects from an older version, e.g. a 20×20 Bessel module next to
a kernel compiled for 32×32, which corrupts exactly the quartet kernel and nothing else).

Changes:
- The Makefile keeps a stamp of compiler, flags and a checksum of all sources and rebuilds
  everything when any of them changes; `pip install` always starts with `make clean`.
- The library reports the compiler and options it was built with (`pyxscat_build_info`).
- `thistle selftest [files]`: build info, then every kernel against sum rules and against the
  fixed-orientation route (13 checks, ~2 s), and for given files the sum rules separately for the
  Python side (tr DS = N, Σ Γ S S = N(N−1)) and the Fortran side.
- The quartet kernel gathers the 2RDM block with explicit loops instead of a vector-subscripted
  section of a pointer array (in case ifx mishandles it in the OpenMP region).
- Regression test: LiF total(0) = elastic(0) = 144 through the interface.

**Diagnosis (with the tester, on the cluster).** The loaded library was ifx 2026.1.1 (conda),
so the "gfortran" run had indeed reused it. Optimisation levels: -O0 and -O1 pass, -O2 and -O3
fail (with fp-model fast or precise), -O3 -fp-model strict passes. Compiling one file at a time
at -O1 with the rest at -O3: only `xs_hermite.f90` cures it. Its recurrence read
`e(t-1, ...)` under `if (t > 0)` inside an OpenMP loop over a private automatic array, and ifx
miscompiled that at -O2/-O3; gfortran 9-14 and a bounds-checked debug build are fine.

**Fix.** The Hermite recurrences (`hermite_table`, and `hermite_pair` used by the pair, transform
and electron-nucleus kernels) now work on an array padded at t = -1 and t = 2 l_max + 2, so
there are no guarded loads at all; `hermite_table` is no longer threaded (it is negligible).
Confirmed on the cluster: ifx 2026.1.1 (conda-forge), Linux x86_64, default flags (-O3
-fp-model precise), 64 threads: `thistle selftest` passes, and so do the LiF sum rules
(I2(0) = Σ Γ S S to 8·10⁻¹¹, one thread vs 64 to 4·10⁻¹⁵).

## Ready for testers; Intel compiler (5 October, evening)

- `make FC=ifx` / `FC=ifx pip install .` (also ifort): `-O3 -fPIC -qopenmp -fp-model precise`
  (value-safe floating point, since the tests compare routes to 1e-12), `-module build`, Intel
  runtime linked statically (`-static-intel`); `EXTRA=` adds flags (e.g. `-xHost`). FFLAGS and
  LDFLAGS from the environment (conda sets them for gfortran) are ignored; on the command line
  they replace the defaults.
- Not run with ifx: no Intel compiler can be installed in this environment. To reduce the risk,
  the sources now compile with `gfortran -std=f2018 -pedantic` without line-length extensions
  (two declarations in the old Bessel module used a variable before declaring it, a GNU
  extension that ifx may reject; two lines were longer than 132 characters).
- Clearer errors: compiler not found (gfortran or oneAPI hint) and library not loadable (Intel
  OpenMP runtime hint).
- `docs/FOR_TESTERS.md`: install (gfortran or ifx), tests, own calculations, what to report.

## h functions (5 October, afternoon)

**What limited us to g.** Only the quartet kernel: the Bessel-derivative tables of the old
module (`rrdj0`, `BesselDeriv`) had fixed 20×20 dimensions, and a quartet needs derivative
orders up to 4 l_max + 2 (+2 for j₂), i.e. 23 for h. Everything written since (pair transforms,
elastic pairs, electron–nucleus, the radial functions) is dimensioned from the basis, and the
Hermite recurrences and the h = 0 closed form were already general.

**Change.** One parameter, `BDIM_MAX = 32` in `xs_bessel_j0.f90` (room up to l = 7).

**Checks** (`tests/test_h_functions.py`, CO with an h shell on C; also CO with cc-pV5Z on C):
- PySCF input, pair transforms, elastic by pairs (l = 0 and 2) and electron–nucleus against the
  Lebedev route: 10⁻¹¹–10⁻¹⁴.
- Quartet kernel: elastic against the pair route 2·10⁻¹⁰, exchange against Lebedev 3·10⁻¹⁰;
  l = 2: 6·10⁻¹² (elastic), 10⁻⁸ (exchange, default screening).
- Cost: CO with cc-pV5Z on C (h, 5 q): elastic quartets 38 s, pairs 3.6 s.

**Not covered.** The Molden format stops at g (no agreed h ordering or normalisation), so
Molden/Molcas/Molpro input cannot carry h shells; PySCF can.

## Radial functions: upward and downward recurrences combined (5 October, afternoon)

**The p0 tables are still used**, for h < `cutoffcentre` = 10⁻¹⁰ bohr, i.e. in practice when the two
pair centres coincide exactly (all four functions on one atom, or the same pair twice). Between
10⁻¹⁰ and 10⁻⁴ bohr (the old threshold) the Bessel route is now used, because it is now exact
there.

**A mistake in the first version of the fix (corrected).** `radial_series` used the upward
recurrence wherever its error amplification looked small by a per-step estimate
Π(2m+1)/x². That estimate is too optimistic: the minimal solution j_b loses to the dominant
y_b at the rate (2b−1)(2b+1)/x² per order, so for orders above x the upward recurrence is
unstable. Measured against 40-digit references: order 10 off by 4·10⁻⁶ at x = 2.4, order 14 by
30 % at x = 3.4, order 18 by 5·10⁴ at x = 3.9. The tests did not catch it because the
independent (Lebedev) checks used s/p/d bases, where orders stay low, and the high-l comparison
was pairs against quartets, which shared the same routine.

**Now:** as in Miller's method, both directions are combined: upward for b ≤ x, backward from a
start chosen so the neglected minimal solution is below 10⁻¹⁷, matched at the junction; series +
downward for x < 1; pure upward only for x > b_max + 1. Accuracy, orders ≤ 18, 10⁻³ ≤ x ≤ 10³:

| x range | old `bessels0rr` | first fix | now |
|---|---|---|---|
| 10⁻³–0.1 | 7·10⁻³ (orders ≤ 3 closed forms) | 10⁻¹⁵ | 10⁻¹⁵ |
| 0.1–1 | 7·10⁻¹¹ | 10⁻¹⁵ | 10⁻¹⁵ |
| 1–10 | 10⁻¹⁵ | up to 5·10⁴ (orders > x) | 2·10⁻¹⁵ |
| 10–300 | 10⁻¹⁴ | 10⁻¹⁵ | 2·10⁻¹⁵ |
| 300–1000 | 10² (Miller start capped at 500) | 10⁻¹⁵ | 2·10⁻¹⁵ |

`tests/test_radial.py` checks every order up to 18 against SciPy/series values for
10⁻⁴ ≤ x ≤ 10³ and specifically the 1 < x < order region. With g functions (CO/cc-pVQZ)
against the Lebedev route: elastic l = 0 7·10⁻¹⁵ (pairs) and 2·10⁻¹⁴ (quartets), exchange
2·10⁻¹⁴; l = 2 3·10⁻¹² (pairs).

Cost of doing it right: benzene/6-31G HF total, 41 q: 4.9 s (earlier executable 5.3 s); elastic
timings in the table above were re-measured.

## j₂ at small q and near-coincident centres: fixed (5 October)

**Cause.** The quartet kernel (as the one-loop code before it) summed the radial functions with
closed forms for j₁, j₂, j₃ (`bessels0rr`), e.g. h₁ = (sin x/x² − cos x/x)(q/h). For small
x = qh these lose 2, 4, 6 powers of x to cancellation. Two symptoms: nearly coincident pair
centres lost digits (hence `cutoffcentre`), and the j₂ moment, which divides by q² after
Hobson's combination, was wrong at small q (36 % at q = 0.01 bohr⁻¹ for water).

**Fix.** `radial_series` (in `xs_radial.f90`, shared with the pair and electron–nucleus
kernels) evaluates Σ_b c_b q^{2b} j_b(qh)/(qh)^b stably for every qh: series + downward
recurrence for x < 1, upward recurrence where safe, Miller in between. The quartet kernel uses
it instead of `bessels0rr`.

**Result.**
- j₂ (exchange, CCSD two-electron, elastic) against the Lebedev projection: 10⁻¹¹ relative at
  q = 0.01 bohr⁻¹, 4·10⁻⁸ at 10⁻³ (before: 36 % and a factor 4700).
- j₀ near coincidence: the quartet kernel equals the exact pair route to 10⁻¹⁴ (before 2·10⁻¹⁰)
  and is exact even with no coincident branch at all, so `cutoffcentre` is now 10⁻¹⁰ (was 10⁻⁴).
  The pending "small-h expansion" is no longer needed.
- Cost: see the entry above (the first version of this fix had an unstable upward recurrence).
- All tests pass, and so do the slow large-basis tests (`PYXSCAT_SLOW=1`).

## Elastic scattering by density-weighted pairs (5 October)

**Question: are we using the optimal routine?** Yes for the general case: `xs_kernel_j0.f90` is
the port of the newest one-loop kernel (`BigMemoryScattering/pipeline_pyrazine/
Total_scattering_one_loop_opt.f90`, 25 Sept., contracted-shell quartets + screening tables),
byte-identical to the copy the regression references were made with. Benzene/6-31G HF total,
41 q, 2 threads: the original executable 5.2 s, ours 3.9 s; agreement 10⁻⁷ (its single-precision π).

**But elastic was going through the total loop**, at the cost of a total. Because D⊗D factorises,
`xs_elastic.f90` folds the density into each primitive-shell pair (one Hermite charge
distribution per pair) and sums over pairs of pairs with the closed-form Cartesian derivatives
of the radial functions (a_{Ti} h^{T−2i} f_{n}). Rigorous q-dependent screening
(|∂^{tuv} j₀| ≤ q^{t+u+v}), j_m(x)/x^m vectorised over q, exact at coincident centres.
Now the default (`elastic(..., method="pairs")`); the quartet route stays available.

| Elastic, 200 q, 2 threads | pairs | quartets |
|---|---|---|
| benzene / 6-31G | 9.8 s | 12.8 s |
| benzene / cc-pVDZ | 15.1 s | 25.6 s |
| benzene / aug-cc-pVDZ | 66 s | 132 s |
| CO / cc-pVQZ (7 q) | 0.57 s | 36 s |

Electron scattering (which needs the X-ray elastic term) and Molden-only runs benefit directly.

**Accuracy findings.**
- Pairs agree with the Lebedev route to 10⁻¹¹ (l = 0 and 2, q from 10⁻³), with the quartet kernel
  to 10⁻⁹, and give I(0) = N² to machine precision (quartets: 2·10⁻¹⁰).
- Same tolerance, tighter screening: 10⁻⁹ instead of 1.5·10⁻⁶ (benzene/aug-cc-pVDZ).
- Eirik's elastic j₂ (and our quartet kernel with his cutoffcentre = 0.01) is off by 2.5·10⁻⁴
  because near-coincident pair centres are treated as coincident; the pair route is exact. The
  test now checks both.
- The quartet kernel's j₂ loses accuracy at very small q (Hobson's cancellation, absolute error
  ~10⁻¹⁰/q²): wrong at q ≤ 0.03 bohr⁻¹. Elastic j₂ is now exact; the exchange and two-electron j₂
  still go through the quartet kernel. Fix planned (direct j₂P₂ derivatives).

`tests/test_elastic_pairs.py` (7 tests).

## Electron scattering (5 October)

**Physics.** First Born approximation, dσ/dΩ = 4/(a₀²q⁴) [I_nn − 2 I_ne + I_ee] (the bracket is
the quantity of the earlier code, modes 5 and 6). The only new integrals are the
electron–nucleus ones, Σ_A Z_A Σ D_μν ∫χ_μχ_ν j₀(q|r − R_A|).

**Kernel (`xs_nuclear.f90`).** Analytic: Hermite expansion at the pair centre, then the
McMurchie–Davidson recursion for the Cartesian derivatives of j₀(q|h|) with the auxiliary
functions f_n = (−q²)ⁿ j_n(qh)/(qh)ⁿ. j_n(x)/xⁿ by series (x < 1), Miller (x ≤ n_max + 1) or
upward recurrence, so nothing special happens at h = 0. j₂P₂ by Hobson. Density-weighted
Hermite coefficients per shell pair, OpenMP over pairs. Benzene/cc-pVDZ, 200 q: 0.7 s (the X-ray
elastic term takes 30 s), so electrons cost the same as X-rays.

**Python.** `pyxscat.electron`: `nucleus_nucleus`, `electron_nucleus`, `decompose`, `scattering`,
`cross_section`, `aligned` (|Z(q) − F(q)|² + inelastic at fixed orientation). Electron IAM in
`pyxscat.iam` (Mott–Bethe amplitudes Z − f; atomic and molecular parts for sM). `pyxscat --probe
electron` in the low-level command.

**Interface.** `probe="electron"` in `scattering`, `iam`, `anisotropy` and `detector`; the
cross section in Å²/sr (bohr²/sr with 1/bohr), NaN at q = 0; `energy_kev` for the relativistic
γ² and the detector wavelength (`th.electron_wavelength`); `th.sm(result, iam)` for sM(s), with
the IAM curve alongside. `thistle run ... --electrons | --probe both --energy 3000`;
`probe = "electron" | "both"` in input files. Electron files are named `<prefix>_electron_*`.

**Validation (`tests/test_electron.py`, 11 tests; 4 more in `test_thistle.py`).**
- I_ne against the Lebedev average of Re[Z*(q)F(q)] from the analytic pair transforms, l = 0
  and 2, Cartesian cc-pVDZ and spherical cc-pVQZ (g functions), q to 12 bohr⁻¹: 10⁻¹⁰.
- Against the old code (`Integrals.f90`, `nuclei_electron_integration`, compiled here): 4·10⁻⁸,
  the single-precision π of the old code. Reference stored in `tests/data`.
- Mott–Bethe for one atom: elastic = (Z − f)², f² the X-ray elastic: 10⁻⁸.
- H atom against the exact [1 − 16/(4+q²)²]²: 2·10⁻⁴ (aug-cc-pV5Z).
- Water at small q: elastic bracket / q² → μ²/3 with PySCF's dipole: 2·10⁻³.
- Fixed orientation averaged on a Lebedev grid = rotational, total/elastic/inelastic and l = 2.

**Observation.** For LiF (SA-CASSCF), sM(s) of the molecule departs strongly from the IAM between
1 and 4 Å⁻¹: electron diffraction weighs the valence density through (Z − f)², so charge
transfer shows far more than in X-ray scattering. PySCF's own SA-CASSCF gives the same curve as
the Molpro wavefunction, so this is physics, not the reader.

**Description.** `docs/THEORY.md` now gives the theory, the algorithms and the validation table
(the basis for the paper's methods section); the README opens with a feature overview.

## A bug in the old elastic j0 (found while testing j2)

`elastic_j0.f90` (all three old folders) calls the j1 Bessel sum `bessels1rr` instead of the j0
one. The elastic j0 result is right only at q = 1 bohr⁻¹: for water HF/6-31G it is 70 % too high
at q = 0.3 and off by 2e4 at q = 5. Details: `benchmarks/legacy/eirik_elastic_j0_bug.md`. Total
scattering and elastic j2 are not affected.

## Open points (6 October)

1. **Absolute Molcas energy.** PySCF's exact integrals give SCF and CASSCF energies about 24 mEh
   below the Molcas output; the input used `RICD`. A rerun without `ricd` settles it (the CI
   conversion is exact either way).
2. **Spherical d/f conventions of OpenMolcas and BAGEL Molden files**: untested (no such file
   yet); inputs in `benchmarks/molden_conventions/`.
3. **Eirik's aligned input** that gave a wrong answer in the earlier code: to compare.
4. **Cost with deeply contracted basis sets.** The quartet kernel works on primitive shells
   (cost ~ N_prim^4): ANO-L on C7H8 (225 primitive shells) ~30 s per elastic curve on 2 cores.
5. **ifx**: tested by one user on one cluster; the `-march=x86-64-v2` flag for the bitwise RDM
   builder is not yet tested with ifx.

(Done since the first list: electron scattering; coincident centres, through the stable radial
functions; partial alignment; h functions; RDMs from long expansions.)

## CAS/CI densities without the n⁴ array (7 October)

Question: does a large CI wavefunction (e.g. FCI) blow up the Python side of the 2RDM? Measured
stage by stage (`benchmarks/ci_rdm_scaling`, 2 cores, 7 GB):

| Stage | Case | Before | Now |
|---|---|---|---|
| CI vector → active RDMs | CAS(12,12), 8.5·10⁵ det. | dense 2.7 s / 1.2 GB; sparse 37 s / 3.1 GB; bitwise 216 s / 0.2 GB | unchanged (already compiled) |
| | CAS(14,14), 1.2·10⁷ det. | dense 39 s / 3.2 GB | unchanged |
| | FCI H₂O/6-31G, 1.7·10⁶ det. | dense 4.9 s / 1.4 GB; sparse 94 s / 3.3 GB; bitwise ≈ 13 min (est.) | unchanged |
| embedding + MO→AO | benzene CAS(6,6)/cc-pVDZ (120 AO) | 1.66 GB array, peak 5.9 GB, 8–10 s | 6 MB of factors, peak 0.14 GB, 0.03 s |
| | benzene CAS(6,6)/cc-pVTZ (300 AO) | 65 GB array (> 200 GB peak): impossible here | 39 MB, peak 0.33 GB, 0.1 s |
| kernel (30 q) | cc-pVDZ | 19.5 s | 28 s (24 or 81 terms: the same) |
| | cc-pVTZ | — | 650 s |

So the CI → RDM step was not the problem (the active RDMs are tiny, n_act⁴); the wall was the
four-index AO array that every CAS/CI route built afterwards (PySCF's `make_rdm12` even over all
MOs, virtuals included). Fixed by the factored embedding (THEORY §4), now the default beyond 60
basis functions for Molpro, Molcas and PySCF CASSCF/CASCI/FCI; `factored=False` restores the
array. Validated against the arrays to 10⁻¹⁴ (total), 10⁻¹⁰ (j₂), against PySCF's own AO array,
and through the existing Molpro/Molcas tests, which now run factored
(`tests/test_factored_cas.py`). Also fixed: picking one root of a state-averaged PySCF
calculation (the state-averaged solver needs `states_make_rdm12`).

Remaining cost driver for very large active spaces: the CI vector itself (CAS(16,16) = 1.7·10⁸
determinants, 1.3 GB per vector, plus the Python list of determinants when read from a file).

## First cluster run (7 October)

Node of partition `lambda`: 2 × Intel Xeon Gold 6326 (Ice Lake, 2 × 16 cores, 2.9 GHz),
gfortran 8.5 (Red Hat), Intel Python 3.12, NumPy 2.5.  Self-test passed; full suite with the slow
cases: 233 passed in 6 min on 32 threads.  Scaling, benzene/cc-pVDZ (120 functions), 100 q:

| threads | elastic (pair kernel) | HF total (quartet kernel) |
|---|---|---|
| 1 | 20.3 s | 41.3 s |
| 8 | 2.61 s (7.8×) | 5.55 s (7.4×) |
| 16 | 1.34 s (15.1×) | 3.55 s (11.6×) |
| 32 | 0.70 s (28.8×) | 2.88 s (14.3×) |

The quartet kernel stalled at ~14×: its OpenMP loop handed out one task per first shell
index, 48 tasks of very unequal size (the largest 6.8 % of the work, so at most 14.7× however many
threads).  A model of that schedule reproduces the measured speed-ups (7.7×, 13.0×, 14.8×).  Fixed:
one task per shell pair (1176 here), largest first (merge sort on a cost estimate); the model now
gives 32.0× on 32 threads.  Results unchanged to rounding (1e-14 relative; full suite passes).
To be re-measured on the cluster.

## j₂ sign convention: i² included (8 October)

At a collaborator's request S₂ now carries the phase i² = −1 of the l = 2 term of the plane-wave
expansion: S₂ = −Σ Γ⟨j₂P₂⟩, and I(q, θ_q) = I₀ + 5⟨P₂⟩P₂(cos θ_q)S₂ (all moments with plus signs;
S₂ = ⟨I P₂⟩ = a₂/5 for perfect alignment).  One place (`core.lsign`) for the kernels; the IAM,
electron (nucleus–nucleus, electron–nucleus) and aligned-ensemble formulas follow.  Every j₂ value
written before this date is −1 times the new one.  Tests and references updated (Eirik's Y₂₀ code:
factor −√(5/4π)).

## MRCI (8 October)

Internally contracted MRCI (Molpro) has no CI vector to export: the external part is expanded in
contracted pair functions and the 2RDM is never formed; only the 1RDM (`dm` / natural orbitals)
is available.  OpenMolcas's MRCI is uncontracted (GUGA direct CI) and prints the complete vector
with `THRPrint = 0.0`: a reader (`interfaces.molcas.read_mrci_output`, `from_mrci`,
`th.load` on the output) converts the CSFs (GUGA case numbers 0/1/2/3 → 0/u/d/2, coupling order
virtual, inactive, active; the determinants reordered to increasing orbital index with the
permutation sign) and builds the RDMs over the correlated orbitals with the determinant builders,
adding the frozen orbitals as core.  Water/cc-pVDZ, CASSCF(6,5) reference, 50 034 CSFs: the MRCI
energy from our RDMs reproduces Molcas's to 8e-8 (the printed six decimals), the natural
occupations to 1e-5; scattering factored vs array 2e-11 (`tests/test_molcas_mrci.py`,
`examples/mrci/molcas_test`).

Same day: the 64-orbital limit is gone (multiword bit strings in the Fortran pair kernel, shared
accumulation with atomic updates when the RDM copies would exceed 2 GB), and
`ci.make_rdm12_external` uses the ≤ 2-electrons-in-the-virtuals structure of MRCI/CISD-like
expansions: grouped pair loop over coupling labels + the vvvv block by matrix products (THEORY
§4).  Verified: water MRCI structured vs full loop 3e-14; the same vector with its virtuals
spread over 120 orbitals (two words, atomic path) 9e-14; a CAS(10,8) vector on scattered
orbitals of a 128-orbital space vs PySCF 2e-14.  The MRCI reader uses it by default, and the
four-index array route (kernel mode 1) is kept beyond ~80 orbitals, where the factored embedding
(eigendecomposition of the active block) is no longer feasible: pyridine MRCI/cc-pVTZ (260
correlated orbitals) means a 36 GB MO array and a 65 GB AO array, for lambda-hmem.  Remaining
cost for very large MRCI vectors: the CSF → determinant expansion and the grouping are in
Python (~10 µs per determinant).  `examples/mrci`
also has the Molpro input for what Molpro can give, and the selected-CI route (PySCF, explicit
determinants, bitwise builder) for large active spaces.

## Before the paper

Must have:
1. ~~Licence~~ (Apache 2.0, done 7 October), public git repository, continuous integration (gfortran Linux/macOS, ifx), a
   release with a DOI (Zenodo).
2. Performance on the cluster: thread scaling, basis-set scaling, comparison with the earlier
   code on the same machine (`benchmarks/cluster_check.py scaling` and a benchmark set).
3. A `paper/` folder: one script per figure and table, from inputs to the final plot.
4. The application cases and the reproduction of published results of the earlier code.
5. Open points 1-3 above.

To do (7 October): j₂ about an arbitrary molecular axis. S₂ is now the moment about the
molecular z axis of the input frame only (the laboratory axis is free); add a `molecular_axis`
option or compute the full l = 2 tensor (five components, then any axis). See THEORY §6.

Could have: MCCI reader, j1 terms, contracted screening in the quartet kernel (open point 4),
BAGEL/NEXMD adapters. (Factored CCSD densities: done.)

## Findings in the existing code

1. **Single-precision π.** `Total_scattering_one_loop_opt.f90` computes π as `acos(-1.0000)`,
   which scales every result by 1 + 8.3e-8.
2. **Old default switch-over distance.** The old default (0.01 bohr) gives errors of about
   4e-6 for ordinary molecules; the new default (1e-4 bohr) gives below 1e-9.
3. **Spin in the old Molcas reader.** The old reader takes `int(float(S))` as M_S, which gives
   the wrong determinant count for doublets and other half-integer spins. The new
   `csf_to_determinants` takes 2M_S explicitly.
