# Python interfaces: thistle-dev `scat_lib` (Patrick Wang) and THISTLE (`thistle`)

9 October 2026, thistle-dev `main` at 1540087. Nothing in thistle-dev was changed; its Fortran
engine was built here with gfortran from a patched scratch copy (below).

## The two designs

| | `scat_lib` (thistle-dev) | `thistle` (this code) |
|---|---|---|
| call | `run_scattering_pyscf(obj, mf, name, backend, engine, type, state, …)` | `th.load(obj).scattering(...)`, `.anisotropy`, `.detector`, … |
| data path | PySCF → Molden file → RDMs → HDF5 deck on disk → `Main.exe` (subprocess) or Julia → text file read back | PySCF → RDMs in memory → Fortran library (ctypes), same process |
| observables per call | one (`total` two-electron term, or `elastic`) | all parts at once (total, elastic, inelastic, Hartree, exchange, cumulant), j₂, fixed orientation, electrons, IAM |
| density representations | MO (`normal`: dense nmo⁴ array) or AO (`zcotr`: dense nao⁴ array) | dense or factored (no n⁴ array), chosen automatically |
| sources | PySCF: CASCI/CASSCF (any fcisolver, DMRG/SHCI), FCI solver, CCSD, any object with `make_rdm12` | PySCF (SCF, CCSD, CISD, CASCI/CASSCF with any fcisolver, FCI solver), Molden, Molpro, OpenMolcas (RASSCF, MRCI), thistle-dev decks |
| basis | Cartesian only (`cart=True` required for d and higher) | Cartesian and spherical |
| build | ifort + MKL (+ HDF5 C) | gfortran or ifx, no MKL, built by pip |

The deck is the useful part of `scat_lib`'s design as a portable input format, and THISTLE already
reads it (`th.load("deck.h5")`, both representations) and writes it (`Molecule.write_deck`).

## Same numbers (`compare_interfaces.py`)

Water, 81 q points to 8 bohr⁻¹. A: `scat_lib` → its engine. B: the same deck → THISTLE kernels.
C: `th.load(obj)` → THISTLE kernels. Relative max differences:

| wavefunction | A–B (cutoffcentre 0.01) | A–B (1e-10) | B–C |
|---|---|---|---|
| CASSCF(4,4)/6-31G | 2–5e-15 | 3–6e-15 | 2e-14 |
| CCSD/6-31G | 3–7e-15 | 4–6e-15 | 2e-14 |
| DMRG-CASCI(6,6)/6-31G (block2, M = 500; ΔE vs CASCI 3e-14) | 3–6e-15 | 5–6e-15 | 2e-14 |
| CASSCF(6,6)/cc-pVDZ (Cartesian) | 3–8e-15 | 1.6e-13 | 2e-14 |
| FCI/STO-3G, roots 1 and 3 (`fci.FCI(mf)`) | 1–4e-15 | 2–3e-15 | 2e-14 |
| CASCI nroots=2, root 2 | 3e-15 | 3e-15 | (vs SA-CASSCF root 2: 9e-9, PySCF's two objects, same through THISTLE alone) |

Every case, both deck representations (MO and AO), both observables. The kernels agree on
identical input, and the two interfaces build the same densities from the same PySCF objects;
`scat_lib`'s `orbital_type` (HF, CASSCF, natural) changes its result by ≤1.3e-15, as it should.

## Where they differ

* **State-averaged CASSCF fails in `scat_lib`, for every root**: it passes one CI vector to the
  state-averaged solver's `make_rdm12`, which expects the list of roots (AssertionError in
  `pyscf.fci.rdm`). THISTLE: `th.load(sa, state=2)`. Workaround there: CASCI with `nroots` on the
  SA orbitals.
* **MO deck route does not scale**: water/cc-pVDZ (25 orbitals, 41 primitives), total, exceeds 7 GB
  in the engine and is killed; the AO route takes 0.3 s.
* **AO route holds n⁴ arrays in Python** (benzene CASSCF(6,6), 81 q, 4 threads, `scaling.py`):

  | basis (n) | | scat_lib (AO) | THISTLE |
  |---|---|---|---|
  | 6-31G (66) | elastic | 29 s, 0.7 GB | 26 s, 0.17 GB (all parts) |
  | | total | 31 s, 0.9 GB | |
  | cc-pVDZ (120) | elastic | 187 s, 5.3 GB | 64 s, 0.5 GB (all parts) |
  | | total | killed (out of memory > 6 GB) | |

* **Standalone FCI solvers** were not read by THISTLE; now `th.load(fci.FCI(mf).run(), mo_coeff=mf)`
  (PySCF keeps only the integrals in the solver, not the orbitals, so they are passed, as
  `scat_lib` does with `mf`).
* `scat_lib` extras not in THISTLE: `fit_utils` (CI coefficients fitted to a reference curve),
  `sine_transform` (pair distribution function), the Julia engine.

## Recommendation

Keep `thistle` as the Python interface: in-process, all observables per call, factored
densities, more sources, no Intel toolchain, and identical numbers. Keep the HDF5 deck as the
exchange format (both directions already work), and port from `scat_lib` the generic
`make_rdm12` fallback, the fitting and PDF utilities, and the Julia engine as an independent
cross-check in CI.

## Building thistle-dev's engine with gfortran

`thistle_dev_gfortran.patch` (casts and one `use ifport`; type mismatches that ifort lets
through: `dphrec`, `combo`, `integer2binary_orbs`, mixed-kind `iand`, `==` on logicals) and
`build_thistle_dev_gfortran.sh` (OpenBLAS, h5py's HDF5 C libraries). With them its own golden
tests pass: 5/5 on both the HDF5 and the `.dat` path, ≤4e-14.

## The kernels themselves

thistle-dev's Fortran is the earlier group code (Eirik's `Integrals.f90` MO engine) with Patrick's
dense-AO kernel (`Total_scattering_j0_groups_fast.f90`) grafted on; its Julia engine is a port of the
same routines. THISTLE's kernels descend from the same code. Differences found:

| | thistle-dev (Fortran and Julia) | THISTLE |
|---|---|---|
| j₀ derivatives | per primitive quartet, polynomial tables (`rrdj0`, `BesselDeriv`) | McMurchie–Davidson recursion seeded by g_n = j_n/xⁿ |
| radial functions, orders 1–3 | closed forms: relative error 5e-8 at qh = 0.1, 8e-2 at 0.01, 2e5 at 1e-3 | series / downward / upward + Miller: ≤ 1e-15 everywhere (`kernels_radial_driver.f90`, against mpmath) |
| radial functions, orders ≥ 4 | Miller started at fixed order 100: 1e-14 up to qh = 50, 1e-3 at 90, O(1) beyond 100 | ≤ 1e-15 up to qh = 10³ |
| coincident centres | `cutoffcentre` = 0.01 bohr → P0 tables (needed because of the closed forms) | exact limit inside the recursion |
| angular momentum | total L ≤ 19 per quartet (s–g) | up to h (l = 5) tested |
| screening | density and coefficient thresholds only (`cutoffz`, `cutoffmd` = 1e-20); every primitive quartet at every q | Gaussian-prefactor and density bounds per shell pair and q |
| density gather | from the n⁴ array per primitive quartet | per contracted-shell quartet; or factored (no n⁴ array) |
| elastic | quartet loop on D⊗D | density-weighted pairs (generalised Debye sum) |
| threads | `OMP_set_num_threads(16)` hard-coded (ignores OMP_NUM_THREADS) | as requested |
| π | `acos(-1.0000)`, single precision: all results × (1 + 8.35e-8) | double |
| MO-route memory | four tensors of N_nonzero(Γ_MO) × n_prim² doubles (water/cc-pVDZ: ~20 GB) | — |
| other observables in the code | j₁, j₂, aligned/detector, electron–nucleus, coherent cross terms (Typec 3–10; gated off in v1) | all of these except j₁, validated |

Consequences, same decks, same `cutoffcentre` so only the radial functions differ:

* CO/cc-pVDZ CASSCF(6,6): agreement 1e-13 to q = 10, 1e-11 at 30, then 2e-8 at q = 100 and 1e-6 at
  250 bohr⁻¹ (thistle-dev's default q range is 0–250). THISTLE is the accurate one (radial test above).
* Butadiene/6-31G* (11 bohr long), elastic: 2e-13 up to q = 8 bohr⁻¹ (15 Å⁻¹), 5e-10 at 10, 2e-6 at
  15. At experimental q the 0.01 cut-off (1–5e-6 at every q) is the larger error; the radial one
  grows with q × molecular size.
* Speed (benzene/6-31G CASSCF(6,6), 81 q, 4 threads each, kernel only on the same deck):
  total 29.5 s vs 8.0 s (dense) / 13.9 s (factored); elastic 29.0 s vs 4.3 s (pairs).
