# THISTLE (this code) and thistle-dev (KirranderGroup/thistle-dev): comparison

Checked 7 October 2026 against thistle-dev at commit 1540087 (9 July 2026, Patrick Wang),
in view of merging the two.

## Same numbers

thistle-dev's five test decks (CH₄ STO-3G total, CH₄ 3-21G elastic, Be 3-21G elastic, and the
two dense-AO variants) read into pyxscat (now `pyxscat.interfaces.deck`, used by `th.load("deck.h5")`; `load_deck.py` was the first version) (no conversion: the deck's MO
coefficients over primitives already carry contraction and normalisation, and its 2RDM has
PySCF's convention, trace N(N-1); our I(0) = N² on every deck):

| Deck | Quantity | Our default (exact) vs golden | Same cutoffcentre (0.01) as the golden |
|---|---|---|---|
| be_elastic | elastic | 3.6e-10 | — |
| ch4_total (+ _zcotr) | two-electron term | 1.6e-5 | **4.1e-13** |
| ch4_elastic (+ _zcotr) | elastic | 7.9e-6 | **4.1e-13** |

after removing the uniform factor 1 + 8.35e-8 of the single-precision π in the legacy Fortran
(their goldens carry it; their Julia engine has the same `fortran_pi` switch). Notes:

* their "total" golden is the two-electron term only (it starts at N(N-1) = 90, not N² = 100);
* the 1e-5 differences are the coincident-centre approximation, `cutoffcentre = 0.01`, which their
  decks and their Python default still use; with the same value the two codes agree to 4e-13;
* beyond q ≈ 50 bohr⁻¹ their goldens hold ~1e-12 rounding noise that our screening sets to zero.

Run: `python compare_goldens.py /path/to/thistle-dev` (with pyxscat importable).

## Their open j₂ question

`docs/FINDINGS_colleague_code_j1j2.md` in thistle-dev reports that Eirik's j₂ (~1e-4 for CH₄)
and their Julia engine's j₂ (~0.1–0.7) disagree, and leaves it open. Our j₂ on their CH₄ deck
(`j2.py`):

| q (bohr⁻¹) | S₂, j₂ kernel | Lebedev projection of the fixed-orientation term | rel. diff |
|---|---|---|---|
| 0.5 | 1.1165e-4 | 1.1165e-4 | 5e-9 |
| 1.0 | 1.9142e-4 | 1.9142e-4 | 8e-9 |
| 2.0 | 6.2251e-5 | 6.2251e-5 | 2e-8 |
| 4.0 | −5.8972e-6 | −5.8972e-6 | 6e-9 |

So Eirik's magnitude is right and the Julia engine's j₂ is wrong (Eirik's Y₂₀ convention is ours ×
√(5/4π) = 0.6308).

## Scope

| | THISTLE (this code) | thistle-dev v1 (beta) |
|---|---|---|
| X-ray, rotational average | total, elastic, inelastic; Hartree/exchange/cumulant | total (two-electron term) and elastic |
| j₂ / partial alignment | yes, validated | gated off (j1/j2 "pending validation") |
| fixed orientation, detector grids | yes (incl. grid.f90's polar grid) | not ported |
| electron scattering (UED/GED), sM(s) | yes | no |
| IAM | ab initio and tables, X-ray and electron | Python modules `iam` (X-ray and electron) and `gas_iam` (xraydb tables) |
| coincident centres | exact limits (default) | `cutoffcentre = 0.01` approximation |
| CCSD | factored density (no n⁴ array) | via PySCF `make_rdm2` (dense) |
| CI vectors → RDMs | Molpro, Molcas (CSFs), dense/bitwise/sparse builders | PySCF objects (CASCI/CASSCF, FCI, DMRG-SCF, SHCI-SCF via `make_rdm12`) |
| input | in memory (PySCF, Molden, Molpro output/punch, Molcas) + TOML input files | the "deck" (HDF5 schema v1, or legacy .dat), written by Python |
| engines | one Fortran library via ctypes (gfortran or ifx, no MKL, built by pip) | Fortran (ifort + MKL + HDF5) and Julia (XScattering.jl, no compile) |
| angular momentum | up to h | Cartesian up to g (d/f/g need `cart=True`) |
| licence | Apache 2.0 | MIT (pyproject) |

## What thistle-dev has that is worth bringing over

1. The **HDF5 deck** as a portable, self-describing input/output format (good for moving a
   calculation to a cluster). Reader: done (`th.load("deck.h5")`, both representations, all five
   test decks to 4e-13 with the deck's cutoffcentre). Writer: to do.
2. The **Julia engine** as an independent implementation for cross-checks in CI (after its j₂ is
   fixed against ours).
3. The generic **PySCF `make_rdm12` route** (DMRG-SCF, SHCI-SCF objects).
4. Their **IAM modules** and logging/banner conventions, to compare with ours.
