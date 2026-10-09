# THISTLE: independent checks, line by line

Run of `validation/independent/run_independent_checks.py`, 9 October 2026: **80 of 80 checks pass.**

Every line compares a THISTLE result with a reference that does not go through THISTLE's kernels. The references are:

- **PySCF analytic Fourier transforms** of basis-function pairs (`ft_aopair`) with PySCF's own density matrices (`make_rdm1`, `make_rdm2`, `mcscf.addons.make_rdm12`), averaged over directions on Lebedev grids (2030 points, degree 77; 5810 points above q = 8 bohr⁻¹). Shares no code with THISTLE.
- **Numerical integration** of the electron density on PySCF's DFT grids (no Gaussian integrals at all).
- **mpmath** at 40 digits; **xraydb** atomic form factors (Waasmaier–Kirfel); **CODATA** constants (scipy).
- **Closed forms and exact limits**: N², N(N−1), N ΣZ, the dipole limit q²μ²/3, the hydrogen atom.
- **Numbers printed by other programs** (Molpro energies, OpenMolcas energies and occupations), the stored output of the earlier production kernels, and thistle-dev's Fortran engine.

Columns: *abs* = max |THISTLE − reference| over the points listed; *rel* = abs / max |reference| (for lines marked *abs*, the tolerance applies to the absolute value). q in bohr⁻¹; Lebedev checks at q = 0.5, 1.5, 3.0 and 5.0 bohr⁻¹ unless stated.

## Summary

| Area | Checks | Pass | Largest deviation (rel or abs) |
|---|---|---|---|
| A. X-ray, rotationally averaged | 24 | 24 | 1.7e-09 |
| B. Exact limits and sum rules | 5 | 5 | 2.6e-10 |
| C. j2 moment (aligned molecules) | 4 | 4 | 2.2e-11 |
| D. Fixed orientation, detector grids | 6 | 6 | 7.6e-10 |
| E. Electron scattering | 7 | 7 | 1.0e-04 |
| F. Exact one-electron results | 3 | 3 | 2.6e-05 |
| G. Independent atom model | 7 | 7 | 5.3e-03 |
| H. Radial functions j_n(x)/x^n | 1 | 1 | 8.2e-15 |
| I. Density matrices | 8 | 8 | 3.4e-10 |
| J. Readers (numbers printed by Molpro and OpenMolcas) | 4 | 4 | 3.4e-06 |
| K. Other implementations | 11 | 11 | 8.3e-08 |

## A. X-ray, rotationally averaged

| ID | Quantity | System | Independent reference | abs | rel / abs | Tolerance | Result |
|---|---|---|---|---|---|---|---|
| A1 | elastic <\|F\|^2>, RHF | water/cc-pVDZ | PySCF ft_aopair + Lebedev 2030 | 1.1e-11 | 1.3e-13 (rel) | 1.0e-09 | PASS |
| A2 | inelastic, RHF (Waller-Hartree) | water/cc-pVDZ | N - 2 sum_ij <\|F_ij\|^2>, PySCF ft | 1.8e-10 | 2.1e-11 (rel) | 1.0e-09 | PASS |
| A3 | total, RHF | water/cc-pVDZ | elastic + Waller-Hartree inelastic, PySCF | 1.8e-10 | 2.1e-12 (rel) | 1.0e-09 | PASS |
| A4 | elastic, RHF (f functions) | CO/cc-pVTZ | PySCF ft_aopair + Lebedev | 5.4e-11 | 3.4e-13 (rel) | 1.0e-09 | PASS |
| A5 | inelastic, RHF (f functions) | CO/cc-pVTZ | Waller-Hartree, PySCF ft | 8.0e-09 | 7.0e-10 (rel) | 1.0e-09 | PASS |
| A6 | elastic, RHF, spherical basis | water/cc-pVDZ (5d) | PySCF ft_aopair on spherical AOs | 1.4e-11 | 1.6e-13 (rel) | 1.0e-09 | PASS |
| A7a | total, CCSD | water/cc-pVDZ | N + <sum Gamma F F*>, PySCF Gamma + ft | 2.9e-10 | 3.4e-12 (rel) | 1.0e-09 | PASS |
| A7b | total, CCSD, factored density | water/cc-pVDZ | same reference | 1.8e-10 | 2.1e-12 (rel) | 1.0e-09 | PASS |
| A7c | elastic, CCSD | water/cc-pVDZ | <\|sum D F\|^2>, PySCF D + ft | 1.8e-11 | 2.1e-13 (rel) | 1.0e-09 | PASS |
| A7d | inelastic, CCSD | water/cc-pVDZ | difference of the two references | 2.7e-10 | 3.2e-11 (rel) | 1.0e-09 | PASS |
| A7e | exchange part, CCSD | water/cc-pVDZ | -1/2 <sum D D F F*>, PySCF | 1.8e-10 | 2.0e-11 (rel) | 1.0e-09 | PASS |
| A7f | cumulant part, CCSD | water/cc-pVDZ | total - N - Hartree - exchange (refs) | 4.5e-10 | 7.2e-10 (rel) | 1.0e-07 | PASS |
| A8a | total, CISD | water/cc-pVDZ | N + <sum Gamma F F*>, PySCF Gamma + ft | 3.5e-10 | 4.1e-12 (rel) | 1.0e-09 | PASS |
| A8b | total, CISD, factored density | water/cc-pVDZ | same reference | 6.9e-11 | 8.0e-13 (rel) | 1.0e-09 | PASS |
| A8c | elastic, CISD | water/cc-pVDZ | <\|sum D F\|^2>, PySCF D + ft | 2.0e-11 | 2.4e-13 (rel) | 1.0e-09 | PASS |
| A8d | inelastic, CISD | water/cc-pVDZ | difference of the two references | 3.3e-10 | 3.9e-11 (rel) | 1.0e-09 | PASS |
| A8e | exchange part, CISD | water/cc-pVDZ | -1/2 <sum D D F F*>, PySCF | 9.2e-11 | 1.0e-11 (rel) | 1.0e-09 | PASS |
| A8f | cumulant part, CISD | water/cc-pVDZ | total - N - Hartree - exchange (refs) | 4.1e-10 | 6.9e-10 (rel) | 1.0e-07 | PASS |
| A9a | total, CASSCF(6,6) | water/cc-pVDZ | N + <sum Gamma F F*>, PySCF Gamma + ft | 2.7e-10 | 3.1e-12 (rel) | 1.0e-09 | PASS |
| A9b | total, CASSCF(6,6), factored density | water/cc-pVDZ | same reference | 4.9e-11 | 5.6e-13 (rel) | 1.0e-09 | PASS |
| A9c | elastic, CASSCF(6,6) | water/cc-pVDZ | <\|sum D F\|^2>, PySCF D + ft | 2.1e-11 | 2.4e-13 (rel) | 1.0e-09 | PASS |
| A9d | inelastic, CASSCF(6,6) | water/cc-pVDZ | difference of the two references | 2.4e-10 | 2.9e-11 (rel) | 1.0e-09 | PASS |
| A9e | exchange part, CASSCF(6,6) | water/cc-pVDZ | -1/2 <sum D D F F*>, PySCF | 2.3e-10 | 2.6e-11 (rel) | 1.0e-09 | PASS |
| A9f | cumulant part, CASSCF(6,6) | water/cc-pVDZ | total - N - Hartree - exchange (refs) | 4.8e-10 | 1.7e-09 (rel) | 1.0e-07 | PASS |

- **A1**: reference converged: Lebedev 2030 vs 3074 differ by 6.6e-16
- **A7f**: small quantity: relative to its own size
- **A8f**: small quantity: relative to its own size
- **A9f**: small quantity: relative to its own size

## B. Exact limits and sum rules

| ID | Quantity | System | Independent reference | abs | rel / abs | Tolerance | Result |
|---|---|---|---|---|---|---|---|
| B1 | I_tot(0) = N^2 | water CCSD | N^2 = 100 | 3.0e-10 | 3.0e-12 (rel) | 1.0e-10 | PASS |
| B2 | I_el(0) = N^2 | water CCSD | N^2 = 100 | 4.7e-11 | 4.7e-13 (rel) | 1.0e-10 | PASS |
| B3 | I_inel(0) = 0 | water CCSD | 0 | 2.6e-10 | 2.6e-10 (abs) | 1.0e-09 | PASS |
| B4 | total at large q (8, 12, 16 bohr^-1) | water CCSD | N + <sum Gamma F F*>, PySCF, Lebedev 5810 | 4.3e-10 | 4.0e-11 (rel) | 1.0e-09 | PASS |
| B5 | sum Gamma S S = N(N-1) | water CCSD (PySCF Gamma) | 90 | 2.6e-13 | 2.8e-15 (rel) | 1.0e-10 | PASS |

- **B4**: two-electron term there: 6.3e-01, 2.9e-01, 1.1e-01 (tight core pairs keep it finite in a Gaussian basis; it does not vanish)

## C. j2 moment (aligned molecules)

| ID | Quantity | System | Independent reference | abs | rel / abs | Tolerance | Result |
|---|---|---|---|---|---|---|---|
| C1a | S2 elastic, RHF | water/cc-pVDZ | <I_el(q n) P2(n_z)>, PySCF ft + Lebedev | 7.4e-13 | 6.0e-12 (rel) | 1.0e-08 | PASS |
| C1b | S2 total, RHF | water/cc-pVDZ | <I_tot(q n) P2(n_z)>, PySCF | 2.7e-12 | 2.2e-11 (rel) | 1.0e-08 | PASS |
| C2a | S2 elastic, CCSD | water/cc-pVDZ | <I_el(q n) P2(n_z)>, PySCF ft + Lebedev | 7.1e-13 | 5.7e-12 (rel) | 1.0e-08 | PASS |
| C2b | S2 total, CCSD | water/cc-pVDZ | <I_tot(q n) P2(n_z)>, PySCF | 2.8e-12 | 2.2e-11 (rel) | 1.0e-08 | PASS |

## D. Fixed orientation, detector grids

| ID | Quantity | System | Independent reference | abs | rel / abs | Tolerance | Result |
|---|---|---|---|---|---|---|---|
| D1 | elastic \|F(q)\|^2, 40 random q vectors | water CCSD | PySCF ft_aopair | 2.8e-14 | 3.6e-16 (rel) | 1.0e-10 | PASS |
| D2 | total, 40 random q vectors | water CCSD | N + sum Gamma F F*, PySCF | 2.8e-14 | 3.5e-16 (rel) | 1.0e-10 | PASS |
| D3 | photon k0 at 10 keV | - | E/(hbar c), CODATA (scipy) | 2.0e-09 | 7.6e-10 (rel) | 1.0e-06 | PASS |
| D4 | \|q\| on the Ewald sphere | 9 directions | 2 k0 sin(2theta/2) | 4.4e-16 | 1.7e-16 (rel) | 1.0e-12 | PASS |
| D5 | polar detector grid (qcoord3) | 10 keV, 60 deg, 49 x 73 | qgrid.dat written by the original grid.f90 | 5.0e-13 | 5.0e-13 (abs) | 1.0e-11 | PASS |
| D6 | aligned ensemble, elastic, <P2> = 0.4 | water RHF, 4 lab q vectors | explicit average over 590 x 24 orientations, PySCF ft | 2.1e-11 | 4.8e-13 (rel) | 1.0e-07 | PASS |

- **D5**: 3529 points; file has 12 decimals

## E. Electron scattering

| ID | Quantity | System | Independent reference | abs | rel / abs | Tolerance | Result |
|---|---|---|---|---|---|---|---|
| E1 | electron-nucleus term I_ne(q) | water RHF | sum_A Z_A int rho j0(q\|r-R_A\|), DFT grid (level 6) | 6.0e-08 | 6.0e-10 (rel) | 1.0e-07 | PASS |
| E2 | I_ne(0) = N sum Z | water RHF | 10 x 10 | 5.7e-14 | 5.7e-16 (rel) | 1.0e-10 | PASS |
| E3 | nucleus-nucleus term | water | sum Z_A Z_B j0(q R_AB) | 2.8e-14 | 2.8e-16 (rel) | 1.0e-14 | PASS |
| E4 | electron bracket, total, CCSD | water/cc-pVDZ | I_nn - 2 I_ne + I_ee, PySCF ft + Lebedev | 2.9e-10 | 5.6e-12 (rel) | 1.0e-09 | PASS |
| E5 | elastic bracket -> q^2 mu^2/3 (q = 2e-3) | water RHF | PySCF SCF dipole | 8.8e-11 | 1.0e-04 (rel) | 3.0e-03 | PASS |
| E6 | electron k0 at 100 keV (relativistic) | - | p/hbar, CODATA (scipy) | 2.4e-06 | 2.7e-08 (rel) | 1.0e-06 | PASS |
| E7 | cross section from the bracket | - | 4 bracket / (a0^2 q^4) (bohr^2, no gamma) | 0.0e+00 | 0.0e+00 (rel) | 1.0e-14 | PASS |

- **E1**: grid integrates N to 6.0e-09
- **E5**: leading term only; next order is O(q^4)
- **E7**: units of the interface: bohr^2/sr for q in bohr^-1

## F. Exact one-electron results

| ID | Quantity | System | Independent reference | abs | rel / abs | Tolerance | Result |
|---|---|---|---|---|---|---|---|
| F1 | H atom elastic F(q)^2 | H/aug-cc-pV5Z UHF | exact [16/(4+q^2)^2]^2 | 2.0e-05 | 2.6e-05 (rel) | 3.0e-04 | PASS |
| F2 | H atom inelastic 1 - F^2 (one electron: I_tot = 1) | H/aug-cc-pV5Z UHF | exact 1 - F(q)^2 | 2.0e-05 | 2.0e-05 (rel) | 3.0e-04 | PASS |
| F3 | H atom electron elastic bracket | H/aug-cc-pV5Z UHF | exact (1 - F)^2 (Mott-Bethe) | 2.1e-05 | 2.3e-05 (rel) | 3.0e-04 | PASS |

- **F1**: UHF energy error 5.2e-06 Eh: basis-set limit

## G. Independent atom model

| ID | Quantity | System | Independent reference | abs | rel / abs | Tolerance | Result |
|---|---|---|---|---|---|---|---|
| G1H | f(q) of the free atom H | HF/aug-cc-pVQZ | xraydb f0 (Waasmaier-Kirfel) | 1.6e-04 | 1.8e-04 (rel) | 1.0e-02 | PASS |
| G1C | f(q) of the free atom C | HF/aug-cc-pVQZ | xraydb f0 (Waasmaier-Kirfel) | 2.9e-02 | 5.2e-03 (rel) | 1.0e-02 | PASS |
| G1N | f(q) of the free atom N | HF/aug-cc-pVQZ | xraydb f0 (Waasmaier-Kirfel) | 3.5e-02 | 5.3e-03 (rel) | 1.0e-02 | PASS |
| G1O | f(q) of the free atom O | HF/aug-cc-pVQZ | xraydb f0 (Waasmaier-Kirfel) | 5.8e-03 | 7.7e-04 (rel) | 1.0e-02 | PASS |
| G1F | f(q) of the free atom F | HF/aug-cc-pVQZ | xraydb f0 (Waasmaier-Kirfel) | 4.6e-03 | 5.4e-04 (rel) | 1.0e-02 | PASS |
| G1Ne | f(q) of the free atom Ne | HF/aug-cc-pVQZ | xraydb f0 (Waasmaier-Kirfel) | 4.1e-03 | 4.3e-04 (rel) | 1.0e-02 | PASS |
| G2 | IAM elastic (Debye sum) | water, xraydb f | sum f_A f_B j0(q R_AB) written here | 3.6e-15 | 4.3e-17 (rel) | 1.0e-12 | PASS |

- **G1H**: method difference (HF in a finite basis vs relativistic HF tables): agreement to ~1% of Z expected
- **G1C**: method difference (HF in a finite basis vs relativistic HF tables): agreement to ~1% of Z expected
- **G1N**: method difference (HF in a finite basis vs relativistic HF tables): agreement to ~1% of Z expected
- **G1O**: method difference (HF in a finite basis vs relativistic HF tables): agreement to ~1% of Z expected
- **G1F**: method difference (HF in a finite basis vs relativistic HF tables): agreement to ~1% of Z expected
- **G1Ne**: method difference (HF in a finite basis vs relativistic HF tables): agreement to ~1% of Z expected

## H. Radial functions j_n(x)/x^n

| ID | Quantity | System | Independent reference | abs | rel / abs | Tolerance | Result |
|---|---|---|---|---|---|---|---|
| H1 | q^2b j_b(x)/x^b, b <= 18, x = 1e-3 .. 1e3 | - | mpmath, 40 digits | 8.2e-15 | 8.2e-15 (abs) | 1.0e-13 | PASS |

- **H1**: pointwise relative error, worst over all orders and x

## I. Density matrices

| ID | Quantity | System | Independent reference | abs | rel / abs | Tolerance | Result |
|---|---|---|---|---|---|---|---|
| I1 | factored 2RDM, CCSD (MO basis) | water/cc-pVDZ | PySCF make_rdm2 | 3.5e-14 | 3.5e-14 (abs) | 1.0e-12 | PASS |
| I2 | factored 2RDM, CISD (MO basis) | water/cc-pVDZ | PySCF make_rdm2 | 3.6e-14 | 3.6e-14 (abs) | 1.0e-12 | PASS |
| I3 | factored 2RDM, CASSCF(6,6) (AO basis) | water/cc-pVDZ | PySCF mcscf.addons.make_rdm12 | 6.2e-15 | 6.2e-15 (abs) | 1.0e-12 | PASS |
| I4d | RDM builder 'dense', CAS(8,8) FCI vector | N2/6-31G | PySCF fci.direct_spin1.make_rdm12 | 3.1e-15 | 3.1e-15 (abs) | 1.0e-12 | PASS |
| I4b | RDM builder 'bitwise', CAS(8,8) FCI vector | N2/6-31G | PySCF fci.direct_spin1.make_rdm12 | 1.3e-14 | 1.3e-14 (abs) | 1.0e-12 | PASS |
| I4s | RDM builder 'sparse', CAS(8,8) FCI vector | N2/6-31G | PySCF fci.direct_spin1.make_rdm12 | 1.3e-14 | 1.3e-14 (abs) | 1.0e-12 | PASS |
| I5 | two-electron term, FCI root 2 via th.load(fci.FCI) | water/STO-3G | same kernel on PySCF's root-2 RDMs | 1.1e-13 | 1.5e-15 (rel) | 1.0e-11 | PASS |
| I6 | total, DMRG-CASCI(6,6) (block2, M = 500), factored | water/cc-pVDZ | PySCF exact CASCI Gamma + ft + Lebedev | 2.9e-08 | 3.4e-10 (rel) | 1.0e-08 | PASS |

- **I1**: max element difference
- **I2**: max element difference
- **I5**: checks the loader (orbitals, root selection) with the kernel held fixed
- **I6**: DMRG - CASCI: energy 5.7e-14 Eh, active 2RDM 2.6e-07 (DMRG convergence of the density, not THISTLE)

## J. Readers (numbers printed by Molpro and OpenMolcas)

| ID | Quantity | System | Independent reference | abs | rel / abs | Tolerance | Result |
|---|---|---|---|---|---|---|---|
| J1 | Molpro MULTI: 4 state energies from our RDMs | LiF SA4-CASSCF(6,4) | energies printed by Molpro | 4.2e-09 | 4.2e-09 (abs) | 2.0e-07 | PASS |
| J2 | OpenMolcas MRCI: energy from our RDMs (12626 CSFs) | water/cc-pVDZ | CI energy printed by OpenMolcas | 7.9e-08 | 7.9e-08 (abs) | 3.0e-07 | PASS |
| J4 | OpenMolcas RASSCF: active natural occupations from our 1RDM | C7H8/ANO-L-MB CASSCF(8,8) | occupations printed by OpenMolcas | 3.4e-06 | 3.4e-06 (abs) | 1.0e-05 | PASS |
| J3 | Molden reader: elastic scattering | CO UHF (PySCF Molden file) | PySCF molden.load + ft_aopair + Lebedev | 6.8e-11 | 5.5e-13 (rel) | 1.0e-09 | PASS |

- **J1**: Hartree; limited by the 8 printed decimals
- **J2**: limited by the 6 printed decimals of the CSF coefficients
- **J4**: CI vector of 1764 CSFs -> determinants; printed to 4 decimals. Energy not usable: Molcas ran with RICD integrals

## K. Other implementations

| ID | Quantity | System | Independent reference | abs | rel / abs | Tolerance | Result |
|---|---|---|---|---|---|---|---|
| K1co_631 | earlier production kernel, HF, its own cut-off | co_631gs | stored output of the earlier code | 5.0e-12 | 2.7e-14 (rel) | 1.0e-09 | PASS |
| K2co_631 | earlier production kernel, CCSD 2-e, its own cut-off | co_631gs | stored output of the earlier code | 1.8e-09 | 1.0e-11 (rel) | 1.0e-09 | PASS |
| K1h2o_63 | earlier production kernel, HF, its own cut-off | h2o_631g | stored output of the earlier code | 2.0e-12 | 2.2e-14 (rel) | 1.0e-09 | PASS |
| K2h2o_63 | earlier production kernel, CCSD 2-e, its own cut-off | h2o_631g | stored output of the earlier code | 1.8e-12 | 2.0e-14 (rel) | 1.0e-09 | PASS |
| K1h2o_cc | earlier production kernel, HF, its own cut-off | h2o_ccpvdz | stored output of the earlier code | 2.5e-12 | 2.8e-14 (rel) | 1.0e-09 | PASS |
| K2h2o_cc | earlier production kernel, CCSD 2-e, its own cut-off | h2o_ccpvdz | stored output of the earlier code | 4.9e-11 | 5.4e-13 (rel) | 1.0e-09 | PASS |
| K1h2o_st | earlier production kernel, HF, its own cut-off | h2o_sto3g | stored output of the earlier code | 1.4e-12 | 1.5e-14 (rel) | 1.0e-09 | PASS |
| K2h2o_st | earlier production kernel, CCSD 2-e, its own cut-off | h2o_sto3g | stored output of the earlier code | 9.8e-13 | 1.1e-14 (rel) | 1.0e-09 | PASS |
| K3total | thistle-dev Fortran engine (MO deck) | CH4/STO-3G | thistle-dev golden output (ifort engine) | 1.1e-13 | 1.3e-15 (rel) | 1.0e-11 | PASS |
| K3zcotr | thistle-dev Fortran engine (AO deck) | CH4/STO-3G | thistle-dev golden output (ifort engine) | 6.5e-13 | 7.3e-15 (rel) | 1.0e-11 | PASS |
| K4 | earlier j2 code (Y20 convention), its own cut-off | water/6-31G | stored output (E. Moen Liane) | 5.8e-09 | 8.3e-08 (rel) | 1.0e-06 | PASS |
- **Three programs** (stored run of `paper/fig_three_programs.py`, not repeated here): OpenMolcas RASSCF read by THISTLE against PySCF CASCI in the same orbitals, water CASSCF(6,5)/cc-pVDZ: total 1.6e-08, elastic 3.7e-08, inelastic 2.9e-07, S2 4.0e-06 (limited by the 8 decimals of the Molden file).

## What the checks found besides the numbers

- **Open-shell SCF objects** are not accepted by `th.load` (`from_scf` handles closed shells only). The hydrogen-atom checks (F1–F3) therefore pass the spin-summed UHF density through `from_density`.
- **Two values of the Bohr radius** are in use: CODATA 2018 (0.529177210903 Å) in `thistle/core.py` and `grids.py`, CODATA 2010 (0.52917721092 Å) in `iam.py` and the Molden reader. They differ by 1.7×10⁻¹¹ relative; harmless, but worth unifying.
- **The OpenMolcas C₇H₈ test file was computed with RICD integrals.** Its printed energy cannot be rebuilt with exact integrals: an exact PySCF CASCI in the same orbitals is 0.024 Eh lower, and our densities reproduce that PySCF energy to 5×10⁻⁶. The reader is checked instead on the natural occupations Molcas prints (J4).
- **Large-q behaviour.** In a Gaussian basis the two-electron term does not vanish at q = 60–80 bohr⁻¹ (3.5×10⁻⁵ for water): tight core pairs keep e^{−q²/4p} near 1. So I_tot → N is not a usable test there; B4 checks the actual values at q = 8–16 bohr⁻¹ against PySCF instead.
- **DMRG.** The block2 density differs from exact CASCI by 2.6×10⁻⁷ (convergence of the density, although the energies agree to 10⁻¹³); THISTLE carries that into the scattering at the 10⁻¹⁰ level.

## Rerun

```
python validation/independent/run_independent_checks.py   # about 5 minutes on 4 cores
python validation/independent/make_report.py              # -> results.md
```
