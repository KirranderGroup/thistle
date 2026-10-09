# MRCI and THISTLE: what can and cannot be extracted

**The CI vector of an internally contracted MRCI cannot be exported as determinants or CSFs.**
Molpro's MRCI is internally contracted: the singly and doubly external parts of the wavefunction
are expanded in *contracted pair functions* (products of the reference states with pair
excitation operators), not in individual configurations.  `print,ci=1` therefore prints the
reference-space (internal) coefficients and the expansion coefficients of those contracted
functions; there is no complete list of determinants with coefficients anywhere in the output,
and `save,civec` stores the vector in the same contracted representation.  That is why the
extraction fails: the object does not exist in the form a CI-to-RDM builder needs.

What Molpro MRCI does provide, and THISTLE can use today:

| Quantity | How | Gives |
|---|---|---|
| MRCI one-particle density of each state | `dm,record` + `natorb` + `put,molden` with `state=` | elastic j0, elastic j2, detector images at the MRCI level (`th.load("pyridine_mrci_no_state1.mld")`, orbitals + occupations) |
| SA-CASSCF 1RDM and 2RDM of each state | `config,det` + `gprint,civector` (or the punch file) | total, inelastic, j2 total: the correlated 2RDM at the CASSCF level |
| MRCI reference-space coefficients | `print,ref=1` | a renormalised CAS vector; not an improvement over CASSCF for the 2RDM |

So the MRCI "endpoint" of the paper is: elastic quantities with the MRCI density, total quantities
with the CASSCF 2RDM.  Molpro has no MRCI 2RDM output (it is never formed: no MRCI analytic
gradients in Molpro).

**Routes to a correlated multireference 2RDM that THISTLE can read**

1. **Selected CI in a large active space (PySCF, `selected_ci` / heat-bath CI; or Dice/SHCI).**
   An explicit determinant expansion, uncontracted, with 1e5-1e7 determinants in 20-40 orbitals:
   exactly what the bitwise RDM builder is for.  `pyridine_sci.py` does this with PySCF and feeds
   the RDMs to THISTLE through the factored embedding; tested on N2 (18 orbitals, 2.7e5
   determinants: RDMs in 26 s, agreeing with PySCF's own to 1e-13).  This is the practical
   "uncontracted MRCI" and the method endpoint I would use.
2. **Uncontracted MRCI in OpenMolcas** (`molcas_test/`): done.  `THRPrint = 0.0` prints the complete
   CSF vector; `th.load("run.log")` reads it (reader `interfaces.molcas.from_mrci`), verified on
   water/cc-pVDZ (energy from our RDMs = Molcas's to 8e-8).  Limit for now: at most 64 correlated
   orbitals (inactive + active + virtual), i.e. small basis sets; a large-basis MRCI needs a
   dedicated RDM builder (to do).  ORCA's uncontracted MRCI (`PrintWF det`) would need its own reader.
3. **CASPT2/NEVPT2 densities**: PySCF's NEVPT2 and Molpro's RS2 give at most the 1RDM
   (and RS2 only with the gradient code); no 2RDM.

Run the Molpro input with `molpro pyridine_mrci.inp`; the geometry is a standard pyridine
geometry in Angstrom (check it).  The MULTI part is what `th.load("pyridine_mrci.out", state=2)`
reads (CASSCF, state 2 = the n-pi* state); the MRCI Molden files give the elastic MRCI signals.
