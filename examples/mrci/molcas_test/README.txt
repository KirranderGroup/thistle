OpenMolcas MRCI test for THISTLE (water/cc-pVDZ, CASSCF(6,5) -> uncontracted MRCI).

On the cluster (OpenMolcas available):
    mkdir mrci_water && cd mrci_water && cp .../water_mrci.input .
    pymolcas water_mrci.input > water_mrci.log 2>&1

Send back (to ~/Downloads on the Mac, or anywhere I can read):
    water_mrci.log                 the complete output (the MRCI CSF printout is what I need)
    water_mrci.rasscf.molden       the RASSCF orbitals (basis of the CSF expansion, after MOTRA)
    water_mrci.rasscf.h5           if written (orbitals and CI vector in HDF5: a cleaner route)
    any CIORB01 / MRCI natural-orbital file the run leaves in the work directory

If MRCI stops because the CI space is too large for the default memory, add
    > export MOLCAS_MEM=4000
before pymolcas.  The expected size here is a few 10^4 CSFs.
