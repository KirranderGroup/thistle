"""Molecular scattering relative to the independent atom model, for a Molpro SA-CASSCF state of
LiF and for water at the CCSD level.

    python examples/compare_with_iam.py
"""
import numpy as np
from pyscf import gto, scf, cc

import pyxscat as xs
from pyxscat import iam
from pyxscat.interfaces import molpro, pyscf as xp

q = np.linspace(0.05, 8.0, 80)

# LiF, Molpro state-averaged CASSCF (test data shipped with the code)
wf = molpro.from_multi("tests/data/molpro/lif_sa4_casscf.out", "tests/data/molpro/lif_sa4_casscf.mld", state=1)
ref = iam.for_wavefunction(wf, q, basis={"F": "aug-cc-pvtz", "Li": "cc-pvtz"})
lif = iam.percent_difference(xs.scattering(wf, q, "total"), ref["total"])

# water, CCSD
mol = gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="aug-cc-pvtz", verbose=0)
wf = xp.from_cc(cc.CCSD(scf.RHF(mol).run()).run())
ref = iam.for_wavefunction(wf, q, basis="aug-cc-pvtz")
h2o = iam.percent_difference(xs.scattering(wf, q, "total"), ref["total"])

np.savetxt("percent_difference_from_iam.dat", np.column_stack([q, lif, h2o]), header="q  LiF_state1  H2O_CCSD  (%)")
print("written percent_difference_from_iam.dat")
