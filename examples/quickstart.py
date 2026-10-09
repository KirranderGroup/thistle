"""THISTLE in a few lines: water at the CCSD level, compared with the independent atom model,
the detector image of one fixed orientation, and the same for electron diffraction.

    python examples/quickstart.py
"""
from pyscf import gto, scf, cc

import thistle as th

# 1. any electronic-structure calculation (here PySCF; a Molden file or a Molcas/Molpro
#    output works the same way: th.load("lif.out", state=2))
mol = gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="aug-cc-pvdz")
water = th.load(cc.CCSD(scf.RHF(mol).run()).run(), name="water CCSD")
print(water)

# 2. rotationally averaged scattering (q in 1/Angstrom unless unit="1/bohr")
res = water.scattering(qmax=4, parts="all")
print(res)
res.save("water_ccsd.csv")
res.plot("water_ccsd.png", curves=["total", "elastic", "inelastic"])

# 3. the same for free atoms, and the difference in percent
th.compare(res, water.iam(qmax=4, basis="aug-cc-pvdz")).plot("water_vs_iam.png")

# 4. a molecule held fixed in the beam (1 Angstrom X-rays, molecular plane facing the beam)
water.detector(wavelength=1.0, pixels=81, orientation=(0, 90, 0)).plot("water_detector.png")

# 5. electron diffraction: the same calls with probe="electron" (3 MeV beam)
e = water.scattering(qmax=10, probe="electron", energy_kev=3000)
e.plot("water_ued.png")                                       # dsigma/dOmega, log scale
th.sm(e, water.iam(qmax=10, probe="electron", basis="aug-cc-pvdz")).plot("water_sm.png")
