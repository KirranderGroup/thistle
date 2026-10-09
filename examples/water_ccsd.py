"""Total, elastic and inelastic scattering of water at the CCSD level, and its split into
mean-field and correlation (cumulant) contributions.

    python examples/water_ccsd.py
"""
import numpy as np
from pyscf import gto, scf, cc

import pyxscat as xs

mol = gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="cc-pvdz", verbose=0)
mf = scf.RHF(mol).run()
mycc = cc.CCSD(mf).run()

q = np.linspace(0.0, 10.0, 101)                      # inverse bohr
wf = xs.interfaces.pyscf.from_cc(mycc)          # factored=True keeps the 2RDM as products of
                                                # one-particle matrices (default beyond 60 functions)
parts = xs.decompose(wf, q)

hf_total = xs.scattering(xs.interfaces.pyscf.from_scf(mf), q, "hf_total")

np.savetxt("water_ccsd.dat",
           np.column_stack([q, parts["total"], parts["elastic"], parts["inelastic"],
                            parts["hartree"], parts["exchange"], parts["cumulant"], hf_total]),
           header="q  total  elastic  inelastic  hartree  exchange  cumulant  HF_total")

print(f"I_total(0) = {parts['total'][0]:.6f}   (N^2 = {mol.nelectron**2})")
print(f"largest correlation effect on the total: {np.abs(parts['total'] - hf_total).max():.4f} "
      f"at q = {q[np.abs(parts['total'] - hf_total).argmax()]:.2f} a.u.")
print("curves written to water_ccsd.dat")
