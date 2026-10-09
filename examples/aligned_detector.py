"""Detector image of a fixed-orientation water molecule (CCSD): total, elastic and the
correlation (cumulant) contribution.

    python examples/aligned_detector.py
"""
import numpy as np
from pyscf import gto, scf, cc

import pyxscat as xs
from pyxscat import aligned

mol = gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692", basis="cc-pvdz", verbose=0)
mycc = cc.CCSD(scf.RHF(mol).run()).run()
wf = xs.interfaces.pyscf.from_cc(mycc)

k0 = 2 * np.pi / (1.0 / 0.529177)                 # 1 Angstrom photons, in bohr^-1
q_lab, on_detector = aligned.detector_q(k0, qmax=3.0, npix=81)
R = aligned.rotation_matrix(0.0, np.pi / 2, 0.0)  # molecular plane turned to face the beam
parts = aligned.decompose(wf, aligned.rotate_q(q_lab[on_detector], R))

np.savez("water_detector.npz", q_lab=q_lab, on_detector=on_detector,
         **{k: v for k, v in parts.items() if k != "qvec"})
print("images written to water_detector.npz:", ", ".join(k for k in parts if k != "qvec"))
