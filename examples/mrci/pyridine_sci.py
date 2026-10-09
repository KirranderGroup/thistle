# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Pyridine: selected CI (heat-bath style, PySCF) in a large active space as the uncontracted
multireference endpoint, and its scattering through THISTLE's bitwise RDM builder and the
factored embedding.  Usage:  python pyridine_sci.py [basis] [norb] [cutoff]
(defaults cc-pvdz, 20 active orbitals around the Fermi level, select_cutoff 1e-3; a few minutes)."""
import sys, time
import numpy as np
from pyscf import gto, scf, mcscf
from pyscf.fci import selected_ci
import pyxscat as xs
from pyxscat import ci
from pyxscat.interfaces import pyscf as xp

basis = sys.argv[1] if len(sys.argv) > 1 else "cc-pvdz"
norb = int(sys.argv[2]) if len(sys.argv) > 2 else 20
cutoff = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
GEOM = """N 0 0 1.4223; C 0 1.1414 0.7179; C 0 -1.1414 0.7179; C 0 1.1955 -0.6775; C 0 -1.1955 -0.6775;
C 0 0 -1.3809; H 0 2.0588 1.3064; H 0 -2.0588 1.3064; H 0 2.1544 -1.2131; H 0 -2.1544 -1.2131; H 0 0 -2.4665"""
mol = gto.M(atom=GEOM, basis=basis, verbose=0)
mf = scf.RHF(mol).run()
nel = 2 * (norb // 2)                                 # active electrons: fill half the active orbitals
na = nb = nel // 2
mc = mcscf.CASCI(mf, norb, nel)
mc.fcisolver = selected_ci.SCI(mol)
mc.fcisolver.select_cutoff = cutoff
mc.fcisolver.ci_coeff_cutoff = cutoff
t = time.time(); mc.kernel()
print(f"SCI/{basis}: CAS({nel},{norb}), E = {mc.e_tot:.8f}, {time.time() - t:.0f} s", flush=True)
sa, sb = mc.ci._strs
dets = [(int(a), int(b)) for a in sa for b in sb]
coeffs = np.asarray(mc.ci).ravel()
print(f"{len(dets)} determinants", flush=True)
t = time.time()
d1, d2 = ci.make_rdm12(norb, na, nb, dets, coeffs, method="bitwise")      # the Fortran pair loop
print(f"RDMs (bitwise builder): {time.time() - t:.0f} s, trace {np.trace(d1):.6f}", flush=True)
D1, F = ci.embed_core(d1, d2, mc.ncore, factored=True)                   # core + active, no n^4 array
wf = xp.from_density(mol, D1, F, mo_coeff=mc.mo_coeff[:, :mc.ncore + norb], label=f"SCI({nel},{norb})")
q = np.linspace(0.0, 8.0, 41) * 0.529177210903                           # 0-8 1/A in bohr^-1
out = np.column_stack([q, xs.scattering(wf, q), xs.scattering(wf, q, "elastic"),
                       np.r_[0.0, xs.scattering(wf, q[1:], l=2)]])
np.savetxt(f"pyridine_sci_{basis}_{norb}.dat", out, header="q(bohr^-1) I_total I_elastic S2_total")
print("I(0) =", out[0, 1], "(N^2 = 1764)")
