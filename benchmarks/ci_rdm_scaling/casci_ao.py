# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""CASCI(6,6) of benzene: time and peak memory of each stage of a total-scattering run,
for the three ways of handing the 2RDM to the kernels.

    python casci_ao.py prepare <basis>            # RHF + CASCI, saved to benzene_<basis>.npz
    python casci_ao.py run <basis> <route>        # route: pyscf | dense | factored

``pyscf``    mcscf.addons.make_rdm12 (all MOs, n_MO^4, then AO n_AO^4), as before;
``dense``    our embedding over core + active orbitals only, then AO n_AO^4;
``factored`` FactoredDM2 (no four-index array anywhere), kernel mode 3.
"""
import json, resource, sys, time
import numpy as np
from pyscf import gto, scf, mcscf

GEOM = """C 0 1.3915 0; C 1.2051 0.69575 0; C 1.2051 -0.69575 0; C 0 -1.3915 0; C -1.2051 -0.69575 0;
C -1.2051 0.69575 0; H 0 2.4715 0; H 2.1404 1.23575 0; H 2.1404 -1.23575 0; H 0 -2.4715 0;
H -2.1404 -1.23575 0; H -2.1404 1.23575 0"""
Q = np.linspace(0.0, 8.0, 30)          # bohr^-1


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6      # GB


def setup(basis):
    mol = gto.M(atom=GEOM, basis=basis, verbose=0)
    mc = mcscf.CASCI(scf.RHF(mol), 6, 6)
    return mol, mc


if sys.argv[1] == "prepare":
    basis = sys.argv[2]
    mol, mc = setup(basis)
    t = time.perf_counter()
    mf = scf.RHF(mol).run()
    mc = mcscf.CASCI(mf, 6, 6)
    mc.kernel()
    np.savez(f"benzene_{basis}.npz", mo=mc.mo_coeff, ci=mc.ci, e=mc.e_tot)
    print(json.dumps({"basis": basis, "nao_sph": mol.nao, "nao_cart": mol.nao_cart(),
                      "e_casci": mc.e_tot, "seconds": round(time.perf_counter() - t, 1)}))
    sys.exit()

basis, route = sys.argv[2], sys.argv[3]
import pyxscat as px
from pyxscat.interfaces import pyscf as xp
from pyxscat.core import FactoredDM2
from pyscf.mcscf import addons
d = np.load(f"benzene_{basis}.npz")
mol, mc = setup(basis)
mc.mo_coeff, mc.ci = d["mo"], d["ci"]
out = {"basis": basis, "route": route, "nao_cart": mol.nao_cart()}
t0 = time.perf_counter()
if route == "pyscf":
    dm1, dm2 = addons.make_rdm12(mc)
    out["rdm_s"] = round(time.perf_counter() - t0, 2); out["rdm_peak_GB"] = round(rss(), 2)
    wf = xp.from_density(mol, dm1, dm2)
    del dm2
else:
    wf = xp.from_casscf(mc, factored=(route == "factored"))
out["to_AO_s"] = round(time.perf_counter() - t0, 2)
out["after_AO_peak_GB"] = round(rss(), 2)
g = wf.dm2
out["dm2_GB"] = round(g.nbytes / 1e9, 3)
if isinstance(g, FactoredDM2):
    out["terms"] = g.nterms
t = time.perf_counter()
total = px.scattering(wf, Q)
out["kernel_s"] = round(time.perf_counter() - t, 1)
out["peak_GB"] = round(rss(), 2)
out["total_q0_q1"] = [round(float(total[0]), 9), float(total[10])]
print(json.dumps(out))
