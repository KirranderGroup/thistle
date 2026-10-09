# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""The earlier production code (BigMemoryScattering, one-loop dense-AO kernel), driven by the
cumulant_scattering package, against THISTLE on the same wavefunctions.

For each wavefunction, cumulant_scattering decomposes the MO 2RDM into Total, Hartree (gamma x gamma),
exchange (gamma x_x gamma) and cumulant (Lambda), writes the inputs and the AO 2RDM (2rdmAO) for
each, and BigMemoryScattering's Main_opt.exe (built with gfortran from its own recipe) computes the
two-electron sum.  THISTLE computes the same parts from the PySCF object.

python run_bigmemory.py BIGMEMORY_EXE CUMULANT_SRC   (build: build_old_codes.sh)   -> bigmemory_results.txt
"""
import os, sys, shutil, subprocess, tempfile, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "/tmp/oldcodes/bigmem/Main_opt.exe")
CUM = os.path.abspath(sys.argv[2] if len(sys.argv) > 2 else "/tmp/oldcodes/cumscat/src")
sys.path.insert(0, CUM)
sys.path.insert(0, os.path.join(HERE, "..", "..", "python"))

import pyxscat as xs
import thistle as th
from pyxscat.interfaces import pyscf as xp
from pyscf import gto, scf, cc, mcscf
from cumulant_scattering.scattering.runner import run_scattering
from cumulant_scattering.rdm.decomposition import decompose_2rdm

PIFAC = (float(np.arccos(np.float32(-1.0))) / np.pi) ** 3     # single-precision pi of the old kernel
QMIN, QMAX, NQ = 1e-10, 10.0, 41
q = np.linspace(QMIN, QMAX, NQ)
OUT = []


def log(s=""):
    print(s, flush=True); OUT.append(s)


def bigmemory(obj, mol, gamma_mo):
    """Run cumulant_scattering -> BigMemoryScattering in a scratch directory; return I_2(q)."""
    d = tempfile.mkdtemp(prefix="bigmem_"); cwd = os.getcwd(); os.chdir(d)
    try:
        t0 = time.time()
        run_scattering(obj, mol, gamma_mo, "out.dat", qmin=QMIN, qmax=QMAX, npoints=NQ, Type=1)
        r = subprocess.run([EXE], capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists("out.dat"):
            raise RuntimeError(f"BigMemoryScattering failed: {r.stdout[-400:]} {r.stderr[-400:]}")
        a = np.loadtxt("out.dat")
        return a[:, 0], a[:, 1] / PIFAC, time.time() - t0
    finally:
        os.chdir(cwd); shutil.rmtree(d, ignore_errors=True)


def ours_parts(wf, cutoffcentre):
    """THISTLE's two-electron parts with the dense kernel on the same (8-fold averaged) densities."""
    kw = dict(screen_tol=0.0, cutoffcentre=cutoffcentre)
    tot = xs.two_electron(wf.basis, q, wf.dm2, **kw)
    hart = xs.elastic(wf.basis, q, wf.dm1, method="quartets", **kw)
    exch = xs.exchange(wf.basis, q, wf.dm1, **kw)
    return {"Total": tot, "H": hart, "X": exch, "Lambda": tot - hart - exch}


WATER = "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"
cases = []
mol = gto.M(atom=WATER, basis="6-31g", cart=True, verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12); mycc = cc.CCSD(mf).run(conv_tol=1e-10)
cases.append(("water/6-31G CCSD", mol, mf, mycc.make_rdm1(), mycc.make_rdm2(), xp.from_cc(mycc, factored=False)))
molc = gto.M(atom=WATER, basis="cc-pvdz", cart=True, verbose=0)
mfc = scf.RHF(molc).run(conv_tol=1e-12); mcc = mcscf.CASSCF(mfc, 6, 6); mcc.conv_tol = 1e-11; mcc.kernel()
from pyscf.mcscf import addons
d1a, d2a = mcc.fcisolver.make_rdm12(mcc.ci, 6, 6)
g1, g2 = addons._make_rdm12_on_mo(d1a, d2a, mcc.ncore, 6, mcc.mo_coeff.shape[1])
cases.append(("water/cc-pVDZ CASSCF(6,6)", molc, mcc, g1, g2, xp.from_casscf(mcc, factored=False)))
molo = gto.M(atom="C 0 0 0; O 0 0 1.128", basis="cc-pvdz", cart=True, verbose=0)
mfo = scf.RHF(molo).run(conv_tol=1e-12); cco = cc.CCSD(mfo).run(conv_tol=1e-10)
cases.append(("CO/cc-pVDZ CCSD", molo, mfo, cco.make_rdm1(), cco.make_rdm2(), xp.from_cc(cco, factored=False)))

log(f"BigMemoryScattering (driven by cumulant_scattering) vs THISTLE; q = {QMIN:g}..{QMAX:g} bohr^-1, {NQ} points")
log("max |old - THISTLE| / max |THISTLE|, per part of the two-electron term")
log(f"{'case':28s} {'part':7s} {'cut-off 0.01':>13s} {'exact':>10s}  {'I(0) old':>12s} {'time':>7s}")
for name, m, obj, gamma, Gamma, wf in cases:
    parts = decompose_2rdm(gamma, Gamma)
    mine01 = ours_parts(wf, 1e-2)
    mineex = ours_parts(wf, 1e-10)
    for key in ("Total", "H", "X", "Lambda"):
        G = parts[key]
        if key == "X":
            G = -0.5 * G                                     # exchange enters as -1/2 gamma x_x gamma
        qq, old, dt = bigmemory(obj, m, G)
        assert np.allclose(qq, q)
        r01 = np.max(np.abs(old - mine01[key])) / np.max(np.abs(mine01[key]))
        rex = np.max(np.abs(old - mineex[key])) / np.max(np.abs(mineex[key]))
        log(f"{name:28s} {key:7s} {r01:13.1e} {rex:10.1e}  {old[0]:12.6f} {dt:6.1f}s")
open(os.path.join(HERE, "bigmemory_results.txt"), "w").write("\n".join(OUT) + "\n")
