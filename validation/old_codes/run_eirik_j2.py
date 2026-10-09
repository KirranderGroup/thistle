# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Eirik's MO-route code (Scattering_eirik/Eirik_j2), driven by the cumulant_scattering package,
against THISTLE on the same wavefunctions.

cumulant_scattering writes basis.dat, MOs.dat and options.dat; the driver adds the MO 2RDM as the
sparse text file 2rdm.txt (i j k l value, chemists' order, 1-based) that Eirik's Main reads with
"readtwordm".  Three observables:

  Type 1  total, j0          (tot_integration)
  Type 7  total, j2          (tot_integration_j2)
  Type 8  elastic, j2        (elastic_integration_j2)

Eirik's j2 has no i^2 = -1 and is normalised with Y20 instead of P2 (both routes), so the j2
comparisons are against -sqrt(5/4pi) THISTLE.  Two j2 routes are left in a debugging state in this
copy, and the second executable ("repaired") is the same code with two lines changed:
  Type 7: the call that fills the coincident-centre table is commented out
          (main_calculation.f90: !CALL set_P0_j2(...)), so every quartet with |P - Q| < cutoffcentre
          contributes zero.  Repaired: the call restored.
  Type 8: elastic_integration_j2 in Integrals.f90 sets f = 0 after every two-centre integral, so only
          coincident-centre quartets are summed (the validated version, elastic_j2.f90, is not in the
          build).  Repaired: that line commented out.
Input reading only (all builds): Main reads three state indices, the number of 2RDM lines is counted
instead of the fixed 4549, and Type 7/8 read the 2RDM file (the switch read2rdm2 existed, set off).  THISTLE is run with the
quartet kernel and cutoffcentre = 0.01 bohr, Eirik's setting.

python run_eirik_j2.py EIRIK_EXE EIRIK_REPAIRED_EXE CUMULANT_SRC   -> eirik_j2_results.txt
"""
import os, sys, shutil, subprocess, tempfile, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "/tmp/oldcodes/eirik/Main.exe")
EXE_P0 = os.path.abspath(sys.argv[2] if len(sys.argv) > 2 else "/tmp/oldcodes/eirik_p0/Main.exe")
CUM = os.path.abspath(sys.argv[3] if len(sys.argv) > 3 else "/tmp/oldcodes/cumscat/src")
sys.path.insert(0, CUM)
sys.path.insert(0, os.path.join(HERE, "..", "..", "python"))

import pyxscat as xs
from pyxscat.interfaces import pyscf as xp
from pyscf import gto, scf, cc, mcscf
from pyscf.mcscf import addons
from cumulant_scattering.scattering.runner import run_scattering
from cumulant_scattering.rdm.decomposition import decompose_2rdm

PIFAC = (float(np.arccos(np.float32(-1.0))) / np.pi) ** 3   # acos(-1.0000) in Variables.f90, (pi/p)^(3/2) twice
Y20 = np.sqrt(5 / (4 * np.pi))
QMIN, QMAX, NQ = 1e-10, 10.0, 100                      # >= 90 points: debug writes read element 90
Q = np.linspace(QMIN, QMAX, NQ)
OUT = []


def log(s=""):
    print(s, flush=True); OUT.append(s)


def eirik(exe, obj, mol, G, typ):
    d = tempfile.mkdtemp(prefix="eirik_"); cwd = os.getcwd(); os.chdir(d)
    try:
        run_scattering(obj, mol, G, "out.dat", qmin=QMIN, qmax=QMAX, npoints=NQ, Type=typ)
        L = open("options.dat").read().split("\n")
        L[L.index("1 1")] = "1 1 1"                    # Eirik's Main reads state1, state2, state3
        open("options.dat", "w").write("\n".join(L))
        open("prueba_am.dat", "w").close()             # debug file opened with STATUS='OLD'
        idx = np.argwhere(np.abs(G) > 1e-14)
        with open("2rdm.txt", "w") as f:
            for a, b, c, e in idx:
                f.write(f"{a+1} {b+1} {c+1} {e+1} {float(G[a, b, c, e])!r}\n")
        t0 = time.time()
        r = subprocess.run([exe], capture_output=True, text=True)
        res = "j2.dat" if typ == 8 else "out.dat"    # elastic_integration_j2 writes j2.dat and stops
        if r.returncode != 0 or not os.path.exists(res):
            raise RuntimeError(f"Eirik's code failed: {r.stderr[-400:]}")
        a = np.loadtxt(res)
        return a[:, 0], a[:, 1], time.time() - t0, len(idx)
    finally:
        os.chdir(cwd); shutil.rmtree(d, ignore_errors=True)


WATER = "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"
cases = []
mol = gto.M(atom=WATER, basis="6-31g", cart=True, verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12); m = cc.CCSD(mf).run(conv_tol=1e-10)
cases.append(("water/6-31G CCSD", mol, mf, m.make_rdm2(), xp.from_cc(m, factored=False)))
molc = gto.M(atom=WATER, basis="cc-pvdz", cart=True, verbose=0)
mfc = scf.RHF(molc).run(conv_tol=1e-12); mc = mcscf.CASSCF(mfc, 6, 6); mc.conv_tol = 1e-11; mc.kernel()
d1a, d2a = mc.fcisolver.make_rdm12(mc.ci, 6, 6)
_, g2 = addons._make_rdm12_on_mo(d1a, d2a, mc.ncore, 6, mc.ncore + 6)   # core + active block
nmo = mc.mo_coeff.shape[1]
g2 = np.pad(g2, [(0, nmo - g2.shape[0])] * 4)                                # zeros for the virtuals
cases.append(("water/cc-pVDZ CASSCF(6,6)", molc, mc, g2, xp.from_casscf(mc, factored=False)))
molh = gto.M(atom="H 0 0 -2.2; H 0.3 0.2 -0.75; H -0.2 0.1 0.75; H 0.1 -0.3 2.2", basis="6-31g", cart=True, verbose=0)
mfh = scf.RHF(molh).run(conv_tol=1e-12); mh = cc.CCSD(mfh).run(conv_tol=1e-10)
cases.append(("H4 chain/6-31G CCSD (s only)", molh, mfh, mh.make_rdm2(), xp.from_cc(mh, factored=False)))

log(f"Eirik's MO-route code (driven by cumulant_scattering) vs THISTLE; q = {QMIN:g}..{QMAX:g} bohr^-1, {NQ} points")
log("max |Eirik - reference| / max |reference|; reference = THISTLE quartets, cutoffcentre 0.01, in Eirik's convention")
log("  as-is: Eirik's kernels unchanged (input reading only); repaired: the two debugging lines undone (see header)")
log(f"{'case':30s} {'observable':18s} {'as-is':>9s} {'repaired':>9s} {'2RDM lines':>11s} {'time':>7s}")
for name, mo, obj, G, wf in cases:
    kw = dict(screen_tol=0.0, cutoffcentre=1e-2)
    for typ, label in ((1, "total j0"), (7, "total j2"), (8, "elastic j2")):
        q, old, dt, nl = eirik(EXE, obj, mo, G, typ)
        old = old / PIFAC
        if typ == 1:
            ref = xs.two_electron(wf.basis, q, wf.dm2, **kw)
        elif typ == 7:
            ref = -Y20 * xs.two_electron(wf.basis, q, wf.dm2, l=2, **kw)
        else:
            ref = -Y20 * xs.elastic(wf.basis, q, wf.dm1, l=2, method="quartets", **kw)
        r0 = np.max(np.abs(old - ref)) / np.max(np.abs(ref))
        r1 = float("nan")
        if typ in (7, 8):
            _, old1, _, _ = eirik(EXE_P0, obj, mo, G, typ)
            r1 = np.max(np.abs(old1 / PIFAC - ref)) / np.max(np.abs(ref))
        log(f"{name:30s} {label:18s} {r0:9.1e} {r1:9.1e} {nl:11d} {dt:6.1f}s")

log()
log("Parts of the two-electron term (cumulant_scattering's decomposition of the same 2RDM), repaired build")
log(f"{'case':30s} {'part':7s} {'total j0':>9s} {'total j2':>9s}")
for name, mo, obj, G, wf in cases[:2]:
    kw = dict(screen_tol=0.0, cutoffcentre=1e-2)
    N = int(round((1 + np.sqrt(1 + 4 * np.einsum("ppqq->", G))) / 2))
    g1 = np.einsum("pqrr->pq", G) / (N - 1)
    parts = decompose_2rdm(g1, G)
    ours = {}
    for lq in (0, 2):
        tot = xs.two_electron(wf.basis, Q, wf.dm2, l=lq, **kw)
        hart = xs.elastic(wf.basis, Q, wf.dm1, l=lq, method="quartets", **kw)
        exch = xs.exchange(wf.basis, Q, wf.dm1, l=lq, **kw)
        ours[lq] = {"Total": tot, "H": hart, "X": exch, "Lambda": tot - hart - exch}
    for key in ("H", "X", "Lambda"):
        Gp = -0.5 * parts[key] if key == "X" else parts[key]
        row = []
        for typ, lq in ((1, 0), (7, 2)):
            qq, old, _, _ = eirik(EXE_P0, obj, mo, Gp, typ)
            assert np.allclose(qq, Q)
            old = old / PIFAC
            ref = ours[lq][key] if lq == 0 else -Y20 * ours[lq][key]
            row.append(np.max(np.abs(old - ref)) / np.max(np.abs(ref)))
        log(f"{name:30s} {key:7s} {row[0]:9.1e} {row[1]:9.1e}")
open(os.path.join(HERE, "eirik_j2_results.txt"), "w").write("\n".join(OUT) + "\n")
