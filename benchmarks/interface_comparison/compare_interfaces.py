# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Python interfaces of thistle-dev (scat_lib, Patrick Wang) and THISTLE, on the same wavefunctions.

For every PySCF wavefunction below, three routes:
  A  scat_lib.run_scattering_pyscf  -> HDF5 deck -> thistle-dev Fortran engine (Main.exe)
  B  scat_lib.write_deck_pyscf      -> the same deck -> THISTLE kernels (th.load("deck.h5"))
  C  th.load(pyscf_object)          -> THISTLE kernels, in memory
A vs B: the kernels on identical input.  B vs C: the two interfaces (RDM extraction, orbitals,
basis handling) into the same kernels.  Both deck representations (MO 'normal', AO 'zcotr'),
both observables (total two-electron term, elastic).

python compare_interfaces.py THISTLE_DEV_CHECKOUT   (Main.exe built in its fortran/)
"""
import os, sys, time, tempfile, shutil, contextlib
import numpy as np

DEV = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "/tmp/pmain")
os.environ["SCAT_DIR"] = os.path.join(DEV, "fortran")
sys.path.insert(0, os.path.join(DEV, "python"))
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "python"))

import thistle as th
import pyxscat as xs
from scat_lib.pyscf_scat import scat_calc as sc
from pyscf import gto, scf, mcscf, cc, fci

PIFAC = (float(np.arccos(np.float32(-1.0))) / np.pi) ** 3    # single-precision pi of the engine
Q = (1e-10, 8.0, 81)                                           # bohr^-1
q = np.linspace(*Q)
OUT = []
def log(s=""):
    print(s, flush=True); OUT.append(s)


@contextlib.contextmanager
def scratch():
    d = tempfile.mkdtemp(prefix="ifc_"); cwd = os.getcwd(); os.chdir(d)
    try:
        yield d
    finally:
        os.chdir(cwd); shutil.rmtree(d, ignore_errors=True)


def rel(a, b):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


def ours(wf, kind, cutoffcentre):
    if kind == "total":
        return xs.two_electron(wf.basis, q, wf.dm2, screen_tol=1e-16, cutoffcentre=cutoffcentre)
    return xs.elastic(wf.basis, q, wf.dm1, screen_tol=1e-16, cutoffcentre=cutoffcentre, method="quartets")


# ------------------------------------------------------------------ a minimal block2 DMRG solver
class DMRGCI:
    """Enough of pyscf.dmrgscf.DMRGCI for CASCI: kernel() keeps the MPS, make_rdm12(state, ...)
    returns spin-summed RDMs in PySCF's convention (what both interfaces call)."""
    def __init__(self, mol, bond_dim=400, nroots=1):
        self.mol, self.M, self.nroots, self.spin = mol, bond_dim, nroots, 0
        self.stdout, self.verbose = mol.stdout, 0
    def kernel(self, h1, h2, norb, nelec, ci0=None, ecore=0, **kw):
        from pyblock2.driver.core import DMRGDriver, SymmetryTypes
        from pyscf import ao2mo
        ne = sum(nelec) if isinstance(nelec, (tuple, list)) else nelec
        d = DMRGDriver(scratch=tempfile.mkdtemp(prefix="b2_"), symm_type=SymmetryTypes.SU2, n_threads=4)
        d.initialize_system(n_sites=norb, n_elec=ne, spin=0)
        g2 = ao2mo.restore(1, h2, norb)
        mpo = d.get_qc_mpo(h1e=h1, g2e=g2, ecore=ecore, iprint=0)
        ket = d.get_random_mps(tag="KET", bond_dim=self.M, nroots=self.nroots)
        sched = [self.M] * 4 + [self.M] * 4 + [self.M] * 12
        noise = [1e-4] * 4 + [1e-5] * 4 + [0] * 12
        thr = [1e-8] * 8 + [1e-12] * 12
        e = d.dmrg(mpo, ket, n_sweeps=20, bond_dims=sched, noises=noise, thrds=thr, iprint=0, tol=1e-12)
        self.driver, self.ket = d, ket
        return (e, 0) if self.nroots == 1 else (e, list(range(self.nroots)))
    def make_rdm12(self, state, norb, nelec, **kw):
        ket = self.ket if self.nroots == 1 else self.driver.split_mps(self.ket, state, "KET%d" % state)
        d1 = self.driver.get_1pdm(ket)
        d2 = self.driver.get_2pdm(ket).transpose(0, 3, 1, 2)     # block2 -> PySCF (pq|rs)
        return d1, d2
    def make_rdm1(self, state, norb, nelec, **kw):
        return self.make_rdm12(state, norb, nelec)[0]
    def spin_square(self, *a, **k):
        return 0.0, 1.0


# ------------------------------------------------------------------ wavefunctions
WATER = "O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692"
cases = []
mol = gto.M(atom=WATER, basis="6-31g", verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12)
mc = mcscf.CASSCF(mf, 4, 4); mc.conv_tol = 1e-11; mc.kernel()
cases.append(("H2O/6-31G CASSCF(4,4)", mc, mf, 0))
sa = mcscf.CASSCF(mf, 4, 4).state_average_([0.5, 0.5]); sa.conv_tol = 1e-11; sa.fcisolver.conv_tol = 1e-14; sa.kernel()
cases.append(("H2O/6-31G SA2-CASSCF(4,4), root 2", sa, mf, 1))
ci2 = mcscf.CASCI(mf, 4, 4); ci2.fcisolver.nroots = 2; ci2.fcisolver.conv_tol = 1e-14; ci2.kernel(sa.mo_coeff)
cases.append(("  same root as CASCI(nroots=2) on SA orbitals", ci2, mf, 1, sa))
cases.append(("H2O/6-31G CCSD", cc.CCSD(mf).run(conv_tol=1e-10), mf, 0))
dm = mcscf.CASCI(mf, 6, 6); dm.fcisolver = DMRGCI(mol, bond_dim=500); dm.kernel()
ref = mcscf.CASCI(mf, 6, 6).run()
cases.append(("H2O/6-31G DMRG-CASCI(6,6), block2", dm, mf, 0))
molc = gto.M(atom=WATER, basis="cc-pvdz", cart=True, verbose=0)
mfc = scf.RHF(molc).run(conv_tol=1e-12)
mcc = mcscf.CASSCF(mfc, 6, 6); mcc.conv_tol = 1e-11; mcc.kernel()
cases.append(("H2O/cc-pVDZ (cart) CASSCF(6,6)", mcc, mfc, 0))
mols = gto.M(atom=WATER, basis="sto-3g", verbose=0)
mfs = scf.RHF(mols).run(conv_tol=1e-12)
fs = fci.FCI(mfs); fs.kernel()
cases.append(("H2O/STO-3G FCI (fci.FCI)", fs, mfs, 0))
fs2 = fci.FCI(mfs); fs2.nroots = 3; fs2.kernel()
cases.append(("H2O/STO-3G FCI, root 3", fs2, mfs, 2))
log(f"DMRG(6,6) M=500 vs CASCI(6,6): dE = {dm.e_tot - ref.e_tot:.1e} Eh")


# ------------------------------------------------------------------ the three routes
log(f"\nq = {Q[0]:g}..{Q[1]:g} bohr^-1, {Q[2]} points; engine output divided by its single-precision pi factor")
log(f"{'case':38s} {'deck':6s} {'obs':7s} | A-B cc=0.01  A-B cc=1e-10 | B-C (exact) | t_A    t_B    t_C")
for name, obj, mfo, state, *mine in cases:
    # route C: our interface (on `mine` when the case feeds scat_lib an equivalent object)
    try:
        t0 = time.time()
        wf_c = th.load(mine[0] if mine else obj, state=state + 1, mo_coeff=mfo).wf
        c_ok = True
    except Exception as e:
        c_ok, c_err = False, f"{type(e).__name__}: {e}"
    for backend in ("normal", "zcotr"):
        for kind in ("total", "elastic"):
            row = f"{name:38s} {('MO' if backend == 'normal' else 'AO'):6s} {kind:7s} |"
            try:
                res = {}
                for cc_ in (1e-2, 1e-10):
                    with scratch():
                        t0 = time.time()
                        out = sc.run_scattering_pyscf(obj, mfo, "scat", backend=backend, type=kind, engine="fortran",
                                                      q_range=Q[:2], q_points=Q[2], cutoffcentre=cc_, state=state)
                        res[cc_] = (np.asarray(out)[:, 1] / PIFAC, time.time() - t0)
                with scratch():
                    path = sc.write_deck_pyscf(obj, mfo, "scat", backend=backend, type=kind, q_range=Q[:2],
                                               q_points=Q[2], cutoffcentre=1e-2, state=state, out_path="deck.h5")
                    wf_b = th.load(path).wf
                t0 = time.time(); b01 = ours(wf_b, kind, 1e-2); tb = time.time() - t0
                b_ex = ours(wf_b, kind, 1e-10)
                row += f"  {rel(res[1e-2][0], b01):9.1e}   {rel(res[1e-10][0], b_ex):9.1e}  |"
            except Exception as e:
                row += f"  scat_lib failed: {type(e).__name__}: {str(e)[:60]}"
                log(row); continue
            if c_ok:
                t0 = time.time(); c_ex = ours(wf_c, kind, 1e-10); tc = time.time() - t0
                row += f"  {rel(b_ex, c_ex):9.1e}  | {res[1e-2][1]:5.1f}s {tb:5.2f}s {tc:5.2f}s"
            else:
                row += f"  ours failed: {c_err[:60]}"
            log(row)

# ------------------------------------------------------------------ scat_lib's orbital_type
log("\nscat_lib orbital_type (route A, AO deck, total), against orbital_type='HF':")
base = None
for ot in ("HF", "CASSCF", "natural"):
    with scratch():
        out = np.asarray(sc.run_scattering_pyscf(mc, mf, "scat", orbital_type=ot, backend="zcotr", type="total",
                                                 q_range=Q[:2], q_points=Q[2]))[:, 1]
    base = out if base is None else base
    log(f"  {ot:8s} {rel(out, base):.1e}")

open(os.path.join(HERE, "compare_interfaces.txt"), "w").write("\n".join(OUT) + "\n")
