import os, sys, subprocess, numpy as np, time
sys.path.insert(0, "/home/claude/pyxscat/python"); sys.path.insert(0, "/tmp/pmain/python")
os.environ["SCAT_DIR"] = "/tmp/pmain/fortran"
import thistle as th, pyxscat as xs
from scat_lib.pyscf_scat import scat_calc as sc
from pyscf import gto, scf, mcscf
BENZ = """C 0.0000 1.3970 0; C 1.2098 0.6985 0; C 1.2098 -0.6985 0; C 0.0000 -1.3970 0; C -1.2098 -0.6985 0;
C -1.2098 0.6985 0; H 0.0000 2.4810 0; H 2.1486 1.2405 0; H 2.1486 -1.2405 0; H 0.0000 -2.4810 0;
H -2.1486 -1.2405 0; H -2.1486 1.2405 0"""
mol = gto.M(atom=BENZ, basis="6-31g", cart=True, verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-10); mc = mcscf.CASSCF(mf, 6, 6).run()
Q = (1e-10, 8.0, 81); q = np.linspace(*Q)
for kind in ("total", "elastic"):
    p = sc.write_deck_pyscf(mc, mf, "bz", backend="zcotr", type=kind, q_range=Q[:2], q_points=Q[2], cutoffcentre=1e-2, out_path=f"bz_{kind}.h5")
    wf = th.load(p).wf
    row = f"benzene/6-31G CAS(6,6) {kind:8s} 4 threads, same deck, kernel only: "
    for exe, lab in (("/tmp/pbuild4/Main.exe", "his"),):
        t0 = time.time(); subprocess.run([exe, p], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True); row += f"{lab} {time.time()-t0:.1f} s; "
    if kind == "total":
        for lab, kw in (("ours (dense, cut-off 0.01)", dict(cutoffcentre=1e-2)), ("ours (dense, exact)", {})):
            t0 = time.time(); xs.two_electron(wf.basis, q, wf.dm2, **kw); row += f"{lab} {time.time()-t0:.1f} s; "
        m = th.load(mc); t0 = time.time(); xs.two_electron(m.wf.basis, q, m.wf.dm2); row += f"ours (factored, exact) {time.time()-t0:.1f} s"
    else:
        for lab, kw in (("ours (quartets, 0.01)", dict(cutoffcentre=1e-2, method="quartets")), ("ours (pairs, exact)", {})):
            t0 = time.time(); xs.elastic(wf.basis, q, wf.dm1, **kw); row += f"{lab} {time.time()-t0:.1f} s; "
    print(row, flush=True)
