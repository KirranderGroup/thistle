import os, sys, subprocess, numpy as np
sys.path.insert(0, "/home/claude/pyxscat/python"); sys.path.insert(0, "/tmp/pmain/python")
os.environ["SCAT_DIR"] = "/tmp/pmain/fortran"
import thistle as th, pyxscat as xs
from scat_lib.pyscf_scat import scat_calc as sc
from pyscf import gto, scf, mcscf
PIFAC = (float(np.arccos(np.float32(-1.0))) / np.pi) ** 3
mol = gto.M(atom="C 0 0 0; O 0 0 1.128", basis="cc-pvdz", cart=True, verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-11); mc = mcscf.CASSCF(mf, 6, 6).run()
for kind in ("elastic", "total"):
    for qmax, npts in ((50.0, 101), (250.0, 251)):
        p = sc.write_deck_pyscf(mc, mf, "co", backend="zcotr", type=kind, q_range=(1e-10, qmax), q_points=npts,
                                cutoffcentre=1e-2, out_path=f"co_{kind}.h5")
        subprocess.run(["/tmp/pmain/fortran/Main.exe", p], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        his = np.loadtxt("co")[:, 1] / PIFAC
        wf = th.load(p).wf; q = np.linspace(1e-10, qmax, npts)
        f = xs.two_electron if kind == "total" else (lambda b, q, d, **k: xs.elastic(b, q, d, method="quartets", **k))
        ours = f(wf.basis, q, wf.dm2 if kind == "total" else wf.dm1, screen_tol=1e-16, cutoffcentre=1e-2)
        r = np.abs(his - ours) / np.abs(ours)
        print(f"CO {kind:8s} q<= {qmax:5.0f}: max rel diff {r.max():.1e} at q = {q[np.argmax(r)]:.1f}; "
              + " ".join(f"q={qq:.0f}:{r[np.argmin(abs(q-qq))]:.0e}" for qq in (10, 30, 60, 100, 150, 200, 250) if qq <= qmax), flush=True)
