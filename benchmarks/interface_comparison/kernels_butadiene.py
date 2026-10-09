import os, sys, subprocess, numpy as np, time
sys.path.insert(0, "/home/claude/pyxscat/python"); sys.path.insert(0, "/tmp/pmain/python")
os.environ["SCAT_DIR"] = "/tmp/pmain/fortran"
import thistle as th, pyxscat as xs
from scat_lib.pyscf_scat import scat_calc as sc
from pyscf import gto, scf
PIFAC = (float(np.arccos(np.float32(-1.0))) / np.pi) ** 3
BUTA = """C -1.8235 -0.4057 0; C -0.6176 0.3486 0; C 0.6176 -0.3486 0; C 1.8235 0.4057 0;
H -2.7618 0.1374 0; H -1.8549 -1.4907 0; H -0.6286 1.4373 0; H 0.6286 -1.4373 0; H 2.7618 -0.1374 0; H 1.8549 1.4907 0"""
mol = gto.M(atom=BUTA, basis="6-31g*", cart=True, verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-10)
print("basis functions", mol.nao, flush=True)
Q = (1e-10, 15.0, 61)
for kind in ("elastic",):
    t0 = time.time()
    p = sc.write_deck_pyscf(mf, mf, "bu", backend="zcotr", type=kind, q_range=Q[:2], q_points=Q[2],
                            cutoffcentre=1e-2, out_path="bu.h5")
    subprocess.run(["/tmp/pmain/fortran/Main.exe", p], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    his = np.loadtxt("bu")[:, 1] / PIFAC; t_his = time.time() - t0
    wf = th.load(p).wf; q = np.linspace(*Q)
    t0 = time.time()
    ours = xs.elastic(wf.basis, q, wf.dm1, screen_tol=1e-16, cutoffcentre=1e-2, method="quartets"); t_q = time.time() - t0
    t0 = time.time(); ours_p = xs.elastic(wf.basis, q, wf.dm1); t_p = time.time() - t0
    r = np.abs(his - ours) / np.abs(ours)
    print(f"butadiene {kind}: his vs ours (same 0.01 cut-off): " + " ".join(
        f"q={qq:.0f}:{r[np.argmin(abs(q-qq))]:.0e}" for qq in (2, 4, 6, 8, 10, 12, 15)), flush=True)
    r2 = np.abs(ours_p - ours) / np.abs(ours_p)
    print(f"  ours exact (pairs) vs ours with the 0.01 cut-off: " + " ".join(
        f"q={qq:.0f}:{r2[np.argmin(abs(q-qq))]:.0e}" for qq in (2, 4, 6, 8, 10, 12, 15)))
    print(f"  time: his deck+engine {t_his:.1f} s; ours quartets {t_q:.1f} s, ours pairs {t_p:.1f} s")
