# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Wall time and memory of the two Python interfaces on benzene CASSCF(6,6), total and elastic.
scat_lib: AO deck ('zcotr', the route that scales) -> thistle-dev Fortran engine (subprocess).
THISTLE:  th.load(mc) -> kernels in process (factored two-particle density beyond 60 functions).
python scaling.py THISTLE_DEV_CHECKOUT BASIS   (one basis per process, so peak RSS is per basis)"""
import os, sys, time, resource, tempfile, shutil
import numpy as np
DEV = os.path.abspath(sys.argv[1]); BASIS = sys.argv[2]; WHO = sys.argv[3]
os.environ["SCAT_DIR"] = os.path.join(DEV, "fortran")
sys.path.insert(0, os.path.join(DEV, "python"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python"))
from pyscf import gto, scf, mcscf

BENZ = """C 0.0000 1.3970 0; C 1.2098 0.6985 0; C 1.2098 -0.6985 0; C 0.0000 -1.3970 0; C -1.2098 -0.6985 0;
C -1.2098 0.6985 0; H 0.0000 2.4810 0; H 2.1486 1.2405 0; H 2.1486 -1.2405 0; H 0.0000 -2.4810 0;
H -2.1486 -1.2405 0; H -2.1486 1.2405 0"""
mol = gto.M(atom=BENZ, basis=BASIS, cart=True, verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-10)
mc = mcscf.CASSCF(mf, 6, 6); mc.kernel()
Q = (1e-10, 8.0, 81)
rss0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
for kind in ("elastic", "total"):
    t0 = time.time()
    if WHO == "scat_lib":
        from scat_lib.pyscf_scat import scat_calc as sc
        d = tempfile.mkdtemp(); cwd = os.getcwd(); os.chdir(d)
        try:
            out = sc.run_scattering_pyscf(mc, mf, "scat", backend="zcotr", type=kind, q_range=Q[:2], q_points=Q[2])
            ok = "ok"
        except Exception as e:
            ok = f"failed ({type(e).__name__})"
        os.chdir(cwd); shutil.rmtree(d, ignore_errors=True)
        child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1e6
    else:
        import thistle as th
        m = th.load(mc)
        r = m.scattering(q=np.linspace(*Q), unit="1/bohr", parts="all" if kind == "total" else "standard")
        ok, child = "ok", 0.0
    dt = time.time() - t0
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    print(f"{WHO:9s} {BASIS:8s} n={mol.nao:4d} {kind:8s} {dt:7.1f} s  peak RSS python {rss:5.2f} GB, engine {child:5.2f} GB  {ok}", flush=True)
